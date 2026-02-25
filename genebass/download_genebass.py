#!/usr/bin/env python3
"""
Download GeneBass pLoF and missense variant data from Google Cloud Storage.

Data source: https://app.genebass.org/downloads
GCS path:    gs://ukbb-exome-public/500k/results/results.mt

Prerequisites:
  - Python 3.8+
  - Java 11 (for Hail)
  - pip install hail

Usage:
  python download_genebass.py                    # Download both pLoF and missense
  python download_genebass.py --annotation pLoF  # Download pLoF only
  python download_genebass.py --annotation missense  # Download missense only
  python download_genebass.py --output-dir ./my_output  # Custom output directory
"""

import argparse
import os
import sys
import time


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download GeneBass pLoF and missense variant results"
    )
    parser.add_argument(
        "--annotation",
        choices=["pLoF", "missense", "both"],
        default="both",
        help="Which annotation type(s) to download (default: both)",
    )
    parser.add_argument(
        "--output-dir",
        default="./genebass_output",
        help="Output directory for CSV files (default: ./genebass_output)",
    )
    parser.add_argument(
        "--format",
        choices=["csv", "tsv", "parquet"],
        default="csv",
        help="Output format (default: csv)",
    )
    parser.add_argument(
        "--mt-path",
        default="gs://ukbb-exome-public/500k/results/results.mt",
        help="GCS path to the GeneBass results MatrixTable",
    )
    parser.add_argument(
        "--describe-only",
        action="store_true",
        help="Just describe the schema and available annotations, don't export",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="If set, export in chunks of this many rows (useful for memory-limited machines)",
    )
    return parser.parse_args()


def init_hail():
    """Initialize Hail with GRCh38 reference genome."""
    try:
        import hail as hl
    except ImportError:
        print("ERROR: Hail is not installed. Install it with:")
        print("  pip install hail")
        print("\nHail requires Java 11. Install Java with:")
        print("  macOS:  brew install openjdk@11")
        print("  Ubuntu: sudo apt install openjdk-11-jdk")
        sys.exit(1)

    # Locate the GCS connector JAR (next to this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gcs_jar = os.path.join(script_dir, "gcs-connector-hadoop3-latest.jar")
    if not os.path.exists(gcs_jar):
        print(f"ERROR: GCS connector JAR not found at {gcs_jar}")
        print("Download it with:")
        print("  curl -L -o gcs-connector-hadoop3-latest.jar \\")
        print("    https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop3-latest.jar")
        sys.exit(1)

    print("Initializing Hail...")
    spark_conf = {
        "spark.jars": gcs_jar,
        "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
        "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",
        "spark.hadoop.google.cloud.auth.null.enable": "true",
        "spark.hadoop.fs.gs.auth.type": "UNAUTHENTICATED",
    }
    hl.init(default_reference="GRCh38", quiet=True, spark_conf=spark_conf)
    return hl


def describe_schema(hl, mt_path):
    """Load and describe the MatrixTable schema."""
    print(f"\nLoading MatrixTable from: {mt_path}")
    mt = hl.read_matrix_table(mt_path)

    print("\n" + "=" * 60)
    print("SCHEMA")
    print("=" * 60)
    mt.describe()

    print("\n" + "=" * 60)
    print("AVAILABLE ANNOTATIONS")
    print("=" * 60)
    annotations = mt.aggregate_rows(hl.agg.collect_as_set(mt.annotation))
    for ann in sorted(annotations):
        print(f"  - {ann}")

    print("\n" + "=" * 60)
    print("DIMENSIONS")
    print("=" * 60)
    print(f"  Rows (gene-annotation pairs): {mt.count_rows()}")
    print(f"  Columns (phenotypes):          {mt.count_cols()}")

    return mt


def export_annotation(hl, mt, annotation, output_dir, fmt, chunk_size=None):
    """Filter MatrixTable to a specific annotation and export."""
    print(f"\n{'=' * 60}")
    print(f"EXPORTING: {annotation}")
    print(f"{'=' * 60}")

    # Filter to the requested annotation
    filtered = mt.filter_rows(mt.annotation == annotation)
    n_rows = filtered.count_rows()
    n_cols = filtered.count_cols()
    print(f"  Rows: {n_rows}, Columns: {n_cols}")
    print(f"  Total entries to export: {n_rows * n_cols:,}")

    # Build output path
    ext = {"csv": "csv", "tsv": "tsv", "parquet": "parquet"}[fmt]
    safe_annotation = annotation.replace("|", "_")
    output_path = os.path.join(output_dir, f"genebass_{safe_annotation}_all.{ext}")

    start = time.time()

    # Use Hail's native export (streams data, avoids loading everything into memory)
    entries = filtered.entries()
    sep = "," if fmt in ("csv", "parquet") else "\t"
    # Use .bgz extension for block-gzip compression (much smaller output)
    compressed_path = output_path + ".bgz"
    hail_output = compressed_path + ".hail_tmp"
    print(f"  Exporting via Hail (streaming to disk, bgzip compressed)...")
    entries.export(hail_output, delimiter=sep)
    # Rename to final path
    os.rename(hail_output, compressed_path)
    output_path = compressed_path
    print(f"  Saved to: {output_path}")

    elapsed = time.time() - start
    file_size = os.path.getsize(output_path)
    print(f"  Done in {elapsed:.1f}s")
    print(f"  File size: {file_size / 1e9:.2f} GB")
    print(f"  Saved to: {output_path}")

    return output_path


def main():
    args = parse_args()
    hl = init_hail()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    if args.describe_only:
        describe_schema(hl, args.mt_path)
        return

    # Load the MatrixTable
    print(f"\nLoading MatrixTable from: {args.mt_path}")
    mt = hl.read_matrix_table(args.mt_path)

    # Show available annotations
    annotations = mt.aggregate_rows(hl.agg.collect_as_set(mt.annotation))
    print(f"Available annotations: {sorted(annotations)}")

    # Determine which annotations to export
    if args.annotation == "both":
        to_export = ["pLoF", "missense|LC"]
    elif args.annotation == "missense":
        to_export = ["missense|LC"]
    else:
        to_export = [args.annotation]

    # Validate requested annotations exist
    for ann in to_export:
        if ann not in annotations:
            # Try case-insensitive match
            match = [a for a in annotations if a.lower() == ann.lower()]
            if match:
                to_export[to_export.index(ann)] = match[0]
            else:
                print(f"WARNING: Annotation '{ann}' not found. Available: {sorted(annotations)}")
                to_export.remove(ann)

    if not to_export:
        print("ERROR: No valid annotations to export.")
        sys.exit(1)

    # Export each annotation
    output_files = []
    for ann in to_export:
        path = export_annotation(hl, mt, ann, args.output_dir, args.format, args.chunk_size)
        output_files.append(path)

    # Summary
    print(f"\n{'=' * 60}")
    print("COMPLETE")
    print(f"{'=' * 60}")
    for f in output_files:
        size = os.path.getsize(f) / 1e9
        print(f"  {f} ({size:.2f} GB)")


if __name__ == "__main__":
    main()
