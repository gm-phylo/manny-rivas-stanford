# GeneBass pLoF & Missense Download

## What This Downloads

From [GeneBass](https://app.genebass.org/downloads), the **gene burden results** MatrixTable containing:
- **pLoF** (predicted Loss of Function) variant associations
- **missense** variant associations

across ~4,500 phenotypes from the UK Biobank 500k exome dataset.

## Quick Start

```bash
# 1. Run the setup script (installs Java, Hail, then downloads)
bash setup_genebass.sh
```

Or manually:

```bash
# 1. Install prerequisites
pip install hail  # Requires Java 11+

# 2. Run the download
python download_genebass.py
```

## Output

Files are saved to `./genebass_output/`:

| File | Description |
|------|-------------|
| `genebass_pLoF_all.csv` | All pLoF gene burden associations |
| `genebass_missense_all.csv` | All missense gene burden associations |

## Options

```bash
# Download only pLoF
python download_genebass.py --annotation pLoF

# Download only missense
python download_genebass.py --annotation missense

# Export as Parquet (smaller, faster)
python download_genebass.py --format parquet

# Export as TSV
python download_genebass.py --format tsv

# Just inspect the schema (no download)
python download_genebass.py --describe-only

# Custom output directory
python download_genebass.py --output-dir /path/to/output
```

## Requirements

- **Python** 3.8+
- **Java** 11, or 17 (for Hail)
- **RAM**: 16 GB+ recommended (the MatrixTable is large)
- **Disk**: ~10-50 GB for the exported CSVs
- **Internet**: Access to Google Cloud Storage (public bucket, no auth needed)

## Alternative: Use Google Cloud Dataproc

For very large exports or if your local machine runs out of memory:

```bash
# Create a Hail-enabled Dataproc cluster
gcloud dataproc clusters create genebass-cluster \
    --region=us-central1 \
    --master-machine-type=n1-standard-8 \
    --worker-machine-type=n1-standard-8 \
    --num-workers=2 \
    --image-version=2.0-debian10 \
    --initialization-actions=gs://hail-common/hailctl/dataproc/0.2/init_notebook.py

# Submit the job
gcloud dataproc jobs submit pyspark download_genebass.py \
    --cluster=genebass-cluster \
    --region=us-central1

# IMPORTANT: Delete cluster when done to avoid charges
gcloud dataproc clusters delete genebass-cluster --region=us-central1
```

## Alternative: Raw gsutil Download

If you just want to copy the raw MatrixTable files locally:

```bash
# Download the full MatrixTable directory (~large)
gsutil -m cp -r gs://ukbb-exome-public/500k/results/results.mt ./results.mt

# Then process locally with Hail
python -c "
import hail as hl
hl.init(default_reference='GRCh38')
mt = hl.read_matrix_table('./results.mt')
mt.describe()
"
```

## Data Source

- **Website**: https://app.genebass.org/downloads
- **GCS Path**: `gs://ukbb-exome-public/500k/results/results.mt`
- **Reference**: Karczewski et al., *Nature* (2022)
