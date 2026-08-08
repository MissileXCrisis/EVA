# config.py

from pathlib import Path

# Base Directory
BASE_DIR = Path(__file__).resolve().parent

# Data Paths
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
RAW_FASTQ_DIR = RAW_DIR / "fastq"
PROCESSED_DIR = DATA_DIR / "processed"

# SRA Target Runs (Pig DMD Satellite Cell Study - Differentiation Stage)
WT_RUNS = [
    "SRR32086846",  # 35WT Differentiation Rep 4
    "SRR32086848",  # 34WT Differentiation Rep 3
    "SRR32086850",  # 33WT Differentiation Rep 2
]

DMD_RUNS = [
    "SRR32086830",  # 17DMD Differentiation Rep 4
    "SRR32086832",  # 16DMD Differentiation Rep 3
    "SRR32086836",  # 14DMD Differentiation Rep 1
]

TARGET_RUNS = WT_RUNS + DMD_RUNS

# Step 1: E.V.A. Ingestion Thresholds
TPM_FOLD_CHANGE_THRESHOLD = 5.0
DELTA_PSI_THRESHOLD = 0.3

# API Endpoints for Steps 2 & 3
ENSEMBL_REST_API = "https://rest.ensembl.org"
ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"