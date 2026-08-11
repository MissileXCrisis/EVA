# src/quantify_and_aggregate.py

import os
import shutil
import subprocess
import sys
from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import (
    CLEANED_FASTQ_DIR,
    CONTROL_RUNS,
    DMD_RUNS,
    PROCESSED_DIR,
    QUANT_DIR,
    SALMON_INDEX_DIR,
    TARGET_RUNS,
)


def purge_file_from_ram(file_path: Path):
    """Tells the Linux kernel to drop a file from RAM cache to manage WSL memory."""
    if file_path.exists():
        try:
            fd = os.open(file_path, os.O_RDONLY)
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            os.close(fd)
        except Exception:
            pass


def run_salmon_quant(sra_id: str) -> Path:
    """Executes Salmon quantification on a single sample."""
    fastq_in = CLEANED_FASTQ_DIR / f"{sra_id}.clean.fastq"
    output_sample_dir = QUANT_DIR / sra_id
    quant_file = output_sample_dir / "quant.sf"

    if not fastq_in.exists():
        print(f"[X] Missing cleaned FASTQ: {fastq_in}")
        return None

    if quant_file.exists():
        print(f"[✓] Quant file already exists for {sra_id}, skipping Salmon execution.")
        return quant_file

    print(f"\n[+] Running Salmon quantification for {sra_id}...")

    cmd = [
        "salmon",
        "quant",
        "-i",
        str(SALMON_INDEX_DIR),
        "-l",
        "A",  # Automatic library type detection
        "-r",
        str(fastq_in),
        "-o",
        str(output_sample_dir),
        "--validateMappings",
        "-p",
        "4",
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"[✓] Salmon completed for {sra_id}.")
    except subprocess.CalledProcessError as e:
        print(f"[X] Salmon failed for {sra_id}:\n{e.stderr}")
        return None
    finally:
        purge_file_from_ram(fastq_in)

    return quant_file


def parse_and_aggregate_quants() -> Path:
    """Merges quant.sf files, calculates WT/DMD mean TPMs, fold changes, and PSI proxy."""
    print("\n[+] Merging sample quantifications into master splicing matrix...")

    dataframes = {}

    # Read each sample's quant.sf file
    for sra_id in TARGET_RUNS:
        quant_file = QUANT_DIR / sra_id / "quant.sf"
        if not quant_file.exists():
            print(f"[X] Quant file missing for {sra_id}. Cannot compute complete matrix.")
            sys.exit(1)

        # Read Name (Transcript ID) and TPM columns
        df = pd.read_csv(quant_file, sep="\t", usecols=["Name", "TPM"])
        # Clean Ensembl transcript IDs (remove version numbers if present, e.g., ENS...1 -> ENS...)
        df["Transcript_ID"] = df["Name"].str.split(".").str[0]
        dataframes[sra_id] = df.set_index("Transcript_ID")["TPM"]

    # Combine all sample TPMs into a single DataFrame
    master_df = pd.DataFrame(dataframes)

    # Calculate Group Means
    master_df["mean_TPM_WT"] = master_df[CONTROL_RUNS].mean(axis=1)
    master_df["mean_TPM_DMD"] = master_df[DMD_RUNS].mean(axis=1)

    # Calculate Fold Change and Fold-Change Ratio (handling division by zero)
    master_df["fold_change"] = (master_df["mean_TPM_DMD"] + 0.01) / (
        master_df["mean_TPM_WT"] + 0.01
    )

    # Calculate Relative Isoform Inclusion / PSI Proxy (Psi = TPM_isoform / Total_TPM_sample)
    total_tpm_dmd = master_df["mean_TPM_DMD"].sum()
    master_df["PSI_DMD"] = (
        master_df["mean_TPM_DMD"] / total_tpm_dmd if total_tpm_dmd > 0 else 0
    )

    # Filter out unexpressed transcripts to keep file lightweight
    filtered_df = master_df[
        (master_df["mean_TPM_WT"] > 0.1) | (master_df["mean_TPM_DMD"] > 0.1)
    ].reset_index()

    output_csv = PROCESSED_DIR / "dmd_splicing_matrix.csv"
    filtered_df.to_csv(output_csv, index=False)

    print(f"[✓] Matrix saved successfully ({len(filtered_df)} transcripts) to: {output_csv}")
    return output_csv


def main():
    if shutil.which("salmon") is None:
        print("ERROR: 'salmon' binary not found in active environment.")
        sys.exit(1)

    if not SALMON_INDEX_DIR.exists():
        print(f"ERROR: Salmon index directory not found at {SALMON_INDEX_DIR}.")
        print("Please build the index first using the salmon index command.")
        sys.exit(1)

    QUANT_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    for run_id in TARGET_RUNS:
        run_salmon_quant(sra_id=run_id)

    parse_and_aggregate_quants()


if __name__ == "__main__":
    main()