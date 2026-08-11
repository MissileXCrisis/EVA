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
    DMD_RUNS,
    PROCESSED_DIR,
    QUANT_DIR,
    SALMON_INDEX_DIR,
    TARGET_RUNS,
    WT_RUNS,
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
    """Executes Salmon quantification on a single cleaned sample."""
    fastq_in = CLEANED_FASTQ_DIR / f"{sra_id}.clean.fastq"
    output_sample_dir = QUANT_DIR / sra_id
    quant_file = output_sample_dir / "quant.sf"

    if not fastq_in.exists():
        print(f"[X] Missing cleaned FASTQ: {fastq_in}")
        return None

    if quant_file.exists():
        print(f"[✓] Quant file already exists for {sra_id}, skipping.")
        return quant_file

    print(f"\n[+] Running Salmon quantification for {sra_id}...")

    cmd = [
        "salmon",
        "quant",
        "-i",
        str(SALMON_INDEX_DIR),
        "-l",
        "A",
        "-r",
        str(fastq_in),
        "-o",
        str(output_sample_dir),
        "--validateMappings",
        "-p",
        "6",
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

    for sra_id in TARGET_RUNS:
        quant_file = QUANT_DIR / sra_id / "quant.sf"
        if not quant_file.exists():
            print(f"[X] Quant file missing for {sra_id}. Cannot compute matrix.")
            sys.exit(1)

        df = pd.read_csv(quant_file, sep="\t", usecols=["Name", "TPM"])
        df["Transcript_ID"] = df["Name"].str.split(".").str[0]
        dataframes[sra_id] = df.set_index("Transcript_ID")["TPM"]

    master_df = pd.DataFrame(dataframes)

    # Compute group averages
    master_df["mean_TPM_WT"] = master_df[WT_RUNS].mean(axis=1)
    master_df["mean_TPM_DMD"] = master_df[DMD_RUNS].mean(axis=1)

    # Compute Fold Change
    master_df["fold_change"] = (master_df["mean_TPM_DMD"] + 0.01) / (
        master_df["mean_TPM_WT"] + 0.01
    )

    # Compute PSI proxy (Isoform Fraction)
    total_tpm_dmd = master_df["mean_TPM_DMD"].sum()
    master_df["PSI_DMD"] = (
        master_df["mean_TPM_DMD"] / total_tpm_dmd if total_tpm_dmd > 0 else 0
    )

    # Filter out unexpressed transcripts
    filtered_df = master_df[
        (master_df["mean_TPM_WT"] > 0.1) | (master_df["mean_TPM_DMD"] > 0.1)
    ].reset_index()

    output_csv = PROCESSED_DIR / "dmd_splicing_matrix.csv"
    filtered_df.to_csv(output_csv, index=False)

    print(f"[✓] Matrix saved ({len(filtered_df)} transcripts) to: {output_csv}")
    return output_csv


def main():
    if shutil.which("salmon") is None:
        print("ERROR: 'salmon' binary not found in active environment.")
        sys.exit(1)

    QUANT_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    for run_id in TARGET_RUNS:
        run_salmon_quant(sra_id=run_id)

    parse_and_aggregate_quants()


if __name__ == "__main__":
    main()