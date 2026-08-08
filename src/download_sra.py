# src/download_sra.py

import shutil
import subprocess
import sys
from pathlib import Path

# Import central settings from config.py
# (Ensure script can resolve parent path if run from project root)
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import RAW_FASTQ_DIR, TARGET_RUNS


def check_tool_availability():
    """Checks if fasterq-dump is installed in the active environment."""
    if shutil.which("fasterq-dump") is None:
        print("ERROR: 'fasterq-dump' (SRA Toolkit) was not found.")
        print("Install via: conda install -c bioconda sra-tools -y")
        sys.exit(1)


def download_run(sra_id: str, output_dir: Path, threads: int = 4):
    """Downloads and extracts FASTQ files for a given SRA run accession."""
    print(f"\n[+] Fetching {sra_id}...")
    cmd = [
        "fasterq-dump",
        "--split-files",
        "--threads",
        str(threads),
        "--outdir",
        str(output_dir),
        sra_id,
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"[✓] Extracted: {sra_id}")
    except subprocess.CalledProcessError as e:
        print(f"[X] Failed {sra_id}:\n{e.stderr}")


def main():
    check_tool_availability()
    RAW_FASTQ_DIR.mkdir(parents=True, exist_ok=True)

    print(
        f"Downloading {len(TARGET_RUNS)} SRA runs to '{RAW_FASTQ_DIR.resolve()}'..."
    )

    for run_id in TARGET_RUNS:
        expected_file = RAW_FASTQ_DIR / f"{run_id}_1.fastq"
        if expected_file.exists():
            print(f"[!] {run_id} already exists. Skipping.")
            continue

        download_run(sra_id=run_id, output_dir=RAW_FASTQ_DIR)


if __name__ == "__main__":
    main()