# src/build_index.py

import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import (
    REFERENCE_DIR,
    SALMON_INDEX_DIR,
    TRANSCRIPTOME_FASTA,
    TRANSCRIPTOME_URL,
)


def download_reference_fasta():
    """Downloads the Sus scrofa cDNA FASTA from Ensembl if not present."""
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    if TRANSCRIPTOME_FASTA.exists():
        print(f"[✓] Reference FASTA already exists at: {TRANSCRIPTOME_FASTA}")
        return

    print(f"[+] Downloading reference cDNA FASTA from Ensembl...")
    print(f"    URL: {TRANSCRIPTOME_URL}")

    try:
        # Stream download with basic progress logging
        def report_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = (downloaded / total_size) * 100 if total_size > 0 else 0
            sys.stdout.write(f"\r    Downloading: {percent:.1f}%")
            sys.stdout.flush()

        urllib.request.urlretrieve(
            TRANSCRIPTOME_URL, TRANSCRIPTOME_FASTA, reporthook=report_progress
        )
        print("\n[✓] Download complete.")
    except Exception as e:
        print(f"\n[X] Failed to download reference FASTA: {e}")
        if TRANSCRIPTOME_FASTA.exists():
            TRANSCRIPTOME_FASTA.unlink()
        sys.exit(1)


def build_salmon_index():
    """Executes salmon index if the index directory is missing or empty."""
    SALMON_INDEX_DIR.mkdir(parents=True, exist_ok=True)

    # Simple check: Salmon index creates specific header/version files in the target folder
    if any(SALMON_INDEX_DIR.iterdir()):
        print(f"[✓] Salmon index already exists at: {SALMON_INDEX_DIR}")
        return

    print("\n[+] Building Salmon index...")
    cmd = [
        "salmon",
        "index",
        "-t",
        str(TRANSCRIPTOME_FASTA),
        "-i",
        str(SALMON_INDEX_DIR),
        "-p",
        "4",
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"[✓] Salmon index built successfully at: {SALMON_INDEX_DIR}")
    except subprocess.CalledProcessError as e:
        print(f"[X] Salmon indexing failed: {e}")
        sys.exit(1)


def main():
    if shutil.which("salmon") is None:
        print("ERROR: 'salmon' binary not found in active conda environment.")
        sys.exit(1)

    download_reference_fasta()
    build_salmon_index()


if __name__ == "__main__":
    main()