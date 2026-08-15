# src/run_qc.py

import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import RAW_FASTQ_DIR, TARGET_RUNS

QC_CLEANED_DIR = RAW_FASTQ_DIR.parent / "fastq_cleaned"
REPORTS_DIR = RAW_FASTQ_DIR.parent / "qc_reports"


def purge_file_from_ram(file_path: Path):
    """Tells the Linux OS kernel to immediately drop a file from RAM page cache."""
    if file_path.exists():
        try:
            fd = os.open(file_path, os.O_RDONLY)
            # POSIX_FADV_DONTNEED (value 4 on Linux) releases the cached memory pages
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            os.close(fd)
        except Exception:
            pass  # Non-Linux environments will safely ignore this call


def run_qc_on_sample(sra_id: str):
    # Keep original uncompressed files
    r1_in = RAW_FASTQ_DIR / f"{sra_id}_3.fastq"
    r1_out = QC_CLEANED_DIR / f"{sra_id}.clean.fastq"

    html_report = REPORTS_DIR / f"{sra_id}_fastp.html"
    json_report = REPORTS_DIR / f"{sra_id}_fastp.json"

    if not r1_in.exists():
        print(f"[X] Missing input file: {r1_in}")
        return

    print(f"\n[+] Running fastp for {sra_id}...")

    cmd = [
        "fastp",
        "--in1",
        str(r1_in),
        "--out1",
        str(r1_out),
        "--html",
        str(html_report),
        "--json",
        str(json_report),
        "--qualified_quality_phred",
        "20",
        "--thread",
        "6",
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"[✓] QC completed for {sra_id}.")
    except subprocess.CalledProcessError as e:
        print(f"[X] QC failed for {sra_id}:\n{e.stderr}")
    finally:
        # Purge both the input file and output file from Linux RAM cache
        purge_file_from_ram(r1_in)
        purge_file_from_ram(r1_out)


def generate_multiqc_summary():
    if shutil.which("multiqc") is not None:
        print("\n[+] Aggregating QC reports with MultiQC...")
        subprocess.run(
            [
                "multiqc",
                str(REPORTS_DIR),
                "-o",
                str(REPORTS_DIR / "multiqc_report"),
                "--force",
            ],
            check=True,
        )
        print(
            f"[✓] MultiQC report generated at: {REPORTS_DIR / 'multiqc_report' / 'multiqc_report.html'}"
        )


def main():
    if shutil.which("fastp") is None:
        print("ERROR: 'fastp' binary not found in active conda environment.")
        sys.exit(1)

    QC_CLEANED_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    for run_id in TARGET_RUNS:
        run_qc_on_sample(sra_id=run_id)

    generate_multiqc_summary()


if __name__ == "__main__":
    main()