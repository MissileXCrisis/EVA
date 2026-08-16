# E.V.A. (Exon & Variant Analyzer)

> A Multiomic Engine Mapping RNA-Seq Splicing & Variants to 3D Protein Biophysics and Targeted Therapeutics

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](#)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)](#)
[![Rendering](https://img.shields.io/badge/WebGL-Mol*-orange.svg)](#)

E.V.A. is an end-to-end bioinformatics and structural biology pipeline that connects transcriptomic dysregulation to 3D protein biophysics and targeted drug design. By integrating differential isoform quantification, bioinformatic de-confounding, fast 3D structure prediction, and spatial biophysical modeling, E.V.A. translates raw sequencing reads into precision therapeutics—including splice-junction-specific antisense oligonucleotides (ASOs) and small-molecule binding pocket coordinates.

## Key Features

* **RNA-Seq Ingestion & Quantification:** Automated SRA retrieval, quality control (`fastp`, `MultiQC`), and pseudoalignment quantification using Salmon.
* **Bioinformatic De-confounding:** Regresses out cell-cycle artifacts and quantifies myogenic clock progression (MSigDB pathways) to isolate true disease-driven splicing events.
* **Translation & Domain Annotation:** Automated retrieval of CDS via the Ensembl REST API, predicting Nonsense-Mediated Decay (NMD) and primary domain gain/loss.
* **Rapid 3D Structural Modeling:** High-throughput structure prediction via the ESMFold API followed by Kabsch structural superposition (BioPython).
* **Spatial Delta & Biophysical Profiling:** Computes backbone dislocation ($\Delta\text{RMSD}$), Solvent Accessible Surface Area ($\Delta\text{SASA}$), and surface charge redistribution.
* **Therapeutic Design Engine:** Automatically generates 20-nt ASO/siRNA target sequences across novel exon-exon junctions and identifies newly formed structural pockets using `fpocket`.

```text
  ____   ____  ____ 
 |  __| \ \ / / / \   |  [E.V.A. v1.0 — Exon & Variant Analyzer]
 |  __|  \ V / / _ \  |  "Aw, what the hell, I don’t got that long a lifespan anyway." - R. Racoon
 |____|   \_/ /     \ |  Initializing ESMFold & RNA-Seq Multiomic Engine...