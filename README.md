# E.V.A. (Exon & Variant Analyzer)

>A Multiomic Engine Mapping RNA-Seq Splicing & Variants to 3D Protein Biophysics and Targeted Therapeutics

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](#)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)](#)
[![Rendering](https://img.shields.io/badge/WebGL-Mol*-orange.svg)](#)

E.V.A. is an end-to-end bioinformatics and structural biology platform that connects transcriptomic dysregulation to 3D protein biophysics and targeted drug design. By integrating differential isoform quantification, bioinformatic de-confounding, fast 3D structure prediction, and spatial biophysical modeling, E.V.A. translates raw sequencing reads into precision therapeutics—including splice-junction-specific antisense oligonucleotides (ASOs) and small-molecule binding pocket coordinates.

Pipeline Architecture
RNA-Seq Ingestion & Quantification (Step 1): Automated SRA metadata querying, raw read filtering (fastp, MultiQC), and pseudoalignment quantification using salmon against reference transcriptomes.Differential Isoform Usage & Bioinformatic De-confounding (Step 2): Unbiased $\Delta\text{IF}$ candidate discovery coupled with cell-state marker profiling and single-sample GSEA (ssGSEA) to separate disease spliceopathy from myogenic differentiation artifacts.cDNA Translation & Functional Annotation (Step 3): Automated CDS retrieval via the Ensembl REST API, frame-shift detection, Nonsense-Mediated Decay (NMD) prediction, and Pfam domain loss/gain mapping.Rapid 3D Structural Modeling & Alignment (Step 4): High-throughput structure prediction via the ESMFold API, Kabsch structural superposition (BioPython.PDB.Superimposer), and spatial metrics computation ($\Delta\text{RMSD}$, $\Delta\text{SASA}$, surface charge redistribution).Therapeutic Target Design Engine (Step 5): Automated generation of 20-nt junction-specific ASO/siRNA target sequences and novel 3D structural pocket identification using fpocket.

Bioinformatic De-confounding Strategy
A core challenge in bulk RNA-seq splicing analysis of dystrophic muscle (DMD) is distinguishing primary pathological splicing dysregulation from secondary isoform shifts driven by active tissue regeneration and cell-type heterogeneity.Diagnostic Marker Panel: Tracks myogenic lineage progression (PAX7, MYF5, MYOD1, MYOG, CKM) and cell cycle activity (CDK1, CDKN1A, TOP2A, PCNA). Elevated expression of MYF5 and MYOG in DMD samples confirms an underlying biological shift toward active myoblast differentiation.Genome-Wide Pathway Scoring (ssGSEA): Single-sample GSEA across MSigDB Hallmark gene sets proves coordinated enrichment of the Myogenesis pathway in DMD replicates (NES ~0.42) while maintaining baseline stability across cell cycle checkpoints (p53, E2F Targets).Target Candidate Discovery: Filters genome-wide differential transcript usage (DTU) targets ($\Delta\text{IF} \approx 1.0$) to isolate high-confidence switching genes—such as MAPKAP1, CEPT1, MCUR1, and MYOT—for downstream 3D structural modeling.

```text
  ____   ____  ____ 
 |  __| \ \ / / / \   |  [E.V.A. v1.0 — Exon & Variant Analyzer]
 |  __|  \ V / / _ \  |  "Aw, what the hell, I don’t got that long a lifespan anyway." - R. Racoon
 |____|   \_/ /     \ |  Initializing ESMFold & RNA-Seq Multiomic Engine...