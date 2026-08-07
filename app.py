# Example: app.py or pipeline.py
from src import ingestion, sequence, modeling, biophysics, therapeutics

def run_eva_pipeline(uploaded_file):
    # Step 1: Ingest & Filter
    filtered_rna = ingestion.parse_rna_seq(uploaded_file)

    # Step 2: Retrieve CDS & Translate
    cds_seq = sequence.fetch_ensembl_cds(filtered_rna.transcript_id)

    # Step 3: 3D Structure Prediction & Superposition
    pdb_canon, pdb_variant = modeling.run_esmfold(cds_seq)
    aligned_pdbs = modeling.align_kabsch(pdb_canon, pdb_variant)

    # Step 4: Spatial Delta Analysis
    rmsd, sasa, charge = biophysics.analyze_deltas(aligned_pdbs)

    # Step 5: Therapeutic Targets
    aso_target = therapeutics.design_aso(filtered_rna.junction)
    pockets = therapeutics.find_pockets(pdb_variant)

    return results