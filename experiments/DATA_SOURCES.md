# Data Sources Documentation

This document records the provenance, release/version identifiers, access metadata, and licensing status for all primary datasets and external cross-references used in the MoleculAR benchmark harness.

---

## 1. Primary Evaluation Database (DB-A)

* **File Path:** `experiments/data/drugrepo_oct6.db` (and root working copy `drugrepo.db`)
* **SHA-256 Checksum:** `D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`
* **Table Schema & Row Counts:**
  * `drugs`: 2,267 rows (columns: `chembl_id`, `drug_name`, `smiles`, `atc_category`, `indication`, `fingerprint_blob`, `hubness_score`, `hubness_percentile`, `generic_base_name`, `trivial_category`)
  * `targets`: 7,068 rows (columns: `chembl_id`, `uniprot_id`, `target_name`)
  * `side_effects`: 119,227 rows (columns: `chembl_id`, `side_effect_name`)
  * `crosswalk`: 2,267 rows (columns: `chembl_id`, `drugbank_id`, `pubchem_cid`)
* **Frozen Candidate Pool (C1):** 2,192 drugs with valid canonical SMILES, Morgan fingerprints (radius 2, 2048-bit), and non-null ChEMBL identifiers.

---

## 2. External Datasets & Provenance

### A. ChEMBL
* **Release:** ChEMBL 37 (`chembl_37_sqlite/chembl_37.db` and `chembl_webresource_client` v0.10.9)
* **Role in Pipeline:** Source of small molecule drug identities, approved/max-phase indications, target mechanisms (`drug_mechanism`), and structural representations.
* **Licence as Stated by Source:** CC BY-SA 3.0 (Creative Commons Attribution-ShareAlike 3.0 Unported), EMBL-EBI.

### B. DrugCentral
* **Release / Version Date:** DrugCentral 2027 Update (Human & Veterinary Drugs & Uses) (Published dump date: 2026-09-25; reference dump: `Drugcentral_2026-09-25.pgdump` / `drugcentral.dump.09022026.sql.gz`)
* **Download Portal:** `https://drugcentral.org/download`
* **Access Endpoint:** Official live public endpoint `unmtid-dbs.net:5433` (database: `drugcentral`, user: `drugman`), accessed 2026-11-10 05:31:46 UTC.
* **Extraction Condition:** `omop_relationship.relationship_name = 'indication'` (approved indications only, excluding off-label uses and contraindications).
* **Raw Extraction File:** `experiments/data/drugcentral_raw_approved_indications.csv` (SHA-256: `9d6eabc60feda0b3354020c8dc1f049bb9ee98fe619e4d87e1226d79988903f9`).
* **Licence as Stated by Source:** CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International), Division of Translational Informatics, University of New Mexico.

### C. SIDER (Side Effect Resource)
* **Release / Version:** SIDER 4.1 (`http://sideeffects.embl.de/`)
* **Downloaded Files:** `meddra_all_se.tsv.gz`, `drug_names.tsv`.
* **Role in Pipeline:** Extraction of MedDRA-coded phenotypic side effect terms mapped to drugs via STITCH/generic name matching.
* **Licence as Stated by Source:** CC BY-NC-SA 4.0 (Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International), EMBL.

### D. UniProt Mapping
* **Source:** Extracted from ChEMBL 37 target sequence components (`target_dictionary` $\to$ `target_components` $\to$ `component_sequences.accession`).
* **Role in Pipeline:** Populates the `uniprot_id` field in the DB-A `targets` table for shared-target Jaccard similarity (M4).
* **Licence as Stated by Source:** CC BY 4.0 (Creative Commons Attribution 4.0 International), UniProt Consortium.

### E. DrugBank / PubChem Crosswalk (in DB-A)
* **DrugBank Cross-References:**
  * **Source:** EMBL-EBI UniChem mapping release table `src1src2.txt.gz` (ChEMBL source 1 to DrugBank source 2).
  * **Licence as Stated by Source:** not confirmed (UniChem cross-reference table published by EMBL-EBI; underlying DrugBank database terms require specific commercial or academic non-commercial licensing, not confirmed for this pipeline).
* **PubChem Cross-References:**
  * **Source:** EMBL-EBI UniChem mapping table `src1src22.txt.gz` (ChEMBL source 1 to PubChem Compound source 22) and PubChem PUG REST API.
  * **Licence as Stated by Source:** Public Domain / Open Access (US National Center for Biotechnology Information / National Library of Medicine).
