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
* **Provenance Note (2026-11-26):** Verified against the local SQLite database `chembl_37/chembl_37_sqlite/chembl_37.db` (`30,480,314,368` bytes; local SHA-256: `4BE13DF3B68E25DCD0BFF44BF094033B5AEBE98F415ACDC8C1CDF380E0C15142`). Its `version` table records `ChEMBL_37` (`creation_date` `2026-05-01 00:00:00.000000`), `INSTALL_sqlite` records `DB Version chembl_37` (`01/05/2026`, SQLite `3.26.0`), and all six core table row counts (`2,921,148` compounds, `3,824,604` compound records, `24,527,044` activities, `1,970,438` assays, `18,552` targets, `101,100` documents) match `chembl_37_release_notes.txt`. The EMBL-EBI download URL for the SQLite release archive is `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz` (`5,764,252,857` bytes via HTTP `HEAD`), whose published SHA-256 in `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/checksums.txt` is `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`. EMBL-EBI does not publish a separate checksum for the extracted `chembl_37.db` file, and the `.tar.gz` archive is not retained locally.

### B. DrugCentral
* **Release / Version Date:** DrugCentral 2027 Update (Human & Veterinary Drugs & Uses) (Published dump date: 2026-09-25; reference dump: `Drugcentral_2026-09-25.pgdump` / `drugcentral.dump.09022026.sql.gz`)
* **Download Portal:** `https://drugcentral.org/download`
* **Access Endpoint:** Official live public endpoint `unmtid-dbs.net:5433` (database: `drugcentral`, user: `drugman`), accessed 2026-11-10 05:31:46 UTC.
* **Extraction Condition:** `omop_relationship.relationship_name = 'indication'` (approved indications only, excluding off-label uses and contraindications).
* **Raw Extraction File:** `experiments/data/drugcentral_raw_approved_indications.csv` (SHA-256: `9d6eabc60feda0b3354020c8dc1f049bb9ee98fe619e4d87e1226d79988903f9`).
* **Licence as Stated by Source:** CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International), Division of Translational Informatics, University of New Mexico.
* **Provenance Note (2026-11-25):** The extraction was run on 2026-11-10 by SQL against the public PostgreSQL instance `unmtid-dbs.net:5433`, not restored from the `.pgdump` file. Verified 2026-11-25 on that instance: `SELECT * FROM dbversion` returns version 56, `dtime` `2026-09-02 00:55:33`. The server is PostgreSQL 16.15. The DrugCentral portal labels this build "DrugCentral 2027 / 2027 Update". No paper describes this release yet, so the paper cites DrugCentral 2023 (Avram et al., Nucleic Acids Res. 2023) as the resource reference. The ATC codes used for P3, P4 and P4b come from the same DrugCentral extract. No separate WHO ATC release version is recorded.

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
