# MoleculAR Evaluation Protocol Specifications

> **File:** `experiments/PROTOCOL.md`  
> **Git Status:** Untracked (`experiments/` is untracked in Git; `PROTOCOL.md` has never been committed to Git).  
> **Status:** LOCKED before running evaluation. Never modify after seeing results.  
> **Random Seed:** 42  

---

## 1. Candidate Pools & Databases

1. **Primary Pool (C1, C2, C3):**
   - Derived from `DB-A` ([`experiments/data/drugrepo_oct6.db`](data/drugrepo_oct6.db), SHA-256 `D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`), which matches the restored root `drugrepo.db` and `backup/drugrepo.db.bak`.
   - Excludes the 65 trivial compounds (`trivial_category IS NOT NULL`).
   - Excludes the 10 inorganic / element / simple-ion compounds:
     `CHEMBL1200345`, `CHEMBL1200458`, `CHEMBL1200539`, `CHEMBL1200691`, `CHEMBL1200696`, `CHEMBL1237066`, `CHEMBL1356`, `CHEMBL2068884`, `CHEMBL2106119`, `CHEMBL3833401`.
   - Primary Pool Size: **$N = 2,192$** drugs.
   - For any query $q$, candidates are all pool drugs except $q$ ($M = N - 1 = 2,191$), unless restricted by protocol (P2b, P3).
   - Uses stored Morgan fingerprints (`fingerprint_blob`, 2048-bit, radius 2) and stored `hubness_percentile`; neither is recomputed.

2. **Database Configurations:**
   - **C1 (DB-A Primary):** Primary pool ($N=2,192$) with `DB-A` targets and side effects.
   - **C2 (DB-R Primary):** Primary pool ($N=2,192$) with `DB-R` ([`experiments/data/drugrepo_oct6_recovered.db`](data/drugrepo_oct6_recovered.db), SHA-256 `66928166E834ABE98AF895E7E8A5234D11EB2D106EF1A034C9801EA027F07A38`), having 314 recovered target rows across 235 drugs.
   - **C3 (DB-A Rule-a-only):** Primary pool ($N=2,192$) on `DB-A`, but gold indications restricted strictly to mapping Rule (a) (full InChIKey matches only; Rules (b) and (c) omitted).
   - **C4 (DB-A with Inorganics):** Full non-trivial pool including the 10 inorganics ($N=2,202$, $M=2,201$) on `DB-A`.

---

## 2. Gold Standard & Indication Exclusions

1. **Gold File:** [`experiments/data/gold_pairs_primary.csv`](data/gold_pairs_primary.csv) (therapeutic pairs only).
2. **Diagnostic Indications Excluded:** 12 diagnostic indications listed in [`experiments/data/excluded_indications.csv`](data/excluded_indications.csv).
3. **Rule (b) Exclusions:** 9 rejected Rule (b) candidate mappings in [`experiments/data/rejected_rule_b.csv`](data/rejected_rule_b.csv):
   - 4 original name mismatches: `hyoscyamine->ATROPINE`, `acetylmethadol->LEVOMETHADYL`, `levmetamfetamine->METHAMPHETAMINE`, `doxribtimine->TELBIVUDINE`.
   - 5 stereo-prefix name differences: `razoxane->DEXRAZOXANE`, `medetomidine->DEXMEDETOMIDINE`, `dexrabeprazole->RABEPRAZOLE`, `levomethadone->METHADONE`, `levomenthol->MENTHOL`.
4. **Gold Counts on DB-A Primary:**
   - Mapped drugs with therapeutic indications: **1,796**
   - Therapeutic pairs by rule: a = **9,565**, b = **291**, c = **58** (Total: **9,914**)

---

## 3. Signals & Scoring

For query $q$ and candidate $c$:
- **Structural Similarity ($S$):** Tanimoto similarity of stored 2048-bit Morgan fingerprints.
- **Target Similarity ($T$):** Jaccard similarity of UniProt ID target sets: $T(q, c) = \frac{|Targets(q) \cap Targets(c)|}{|Targets(q) \cup Targets(c)|}$. If either set is empty, $T = 0.0$.
- **Side-Effect Similarity ($E$):** Jaccard similarity of MedDRA side effect concept names: $E(q, c) = \frac{|SE(q) \cap SE(c)|}{|SE(q) \cup SE(c)|}$. If either set is empty, $E = 0.0$.

---

## 4. Ranking & Tie Handling

For a query $q$, each candidate $c$ in the candidate pool receives score $\text{score}(c)$.
Rank of candidate $c$ uses average rank for ties:
$$r(c) = 1 + \sum_{c' \neq c} \mathbb{I}(\text{score}(c') > \text{score}(c)) + 0.5 \sum_{c' \neq c} \mathbb{I}(\text{score}(c') == \text{score}(c))$$

---

## 5. Per-Query Metrics

For query $q$ with relevant drug set $Rel(q)$ ($|Rel(q)| = R > 0$) and candidate set size $M$:
1. **Minimum Rank of First Relevant Drug:**
   $$r_{\min} = \min_{d \in Rel(q)} r(d)$$
2. **Hit@k:**
   $$\text{hit@}k = \begin{cases} 1 & \text{if } r_{\min} \le k \\ 0 & \text{otherwise} \end{cases}$$
   Evaluated for $k \in \{10, 20\}$.
3. **Mean Reciprocal Rank (MRR):**
   $$\text{MRR} = \frac{1}{r_{\min}}$$
4. **Precision@10:**
   $$\text{precision@}10 = \frac{|\{d \in Rel(q) \mid r(d) \le 10\}|}{10}$$
5. **Percentile Rank of First Relevant Drug:**
   $$\text{percentile\_rank} = \frac{r_{\min}}{M}$$
   (Lower is better; rank 1 out of $M$ gives $1/M$).
6. **Enrichment:**
   $$\text{Enrichment} = \frac{\overline{\text{hit@}10}}{\overline{\text{rand\_hit@}10}}$$
   Where analytical random hit@10 for candidate size $M$ and relevant size $R$ is:
   $$\text{rand\_hit@}10(q) = 1 - \frac{\binom{M - R}{10}}{\binom{M}{10}} \quad (\text{for } M \ge 10, M - R \ge 10)$$

---

## 6. Protocols

1. **P1 (All Relevant):**
   - Relevant $Rel_{P1}(q)$: all other pool drugs sharing $\ge 1$ gold therapeutic indication with $q$.
   - Candidates: all pool drugs except $q$ ($M = 2,191$).
   - Drop queries where $|Rel_{P1}(q)| = 0$.
2. **P2a ($\theta \in \{0.3, 0.4, 0.5\}$ - Analogue Neutralization):**
   - Relevant $Rel_{P2a,\theta}(q) = Rel_{P1}(q) \setminus \{c \mid S(q, c) > \theta\}$.
   - Structural analogues ($S > \theta$) **remain** in the candidate list as non-relevant competitors.
   - Candidates: all pool drugs except $q$ ($M = 2,191$).
   - Drop queries where $|Rel_{P2a,\theta}(q)| = 0$.
3. **P2b ($\theta = 0.4$ - Analogue Removal):**
   - Relevant $Rel_{P2b}(q) = Rel_{P2a, 0.4}(q)$.
   - All analogues with $S(q, c) > 0.4$ are **removed** from the candidate list for $q$.
   - Candidate pool size shrinks per query: $M_q = 2,191 - |\{c \mid S(q, c) > 0.4, c \neq q\}|$.
   - Analytical random expectation uses reduced $M_q$.
   - Drop queries where $|Rel_{P2b}(q)| = 0$.
4. **P3 (Cross-Class Relevant, ATC Level-3 Divergence - Primary):**
   - Candidates restricted to pool drugs that possess $\ge 1$ DrugCentral ATC level-3 code ($M = N_{\text{ATC}} - 1 = 1,899$).
   - Query must possess $\ge 1$ ATC level-3 code.
   - Relevant $Rel_{P3}(q) = \{c \in \text{pool with ATC} \mid c \in Rel_{P1}(q) \text{ and } ATC_3(q) \cap ATC_3(c) = \emptyset\}$.
   - Drop queries where $|Rel_{P3}(q)| = 0$ (or query lacks ATC code).
5. **P3-all (Cross-Class Relevant, Full Candidate Pool):**
   - Same relevant set as P3.
   - Candidate pool is the full primary pool ($M = 2,191$).

---

## 6.1 Duplicate-Neutral Protocols (-dn)

1. **Same-Drug Duplicate Definition:**
   For any query $q$, candidate $c \neq q$ is a same-drug duplicate iff:
   - Same `generic_base_name` (case-insensitive); OR
   - Same first-block InChIKey; OR
   - Same desalted canonical SMILES (using RDKit `LargestFragmentChooser`).
   Let $D(q) = \{c \in \text{pool} \setminus \{q\} \mid c \text{ is a duplicate of } q\}$.

2. **Protocol Variants:**
   - **P1-dn:**
     - Candidates: $\text{pool} \setminus (\{q\} \cup D(q))$, $M_q = 2,191 - |D(q)|$.
     - Relevant: $Rel_{P1\text{-dn}}(q) = Rel_{P1}(q) \setminus D(q)$.
     - Drop queries where $|Rel_{P1\text{-dn}}(q)| = 0$.
     - Analytical random hit@10 uses per-query $M_q$ and $R_q = |Rel_{P1\text{-dn}}(q)|$.
   - **P2a-dn ($\theta = 0.4$):**
     - Candidates: $\text{pool} \setminus (\{q\} \cup D(q))$, $M_q = 2,191 - |D(q)|$.
     - Relevant: $Rel_{P2a\text{-dn}}(q) = Rel_{P1}(q) \setminus (\{c \mid S(q, c) > 0.4\} \cup D(q))$.
     - Analogues ($S > 0.4$) not in $D(q)$ remain as non-relevant competitors in candidates; duplicates $D(q)$ are removed from both candidates and relevant set.
     - Drop queries where $|Rel_{P2a\text{-dn}}(q)| = 0$.
     - Analytical random hit@10 uses per-query $M_q$ and $R_q$.
   - **P2b-dn ($\theta = 0.4$):**
     - Candidates: $\text{pool} \setminus (\{q\} \cup \{c \mid S(q, c) > 0.4\} \cup D(q))$.
     - Relevant: $Rel_{P2b\text{-dn}}(q) = Rel_{P1}(q) \setminus (\{c \mid S(q, c) > 0.4\} \cup D(q))$.
     - Both analogues ($S > 0.4$) and duplicates $D(q)$ are removed from candidate list and relevant set.
     - Drop queries where $|Rel_{P2b\text{-dn}}(q)| = 0$.
     - Analytical random hit@10 uses per-query $M_q$ and $R_q$.
   - **P3-dn:**
     - Candidates: $(\text{pool with ATC}) \setminus (\{q\} \cup D(q))$.
     - Relevant: $Rel_{P3\text{-dn}}(q) = Rel_{P3}(q) \setminus D(q)$.
     - Candidate pool size: $M_q = (N_{\text{ATC}} - 1) - |D(q) \cap \text{pool with ATC}|$.
     - Drop queries where $|Rel_{P3\text{-dn}}(q)| = 0$.
     - Analytical random hit@10 uses per-query $M_q$ and $R_q$.

---

## 7. Methods Evaluated

- **M1 (Random Baseline):**
  - Analytical hit@10 used for expectation.
  - 100 seeded permutations (seed 42) per query used for empirical MRR and percentile rank.
- **M2 ($S$):** Fingerprint Tanimoto similarity.
- **M3 ($T$):** Target UniProt Jaccard similarity.
- **M4 ($E$):** Side-effect concept Jaccard similarity.
- **M5 ($S + T$):** Composite score $\frac{0.4 S + 0.4 T}{0.8} = 0.5 S + 0.5 T$.
- **M6 ($S + T + E$ Unpenalized):** Composite score $0.4 S + 0.4 T + 0.2 E$.
- **M7 (Production Score):**
  $$\text{score}(q, c) = \frac{0.4 S(q, c) + 0.4 T(q, c) + 0.2 E(q, c)}{1 + (\text{hubness\_percentile}(c) / 100)^2}$$
- **M8 (Indication Popularity Baseline):**
  - Rank candidates by their total number of therapeutic gold indications descending. Independent of query $q$.
- **M9 (Available-Signals Reweighted):**
  - A signal is "available" for pair $(q, c)$ iff both drugs have non-empty sets for that modality ($S$ is always available; $T$ available if both have $\ge 1$ target; $E$ available if both have $\ge 1$ side effect).
  - Score:
    $$\text{score}_{\text{avail}}(q, c) = \frac{\sum_{i \in \text{avail}} w_i s_i(q, c)}{\sum_{i \in \text{avail}} w_i}$$
    with $(w_S, w_T, w_E) = (0.4, 0.4, 0.2)$, then penalized by the same hubness penalty:
    $$\text{score}(q, c) = \frac{\text{score}_{\text{avail}}(q, c)}{1 + (\text{hubness\_percentile}(c) / 100)^2}$$

---

## 8. Statistical Inference & Bootstrap

- Mean metric across queries reported with 95% confidence intervals from **1,000 bootstrap resamples** of queries (percentile method: 2.5th and 97.5th percentiles; seed 42).
- For **Enrichment**, the ratio is computed inside each bootstrap resample:
  $$\text{Enrichment}^* = \frac{\frac{1}{|Q^*|} \sum_{q \in Q^*} \text{hit@}10(q)}{\frac{1}{|Q^*|} \sum_{q \in Q^*} \text{rand\_hit@}10(q)}$$
- Stratified reporting for C1 across relevant set size strata:
  - Small: $1 \le R \le 10$
  - Medium: $11 \le R \le 50$
  - Large: $R > 50$

---

## 8.1 Paired Comparisons

On C1 for `P1-dn`, `P2a-dn (0.4)`, `P2b-dn (0.4)`, and `P3-dn`:
- Compute paired bootstrap differences using the same 1,000 resamples of queries (seed 42).
- Differences evaluated for metric $\in \{\text{MRR}, \text{hit@}10\}$:
  - M6 - M5
  - M7 - M6
  - M7 - M5
  - M9 - M7
  - M7 - M2
  - M7 - M8
  - M6 - M2
- 95% confidence intervals computed via empirical 2.5th and 97.5th percentiles of resample mean differences.
- Saved to `experiments/results/paired_differences.csv`.

---

## 8.2 Weight Sensitivity (20 Random Splits)

- Grid search over 66 triples $(w_S, w_T, w_E)$ with step size 0.1 summing to 1.0.
- Scoring form: M7 score $(w_S S + w_T T + w_E E) / (1 + (p/100)^2)$.
- Objectives: MRR on `P2b-dn (0.4)` and MRR on `P1-dn`.
- For each objective, run 20 random 50/50 splits of query drugs using seeds 0 to 19.
- For each split:
  - Best triple on tune set by MRR.
  - Test MRR and test hit@10 for the best tune triple.
  - Test MRR and test hit@10 for default triple (0.4, 0.4, 0.2).
  - Default triple's rank among all 66 triples on tune set.
- Report selection frequency across 20 splits, and mean/standard deviation of test difference (best minus default).
- Saved to `experiments/results/weight_stability_p2b_dn.csv` and `experiments/results/weight_stability_p1_dn.csv`.

---

## 8.3 Hub Penalty Exponent Ablation

- Score formulation:
  $$\text{score}(q, c) = \frac{0.4 S(q, c) + 0.4 T(q, c) + 0.2 E(q, c)}{1 + (p(c) / 100)^k}$$
  for $k \in \{0, 1, 2, 3, 4\}$, where $p(c)$ is the stored `hubness_percentile`.
- Protocols evaluated: `P1-dn`, `P2a-dn (0.4)`, `P2b-dn (0.4)`, `P3-dn`.
- Performance metrics: MRR, hit@10, and enrichment with 95% bootstrap CIs (1,000 resamples, seed 42).
- Recurrence metrics across all queries:
  - Mean hubness percentile of top-10 candidates.
  - Maximum share of top-10 appearances by any single candidate.
  - Top-1% candidates' share of total top-10 appearances.
- Saved to `experiments/results/hub_penalty_ablation.csv`.

---

## 8.4 Failure Analysis

- For M7 on `P2b-dn (0.4)` and `P3-dn`:
  - List 20 queries with best first-relevant rank ($r_{\min}$) and 20 queries with worst first-relevant rank.
  - Columns: query ChEMBL ID, query name, number of relevant candidates ($|Rel(q)|$), first relevant candidate ChEMBL ID, first relevant candidate name, first relevant candidate rank, similarity signals ($S, T, E$), whether query and candidate have targets, whether query and candidate have side effects.
- Saved to `experiments/results/failure_analysis_p2b_dn.csv` and `experiments/results/failure_analysis_p3_dn.csv`.
