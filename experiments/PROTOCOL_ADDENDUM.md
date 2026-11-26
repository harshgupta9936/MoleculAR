# PROTOCOL ADDENDUM: Specifications for Protocols P4 and P4b

**Date:** 2026-11-14  
**Status:** APPROVED (Approved for evaluation in Task N10)  
**Parent Document:** [`experiments/PROTOCOL.md`](PROTOCOL.md) (Frozen - unmodified)

---

## 1. Motivation & Purpose

Protocol P3 enforces cross-class drug repurposing by requiring that candidate and query drugs share no ATC Level-3 (pharmacological subgroup) code. However, drugs within the same anatomical or organ system (ATC Level 1) may still share broad physiological or target class overlap.

**Protocol P4** and **Protocol P4b** establish a stricter anatomical/organ-level cross-class barrier:
A candidate drug is considered relevant to query drug $q$ if and only if they share a validated therapeutic indication in the gold standard but have **NO ATC Level-1 anatomical main group letter in common**.

- **Protocol P4 (Full Anatomical Candidate Pool):** Candidates comprise all pool drugs with $\ge 1$ ATC code ($M_q \approx 1,899$).
- **Protocol P4b (Anatomical Removal):** Candidates sharing any ATC Level-1 letter with query $q$ are removed from the candidate pool ($M_q \approx 1,636$). This mirrors the transition between Protocol P2a (analogue neutralization) and Protocol P2b (analogue removal).

This evaluates true physiological cross-domain drug repurposing (e.g., repurposing a dermatological agent for nervous system disease, or an antineoplastic agent for cardiovascular disease).

---

## 2. Formal Protocol Specifications

### 2.1 Eligibility and Pools
- **Drug Universe:** Primary pool $P$ consisting of $2,192$ non-trivial drugs from DB-A (excluding the 10 inorganic compounds defined in PROTOCOL.md Section 1.3).
- **ATC Eligibility:** Drugs must have at least one valid WHO ATC code mapped via DrugCentral (`drugcentral_raw_atc.csv`), excluding rejected Rule (b) structures.
  - Pool drugs with $\ge 1$ ATC code: **1,900** drugs ($P_{\text{ATC}}$).

### 2.2 Candidate Pools
- **Protocol P4 ($C_{q, 4}$):**
  - Restricted to pool drugs having at least one ATC code ($c \in P_{\text{ATC}}$).
  - Excludes the query drug itself ($c \ne q$).
  - **Duplicate-Neutral:** Excludes any same-drug duplicate of $q$ (matching generic base name, matching first-block InChIKey, or matching desalted canonical SMILES).
  - Effective candidate pool size: $M_{q, 4} = |P_{\text{ATC}}| - 1 - |\text{Duplicates}(q) \cap P_{\text{ATC}}|$.
- **Protocol P4b ($C_{q, 4b}$):**
  - Restricted to pool drugs having at least one ATC code ($c \in P_{\text{ATC}}$) that share **NO ATC Level-1 letter with query $q$**:
    $$C_{q, 4b} = \{ c \in P_{\text{ATC}} \setminus (\{q\} \cup \text{Duplicates}(q)) : \text{ATC}_1(q) \cap \text{ATC}_1(c) = \emptyset \}$$
  - Effective candidate pool size: $M_{q, 4b} < M_{q, 4}$.

### 2.3 Relevance Definition
A candidate drug $c$ is **relevant** to query $q$ if and only if:
1. $c$ and $q$ share at least one therapeutic indication in the primary gold standard (`gold_pairs_primary.csv`):
   $$\text{Ind}(q) \cap \text{Ind}(c) \ne \emptyset$$
2. $c$ and $q$ have **no ATC Level-1 letter in common**:
   $$\text{ATC}_1(q) \cap \text{ATC}_1(c) = \emptyset$$
   where $\text{ATC}_1(d) = \{ \text{code}[0] : \text{code} \in \text{ATC}(d) \}$.
3. $c$ is not a same-drug duplicate of $q$ ($c \notin \text{Duplicates}(q)$).

*(Note: Because all candidates in $C_{q, 4b}$ already have disjoint ATC Level-1 sets from $q$, the relevant set $\text{Rel}_{4b}(q)$ is identical to $\text{Rel}_4(q)$).*

### 2.4 Query Set
A pool drug $q \in P_{\text{ATC}}$ is an eligible **query** if and only if it has at least one relevant candidate under the condition above:
$$|\text{Rel}_4(q)| \ge 1$$

---

## 3. Dataset & Protocol Parameters (Audit Counts)

| Protocol Metric | Protocol P4 | Protocol P4b |
| :--- | :---: | :---: |
| Number of Qualifying Queries ($N$) | **965** | **965** |
| Total Relevant Pairs | **10,464** | **10,464** |
| Mean Candidate Pool Size ($M$) | **1898.95** (1897 – 1899) | **1636.40** (733 – 1894) |
| Mean Relevant Count per Query ($R$) | **10.84** (1 – 177) | **10.84** (1 – 177) |
| Analytical Random Hit@10 Expectation | **0.052461** (5.25%) | **0.061161** (6.12%) |

### 3.1 Ties & Precision
- **Rank tie-handling:** Exact mid-rank average ranking:
  $$r(c) = 1 + N_{>} + 0.5 \times N_{=}$$
- **Precision:** IEEE-754 float64.
- **Bootstrap:** 1,000 resamples (seed 42) for 95% confidence intervals.
