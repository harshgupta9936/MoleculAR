# MoleculAR Second-Round Findings & Analyses (Duplicate-Neutral)

> **Generated strictly from on-disk CSV files by:** `experiments/generate_findings_tables.py`  
> **Compliance:** Fulfills Requirements N1 and N2 (no hand-typed numbers; source CSV & SHA-256 printed under every table).

---

## 1. Paired Bootstrap Comparisons (M6 vs M5, M7 vs M6, M7 vs M5, M9 vs M7, M7 vs M2, M7 vs M8, M6 vs M2)

Paired bootstrap differences evaluated on Configuration C1 across duplicate-neutral protocols (`P1-dn`, `P2a-dn (0.4)`, `P2b-dn (0.4)`, `P3-dn`) with 1,000 paired resamples (seed 42).

### Table 1.P1-dn: Paired Differences on P1-dn

| Comparison | Metric | Mean Difference | 95% Bootstrap CI | Significant (CI excludes 0)? |
| :--- | :---: | :---: | :---: | :---: |
| M6 - M5 | MRR | +0.0211 | [+0.0148, +0.0276] | Yes |
| M6 - M5 | Hit@10 | +0.0175 | [+0.0108, +0.0247] | Yes |
| M7 - M6 | MRR | -0.0221 | [-0.0325, -0.0120] | Yes |
| M7 - M6 | Hit@10 | -0.0030 | [-0.0121, +0.0060] | No |
| M7 - M5 | MRR | -0.0010 | [-0.0117, +0.0103] | No |
| M7 - M5 | Hit@10 | +0.0145 | [+0.0036, +0.0259] | Yes |
| M9 - M7 | MRR | -0.0327 | [-0.0401, -0.0251] | Yes |
| M9 - M7 | Hit@10 | -0.0217 | [-0.0313, -0.0127] | Yes |
| M7 - M2 | MRR | +0.0987 | [+0.0807, +0.1157] | Yes |
| M7 - M2 | Hit@10 | +0.1085 | [+0.0898, +0.1278] | Yes |
| M7 - M8 | MRR | +0.4330 | [+0.4087, +0.4552] | Yes |
| M7 - M8 | Hit@10 | +0.4858 | [+0.4563, +0.5142] | Yes |
| M6 - M2 | MRR | +0.1208 | [+0.1046, +0.1368] | Yes |
| M6 - M2 | Hit@10 | +0.1115 | [+0.0952, +0.1290] | Yes |

> **Source File:** [`experiments\results\paired_differences.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/paired_differences.csv)  
> **SHA-256:** `ea8cae01af6419699f0b8b86268c4cc90d157b227bcbe4a7610635fbde83e93b`

### Table 1.P2a-dn (0.4): Paired Differences on P2a-dn (0.4)

| Comparison | Metric | Mean Difference | 95% Bootstrap CI | Significant (CI excludes 0)? |
| :--- | :---: | :---: | :---: | :---: |
| M6 - M5 | MRR | +0.0236 | [+0.0180, +0.0292] | Yes |
| M6 - M5 | Hit@10 | +0.0401 | [+0.0292, +0.0517] | Yes |
| M7 - M6 | MRR | +0.0320 | [+0.0206, +0.0430] | Yes |
| M7 - M6 | Hit@10 | +0.0109 | [-0.0018, +0.0249] | No |
| M7 - M5 | MRR | +0.0557 | [+0.0433, +0.0672] | Yes |
| M7 - M5 | Hit@10 | +0.0511 | [+0.0353, +0.0675] | Yes |
| M9 - M7 | MRR | -0.0297 | [-0.0369, -0.0232] | Yes |
| M9 - M7 | Hit@10 | -0.0304 | [-0.0419, -0.0182] | Yes |
| M7 - M2 | MRR | +0.2002 | [+0.1831, +0.2164] | Yes |
| M7 - M2 | Hit@10 | +0.2286 | [+0.2055, +0.2523] | Yes |
| M7 - M8 | MRR | +0.2206 | [+0.1974, +0.2452] | Yes |
| M7 - M8 | Hit@10 | +0.3830 | [+0.3538, +0.4146] | Yes |
| M6 - M2 | MRR | +0.1681 | [+0.1531, +0.1822] | Yes |
| M6 - M2 | Hit@10 | +0.2176 | [+0.1969, +0.2383] | Yes |

> **Source File:** [`experiments\results\paired_differences.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/paired_differences.csv)  
> **SHA-256:** `ea8cae01af6419699f0b8b86268c4cc90d157b227bcbe4a7610635fbde83e93b`

### Table 1.P2b-dn (0.4): Paired Differences on P2b-dn (0.4)

| Comparison | Metric | Mean Difference | 95% Bootstrap CI | Significant (CI excludes 0)? |
| :--- | :---: | :---: | :---: | :---: |
| M6 - M5 | MRR | +0.0280 | [+0.0212, +0.0355] | Yes |
| M6 - M5 | Hit@10 | +0.0243 | [+0.0146, +0.0340] | Yes |
| M7 - M6 | MRR | -0.0231 | [-0.0358, -0.0114] | Yes |
| M7 - M6 | Hit@10 | -0.0036 | [-0.0140, +0.0073] | No |
| M7 - M5 | MRR | +0.0049 | [-0.0085, +0.0177] | No |
| M7 - M5 | Hit@10 | +0.0207 | [+0.0079, +0.0340] | Yes |
| M9 - M7 | MRR | -0.0325 | [-0.0410, -0.0250] | Yes |
| M9 - M7 | Hit@10 | -0.0371 | [-0.0486, -0.0255] | Yes |
| M7 - M2 | MRR | +0.1626 | [+0.1437, +0.1816] | Yes |
| M7 - M2 | Hit@10 | +0.1447 | [+0.1246, +0.1666] | Yes |
| M7 - M8 | MRR | +0.3615 | [+0.3380, +0.3860] | Yes |
| M7 - M8 | Hit@10 | +0.4225 | [+0.3939, +0.4511] | Yes |
| M6 - M2 | MRR | +0.1857 | [+0.1668, +0.2044] | Yes |
| M6 - M2 | Hit@10 | +0.1483 | [+0.1295, +0.1672] | Yes |

> **Source File:** [`experiments\results\paired_differences.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/paired_differences.csv)  
> **SHA-256:** `ea8cae01af6419699f0b8b86268c4cc90d157b227bcbe4a7610635fbde83e93b`

### Table 1.P3-dn: Paired Differences on P3-dn

| Comparison | Metric | Mean Difference | 95% Bootstrap CI | Significant (CI excludes 0)? |
| :--- | :---: | :---: | :---: | :---: |
| M6 - M5 | MRR | +0.0068 | [+0.0031, +0.0103] | Yes |
| M6 - M5 | Hit@10 | +0.0210 | [+0.0097, +0.0331] | Yes |
| M7 - M6 | MRR | +0.0122 | [+0.0019, +0.0226] | Yes |
| M7 - M6 | Hit@10 | -0.0065 | [-0.0226, +0.0105] | No |
| M7 - M5 | MRR | +0.0189 | [+0.0087, +0.0291] | Yes |
| M7 - M5 | Hit@10 | +0.0145 | [-0.0024, +0.0323] | No |
| M9 - M7 | MRR | -0.0109 | [-0.0167, -0.0045] | Yes |
| M9 - M7 | Hit@10 | -0.0250 | [-0.0387, -0.0105] | Yes |
| M7 - M2 | MRR | +0.0351 | [+0.0229, +0.0472] | Yes |
| M7 - M2 | Hit@10 | +0.0540 | [+0.0306, +0.0774] | Yes |
| M7 - M8 | MRR | -0.0325 | [-0.0531, -0.0135] | Yes |
| M7 - M8 | Hit@10 | -0.0460 | [-0.0791, -0.0129] | Yes |
| M6 - M2 | MRR | +0.0229 | [+0.0136, +0.0328] | Yes |
| M6 - M2 | Hit@10 | +0.0605 | [+0.0403, +0.0823] | Yes |

> **Source File:** [`experiments\results\paired_differences.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/paired_differences.csv)  
> **SHA-256:** `ea8cae01af6419699f0b8b86268c4cc90d157b227bcbe4a7610635fbde83e93b`

---

## 2. Weight Stability Analysis (20 Random 50/50 Splits, Seeds 0–19)

Evaluates whether optimal modality weights $(w_S, w_T, w_E)$ selected on a 50% tuning split over the 66-simplex grid generalize to an independent 50% test split, compared against production default weights $(0.4, 0.4, 0.2)$.

### Table 2.P2b-dn (0.4): Weight Stability on P2b-dn (0.4)

| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Tune MRR | Best Test MRR | Default Test MRR | Test MRR Diff (Best - Def) | Default Tune MRR | Default Tune Rank (/66) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0 | (0.3, 0.2, 0.5) | 0.5479 | 0.5391 | 0.5226 | +0.0164 | 0.5106 | 24 |
| 1 | (0.3, 0.2, 0.5) | 0.5318 | 0.5551 | 0.5373 | +0.0178 | 0.4959 | 26 |
| 2 | (0.3, 0.2, 0.5) | 0.5264 | 0.5605 | 0.5301 | +0.0304 | 0.5031 | 23 |
| 3 | (0.2, 0.2, 0.6) | 0.5615 | 0.5200 | 0.4957 | +0.0243 | 0.5376 | 23 |
| 4 | (0.2, 0.3, 0.5) | 0.5225 | 0.5560 | 0.5404 | +0.0156 | 0.4929 | 26 |
| 5 | (0.3, 0.2, 0.5) | 0.5591 | 0.5279 | 0.5010 | +0.0269 | 0.5323 | 23 |
| 6 | (0.4, 0.2, 0.4) | 0.5439 | 0.5346 | 0.5150 | +0.0195 | 0.5182 | 20 |
| 7 | (0.3, 0.2, 0.5) | 0.5366 | 0.5503 | 0.5202 | +0.0301 | 0.5131 | 25 |
| 8 | (0.3, 0.2, 0.5) | 0.5386 | 0.5484 | 0.5266 | +0.0217 | 0.5066 | 24 |
| 9 | (0.2, 0.3, 0.5) | 0.5396 | 0.5389 | 0.5185 | +0.0204 | 0.5147 | 23 |
| 10 | (0.2, 0.2, 0.6) | 0.5146 | 0.5668 | 0.5415 | +0.0252 | 0.4917 | 21 |
| 11 | (0.1, 0.3, 0.6) | 0.5450 | 0.5278 | 0.5175 | +0.0102 | 0.5157 | 24 |
| 12 | (0.3, 0.2, 0.5) | 0.5436 | 0.5434 | 0.5222 | +0.0211 | 0.5110 | 24 |
| 13 | (0.3, 0.2, 0.5) | 0.5461 | 0.5409 | 0.5189 | +0.0220 | 0.5144 | 24 |
| 14 | (0.3, 0.2, 0.5) | 0.5385 | 0.5484 | 0.5204 | +0.0280 | 0.5128 | 18 |
| 15 | (0.2, 0.2, 0.6) | 0.5553 | 0.5261 | 0.5103 | +0.0159 | 0.5230 | 22 |
| 16 | (0.1, 0.3, 0.6) | 0.5338 | 0.5389 | 0.5247 | +0.0142 | 0.5085 | 24 |
| 17 | (0.3, 0.2, 0.5) | 0.5468 | 0.5401 | 0.5215 | +0.0186 | 0.5117 | 26 |
| 18 | (0.3, 0.2, 0.5) | 0.5445 | 0.5425 | 0.5117 | +0.0307 | 0.5215 | 20 |
| 19 | (0.3, 0.2, 0.5) | 0.5485 | 0.5385 | 0.5127 | +0.0258 | 0.5206 | 23 |
| **Mean** | — | — | **0.5422** | **0.5204** | **+0.0218 (std: 0.0057)** | — | — |

> **Source File:** [`experiments\results\weight_stability_p2b_dn.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/weight_stability_p2b_dn.csv)  
> **SHA-256:** `952206ddd6e9edeafc9e8aaa89fb7d9dad94190a171407955d1690c7ce564209`

### Table 2.P1-dn: Weight Stability on P1-dn

| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Tune MRR | Best Test MRR | Default Test MRR | Test MRR Diff (Best - Def) | Default Tune MRR | Default Tune Rank (/66) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0 | (0.3, 0.2, 0.5) | 0.6076 | 0.5985 | 0.5802 | +0.0183 | 0.5938 | 17 |
| 1 | (0.3, 0.2, 0.5) | 0.5935 | 0.6126 | 0.6015 | +0.0110 | 0.5725 | 19 |
| 2 | (0.3, 0.2, 0.5) | 0.5993 | 0.6068 | 0.5932 | +0.0135 | 0.5808 | 17 |
| 3 | (0.2, 0.2, 0.6) | 0.6024 | 0.6046 | 0.5827 | +0.0219 | 0.5913 | 14 |
| 4 | (0.2, 0.2, 0.6) | 0.5872 | 0.6198 | 0.6002 | +0.0196 | 0.5738 | 17 |
| 5 | (0.4, 0.2, 0.4) | 0.5923 | 0.6104 | 0.5995 | +0.0109 | 0.5746 | 20 |
| 6 | (0.2, 0.2, 0.6) | 0.5992 | 0.6078 | 0.5896 | +0.0182 | 0.5844 | 16 |
| 7 | (0.2, 0.2, 0.6) | 0.5886 | 0.6184 | 0.6049 | +0.0134 | 0.5691 | 17 |
| 8 | (0.3, 0.2, 0.5) | 0.6008 | 0.6052 | 0.5930 | +0.0122 | 0.5811 | 18 |
| 9 | (0.2, 0.2, 0.6) | 0.6040 | 0.6030 | 0.5958 | +0.0072 | 0.5783 | 21 |
| 10 | (0.2, 0.2, 0.6) | 0.6004 | 0.6066 | 0.5866 | +0.0201 | 0.5875 | 15 |
| 11 | (0.2, 0.2, 0.6) | 0.6156 | 0.5914 | 0.5771 | +0.0142 | 0.5969 | 19 |
| 12 | (0.3, 0.1, 0.6) | 0.6025 | 0.5866 | 0.5883 | -0.0018 | 0.5857 | 15 |
| 13 | (0.2, 0.2, 0.6) | 0.6071 | 0.5999 | 0.5799 | +0.0200 | 0.5941 | 13 |
| 14 | (0.4, 0.1, 0.5) | 0.6168 | 0.5764 | 0.5818 | -0.0054 | 0.5922 | 19 |
| 15 | (0.5, 0.2, 0.3) | 0.6211 | 0.5781 | 0.5682 | +0.0099 | 0.6059 | 13 |
| 16 | (0.3, 0.2, 0.5) | 0.5967 | 0.6093 | 0.5926 | +0.0168 | 0.5815 | 16 |
| 17 | (0.3, 0.2, 0.5) | 0.5966 | 0.6094 | 0.5967 | +0.0127 | 0.5773 | 14 |
| 18 | (0.2, 0.2, 0.6) | 0.6101 | 0.5968 | 0.5847 | +0.0121 | 0.5893 | 17 |
| 19 | (0.2, 0.2, 0.6) | 0.6064 | 0.6005 | 0.5846 | +0.0159 | 0.5894 | 18 |
| **Mean** | — | — | **0.6021** | **0.5891** | **+0.0130 (std: 0.0067)** | — | — |

> **Source File:** [`experiments\results\weight_stability_p1_dn.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/weight_stability_p1_dn.csv)  
> **SHA-256:** `5a9f5ba5a7d79b104defccc8f3dc7b2a2da52acf622f51afe6ddaf1b484aff6b`

---

## 3. Hub Penalty Exponent Ablation and Hubness Recurrence Statistics

Evaluates the impact of hub penalty exponent $k \in \{0, 1, 2, 3, 4\}$ on ranking performance and hub drug dominance.

### Table 3.1: Hub Penalty Exponent Ablation across Duplicate-Neutral Protocols

| Protocol | $k$ | N | MRR (95% CI) | Hit@10 (95% CI) | Enrichment (95% CI) | Mean Top-10 Hubness | Max Candidate Appearances | Max Candidate Share | Top 1% Candidate Share |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| P1-dn | 0 | 1659 | 0.6091 [0.5875, 0.6281] | 0.7571 [0.7366, 0.7764] | 5.22 [4.98, 5.49] | 58.00 | 41 | 0.25% | 4.09% |
| P1-dn | 1 | 1659 | 0.5914 [0.5694, 0.6120] | 0.7595 [0.7378, 0.7788] | 5.24 [4.99, 5.51] | 41.24 | 61 | 0.37% | 5.63% |
| P1-dn | 2 | 1659 | 0.5870 [0.5653, 0.6065] | 0.7541 [0.7312, 0.7740] | 5.20 [4.95, 5.47] | 41.08 | 55 | 0.33% | 5.54% |
| P1-dn | 3 | 1659 | 0.5899 [0.5676, 0.6088] | 0.7547 [0.7330, 0.7740] | 5.21 [4.95, 5.46] | 42.66 | 48 | 0.29% | 5.24% |
| P1-dn | 4 | 1659 | 0.5903 [0.5689, 0.6086] | 0.7559 [0.7336, 0.7752] | 5.22 [4.96, 5.47] | 44.27 | 48 | 0.29% | 4.96% |
| P2a-dn (0.4) | 0 | 1645 | 0.3426 [0.3252, 0.3607] | 0.6426 [0.6176, 0.6657] | 4.57 [4.33, 4.84] | 58.11 | 40 | 0.24% | 4.10% |
| P2a-dn (0.4) | 1 | 1645 | 0.3740 [0.3557, 0.3911] | 0.6553 [0.6316, 0.6784] | 4.66 [4.42, 4.92] | 41.30 | 61 | 0.37% | 5.65% |
| P2a-dn (0.4) | 2 | 1645 | 0.3746 [0.3551, 0.3914] | 0.6535 [0.6292, 0.6766] | 4.64 [4.42, 4.92] | 41.14 | 55 | 0.33% | 5.57% |
| P2a-dn (0.4) | 3 | 1645 | 0.3701 [0.3514, 0.3871] | 0.6535 [0.6292, 0.6766] | 4.64 [4.41, 4.91] | 42.72 | 47 | 0.29% | 5.26% |
| P2a-dn (0.4) | 4 | 1645 | 0.3652 [0.3465, 0.3816] | 0.6517 [0.6273, 0.6748] | 4.63 [4.40, 4.90] | 44.34 | 47 | 0.29% | 4.98% |
| P2b-dn (0.4) | 0 | 1645 | 0.5397 [0.5183, 0.5607] | 0.6967 [0.6736, 0.7191] | 4.94 [4.71, 5.21] | 58.66 | 42 | 0.26% | 4.47% |
| P2b-dn (0.4) | 1 | 1645 | 0.5221 [0.5013, 0.5430] | 0.6979 [0.6760, 0.7210] | 4.95 [4.73, 5.21] | 38.97 | 65 | 0.40% | 5.99% |
| P2b-dn (0.4) | 2 | 1645 | 0.5166 [0.4957, 0.5364] | 0.6930 [0.6699, 0.7155] | 4.91 [4.69, 5.18] | 39.04 | 56 | 0.34% | 5.85% |
| P2b-dn (0.4) | 3 | 1645 | 0.5204 [0.4988, 0.5399] | 0.6948 [0.6729, 0.7173] | 4.93 [4.70, 5.20] | 40.95 | 55 | 0.33% | 5.68% |
| P2b-dn (0.4) | 4 | 1645 | 0.5231 [0.5018, 0.5430] | 0.6930 [0.6705, 0.7155] | 4.91 [4.68, 5.19] | 42.77 | 52 | 0.32% | 5.34% |
| P3-dn | 0 | 1240 | 0.1101 [0.0990, 0.1222] | 0.2540 [0.2290, 0.2774] | 2.09 [1.88, 2.29] | 58.01 | 38 | 0.31% | 5.00% |
| P3-dn | 1 | 1240 | 0.1221 [0.1079, 0.1351] | 0.2460 [0.2218, 0.2702] | 2.02 [1.82, 2.24] | 41.22 | 59 | 0.48% | 6.69% |
| P3-dn | 2 | 1240 | 0.1223 [0.1088, 0.1347] | 0.2476 [0.2234, 0.2726] | 2.03 [1.83, 2.25] | 41.03 | 52 | 0.42% | 6.65% |
| P3-dn | 3 | 1240 | 0.1193 [0.1064, 0.1322] | 0.2460 [0.2226, 0.2686] | 2.02 [1.82, 2.23] | 42.51 | 43 | 0.35% | 6.31% |
| P3-dn | 4 | 1240 | 0.1192 [0.1058, 0.1316] | 0.2492 [0.2250, 0.2726] | 2.05 [1.85, 2.25] | 44.16 | 41 | 0.33% | 5.97% |

> **Source File:** [`experiments\results\hub_penalty_ablation.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/hub_penalty_ablation.csv)  
> **SHA-256:** `3e50d709443e304476bd943d7549f29f8acc46a6bc28b0ee231b9789aeb9f8b1`

---

## 4. Failure Analysis (Top 20 Best and Worst Queries for Method M7)

### Table 4.P2b-dn (0.4).1: Top 20 Best Performing Queries for M7 on P2b-dn (0.4)

| Query Drug | N (Rel) | First Rel Drug | Rank ($r_{\min}$) | $S$ | $T$ | $E$ | Targets (Q / Cand)? | SEs (Q / Cand)? |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| CHEMBL1650 (CORTISONE) | 297 | CHEMBL131 (PREDNISOLONE) | 1.0 | 0.382 | 1.000 | 0.377 | Y / Y | Y / Y |
| CHEMBL389621 (HYDROCORTISONE) | 288 | CHEMBL1200386 (PREDNICARBATE) | 1.0 | 0.395 | 1.000 | 0.500 | Y / Y | Y / Y |
| CHEMBL384467 (DEXAMETHASONE) | 278 | CHEMBL635 (PREDNISONE) | 1.0 | 0.392 | 1.000 | 0.521 | Y / Y | Y / Y |
| CHEMBL131 (PREDNISOLONE) | 257 | CHEMBL1200600 (FLUOROMETHOLONE) | 1.0 | 0.378 | 1.000 | 0.000 | Y / Y | Y / N |
| CHEMBL1451 (TRIAMCINOLONE) | 221 | CHEMBL1512 (FLUNISOLIDE) | 1.0 | 0.321 | 1.000 | 0.256 | Y / Y | Y / Y |
| CHEMBL1504 (TRIAMCINOLONE ACETONIDE) | 218 | CHEMBL1201064 (FLUOROMETHOLONE) | 1.0 | 0.395 | 1.000 | 0.000 | Y / Y | N / Y |
| CHEMBL632 (BETAMETHASONE) | 217 | CHEMBL1201109 (DESONIDE) | 1.0 | 0.388 | 1.000 | 0.299 | Y / Y | Y / Y |
| CHEMBL650 (METHYLPREDNISOLONE) | 217 | CHEMBL384467 (DEXAMETHASONE) | 1.0 | 0.377 | 1.000 | 0.543 | Y / Y | Y / Y |
| CHEMBL635 (PREDNISONE) | 206 | CHEMBL384467 (DEXAMETHASONE) | 1.0 | 0.392 | 1.000 | 0.521 | Y / Y | Y / Y |
| CHEMBL1091 (HYDROCORTISONE) | 204 | CHEMBL1201173 (MEDRYSONE) | 1.0 | 0.387 | 1.000 | 0.000 | Y / Y | N / Y |
| CHEMBL1200449 (TRIAMCINOLONE DIACETATE) | 199 | CHEMBL989 (FLUOCINOLONE ACETONIDE) | 1.0 | 0.388 | 1.000 | 0.000 | Y / Y | N / N |
| CHEMBL1200538 (BETAMETHASONE) | 181 | CHEMBL1473 (FLUTICASONE PROPIONATE) | 1.0 | 0.400 | 1.000 | 0.000 | Y / Y | N / N |
| CHEMBL1200532 (METHDILAZINE) | 155 | CHEMBL3989885 (TRIMEPRAZINE) | 1.0 | 0.391 | 1.000 | 0.000 | Y / Y | N / N |
| CHEMBL3989885 (TRIMEPRAZINE) | 155 | CHEMBL1255739 (METHAPYRILENE) | 1.0 | 0.196 | 1.000 | 0.000 | Y / Y | N / N |
| CHEMBL428647 (PACLITAXEL) | 155 | CHEMBL1683544 (ERIBULIN) | 1.0 | 0.097 | 1.000 | 0.269 | Y / Y | Y / Y |
| CHEMBL9 (NORFLOXACIN) | 146 | CHEMBL1200779 (TROVAFLOXACIN) | 1.0 | 0.333 | 1.000 | 0.278 | Y / Y | Y / Y |
| CHEMBL4 (OFLOXACIN) | 145 | CHEMBL1202 (CIPROFLOXACIN) | 1.0 | 0.386 | 0.500 | 0.423 | Y / Y | Y / Y |
| CHEMBL1201109 (DESONIDE) | 142 | CHEMBL1766 (DESOXIMETASONE) | 1.0 | 0.341 | 1.000 | 0.313 | Y / Y | Y / Y |
| CHEMBL1200806 (NAPROXEN) | 141 | CHEMBL898 (DIFLUNISAL) | 1.0 | 0.157 | 1.000 | 0.202 | Y / Y | Y / Y |
| CHEMBL1201071 (TETRACYCLINE PHOSPHATE COMPLEX) | 140 | CHEMBL1200506 (ERYTHROMYCIN) | 1.0 | 0.121 | 1.000 | 0.000 | Y / Y | N / Y |

> **Source File:** [`experiments\results\failure_analysis_p2b_dn.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/failure_analysis_p2b_dn.csv)  
> **SHA-256:** `bd8197e28ecbd4f3221c05afd41e20740d873b38d357e9b7a19f55d4b2d6cc7f`

### Table 4.P2b-dn (0.4).2: Top 20 Worst Performing Queries for M7 on P2b-dn (0.4)

| Query Drug | N (Rel) | First Rel Drug | Rank ($r_{\min}$) | $S$ | $T$ | $E$ | Targets (Q / Cand)? | SEs (Q / Cand)? |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| CHEMBL1200627 (PROPIOLACTONE) | 1 | CHEMBL582 (ISOPROPYL ALCOHOL) | 2172.0 | 0.000 | 0.000 | 0.000 | N / N | N / N |
| CHEMBL4101807 (VOXELOTOR) | 1 | CHEMBL467 (HYDROXYUREA) | 2164.0 | 0.034 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL582 (ISOPROPYL ALCOHOL) | 1 | CHEMBL1200627 (PROPIOLACTONE) | 2128.0 | 0.000 | 0.000 | 0.000 | N / N | N / N |
| CHEMBL1200358 (SODIUM PHENYLACETATE) | 1 | CHEMBL296306 (LACTULOSE) | 2123.0 | 0.041 | 0.000 | 0.000 | N / N | N / Y |
| CHEMBL426559 (LAROPIPRANT) | 1 | CHEMBL247951 (PEMAFIBRATE) | 2011.0 | 0.104 | 0.000 | 0.000 | Y / Y | N / N |
| CHEMBL247951 (PEMAFIBRATE) | 1 | CHEMBL426559 (LAROPIPRANT) | 1976.0 | 0.104 | 0.000 | 0.000 | Y / Y | N / N |
| CHEMBL1251 (GANIRELIX) | 1 | CHEMBL2355051 (CLOMIPHENE) | 1944.0 | 0.078 | 0.000 | 0.000 | Y / N | N / Y |
| CHEMBL3545185 (SELINEXOR) | 1 | CHEMBL451887 (CARFILZOMIB) | 1811.0 | 0.077 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL272427 (TAURURSODIOL) | 3 | CHEMBL1469 (PHENYLBUTANOIC ACID) | 1806.0 | 0.069 | 0.000 | 0.000 | N / N | N / N |
| CHEMBL2105680 (RIFAMYCIN) | 1 | CHEMBL1707 (LOPERAMIDE) | 1697.0 | 0.078 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL106258 (AMINOREX) | 1 | CHEMBL4084119 (LIRAGLUTIDE) | 1621.0 | 0.051 | 0.000 | 0.000 | N / Y | N / Y |
| CHEMBL2106155 (ENCAINIDE) | 1 | CHEMBL1200813 (ACEBUTOLOL) | 1548.0 | 0.145 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL1096885 (VALRUBICIN) | 1 | CHEMBL671 (THIOTEPA) | 1542.0 | 0.010 | 0.000 | 0.123 | Y / N | Y / Y |
| CHEMBL1256391 (PIRFENIDONE) | 1 | CHEMBL3039504 (NINTEDANIB ESYLATE) | 1531.0 | 0.128 | 0.000 | 0.000 | N / Y | Y / N |
| CHEMBL744 (RILUZOLE) | 3 | CHEMBL290916 (EDARAVONE) | 1530.0 | 0.075 | 0.000 | 0.000 | Y / N | Y / N |
| CHEMBL1863513 (INGENOL MEBUTATE) | 6 | CHEMBL1200582 (AMINOLEVULINIC ACID) | 1517.0 | 0.056 | 0.000 | 0.000 | N / N | N / N |
| CHEMBL1200729 (CHLORAMPHENICOL) | 3 | CHEMBL4 (OFLOXACIN) | 1516.0 | 0.078 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL2104415 (MOXIDECTIN) | 1 | CHEMBL684 (DIETHYLCARBAMAZINE) | 1500.0 | 0.047 | 0.000 | 0.000 | N / N | N / N |
| CHEMBL2311030 (ERGOLOID MESYLATES) | 1 | CHEMBL71752 (VINPOCETINE) | 1469.0 | 0.128 | 0.000 | 0.000 | Y / N | N / N |
| CHEMBL544 (MEQUINOL) | 1 | CHEMBL38 (TRETINOIN) | 1446.0 | 0.057 | 0.000 | 0.000 | Y / Y | N / N |

> **Source File:** [`experiments\results\failure_analysis_p2b_dn.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/failure_analysis_p2b_dn.csv)  
> **SHA-256:** `bd8197e28ecbd4f3221c05afd41e20740d873b38d357e9b7a19f55d4b2d6cc7f`

### Table 4.P3-dn.1: Top 20 Best Performing Queries for M7 on P3-dn

| Query Drug | N (Rel) | First Rel Drug | Rank ($r_{\min}$) | $S$ | $T$ | $E$ | Targets (Q / Cand)? | SEs (Q / Cand)? |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| CHEMBL1200467 (HYDROXYZINE) | 124 | CHEMBL1172 (DESLORATADINE) | 1.0 | 0.129 | 1.000 | 0.113 | Y / Y | Y / Y |
| CHEMBL1201109 (DESONIDE) | 113 | CHEMBL1512 (FLUNISOLIDE) | 1.0 | 0.721 | 1.000 | 0.082 | Y / Y | Y / Y |
| CHEMBL1200806 (NAPROXEN) | 96 | CHEMBL898 (DIFLUNISAL) | 1.0 | 0.157 | 1.000 | 0.202 | Y / Y | Y / Y |
| CHEMBL1200820 (PIPERACILLIN) | 92 | CHEMBL127 (MEROPENEM) | 1.0 | 0.144 | 1.000 | 0.254 | Y / Y | Y / Y |
| CHEMBL1200855 (TICARCILLIN DISODIUM) | 90 | CHEMBL927 (CEFDINIR) | 1.0 | 0.189 | 1.000 | 0.000 | Y / Y | N / Y |
| CHEMBL1082 (AMOXICILLIN) | 88 | CHEMBL927 (CEFDINIR) | 1.0 | 0.193 | 1.000 | 0.306 | Y / Y | Y / Y |
| CHEMBL1200326 (NICARDIPINE) | 88 | CHEMBL1280 (VERAPAMIL) | 1.0 | 0.154 | 1.000 | 0.293 | Y / Y | Y / Y |
| CHEMBL1201080 (TRIMETHOPRIM) | 87 | CHEMBL119 (TRIMETREXATE) | 1.0 | 0.393 | 0.000 | 0.093 | Y / N | Y / Y |
| CHEMBL1055 (CHLORTHALIDONE) | 81 | CHEMBL1763 (HYDROFLUMETHIAZIDE) | 1.0 | 0.182 | 1.000 | 0.615 | Y / Y | Y / Y |
| CHEMBL193 (NIFEDIPINE) | 81 | CHEMBL1280 (VERAPAMIL) | 1.0 | 0.091 | 1.000 | 0.252 | Y / Y | Y / Y |
| CHEMBL1237132 (CLEVIDIPINE) | 78 | CHEMBL1280 (VERAPAMIL) | 1.0 | 0.105 | 1.000 | 0.069 | Y / Y | Y / Y |
| CHEMBL1648 (ISRADIPINE) | 78 | CHEMBL1280 (VERAPAMIL) | 1.0 | 0.106 | 1.000 | 0.283 | Y / Y | Y / Y |
| CHEMBL1726 (NISOLDIPINE) | 78 | CHEMBL1280 (VERAPAMIL) | 1.0 | 0.119 | 1.000 | 0.235 | Y / Y | Y / Y |
| CHEMBL1200965 (BACAMPICILLIN) | 77 | CHEMBL3544973 (CEFTAROLINE FOSAMIL) | 1.0 | 0.197 | 1.000 | 0.000 | Y / Y | N / N |
| CHEMBL1558 (PRAZOSIN) | 77 | CHEMBL3989562 (TERAZOSIN) | 1.0 | 0.585 | 1.000 | 0.343 | Y / Y | Y / Y |
| CHEMBL1697708 (MEZLOCILLIN) | 77 | CHEMBL127 (MEROPENEM) | 1.0 | 0.143 | 1.000 | 0.000 | Y / Y | N / Y |
| CHEMBL1639 (ALISKIREN) | 72 | CHEMBL1280 (VERAPAMIL) | 1.0 | 0.236 | 0.000 | 0.164 | N / Y | Y / Y |
| CHEMBL1200852 (PENICILLIN V) | 68 | CHEMBL3544973 (CEFTAROLINE FOSAMIL) | 1.0 | 0.189 | 1.000 | 0.000 | Y / Y | N / N |
| CHEMBL898 (DIFLUNISAL) | 68 | CHEMBL15770 (SULINDAC) | 1.0 | 0.246 | 1.000 | 0.497 | Y / Y | Y / Y |
| CHEMBL126 (LINEZOLID) | 64 | CHEMBL1741 (CLARITHROMYCIN) | 1.0 | 0.088 | 1.000 | 0.253 | Y / Y | Y / Y |

> **Source File:** [`experiments\results\failure_analysis_p3_dn.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/failure_analysis_p3_dn.csv)  
> **SHA-256:** `9e66d508445a8aa5e94548771d15ec6c58d9853c343418659a1db573365b13bf`

### Table 4.P3-dn.2: Top 20 Worst Performing Queries for M7 on P3-dn

| Query Drug | N (Rel) | First Rel Drug | Rank ($r_{\min}$) | $S$ | $T$ | $E$ | Targets (Q / Cand)? | SEs (Q / Cand)? |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| CHEMBL1200485 (SORAFENIB) | 1 | CHEMBL3989777 (TRIENTINE) | 1891.0 | 0.016 | 0.000 | 0.000 | Y / N | Y / N |
| CHEMBL2103868 (CABOZANTINIB S-MALATE) | 1 | CHEMBL3989777 (TRIENTINE) | 1891.0 | 0.016 | 0.000 | 0.000 | Y / N | N / N |
| CHEMBL5314934 (VANZACAFTOR) | 1 | CHEMBL689 (MANNITOL) | 1891.0 | 0.011 | 0.000 | 0.000 | Y / N | N / Y |
| CHEMBL1289926 (AXITINIB) | 1 | CHEMBL3989777 (TRIENTINE) | 1885.0 | 0.016 | 0.000 | 0.000 | Y / N | Y / N |
| CHEMBL1289494 (TIVOZANIB) | 1 | CHEMBL3989777 (TRIENTINE) | 1884.0 | 0.015 | 0.000 | 0.000 | Y / N | N / N |
| CHEMBL413 (SIROLIMUS) | 1 | CHEMBL3989777 (TRIENTINE) | 1881.0 | 0.018 | 0.000 | 0.000 | Y / N | N / N |
| CHEMBL1677 (TACRINE) | 1 | CHEMBL1149 (LEVOCARNITINE) | 1879.5 | 0.000 | 0.000 | 0.000 | Y / N | N / N |
| CHEMBL4101807 (VOXELOTOR) | 1 | CHEMBL467 (HYDROXYUREA) | 1877.0 | 0.034 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL4297603 (DEUTIVACAFTOR) | 1 | CHEMBL689 (MANNITOL) | 1876.0 | 0.017 | 0.000 | 0.000 | Y / N | N / Y |
| CHEMBL1289601 (LENVATINIB) | 1 | CHEMBL3989777 (TRIENTINE) | 1875.0 | 0.029 | 0.000 | 0.000 | N / N | N / N |
| CHEMBL4298128 (ELEXACAFTOR) | 1 | CHEMBL689 (MANNITOL) | 1873.0 | 0.037 | 0.000 | 0.000 | Y / N | N / Y |
| CHEMBL1097279 (CANGRELOR TETRASODIUM) | 1 | CHEMBL1200379 (APRACLONIDINE) | 1866.0 | 0.030 | 0.000 | 0.000 | Y / Y | N / Y |
| CHEMBL1201749 (DIFLUPREDNATE) | 1 | CHEMBL1200331 (FOSFOMYCIN TROMETHAMINE) | 1862.0 | 0.027 | 0.000 | 0.000 | Y / Y | N / N |
| CHEMBL1678 (DONEPEZIL) | 1 | CHEMBL1149 (LEVOCARNITINE) | 1848.0 | 0.047 | 0.000 | 0.000 | Y / N | Y / N |
| CHEMBL1200358 (SODIUM PHENYLACETATE) | 1 | CHEMBL296306 (LACTULOSE) | 1844.0 | 0.041 | 0.000 | 0.000 | N / N | N / Y |
| CHEMBL1200448 (TILUDRONATE DISODIUM) | 1 | CHEMBL3989767 (CALCITONIN SALMON) | 1836.0 | 0.034 | 0.000 | 0.000 | N / Y | N / N |
| CHEMBL1201733 (PAZOPANIB) | 1 | CHEMBL3989777 (TRIENTINE) | 1812.0 | 0.030 | 0.000 | 0.000 | Y / N | Y / N |
| CHEMBL1555 (GALANTAMINE) | 1 | CHEMBL1149 (LEVOCARNITINE) | 1779.0 | 0.047 | 0.000 | 0.000 | Y / N | Y / N |
| CHEMBL1220 (TINIDAZOLE) | 1 | CHEMBL2206196 (PAROMOMYCIN) | 1765.0 | 0.025 | 0.000 | 0.031 | N / Y | Y / Y |
| CHEMBL1201110 (VALACYCLOVIR) | 1 | CHEMBL1200453 (DOCOSANOL) | 1758.0 | 0.051 | 0.000 | 0.000 | Y / Y | N / N |

> **Source File:** [`experiments\results\failure_analysis_p3_dn.csv`](file:///C:/Users/hg340/Documents/MoleculAR/experiments/results/failure_analysis_p3_dn.csv)  
> **SHA-256:** `9e66d508445a8aa5e94548771d15ec6c58d9853c343418659a1db573365b13bf`
