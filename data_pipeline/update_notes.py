import re

with open("methodology_saga_notes.md", "r", encoding="utf-8") as f:
    text = f.read()

old_text_pattern = re.compile(r"## DrugBank Open-Access Expansion.*?hits this ceiling perfectly\.", re.DOTALL)

new_text = """## Phase 2: Open-Access Database Expansion (Deferred)
*   **DrugBank Assessment:** Investigated incorporating DrugBank's open-access dataset to expand the pipeline. However, the open-access tier is severely license-limited and contributed zero genuinely new, actionable records with usable structures. It was only useful as a provenance-tag reference.
*   **PubChem Merge Assessment:** A parallel attempt to merge FDA-approved structures from PubChem yielded a potential addition of only 23 genuinely new drugs (a ~1% dataset growth).
*   **Deferred Decision:** Given the marginal dataset growth (~1%), finalizing the 23-drug insert does not justify the cost of re-validating the full algorithmic suite (re-running the strict 8-case benchmark, re-calculating all hub-bias percentile distributions, and executing a full latency test) at this specific stage of the project.
*   **Baseline Lock:** The PubChem merge is explicitly documented as a deferred decision. The dataset remains strictly locked at 2,267 unique drugs as the stable, validated baseline for all subsequent ML and pipeline development."""

text = old_text_pattern.sub(new_text, text)

with open("methodology_saga_notes.md", "w", encoding="utf-8") as f:
    f.write(text)

print("Updated methodology_saga_notes.md")
