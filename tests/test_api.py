"""
FastAPI integration tests for MoleculAR API endpoints.
Tests /api/rank, /api/target-search, and /api/structure.
"""

import sys
from pathlib import Path

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pytest
from fastapi.testclient import TestClient

from server import app
from pipeline import drug_features, get_repurposing_candidates

client = TestClient(app)

QUERY_DRUGS_8 = [
    "SILDENAFIL",
    "MINOXIDIL",
    "METFORMIN",
    "THALIDOMIDE",
    "DULOXETINE",
    "GABAPENTIN",
    "BUPROPION",
    "RALOXIFENE",
]


def test_rank_sildenafil_fields_and_types():
    """(a) /api/rank for Sildenafil returns all new fields with the right types."""
    resp = client.post("/api/rank", json={"drug_name": "SILDENAFIL", "top_n": 10})
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # Envelope fields and types
    assert data["query_drug"] == "SILDENAFIL"
    assert isinstance(data["query_hubness_percentile"], (float, int))
    assert data["query_has_targets"] is True
    assert data["query_has_side_effects"] is True
    assert isinstance(data["results"], list)
    assert len(data["results"]) == 10

    # Per-candidate fields and types
    for cand in data["results"]:
        assert isinstance(cand["candidate"], str)
        assert isinstance(cand["score"], (float, int))
        assert isinstance(cand["struct_sim"], (float, int))
        assert isinstance(cand["target_sim"], (float, int))
        assert isinstance(cand["side_effect_sim"], (float, int))
        assert isinstance(cand["shared_targets"], list)
        assert isinstance(cand["shared_side_effects"], list)
        assert isinstance(cand["indication"], str)
        # Additive fields
        assert cand["hubness_percentile"] is None or isinstance(cand["hubness_percentile"], (float, int))
        assert isinstance(cand["has_targets"], bool)
        assert isinstance(cand["has_side_effects"], bool)
        assert isinstance(cand["is_duplicate_ingredient"], bool)


def test_rank_order_and_scores_equal_pipeline():
    """(b) The returned order and scores equal a direct call to get_repurposing_candidates."""
    resp = client.post("/api/rank", json={"drug_name": "SILDENAFIL", "top_n": 10})
    assert resp.status_code == 200
    api_results = resp.json()["results"]

    direct_results = get_repurposing_candidates("SILDENAFIL", top_n=10)
    assert len(api_results) == len(direct_results)

    for api_cand, direct_cand in zip(api_results, direct_results):
        assert api_cand["candidate"] == direct_cand["candidate"]
        assert pytest.approx(api_cand["score"], abs=1e-9) == direct_cand["score"]
        assert pytest.approx(api_cand["struct_sim"], abs=1e-9) == direct_cand["struct_sim"]
        assert pytest.approx(api_cand["target_sim"], abs=1e-9) == direct_cand["target_sim"]
        assert pytest.approx(api_cand["side_effect_sim"], abs=1e-9) == direct_cand["side_effect_sim"]


def test_drug_no_hubness_returns_null():
    """(c) A drug with no hubness value returns null from API, and candidate score equals raw score without penalty."""
    from pipeline import resolve_chembl_id

    q_cid = resolve_chembl_id("SILDENAFIL")
    cand_cid = resolve_chembl_id("VARDENAFIL")

    orig_q_hub = drug_features[q_cid].get("hubness_percentile")
    orig_cand_hub = drug_features[cand_cid].get("hubness_percentile")

    try:
        # Temporarily set hubness_percentile to None for both query and candidate
        drug_features[q_cid]["hubness_percentile"] = None
        drug_features[cand_cid]["hubness_percentile"] = None

        resp = client.post("/api/rank", json={"drug_name": "SILDENAFIL", "top_n": 5})
        assert resp.status_code == 200
        data = resp.json()

        # (1) Query drug with no hubness returns null in JSON
        assert data["query_hubness_percentile"] is None

        # (2) Candidate drug with no hubness returns null in JSON and score equals raw score
        vardenafil_cand = next((c for c in data["results"] if c["candidate"] == "VARDENAFIL"), None)
        assert vardenafil_cand is not None, "VARDENAFIL should be in top results"
        assert vardenafil_cand["hubness_percentile"] is None

        # Raw score: 0.4 * struct + 0.4 * target + 0.2 * side_effect
        raw_score = (
            0.4 * vardenafil_cand["struct_sim"]
            + 0.4 * vardenafil_cand["target_sim"]
            + 0.2 * vardenafil_cand["side_effect_sim"]
        )
        assert pytest.approx(vardenafil_cand["score"], abs=1e-9) == raw_score
    finally:
        # Always restore original hubness values even on failure
        drug_features[q_cid]["hubness_percentile"] = orig_q_hub
        drug_features[cand_cid]["hubness_percentile"] = orig_cand_hub


def test_structure_known_and_unknown():
    """(d) /api/structure returns 200 image/svg+xml for a known drug and 404 for an unknown one."""
    # Known drug
    resp_known = client.get("/api/structure", params={"name": "SILDENAFIL"})
    assert resp_known.status_code == 200
    assert "image/svg+xml" in resp_known.headers["content-type"]
    assert resp_known.content.startswith(b"<?xml") or b"<svg" in resp_known.content
    assert b"</svg>" in resp_known.content

    # Unknown drug
    resp_unknown = client.get("/api/structure", params={"name": "NONEXISTENT_DRUG_XYZ_123"})
    assert resp_unknown.status_code == 404


def test_structure_every_candidate_for_8_query_drugs():
    """(e) Every candidate name returned for the 8 query drugs returns 200 from /api/structure."""
    for drug in QUERY_DRUGS_8:
        rank_resp = client.post("/api/rank", json={"drug_name": drug, "top_n": 10})
        assert rank_resp.status_code == 200, f"Failed rank for {drug}"
        candidates = rank_resp.json()["results"]
        assert len(candidates) > 0

        for cand in candidates:
            cname = cand["candidate"]
            struct_resp = client.get("/api/structure", params={"name": cname})
            assert struct_resp.status_code == 200, f"Structure 404 for candidate {cname} from query {drug}"
            assert "image/svg+xml" in struct_resp.headers["content-type"]
            assert b"<svg" in struct_resp.content


def test_target_search_known_and_unknown():
    """(f) /api/target-search for P35354 returns results and for a made-up ID returns empty list with query_has_targets false."""
    # Known target (PTGS2 / COX-2)
    resp_known = client.post("/api/target-search", json={"target_id": "P35354", "top_n": 10})
    assert resp_known.status_code == 200
    data_known = resp_known.json()
    assert data_known["query_target"] == "P35354"
    assert data_known["query_has_targets"] is True
    assert data_known["query_has_side_effects"] is False
    assert len(data_known["results"]) > 0

    # Made-up target
    resp_made_up = client.post("/api/target-search", json={"target_id": "NONEXISTENT_TARGET_9999", "top_n": 10})
    assert resp_made_up.status_code == 200
    data_made_up = resp_made_up.json()
    assert data_made_up["query_target"] == "NONEXISTENT_TARGET_9999"
    assert data_made_up["query_has_targets"] is False
    assert data_made_up["query_has_side_effects"] is False
    assert data_made_up["results"] == []


def test_duplicate_pairs_count_and_detection():
    """(g) duplicate_pairs.json has exactly 107 pairs and API sets is_duplicate_ingredient accurately."""
    import json
    dup_json_path = REPO_ROOT / "duplicate_pairs.json"
    assert dup_json_path.exists(), "duplicate_pairs.json must exist"

    with open(dup_json_path, "r", encoding="utf-8") as f:
        dup_data = json.load(f)

    # Exactly 107 pairs recorded
    assert dup_data["metadata"]["count"] == 107
    assert len(dup_data["pairs"]) == 107

    # Test Cetirizine rank: Levocetirizine Dihydrochloride is a known duplicate pair (first_block_inchikey)
    resp = client.post("/api/rank", json={"drug_name": "CETIRIZINE", "top_n": 10})
    assert resp.status_code == 200
    results = resp.json()["results"]

    levo = next((c for c in results if c["candidate"] == "LEVOCETIRIZINE DIHYDROCHLORIDE"), None)
    assert levo is not None, "LEVOCETIRIZINE DIHYDROCHLORIDE should be in top results for CETIRIZINE"
    assert levo["is_duplicate_ingredient"] is True

    hydroxyzine = next((c for c in results if c["candidate"] == "HYDROXYZINE"), None)
    assert hydroxyzine is not None
    assert hydroxyzine["is_duplicate_ingredient"] is False

    # Target search candidates should all have is_duplicate_ingredient False
    t_resp = client.post("/api/target-search", json={"target_id": "P35354", "top_n": 5})
    assert t_resp.status_code == 200
    for cand in t_resp.json()["results"]:
        assert cand["is_duplicate_ingredient"] is False


def test_health_ok_returns_200():
    """(h) /api/health returns 200 OK with database hash verification and table row counts."""
    resp = client.get("/api/health")
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "ok"
    assert data["database"]["sha256_verified"] is True
    assert data["database"]["expected_sha256"] == "D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7"
    assert data["database"]["row_counts"] == {
        "drugs": 2267,
        "targets": 7068,
        "side_effects": 119227,
        "crosswalk": 2267,
    }


def test_health_mismatch_returns_503(monkeypatch):
    """(i) /api/health returns HTTP 503 with the same JSON schema when database hash verification fails."""
    import server as srv

    patched_payload = {
        "status": "sha256_mismatch",
        "database": {
            "expected_sha256": srv.EXPECTED_DB_SHA256,
            "sha256_verified": False,
            "row_counts": srv._HEALTH_PAYLOAD["database"]["row_counts"],
        },
    }
    monkeypatch.setattr(srv, "_HEALTH_PAYLOAD", patched_payload)

    resp = client.get("/api/health")
    assert resp.status_code == 503
    assert resp.json() == patched_payload


