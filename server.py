import sqlite3
from functools import lru_cache
from typing import List, Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response
from pydantic import BaseModel, Field

# Import the existing production pipeline
from pipeline import (
    DB_PATH,
    all_target_ids,
    drug_features,
    get_repurposing_candidates,
    resolve_chembl_id,
)

app = FastAPI(title="MoleculAR API", version="1.0.0")

# Configure CORS: Allow only localhost origins with credentials
ALLOWED_ORIGINS = [
    "http://localhost:5173",
    "http://127.0.0.1:5173",
    "http://localhost:3000",
    "http://127.0.0.1:3000",
    "http://localhost:8000",
    "http://127.0.0.1:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------------------------
# Request / Response schemas
# ---------------------------------------------------------------------------

class RankRequest(BaseModel):
    drug_name: str = Field(..., description="The name of the query drug (e.g., 'SILDENAFIL')")
    top_n: int = Field(10, description="Number of top candidates to return")


class CandidateResult(BaseModel):
    candidate: str
    score: float
    struct_sim: float
    target_sim: float
    side_effect_sim: float
    shared_targets: List[str]
    shared_side_effects: List[str]
    indication: str
    # Additive fields (pass-through from drug_features; do not affect score/ordering)
    hubness_percentile: Optional[float] = None
    has_targets: bool
    has_side_effects: bool
    is_duplicate_ingredient: bool = False


class RankResponse(BaseModel):
    query_drug: str
    query_hubness_percentile: Optional[float] = None
    query_has_targets: bool
    query_has_side_effects: bool
    results: List[CandidateResult]


class TargetRequest(BaseModel):
    target_id: str = Field(..., description="The UniProt ID of the target protein (e.g., 'P35354')")
    top_n: int = Field(10, description="Number of top candidates to return")


class TargetResponse(BaseModel):
    query_target: str
    query_has_targets: bool
    query_has_side_effects: bool
    results: List[CandidateResult]


# ---------------------------------------------------------------------------
# Endpoints: rank and target-search (unchanged behaviour, additive envelope)
# ---------------------------------------------------------------------------

@app.post("/api/rank", response_model=RankResponse)
def rank_candidates(request: RankRequest):
    """
    Rank repurposing candidates for a given drug.
    Applies hand-weights 0.4/0.4/0.2 and hub-percentile correction (unchanged).
    Additive envelope fields: query_hubness_percentile, query_has_targets,
    query_has_side_effects. Per-candidate: hubness_percentile, has_targets,
    has_side_effects.
    """
    query = request.drug_name.upper()

    if query not in drug_features:
        raise HTTPException(status_code=404, detail=f"Drug '{query}' not found in the dataset.")

    try:
        candidates = get_repurposing_candidates(query, top_n=request.top_n)

        # Envelope: query drug properties (read-only from drug_features cache)
        q_cid = resolve_chembl_id(query)
        q_feats = drug_features[q_cid]

        return RankResponse(
            query_drug=query,
            query_hubness_percentile=q_feats.get("hubness_percentile"),
            query_has_targets=len(q_feats.get("targets", set())) > 0,
            query_has_side_effects=len(q_feats.get("side_effects", set())) > 0,
            results=candidates,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/target-search", response_model=TargetResponse)
def target_search(request: TargetRequest):
    """
    Rank repurposing candidates for a given UniProt target.
    Returns 200 OK with an empty array if the target is not found in any drug.
    query_has_targets is True only if at least one drug in the database has that target.
    """
    from pipeline import get_target_candidates

    query = request.target_id.upper()

    try:
        candidates = get_target_candidates(query, top_n=request.top_n)
        return TargetResponse(
            query_target=query,
            query_has_targets=query in all_target_ids,
            query_has_side_effects=False,
            results=candidates,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ---------------------------------------------------------------------------
# Structure endpoint: GET /api/structure?name=<drug_name>
# Read-only. Looks up SMILES from DB by chembl_id. Returns image/svg+xml.
# Results are cached in memory. Never accepts SMILES from the client.
# ---------------------------------------------------------------------------

# In-memory SVG cache: chembl_id -> SVG bytes
_structure_svg_cache: dict[str, bytes] = {}


@lru_cache(maxsize=1024)
def _lookup_smiles(chembl_id: str) -> str | None:
    """
    Look up canonical SMILES from drugrepo.db by chembl_id.
    Cached via lru_cache; read-only via URI mode=ro.
    """
    conn = sqlite3.connect(f"file:{DB_PATH.as_posix()}?mode=ro", uri=True)
    cur = conn.cursor()
    cur.execute(
        "SELECT smiles FROM drugs WHERE chembl_id = ?",
        (chembl_id,),
    )
    row = cur.fetchone()
    conn.close()
    return row[0] if row else None


@app.get("/api/structure")
def get_structure(name: str = Query(..., max_length=100, description="Drug name (e.g., SILDENAFIL)")):
    """
    Return a 300×200 SVG 2D structure drawing for the named drug.
    Resolves the name via resolve_chembl_id and looks up SMILES by chembl_id in drugrepo.db.
    Returns 404 if not found.
    SVG has a white background and is safe to embed via <img src="...">.
    Results are cached in memory for the server lifetime.
    Never accepts SMILES from the client.
    """
    name_clean = name.strip()
    cid = resolve_chembl_id(name_clean)
    if not cid:
        raise HTTPException(status_code=404, detail=f"Drug '{name}' not found in the database.")

    if cid in _structure_svg_cache:
        return Response(content=_structure_svg_cache[cid], media_type="image/svg+xml")

    smiles = _lookup_smiles(cid)
    if not smiles:
        raise HTTPException(status_code=404, detail=f"Drug '{name}' not found in the database.")

    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(
            status_code=404,
            detail=f"Could not parse a valid structure for '{name}'.",
        )

    # 300×200 px SVG with white background (readable in both light and dark mode)
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg_text = drawer.GetDrawingText()
    svg_bytes = svg_text.encode("utf-8")

    _structure_svg_cache[cid] = svg_bytes
    return Response(content=svg_bytes, media_type="image/svg+xml")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
