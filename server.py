from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List

# Import the existing production pipeline
from pipeline import get_repurposing_candidates, drug_features

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

# Request Schema
class RankRequest(BaseModel):
    drug_name: str = Field(..., description="The name of the query drug (e.g., 'SILDENAFIL')")
    top_n: int = Field(10, description="Number of top candidates to return")

# Response Schema components
class CandidateResult(BaseModel):
    candidate: str
    score: float
    struct_sim: float
    target_sim: float
    side_effect_sim: float
    shared_targets: List[str]
    shared_side_effects: List[str]
    indication: str

class RankResponse(BaseModel):
    query_drug: str
    results: List[CandidateResult]

@app.post("/api/rank", response_model=RankResponse)
def rank_candidates(request: RankRequest):
    """
    Rank repurposing candidates for a given drug using the finalized Phase 4.4 Hand-Weighted Pipeline.
    Includes percentile-rank hub bias correction.
    """
    query = request.drug_name.upper()
    
    # Verify drug exists in dataset
    if query not in drug_features:
        raise HTTPException(status_code=404, detail=f"Drug '{query}' not found in the dataset.")
        
    try:
        # get_repurposing_candidates applies hand-weights (0.4/0.4/0.2) and hub-correction by default
        candidates = get_repurposing_candidates(query, top_n=request.top_n)
        
        return RankResponse(
            query_drug=query,
            results=candidates
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

class TargetRequest(BaseModel):
    target_id: str = Field(..., description="The UniProt ID of the target protein (e.g., 'P35354')")
    top_n: int = Field(10, description="Number of top candidates to return")

class TargetResponse(BaseModel):
    query_target: str
    results: List[CandidateResult]

@app.post("/api/target-search", response_model=TargetResponse)
def target_search(request: TargetRequest):
    """
    Rank repurposing candidates for a given target using a virtual drug feature vector.
    Applies the same hub-correction penalty as Mode A.
    Returns 200 OK with an empty array if the target is not found in any drug.
    """
    from pipeline import get_target_candidates
    
    query = request.target_id.upper()
    
    try:
        candidates = get_target_candidates(query, top_n=request.top_n)
        return TargetResponse(
            query_target=query,
            results=candidates
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
