# backend/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from api.run_beam import router as beam_router  # import the router
from api.BeamDeflection import solve_beam  # relative import works because of __init__.py


app = FastAPI(title="Beam Deformation API")

# CORS for local frontend
# CORS setup...
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(beam_router, prefix="/api")  # <- include your run_beam router

@app.post("/api/beam")
async def beam_solver(params: dict):
    L = params.get("length", 1.0)
    E = params.get("E", 2.1e11)
    I = params.get("I", 1e-6)
    load_pos = params.get("loadPos", L / 2)
    load_val = params.get("loadVal", -1000)
    constraintA = params.get("constraintA", "PIN")
    constraintB = params.get("constraintB", "PIN")

    xs, ys = solve_beam(L, E, I, [(load_pos, load_val)], [constraintA, constraintB])
    return {"x": xs, "y": ys}

@app.get("/")
def root():
    return {"message": "Beam deformation API is running!"}
