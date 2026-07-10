# app/main.py

from pathlib import Path # defines frontend_path
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

# Importing router objects for each table router 
from app.routers.samples import router as samples_router # Import router object directly
from app.routers.concentration import router as concentration_router
from app.routers.extractions import router as extractions_router
from app.routers.assay import router as assay_router
from app.routers.cbatch import router as cbatch_router
from app.routers.ebatch import router as ebatch_router
from app.routers.abatch import router as abatch_router
from app.routers.results import router as results_router
from app.routers.location import router as location_router
from app.routers.county import router as county_router
from app.routers.wwtp import router as wwtp_router

# Create FastAPI app
app = FastAPI(
    title="PATCHR API",
    description="Backend API for PATCHR database"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins = ["*"], # restrict later
    allow_credentials = True,
    allow_methods = ["*"],
    allow_headers = ["*"],
)

# Root endpoint, verifies API is running
@app.get("/")
def root():
    return {"status": "PATCHR API is running"}

# Include routers 
app.include_router(samples_router)
app.include_router(concentration_router)
app.include_router(extractions_router)
app.include_router(assay_router)
app.include_router(cbatch_router)
app.include_router(ebatch_router)
app.include_router(abatch_router)
app.include_router(results_router)
app.include_router(location_router)
app.include_router(county_router)
app.include_router(wwtp_router)

# Mount the frontend directory
frontend_path = Path(__file__).parent / "static"

# Serve frontend files at /static URL path
# For example, index.html will be accessible at https://127.0.0.1:8000/static/index.html
app.mount("/static", StaticFiles(directory=frontend_path), name="static")