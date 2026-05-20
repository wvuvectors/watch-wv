# app/main.py

from pathlib import Path # defines frontend_path
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

# Importing router objects for each table router 
from app.routers.samples import router as samples_router # Import router object directly
from app.routers.concentration import router as concentration_router

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

# Mount the frontend directory
frontend_path = Path(__file__).parent / "static"

# Serve frontend files at /static URL path
# For example, index.html will be accessible at https://127.0.0.1:8000/static/index.html
app.mount("/static", StaticFiles(directory=frontend_path), name="static")