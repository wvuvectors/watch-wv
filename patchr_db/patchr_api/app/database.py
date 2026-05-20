# app/database.py

# Imports
from sqlalchemy import create_engine  # Core SQLAlchemy engine for DB connection
from sqlalchemy.orm import sessionmaker, declarative_base  # ORM session factory and base class
from sqlalchemy.exc import SQLAlchemyError
from typing import Generator 
from app.core.config import settings  # Import general settings

# Build the full database URL from settings
DB_URL = (
    f"mysql+pymysql://{settings.DB_USER}:{settings.DB_PASSWORD}"
    f"@{settings.DB_HOST}:{settings.DB_PORT}/{settings.DB_NAME}"
)

# Create SQLAlchemy engine 
engine = create_engine(
    DB_URL,
    pool_pre_ping=True,  # Automatically checks connections are alive before using them
    future=True  # Enables SQLAlchemy 2.0 style usage (recommended)
)

# Creates db sessions used in dependency injection in routers
SessionLocal = sessionmaker(
    bind=engine,
    autocommit=False,  # Sessions require explicit commit
    autoflush=False    # Sessions won’t auto-flush changes until commit
)

# Declarative base: base class for all ORM models
Base = declarative_base()

# Dependency function for FastAPI, provides db session to FastAPI route handlers via Depends
def get_db() -> Generator:
    """
    Usage in routers:
        @router.get("/samples/")
        def get_samples(db: Session = Depends(get_db)):
            ...
    """
    db = SessionLocal()
    try:
        yield db  # Provide the session to the route
    except SQLAlchemyError:
        db.rollback()  # Rollback if any error occurs during the session
        raise
    finally:
        db.close()  # Always close the session
