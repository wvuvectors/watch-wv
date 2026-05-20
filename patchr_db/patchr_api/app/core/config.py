# app/core/config.py
# Handles app-wide settings, contains db credentials and env-dependent settings

from pydantic_settings import BaseSettings  # Pydantic v2 settings management

# App settings
class Settings(BaseSettings):
    # Database credentials
    DB_USER: str = "labuser"
    DB_PASSWORD: str = "password123"
    DB_HOST: str = "localhost"
    DB_PORT: int = 3306
    DB_NAME: str = "patchr_db"
    
    FRONTEND_PATH: str = "patchr_frontend"
    
    # Load values from .env file if present
    class Config:
        env_file = ".env"

# Global settings object for importing    
settings = Settings()