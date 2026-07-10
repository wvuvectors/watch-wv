# app/schemas/county.py

from pydantic import BaseModel
from typing import Optional

class CountySchema(BaseModel):
    county_id: str
    county_labcode: Optional[str]
    county_fips: Optional[str]
    county_name: Optional[str]
    county_population: Optional[float]
    
    class Config:
        from_attributes = True