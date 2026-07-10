# app/schemas/wwtp.py

from pydantic import BaseModel
from typing import Optional

class WWTPSchema(BaseModel):
    wwtp_id: str
    wwtp_site_id: Optional[str]
    wwtp_common_name: Optional[str]
    wwtp_authority_name: Optional[str]
    wwtp_counties_served: Optional[str]
    wwtp_epaid_id: Optional[str]
    wwtp_cwns_id: Optional[str]
    wwtp_capacity_mgd: Optional[float]
    wwtp_population_served: Optional[float]
    
    class Config:
        from_attributes = True