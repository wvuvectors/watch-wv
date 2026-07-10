# app/models/wwtp.py

from sqlalchemy import Column, String, Float
from app.database import Base

class WWTP(Base):
    __tablename__ = "wwtp"
    
    wwtp_id = Column(String(100), primary_key=True)
    wwtp_site_id = Column(String(30), nullable=True)
    wwtp_common_name = Column(String(50), nullable=True)
    wwtp_authority_name = Column(String(100), nullable=True)
    wwtp_counties_served = Column(String(20), nullable=True)
    wwtp_epaid_id = Column(String(20), nullable=True)
    wwtp_cwns_id = Column(String(20), nullable=True)
    wwtp_capacity_mgd = Column(Float, nullable=True)
    wwtp_population_served = Column(Float, nullable=True)