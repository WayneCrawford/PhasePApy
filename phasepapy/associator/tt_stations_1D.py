"""
Define traveltime tables:
:TTtable1D: id, distance_km, delta_degrees, p_traveltime, s_traveltime,
            s_p_time
:Station1D: id, station, network, location, latitide, longitude, elevation,
            starttime, endtime

Station is same for 1D and 3D

"""
from sqlalchemy import (Column, Integer, Float, DateTime, String)
from sqlalchemy.ext.declarative import declarative_base

BaseTT1D = declarative_base()


class TTtable1D(BaseTT1D):
    """
    1D travel-time table
    """
    __tablename__ = 'traveltimes'
    id = Column(Integer, primary_key=True)
    d_km = Column(Float)     # Distance in kilometers
    delta = Column(Float)    # Distance in degrees
    p_tt = Column(Float)     # P travel-time
    s_tt = Column(Float)     # S travel-time
    s_p = Column(Float)      # S-P time

    def __init__(self, d_km, delta, p_tt, s_tt, s_p):
        """
        Create an entry in the travel-time table for later lookup.

        The travel-time table is stored in a sqlalchemy data table
        d_km=Column(Float) # Distance in kilometers
        delta=Column(Float) # Distance in degrees
        p_tt=Column(Float) # P travel-time
        s_tt=Column(Float) # S travel-time
        s_p=Column(Float) # S-P time
        """
        self.d_km = d_km
        self.delta = delta
        self.p_tt = p_tt
        self.s_tt = s_tt
        self.s_p = s_p

    def __repr__(self):
        return "TTtable1D({:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f})".format(
               self.d_km, self.delta, self.p_tt, self.s_tt, self.s_p)


class Station1D(BaseTT1D):
    __tablename__ = "stations"
    id = Column(Integer, primary_key=True)
    sta = Column(String(5))
    net = Column(String(2))
    loc = Column(String(2))
    latitude = Column(Float)
    longitude = Column(Float)
    elevation = Column(Float)
    starttime = Column(DateTime)
    endtime = Column(DateTime)

    def __init__(self, sta, net, loc, latitude, longitude, elevation):
        self.sta = sta
        self.net = net
        self.loc = loc
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.starttime = None
        self.endtime = None

    def __repr__(self):
        return "Station1D('{}', '{}', '{}', {}, {}, {}): <{} {}>".format(
               self.sta, self.net, self.loc, self.latitude, self.longitude,
               self.elevation, self.starttime, self.endtime)
