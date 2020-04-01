"""
Define traveltime tables for sqlalchemy

Station is same for 1D and 3D
"""
from sqlalchemy import (Column, String, Float, DateTime, Integer)
from sqlalchemy.ext.declarative import declarative_base

BaseTT3D = declarative_base()


class TTtable3D(BaseTT3D):
    """
    3D travel time table class and sqlalchemy table

    Columns are:
        id:    primary key
        sta:   station code
        sgid:  source grid ID
        d_km:  Distance in kilometers
        delta: Distance in degrees
        p:     P travel-time
        s:     S travel-time
        s_p:   S-P time
        pn:    Pn travel-time
        sn:    Sn travel-time
        sn_pn: Sn-Pn time

    What happens if there is no Pn or Sn phase for a given distance?
    """
    __tablename__ = 'traveltimes'
    id = Column(Integer, primary_key=True)
    sta = Column(String(5))
    sgid = Column(Integer)
    d_km = Column(Float)
    delta = Column(Float)
    p = Column(Float)
    s = Column(Float)
    s_p = Column(Float)
    pn = Column(Float)
    sn = Column(Float)
    sn_pn = Column(Float)

    def __init__(self, sta, sgid, d_km, delta, p_tt, s_tt, s_p,
                 pn_tt, sn_tt, sn_pn):
        # Does there need to be a different row for each station if all
        # we want is to allow varying depths?
        self.d_km = d_km
        self.sta = sta
        self.sgid = sgid
        self.delta = delta
        self.p = p_tt
        self.s = s_tt
        self.s_p = s_p
        self.pn = pn_tt
        self.sn = sn_tt
        self.sn_pn = sn_pn

    def __repr__(self):
        fmt = "TTtable3D('{}', {:d}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, "
        fmt += "{:.2f}, {:.2f}, {:.2f}, {:.2f})"
        return fmt.format(self.sta, self.sgid, self.d_km, self.delta,
                          self.p, self.s, self.s_p,
                          self.pn, self.sn, self.sn_pn)


class Station3D(BaseTT3D):
    """
    3D Station class and sqlalchemy table

    Columns are:
        id:        primary key
        sta:       station code
        net:       network code
        loc:       location code
        latitude:  station latitude
        longitude: station longitude
        elevation: station elevation
        starttime: station starttime
        endtime:   station endtime
    """
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
        return "Station3D('{}', '{}', '{}', {}, {}, {}): <{} {}>".format(
               self.sta, self.net, self.loc, self.latitude, self.longitude,
               self.elevation, self.starttime, self.endtime)


class SourceGrids(BaseTT3D):
    """
    3D traveltime table point class and sqlalchemy table

    Columns are:
        id:    primary key
        latitude:
        longitude:
        depth:  (km)
    """
    __tablename__ = "sourcegrids"
    id = Column(Integer, primary_key=True)
    latitude = Column(Float)
    longitude = Column(Float)
    depth = Column(Float)

    def __init__(self, latitude, longitude, depth):
        """
        :param latitude:
        :param longitude:
        :param depth: km
        """
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth

    def __repr__(self):
        return "SourceGrids({}, {}, {})".format(self.latitude,
                                                self.longitude, self.depth)
