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
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

from .misc import isoformat_digits

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

    def __str__(self, table_format=False, table_header=False):
        if table_header:
            fmt = " {:>5s} | {:^7s} | {:^7s} | {:^6s} | {:^6s} | {:^6s}"
            s = fmt.format('id', 'd_km', 'delta', 'p_tt', 's_tt', 's_p')
        elif table_format:
            fmt = " {:5d} | {:7.2f} | {:7.4f} | {:6.2f} | {:6.2f} | {:6.2f}"
            s = fmt.format(self.id, self.d_km, self.delta, self.p_tt,
                          self.s_tt, self.s_p)
        else:
            s= self.__repr__()
        return s


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

    def __str__(self, time_digits=0, table_format=False, table_header=False,
                include_times=False):
        if table_header:
            fmt = " {:>5s} | {:^5s} | {:^5s} | {:^3s} | {:^8s} | {:^9s} | {:^9s}"
            s = fmt.format('id', 'sta', 'net', 'loc', 'latitude', 'longitude',
                           'elevation')
            if include_times:
                s += " | {:^19s} | {:^19s}".format('starttime', 'endtime')
        elif table_format:
            fmt = " {:5d} | {:^5s} | {:^5s} | {:^3s} | {:8.2f} | {:9.2f} | {:9.0f}"
            s = fmt.format(self.id, self.sta, self.net, self.loc,
                          self.latitude, self.longitude, self.elevation)
            if include_times:
                s += " | {:^19s} | {:^19s}".format(
                    isoformat_digits(self.starttime, time_digits),
                    isoformat_digits(self.endtime, time_digits))
                    
        else:
            s= "Station1D('{}', '{}', '{}', {}, {}, {})".format(
               self.sta, self.net, self.loc, self.latitude, self.longitude,
               self.elevation)
            if include_times:
                s += ": <{}, {}>".format(
                    isoformat_digits(self.starttime, time_digits),
                    isoformat_digits(self.endtime, time_digits))
        return s


def make_tt1D_session(url):
    """Make a ORM session from a sqlalchemy database URL"""
    engine = create_engine(url, echo=False)
    BaseTT1D.metadata.create_all(engine)
    session = sessionmaker(bind=engine)
    return session()


def str_station_table(session, order_by='sta'):
    """
    Print Station1D table from traveltime_station session

    :param session: sqlalchemy session containing Station1D method
    :param order_by: 'sta' or 'id'
    """
    basis = Station1D
    order_basis = basis.sta
    if order_by == 'id':
        order_basis = basis.id
    if session.query(basis).count():
        print('Stations are:')
        first_time = True
        for obj in session.query(basis).\
                order_by(order_basis):
            if first_time:
                print(obj.__str__(table_header=True))
                first_time = False
            print(obj.__str__(table_format=True))
    else:
        print('No Stations')


def str_traveltime_table(session, order_by='d_km'):
    """
    Print Traveltime table from traveltime_station session

    :param order_by: 'd_km' or 'id'
    """
    basis = TTtable1D
    order_basis = basis.d_km
    if order_by == 'id':
        order_basis = basis.id
    if session.query(basis).count():
        print('TravelTimes are:')
        first_time = True
        for obj in session.query(basis).\
                order_by(order_basis):
            if first_time:
                print(obj.__str__(table_header=True))
                first_time = False
            print(obj.__str__(table_format=True))
    else:
        print('No TravelTimes')
