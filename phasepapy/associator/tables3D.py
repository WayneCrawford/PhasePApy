"""
Define database tables
"""
from sqlalchemy import Column, String, Float, Integer, Boolean, DateTime
from sqlalchemy.ext.declarative import declarative_base

from ..phasepicker.scnl import SCNL

# Declare mapping
Base = declarative_base()


class Pick(Base):
    """
    Pick class and sqlalchemy table type

    Attributes/Columns are:
        id: table primary key
        sta: station code
        chan: channel code
        net: network code
        loc: location code
        time: pick time (DateTime)
        snr: Pick signal to noise ratio
        phase: Pick phase name ('P' or 'S')
        uncert: Pick uncertainty (seconds?)
        polarity: 'C' or 'D'????
        locate_flag: ??? (Boolean)
        assoc_id: id of linked Associated()
        modified_id: id of linked PickModified()
        t_create: creation time of this pick (DateTime)
    """
    __tablename__ = "picks"
    id = Column(Integer, primary_key=True)
    sta = Column(String(5))
    chan = Column(String(3))
    net = Column(String(2))
    loc = Column(String(2))
    time = Column(DateTime)
    snr = Column(Float)
    phase = Column(String(1))
    uncert = Column(Float)
    polarity = Column(String(1))
    locate_flag = Column(Boolean)
    assoc_id = Column(Integer)
    modified_id = Column(Integer)
    t_create = Column(DateTime)

    def __init__(self, scnl, picktime, polarity, snr, uncert, t_create):
        """
        :param scnl: station, component, network, location codes
        :type scnl: class SCNL
        :param picktime: phase pick time (DateTime?)
        :param polarity: ???
        :param snr: Signal to nois ratio
        :param uncert: Pick uncertainty? (seconds?)
        :param t_create: time that pick was created (DateTime?)
        """
        self.sta = scnl.station
        self.chan = scnl.channel
        self.net = scnl.network
        self.loc = scnl.location
        self.time = picktime
        self.snr = snr
        self.uncert = uncert
        self.modified_id = None
        self.phase = None
        self.polarity = polarity
        self.locate_flag = None
        self.assoc_id = None
        self.t_create = t_create

    def __repr__(self):
        return 'Pick({}, "{}", "{}", "{}", "{}")'.format(
            SCNL(self.sta, self.chan, self.net, self.loc),
            self.time.isoformat("T"), self.phase,
            self.modified_id, self.assoc_id)


class PickModified(Base):
    """
    PickModified class and sqlalchemy table type
    """
    __tablename__ = "picks_modified"
    id = Column(Integer, primary_key=True)
    sta = Column(String(5))
    chan = Column(String(3))
    net = Column(String(2))
    loc = Column(String(2))
    time = Column(DateTime)
    phase = Column(String(1))
    error = Column(Float)
    locate_flag = Column(Boolean)
    assoc_id = Column(Integer)

    def __init__(self, sta, chan, net, loc, picktime, phase, error, assoc_id):
        """
        :param sta: station name
        :param chan: channel name
        :param net: network name
        :param loc: location code
        :param picktime: phase pick time (DateTime?)
        :paraam phase: "P" or "S"?  (Pn and Sn too?)
        :param error: Pick uncertainty? (seconds?)
        :param assoc_id): ID of linked Associated()
        """
        self.sta = sta
        self.chan = chan
        self.net = net
        self.loc = loc
        self.time = picktime
        self.phase = phase
        self.error = error
        self.locate_flag = None
        self.assoc_id = assoc_id

    def __repr__(self):
        return 'PickModified({}, {}, {}, "{}", {}, {}, {}, {})'.format(
            self.sta, self.chan, self.net, self.loc, self.time.isoformat("T"),
            self.phase, self.error, self.assoc_id)


class Candidate(Base):
    """
    Candidate class and sqlalchemy table type
    """
    __tablename__ = "candidate"
    id = Column(Integer, primary_key=True)
    ot = Column(DateTime)
    sta = Column(String(5))
    d_km = Column(Float)
    delta = Column(Float)
    weight = Column(Float)
    # P and S travel times are not completely necessary because they can be
    # calculated, but simpler to save here
    tp = Column(DateTime)
    p_modified_id = Column(Integer)    # modified pick ID
    ts = Column(DateTime)
    s_modified_id = Column(Integer)    # modified pick ID
    locate_flag = Column(Boolean)
    assoc_id = Column(Integer)

    def __init__(self, ot, sta, d_km, delta, tp, p_modified_id, ts,
                 s_modified_id):
        """
        :param ot: origin time (DateTime?)
        :param sta: station name?
        :param d_km: distance in km station-origin??
        :param delta:
        :param tp: P-wave traveltime?
        :param p_modified_id: id of corresponding PickModified?
        :param ts: S-wave traveltime?
        :param s_modified_id: id of corresponding PickModified?
        """
        self.ot = ot
        self.sta = sta
        self.d_km = d_km
        self.delta = delta
        self.weight = None
        self.tp = tp
        self.ts = ts
        self.p_modified_id = p_modified_id
        self.s_modified_id = s_modified_id
        self.locate_flag = None
        self.assoc_id = None

    def __repr__(self):
        s = 'Candidate("{}", "{}", {:.2f}, {:.2f}, '.format(
            self.ot.isoformat("T"), self.sta, self.d_km, self.delta)
        s += '{:g}, {:d}, {:g}, {:d})'.format(
            self.tp, self.p_modified_id, self.ts, self.s_modified_id)
        return s

    # def __str__(self):
    #     return "Candidate Event <%s %s %.2f %.2f %d %d>" %
    #            self.ot.isoformat("T"), self.sta, self.d_km, self.delta,
    #            self.p_modified_id, self.s_modified_id)
    def set_assoc_id(self, assoc_id, session, FT):
        self.assoc_id = assoc_id
        self.locate_flag = FT
        # Assign phases to modified picks

        # Actually only one pick_p and pick_s
        pick_p = session.query(PickModified).\
            filter(PickModified.id == self.p_modified_id)
        for pick in pick_p:
            pick.phase = 'P'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT

        pick_s = session.query(PickModified).\
            filter(PickModified.id == self.s_modified_id)
        for pick in pick_s:
            pick.phase = 'S'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT

        # Assign the phases to picks contribute to a modified picks
        picks_p = session.query(Pick).\
            filter(Pick.modified_id == self.p_modified_id).all()
        for pick in picks_p:
            pick.phase = 'P'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT

        picks_s = session.query(Pick).\
            filter(Pick.modified_id == self.s_modified_id).all()
        for pick in picks_s:
            pick.phase = 'S'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT


class Associated(Base):
    """
    Associated class and sqlalchemy table type
    """
    __tablename__ = "associated"
    id = Column(Integer, primary_key=True)
    ot = Column(DateTime)
    ot_uncert = Column(Float)
    latitude = Column(Float)
    longitude = Column(Float)
    depth = Column(Float)
    rms = Column(Float)
    nsta = Column(Integer)
    t_create = Column(DateTime)
    t_update = Column(DateTime)

    def __init__(self, ot, ot_uncert, latitude, longitude, depth, rms, nsta,
                 t_create, t_update):
        """
        :param ot: origin time (DateTime?)
        :param ot_uncert: origin time uncertainty (seconds?)
        :param latitude: origin latitude
        :param longitude: origin longitude
        :param loc_uncert: location uncertainty (km?)
        :param nsta: number of stations (fitting the origin?)
        :param t_create: time this Association was created? (DateTime?)
        :param t_update:
        """
        self.ot = ot
        self.ot_uncert = ot_uncert
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.rms = rms
        self.nsta = nsta
        self.t_create = t_create
        self.t_update = t_update

    def __repr__(self):
        fmt = "Associated({}, {:.2f}, {:.3f}, {:.3f}, {:.3f}, {:d}, {}, {})"
        return fmt.format(self.ot.isoformat("T"), self.ot_uncert,
                          self.latitude, self.longitude, self.loc_uncert,
                          self.nsta, self.t_create.isoformat("T"),
                          self.t_update.isoformat("T"))
