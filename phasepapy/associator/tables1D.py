"""
Define database tables:
:Pick: id, station, channel, network, location, time, SNR, phase,
       uncertainty, polarity, locate_flag, assoc_id, modified_id
       creation_time
:PickModified: id, station, channel, network, location, picktime,
               phase, error, locate_flag, associate_id
:Candidate: id, origin_time, station, distance_km, delta, weight,
            traveltime_P, modified_P_id, traveltime_S, modified_S_id,
            locate_flag, associate_id
:Associated: id, origin_time, origin_time_uncertainty,
             latitude, longitude, location_uncertainty,
             num_stations, creation_time, update_time
"""

from sqlalchemy import (Column, Integer, String, DateTime, Float,
                        Boolean)
from sqlalchemy.ext.declarative import declarative_base
from ..phasepicker.scnl import SCNL
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
        self.polarity = polarity
        self.snr = snr
        self.uncert = uncert
        self.modified_id = None
        self.phase = None
        self.locate_flag = None
        self.assoc_id = None
        self.t_create = t_create

    def __repr__(self):
        s = "Pick({}, '{}', '{}', '{}', '{}', '{}')".format(
            SCNL([self.sta, self.chan, self.net, self.loc]).__repr__(),
            self.time.isoformat(),
            self.polarity, self.snr, self.uncert,
            self.t_create.isoformat())
        if self.phase:
            s += ' + phase="{}"'.format(self.phase)
        if self.modified_id:
            s += ' + modified_id="{}"'.format(self.modified_id)
        if self.assoc_id:
            s += ' + assoc_id="{}" '.format(self.assoc_id)
        return s

    def __str__(self, time_digits=2, table_format=False, table_header=False,
                include_create_update=False):
        """
        :param time_digits: number of digits after seconds in times
        :param table_format: output format for a table_format
        :param table_header: return table_format header (ignores table_format)
        :param include_create_update: include creation and update times
        """
        assoc_wd, modified_wd, phase_wd = '-', '-', '-'
        if self.assoc_id:
            assoc_wd = f"{self.assoc_id:d}"
        if self.phase:
            phase_wd=self.phase
        if self.modified_id:
            modified_wd = f"{self.modified_id:d}"
        if table_header:
            fmt = " {:>5s} | {:^22s} | {:^14s} | {:^8s} | {:^8s} | {:^8s} | {:^5s} | {:^8s} | {:^8s}"
            s = fmt.format('id', 'time', 'SCNL', 'polarity', 'snr',  'uncert',
                           'phase', 'assoc_id', 'modified_id')
        elif table_format:
            fmt = " {:5d} | {:22s} | {:^14s} | {:8s} | {:8.2f} | {:8.2f} | {:5s} | {:8s} | {:8s}"
            s = fmt.format(self.id, _isoformat_digits(self.time, time_digits),
                          SCNL([self.sta, self.chan, self.net, self.loc]).__str__(),
                          self.polarity, self.snr, self.uncert, phase_wd,
                          assoc_wd, modified_wd)
        else: 
            s = "Pick({}, '{}', '{}', '{}', '{}', '{}')".format(
                SCNL([self.sta, self.chan, self.net, self.loc]),
                _isoformat_digits(self.time, time_digits),
                self.polarity, self.snr, self.uncert,
                _isoformat_digits(self.t_create, time_digits))
            if self.phase:
                s += ' + phase="{}"'.format(self.phase)
            if self.modified_id:
                s += ' + modified_id="{}"'.format(self.modified_id)
            if self.assoc_id:
                s += ' + assoc_id="{}" '.format(self.assoc_id)
        return s


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

    def __init__(self, sta, chan, net, loc, picktime, phase,
                 uncert, assoc_id):
        """
        :param sta: station name
        :param chan: channel name
        :param net: network name
        :param loc: location code
        :param picktime: phase pick time (DateTime?)
        :paraam phase: "P" or "S"?  (Pn and Sn too?)
        :param uncert: Pick uncertainty? (seconds?)
        :param assoc_id): ID of linked Associated()
        """
        self.sta = sta
        self.chan = chan
        self.net = net
        self.loc = loc
        self.time = picktime
        self.phase = phase
        self.error = uncert
        self.locate_flag = None
        self.assoc_id = assoc_id

    def __repr__(self):
        fmt = "PickModified('{}', '{}', '{}', '{}', '{}',"
        fmt += " '{}', '{}', {}, '{}')"
        return fmt.format(self.sta, self.chan, self.net, self.loc,
                          self.time.isoformat(), self.phase, self.error,
                          self.assoc_id)

    def __str__(self, time_digits=2, table_format=False, table_header=False):
        """
        :param time_digits: number of digits after seconds in times
        :param table_format: output format for a table_format
        :param table_header: return table_format header (ignores table_format)
        """
        assoc_wd, phase_wd = '-', '-'
        if self.assoc_id:
            assoc_wd = f"{self.assoc_id:d}"
        if self.phase:
            phase_wd=self.phase
        if table_header:
            fmt = " {:>5s} | {:^22s} | {:^14s} | {:^8s} | {:^5s} | {:^8s}"
            s = fmt.format('id', 'time', 'SCNL', 'uncert', 'phase', 'assoc_id')
        elif table_format:
            fmt = " {:5d} | {:22s} | {:^14s} | {:8.2f} | {:5s} | {:8s}"
            s = fmt.format(self.id, _isoformat_digits(self.time, time_digits),
                           SCNL([self.sta, self.chan, self.net,
                                self.loc]).__str__(),
                           self.error, phase_wd, assoc_wd)
        else: 
            fmt = "PickModified('{}', '{}', '{}', '{}', '{}',"
            fmt += " '{}', '{}', {}, '{}')"
            return fmt.format(self.sta, self.chan, self.net, self.loc,
                              _isoformat_digits(self.time, time_digits),
                              self.phase, self.error, self.assoc_id)
            s = self.__repr__()
        return s


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
    # P and S travel times are not completely necessary
    # because they can be calculated, but simpler to save
    tp = Column(DateTime)
    p_modified_id = Column(Integer)  # modified pick ID
    ts = Column(DateTime)
    s_modified_id = Column(Integer)  # modified pick ID
    locate_flag = Column(Boolean)
    assoc_id = Column(Integer)

    def __init__(self, ot, sta, d_km, delta, tp, p_modified_id,
                 ts, s_modified_id):
        """
        :param ot: origin time (DateTime?)
        :param sta: station name?
        :param d_km: distance in km station-origin??
        :param delta: distance in degrees?
        :param tp: P-pick time?
        :param p_modified_id: id of corresponding PickModified?
        :param ts: S-pick time?
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
        fmt = "Candidate('{}', '{}', {:g}, {:g}, '{}', {:d}, '{}', {:d})"
        return fmt.format(self.ot.isoformat(), self.sta, self.d_km, self.delta,
                          self.tp.isoformat(), self.p_modified_id,
                          self.ts.isoformat(), self.s_modified_id)

    def __str__(self, time_digits=2, table_format=False, table_header=False,
                include_create_update=False):
        """
        :param time_digits: number of digits after seconds in times
        :param table_format: output format for a table_format
        :param table_header: return table_format header (ignores table_format)
        :param include_create_update: include creation and update times
        """
        if table_header:
            fmt = " {:>5s} | {:^22s} | {:^9s} | {:^8s} | {:^8s} | {:^8s} | {:^8s} | {:^8s} "
            s = fmt.format('id', 'origin time (ot)', 'station', 'd_km', 'delta',
                                  'tp-ot', 'ts-tp', 'assoc_id')
        elif table_format:
            if self.assoc_id:
                assoc_wd = f"{self.assoc_id:d}"
            else:
                assoc_wd = '-'
            fmt = " {:5d} | {:22s} | {:^9s} | {:8.3f} | {:8.3f} | {:8.3f} | {:8.3f} | {:8s}"
            s = fmt.format(self.id, _isoformat_digits(self.ot, time_digits),
                          self.sta, self.d_km, self.delta,
                          (self.tp - self.ot).total_seconds(),
                          (self.ts - self.tp).total_seconds(),
                          assoc_wd)
        else: 
            fmt = "Candidate(ot={}, sta='{}', d_km={:.2f}, delta={:.2f},"\
                " tp=ot+{:." + str(time_digits) + "f}s,"\
                " ts=tp+{:." + str(time_digits) + "f}s)"
            return fmt.format(_isoformat_digits(self.ot, time_digits), self.sta,
                              self.d_km, self.delta,
                              (self.tp - self.ot).total_seconds(),
                              (self.ts - self.tp).total_seconds())
        return s


    def set_assoc_id(self, assoc_id, session, FT):
        self.assoc_id = assoc_id
        self.locate_flag = FT
        # Assign phases to modified picks

        # Actually only one pick_p and pick_s
        pick_p = session.query(PickModified).filter(
                    PickModified.id == self.p_modified_id)
        for pick in pick_p:
            pick.phase = 'P'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT

        pick_s = session.query(PickModified).filter(
                    PickModified.id == self.s_modified_id)
        for pick in pick_s:
            pick.phase = 'S'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT

        # Assign the phases to picks contribute to a modified picks
        picks_p = session.query(Pick).filter(
                    Pick.modified_id == self.p_modified_id).all()
        for pick in picks_p:
            pick.phase = 'P'
            pick.assoc_id = assoc_id
            pick.locate_flag = FT

        picks_s = session.query(Pick).filter(
                    Pick.modified_id == self.s_modified_id).all()
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
    loc_uncert = Column(Float)
    nsta = Column(Integer)
    t_create = Column(DateTime)
    t_update = Column(DateTime)

    def __init__(self, ot, ot_uncert, latitude, longitude, loc_uncert,
                 nsta, t_create, t_update):
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
        self.loc_uncert = loc_uncert
        self.nsta = nsta
        self.t_create = t_create
        self.t_update = t_update

    def __repr__(self):
        fmt = "Associated({}, {:.2f}, {:.3f}, {:.3f}, {:.3f}, {:d}, {}, {})"
        return fmt.format(self.ot.isoformat(),
                          self.ot_uncert, self.latitude, self.longitude,
                          self.loc_uncert, self.nsta,
                          self.t_create.isoformat(),
                          self.t_update.isoformat())


    def __str__(self, time_digits=2, table_format=False, table_header=False,
                include_create_update=False):
        """
        :param time_digits: number of digits after seconds in times
        :param table_format: output format for a table_format
        :param table_header: return table_format header (ignores table_format)
        :param include_create_update: include creation and update times
        """
        if table_header:
            fmt = " {:>5s} | {:^22s} | {:^9s} | {:^8s} | {:^8s} | {:^10s} | {:^4s}"
            s = fmt.format('id', 'orig_time', 'ot_uncert', 'lat',
                                  'lon', 'loc_uncert', 'nsta')
            if include_create_update:
                fmt = " | {:^15s} | {:^15s}"
                s += fmt.format('t_create', 't_update')
        elif table_format:
            fmt = " {:5d} | {:22s} | {:9.2f} | {:8.3f} | {:8.3f} | {:10.3f} | {:4d}"
            s = fmt.format(self.id, _isoformat_digits(self.ot, time_digits),
                          self.ot_uncert, self.latitude, self.longitude,
                          self.loc_uncert, self.nsta)
            if include_create_update:
                s += " | {:15s} | {:15s})".format(
                    _isoformat_digits(self.t_create, 0),
                    _isoformat_digits(self.t_update, 0))
        else: 
            fmt = "Associated({}, {:.2f}, {:.3f}, {:.3f}, {:.3f}, {:d}"
            s = fmt.format(_isoformat_digits(self.ot, time_digits),
                          self.ot_uncert, self.latitude, self.longitude,
                          self.loc_uncert, self.nsta)
            if include_create_update:
                s += ", {}, {})".format(
                    _isoformat_digits(self.t_create, 0),
                    _isoformat_digits(self.t_update, 0))
            else:
                s += ")"
        return s


def _isoformat_digits(time, digits):
    """
    return datetime as isoformat with specified digits after decimal
    
    :param digits: 0-6
    """
    s = time.strftime('%Y-%m-%dT%H:%M:%S')
    digits = int(digits)
    if digits <= 0:
        return s
    if digits > 6:
        digits = 6
    fmt='.{:0' + str(digits) + 'd}'
    s += fmt.format(int(time.microsecond * 10**(digits-6)))
    return s