"""
Find traveltime table rows closest to requested distance or S-P offset
"""
from .tt_stations_1D import TTtable1D


def tt_km(session, d_km):
    """
    Return the closest traveltime table row to the requested distance

    :param session: the database connection
    :param d_km: requested distance in km
    :returns: TTtable_object, km_difference the table row, the difference
              betwrrn that row's distance and the requested distance
    """

    min = session.query(TTtable1D).filter(TTtable1D.d_km <= d_km).\
        order_by(TTtable1D.d_km.desc()).first()
    max = session.query(TTtable1D).filter(TTtable1D.d_km >= d_km).\
        order_by(TTtable1D.d_km).first()
    if abs(min.d_km - d_km) <= abs(max.d_km - d_km):
        return min, abs(min.d_km - d_km)
    else:
        return max, abs(max.d_km - d_km)


def tt_s_p(session, s_p):
    """
    Return the closest traveltime table row to a requested S-P offset

    :param session: the database connection
    :param s_p: requested S-P arrival time offset
    :returns: TTtable_object, s_p_difference, the table row, the difference
              between that row's s_p offset and the requested value
    """
    min = session.query(TTtable1D).filter(TTtable1D.s_p <= s_p).\
        order_by(TTtable1D.s_p.desc()).first()
    max = session.query(TTtable1D).filter(TTtable1D.s_p >= s_p).\
        order_by(TTtable1D.s_p).first()
    if abs(min.s_p - s_p) <= abs(max.s_p - s_p):
        return min, abs(min.s_p - s_p)
    else:
        return max, abs(max.s_p - s_p)
