"""
Find traveltime table rows closest to requested distance or S-P offset
"""
from .tt_stations_3D import TTtable3D


def tt_km(session, d_km):
    """
    Return the closest traveltime table row to the requested distance

    :param session: the database connection
    :param d_km: requested distance in km
    :returns: TTtable_object, km_difference the table row, the difference
              between that row's distance and the requested distance
    """

    min = session.query(TTtable3D).filter(TTtable3D.d_km <= d_km).\
        order_by(TTtable3D.d_km.desc()).first()
    max = session.query(TTtable3D).filter(TTtable3D.d_km >= d_km).\
        order_by(TTtable3D.d_km).first()
    if abs(min.d_km - d_km) <= abs(max.d_km - d_km):
        return min, abs(min.d_km - d_km)
    else:
        return max, abs(max.d_km - d_km)


def tt_s_p(session, s_p, uncert):
    """
    Return the closest traveltime table row to a requested S-P offset

    :param session: the database connection
    :param s_p: requested S-P arrival time offset
    :param uncert: allowed uncertainty?
    :returns TTtable_object: the traveltime table row
    """

    nodes = session.query(TTtable3D).filter(TTtable3D.s_p <= s_p + uncert).\
        filter(TTtable3D.s_p >= s_p - uncert).all()
    return nodes
