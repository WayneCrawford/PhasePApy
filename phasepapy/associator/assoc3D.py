"""
Create associations from Picks and Traveltime tables
"""
import numpy as np
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, PyramidSearching, SourceGrids

from .tables3D import Base, Pick, PickModified, Candidate, Associated
from .tt_stations_3D import BaseTT3D, TTtable3D
from .func3D import tt_km, tt_s_p
# from .search import *
from datetime import timedelta, datetime
from operator import itemgetter
from itertools import combinations
import logging
# import time       "from datetime import *" will import time


class LocalAssociator():
    """
    The 3D Associator associate picks with travel time curve of 3D velocity.
    """
    def __init__(self, db_assoc, db_tt, max_km=350, aggregation=1,
                 aggr_norm='L2', assoc_ot_uncert=3, nsta_declare=3,
                 nt=31, np=41, nr=5):
        """
        :param db_assoc: associator database
        :param db_tt: travel time table database
        :param max_km: maximum distance of S-P interval in distance
        :param aggregation: the coefficient multiplied to minimum travel time
        :param aggr_norm: L2: median; L1: mean
        :param assoc_ot_uncert: origin time uncertainty of candidates
        :param nsta_declare: minimum number of stations to declare an event
        :param nt, np, nr: node geometry
        """

        engine_associator = create_engine(db_assoc, echo=False)
        engine_tt_stations_3D = create_engine(db_tt, echo=False)
        Base.metadata.create_all(engine_associator)
        BaseTT3D.metadata.create_all(engine_tt_stations_3D)
        Session1 = sessionmaker(bind=engine_associator)  # events table
        Session2 = sessionmaker(bind=engine_tt_stations_3D)
        self.assoc_db = Session1()
        self.tt_stations_db_3D = Session2()
        self.max_km = max_km
        # Set our maximum travel_time from max distance
        tmp, d_diff = tt_km(self.tt_stations_db_3D, self.max_km)
        self.max_tt = tmp.s
        self.max_s_p = tmp.s_p
        self.min_s_p = self.tt_stations_db_3D.\
            query(TTtable3D.s_p).order_by(TTtable3D.s_p).first()[0]
        self.aggregation = aggregation
        self.aggr_window = self.aggregation * self.min_s_p
        self.aggr_norm = aggr_norm
        self.assoc_ot_uncert = assoc_ot_uncert
        self.nsta_declare = nsta_declare
        self.nt = nt
        self.np = np
        self.nr = nr

    def id_candidate_events(self):
        """
        Create a set of possible candidate events from our picks table.

        Where session is the connection to the sqlalchemy database.
        This method simply takes all picks with time differences less than
        our maximum S-P times for each station and generates a list of
        candidate events.
        """
        # now1 = time.time()
        #############
        # Get all stations with unnassoiated picks
        stations = self.assoc_db.query(Pick.sta).\
            filter(Pick.assoc_id == None).distinct().all()

        for sta, in stations:    # comma is needed because tupe
            picks = self.assoc_db.query(Pick).filter(Pick.sta == sta).\
                    filter(Pick.assoc_id == None).order_by(Pick.time).all()
            # Condense picktimes that are within our pick uncertainty value
            # picktimes are python datetime objects
            if stations.index((sta,)) == 0:  # stupid tuple
                counter0 = 0
                picktimes_new, counter = pick_cluster(self.assoc_db, picks,
                                                      self.aggr_window,
                                                      self.aggr_norm,
                                                      counter0)
            else:
                picktimes_new, counter = pick_cluster(self.assoc_db, picks,
                                                      self.aggr_window,
                                                      self.aggr_norm,
                                                      counter)
            picks_modified = self.assoc_db.query(PickModified).\
                filter(PickModified.sta == sta).\
                filter(PickModified.assoc_id == None).\
                order_by(PickModified.time).all()

            # Generate all possible candidate events
            for i in range(0, len(picks_modified) - 1):
                for j in range(i + 1, len(picks_modified)):
                    s_p = (picks_modified[j].time
                           - picks_modified[i].time).total_seconds()
                    # print(s_p)
                    if s_p <= self.max_s_p and s_p >= self.min_s_p:
                        nodes = tt_s_p(self.tt_stations_db_3D, s_p, .1)
                        ots = []
                        d_kms = []
                        deltas = []
                        for node in nodes:
                            ots.append(picks_modified[i].time
                                       - timedelta(seconds=node.p))
                            d_kms.append(node.d_km)
                            deltas.append(node.delta)
                        # L1: mean, L2: median
                        ot, ot_uncert = datetime_statistics(ots, norm='L2')
                        d_km = np.median(d_kms)
                        delta = np.median(deltas)
                        # print 'ot:',ot, d_km, delta
                        # print 'length of nodes:',len(nodes)

                        # ot=picks_modified[i].time-timedelta(seconds=tt.p_tt)
                        new_candidate = Candidate(ot, sta, d_km, delta,
                                                  picks_modified[i].time,
                                                  picks_modified[i].id,
                                                  picks_modified[j].time,
                                                  picks_modified[j].id)
                        self.assoc_db.add(new_candidate)
                        self.assoc_db.commit()
        # print 'id_candidate time in seconds: ',time.time()-now1

    def associate_candidates(self):
        """
        Associate all possible candidate events

        Associate by comparing projected origin-times. Does not yet deal with
        the condition that more picks and candidate events
        could arrive while we do our associations.
        """
        # now2 = time.time()
        dt_ot = timedelta(seconds=self.assoc_ot_uncert)

        # Query all candidate ots
        candidate_ots = self.assoc_db.query(Candidate).\
            filter(Candidate.assoc_id == None).\
            order_by(Candidate.ot).all()
        L_ots = len(candidate_ots)
        # print L_ots
        Array = []
        for i in range(L_ots):
            cluster = self.assoc_db.query(Candidate).\
                filter(Candidate.assoc_id == None).\
                filter(Candidate.ot >= candidate_ots[i].ot).\
                filter(Candidate.ot < (candidate_ots[i].ot + dt_ot)).\
                order_by(Candidate.ot).all()
            cluster_sta = self.assoc_db.query(Candidate.sta).\
                filter(Candidate.assoc_id == None).\
                filter(Candidate.ot >= candidate_ots[i].ot).\
                filter(Candidate.ot < (candidate_ots[i].ot + dt_ot)).\
                order_by(Candidate.ot).all()
            l_cluster = len(set(cluster_sta))
            Array.append((i, l_cluster, len(cluster)))
        # print Array
        # sort Array by l_cluster, notice Array has been changed
        Array.sort(key=itemgetter(1), reverse=True)
        # print Array
        # print 'candidates_ots:', time.time()-now2

        for i in range(len(Array)):
            index = Array[i][0]
            if Array[i][1] >= self.nsta_declare:
                matches = self.assoc_db.query(Candidate).\
                    filter(Candidate.assoc_id == None).\
                    filter(Candidate.ot >= candidate_ots[index].ot).\
                    filter(Candidate.ot < (candidate_ots[index].ot + dt_ot)).\
                    order_by(Candidate.ot).all()

                # remove candidates whose modified picks have been associated
                picks_associated_id = list(set(
                    self.assoc_db.query(PickModified.id).
                    filter(PickModified.assoc_id != None).all()))
                index_matches = []
                for id, in picks_associated_id:
                    for i, match in enumerate(matches):
                        if ((match.p_modified_id == id) or (match.s_modified_id
                                                            == id)):
                            index_matches.append(i)
                # delete from the end
                if index_matches:
                    for j in sorted(set(index_matches), reverse=True):
                        del matches[j]

                # 3D Associator
                # now = time.time()
                tt = []
                for match in matches:
                    # print 'sta:',match.sta
                    match_p = match.tp
                    match_s = match.ts
                    match_ot = match.ot
                    match_ttp = (match_p - match_ot).total_seconds()
                    match_tts = (match_s - match_ot).total_seconds()
                    tt.append([match.sta, match_ttp, match_tts, match.ot,
                               match])
                # print matches
                # print tt

                cb = self.comb(tt)
                # print 'cb',cb

                rms_sort = []
                for i in range(len(cb)):
                    tt_cb = cb[i]
                    # self.nsta_declare must be >= 3
                    if len(tt_cb) >= self.nsta_declare:
                        tt_new, sourcegrid, rms =\
                            PyramidSearching(self.tt_stations_db_3D, self.nt,
                                             self.np, self.nr, tt_cb)
                        rms_sort.append((tt_new, sourcegrid, rms, i))

                # rms_sort can be empty if all tt_cb have length < 3
                if rms_sort:
                    rms_sort.sort(key=itemgetter(2))
                    tt_new, sourcegrid, rms, index = rms_sort[0]
                    lat, lon, dep = self.tt_stations_db_3D.query(
                        SourceGrids.latitude,
                        SourceGrids.longitude,
                        SourceGrids.depth).\
                        filter(SourceGrids.id == sourcegrid).first()

                    nsta = len(tt_new)

                    all_ots = []
                    for i in range(nsta):
                        all_ots.append(tt_new[i][3])
                    origintime, ot_unc = datetime_statistics(all_ots)
                    # 3D Associator uses pick rms instead of loc_uncert
                    t_create = datetime.utcnow()
                    t_update = datetime.utcnow()
                    new_event = Associated(origintime, round(ot_unc, 3),
                                           lat, lon, dep, round(rms, 3),
                                           nsta, t_create, t_update)
                    self.assoc_db.add(new_event)
                    self.assoc_db.flush()
                    self.assoc_db.refresh(new_event)
                    self.assoc_db.commit()
                    event_id = new_event.id
                    logging.info(str(['event_id:', event_id]))
                    logging.info(str(['ot:', origintime,
                                      'ot_uncert:', round(ot_unc, 3),
                                      'loc:', lat, lon,
                                      'rms:', round(rms, 3),
                                      'nsta:', nsta]))

                    for tt_tuple in cb[index]:
                        match = tt_tuple[4]
                        match.set_assoc_id(event_id, self.assoc_db, True)
                    self.assoc_db.commit()
                # print 'time: ', time.time() - now
                # 3D Associator
            else:
                break

    def comb(self, tt):
        """
        Create the combinations from different stations

        :param tt: ???
        :returns: combinations of different stations
        """
        L = len(set([item[0] for item in tt]))   # length of the set(sta)
        # combinations of the array, some are repeated (sta1, sta1, sta2,...)
        cb = list(combinations((tt), L))

        # remove those combinations of repeated station
        index = []
        for i in range(len(cb)):
            temp = []
            for j in range(L):
                temp.append(cb[i][j][0])
            if len(set(temp)) < L:
                index.append(i)
        index.reverse()
        for i in index:
            del cb[i]

        # only return combinations of different stations
        return cb


def datetime_statistics(dt_list, norm='L2'):
    """
    Calculate the mean and standard deviations of a list of datetimes

    :param dt_list: list of datetimes
    :param norm: type of norm to use (L1 or L2)
    :returns mean, std:  mean and standard deviation (seconds)
    """
    offsets = []
    for dt in dt_list:
        offsets.append((dt - dt_list[0]).total_seconds())
    if norm == 'L1':
        mean_offsets = np.mean(offsets)
    elif norm == 'L2':
        mean_offsets = np.median(offsets)
    std_offsets = np.std(offsets)
    return dt_list[0] + timedelta(seconds=mean_offsets), std_offsets


def pick_cluster(session, picks, pickwindow, aggr_norm, counter):
    """
    Clean up very close picks on different channels of same station

    :param session:
    :param picks:
    :param pickwindow:
    :param aggr_norm:
    :param counter:
    """
    #               |    |                     /\
    #               |    |                    /  \          /\
    #               |    | /\      /\        /    \        /  \      /\
    #         ______|/\__|/  \    /  \      /      \      /    \    /  \  /\___
    #               |    |    \  /    \    /        \    /      \  /    \/
    #               |    |     \/      \  /          \  /        \/
    #               |    |              \/            \/
    # pickwindow:    ----
    # STA1 E     ---|----|--------------------|--------------
    # STA1 N     ----|-------------------------|-------------
    # STA1 Z     -----|-------------------------|------------
    # stack      ---|||--|--------------------|||------------
    # cluster STA1 --|---|---------------------|-------------
    # Chen highly recommends using norm=='L2' to lower the effect of outliers
    # Better to set pickwindow==t_up, t_up is to clean closed picks
    # ARGUE: whether to only take the median or mean of the picks from
    #        different stations? won't count the followings after first one

    picks_new = []
    # only one pick in picks
    if len(picks) == 1:
        cluster = []
        cluster.append(picks[0])
        cluster_time = []
        cluster_time.append(picks[0].time)
        picks[0].modified_id = 1 + counter  # assign modified id to picks
        counter += 1
        pickave, pickstd = datetime_statistics(cluster_time, aggr_norm)
        # append the row to the picks_new, not only the pick time
        picks_new.append(picks[0])
        pick_modified = PickModified(picks[0].sta, picks[0].chan,
                                     picks[0].net, picks[0].loc, picks[0].time,
                                     picks[0].phase, round(pickstd, 3),
                                     picks[0].assoc_id)
        session.add(pick_modified)
        session.commit()
    # more than one pick in picks
    else:
        j = 0
        counter += 1
        while True:
            i = j
            cluster = []
            cluster.append(picks[i])
            cluster_time = []
            cluster_time.append(picks[i].time)
            channel = []
            channel.append(picks[i].chan)
            picks[i].modified_id = counter
            while True:
                # cluster picks of different channels; notice that for the
                # closed picks on the same station, those picks behind the
                # first pick could be separated lonely or separated cluster
                if ((picks[i+1].chan not in channel)
                    and ((picks[i+1].time-picks[i].time).total_seconds()
                         < pickwindow)):
                    cluster.append(picks[i+1])
                    cluster_time.append(picks[i+1].time)
                    channel.append(picks[i+1].chan)
                    picks[i+1].modified_id = counter  # assign modified id
                    i += 1
                    # make sure do not go over the range limit because j=i+1
                    # below, jump out inner while loop
                    if i == len(picks) - 1:
                        break
                # elif is dealing with the exactly same picks, probably from
                # processing same stream twice
                elif (picks[i+1].sta == picks[i].sta
                      and picks[i+1].chan == picks[i].chan
                      and picks[i+1].time == picks[i].time):
                    # and picks[i+1].snr == picks[i].snr
                    # and picks[i+1].phase == picks[i].phase
                    # and picks[i+1].uncert == picks[i].uncert:
                    cluster.append(picks[i+1])
                    cluster_time.append(picks[i+1].time)
                    channel.append(picks[i+1].chan)
                    picks[i+1].modified_id = counter  # assign modified id
                    i += 1
                    # make sure do not go over the range limit because j=i+1
                    # below, jump out inner while loop
                    if i == len(picks) - 1:
                        break
                else:
                    break
            pickave, pickstd = datetime_statistics(cluster_time, aggr_norm)

            # append whole rows to the picks_new, not only the pick time
            for pick in cluster:
                if (pick.time - pickave).total_seconds() >= 0:
                    break
            picks_new.append(pick)
            pick_modified = PickModified(pick.sta, pick.chan, pick.net,
                                         pick.loc, pick.time, pick.phase,
                                         round(pickstd, 3), pick.assoc_id)
            session.add(pick_modified)
            session.commit()
            # next cluster
            j = i + 1
            counter += 1

            # jump outer while loop and compare last two picks. For the
            # situation that last one == ungrouped, use if statement to add
            # in picks_new
            if j >= len(picks) - 1:
                if ((picks[-1].time
                     - picks[-2].time).total_seconds() > pickwindow):
                    picks_new.append(picks[-1])
                    picks[-1].modified_id = counter
                    pick_modified = PickModified(picks[-1].sta, picks[-1].chan,
                                                 picks[-1].net, picks[-1].loc,
                                                 picks[-1].time,
                                                 picks[-1].phase,
                                                 round(pickstd, 3),
                                                 picks[-1].assoc_id)
                    session.add(pick_modified)
                    session.commit()
                else:
                    if picks[-1] in cluster:
                        counter -= 1
                    else:
                        picks[-1].modified_id = counter
                        pick_modified = PickModified(picks[-1].sta,
                                                     picks[-1].chan,
                                                     picks[-1].net,
                                                     picks[-1].loc,
                                                     picks[-1].time,
                                                     picks[-1].phase,
                                                     round(pickstd, 3),
                                                     picks[-1].assoc_id)
                        session.add(pick_modified)
                        session.commit()
                break

    return picks_new, counter
