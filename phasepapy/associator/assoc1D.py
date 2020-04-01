"""
Create associations from Picks and Traveltime tables
"""
from .tables1D import Base, Pick, PickModified, Candidate, Associated
from .tt_stations_1D import BaseTT1D, Station1D, TTtable1D
from .func1D import tt_s_p, tt_km
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import numpy as np
from scipy.optimize import fmin
from obspy.geodetics import locations2degrees, gps2dist_azimuth
from datetime import datetime, timedelta
from operator import itemgetter
from itertools import combinations
import logging
import warnings


class LocalAssociator():
    """
    Associate picks with 1D velocity travel-times and fixed hypocenter depth.
    """
    def __init__(self, db_assoc, db_tt, max_km=350, aggregation=1,
                 aggr_norm='L2', assoc_ot_uncert=3, nsta_declare=3,
                 cutoff_outlier=30, loc_uncert_thresh=0.2):
        """
        :param db_assoc: associator database
        :param db_tt: travel time table database
        :param max_km: maximum distance of S-P interval in distance
        :param aggregation: the coefficient multiplied to minimum travel time
        :param aggr_norm: L2: median; L1: mean
        :param assoc_ot_uncert:  # of seconds between predicted origin times
                                 to associate candidate events
        :param nsta_declare: minimum station number to declare a earthquake
        :param cutoff_outlier: the outlier cut off distance in km
        :param loc_uncert_thresh: location uncertainty in degree
        """
        engine_associator = create_engine(db_assoc, echo=False)
        engine_tt_stations_1D = create_engine(db_tt, echo=False)
        # Base and BaseTT1D are imported from tables1D
        Base.metadata.create_all(engine_associator)
        BaseTT1D.metadata.create_all(engine_tt_stations_1D)
        Session1 = sessionmaker(bind=engine_associator)  # events table
        Session2 = sessionmaker(bind=engine_tt_stations_1D)  # traveltime table
        self.assoc_db = Session1()
        self.tt_stations_db_1D = Session2()

        self.max_km = max_km
        # Set max travel_time from max distance
        tmp, d_diff = tt_km(self.tt_stations_db_1D, self.max_km)
        self.max_tt = tmp.s_tt
        self.max_s_p = tmp.s_p
        self.min_s_p = self.tt_stations_db_1D.query(TTtable1D.s_p).\
            filter(TTtable1D.d_km == 0.0).first()[0]
        self.aggregation = aggregation
        self.aggr_window = self.aggregation * self.min_s_p
        self.aggr_norm = aggr_norm
        self.assoc_ot_uncert = assoc_ot_uncert
        self.nsta_declare = nsta_declare
        self.cutoff_outlier = cutoff_outlier
        self.loc_uncert_thresh = loc_uncert_thresh

    def __repr__(self):
        s = 'LocalAssociator({}, {}, {:g}, {:g}, "{}", '.format(
            self.assoc_db, self.tt_stations_db_1D, self.max_km,
            self.aggregation, self.aggr_norm)
        s += '{:g}, {:d}, {:g}, {:g})'.format(
            self.assoc_ot_uncert, self.nsta_declare, self.cutoff_outlier,
            self.loc_uncert_thresh)
        return s

    def num_associated(self):
        return self.assoc_db.query(Associated).count()

    def num_candidates(self):
        return self.assoc_db.query(Candidate).count()

    def num_picks(self):
        return self.assoc_db.query(Pick).count()

    def num_pickmodifieds(self):
        return self.assoc_db.query(PickModified).count()

    def print_counts(self):
        """
        Print the number of rows in the associator tables
        """
        print(f'{self.num_picks():2d} Picks, '
              f'{self.num_pickmodifieds():2d} PickModifieds, '
              f'{self.num_candidates():2d} Candidates, '
              f'{self.num_associated():2d} Associated')
              # f'{self.assoc_db.query(Associated).count():2d} Associated')

    def print_tables(self, order_by='time'):
        """
        Print all associator database tables
        
        :param order_by: 'time' or 'id'
        """
        for type in ['Pick', 'PickModified', 'Candidate', 'Associated']:
            self.print_table(type, order_by)

    def print_table(self, type, order_by='time'):
        """
        Print database table
        
        :param type: 'Pick', 'PickModified', 'Candidate' or 'Associated'
        :param order_by: 'time' or 'id'
        """
        if type == 'Pick':
            basis = Pick
            order_basis = basis.time
        elif type == 'PickModified':
            basis = PickModified
            order_basis = basis.time
        elif type == 'Candidate':
            basis = Candidate
            order_basis = basis.ot
        elif type == 'Associated':
            basis = Associated
            order_basis = basis.ot
        else:
            print(f'{__name__}.LocalAssociator.print_table(): Invalid table type: "{type}"')
            return
        if order_by == 'id':
            order_basis = basis.id
        if self.assoc_db.query(basis).count():
            print(type + 's are:')
            first_time = True
            for obj in self.assoc_db.query(basis).\
                    order_by(order_basis):
                if first_time:
                     print(obj.__str__(table_header=True))
                     first_time = False
                print(obj.__str__(table_format=True))
        else:
            print('No ' + type + 's')

    def print_associated(self, type, order_by='time'):
        """
        Print associated rows in a database table
        
        :param type: 'Pick', 'PickModified', or 'Candidate', will print
                     Associated otherwise
        :param order_by: 'time', 'id' or 'station'
        """
        if type == 'Pick':
            basis = Pick
            order_basis = basis.time
        elif type == 'PickModified':
            basis = PickModified
            order_basis = basis.time
        elif type == 'Candidate':
            basis = Candidate
            order_basis = basis.ot
        else:
            self.print_table('Associated')
            return
        if order_by == 'id':
            order_basis = basis.id
        if order_by == 'station':
            order_basis = basis.sta
        if self.assoc_db.query(basis).filter(basis.assoc_id != None).count():
            print('Associated ' + type + 's are:')
            first_time = True
            for obj in self.assoc_db.query(basis).\
                    filter(basis.assoc_id != None).\
                    order_by(order_basis):
                if first_time:
                     print(obj.__str__(table_header=True))
                     first_time = False
                print(obj.__str__(table_format=True))
        else:
            print('No Associated' + type + 's')

    def id_candidate_events(self, verbose=False):
        """
        Create candidate events based on possible S-P offsets.

        - Get all stations with unassociated picks
        - Condense picktimes within our pick uncertainty value
        - Generate all possible candidate events
        """
        # Get all stations with unassociated picks
        stations = self.assoc_db.query(Pick.sta).\
            filter(Pick.assoc_id == None).distinct().all()
        if verbose:
            print("id_candidate_events: "
                  f"{len(stations):d} stations with unassociated picks")

        counter = 0  # WCC added and got rid of counter0 below
        for sta, in stations:  # the comma is needed
            picks = self.assoc_db.query(Pick).\
                filter(Pick.sta == sta).\
                filter(Pick.assoc_id == None).order_by(Pick.time).all()
            # Condense picktimes that are within our pick uncertainty value
            # picktimes are python datetime objects
            # if stations.index((sta, )) == 0:  # stupid tuple
            #     counter0 = 0
            #     picktimes_new, counter = pick_cluster(
            #             self.assoc_db, picks, self.aggr_window,
            #             self.aggr_norm, counter0)
            # else:
            #     picktimes_new, counter = pick_cluster(
            #             self.assoc_db, picks, self.aggr_window,
            #             self.aggr_norm, counter)
            temp, counter = pick_cluster(self.assoc_db, picks,
                                         self.aggr_window,
                                         self.aggr_norm, counter)
                    
            picks_modified = self.assoc_db.query(PickModified).\
                filter(PickModified.sta == sta).\
                filter(PickModified.assoc_id == None).\
                order_by(PickModified.time).all()
            if verbose:
                print('\tFound {:d} picks for station {}'.format(
                    len(picks_modified), sta))
            # Generate all possible candidate events
            for i in range(0, len(picks_modified) - 1):
                for j in range(i + 1, len(picks_modified)):
                    s_p = (picks_modified[j].time -
                           picks_modified[i].time).total_seconds()
                    if s_p <= self.max_s_p and s_p >= self.min_s_p:
                        tt, tt_uncert = tt_s_p(self.tt_stations_db_1D, s_p)
                        if verbose:
                            print('\t\tFound an S-P candidate: {:g} seconds'.format(s_p))
                        ot = picks_modified[i].time -\
                            timedelta(seconds=tt.p_tt)
                        new_candidate = Candidate(ot, sta, tt.d_km, tt.delta,
                                                  picks_modified[i].time,
                                                  picks_modified[i].id,
                                                  picks_modified[j].time,
                                                  picks_modified[j].id)
                        self.assoc_db.add(new_candidate)
                        self.assoc_db.commit()
        # print 'create candidate time in seconds: ',time.time()-now1, 's'

    def associate_candidates(self):
        """
        Associate all possible candidate events

        by comparing the projected origin-times
        A cluster is all events whose origin times fall within
        self.assoc_ot_uncert of each other
        Association is declared if there are > self.nsta_declare stations in a cluster
        """
        # now2 = time.time()
        Array = self.get_clusters(self.assoc_db, self.assoc_ot_uncert, True)
        dt_ot = timedelta(seconds=self.assoc_ot_uncert)
        #
        # # Query all candidate origin times
        # candidate_ots = self.assoc_db.query(Candidate).\
        #     filter(Candidate.assoc_id == None).\
        #     order_by(Candidate.ot).all()
        # L_ots = len(candidate_ots)
        # # print L_ots
        # Array = []
        # for i in range(L_ots):
        #     cluster = self.assoc_db.query(Candidate).\
        #         filter(Candidate.assoc_id == None).\
        #         filter(Candidate.ot >= candidate_ots[i].ot).\
        #         filter(Candidate.ot < (candidate_ots[i].ot + dt_ot)).\
        #         order_by(Candidate.ot).all()
        #     cluster_sta = self.assoc_db.query(Candidate.sta).\
        #         filter(Candidate.assoc_id == None).\
        #         filter(Candidate.ot >= candidate_ots[i].ot).\
        #         filter(Candidate.ot < (candidate_ots[i].ot + dt_ot)).\
        #         order_by(Candidate.ot).all()
        #     l_cluster = len(set(cluster_sta))
        #     Array.append((i, l_cluster, len(cluster)))
        # print Array
        # sort Array by number of stations in cluster
        Array.sort(key=itemgetter(1), reverse=True)
        # print Array

        # print 'cluster analysis time:', time.time()-now2, 's'

        for i in range(len(Array)):
            candidate_ot = Array[i][0]
            # index = Array[i][0]
            if Array[i][1] >= self.nsta_declare:
                candis = self.assoc_db.query(Candidate).\
                    filter(Candidate.assoc_id == None).\
                    filter(Candidate.ot >= candidate_ot).\
                    filter(Candidate.ot < (candidate_ot + dt_ot)).\
                    order_by(Candidate.ot).all()
                # candis = self.assoc_db.query(Candidate).\
                #     filter(Candidate.assoc_id == None).\
                #     filter(Candidate.ot >= candidate_ots[index].ot).\
                #     filter(Candidate.ot < (candidate_ots[index].ot + dt_ot)).\
                #     order_by(Candidate.ot).all()

                # remove the candidates with associated modified picks 
                picks_associated_id = list(set(
                    self.assoc_db.query(PickModified.id).
                    filter(PickModified.assoc_id != None).all()))
                index_candis = []
                for id, in picks_associated_id:
                    for i, candi in enumerate(candis):
                        if (candi.p_modified_id == id) or\
                           (candi.s_modified_id == id):
                            index_candis.append(i)
                # delete from the end
                if index_candis:
                    for j in sorted(set(index_candis), reverse=True):
                        del candis[j]
                # print 'candis', candis

                # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                # 1D Associator
                # store all necessary parameter in lists
                radius = []
                for i, candi in enumerate(candis):
                    # pass in the radius for map plotting
                    lon, lat = self.tt_stations_db_1D.\
                               query(Station1D.longitude, Station1D.latitude).\
                               filter(Station1D.sta == candi.sta).first()
                    radius.append((candi.sta, lon, lat, candi.d_km,
                                   candi.delta, i))
                cb = self.comb(radius)
                # print 'cb',cb

                rms_sort = []
                for i in range(len(cb)):
                    radius_cb = cb[i]
                    # self.nsta_declare has to be greater than or equal to 3
                    if len(radius_cb) >= self.nsta_declare:
                        # Set disp=True to print convergence messages.
                        location = fmin(locating, [lon, lat], radius_cb,
                                        disp=False)
                        residual_minimum = residuals_minimum(location,
                                                             radius_cb)
                        rms_sort.append((location, residual_minimum, i))

                # It is possible to have empty rms_sort
                if rms_sort:
                    rms_sort.sort(key=itemgetter(1))
                    # loc is the location before outlier cutoff
                    loc, rms, index = rms_sort[0]
                    lon = loc[0]
                    lat = loc[1]
                    # matches is one of combination of
                    # radius.append([candi.sta, lon, lat, candi.d_km,
                    #                candi.delta, i])
                    matches = cb[index]
                    # print 'location: ', lat, lon, rms
                    # print 'matches',matches

                    # cut off outlier
                    MISMATCHES = []
                    # MATCHES_nol are the matches with no outliers,
                    # MISMATCHES is the outliers, which are not for locating
                    MATCHES_nol, mismatches = outlier_cutoff(
                            matches, loc, self.cutoff_outlier)
                    if mismatches:
                        MISMATCHES.append(mismatches[0])
                    while mismatches:
                        loc = fmin(locating, [lon, lat], MATCHES_nol,
                                   disp=False)
                        MATCHES_nol, mismatches = outlier_cutoff(
                                MATCHES_nol, loc, self.cutoff_outlier)
                        if mismatches:
                            MISMATCHES.append(mismatches[0])
                    # print(f"MATCHES_nol: {MATCHES_nol:d}, "
                    #       f"MISMATCHES: {MISMATCHES:d}")

                    # declare event when nsta and RMS are under control
                    nsta = len(MATCHES_nol)
                    if nsta >= self.nsta_declare:
                        LOC = fmin(locating, (lon, lat), MATCHES_nol,
                                   disp=False)
                        LON = round(LOC[0], 3)
                        LAT = round(LOC[1], 3)
                        OTS = []
                        for i in range(nsta):
                            OTS.append(candis[MATCHES_nol[i][5]].ot)
                        origintime, ot_unc = datetime_statistics(OTS)
                        RMS = residuals_minimum(LOC, MATCHES_nol)
                        t_create = datetime.utcnow()
                        t_update = datetime.utcnow()
                        if RMS <= self.loc_uncert_thresh:
                            new_event = Associated(origintime,
                                                   round(ot_unc, 3),
                                                   LAT, LON, round(RMS, 3),
                                                   nsta, t_create, t_update)
                            self.assoc_db.add(new_event)
                            self.assoc_db.flush()
                            self.assoc_db.refresh(new_event)
                            self.assoc_db.commit()
                            event_id = new_event.id

                            logging.info('event_id: ' + str(event_id))
                            logging.info(str(['ot:', origintime,
                                              'ot_uncert:', ot_unc,
                                              'loc:', LAT, LON,
                                              'loc_uncert:', RMS,
                                              'nsta:', nsta]))

                            # Associate candidates, picks with the identified
                            # event
                            for candi in MATCHES_nol:
                                candis[candi[5]].set_assoc_id(
                                        event_id, self.assoc_db, True)
                            self.assoc_db.commit()

                            # Associate candidates from outliers if the d_km
                            # intersect loc_uncert
                            if MISMATCHES:
                                for i in range(len(MISMATCHES)):
                                    d = gps2dist_azimuth(
                                        LAT, LON, MISMATCHES[i][2],
                                        MISMATCHES[i][1])[0] / 1000
                                    r = MISMATCHES[i][3]
                                    uncert_km = RMS * np.pi / 180.0 * 6371
                                    if abs(d - r) <= uncert_km:
                                        candis[MISMATCHES[i][5]].\
                                            set_assoc_id(event_id,
                                                         self.assoc_db, False)
                            self.assoc_db.commit()
                # 1D Associator
                # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            else:
                break

    def single_phase(self, verbose=False):
        """
        Associate stations with only S or P pick (not both)

        Does this require existing associated events and/or candidates
        based on S-P?
        """
        if verbose:
            print(f'In {__name__}.LocalAssociator.single_phase()')
        events = self.assoc_db.query(Associated).all()
        for event in events:
            event_id = event.id
            ot = event.ot
            if verbose:
                print(f'\tTesting event {event_id} at time {ot}')

            # Make a list of stations that are already associated
            sta_assoc = []
            for sta, in self.assoc_db.query(PickModified.sta).\
                    filter(PickModified.assoc_id == event_id).\
                    distinct().all():
                sta_assoc.append(sta)

            # associate single phase
            for sta, in self.assoc_db.query(PickModified.sta).\
                    filter(PickModified.assoc_id == None).\
                    filter(PickModified.time > ot).\
                    filter(PickModified.time <= (ot +
                           timedelta(seconds=self.max_tt))).\
                    distinct().all():
                station = self.tt_stations_db_1D.query(Station1D).\
                          filter(Station1D.sta == sta).first()
                d_km = gps2dist_azimuth(
                    event.latitude, event.longitude,
                    station.latitude, station.longitude)[0]/1000.

                # only associate single phase from stations without p-s pairs
                if (d_km < self.max_km) and (sta not in sta_assoc):
                    if verbose:
                        print(f'\t\tTesting Station {sta} at {d_km:.4g} km...',
                              end='')
                    tt, d_diff = tt_km(self.tt_stations_db_1D, d_km)

                    picks_p = self.assoc_db.query(PickModified).\
                        filter(PickModified.sta == sta).\
                        filter(PickModified.time >=
                               (ot + timedelta(seconds=tt.p_tt -
                                0.5 * self.aggr_window))).\
                        filter(PickModified.time <=
                               (ot + timedelta(seconds=tt.p_tt +
                                0.5 * self.aggr_window))).\
                        all()
                    # print('picks_p: ', picks_p, 'tt.p_tt: ', tt.p_tt)
                    # if there is more than one modified pick in the
                    # aggr_window range, only associate the first one
                    if picks_p:
                        if verbose:
                            print('P found...', end='')
                        modi_pick = picks_p[0]  # the first modified pick
                        modi_pick.phase = 'P'
                        modi_pick.assoc_id = event.id
                        modi_pick.locate_flag = False
                        # Associated all the picks contribute to this single
                        # modified picks with assoc_id and phase
                        picks = self.assoc_db.query(Pick).\
                            filter(Pick.modified_id == modi_pick.id).all()
                        for pick in picks:
                            pick.phase = 'P'
                            pick.assoc_id = event.id
                            pick.locate_flag = False

                    picks_s = self.assoc_db.query(PickModified).\
                        filter(PickModified.sta == sta).\
                        filter(PickModified.time >=
                               (ot + timedelta(seconds=tt.s_tt -
                                (0.5 * self.aggr_window)))).\
                        filter(PickModified.time <=
                               (ot + timedelta(seconds=tt.s_tt +
                                (0.5 * self.aggr_window)))).\
                        all()
                    # if there are more than one modified pick in the
                    # aggr_window range, only associate the first modified pick
                    if picks_s:
                        if verbose:
                            print('S found...', end='')
                        modi_pick = picks_s[0]  # the first modified pick
                        modi_pick.phase = 'S'
                        modi_pick.assoc_id = event.id
                        modi_pick.locate_flag = None
                        # Associated all the picks contribute to this single
                        # modified picks with assoc_id and phase
                        picks = self.assoc_db.query(Pick).\
                            filter(Pick.modified_id == modi_pick.id).all()
                        for pick in picks:
                            pick.phase = 'S'
                            pick.assoc_id = event.id
                            pick.locate_flag = None
                    if verbose:
                        print()
                    if picks_p or picks_s:
                        event.nsta += 1
            self.assoc_db.commit()

    def comb(self, tt):
        """
        Create the combinations from different stations
        """
        L = len(set([item[0] for item in tt]))   # length of the set(sta)
        # combinations of the array, some repeated, like (sta1, sta1, sta2,..)
        cb = list(combinations((tt), L))

        # remove those combinations of repeated station
        index = []
        for i in range(len(cb)):
            temp = []
            for j in range(L):
                temp.append(cb[i][j][0])
            lt = len(set(temp))
            if lt < L:
                index.append(i)
        index.reverse()
        for i in index:
            del cb[i]

        # only return combinations of different stations
        return cb

    @staticmethod
    def get_clusters(assoc_db, assoc_ot_uncert, only_unassociated):
        """
        Identify and return clusters among Candidates in the database
        
        :param assoc_db: associater database session
        :param assoc_ot_uncert: max difference in origin time within a cluster
        :param only_unassociated: only test unassociated Candidates
        :returns: Array of clusters: row = [origin_time, n_sta, n_picks]
        """
        dt_ot = timedelta(assoc_ot_uncert)
        candidate_base = assoc_db.query(Candidate)
        if only_unassociated:
            candidate_base = candidate_base.filter(Candidate.assoc_id == None)
        # Get ordered candidate origin times
        candidate_ots = [c.ot for c in candidate_base.order_by(Candidate.ot).all()]
        Array = []
        for candidate_ot in candidate_ots:
            cluster = candidate_base.\
                filter(Candidate.ot >= candidate_ot).\
                filter(Candidate.ot < (candidate_ot + dt_ot)).\
                order_by(Candidate.ot).all()
            cluster_stas = set([c.sta for c in cluster])
            Array.append((candidate_ot, len(cluster_stas), len(cluster)))
        return Array

def datetime_statistics(dt_list, norm='L2'):
    """
    Mean and standard deviations in seconds of a list of datetime values
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


def pick_cluster(session, picks, pickwindow, pickaveraging_norm, counter):
    """
    Cluster picks from different components on the same station.

    :param session:
    :param picks:
    :param pickwindow:
    :param pickaveraging_norm:
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

        picks[0].modified_id = 1+counter  # assign modified id to picks
        counter += 1
        pickave, pickstd = datetime_statistics(cluster_time,
                                               pickaveraging_norm)
        # append the row to the picks_new, not only the pick time
        picks_new.append(picks[0])
        pick_modified = PickModified(picks[0].sta, picks[0].chan,
                                     picks[0].net, picks[0].loc,
                                     picks[0].time, picks[0].phase,
                                     round(pickstd, 3), picks[0].assoc_id)
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
                if (picks[i + 1].chan not in channel) and\
                   ((picks[i + 1].time-picks[i].time).total_seconds() <
                   pickwindow):
                    cluster.append(picks[i + 1])
                    cluster_time.append(picks[i + 1].time)
                    channel.append(picks[i + 1].chan)
                    # assign modified id to picks
                    picks[i + 1].modified_id = counter
                    i += 1
                    # make sure do not go over the range limit because
                    # j=i+1 below, jump out inner while loop
                    if i == len(picks) - 1:
                        break
                # elif is dealing with the exactly same picks, probably from
                # processing same stream twice
                elif (picks[i + 1].sta == picks[i].sta) and\
                     (picks[i+1].chan == picks[i].chan) and\
                     (picks[i+1].time == picks[i].time):
                    # and picks[i+1].snr==picks[i].snr\
                    # and picks[i+1].phase==picks[i].phase\
                    # and picks[i+1].uncert==picks[i].uncert:
                    cluster.append(picks[i + 1])
                    cluster_time.append(picks[i + 1].time)
                    channel.append(picks[i + 1].chan)
                    # assign modified id to picks
                    picks[i + 1].modified_id = counter
                    i += 1
                    # make sure do not go over the range limit because j=i+1
                    # below, jump out inner while loop
                    if i == len(picks) - 1:
                        break
                else:
                    break
            pickave, pickstd = datetime_statistics(cluster_time,
                                                   pickaveraging_norm)

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
            # situation that last one is ungrouped, use if statement to add
            # in picks_new
            if j >= len(picks) - 1:
                if (picks[-1].time - picks[-2].time).total_seconds() >\
                        pickwindow:
                    picks_new.append(picks[-1])
                    # assign modified id to picks
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


def locating(guess, *args):
    """"
    Locating function

    Returns RMS of the distances between the guess and the args
    Sums the distance differences between the iterating guess distance and
    circle radius

    :guess: (lon, lat)
    :*args: list of (sta, lon, lat, d_km, delta) tuples
    """
    residuals = 0
    for arg in args:
        residuals += (residual(guess, arg))**2
    # # The following is the original version (does same thing)
    # L = len(args)
    # i = 0
    # while True:
    #     # gps2DistAzimuth(lat1, lon1, lat2, lon2)
    #     # Returns: (Great circle distance in m,
    #     #           azimuth A->B in degrees,
    #     #           azimuth B->A in degrees)
    #     residuals += (gps2dist_azimuth(
    #         guess[1], guess[0], args[i][2], args[i][1])[0] / 1000 *
    #         180 / (np.pi * 6371) - args[i][4])**2
    #     # np.sqrt((guess[0]-args[i][1])**2+(guess[1]-args[i][2])**2)
    #     #         - args[i][4])**2
    #     if i == L - 1:
    #         break
    #     else:
    #         i += 1
    return np.sqrt(residuals / len(args))


def residuals_minimum(location, args):
    """
    Return the minimum residual.

    The only difference from locating() is that there is no * before args.
    """
    residuals = 0
    for arg in args:
        residuals += (residual(location, arg))**2
    # # The following is the original version
    # L = len(args)
    # i = 0
    # while True:
    #     residuals += (gps2dist_azimuth(
    #         location[1], location[0], args[i][2], args[i][1])[0] / 1000 *
    #         180 / (np.pi * 6371) - args[i][4])**2
    #     if i == L - 1:
    #         break
    #     else:
    #         i += 1
    return np.sqrt(residuals / len(args))


def residual(location, args):
    """
    Return the residual

    location is (lon_median, lat_median)
    args format is (sta, lon, lat, d_km, delta)

    uses obspy.geodetics.locations2degrees
    """
    # x = gps2dist_azimuth(location[1], location[0], args[2], args[1])[0]\
    #     / 1000 * 180 / (np.pi * 6371) - args[4]
    # return x

    # offset_m = gps2dist_azimuth(location[1], location[0], args[2],
    #                             args[1])[0]
    # offset_deg = np.degrees(offset_m / 6371000)

    offset_deg = locations2degrees(location[1], location[0],
                                   args[2], args[1])
    return offset_deg - args[4]


def outlier_cutoff(matches, location, cutoff_outlier):
    """
    Cuts of outliers, or defines the outlier cutoff?
    """
    # 'tuple' object has no attribute 'remove', the matches passed in is
    # tuple, has to change to list
    matches = list(matches)

    res = []
    for n in range(len(matches)):
        x = residual(location, matches[n])
        res.append(x**2)  # (di - ri)**2
    m = max(res)
    index = [i for i, j in enumerate(res) if j == m]
    mismatch = []
    # if the min and max have the same absolute value, there will be two
    # indexes in index list
    for i in index:
        # print 'res',abs(6371*res[i]**0.5*np.pi/180.)
        if abs(6371 * res[i]**0.5 * np.pi / 180.) > cutoff_outlier:
            mismatch.append(matches[i])
            matches.remove(matches[i])

    # has to return tuple to locate
    return tuple(matches), tuple(mismatch)
