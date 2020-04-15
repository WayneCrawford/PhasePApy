"""
Create associations from Picks and Traveltime tables
"""
import sys
import logging
from datetime import datetime, timedelta
from operator import itemgetter
from itertools import combinations
from .tables1D import (Pick, PickModified, Candidate, Associated,
                       make_assoc_session)
from .tt_stations_1D import Station1D, TTtable1D, make_tt1D_session
# from .func1D import tt_s_p, tt_km
import numpy as np
from scipy.optimize import fmin
from obspy.geodetics import locations2degrees, gps2dist_azimuth
# Added for plotting
from obspy.core import Stream, read, UTCDateTime
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt

class LocalAssociator():
    """
    Associate picks with 1D velocity travel-times and fixed hypocenter depth.
    """
    def __init__(self, assoc_db_url, tt_db_url, max_km=350, aggregation=1,
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
        self.assoc_db = make_assoc_session(assoc_db_url)
        self.tt_stations_db_1D = make_tt1D_session(tt_db_url)

        self.max_km = max_km
        # Set max travel_time from max distance
        # tmp, d_diff = tt_km(self.tt_stations_db_1D, self.max_km)
        tmp = _tt_km(self.tt_stations_db_1D, self.max_km)
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
        s = 'LocalAssociator("{}", "{}", {:g}, {:g}, "{}", '.format(
            self.assoc_db.get_bind().url,
            self.tt_stations_db_1D.get_bind().url,
            self.max_km, self.aggregation, self.aggr_norm)
        s += '{:g}, {:d}, {:g}, {:g})'.format(
            self.assoc_ot_uncert, self.nsta_declare, self.cutoff_outlier,
            self.loc_uncert_thresh)
        return s

    def __str__(self):
        return self.str_associator_counts()

    def id_candidate_events(self, verbose=False):
        """
        Create candidate events based on S-P offsets.

        - Get all stations with unassociated picks
        - Condense picktimes within our pick uncertainty value
        - Generate all possible candidate events
        """
        # Get all stations with unassociated picks
        stations = [s.sta for s in self.get_stations('unassociated')]
        if verbose:
            print("id_candidate_events: "
                  f"{len(stations):d} stations with unassociated picks")

        counter = 0  # WCC added and got rid of counter0 below
        for sta in stations:
            picks = self.get_picks('unassociated', limit_list=[sta],
                                   order_by='time')
            # Condense picktimes that are within our pick uncertainty value
            temp, counter = _pick_cluster(self.assoc_db, picks,
                                          self.aggr_window, self.aggr_norm,
                                          counter)

            picks_modified = self.get_pick_modifieds('unassociated',
                                                     limit_list=[sta],
                                                     order_by='time')
            if verbose:
                print('\tFound {:d} picks for station {}'.format(
                    len(picks_modified), sta))
            # Generate all possible candidate events
            for i in range(0, len(picks_modified) - 1):
                for j in range(i + 1, len(picks_modified)):
                    s_p = (picks_modified[j].time
                           - picks_modified[i].time).total_seconds()
                    if s_p <= self.max_s_p and s_p >= self.min_s_p:
                        tt = _tt_s_p(self.tt_stations_db_1D, s_p)
                        if verbose:
                            print('\t\tFound an S-P candidate: {:g} seconds'.
                                  format(s_p))
                        ot = picks_modified[i].time -\
                            timedelta(seconds=tt.p_tt)
                        new_candidate = Candidate(ot, sta, tt.d_km, tt.delta,
                                                  picks_modified[i].time,
                                                  picks_modified[i].id,
                                                  picks_modified[j].time,
                                                  picks_modified[j].id)
                        self.assoc_db.add(new_candidate)
                        self.assoc_db.commit()
        # print('create candidate time in seconds: ',time.time()-now1, 's')

    def associate_candidates(self, verbose=False):
        """
        Associate candidate events

        A cluster is all events whose S-P-projected origin times fall within
            self.assoc_ot_uncert of each other
        Association is declared if there are > self.nsta_declare stations
            in a cluster
        """
        if verbose:
            print(f"{__name__}.assoc1D.associate_candidates()")
        # now2 = time.time()
        clusters = self.get_clusters(self.assoc_db, self.assoc_ot_uncert, True)
        dt_ot = timedelta(seconds=self.assoc_ot_uncert)
        clusters.sort(key= lambda i: i['n_stations'], reverse=True)

        if verbose:
            print(f"  {len(clusters)} origin time clusters found")
        i_cluster = 0
        for cluster in clusters:
            i_cluster += 1
            candidate_ot = cluster['origin_time']
            if cluster['n_stations'] >= self.nsta_declare:
                candis = self._get_goodtimes_nonassoc(candidate_ot, dt_ot)
                rms_sort = self._get_location_combinations(candis)

                if rms_sort:
                    rms_sort.sort(key=itemgetter(1))
                    loc, rms, matches = rms_sort[0]
                    lon, lat = loc[0], loc[1]
                    # remove outliers
                    # MATCHES_nol are matches with no outliers,
                    MATCHES_nol, mismatches = _outlier_cutoff(
                            matches, loc, self.cutoff_outlier)
                    MISMATCHES = []   # outliers, not for locating
                    while mismatches:
                        MISMATCHES.append(mismatches[0])
                        loc = fmin(_locating, [lon, lat], MATCHES_nol,
                                   disp=False)
                        MATCHES_nol, mismatches = _outlier_cutoff(
                                MATCHES_nol, loc, self.cutoff_outlier)
  
                    # declare event when RMS under control
                    nsta = len(MATCHES_nol)
                    if nsta >= self.nsta_declare:
                        LOC = fmin(_locating, (lon, lat), MATCHES_nol,
                                   disp=False)
                        LON, LAT = round(LOC[0], 3), round(LOC[1], 3)
                        OTS = [candis[c[5]].ot for c in MATCHES_nol]
                        origintime, ot_unc = _datetime_statistics(OTS)
                        RMS = _residuals_minimum(LOC, MATCHES_nol)
                        if RMS <= self.loc_uncert_thresh:
                            if verbose:
                                print(f"  Associated cluster {i_cluster}")
                            event_id = self._save_event(origintime, ot_unc,
                                                        LAT, LON, RMS, nsta)
                            # Associate candidates and picks with the event
                            for candi in MATCHES_nol:
                                candis[candi[5]].set_assoc_id(
                                    event_id, self.assoc_db, True)
                            self.assoc_db.commit()

                            # Associate candidates from outliers if d_km
                            # < loc_uncert
                            for MISMATCH in MISMATCHES:
                                d = gps2dist_azimuth(LAT, LON, MISMATCH[2],
                                                     MISMATCH[1])[0] / 1000
                                r = MISMATCH[3]
                                uncert_km = RMS * np.pi / 180.0 * 6371
                                if abs(d - r) <= uncert_km:
                                    candis[MISMATCH[5]].set_assoc_id(
                                        event_id, self.assoc_db, False)
                            self.assoc_db.commit()
            else:
                break
        if verbose:
            print("  Tested {} clusters with more than {} stations".format(
                i, self.nsta_declare))

    def _save_event(self, origintime, ot_unc, LAT, LON, RMS, nsta):
        """Save a new associated event to the database"""
        new_event = Associated(origintime, round(ot_unc, 3), LAT, LON,
                               round(RMS, 3), nsta, datetime.utcnow(),
                               datetime.utcnow())
        self.assoc_db.add(new_event)
        self.assoc_db.flush()
        self.assoc_db.refresh(new_event)
        self.assoc_db.commit()
        event_id = new_event.id
        logging.info('event_id: ' + str(event_id))
        logging.info(str(['ot:', origintime, 'ot_uncert:', ot_unc,
                          'loc:', LAT, LON, 'loc_uncert:', RMS,
                          'nsta:', nsta]))
        return event_id

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
                    # tt, d_diff = _tt_km(self.tt_stations_db_1D, d_km)
                    tt = _tt_km(self.tt_stations_db_1D, d_km)
                    p_predict = ot + timedelta(seconds=tt.p_tt)
                    s_predict = ot + timedelta(seconds=tt.s_tt)
                    half_win = timedelta(seconds=0.5 * self.aggr_window)

                    modi_pick_p = self._find_pick(sta, p_predict, half_win, 'P',
                                             event.id, verbose=verbose)
                    modi_pick_s = self._find_pick(sta, s_predict, half_win, 'S',
                                             event.id, verbose=verbose)
                    if verbose:
                        print()
                    if modi_pick_p or modi_pick_s:
                    #    self.assoc_db.add(candi)
                        event.nsta += 1
            self.assoc_db.commit()

    def count_associated(self):
        return self.assoc_db.query(Associated).count()

    def get_associated(self, assoc_id=None):
        """
        Get list of Associated events
        
        :param assoc_id: limit to this assoc_id
        """
        if assoc_id:
            return self.assoc_db.query(Associated).\
                filter(Associated.id==assoc_id).all()
        return self.assoc_db.query(Associated).all()

    def count_stations(self, type='all', limit_list=None, assoc_id=None):
        """
        Count the number of stations fitting the given criteria
        
        See get_stations for explanation of arguments
        """
        sel = self._select_stations(type, limit_list, assoc_id)
        return sel.count()
                       
    def get_stations(self, type='all', limit_list=None, assoc_id=None):
        """
        Return a list of Station1Ds fitting the given criteria
        
        :param type: type of stations to look for: "all", "associated",
                     "unassociated", "matched", "mismatched", "single-phase".
        :param limit_list: list of Candidates or Associateds to limit to
                           (must all be same type)
        :param assoc_id: assoc_id to limit to
        """
        sel = self._select_stations(type, limit_list, assoc_id)
        return sel.all()
                       
    def count_picks(self, type='all', limit_list=None, phase=None,
                    assoc_id=None):
        """
        Return number of Picks fitting given criteria
        
        See get_picks for explanation of arguments
        """
        sel = self._select_picks(type, limit_list, phase, assoc_id)
        #return eval(cmd + '.count()')
        return sel.count()

    def get_picks(self, type='all', limit_list=None, phase=None,
                    assoc_id=None, order_by=None):
        """
        Return list of Pick
        
        :param type: type of picks to look for: "all", "associated",
                     "unassociated","matched", "mismatched", "single-phase".
        :param limit_list: list of Candidates, Station1Ds, Associateds, or
                           station names to limit to (must all be same type)
        :param phase: limit to the given phase ('P' or 'S')
        :param assoc_id: assoc_id to limit to
        :param order_by: None or "time"
        """
        sel = self._select_picks(type, limit_list, phase, assoc_id, order_by)
        return sel.all()

    def count_pick_modifieds(self, type='all', limit_list=None, phase=None,
                    assoc_id=None):
        """
        Return number of PickModifieds fitting given criteria
        
        See get_picks for explanation of arguments
        """
        sel = self._select_pick_modifieds(type, limit_list, phase, assoc_id)
        #return eval(cmd + '.count()')
        return sel.count()
                       
    def get_pick_modifieds(self, type='all', limit_list=None, phase=None,
                    assoc_id=None, order_by=None):
        """
        Return list of PickModified
        
        See get_picks for explanation of arguments
        """
        sel = self._select_pick_modifieds(type, limit_list, phase, assoc_id)
        return sel.all()
                       
    def count_candidates(self, type='all', limit_list=None, assoc_id=None):
        """
        Return acces to selected stations from tt_stations database
        
        :param type: type of candidates to look for: "all", "associated",
                     "unassociated", "matched", "mismatched"
        :param limit_list: list of Stations, Associates or station names to
                           limit to (must all be same type)
        :param assoc_id: limit to given assoc_id
        """
        sel = self._select_candidates(type, limit_list, assoc_id)
        return sel.count()

    def get_candidates(self, type='all', limit_list=None, assoc_id=None):
        """
        Return acces to selected stations from tt_stations database
        
        :param type: type of candidates to look for: "all", "associated",
                     "unassociated", "matched", "mismatched"
        :param limit_list: list of Stations, Associates or station names to
                           limit to (must all be same type)
        :param assoc_id: limit to given assoc_id
        """
        sel = self._select_candidates(type, limit_list, assoc_id)
        return sel.all()

    def str_associator_counts(self):
        """
        Print a count of important Associator metrics
        """
        s = '{:2d} Picks, {:2d} PickModifieds, '.format(
            self.count_picks(), self.count_pick_modifieds())
        s += '{:2d} Candidates, {:2d} Associated'.format(
            self.count_candidates(), self.count_associated())
        if self.count_associated:
            s += ' ('
            first_time = True
            for assoc_id in [a.id for a in
                             self.assoc_db.query(Associated).all()]:
                pick_mods = self.assoc_db.query(PickModified).\
                    filter(PickModified.assoc_id == assoc_id)
                if first_time:
                    first_time = False
                else:
                    s += '; '
                s += 'nSta={:d},nP={:d},nS={:d}'.format(
                    len(set([p.sta for p in pick_mods.all()])),
                    pick_mods.filter(PickModified.phase == 'P').count(),
                    pick_mods.filter(PickModified.phase == 'S').count())
            s += ')'
        return s

    def str_associator_tables(self, order_by='time'):
        """
        Associator database tables string

        :param order_by: 'time' or 'id'
        """
        s = ''
        for type in ['Pick', 'PickModified', 'Candidate', 'Associated']:
            s+= self.set_associator_table(type, order_by) + '\n'
        return s

    def str_associator_table(self, type, order_by='time'):
        """
        Database table string

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
            return f'{__name__}.LocalAssociator.print_table(): ' +\
                   f'Invalid table type: "{type}"'
            return
        if order_by == 'id':
            order_basis = basis.id
        if self.assoc_db.query(basis).count():
            s = type + 's are:\n'
            first_time = True
            for obj in self.assoc_db.query(basis).\
                    order_by(order_basis):
                if first_time:
                    s += obj.__str__(table_header=True)
                    first_time = False
                s += '\n' + print(obj.__str__(table_format=True))
            return s
        else:
            return 'No ' + type + 's'

    def str_station_table(self, order_by='sta'):
        """
        Return database table

        :param order_by: 'sta' or 'id'
        """
        basis = Station1D
        order_basis = basis.sta
        if order_by == 'id':
            order_basis = basis.id
        if self.tt_stations_db_1D.query(basis).count():
            s = 'Stations are:\n'
            first_time = True
            for obj in self.tt_stations_db_1D.query(basis).\
                    order_by(order_basis):
                if first_time:
                    s += obj.__str__(table_header=True)
                    first_time = False
                s += '\n' + obj.__str__(table_format=True)
        else:
            return 'No Stations'


    def str_associated_candidate_stations(self):
        s=''
        for assoc in self.assoc_db.query(Associated).all():
            stas = self.assoc_db.query(Candidate.sta).\
                filter(Candidate.assoc_id == assoc.id).distinct().all()
            s += f'{assoc.id:2d}: {[s[0] for s in stas]}\n'
        return s


    def str_traveltime_table(self, order_by='d_km'):
        """
        Return database table

        :param order_by: 'd_km' or 'id'
        """
        basis = TTtable1D
        order_basis = basis.d_km
        if order_by == 'id':
            order_basis = basis.id
        if self.tt_stations_db_1D.query(basis).count():
            s = 'TravelTimes are:\n'
            first_time = True
            for obj in self.tt_stations_db_1D.query(basis).\
                    order_by(order_basis):
                if first_time:
                    s += obj.__str__(table_header=True) + '\n'
                    first_time = False
                s += obj.__str__(table_format=True) + '\n'
            return s
        else:
            return 'No TravelTimes'

    def str_associated(self, type, order_by='time'):
        """
        Print associated rows in an associator database table

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
            return(self.str_associator_table('Associated'))
        if order_by == 'id':
            order_basis = basis.id
        if order_by == 'station':
            order_basis = basis.sta
        if self.assoc_db.query(basis).filter(basis.assoc_id != None).count():
            s = 'Associated ' + type + 's are:\n'
            first_time = True
            for obj in self.assoc_db.query(basis).\
                    filter(basis.assoc_id != None).\
                    order_by(order_basis):
                if first_time:
                    s += obj.__str__(table_header=True) + '\n'
                    first_time = False
                s += obj.__str__(table_format=True) + '\n'
            return s
        else:
            return 'No Associated' + type + 's'

    def plot_section(self, assoc_id, files, seconds_before=5,
                     seconds_after=100, channel='Z', scale_km=2,
                     outfile=None):
        """
        Plot a record section, with picks and waveforms ordered by distance

        :param assoc_id: ID number for the association to plot
        :param files: list of waveform files
        :param seconds_before: seconds before origin_time to start plot
        :param seconds_after: seconds after origin_time to end plot
        :param channel: channel to plot (default='Z')
        :param scale_km: maximum wiggle will be this many km on plot
        :param outfile: save plot to this filename
        """
        assert files, 'empty wavefile list'
        # Read the waveforms
        ST = Stream()
        for file in files:
            ST += read(file)
        
        # Get selection criteria
        sta_list = [s.sta for s in self.get_stations('associated')]
        sta_list.extend([s.sta for s in self.get_stations('single-phase')])
        eve = self.get_associated(assoc_id=assoc_id)[0]
        ot_utc = UTCDateTime(eve.ot)    # obspy datetime format
 
        # Select and trim/clean traces
        ST_new = _select_traces(ST, sta_list, channel, ot_utc)
        ST_new.trim(ot_utc - seconds_before, ot_utc + seconds_after)
        ST_new.detrend('demean')
        # ST_new.filter('bandpass', freqmin=0.1, freqmax=100)

        # Get plot data
        sta_data = self._make_section_arrays(ST_new, scale_km, channel,
                                             seconds_before, assoc_id, eve)

        # Set up axis
        station_ticks = [x['tickloc'] for x in sta_data.values()]
        x_max = 1.1*max(station_ticks) - 0.1*min(station_ticks)
        fig = plt.figure(figsize=(15, 8))
        ax1 = fig.add_subplot(111)
        plt.ylim(-seconds_before, seconds_after)
        plt.xlim(0, x_max)

        # Plot data and picks
        for sta in sta_data.values():
            offset = sta['tickloc']
            # dots where picks cross the waveform
            for circle in sta['circles']:
                ax1.plot(circle['x'], circle['y'], 'o', c='gray')
            # data
            ax1.plot(sta['seg_data'][:,0] + offset, sta['seg_data'][:,1],
                     linewidth=.25, color='blue')
            # All modified picks
            for seg_pick in sta['seg_all_pickmodifieds']:
                ax1.plot(seg_pick[:,0]*.3 + offset, seg_pick[:,1],
                         linewidth=2, color='black') 
            # All initial picks
            for seg_pick in sta['seg_all_picks']:
                ax1.plot(seg_pick[:,0]*.3 + offset, seg_pick[:,1],
                         linewidth=1, color='black') 
            # Associated picks
            for seg_pick in sta['seg_assoc_picks']:
                ax1.plot(seg_pick[:,0] + offset, seg_pick[:,1],
                         linewidth=2, color='black') 
        
        # Set axis labels and ticks
        ax1.set_xticks(station_ticks)
        ax1.set_xticklabels([x for x in sta_data.keys()])
        ax1.invert_yaxis()
        ax1.xaxis.tick_top()
        plt.setp(plt.xticks()[1], rotation=45)
        ax2 = ax1.secondary_xaxis('bottom')
        ax2.set_xlabel('Offset(km)')
        plt.ylabel('Record Length (s)', fontsize=18)
        plt.title('{}-channel Section Plot of Event {:d} at {}'.format(
            channel, assoc_id, ot_utc))
        # plt.tight_layout()

        # Plot travel-time curves
        ttable = self.tt_stations_db_1D.query(TTtable1D).\
            filter(TTtable1D.d_km <= x_max).all()
        D = [x.d_km for x in ttable]
        ax1.plot(D, [x.p_tt for x in ttable], 'b--',
                 D, [x.s_tt for x in ttable], 'r--', linewidth=2)

        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()        

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
        candidate_ots = [c.ot for c in
                         candidate_base.order_by(Candidate.ot).all()]
        Array = []
        for candidate_ot in candidate_ots:
            cluster = candidate_base.\
                filter(Candidate.ot >= candidate_ot).\
                filter(Candidate.ot < (candidate_ot + dt_ot)).\
                order_by(Candidate.ot).all()
            cluster_stas = set([c.sta for c in cluster])
            # Array.append((candidate_ot, len(cluster_stas), len(cluster)))
            Array.append({'origin_time': candidate_ot,
                          'n_stations': len(cluster_stas),
                          'n_candidates': len(cluster)})
        return Array

    def _make_section_arrays(self, ST_new, scale_km, channel,
                             seconds_before, assoc_id, eve):
        """
        Make arrays used in section plots
        
        First step in cleaning up code
        """
        eq_lat, eq_lon = eve.latitude, eve.longitude
        sta_data={}
        for tr in ST_new:
            sta = tr.stats.station
            data = tr.data / (tr.data.max() - tr.data.min()) * scale_km
            
            # due to float point arithmetic issue, cannot use
            # "t=np.arange(0, tr.stats.npts / tr.stats.sampling_rate,
            #              tr.stats.delta)"
            t = np.arange(0, round(tr.stats.npts / tr.stats.sampling_rate
                                   / tr.stats.delta)) * tr.stats.delta
            t -= seconds_before
            seg=np.hstack((data[:, np.newaxis], t[:, np.newaxis]))
            lon, lat = self.tt_stations_db_1D.\
                query(Station1D.longitude, Station1D.latitude).\
                filter(Station1D.sta == tr.stats.station).first()
            distance = int(gps2dist_azimuth(lat, lon, eq_lat,
                                            eq_lon)[0] / 1000.)
            picks_p = None
            if channel == 'Z3':
                picks_p = self.assoc_db.query(Pick.time).\
                    filter(Pick.assoc_id == assoc_id).\
                    filter(Pick.sta == tr.stats.station).\
                    filter(Pick.chan == tr.stats.channel).\
                    filter(Pick.phase == 'P').all()
            if not picks_p:
                # picks_p = [x.time for x in self.get_pick_modifieds(
                #     assoc_id=assoc_id, phase='P', limit_list=[sta])]
                picks_p = self.assoc_db.query(PickModified.time).\
                    filter(PickModified.assoc_id == assoc_id).\
                    filter(PickModified.sta == tr.stats.station).\
                    filter(PickModified.phase == 'P').all()
            picks_s = self.assoc_db.query(PickModified.time).\
                filter(PickModified.assoc_id == assoc_id).\
                filter(PickModified.sta == tr.stats.station).\
                filter(PickModified.phase == 'S').all()
            picks_assoc = picks_p + picks_s
            picks_all = [x.time for x in self.get_picks('all', limit_list=[sta])]
            pickmodifieds_all = [x.time for x in self.get_pick_modifieds('all', limit_list=[sta])]
            seg_assoc_picks = []
            seg_all_picks = []
            seg_all_pickmodifieds = []
            circles=[]
            tickloc_picks2=[]
            for pick, in picks_assoc:
                index = int((pick-eve.ot + timedelta(seconds=seconds_before)).
                            total_seconds()/tr.stats.delta)
                if index >= 0 and index < len(data):
                    circles.append({'x': distance + data[index], 'y': t[index]})
                seg_assoc_picks.append(self._pick_seg(tr, t, data, pick, eve.ot, seconds_before))
            for pick in picks_all:
                pick_seg = self._pick_seg(tr, t, data, pick, eve.ot, seconds_before)
                if pick_seg is not None:
                    seg_all_picks.append(pick_seg)
            for pick in pickmodifieds_all:
                pick_seg = self._pick_seg(tr, t, data, pick, eve.ot, seconds_before)
                if pick_seg is not None:
                    seg_all_pickmodifieds.append(pick_seg)
            sta_data[sta] = {'seg_data': seg, 'tickloc': distance,
                             'circles': circles,
                             'seg_assoc_picks': seg_assoc_picks,
                             'seg_all_picks': seg_all_picks,
                             'seg_all_pickmodifieds': seg_all_pickmodifieds}
                                              

        return (sta_data)
    
    @staticmethod
    def _pick_seg(tr, t, data, picktime, orig_time, seconds_before):
        """ Make a bar segment at a given time"""
        index = int((picktime - orig_time + timedelta(seconds=seconds_before)).
                    total_seconds()/tr.stats.delta)
        if index < 0 or index >= len(data):
            return None
        # The following COULD BE ONE LINE?
        tp = np.array([t[index], t[index]])
        dp = np.array([data.min(), data.max()])
        return np.hstack((dp[:, np.newaxis], tp[:, np.newaxis]))
        
    def _get_goodtimes_nonassoc(self, target_ot, dt_ot):
        """Get Candidates with origin time matching target and without assocd picks"""
        candis = self.assoc_db.query(Candidate).\
            filter(Candidate.assoc_id == None).\
            filter(Candidate.ot >= target_ot).\
            filter(Candidate.ot < (target_ot + dt_ot)).\
            order_by(Candidate.ot).all()

        # remove the candidates with associated modified picks
        picks_associated_id = list(set(
            self.assoc_db.query(PickModified.id).
            filter(PickModified.assoc_id != None).all()))
        index_candis = []
        for id, in picks_associated_id:
            for j, candi in enumerate(candis):
                if (candi.p_modified_id == id) or\
                   (candi.s_modified_id == id):
                    index_candis.append(j)
        # delete from the end
        if index_candis:
            for j in sorted(set(index_candis), reverse=True):
                del candis[j]
        # print 'candis', candis
        return candis

    def _get_location_combinations(self, candis):
        """Get possible event locations and corresponding stations"""
        radius = []
        for i, candi in enumerate(candis):
            # pass in the radius for map plotting
            lon, lat = self.tt_stations_db_1D.\
                       query(Station1D.longitude, Station1D.latitude).\
                       filter(Station1D.sta == candi.sta).first()
            radius.append((candi.sta, lon, lat, candi.d_km,
                           candi.delta, i))
        cb = self._comb(radius)
        # print 'cb',cb

        rms_sort = []
        for i in range(len(cb)):
            radius_cb = cb[i]
            # self.nsta_declare has to be greater than or equal to 3
            if len(radius_cb) >= self.nsta_declare:
                # Set disp=True to print convergence messages.
                location = fmin(_locating, [lon, lat], radius_cb,
                                disp=False)
                residual_minimum = _residuals_minimum(location,
                                                     radius_cb)
                # rms_sort.append((location, residual_minimum, i))
                rms_sort.append((location, residual_minimum, cb[i]))
        return rms_sort
        #return cb, rms_sort

    def _find_pick(self, sta, t_predict, half_win, phase, event_id, verbose=False):
        """
        Find a pick matching the predicted time
        
        Associates the pickmodified and its parent picks()
        """
        picks_valid = self.assoc_db.query(PickModified).\
            filter(PickModified.sta == sta).\
            filter(PickModified.time >= (t_predict - half_win)).\
            filter(PickModified.time <= (t_predict + half_win)).\
            all()
        # if > 1 modified pick in range, associate the first one
        if picks_valid:
            if verbose:
                print(f'{phase} found...', end='')
            modi_pick = picks_valid[0]  # the first modified pick
            modi_pick.phase = phase
            modi_pick.assoc_id = event_id
            modi_pick.locate_flag = None
            # Associate all picks contributing to this
            # modified pick
            picks = self.assoc_db.query(Pick).\
                filter(Pick.modified_id == modi_pick.id).all()
            for pick in picks:
                pick.phase = phase
                pick.assoc_id = event_id
                pick.locate_flag = None
            return modi_pick
        return None

    def _comb(self, tt):
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
    def _verify_limit_list(limit_list):
        if limit_list:
            list_type = type(limit_list[0])
            assert np.all([type(l)==list_type for l in limit_list]),\
                f'list objects are not all type {list_type}'
            return list_type
        return None

    def _select_stations(self, type='all', limit_list=None, assoc_id=None):
        """
        Return a list of Station1Ds fitting the given criteria
        
        :param type: type of stations to look for:
        :param limit_list: list of elements to limit to
        :param assoc_id: assoc_id to limit to
        """
        if type == 'unassociated':
            picks = self.get_picks('associated', limit_list, assoc_id)
            assoc_sta_codes = list(set([p.sta for p in picks]))
            all_sta_codes = list(set([s.sta for s in self._select_stations('all')]))
            sta_codes = [c for c in all_sta_codes if c not in assoc_sta_codes]
        else:
            picks = self.get_picks(type, limit_list, assoc_id)       
            sta_codes = list(set([p.sta for p in picks]))
        sel = self.tt_stations_db_1D.query(Station1D).\
            filter(Station1D.sta.in_(sta_codes))
        return sel

    def _select_candidates(self, type='all', limit_list=None, assoc_id=None):
        """
        Return acces to selected stations from tt_stations database
        
        :param type: type of stations to look for
        :param limit_list: list of elements to limit to
        :param assoc_id: limit to given assoc_id
        """
        basis = 'Candidate'
        list_type = self._verify_limit_list(limit_list)

        # Build query + filter command
        cmd = f'self.assoc_db.query({basis})'
        if type == "associated":
            cmd += f'.filter({basis}.assoc_id != None)'
        if type == "unassociated":
            cmd += f'.filter({basis}.assoc_id != None)'
        elif type == "matched":
            cmd += f'.filter({basis}.assoc_id != None)'
            cmd += f'.filter({basis}.locate_flag == True)'
        elif type == "mismatched":
            cmd += f'.filter({basis}.assoc_id != None)'
            cmd += f'.filter({basis}.locate_flag == False)'

        if assoc_id:
            cmd += f'.filter({basis}.assoc_id == {assoc_id})'        

        if limit_list:
            if list_type is str:
                cmd += f'.filter({basis}.sta.in_(limit_list))'
            elif list_type is Associated:
                cmd += f'.filter({basis}.assoc_id.in_([l.id for l in limit_list])'
            elif list_type is Station1D:
                cmd += f'.filter({basis}.sta.in_([l.sta for l in limit_list])'
        
        return eval(cmd)
                       
    def _select_picks_common(self, basis, type='all', limit_list=None,
                                    phase=None, assoc_id=None, order_by=None):
        """
        Return acces to selected Picks or Pickmodifieds
        
        :param basis: "Pick" or "PickModified"
        :param type: type of stations to look for
        :param limit_list: list of elements to limit to
        :param assoc_id: limit to given assoc_id
        :param order_by: order list by: "time" or None
        """
        assert basis in ['Pick', 'PickModified'], f'invalid basis: {basis}'
        list_type = self._verify_limit_list(limit_list)

        # Build query + filter command
        cmd = f'self.assoc_db.query({basis})'
        if type == "associated":
            cmd += f'.filter({basis}.assoc_id != None)'
        if type == "unassociated":
            cmd += f'.filter({basis}.assoc_id == None)'
        elif type == "matched":
            cmd += f'.filter({basis}.assoc_id != None)'
            cmd += f'.filter({basis}.locate_flag == True)'
        elif type == "mismatched":
            cmd += f'.filter({basis}.assoc_id != None)'
            cmd += f'.filter({basis}.locate_flag == False)'
        elif type == "single-phase":
            # Associated PickModifieds that are not part of Candidates
            lt = self.assoc_db.query(
                Candidate.p_modified_id,Candidate.s_modified_id).\
                filter(Candidate.assoc_id != None).all()
            cand_pick_ids = [item for t in lt for item in t]
            cmd += f'.filter({basis}.assoc_id != None)'            
            cmd += f'.filter({basis}.id.notin_(cand_pick_ids))'            
            
        if phase:
            cmd += f'.filter({basis}.phase == "{phase}")'

        if assoc_id:
            cmd += f'.filter({basis}.assoc_id == {assoc_id})'        

        if limit_list:
            if list_type is str:
                cmd += f'.filter({basis}.sta.in_(limit_list))'
            elif list_type is Associated:
                cmd += f'.filter({basis}.assoc_id.in_([l.id for l in limit_list])'
            elif list_type is Station:
                cmd += f'.filter({basis}.sta.in_([l.sta for l in limit_list])'
            elif list_type is Candidate:
                cmd += f'.filter({basis}.id.in_([v for l in limit_list for m in [l.p_modfied_id, l.s_modfied_id] for v in m.values()]))'            
            else:
                f'Unsupported type "{list_type}" for limit_list'

        if order_by =='time':
            cmd += f'.order_by({basis}.time)'
        return eval(cmd)

    def _select_picks(self, type='all', limit_list=None, phase=None,
                    assoc_id=None, order_by=None, basis="Pick"):
        """
        Return acces to selected Picks
        
        :param type: type of stations to look for
        :param limit_list: list of elements to limit to
        :param assoc_id: limit to given assoc_id
        :param order_by: order list by: "time" or None
        """
        return self._select_picks_common("Pick", type, limit_list, phase,
                                    assoc_id, order_by)
                       
    def _select_pick_modifieds(self, type='all', limit_list=None, phase=None,
                    assoc_id=None, order_by=None):
        """
        Return acces to selected Pickmodifieds
        
        :param type: type of stations to look for
        :param limit_list: list of elements to limit to
        :param assoc_id: limit to given assoc_id
        :param order_by: order list by: "time" or None
        """
        return self._select_picks_common("PickModified", type, limit_list, phase,
                                    assoc_id, order_by)
                       

def _select_traces(ST, sta_list, channel, ot_utc):
    """
    Choose stream traces to plot
    
    :param ST: input Stream
    :param sta_list: list of stations to plot
    :param channel: channel code to select
    :param ot_utc: origin time of event
    :type ot_utc: UTCDateTime
    :returns: output Stream
    """
    # in case some stations use channel code like BH1, BH2
    # or BH3, make a universal search string:
    if channel == 'E' or channel == 'e':
        Chan = 'E1'
    elif channel == 'N' or channel == 'n':
        Chan = 'N2'
    elif channel == 'Z' or channel == 'z':
        Chan = 'Z3'
    else:
        print('Please input component E, e, N, n, Z, or z,'
              ' the default is Z')

    # Calculate distance from headers lat/lon
    ST_new = Stream()
    for tr in ST:
        stat = tr.stats
        if stat.channel[-1] in Chan and stat.station in sta_list:
            if ((stat.starttime < ot_utc) and (stat.endtime > ot_utc)):
                ST_new += tr

    # Remove duplicate traces (same station). VERY INEFFICIENT AND OVERLONG
    while True:
        stas = [tr.stats.station for tr in ST_new]
        dup = list(set([x for x in stas if stas.count(x) > 1]))
        if not dup:
            break
        index = [i for (i, j) in enumerate(stas)
                 if j == dup[-1]]
        i = 0
        while True:
            # Remove the shorter trace
            if (ST_new[index[i]].stats.npts
                    < ST_new[index[i+1]].stats.npts):
                del ST_new[index[i]]
                break
            elif (ST_new[index[i]].stats.npts
                  >= ST_new[index[i+1]].stats.npts):
                del ST_new[index[i+1]]
                break
    return ST_new


def _datetime_statistics(dt_list, norm='L2'):
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


def _pick_cluster(session, picks, pickwindow, pickaveraging_norm, counter):
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
        pickave, pickstd = _datetime_statistics(cluster_time,
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
            pickave, pickstd = _datetime_statistics(cluster_time,
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


def _locating(guess, *args):
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
        residuals += (_residual(guess, arg))**2
    return np.sqrt(residuals / len(args))


def _residuals_minimum(location, args):
    """
    Return the minimum residual.

    The only difference from locating() is that there is no * before args.
    """
    residuals = 0
    for arg in args:
        residuals += (_residual(location, arg))**2
    return np.sqrt(residuals / len(args))


def _residual(location, args):
    """
    Return the residual

    location is (lon_median, lat_median)
    args format is (sta, lon, lat, d_km, delta)

    uses obspy.geodetics.locations2degrees
    """
    offset_deg = locations2degrees(location[1], location[0],
                                   args[2], args[1])
    return offset_deg - args[4]


def _outlier_cutoff(matches, location, cutoff_outlier):
    """
    Cuts of outliers, or defines the outlier cutoff?
    """
    # 'tuple' object has no attribute 'remove', the matches passed in is
    # tuple, has to change to list
    matches = list(matches)

    res = []
    for n in range(len(matches)):
        x = _residual(location, matches[n])
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
    
def _tt_km(session, d_km):
    """
    Return a traveltime table row interpolated to the requested distance

    :param session: the database connection
    :param d_km: requested distance in km
    :returns: TTtable_object
    """

    closer = session.query(TTtable1D).filter(TTtable1D.d_km <= d_km).\
        order_by(TTtable1D.d_km.desc()).first()
    farther = session.query(TTtable1D).filter(TTtable1D.d_km >= d_km).\
        order_by(TTtable1D.d_km).first()
    if farther == closer:
        return closer
    alpha = (d_km - closer.d_km)/(farther.d_km - closer.d_km)
    interpd = TTtable1D(closer.d_km*(1 - alpha) + farther.d_km*alpha,
                        closer.delta*(1 - alpha) + farther.delta*alpha,
                        closer.p_tt*(1 - alpha) + farther.p_tt*alpha,
                        closer.s_tt*(1 - alpha) + farther.s_tt*alpha,        
                        closer.s_p*(1 - alpha) + farther.s_p*alpha)
    return interpd


def _tt_s_p(session, s_p):
    """
    Return a traveltime table row interpolated to the requested S-P offset

    :param session: the database connection
    :param s_p: requested S-P arrival time offset
    :returns: TTtable_object
    """

    closer = session.query(TTtable1D).filter(TTtable1D.s_p <= s_p).\
        order_by(TTtable1D.s_p.desc()).first()
    farther = session.query(TTtable1D).filter(TTtable1D.s_p >= s_p).\
        order_by(TTtable1D.s_p).first()
    if farther == closer:
        return closer
    alpha = (s_p - closer.s_p)/(farther.s_p - closer.s_p)
    interpd = TTtable1D(closer.d_km*(1 - alpha) + farther.d_km*alpha,
                        closer.delta*(1 - alpha) + farther.delta*alpha,
                        closer.p_tt*(1 - alpha) + farther.p_tt*alpha,
                        closer.s_tt*(1 - alpha) + farther.s_tt*alpha,        
                        closer.s_p*(1 - alpha) + farther.s_p*alpha)
    return interpd


