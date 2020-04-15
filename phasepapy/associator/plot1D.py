"""
Plot class for 1D tables
"""
from datetime import timedelta
import warnings

from obspy.core import UTCDateTime, Stream, read
from obspy.geodetics import gps2dist_azimuth
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
# from matplotlib.colors import colorConverter
from mpl_toolkits.basemap import Basemap
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

from .tables1D import Pick, PickModified, Candidate, Associated
from .tt_stations_1D import TTtable1D, Station1D
from .assoc1D import LocalAssociator
from .misc import isoformat_digits


class Plot():
    def __init__(self, db_assoc, db_tt):
        """ Define the travel time and associator database """
        engine_assoc = create_engine(db_assoc, echo=False)
        # create a configuration file including paths
        engine_tt_stations = create_engine(db_tt, echo=False)
        Session1 = sessionmaker(bind=engine_assoc)  # events table
        Session2 = sessionmaker(bind=engine_tt_stations)  # traveltime table
        self.assoc_db = Session1()
        self.tt_stations_db_1D = Session2()

    def cluster_plot(self, assoc_ot_uncert=3, outfile=None):
        """
        Plot clustering nodes as a function of event time
        
        :param assoc_ot_uncert: the maximum allowed offset between origin
                                times within a cluster
        """
        matplotlib.rcParams["axes.labelsize"] = "large"
        matplotlib.rcParams["axes.linewidth"] = 2.0
        matplotlib.rcParams["xtick.major.size"] = 8
        matplotlib.rcParams["ytick.major.size"] = 8
        matplotlib.rcParams["ytick.minor.size"] = 5
        matplotlib.rcParams["xtick.labelsize"] = "large"
        matplotlib.rcParams["ytick.labelsize"] = "large"

        clusters = LocalAssociator.get_clusters(self.assoc_db, assoc_ot_uncert,
                                             False)

        x2 = np.array([x['origin_time'] for x in clusters])
        y1 = np.array([x['n_stations'] for x in clusters])
        y2 = np.array([x['n_candidates'] for x in clusters])
        plt.figure(figsize=(15, 6))
        # p1=plt.plot(x2,y1,'ro-')
        p1 = plt.plot(x2, y1, 'o-', c='gray')
        plt.ylim(0, max(y1) + 1)
        plt.xlim(min(x2) - timedelta(seconds=60),
                 max(x2) + timedelta(seconds=60))
        plt.xlabel('Time (s)', fontsize=20)
        plt.ylabel('Count', fontsize=20)
        plt.legend(p1, ['count'])
        plt.tight_layout()
        plt.title('Unique Candidates Cluster Analysis', fontsize=20)
        plt.figure(figsize=(15, 6))
        # p2=plt.plot(x2,y2,'bo-')
        p2 = plt.plot(x2, y2, 'o-', c='gray')
        plt.ylim(0, max(y2) + 1)
        plt.xlim(min(x2) - timedelta(seconds=60),
                 max(x2) + timedelta(seconds=60))
        plt.xlabel('Time (s)', fontsize=20)
        plt.legend(p2, ['count'])
        plt.ylabel('Count', fontsize=20)
        plt.tight_layout()
        plt.title('Total Candidates Cluster Analysis', fontsize=20)
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()

    def event_plot(self, assoc_id, west=-104.5, east=-94, south=33.5,
                   north=37.5, deltalon=None, deltalat=None, outfile=None,
                   plot_residuals=False, max_dist=None, max_s_p=None,
                   plot_legend=True, plot_create_update_times=False):
        """
        Plot the circles, stations, location and residuals on one map

        :param assoc_id:  event number after the event been associated
        :param west: west bound of map, None autocalculates
        :param east: east bound of map, None autocalculates
        :param south: south bound of map, None autocalculates
        :param north: north bound of map, None autocalculates
        :param deltalon: logitude tick spacing, None autocalculates
        :param deltalat: latitude tick spacing, None autocalculates
        :param plot_residuals:  plot residuals on the map
        :param max_dist:  maximum distance on residuals plot
        :param max_s_p:  maximum S-P time on residuals plot
        :param plot_legend:  plot a map legend
        :param plot_legend:  plot event associator create and update times
        """
        # Set up plot
        from itertools import cycle
        plot_colors = cycle(['r', 'g', 'c', 'm', 'y', 'k', 'b'])
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)        
        # Draw basemap and boundaries
        bounds = self._get_bounds(assoc_id, west, east, south, north, deltalon, deltalat)
        m = Basemap(projection='merc', llcrnrlat=bounds['south'], urcrnrlat=bounds['north'],
                    llcrnrlon=bounds['west'], urcrnrlon=bounds['east'], lat_ts=0, resolution='c')
        m.drawcountries()
        m.drawstates()
        m.fillcontinents(color='white', lake_color='aqua', zorder=1)
        # draw parallels, meridians and structures.
        m.drawparallels(np.arange(bounds['south'], bounds['north'],
                                  bounds['deltalat']),
                        labels=[1, 0, 0, 0], color='gray', dashes=[1, 1e-5],
                        labelstyle='+/-', linewidth=0.1)
        m.drawmeridians(np.arange(bounds['west'], bounds['east'],
                                  bounds['deltalon']),
                        labels=[0, 0, 0, 1], color='gray', dashes=[1, 1e-5],
                        labelstyle='+/-', linewidth=0.1)
        m.drawmapboundary(fill_color='blue')

        # Get lists of Candidates and single-phase picks
        matches = self.assoc_db.query(Candidate).\
            filter(Candidate.assoc_id == assoc_id).\
            filter(Candidate.locate_flag == True).all()
        mismatches = self.assoc_db.query(Candidate).\
            filter(Candidate.assoc_id == assoc_id).\
            filter(Candidate.locate_flag == False).all()        
        single_phases = self.assoc_db.query(PickModified).\
            filter(PickModified.assoc_id == assoc_id).\
            filter(PickModified.locate_flag == None).all()

        # Get lists of stations
        match_stations = [self.tt_stations_db_1D.query(Station1D).
                          filter(Station1D.sta == x.sta).first()
                          for x in matches]
        mismatch_stations = [self.tt_stations_db_1D.query(Station1D).
                             filter(Station1D.sta == x.sta).first()
                             for x in mismatches]
        single_phase_stations = [self.tt_stations_db_1D.query(Station1D).
                                 filter(Station1D.sta == x.sta).first()
                                 for x in single_phases]
        
        # Make S-P lists        
        s_p_list=[{'time': (c.ts - c.tp).total_seconds(),
                   'latitude': s.latitude,
                   'longitude': s.longitude,
                   'radius': c.d_km,
                   'color': plot_colors.__next__(),
                   'station': c.sta} for c, s in zip(matches, match_stations)]
        s_p_list2=[{'time': (c.ts - c.tp).total_seconds(),
                    'latitude': s.latitude,
                    'longitude': s.longitude,
                    'radius': c.d_km,
                    'color': plot_colors.__next__(),
                    'station': c.sta} for c, s in zip(mismatches, mismatch_stations)]
        
        # Plot match S-P distance circles
        for s_p in s_p_list:
            line = _draw_circle(m, s_p['longitude'], s_p['latitude'],
                                s_p['radius'], lw=2., color=s_p['color'])
        legend_list = {"S-P Pair for locating": line}

        # Plot mismatch S-P distance circles
        for s_p in s_p_list2:
            line = _draw_circle(m, s_p['longitude'], s_p['latitude'],
                                s_p['radius'], lw=2., color=s_p['color'])
        if mismatches:
            legend_list["S-P Pair not for locating"] = line
            
        s_p_list.extend(s_p_list2)

        # Plot all S-P stations
        xs, ys = m([s['longitude'] for s in s_p_list],
                   [s['latitude'] for s in s_p_list])
        m.scatter(xs, ys, marker='^', s=100, c=[s['color'] for s in s_p_list], zorder=3)
        for x, y, sta in zip(xs, ys, [s['station'] for s in s_p_list]):
            plt.text(x, y + 1, sta, fontweight='bold', ha='center',
                     color='k')

        # plot event location and uncertainty
        event = self.assoc_db.query(Associated).\
            filter(Associated.id == assoc_id).first()
        uncert = event.loc_uncert * 6371 * np.pi / 180.
        x, y = m(event.longitude, event.latitude)
        m.scatter(x, y, marker='*', s=100, c='k', zorder=3)
        _draw_circle(m, event.longitude, event.latitude, uncert, color='k', lw=2.)
        plt.title('Event {:d} at Origin Time: {}'.format(event.id,
            isoformat_digits(event.ot,2)), fontsize=18)

        if plot_create_update_times:
            t_create = event.t_create
            t_update = event.t_update
            x_text, y_text = m(bounds['west']
                               + 0.2 * (bounds['east'] - bounds['west']),
                               bounds['south']
                               + 0.05 * (bounds['north'] - bounds['south']))
            plt.text(x_text, y_text,
                     f'Create Time: {t_create}\nUpdate Time:{t_update}',
                     horizontalalignment='center', fontsize=12,
                     fontweight='bold')

        # plot single phase stations and their circles
        for single_phase in single_phases:
            lon, lat = self.tt_stations_db_1D.\
                query(Station1D.longitude, Station1D.latitude).\
                filter(Station1D.sta == single_phase.sta).first()
            x, y = m(lon, lat)
            pick_time = single_phase.time
            travel_time = (pick_time - event.ot).total_seconds()
            # ANY DIFFERENCE BETWEEN THE NEXT TWO BESIDES ZORDER?
            zorder=4
            if single_phase.phase == 'P':
                zorder=5
             # tt, tt_uncert = self._distance_singlephase(single_phase.phase,
            #                                           travel_time)
            tt = self._distance_singlephase(single_phase.phase, travel_time)
            line = _draw_circle(m, lon, lat, tt.d_km, ls=':', lw=2.,
                                c='gray', zorder=4)
            m.scatter(x, y, marker='^', s=100, c='gray', zorder=zorder)
            plt.text(x, y+1, single_phase.sta, fontweight='bold',
                     ha='center')
#             if single_phase.phase == 'S':
#                 zorder=4
#                 # tt, tt_uncert = self._distance_singlephase(single_phase.phase,
#                 #                                           travel_time)
#                 tt = self._distance_singlephase(single_phase.phase,
#                                                 travel_time)
#                 d_km = tt.d_km
#                 line = _draw_circle(m, lon, lat, d_km, ls=':', lw=2.,
#                                    c='gray', zorder=4)
#                 m.scatter(x, y, marker='^', s=100, c='gray', zorder=zorder)
#                 plt.text(x, y+1, single_phase.sta, fontweight='bold',
#                          ha='center')
        if single_phases:
            legend_list["Single Phase"] = line

        if plot_legend:
            legend = plt.legend(legend_list.values(), legend_list.keys(),
                                loc='upper left')

        if plot_residuals:
            subpos = [0.05, 0.28, 0.35, 0.6]
            subax = add_subplot_axes(ax, subpos)
            self._plot_residuals(subax, s_p_list, max_dist, max_s_p)
        
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
            
    def section_plot(self, assoc_id, files, seconds_ahead=5,
                     record_length=100, channel='Z', scale_factor=2,
                     outfile=None):
        """
        Plot a record section, with picks and waveforms ordered by distance

        :param assoc_id: ID number for the association to plot
        :param files: list of waveform files
        :param seconds_ahead: time before origin_time to start plot
        :param record_length: time after origin_time to end plot
        :param channel: channel to plot (default='Z')
        :param scale_factor: factor to scale wiggles by
        :param outfile: save plot to this filename
        """
        assert files, 'empty wavefile list'
        station = self.assoc_db.query(Candidate.sta).\
            filter(Candidate.assoc_id == assoc_id).all()
        sta_list = []
        for sta, in station:
            sta_list.append(str(sta))
        station_single = self.assoc_db.query(Pick.sta).\
            filter(Pick.assoc_id == assoc_id).\
            filter(Pick.locate_flag == None).all()
        for sta, in station_single:
            sta_list.append(str(sta))
        # print sta_list

        eve = self.assoc_db.query(Associated).\
            filter(Associated.id == assoc_id).first()
        # Earthquakes' epicenter
        eq_lat = eve.latitude
        eq_lon = eve.longitude

        # Read the waveforms
        ST = Stream()
        for file in files:
            st = read(file)
            ST += st

        # in case of some seismometer use channel code like BH1, BH2
        # or BH3, resign the channel code as:
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
            if tr.stats.channel[2] in Chan and tr.stats.station in sta_list:
                if ((tr.stats.starttime.datetime < eve.ot)
                        and (tr.stats.endtime.datetime > eve.ot)):
                    tr.trim(UTCDateTime(eve.ot
                                        - timedelta(seconds=seconds_ahead)),
                            UTCDateTime(eve.ot
                                        + timedelta(seconds=record_length)))
                    ST_new += tr
        # print ST_new.__str__(extended=True)

        # remove traces from same station, but different # samples
        # X31A..B Z | 2011-01-18T01:40:30.325Z - 18T01:41:44.975Z | 2987 samps
        # WMOK..B Z | 2011-01-18T01:40:30.325Z - 18T01:41:44.975Z | 2987 samps
        # X31A..B Z | 2011-01-18T01:40:30.325Z - 18T01:42:10.325Z | 4001 samps
        # WMOK..B Z | 2011-01-18T01:40:30.325Z - 18T01:42:10.325Z | 4001 samps

        # Remove traces from same station with different # samples
        while True:
            ST_new_sta = []
            for tr in ST_new:
                ST_new_sta.append(tr.stats.station)
            duplicate = list(set([tr for tr in ST_new_sta
                                  if ST_new_sta.count(tr) > 1]))
            if not duplicate:
                break
            index = [i for (i, j) in enumerate(ST_new_sta)
                     if j == duplicate[-1]]
            i = 0
            while True:
                if (ST_new[index[i]].stats.npts
                        < ST_new[index[i+1]].stats.npts):
                    del ST_new[index[i]]
                    break
                elif (ST_new[index[i]].stats.npts
                      >= ST_new[index[i+1]].stats.npts):
                    del ST_new[index[i+1]]
                    break
        # print ST_new.__str__(extended=True)

        ST_new.detrend('demean')
        # ST_new.filter('bandpass', freqmin=0.1, freqmax=100)

        segs, ticklocs, sta, circle_x, circle_y, segs_picks, ticklocs_picks, \
            data_picks = self._make_section_arrays(ST_new, scale_factor,
                                                   channel, seconds_ahead,
                                                   record_length, assoc_id,
                                                   eve)
        tick_min, tick_max = min(ticklocs), max(ticklocs)
        offsets = np.zeros((len(ST_new), 2), dtype=float)
        offsets[:, 0] = ticklocs
        offsets_picks = np.zeros((len(segs_picks), 2), dtype=float)
        offsets_picks[:, 0] = ticklocs_picks
        lines = LineCollection(segs, offsets=offsets, transOffset=None,
                               linewidths=.25, color='gray')
        lines_picks = LineCollection(segs_picks, offsets=offsets_picks,
                                     transOffset=None, linewidths=1,
                                     color='k')

        # Set up axis
        fig = plt.figure(figsize=(15, 8))
        ax1 = fig.add_subplot(111)
        x1 = tick_max + (tick_max - tick_min) * 0.1
        plt.ylim(0, record_length)
        plt.xlim(0, x1)

        # Plot data and picks
        ax1.plot(circle_x, circle_y, 'o', c='gray')  # dots where picks cross the waveform
        ax1.add_collection(lines)        # Plot the data
        ax1.add_collection(lines_picks)  # Plot the picks
        
        # Set axis labels and ticks
        ax1.set_xticks(ticklocs)
        ax1.set_xticklabels(sta)
        ax1.invert_yaxis()
        ax1.xaxis.tick_top()
        plt.setp(plt.xticks()[1], rotation=45)
        ax2 = ax1.secondary_xaxis('bottom')
        ax2.set_xlabel('Offset(km)')
        plt.ylabel('Record Length (s)', fontsize=18)
        plt.title(f'{channel}-channel Section Plot of Event at {tr.stats.starttime}')
        # plt.tight_layout()

        # Plot travel-time curves
        ttable = self.tt_stations_db_1D.query(TTtable1D).\
            filter(TTtable1D.d_km <= x1).all()
        D = [x.d_km for x in ttable]
        P = [x.p_tt + seconds_ahead for x in ttable]
        S = [x.s_tt + seconds_ahead for x in ttable]
        ax1.plot(D, P, 'b--', D, S, 'r--', linewidth=2)

        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()

    def _get_bounds(self, assoc_id, west, east, south, north, deltalon,
                    deltalat):
        """Set any non-specified map bounds, return as dict"""
        # Get event and all candidate station positions
        eve = self.assoc_db.query(Associated).\
            filter(Associated.id == assoc_id).first()
        stations = []
        for candidate in self.assoc_db.query(PickModified).\
                filter(PickModified.assoc_id == assoc_id).all():
            stations.append(self.tt_stations_db_1D.
                            query(Station1D).
                            filter(Station1D.sta == candidate.sta).first())
        lons = [eve.longitude] + [x.longitude for x in stations]
        lats = [eve.latitude] + [x.latitude for x in stations]
        lat_range = max(lats) - min(lats)
        lon_range = max(lons) - min(lons)
        if not west:
            west = min(lons) - lon_range/20
        if not east:
            east = max(lons) + lon_range/20
        if not south:
            south = min(lats) - lat_range/20
        if not north:
            north = max(lats) + lat_range/20
        
        # Adjust to avoid too-narrow plots
        lat_range = north - south
        lon_range = east - west
        if (lat_range) / (lon_range) < 0.6:
            north += (lon_range - lat_range)/3
            south -= (lon_range - lat_range)/3
        elif (lon_range) / (lat_range) < 0.6:
            east += (lat_range - lon_range)/3
            west -= (lat_range - lon_range)/3
        
        # Calculate tick spacing  
        if not deltalon:
            deltalon = self._tick_spacing(east - west)
        if not deltalat:
            deltalat = self._tick_spacing(north - south)
        # Adjust bounds to start at deltalat, deltalon intervals
        west, east = self._adjust_bounds(west, east, deltalon)
        south, north = self._adjust_bounds(south, north, deltalat)

        return dict(west=west, east=east, south=south, north=north,
                    deltalon=deltalon, deltalat=deltalat)


    @staticmethod
    def _adjust_bounds(below, above, step):
        below -= (below % step)
        if (above % step):
            above += step - (above % step)
        return below, above

    @staticmethod
    def _tick_spacing(ax_range, max_ticks=8):
        """
        Calculate the tick spacing corresponding to an axis range
        
        Try to have at least 3 and at most 10 ticks
        """
        # find a factor of ten that makes at least three ticks
        spacing = 10**(np.floor(np.log10(ax_range / 3)))
        # if it makes more than 10 ticks, multiply spacing by 5
        if ax_range / spacing > max_ticks:
            spacing *= 2  # spacing starts with "2"
            if ax_range / spacing > max_ticks:
                spacing *= 2.5 # spacing starts with "5"
        assert ax_range / spacing <= max_ticks, 'Error! too many ticks'
        return spacing
        
    def _plot_residuals(self, subax, s_p_list, max_dist=None, max_s_p=None):
        """Plot distribution of residuals
        
        :param s_p_list: list of dicts with keys 'time', 'radius' and 'color'
        :param max_distance: maximum distance to plot (None calculates from data)
        :param max_s_p: maximum S-P time to plot (None calculates from data)
        """
        Dist = self.tt_stations_db_1D.query(TTtable1D.d_km).all()
        SP_intv = self.tt_stations_db_1D.query(TTtable1D.s_p).all()
        S_P = [x[0] for x in SP_intv]
        D_S_P = [x[0] for x in Dist]
        # D_S_P = []
        # for dist, in Dist:
        #     D_S_P.append(dist)
        # for s_p, in SP_intv:
        #     S_P.append(s_p)

        subax.plot(D_S_P, S_P, 'k-', linewidth=2, zorder=1)

        for s_p in s_p_list:
        # for i in range(len(radius_rainbow)):
            subax.scatter(s_p['radius'], s_p['time'], s=50,
                          color=s_p['color'], zorder=2)
        if not max_dist:
            max_dist = max(D_S_P)
        if not max_s_p:
            max_s_p = max(S_P)

        plt.xlim([0, max_dist])
        plt.ylim([0, max_s_p])
        plt.xlabel('Distance (km)')
        plt.ylabel('S-P (s)')

    def _make_section_arrays(self, ST_new, scale_factor, channel,
                             seconds_ahead, record_length, assoc_id, eve):
        """
        Make arrays used in section plots
        
        First step to cleaning up code
        """
        eq_lat = eve.latitude
        eq_lon = eve.longitude
        segs = []
        ticklocs = []
        sta = []
        circle_x = []
        circle_y = []
        segs_picks = []
        ticklocs_picks = []
        for tr in ST_new:
            data = tr.data / (tr.data.max() - tr.data.min()) * scale_factor
            # due to float point arithmetic issue, cannot use
            # "t=np.arange(0, tr.stats.npts / tr.stats.sampling_rate,
            #              tr.stats.delta)"
            t = np.arange(0, round(tr.stats.npts / tr.stats.sampling_rate
                                   / tr.stats.delta)) * tr.stats.delta
            segs.append(np.hstack((data[:, np.newaxis], t[:, np.newaxis])))
            lon, lat = self.tt_stations_db_1D.\
                query(Station1D.longitude, Station1D.latitude).\
                filter(Station1D.sta == tr.stats.station).first()
            # gps2DistAzimuth returns in meters, convert to km by /1000
            distance = int(gps2dist_azimuth(lat, lon, eq_lat,
                                            eq_lon)[0] / 1000.)
            # distance = self.assoc_db.query(Candidate.d_km).\
            #   filter(Candidate.assoc_id==assoc_id).\
            #   filter(Candidate.sta==tr.stats.station).first()[0]
            # print distance,tr.stats.station
            ticklocs.append(distance)
            sta.append(tr.stats.station)
            # DOT plot where picks are picked, notice that for vertical trace
            # plot p is queried from Pick table, s from PickModified table
            # horizontal trace plot p and s queried from PickModified table
            if channel == 'Z3':
                picks_p = self.assoc_db.query(Pick.time).\
                    filter(Pick.assoc_id == assoc_id).\
                    filter(Pick.sta == tr.stats.station).\
                    filter(Pick.chan == tr.stats.channel).\
                    filter(Pick.phase == 'P').all()
                if not picks_p:
                    picks_p = self.assoc_db.query(PickModified.time).\
                        filter(PickModified.assoc_id == assoc_id).\
                        filter(PickModified.sta == tr.stats.station).\
                        filter(PickModified.phase == 'P').all()
                picks_s = self.assoc_db.query(PickModified.time).\
                    filter(PickModified.assoc_id == assoc_id).\
                    filter(PickModified.sta == tr.stats.station).\
                    filter(PickModified.phase == 'S').all()
            else:
                picks_p = self.assoc_db.query(PickModified.time).\
                    filter(PickModified.assoc_id == assoc_id).\
                    filter(PickModified.sta == tr.stats.station).\
                    filter(PickModified.phase == 'P').all()
                picks_s = self.assoc_db.query(PickModified.time).\
                    filter(PickModified.assoc_id == assoc_id).\
                    filter(PickModified.sta == tr.stats.station).\
                    filter(PickModified.phase == 'S').all()
            picks = picks_p + picks_s
            for pick, in picks:
                index = int((pick-eve.ot + timedelta(seconds=seconds_ahead)).
                            total_seconds()/tr.stats.delta)
                circle_x.append(distance + data[index])
                circle_y.append(t[index])
                # BAR plot where picks are picked
                t_picks = np.array([t[index], t[index]])
                data_picks = np.array([data.min(), data.max()])
                segs_picks.append(np.hstack((data_picks[:, np.newaxis],
                                             t_picks[:, np.newaxis])))
                ticklocs_picks.append(distance)

        return (segs, ticklocs, sta, circle_x, circle_y, segs_picks,
                ticklocs_picks, data_picks)


    def _distance_singlephase(self, phase, time):
        """
        Return the travel-time table value to a single-phase travel time

        Modified to return interpolated time rather than closest table value
        :param phase: phase to look for
        :param time: travel time (seconds)
        :returns: tt_object
        """
        if phase == 'P':
            closer = self.tt_stations_db_1D.query(TTtable1D).\
                filter(TTtable1D.p_tt <= time).\
                order_by(TTtable1D.p_tt.desc()).first()
            farther = self.tt_stations_db_1D.query(TTtable1D).\
                filter(TTtable1D.p_tt >= time).\
                order_by(TTtable1D.p_tt).first()
            if farther == closer:
                return closer
            alpha = (time - closer.p_tt)/(farther.p_tt - closer.p_tt)
            interpd = TTtable1D(closer.d_km*(1 - alpha) + farther.d_km*alpha,
                                closer.delta*(1 - alpha) + farther.delta*alpha,
                                closer.p_tt*(1 - alpha) + farther.p_tt*alpha,
                                closer.s_tt*(1 - alpha) + farther.s_tt*alpha,        
                                closer.s_p*(1 - alpha) + farther.s_p*alpha)
        elif phase == 'S':
            closer = self.tt_stations_db_1D.query(TTtable1D).\
                filter(TTtable1D.s_tt <= time).\
                order_by(TTtable1D.s_tt.desc()).first()
            farther = self.tt_stations_db_1D.query(TTtable1D).\
                filter(TTtable1D.s_tt >= time).\
                order_by(TTtable1D.s_tt).first()
            if farther == closer:
                return closer
            alpha = (time - closer.s_tt)/(farther.s_tt - closer.s_tt)
            interpd = TTtable1D(closer.d_km*(1 - alpha) + farther.d_km*alpha,
                                closer.delta*(1 - alpha) + farther.delta*alpha,
                                closer.p_tt*(1 - alpha) + farther.p_tt*alpha,
                                closer.s_tt*(1 - alpha) + farther.s_tt*alpha,        
                                closer.s_p*(1 - alpha) + farther.s_p*alpha)
        return interpd


def add_subplot_axes(ax, rect, axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position = ax.transAxes.transform(rect[0: 2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[2]
    subax = fig.add_axes([x, y, width, height])  #, axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def shoot(lon, lat, azimuth, maxdist=None):
    """
    Shooter Function

    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    :param lat:
    :param lon:
    :param azimuth:
    :param maxdist:
    :returns: glon2, glat2, baz
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS = 0.00000000005
    if ((np.abs(np.cos(glat1)) < EPS) and not (np.abs(np.sin(faz)) < EPS)):
        warnings.warn("Only N-S courses are meaningful, starting at a pole!")

    # a=6378.13/1.852
    a = 6371 / 1.852
    f = 1 / 298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf == 0):
        b = 0.
    else:
        b = 2. * np.arctan2(tu, cf)
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs(y - c) > EPS):
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
             d / 4. - cz) * sy * d + tu
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2 * np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180. / np.pi
    glat2 *= 180. / np.pi
    baz *= 180. / np.pi

    return (glon2, glat2, baz)


def _draw_circle(m, centerlon, centerlat, radius, *args, **kwargs):
    """
    Draw a circle?

    http://www.geophysique.be/2011/02/20/matplotlib-basemap-tutorial
         -09-drawing-circles/
    :param m:
    :param centerlon:
    :param centerlat:
    :param radius:
    """
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    # ~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X, Y = m(X, Y)
    plt.plot(X, Y, **kwargs)
    return plt.plot(X, Y, **kwargs)[0]  # just return for legend

#     def _distance_singlephase(self, phase, time):
#         """
#         Return the travel-time table value to a single-phase travel time
# 
#         :param phase: phase to look for
#         :param time: travel time (seconds)
#         :returns: closest tt_object, s_p_difference
#         """
#         if phase == 'P':
#             min = self.tt_stations_db_1D.query(TTtable1D).\
#                 filter(TTtable1D.p_tt <= time).\
#                 order_by(TTtable1D.p_tt.desc()).first()
#             max = self.tt_stations_db_1D.query(TTtable1D).\
#                 filter(TTtable1D.p_tt >= time).\
#                 order_by(TTtable1D.p_tt).first()
#             if abs(min.p_tt - time) <= abs(max.p_tt - time):
#                 return min, abs(min.p_tt - time)
#             else:
#                 return max, abs(max.p_tt - time)
# 
#         if phase == 'S':
#             min = self.tt_stations_db_1D.query(TTtable1D).\
#                 filter(TTtable1D.s_tt <= time).\
#                 order_by(TTtable1D.s_tt.desc()).first()
#             max = self.tt_stations_db_1D.query(TTtable1D).\
#                 filter(TTtable1D.s_tt >= time).\
#                 order_by(TTtable1D.s_tt).first()
#             if abs(min.s_tt - time) <= abs(max.s_tt - time):
#                 return min, abs(min.s_tt - time)
#             else:
#                 return max, abs(max.s_tt - time)
# 
