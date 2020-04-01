#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib
import copy
import numpy as np
from obspy.core import UTCDateTime
from .scnl import SCNL
from .cf_kt import Kurtosis
from .util import rms, rolling_window


class KTPicker():
    """
    KTPicker is designed based on kurtosis.
    """

    def __init__(self, t_win=1, t_ma=10, nsigma=6, t_up=0.2, nr_len=2,
                 nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3):
        """
        :param t_win: Time (sec) of moving window to calculate kurtosis
        :param t_ma: Time (sec) of  moving window for dynamic threshold
        :param nsigma: Threshold ratio (x t_ma window rms) for potential pick
        :param t_up: Minimum spacing (sec) between potential picks
        :param nr_len: Signal-to-noise ratio window length (secs)
        :param nr_coeff: Signal-to-noise ratio required to accept potential
                         pick
        :param pol_len: Polarity: # of samples used to calculate pre-pick STD
        :param pol_coeff: Polarity: Threshold ratio between first extreme and
                          STD
        :param uncert_coeff: CF threshold ratio for identifying picks?

        :uncert_len: Window length in time to calculate the rms of the CF
                     before the picks. Fixed equal to t_ma
        """
        self.t_win = t_win
        self.t_ma = t_ma
        self.nsigma = nsigma
        self.t_up = t_up
        self.nr_len = nr_len
        self.nr_coeff = nr_coeff
        self.pol_len = pol_len
        self.pol_coeff = pol_coeff
        self.uncert_coeff = uncert_coeff
        self.uncert_len = self.t_ma

    def picks(self, tr):
        """
        Make picks, polarity, snr, and uncertainty.
        """
        # tr = trace.detrend('linear')
        summary = KTSummary(self, tr)
        # threshold
        # threshold = summary.threshold
        # picks
        scnl, picks, trigger, snr = summary.pick_ident()
        # uncertainty
        uncertainty = summary.uncertainty()
        # polarity
        polarity = summary.polarity()
        return scnl, picks, polarity, snr, uncertainty


class KTSummary():
    """
    Class to calculate CF

    Also calculates threshold level, cleans the false picks, determines
    uncertainty, polarity and plots CF.
    """

    def __init__(self, picker, trace):
        self.picker = picker
        self.tr = trace
        self.stats = self.tr.stats
        self.t_win = self.picker.t_win
        self.cf = Kurtosis(self.tr, self.t_win)
        self.FC = self.cf._statistics()
        self.summary = self.FC
        # Create the threshold from our picker config
        self.thres = self.threshold()
        self.uncert = self.uncertainty()
        self.pol = self.polarity()

    def threshold(self):
        """
        Control the threshold level with nsigma.
        """
        dt = self.stats.delta
        npts_Tma = int(round(self.picker.t_ma / dt, 0))
        LEN = self.stats.npts
        # print "npts_Tma: ",npts_Tma
        threshold = np.zeros(LEN)
        # threshold[0:npts_Tma] = 1
        threshold[npts_Tma:LEN] = rms(rolling_window(self.summary[0:LEN - 1],
                                                     npts_Tma),
                                      -1) * self.picker.nsigma
        return threshold

    def pick_ident(self):
        """
        Clean false picks and Make picks.
        """
        scnl = SCNL([self.stats.station, self.stats.channel,
                     self.stats.network, self.stats.location])
        dt = self.stats.delta
        npts_Tma = int(round(self.picker.t_ma/dt, 0))
        LEN = self.stats.npts

        # trigger the earthquakes
        trigger_ptnl_index = np.where(self.summary[npts_Tma:LEN] >
                                      self.thres[npts_Tma:LEN])
        trigger_ptnl_index = trigger_ptnl_index + np.array(npts_Tma)
        t = np.arange(0, self.stats.npts / self.stats.sampling_rate, dt)
        trigger_ptnl = t[trigger_ptnl_index][0]

        # clean close picks
        window_t_up = int(round(self.picker.t_up / dt, 0))
        trigger_remove1_index = []
        for i in range(0, len(trigger_ptnl) - 1):  # second from last
            # avoid consecutive picking
            if (trigger_ptnl[i+1] - trigger_ptnl[i]) <= window_t_up * dt:
                trigger_remove1_index.append(i + 1)
        # delete close picks
        trigger_ptnl = np.delete(trigger_ptnl, trigger_remove1_index)

        # clean_filtering
        trigger_remove2_index = []
        N = self.picker.nr_coeff
        filter_length = self.picker.nr_len
        for i in range(len(trigger_ptnl)):
            # determine filter_length for each pick:
            r, R = self.winlen(i, trigger_ptnl, filter_length, t, dt)
            M = min(r, R)
            trig_addr = int(round(trigger_ptnl[i] / dt, 0))
            # if N * np.std(self.tr.data[int(round(trigger_ptnl[i] / dt, 0))\
            # - M:int(round(trigger_ptnl[i] / dt, 0))]) >=\
            # np.std(self.tr[int(round(trigger_ptnl[i] / dt, 0)):\
            # int(round(trigger_ptnl[i] / dt, 0)) + M]):
            if N * np.std(self.tr.data[trig_addr - M:trig_addr]) >=\
               np.std(self.tr[trig_addr:trig_addr + M]):
                trigger_remove2_index.append(i)
        # delete fake picks
        trigger_ptnl = np.delete(trigger_ptnl, trigger_remove2_index)

        # assign potential picks to trigger
        trigger = trigger_ptnl

        # Need to be careful copy list or array, since if just use
        # A=B, when any element in B changes, A will change as well
        # roll backward for picking
        picks = []
        for i in range(len(trigger)):
            index = int(round(trigger[i] / dt, 0))
            while True:
                if self.summary[index] > self.summary[index - 1]:
                    index -= 1
                else:
                    break
            picks.append(UTCDateTime(self.tr.stats.starttime +
                                     round(t[index], 3)))

        # Need to be careful copy list or array, since if just copy
        # like A=B, when any element in B changes, A will change as well
        # roll forward for maximum signal values
        maxes = copy.deepcopy(trigger)
        for i in range(len(trigger)):
            index = int(round(trigger[i] / dt, 0))
            while True:
                if self.summary[index] < self.summary[index + 1]:
                    index += 1
                else:
                    break
            maxes[i] = round(self.summary[index], 3)

        # Need to be careful copying lists or arrays: if we use A=B
        # A will change when any element in B changes
        # Signal noise ratio: SNR
        SNR = copy.deepcopy(trigger)

        for i in range(len(picks)):
            index = int(round(trigger[i] / dt, 0))
            noise = rms(self.summary[index-npts_Tma:index])
            SNR[i] = round(maxes[i] / noise, 1)

        return scnl, picks, trigger, SNR

    def uncertainty(self):
        """
        Uncertainty based on the noise level of CF.
        """
        scnl, picks, trigger, SNR = self.pick_ident()
        dt = self.stats.delta
        npts_Tma = int(round(self.picker.uncert_coeff / dt, 0))
        t = np.arange(0, self.stats.npts / self.stats.sampling_rate, dt)
        pick_uncert = copy.deepcopy(trigger)
        
        # print(scnl, len(picks), len(trigger), len(SNR), flush=True)

        for i in range(len(trigger)):
            # print(i, flush=True)
            r, R = self.winlen(i, trigger, npts_Tma, t, dt)
            index0 = int(round((picks[i] - self.tr.stats.starttime) / dt, 0))
            index = int(round(trigger[i] / dt, 0))
            uncert_level = self.picker.uncert_coeff * rms(
                                self.summary[index0 - npts_Tma: index0])
            while True:
                if self.summary[index] > uncert_level and\
                   self.summary[index] > self.summary[index - 1]:
                    index -= 1
                else:
                    break
            pick_uncert[i] = round(t[index] - 
                                   (picks[i] - self.tr.stats.starttime),
                                   3)
        return pick_uncert

    def polarity(self):
        """
        Determine polarity for declared picks.
        """
        dt = self.stats.delta
        # t = np.arange(0, self.stats.npts / self.stats.sampling_rate, dt)
        pol = []
        scnl, picks, trigger, snr = self.pick_ident()
        for i in range(len(picks)):
            index0 = int(round((picks[i] - self.tr.stats.starttime)/dt, 0))
            index = index0

            # roll forward index+=1
            while True:
                if index >= self.stats.npts - 1 - 2:
                    break
                elif (self.tr[index+1] - self.tr[index]) *\
                     (self.tr[index+2] - self.tr[index+1]) > 0:
                    index += 1
                else:
                    break

            # if (self.tr[index+1] - self.tr[index0]) > 0 and\
            #         abs(self.tr[index+1] - self.tr[index0]) >\
            #         (self.picker.pol_coeff *
            #          np.std(self.tr[index0 - self.picker.pol_len: index0])):
            #     polarity = 'C'
            # elif ((self.tr[index+1] - self.tr[index0]) < 0) and\
            #      abs(self.tr[index+1] - self.tr[index0]) >\
            #      (self.picker.pol_coeff *
            #       np.std(self.tr[index0 - self.picker.pol_len: index0])):
            #     polarity = 'D'

            # notice index+1, rolling stop one point before extreme, compare
            # with std to avoid very small
            delta_peak = self.tr[index+1] - self.tr[index0]
            if abs(delta_peak) > (self.picker.pol_coeff *
                                  np.std(self.tr[index0 - self.picker.pol_len:
                                                 index0])):
                if delta_peak > 0:
                    polarity = 'C'
                else:
                    polarity = 'D'
            else:
                polarity = ''
            pol.append(polarity)
        return pol

    def winlen(self, index, trigger_ptnl, filter_length, t, dt):
        """
        Determine the filter window length.

        If the time difference between two picks is less
        than window length, use the picks interval as window.
        """
        i = index
        if len(trigger_ptnl) == 1:
            # print 'A'
            if trigger_ptnl[i] <= filter_length:
                r = int(round(trigger_ptnl[i] / dt, 0))
            if trigger_ptnl[i] > filter_length:
                r = int(round(filter_length / dt, 0))
            if trigger_ptnl[i] + filter_length >= t[-1]:
                R = int(round((t[-1] - trigger_ptnl[i]) / dt, 0))
            if trigger_ptnl[i] + filter_length < t[-1]:
                R = int(round(filter_length / dt, 0))
        elif len(trigger_ptnl) > 1:
            # print 'B'
            if i == 0:
                if trigger_ptnl[i] <= filter_length:
                    r = int(round(trigger_ptnl[i] / dt, 0))  # print 'a'
                if trigger_ptnl[i] > filter_length:
                    r = int(round(filter_length / dt, 0))  # print 'b'
                if (trigger_ptnl[i+1] - trigger_ptnl[i]) <= filter_length:
                    R = int(round((trigger_ptnl[i+1] -
                                   trigger_ptnl[i]) / dt, 0))  # print 'c'
                if (trigger_ptnl[i+1] - trigger_ptnl[i]) > filter_length:
                    R = int(round(filter_length / dt, 0))  # print 'd'
            elif i > 0 and i < len(trigger_ptnl) - 1:
                if trigger_ptnl[i] - trigger_ptnl[i-1] <= filter_length:
                    r = int(round((trigger_ptnl[i] - trigger_ptnl[i-1]) / dt,
                                  0))  # print 'e'
                if trigger_ptnl[i] - trigger_ptnl[i-1] > filter_length:
                    r = int(round(filter_length / dt, 0))  # print 'f'
                if trigger_ptnl[i+1] - trigger_ptnl[i] <= filter_length:
                    R = int(round((trigger_ptnl[i + 1] - trigger_ptnl[i]) / dt,
                                  0))
                    # print 'g'
                if trigger_ptnl[i + 1] - trigger_ptnl[i] > filter_length:
                    R = int(round(filter_length / dt, 0))
                    # print 'h'
            elif i == len(trigger_ptnl) - 1:
                if trigger_ptnl[i] - trigger_ptnl[i - 1] <= filter_length:
                    r = int(round((trigger_ptnl[i] - trigger_ptnl[i - 1]) / dt,
                                  0))
                    # print 'i'
                if trigger_ptnl[i] - trigger_ptnl[i - 1] > filter_length:
                    r = int(round(filter_length / dt, 0))
                    # print 'j'
                if trigger_ptnl[i] + filter_length > t[-1]:
                    R = int(round((t[-1] - trigger_ptnl[i]) / dt, 0))
                    # print 'k'
                if trigger_ptnl[i] + filter_length <= t[-1]:
                    R = int(round(filter_length / dt, 0))
                    # print 'l'
        return r, R

    def plot_picks(self):
        """
        Plot picks and waveform.
        """
        matplotlib.rcParams["axes.labelsize"] = "large"
        matplotlib.rcParams["axes.linewidth"] = 2.0
        matplotlib.rcParams["xtick.major.size"] = 8
        matplotlib.rcParams["ytick.major.size"] = 8
        matplotlib.rcParams["ytick.minor.size"] = 5
        matplotlib.rcParams["xtick.labelsize"] = "large"
        matplotlib.rcParams["ytick.labelsize"] = "large"

        dt = self.stats.delta
        t = np.arange(0, self.stats.npts/self.stats.sampling_rate, dt)
        scnl, picks, trigger, snr = self.pick_ident()
        # fig = plt.figure(figsize=(10, 5))
        plt.figure(figsize=(10, 5))
        plt.plot(t, self.tr, c='gray')
        for i in range(len(picks)):
            plt.plot([(picks[i] - self.tr.stats.starttime),
                      (picks[i] - self.tr.stats.starttime)],
                     [min(self.tr), max(self.tr)], 'k--')
            plt.text((picks[i] - self.tr.stats.starttime),
                     max(self.tr) - 0.3 * (max(self.tr) - min(self.tr)),
                     '%s' % (self.pol[i]), color='black')
        plt.xlabel('Time (s)')
        plt.show()

    def plot_summary(self):
        """
        Plot CF.
        """
        matplotlib.rcParams["axes.labelsize"] = "large"
        matplotlib.rcParams["axes.linewidth"] = 2.0
        matplotlib.rcParams["xtick.major.size"] = 8
        matplotlib.rcParams["ytick.major.size"] = 8
        matplotlib.rcParams["ytick.minor.size"] = 5
        matplotlib.rcParams["xtick.labelsize"] = "large"
        matplotlib.rcParams["ytick.labelsize"] = "large"

        # fig = plt.figure(figsize=(12, 8))
        plt.figure(figsize=(12, 8))
        dt = self.stats.delta
        t = np.arange(0, self.stats.npts / self.stats.sampling_rate, dt)

        # Plot raw data
        ax = plt.subplot(2, 1, 1)
        ax.plot(t, self.tr, c='gray')
        plt.ylabel('Raw Data')

        # Plot summary
        ax1 = plt.subplot(2, 1, 2)
        ax1.plot(t, self.summary / np.amax(self.summary), c='k')
        ax1.plot(t, self.thres / np.amax(self.summary), '--', linewidth=2.0,
                 c='k')
        plt.ylabel('Characteristic Function')
        ax1.legend(('Normalized CF', 'Threshold', 'Picks'), 'upper right',
                   shadow=True, fancybox=True)

        # scnl, picks, trigger, snr=self.pick_ident()
        # for i in range(len(picks)):
        #   ax.plot([(picks[i]-self.tr.stats.starttime),
        #            (picks[i]-self.tr.stats.starttime)],
        #           [min(self.tr.data),max(self.tr.data)],'r--')
        #   ax.text((picks[i]-self.tr.stats.starttime),0.5,
        #           '%s' % (self.pol[i]),color='red')

        plt.xlabel('Time (s)')
        # plt.title('CF')
        plt.tight_layout()
        plt.show()
