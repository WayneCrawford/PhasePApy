#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import os
import glob
import unittest
import inspect
import sys
import difflib
# from pprint import pprint
# import xml.etree.ElementTree as ET
# from CompareXMLTree import XmlTree
# from obsinfo.network.network import _make_stationXML_script
from obspy import Catalog
from obspy.io.nordic.core import read_nordic

# sys.path.append("/home/soumya/programming/PhasePApy/")
sys.path.append("../../")
from phasepapy.associator import tables1D, assoc1D, plot1D
from phasepapy.associator.tables1D import make_assoc_session
from phasepapy.associator.tt_stations_1D import (make_tt1D_session)


class TestADDONSMethods(unittest.TestCase):
    """
    Test suite for phasepapy operations.
    """
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(inspect.getfile(
            inspect.currentframe())))
        self.testing_path = os.path.join(self.path, "data")

    def assertTextFilesEqual(self, first, second, msg=None):
        with open(first) as f:
            str_a = f.read()
        with open(second) as f:
            str_b = f.read()

        if str_a != str_b:
            first_lines = str_a.splitlines(True)
            second_lines = str_b.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=first, tofile=second)
            message = ''.join(delta)

            if msg:
                message += " : " + msg

            self.fail("Multi-line strings are unequal:\n" + message)

    def test_fullrun_NORDIC(self):
        """
        Test Read in catalog from NORDIC and run with a tt_stations_1D file
        """
        db_tt = 'sqlite:///' + os.path.join(self.testing_path, 'tt_lsv_1D.db')
        assoc_params = dict(max_km=80, aggregation=1, aggr_norm='L2',
                            assoc_ot_uncert=1, nsta_declare=2,
                            cutoff_outlier=10, loc_uncert_thresh=0.1)
        catalog_file = os.path.join(self.testing_path, 'test_catalog.nordic')
        db_assoc_file = 'assoc_1D.db'
        db_assoc_url = 'sqlite:///' + db_assoc_file
        events, wavefiles = read_nordic(catalog_file, True)
        txt = ''
        for event, wavefile in zip(events, wavefiles):
            if os.path.exists(db_assoc_file):
                os.remove(db_assoc_file)
            dbsession = make_assoc_session(db_assoc_url)
            for pick in event.picks:
                my_pick = tables1D.Pick.from_obspy(pick)
                dbsession.add(my_pick)
                dbsession.commit()

            assoc = assoc1D.LocalAssociator(db_assoc_url, db_tt, **assoc_params)
            assoc.id_candidate_events()
            assoc.associate_candidates()
            if assoc.count_associated():
                assoc.single_phase()
            txt += str(assoc) + '\n'
        with open('temp.txt', 'w') as f:
            f.write(txt)
        self.assertTextFilesEqual('temp.txt',
                                  os.path.join(self.testing_path,
                                               'test_catalog_out.txt'))
        os.remove('temp.txt')
        os.remove('assoc_1D.db')

    def test_count_picks(self):
        """
        Test counting of picks in an associator database
        """
        db_tt = 'sqlite:///' + os.path.join(self.testing_path, 'tt_lsv_1D.db')
        db_tables = 'sqlite:///' + os.path.join(self.testing_path,
                                                'assoc7_example_1D.db')
        assoc = assoc1D.LocalAssociator(db_tables, db_tt, max_km=80)
        self.assertEqual(assoc.count_picks(), 15)
        self.assertEqual(assoc.count_pick_modifieds(), 12)
        self.assertEqual(assoc.count_pick_modifieds('associated'), 7)
        self.assertEqual(assoc.count_pick_modifieds('associated', phase='P'), 3)
        self.assertEqual(assoc.count_pick_modifieds('associated', phase='S'), 4)
        self.assertEqual(assoc.count_pick_modifieds('matched'), 4)
        self.assertEqual(assoc.count_pick_modifieds('mismatched'), 0)
        self.assertEqual(assoc.count_pick_modifieds('single-phase'), 3)


def suite():
    return unittest.makeSuite(TestADDONSMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
