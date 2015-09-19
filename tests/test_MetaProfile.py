#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_MetaProfile
----------------------------------

Tests for `MetaProfile` module.
"""

import os
import unittest
import time
import sys
sys.path.append("..")
from MetaProfile import utils as metafun
from MetaProfile import Signal, Window, Profile, MetaProfiler, get_profiler, get_signals
os.system("taskset -p 0xff %d" % os.getpid())


class TestMetaprofile(unittest.TestCase):

    def setUp(self):
        pass

    def test_something(self):
        pass

    def tearDown(self):
        pass

signal_list = ({'name': "Test-1",
                'filepath': "test1.bed"},
               {'name': "Test-2",
                'filepath': "test2.bed"},
               {'name': "Test-3",
                'filepath': "test1.bed"},
               {'name': "Test-4",
                'filepath': "test2.bed"},
               {'name': "Test-5",
                'filepath': "test3.bed"},
               {'name': "Test-6",
                'filepath': "test2.bed"},
               {'name': "Test-7",
                'filepath': "test3.bed"})
windows_list = ({'name': "PolyA-signal",
                 'filepath': "test_windows.bed",
                 'pseudocount': 0},)

# print "########################################################################################################"
# start_time = time.time()
# start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
# sys.stderr.write("############## Started script on %s ##############\n" % start_date)

# pc = get_profiler(signal_list, windows_list)

# sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
# print "########################################################################################################"
# start_time = time.time()
# start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
# sys.stderr.write("############## Started script on %s ##############\n" % start_date)

# pc = get_profiler(signal_list, windows_list, parallel=False)

# sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
# print "########################################################################################################"

print "########################################################################################################"
start_time = time.time()
start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
sys.stderr.write("############## Started script on %s ##############\n" % start_date)
signal_list = ({'name': "CPSF-160-60min",
                'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_60/mapped_sequences.all.bed"},
               {'name': "CPSF-160-120min",
                'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_120/mapped_sequences.all.bed"},
               {'name': "CPSF-160-30min",
                'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_30/mapped_sequences.all.bed"},
               {'name': "CPSF-160-240min",
                'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_240/mapped_sequences.all.bed"})

windows_list = ({'name': "PolyA-signal",
                 'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/AdditionalData/PolyAsignal/PolyAsignal_from_asymetric.bed",
                 'pseudocount': 0})
pc = get_profiler(signal_list, windows_list)
sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))

print "########################################################################################################"
start_time = time.time()
start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
sys.stderr.write("############## Started script on %s ##############\n" % start_date)

signal60 = Signal("CPSF-160-60min",
                       "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_60/mapped_sequences.all.bed")
signal120 = Signal("CPSF-160-120min",
                       "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_120/mapped_sequences.all.bed")
signal30 = Signal("CPSF-160-30min",
                       "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_30/mapped_sequences.all.bed")
signal240 = Signal("CPSF-160-240min",
                       "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_240/mapped_sequences.all.bed")

windows = Window("PolyA-signal",
                        "/import/bc2/home/zavolan/gumiennr/Ule/Ule/AdditionalData/PolyAsignal/PolyAsignal_from_asymetric.bed",
                       pseudocount=0)

pc = MetaProfiler([signal30, signal60, signal120, signal240], [windows])
pc.create_profiles(normalize_to_library=True)

sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
print "########################################################################################################"
