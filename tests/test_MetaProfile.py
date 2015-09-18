#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_MetaProfile
----------------------------------

Tests for `MetaProfile` module.
"""

import unittest
from MetaProfile import utils as metafun
from MetaProfile import Signal, Window, Profile, MetaProfiler, get_profiler, get_signals


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
                'filepath': "test3.bed"})
windows_list = ({'name': "PolyA-signal",
                 'filepath': "test_windows.bed",
                 'pseudocount': 0},)

pc = get_profiler(signal_list, windows_list)
