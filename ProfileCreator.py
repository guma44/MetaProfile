import os
import sys
import numpy as np
from pandas import DataFrame
import matplotlib.gridspec as gridspec
import pylab as pl
import seaborn as sns
from Bio import SeqIO, Seq, motifs
from Bio.Alphabet import IUPAC
import HTSeq


class WrongFiletypeException(Exception): pass
class NotSuchAFileException(Exception): pass
class NoSignalsException(Exception): pass
class NoProfilesException(Exception): pass
class NotSuchAProfileException(Exception): pass
class NotSuchAProfileTypeException(Exception): pass

class ProfileCreator():

    """A class that handles creation of the genomic profile"""

    def __init__(self, windows_path, windows_filetype, window_length=200, pseudocount=1):
        """Initialize ProfileCreator

        :param signals_path: @todo
        :param signals_filetype: @todo
        :param windows_path: @todo
        :param windows_filetype: @todo
        :param window_length: @todo
        :param pseudocount: @todo

        """

        self.window_length = window_length
        self.pseudocount = pseudocount
        # self.set_signals(signals_path, signals_filetype)
        self.set_windows(windows_path, windows_filetype)
        self.signals = {}
        self.profiles = {}

    def create_genomic_signals(self, path_to_file, filetype, stranded=True):
        """Prepares coverage as a HTSeq.GenomicArray

        :param filepath: path to file
        :param filetype: type of the file (can be bed etc.)
        :returns: 2-tuple of HTSeq.GenomicArray, library_size

        """
        sys.stderr.write("Creating the signal. It may take few minutes...\n")
        coverage = HTSeq.GenomicArray("auto", stranded=stranded, typecode="d")
        library_size = 0
        if filetype.upper() == "BED":
            for line in HTSeq.BED_Reader(path_to_file):
                coverage[line.iv] += 1
                library_size += 1
        elif filetype.upper() == "GFF" or filetype.upper() == "GTF":
            for line in HTSeq.GFF_Reader(path_to_file):
                coverage[line.iv] += 1
                library_size += 1
        elif filetype.upper() == "SAM":
            for line in HTSeq.SAM_Reader(path_to_file):
                coverage[line.iv] += 1
                library_size += 1
        elif filetype.upper() == "BAM":
            for line in HTSeq.BAM_Reader(path_to_file):
                coverage[line.iv] += 1
                library_size += 1
        return coverage, library_size

    def get_windows_iterator(self):
        """Return an generator over intervals from the file including extension

        :param filepath: path to file
        :param filetype: type of the file (can be bed etc.)
        :param window_length: length of the window to add on both sides
        :returns: generator over expanded GenomicIntervals

        """
        filepath = self.windows_path
        filetype = self.windows_filetype
        window_length = self.window_length
        class Iterable(object):
            def __iter__(self):
                if filetype.upper() == "BED":
                    for line in HTSeq.BED_Reader(filepath):
                        line.iv.start -= window_length
                        line.iv.end += window_length
                        yield line.iv
                elif filetype.upper() == "GFF" or filetype.upper() == "GTF":
                    for line in HTSeq.GFF_Reader(filepath):
                        line.iv.start -= window_length
                        line.iv.end += window_length
                        yield line.iv
                elif filetype.upper() == "SAM":
                    for line in HTSeq.SAM_Reader(filepath):
                        line.iv.start -= window_length
                        line.iv.end += window_length
                        yield line.iv
                elif filetype.upper() == "BAM":
                    for line in HTSeq.BAM_Reader(filepath):
                        line.iv.start -= window_length
                        line.iv.end += window_length
                        yield line.iv
        self.windows = Iterable()

    def create_profiles(self):
        """Prepare array of profiles for each each signal. For each window a pseudocount
        is added.

        :param windows: iterator that returns a windows of equal length as HTSeq.GenomicInterval
        :param signals: HTSeq.GenomicArray of signals over chromosomes
        :param pseudocount: add a pseudocount to the window to avoid dividing by 0 in some operations,
                            defaults to 1
        :returns: numpy array of signals count for each window

        """
        if not self.signals:
            raise NoSignalsException("No signals detected. Add some genomic signals before using this function.\n")

        for signal_name, (signals_coverage, library_size) in self.signals.iteritems():
            sys.stderr.write("Calculating profile for %s" % signal_name)
            profile = []
            window_lengths = set()
            for window in self.windows:
                window_lengths.add(window.length)
                assert len(window_lengths) == 1, "Windows has at least two different lengths"
                try:
                    wincvg = np.fromiter(signals_coverage[window], dtype='i', count=window.length)
                except IndexError:
                    sys.stderr.write("Wrong window: %s\n" % str(window))
                    continue
                if window.strand == "+":
                    profile.append(wincvg + self.pseudocount)
                else:
                    profile.append(wincvg[::-1] + self.pseudocount)
            self.profiles[signal_name] = {"raw": np.asarray(profile)}

    def normalize_to_library(self, metric='mean', aggregate=True):
        """Normalize profile to the library size

        :param profile: numpy array of windows, could be from make_profile_array
        :param library_size: size of the library used to generate signals
        :param metric: metric to use when aggregating profile, can be sum or mean (average),
                       defaults to mean
        :param aggregate: if True the data will be aggregated into 1-D vector according to
               metric, if False metric will be ignored, defaults to True
        :returns: @todo

        """
        if not self.profiles:
            raise NoProfilesException("No profiles detected. Create raw profiles with create_profiles before using this function.\n")
        for profile_name, profile in self.profiles.iteritems():
            self.profiles[profile_name]['normalized_to_gene'] =  np.asarray(map(lambda i: i/float(self.signals[profile_name][1]), profile))

    def normalize_to_gene(self):
        """Normalize profile to the library size and to the gene count i.e.
        normalize each window to the total count in it by dividing each
        entry by the sum of the counts.

        :param profile: numpy array of windows, could be from make_profile_array
        :param metric: metric to use when aggregating profile, can be sum or mean (average),
                       defaults to mean
        :param aggregate: if True the data will be aggregated into 1-D vector according to
               metric, if False metric will be ignored, defaults to True
        :returns: normalized profile as np.array

        """
        if not self.profiles:
            raise NoProfilesException("No profiles detected. Create raw profiles with create_profiles before using this function.\n")
        for profile_name, profile in self.profiles.iteritems():
            adjusted_profile = np.asarray(map(lambda i: i/float(i.sum()), profile['raw']))
            self.profiles[profile_name]['normalized_to_gene'] = adjusted_profile

    def get_aggregated_profile(self, profile_name, profile_type, metric='mean'):
        """Get the profile that is aggregated i.e. summarized along
        all windows.

        :param profile_name: @todo
        :param profile_type: @todo
        :param metric: @todo
        :returns: @todo

        """
        if profile_name not in self.profiles:
            raise NotSuchAProfileException("No profile named %s. Available profiles are: %s" % (profile_name, ", ".join(map(str, self.profiles.keys()))))
        if profile_type not in self.profiles[profile_name]:
            raise NotSuchAProfileTypeException("No profile type %s. Available types are: %s" % (profile_type, ", ".join(map(str, self.profiles[profile_name].keys()))))
        if metric.upper() == "MEAN" or metric.upper() == "AVERAGE":
            return self.profiles[profile_name][profile_type].mean(axis=0)
        elif metric.upper() == "SUM":
            return self.profiles[profile_name][profile_type].sum(axis=0)
        else:
            raise NotSuchAMetricException("%s is not valid metric. Available are sum or mean.")

    def set_window_length(self, window_length):
        """
         Set window length
        """
        self.window_length = window_length

    def set_windows_path(self, windows_path):
        """
         Set path to windows file
        """
        if not os.path.isfile(windows_path):
            raise NotSuchAFileException("No file named %s" % signals_path)
        self.windows_path = windows_path

    def set_windows_filetype(self, windows_filetype):
        """
         Set windows type of file (BED, GFF, BAM, SAM)
        """
        if windows_filetype.upper() not in ["BED", "GFF", "GTF", "BAM", "SAM"]:
            raise WrongFiletypeException("Wrong filetype, allowed: %s" % "(BED, GTF, GFF, BAM, SAM)")
        self.windows_filetype = windows_filetype

    def set_windows(self, filepath, filetype):
        """Set windows filepath and filetype

        :param filepath: @todo
        :param filetype: @todo
        :returns: @todo

        """
        self.set_windows_filetype(filetype)
        self.set_windows_path(filepath)
        self.get_windows_iterator()

    def add_signal(self, name, path, filetype, force=False):
        """Add a signal to ProfileCreator. Can be used many times
        to add more signals.

        :param name: @todo
        :param path: @todo
        :param filetype: @todo
        :param force: @todo
        :returns: @todo

        """
        if not os.path.isfile(path):
            raise NotSuchAFileException("No file named %s" % signals_path)
        if filetype.upper() not in ["BED", "GFF", "GTF", "BAM", "SAM"]:
            raise WrongFiletypeException("Wrong filetype, allowed: %s" % "(BED, GTF, GFF, BAM, SAM)")
        if name in self.signals and not force:
            raise KeyError("%s already is already present in signals. Use force=True to override." % name)
        self.signals[name] = self.create_genomic_signals(path, filetype)
