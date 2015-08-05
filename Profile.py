import os
import sys
from numpy import fromiter, asarray

class WrongFiletypeException(Exception): pass
class NoSuchAFileException(Exception): pass
class NoSuchAProfileException(Exception): pass
class NoSuchAProfileTypeException(Exception): pass

class Profile(object):

    """Profile of the ..."""

    def __init__(self, signal, window):
        """Initialize profile object

        :param signals: @todo
        :param windows: @todo

        """
        #TODO assert for signal and windows type
        self.signal = signal
        self.window = window
        self.name = "%s on %s" % (signal.name, window.name)
        self.profile_raw = None
        self.profile_normalized_to_gene = None
        self.profile_normalized_to_library = None
        self.__create_profile()

    def __create_profile(self):
        """Prepare profile with given signal and windows. For each window a pseudocount
        is added.
        """
        sys.stderr.write("Calculating profile for %s on %s\n" % (self.signal.name, self.window.name))
        profile = []
        window_lengths = set()
        for window in self.window.windows:
            window_lengths.add(window.length)
            assert len(window_lengths) == 1, "Windows has at least two different lengths"
            try:
                wincvg = fromiter(self.signal.coverage[window], dtype='i', count=window.length)
            except IndexError:
                sys.stderr.write("Wrong window: %s\n" % str(window))
                continue
            if window.strand == "+":
                profile.append(wincvg + self.window.pseudocount)
            else:
                profile.append(wincvg[::-1] + self.window.pseudocount)

        self.profile_raw = asarray(profile)

    def normalize_to_library(self):
        """Normalize profile to the library size
        """
        self.profile_normalized_to_gene = asarray(map(lambda i: i/float(self.signal.library_size), self.profile_raw))

    def normalize_to_gene(self):
        """Normalize profile to the library size and to the gene count i.e.
        normalize each window to the total count in it by dividing each
        entry by the sum of the counts.
        """
        adjusted_profile = asarray(map(lambda i: i/float(i.sum()), self.profile_raw))
        self.profile_normalized_to_gene = adjusted_profile

    def get_aggregated_profile(self, profile_type, metric='mean'):
        """Get the profile that is aggregated i.e. summarized along
        all windows.
        :param profile_type: @todo
        :param metric: @todo
        :returns: @todo

        """
        if profile_type == "raw":
            return self.__agregate(self.profile_raw, metric)
        elif profile_type == "normalized_to_gene":
            if self.profile_normalized_to_gene is not None:
                return self.__agregate(self.profile_normalized_to_gene, metric)
            else:
                raise NoSuchAProfileTypeException("Profile is not normalized by default. "
                                            "First run Profile.normalize_to_gene to normalize it")
        elif profile_type == "normalized_to_library":
            if self.profile_normalized_to_gene is not None:
                return self.__agregate(self.profile_normalized_to_library, metric)
            else:
                raise NoSuchAProfileTypeException("Profile is not normalized by default. "
                                            "First run Profile.normalize_to_library to normalize it")
        else:
            raise NoSuchAProfileTypeException("No profile type %s. " % profile_type + \
                    "Available types are: raw, normalized_to_gene and normalized_to_library")

    def __agregate(self, profile, metric):
        if metric.upper() == "MEAN" or metric.upper() == "AVERAGE":
            return profile.mean(axis=0)
        elif metric.upper() == "SUM":
            return profile.sum(axis=0)
        else:
            raise NoSuchAMetricException("%s is not valid metric. Available are sum or mean.")

    def __repr__(self):
        return "Profile of %s on %s" % (self.signal.name, self.window.name)
