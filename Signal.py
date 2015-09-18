import os
from sys import stderr
import HTSeq
from .utils import infer_extension

class WrongFiletypeException(Exception): pass
class NotSuchAFileException(Exception): pass
class NoSignalsException(Exception): pass


def get_signals(signals_list, number_of_processes=10):
    """Get Signals based on the provided list. Use multiprocessing.

    :param signals_list: list of dicts -- a list with dictionaries of signals to take into account [{'argname': arg_value}]
    :param number_of_processes: int -- number of parallel processes
    :returns: list of Signals

    Example:
        >>> from MetaPy import Signal
        >>> signals_list = ({'name': "input1",
        >>>                  'filepath': "/path/to/input1.bed"},
        >>>                 {'name': "input2",
        >>>                  'filepath': "/path/to/input2.bed"})
        >>> signals = Signal.get_signals(signals_list)

    """
    pool = multiprocessing.Pool(3)
    results = pool.map(helper_func, job_kwargs)
    pool.close()
    pool.join()
    return results

class Signal(object):

    """Genomic signal representation as HTSeq.GenomicArray"""

    def __init__(self, name, filepath, filetype=None, stranded=True, func=None):
        """
        Create genomic Signal

        :param name: str -- name of the signal
        :param filepath: str -- path to file with signals
        :param filetype: str -- type of the file (can be one of the following: BED, GTF, GFF, BAM, SAM or OTHER),
                                if the OTHER is specified a func parameter has to be provided
        :param stranded: bool -- specify whether the input is stranded
        :param func: function -- iterator function that takes filepath as argument and in each iteration it returns an
                                 object that has iv attribute that is an HTSeq.GenomicInterval
        """
        self.name = name
        if not os.path.isfile(filepath):
            raise NotSuchAFileException("No file named %s" % filepath)
        if filetype is None:
            filetype = infer_extension(filepath)
        if filetype.upper() not in ["BED", "GFF", "GTF", "BAM", "SAM", "OTHER"]:
            raise WrongFiletypeException("Wrong file type, allowed: %s" % "(BED, GTF, GFF, BAM, SAM, OTHER)")
        if filetype.upper() == "OTHER":
            if func is not None:
                raise ValueError("If the file type is OTHER you have to provide func parameter!")
            if not hasattr(func, '__call__'):
                raise TypeError("The func argument shall be a function!")
        self.filepath = filepath
        self.filetype = filetype
        self.__create_genomic_signals(stranded=stranded, func=func)

    @classmethod
    def get_signal(cls, name, filepath, filetype=None, stranded=True, func=None):
        """
        Create genomic Signal and return it. Can be used in threads and multiprocessing.

        :param name: str -- name of the signal
        :param filepath: str -- path to file with signals
        :param filetype: str -- type of the file (can be one of the following: BED, GTF, GFF, BAM, SAM or OTHER),
                                if the OTHER is specified a func parameter has to be provided
        :param stranded: bool -- specify whether the input is stranded
        :param func: function -- iterator function that takes filepath as argument and in each iteration it returns an
                                 object that has iv attribute that is an HTSeq.GenomicInterval
        """
        signal = cls(name=name, filepath=filepath, filetype=filetype, stranded=stranded, func=func)
        return signal

    def __create_genomic_signals(self, stranded=True, func=None):
        """Prepares coverage as a HTSeq.GenomicArray

        :param filepath: path to file
        :param filetype: type of the file (can be bed etc.)
        """
        stderr.write("Creating %s signal. It may take few minutes...\n" % self.name)
        self.coverage = HTSeq.GenomicArray("auto", stranded=stranded, typecode="d")
        self.library_size = 0
        if self.filetype.upper() == "BED":
            for line in HTSeq.BED_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif self.filetype.upper() == "GFF" or self.filetype.upper() == "GTF":
            for line in HTSeq.GFF_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif self.filetype.upper() == "SAM":
            for line in HTSeq.SAM_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif self.filetype.upper() == "BAM":
            for line in HTSeq.BAM_Reader(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        elif self.filetype.upper() == "OTHER":
            for line in func(self.filepath):
                self.coverage[line.iv] += 1
                self.library_size += 1
        else:
            assert False, "I should not be here!"

    def __repr__(self):
        return "<Signal: %s>" % self.name
