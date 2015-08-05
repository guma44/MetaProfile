import os
from sys import stderr
import HTSeq
from functions import infer_extension

class WrongFiletypeException(Exception): pass
class NotSuchAFileException(Exception): pass
class NoSignalsException(Exception): pass

class Signal(object):

    """Genomic signal representation as HTSeq.GenomicArray"""

    def __init__(self, name, filepath, filetype=None, stranded=True):
        self.name = name
        if not os.path.isfile(filepath):
            raise NotSuchAFileException("No file named %s" % filepath)
        if filetype is None:
            filetype = infer_extension(filepath)
        if filetype.upper() not in ["BED", "GFF", "GTF", "BAM", "SAM"]:
            raise WrongFiletypeException("Wrong file type, allowed: %s" % "(BED, GTF, GFF, BAM, SAM)")
        self.filepath = filepath
        self.filetype = filetype
        self.__create_genomic_signals(stranded=stranded)

    def __create_genomic_signals(self, stranded=True):
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
