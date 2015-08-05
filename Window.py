import os
import HTSeq
from functions import infer_extension

class Window(object):

    """Representation of the windows"""

    def __init__(self, name, filepath, filetype=None, window_length=100, pseudocount=1):
        """Initialization of Window

        :param filepath: @todo
        :param filetype: @todo
        :param window_length: @todo
        :param pseudocount: @todo

        """
        self.name = name
        if not os.path.isfile(filepath):
            raise NotSuchAFileException("No file named %s" % filepath)
        if filetype is None:
            filetype = infer_extension(filepath)
        if filetype.upper() not in ["BED", "GFF", "GTF", "BAM", "SAM"]:
            raise WrongFiletypeException("Wrong file type, allowed: %s" % "(BED, GTF, GFF, BAM, SAM)")
        self.filepath = filepath
        self.filetype = filetype
        if type(window_length) != type(0):
            raise Exception("Window length has to be an integer")
        self.window_length = window_length
        if type(pseudocount) != type(0):
            raise Exception("Pseudocount has to be an integer")
        self.pseudocount = pseudocount
        self.__get_windows_iterator()

    def set_window_length(self, window_length):
        """Set window length
        """
        self.window_length = window_length

    def __get_windows_iterator(self):
        """Return an generator over intervals from the file including extension

        :param filepath: path to file
        :param filetype: type of the file (can be bed etc.)
        :param window_length: length of the window to add on both sides
        :returns: generator over expanded GenomicIntervals

        """
        filepath = self.filepath
        filetype = self.filetype
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
