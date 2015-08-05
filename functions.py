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

class NotSuchAMetricException(Exception): pass
class NoExtensionException(Exception): pass

def get_windows(filepath, filetype, window_length):
    """Return an generator over intervals from the file including extension

    :param filepath: path to file
    :param filetype: type of the file (can be bed etc.)
    :param window_length: length of the window to add on both sides
    :returns: generator over expanded GenomicIntervals

    """
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
    return Iterable()

def get_genomic_signals(path_to_file, filetype, stranded=True):
    """Prepares coverage as a HTSeq.GenomicArray

    :param filepath: path to file
    :param filetype: type of the file (can be bed etc.)
    :returns: 2-tuple of HTSeq.GenomicArray, library_size

    """
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


def make_profile_array(windows, signals, pseudocount=1):
    """Prepare array of profiles for each window. For each window a pseudocount
    is added.

    :param windows: iterator that returns a windows of equal length as HTSeq.GenomicInterval
    :param signals: HTSeq.GenomicArray of signals over chromosomes
    :param pseudocount: add a pseudocount to the window to avoid dividing by 0 in some operations,
                        defaults to 1
    :returns: numpy array of signals count for each window

    """
    profile = []
    window_lengths = set()
    for window in windows:
        window_lengths.add(window.length)
        assert len(window_lengths) == 1, "Windows has at least two different lengths"
        try:
            wincvg = np.fromiter(signals[window], dtype='i', count=window.length)
        except IndexError:
            sys.stderr.write("Wrong window: %s\n" % str(window))
            continue
        if window.strand == "+":
            profile.append(wincvg + pseudocount)
        else:
            profile.append(wincvg[::-1] + pseudocount)
    return np.asarray(profile)


def normalize_to_library(profile, library_size, metric='mean', aggregate=True):
    """Normalize profile to the library size

    :param profile: numpy array of windows, could be from make_profile_array
    :param library_size: size of the library used to generate signals
    :param metric: metric to use when aggregating profile, can be sum or mean (average),
                   defaults to mean
    :param aggregate: if True the data will be aggregated into 1-D vector according to
           metric, if False metric will be ignored, defaults to True
    :returns: @todo

    """
    if aggregate:
        if metric.upper() == "MEAN" or metric.upper() == "AVERAGE":
            return profile.mean(axis=0)/float(library_size)
        elif metric.upper() == "SUM":
            return profile.sum(axis=0)/float(library_size)
        else:
            raise NotSuchAMetricException("%s is not valid metric. Available are sum or mean.")
    else:
        return np.asarray(map(lambda i: i/float(library_size), profile))


def normalize_to_gene(profile, metric='mean', aggregate=True):
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
    adjusted_profile = np.asarray(map(lambda i: i/float(i.sum()), profile))
    if aggregate:
        if metric.upper() == "MEAN" or metric.upper() == "AVERAGE":
            return adjusted_profile.mean(axis=0)
        elif metric.upper() == "SUM":
            return adjusted_profile.sum(axis=0)
        else:
            raise NotSuchAMetricException("%s is not valid metric. Available are sum or mean.")
    else:
        return adjusted_profile


def create_nucleotide_frequency_from_fasta(path_to_fasta, sequence_length=None, verbose=False):
    """Create a pandas.DataFrame with the frequency of each nucleotide
    on each position. Sequences are converted to uppercase DNA and saved
    as IUPAC.ambiguous_dna, however only the results for A, C, T and G are
    returned.

    :param path_to_fasta: path to the fasta file
    :sequence_length: assumed length of the input sequences. If not provided length of the first
                      sequence will be chosen. It will be check and any sequece of different 
                      length will be ignored
    :returns: 2-tuple: pandas.DataFrame, Bio.motifs.Motif

    """
    seqs = []
    number_of_ignored = 0
    seq_counter = 0
    if not sequence_length:
        for rec in SeqIO.parse(path_to_fasta, 'fasta'):
            sequence_length = len(rec.seq)
            break
        sys.stderr.write("Presumed sequence length was not provided. First encounterd sequence length" +
                         " will be used i.e. %i\n" % sequence_length)
    for rec in SeqIO.parse(path_to_fasta, 'fasta'):
        seq_counter += 1
        if len(rec.seq) != sequence_length:
            if verbose:
                sys.stderr.write("%s has wrong length (%i)." % (rec.id, len(rec.seq)) +
                                 " It will be ignored.\n")
            number_of_ignored += 1
        else:
            seqs.append(Seq.Seq(str(rec.seq).upper().replace("U", "T"), IUPAC.ambiguous_dna))
    motifs_obj = motifs.create(seqs, IUPAC.ambiguous_dna)
    frequency_df = DataFrame(motifs_obj.pwm)[["A", "C", "T", "G"]]
    sys.stderr.write("%i sequences out of %i was ignored because of the length issue.\n" % (number_of_ignored, seq_counter))
    return frequency_df, motifs_obj


def plot_line_profile(profile, x=None, ax=None, smooth_by=1, **kwargs):
    if ax is None:
        fig, ax = pl.subplots()
    if profile.ndim == 1:
        if x is None:
            x = range(-(len(profile)//2), len(profile)//2 + len(profile)%2)
        assert len(profile) == len(x), "x vector and profile length are different"
        smoothed = pd.rolling_mean(profile, window=smooth_by, center=True)
        ax.plot(x, smoothed, **kwargs)
        return ax
    elif profile.ndim == 2:
        if x is None:
            x = range(-(profile.shape[1]//2), profile.shape[1]//2 + profile.shape[1]%2)
        assert profile.shape[1] == len(x), "x vector and profile length are different"
        smoothed = pd.rolling_mean(profile, window=smooth_by, center=True, axis=1)
        sns.tsplot(smoothed, ax=ax, **kwargs)
        return ax
    else:
        raise Exception("Input profile has wrong dimentionality")


def clean_axis(ax):
    """
    Remove ticks and tick labels

    Parameters
    ----------
    ax : matplotlib.axes
        Axes to modify
    """
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


def plot_heatmap(df, heatmapAX=None, xlabel='', ylabel='', title='', percentiles=(1, 99), sort_by=None, **imshow_kws):
    if heatmapAX is None:
        fig, heatmapAX = pl.subplots()
    vmin = np.percentile(df, percentiles[0])
    vmax = np.percentile(df, percentiles[1])

    # heatmap
    if sort_by is not None:
        assert df.shape[0] == sort_by.shape[0], "Sort values and profile has different shape"
        ind = sort_by.argsort()
    else:
        ind = np.arange(df.shape[0])
        assert df.shape[0] == len(ind)

    axi = heatmapAX.imshow(df[ind, :], interpolation='none',
            aspect='auto', origin='lower', vmin=vmin, vmax=vmax, **imshow_kws)
    heatmapAX.yaxis.set_label_position('right')
    heatmapAX.set_ylabel(ylabel)
    heatmapAX.set_xlabel(xlabel)
    clean_axis(heatmapAX)

    heatmapAX.set_title(title, fontsize=23)
    heatmapAX.set_ylabel(ylabel)
    heatmapAX.set_xlabel(xlabel)

def create_grid(figsize=(8, 12), strip=False, height_ratios=(1, 4), width_ratios=(1, 4), subplot_params=None):
    if subplot_params is None:
        subplot_params = {}
    fig = pl.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        len(height_ratios),
        len(width_ratios),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        wspace=0.1, hspace=0.1,
        **subplot_params)
    fig.array_axes = plt.subplot(gs[1, 0:])
    fig.line_axes = pl.subplot(gs[0, 0:], sharex=fig.array_axes)
    fig.gs = gs
    return fig

def infer_extension(path):
    """Infer extension of a file based on its path

    :param path: @todo
    :returns: @todo

    """
    filename, extension = os.path.splitext(path)
    if extension == '':
        raise NoExtensionException("File has no extension, please provide file type manually")
    else:
        return extension[1:]
