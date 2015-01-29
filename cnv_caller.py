import pysam
import numpy
import sys
import os

from optparse import OptionParser as opt
from subprocess import call

# Set the options that need to be set
prsr = opt()
prsr.add_option("-w", "--windowSize", dest="winsize", metavar="INT", default=500,
                help="Windowsize (bp) to be used to calculate log2ratio [Default:%default]")
prsr.add_option("-m", "--mappingQuality", dest="mapq", metavar="INT", default=0,
                help="Mapping quality cutoff for reads to be used in the calculation [Default:%default]")
prsr.add_option("-f", "--file", dest="bam", metavar="FILE",
                help="Input bam file to be analyzed, should be sorted and indexed")
prsr.add_option("-o", "--ouput", dest="path", metavar="PATH", default=os.getcwd(), help="Output path")
prsr.add_option("-l", "--name-list", dest="order", metavar="FILE",
                default=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "chr.list"),
                help="List of bam headers in order as they should be plotted, [Default:%default]")
prsr.add_option("-a", "--plot", dest="plot", metavar="BOOLEAN", default=True,
                help="Specify if plotting should be done using DNAcopy [Default:%default]")
prsr.add_option("-r", "--reference", dest="ref", metavar="FILE", help="Bam file to be used as refernce / control")

# Get options
(options, args) = prsr.parse_args()
DEVNULL = open(os.devnull, 'wb')


def checkFile(test_file):
    if not os.path.isfile(test_file):
        quit("Could not find the file:" + test_file)


def checkValidArgs(options):
    if options.bam == None:
        quit("ERROR: No BAM file submitted")
    checkFile(options.bam)


def setAbsPath(options):
    options.path = os.path.abspath(options.path)
    return options.path


def getNames(bam):
    header_dictSQ = bam.header.get('SQ')
    numberOfHeaders = len(header_dictSQ)
    #NOTE: Appending takes quite a while, if length is known, better create it and insert values
    NAMES = [None] * numberOfHeaders
    LENGTHS = [None] * numberOfHeaders
    i = 0
    while i < numberOfHeaders:
        NAMES[i] = header_dictSQ[i].get('SN')
        LENGTHS[i] = header_dictSQ[i].get('LN')
        i = i + 1
    return NAMES, LENGTHS


def median(valueList):
    #TODO: These numpy functions take any iterable and outputs an array, no need to convert it first
    return numpy.median(valueList)


def _pileUpBam(bam, name, ln):
    #TODO: Verify that the pile has a length
    bamPile = bam.pileup(name, 0, ln)
    result = [None] * len(bamPile)
    for i, pile in enumerate(bamPile):
        result[i] = pile.n
    return result

def getMedianCov(bam, NAMES, LENGTH):
    medians_dict = {}
    for name, ln in zip(NAMES, LENGTH):

        medians_dict[name] = median(_pileUpBam(bam, name, ln))

    return medians_dict


def getNormalizer(bam, ref, NAMES, LENGTH):
    REF = []
    BAM = []
    for name, ln in zip(NAMES, LENGTH):
        BAM += _pileUpBam(bam, name, ln)
        REF += _pileUpBam(ref, name, ln)

    normalizer = float(sum(REF)) / float(sum(BAM))
    return normalizer


class RWriter(object):
    def __init__(self, name_list):
        self.data = ()
        self.fh = None
        self.name = ()
        self.head = ()
        self.mid = ()
        self.merge = ()
        self.main = ()
        self.name_list = open(name_list, "r")
        self.counter = 0
        self.path = ()
        for name in self.name_list:
            self.counter = self.counter + 1
        self.name_list.seek(0)

    def getName(self):
        return self.name

    def setPath(self, path):
        self.path = path

    def writeHeader(self):
        HEAD = ["#!/usr/bin/env Rscript", "library(DNAcopy)"]
        self.head = '\n'.join(HEAD) + '\n'

    def writeMid(self):
        MID = []
        n = 1
        for name in self.name_list:
            name = name.rstrip()
            MID.append("c%s <- read.table(\"%s\")" % (n, os.path.join(self.path, name)))
            n = n + 1
        n = 1
        self.name_list.seek(0)
        for name in self.name_list:
            name = name.rstrip()
            MID.append("c%s$chr <- (\"%s\")" % (n, name))
            n = n + 1
        self.mid = '\n'.join(MID) + '\n'
        self.name_list.seek(0)

    def writeMerger(self):
        MERGE = []
        MERGE.append("merged <- rbind(")
        n = 1
        for name in self.name_list:
            name = name.rstrip()
            while n <= self.counter:
                if n == 1:
                    MERGE.append("c")
                    MERGE.append(str(n))
                else:
                    MERGE.append(",c")
                    MERGE.append(str(n))
                n = n + 1
        MERGE.append(")")
        MERGE.append("\npdf(\"%s\")" % os.path.join(self.path, "cnv_report.pdf"))
        self.merge = ''.join(MERGE) + '\n'
        self.name_list.seek(0)

    def writeMain(self):
        MAIN = []
        MAIN.append("CNA.object <- CNA(merged$V1, merged$chr, merged$V2, data.type=(\"logratio\"), presorted=TRUE)")
        MAIN.append("CNA.smooth <- smooth.CNA(CNA.object)")
        MAIN.append("CNA.segm <- segment(CNA.smooth)")
        MAIN.append("plot(CNA.segm, plot.type=\"w\")")
        MAIN.append("plot(CNA.segm, plot.type=\"s\")")
        MAIN.append("plot(CNA.segm, plot.type=\"p\")")
        n = 1
        for name in self.name_list:
            MAIN.append(
                "c%s.object <- CNA(c%s$V1, c%s$chr, c%s$V2, data.type=(\"logratio\"), presorted=TRUE)" % (n, n, n, n))
            MAIN.append("c%s.smooth <- smooth.CNA(c%s.object)" % (n, n))
            MAIN.append("c%s.segm <- segment(c%s.smooth)" % (n, n))
            MAIN.append("plot(c%s.segm, plot.type=\"s\")" % n)
            n = n + 1
        MAIN.append("dev.off()")
        self.main = '\n'.join(MAIN)

    def setName(self, name):
        self.name = name

    def assembler(self):
        self.writeHeader()
        self.writeMid()
        self.writeMerger()
        self.writeMain()
        self.writer(self.name)

    def writer(self, name):
        self.fh = open(name, "w")
        self.fh.write(self.head)
        self.fh.write(self.mid)
        self.fh.write(self.merge)
        self.fh.write(self.main)
        self.fh.close()


class CovScanner(object):
    def __init__(self):
        self.window_size = -1
        self.pos_start = 0
        self.pos_end = -1
        self.bamfile = None
        self.ref = None
        self.name = ""
        self.mapq_cutoff = -1
        self.RATIOS = []
        self.POS = []
        self.medians = None
        self.normalizer = None

    def setWindowSize(self, win):
        self.window_size = win
        self.pos_end = win

    def setChrMedians(self, medians):
        self.medians = medians

    def setSamFile(self, bam, ref):
        self.bamfile = bam
        self.ref = ref

    def setNormalizer(self, normalizer):
        self.normalizer = normalizer

    def setMapq(self, mapq):
        self.mapq_cutoff = int(mapq)

    def setName(self, name):
        self.name = name

    def move(self, ln):
        while self.pos_end < ln:
            window_mean = self.windowMean(self.bamfile)
            self.POS.append(self.pos_end)
            self.RATIOS.append(self.getLogRatios(window_mean))
            self.pos_start = self.pos_start + self.window_size
            self.pos_end = self.pos_end + self.window_size
        return self.RATIOS, self.POS

    def getLogRatios(self, window_mean):
        if self.ref is None:
            logRatio = window_mean / self.medians[self.name]
        else:
            win_mean_ref = self.windowMean(self.ref)
            if win_mean_ref == 0:
                win_mean_ref = 1
            logRatio = (window_mean * self.normalizer) / win_mean_ref

        if logRatio != 0:
            logRatio = numpy.log2(logRatio)

        return logRatio

    def windowMean(self, bam):
        bam = pysam.AlignmentFile(bam, "rb")
        bamPile = bam.pileup(self.name, self.pos_start, self.pos_end, truncate=True)
        WINDOW_COVS = numpy.zeros((len(bamPile),))

        for i, pile in enumerate(bamPile):
            cov = 0
            for reads in pile.pileups:
                if reads.alignment.mapq >= self.mapq_cutoff:
                    cov = cov + 1
            WINDOW_COVS[i] = cov

        bam.close()

        if WINDOW_COVS.sum() > 0:
            return numpy.mean(WINDOW_COVS)
        else:
            return 0


class FilePrinter(object):
    def __init__(self, path):
        self._path = path
        self._fh = None

    def __enter__(self):

        self._fh = open(self._path, "w")
        return self

    def __exit__(self, *args):

        self._fh.close()
        self._fh = None

    def printToFile(self, *data):

        if (self._fh is None):
            raise Exception("Only use within with statement blocks")

        line = '\t'.join(map(str, data)) + "\n"
        self._fh.write(line)

    def printListToFile(self, LIST):

        #TODO: If line always is a string, this is OK else it would be correct to unpack it so that each item of
        #the line is sent separately using
        #self.printToFile(*line)
        for line in LIST:
            self.printToFile(line)

    def printZipListToFile(self, *listOfLists):

        """Zips all parameter lists and prints them

        :param listOfLists: All lists sent to function to be printed
        """
        for line in zip(*listOfLists):
            self.printToFile(*line)


""" START """

if __name__ == "__main__":

    options.path = setAbsPath(options)
    normalizer = 0

    if not os.path.exists(options.path):
        os.makedirs(options.path)

    bam = pysam.AlignmentFile(options.bam, "rb")
    NAMES, LENGTH = getNames(bam)

    if bool(options.plot) == True:
        rscripter = RWriter(options.order)
        rscripter.setName(os.path.join(options.path, "run_DNAcopy.r"))
        rscripter.setPath(options.path)
        rscripter.assembler()

    if options.ref == None:
        chr_medians = getMedianCov(bam, NAMES, LENGTH)
    else:
        ref = pysam.AlignmentFile(options.ref, "rb")
        normalizer = getNormalizer(bam, ref, NAMES, LENGTH)

    for name, ln in zip(NAMES, LENGTH):
        with FilePrinter(os.path.join(options.path, name)) as out:
            scan = CovScanner()
            if options.ref == None:
                scan.setChrMedians(chr_medians)
            scan.setSamFile(options.bam, options.ref)
            scan.setMapq(options.mapq)
            scan.setWindowSize(options.winsize)
            scan.setNormalizer(normalizer)
            scan.setName(name)
            WINDOW_RATIOS, POS = scan.move(ln)
            out.printZipListToFile(WINDOW_RATIOS, POS)

    if bool(options.plot) == True:
        cmd = ["Rscript", rscripter.getName()]
        call(cmd, stdout=DEVNULL, stderr=DEVNULL)
        rm_cmd = ["rm", rscripter.getName()]
        call(rm_cmd)

    bam.close()
    DEVNULL.close()
