#!/usr/bin/env python
# encoding: utf-8
#
### CHANGELOG ###
# 1.2: Add support for InDel detection method.
#      - The method used is depth analysis: by segmenting genome's depth, we can spot InDels.
#      - In order to do that, we launch SegmentIt, and gather the segments.
#      - These segments will point to the InDels.
# 1.1: Added flag for multiple indels at position
# 1.1: Added possibility for non-ref variants (ex: ref=A, variant=C/G)
#
### LICENCE: GNU GPLv3

from operator import itemgetter, attrgetter
from argparse import ArgumentParser
from collections import defaultdict, OrderedDict
import gzip
import sys
from decimal import *
from scipy.stats import fisher_exact as fisher
from array import array
import numpy as np
from time import clock
from matplotlib import pyplot as plt
from Bio.Statistics.lowess import lowess

from GenoPy.Genomic import *
from GenoPy.Variants import *
from General.IO import *
from GenoPy.MPileup import *
from GenoPy.InDels import *
from GenoPy.DepthProfile import DepthProfile
from GenoPy.Bins import *
from GenoPy.VCF import *

DEBUG=True
VERSION="1.1.3"
DEVELOPMENT=True
UNSTABLE=True


###############################################################################
###############################  Main   #######################################
###############################################################################

print >>sys.stdout, "Command line:", sys.argv
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input", type=str, required=True, help="Please input the VCF file outputed by IonTorrent VariantCaller")
parser.add_argument("--minFreq", dest="minFreq", type=float, required=False, default=0.01, help="Please input the minimum frequency of the variants you want reported. [0.01]")
parser.add_argument("--minDepth", dest="minDepth", type=int, required=False, default=100, help="Please input minimum depth for variant calling. [100]")
parser.add_argument("--minReadsAlt", dest="minReadsAlt", type=int, required=False, default=5, help="Minimum reads supporting ALT to call a variant [5]")
parser.add_argument("--strandBias", dest="strandBias", type=int, required=False, default=90, help="This is supposed to be the threshold for maximum strand bias: if more than <value> percents reads are on one strand, the variant will be discarded. [90]")
parser.add_argument("--fisherBias", dest="fisherBias", type=float, required=False, default=60, help="This is another indicator for strand bias. A Fisher Exact Test is computer for reference and alternate contigency tables. The p-value is then phred scaled. High values indicate more bias. [60]")
parser.add_argument("--outputTrash", dest="outputTrash", default="Trash", type=str, required=False, help="This option allows you to output trashed variants. Please give either a path to a file, or a '-' to output on stderr")
parser.add_argument("--quality", dest="quality", type=int, required=False, default=10, help="This allows to filter on quality. Quality shows the distance between called allele and other alleles at this position [10]")
parser.add_argument("--output", "-o", dest="output", type=str, required=False, default="GoodVariants", help="This option allows you to output called variants. Please give either a path to a file, or a '-' to output on stdout")
parser.add_argument("--statistics", "-s", dest="statistics", type=str, required=False, default="Statistics")
parser.add_argument("--localStrandBias", dest="localStrandBias", required=False, default=100, type=int, help="Local strand bias. This checks for desequilibrium between forward and reverse strand in alternate allele [100 (deactivated)]")
parser.add_argument("--hetTreshold", dest="hetTreshold", required=False, default=20, type=int, help="Treshold on variant frequency to call a heterozygous genotype [20]")
parser.add_argument("--homTreshold", dest="homTreshold", required=False, default=75, type=int, help="Treshold on variant frequency to call a homozygous genotype [75]")
#parser.add_argument("--snpCluster", dest="snpCluster", required=False, default=3, type=int, help="How many consecutive variants to filter out as snpCluster ? [3]")
#parser.add_argument("--snpClusterWindowSize", dest="snpClusterWindowSize", required=False, default=20, type=int, help="Size of the window used for snpCluster filter computation [20bp]")
parser.add_argument("--minBinSize", dest="minBinSize", required=False, default=5, help="Minimum size of any Bin to be called. A bin will not be called if there is less than this number of positions. Put 0 to deactivate. [5]")
parser.add_argument("--maxBinSize", dest="maxBinSize", required=False, default=50, help="Maximum size of any Bin to be called. A new bin will be created each time this size is reached. [50]")
parser.add_argument("--outputBins", dest="outputBins", required=False, default="Bins", help="Bin file filename")
parser.add_argument("--depthProfile", dest="depthProfile", required=False, default="DepthProfile", help="DepthProfile file filename")
parser.add_argument("--regions", "-r", dest="regions", required=True, help="Please give as input a BED file describing the regions you sequenced")
parser.add_argument("--loglevel", dest="ll", required=False, default=1, help="Loglevel: 1: normally verbose; 5: debug mode [1]")
parser.add_argument("--refGenome", "-R", dest="refGenome", required=True, help="Please give as input the Fasta file of the genome you used to align your sequences")
parser.add_argument("--nofilter", action='store_true', dest="nofilter", required=False, help="Disables all filters. Please note that this cannot be overwritten by parameters you would enter.")
parser.add_argument("--clinics", action='store_true', dest="clinics", required=False, help="Presets all the filters for clinical uses. Please not that this cannot be overwritten by parameters you would enter.")
parser.add_argument("--minDeltaDepthToCallIndels", "-D", dest="minDeltaDepth", default=200, required=False, help="This is the minimum delta depth between 2 positions used in order to call for an InDel event [200]")
parser.add_argument("--maxDeltaDeltaDepthToCallTrueIndel", "-F", dest="maxDeltaDelta", default=300, required=False, help="This is the maximum delta between two delta depth not to filter the event as a false positive [300]")
parser.add_argument("--bigIndelsFile", dest="bigIndelsFile", default="bigIndels", required=False, help="BigIndel tabulated file filename")
parser.add_argument("--outputPrefix", dest="prefix", default=None, required=False, help="This is a prefix that will be prepended to all default filenames. When using this option, you will get every output by default")
parser.add_argument("--toleranceToDetectFalseStopInBigIndels", "-T", dest="tolerance", default=0.15, type=float, required=False)
parser.add_argument("--toleranceFalsePositiveEnd", "-C", dest="toleranceFalsePositiveEnd", default=2, type=int, required=False)
parser.add_argument("--indelWindowLength", dest="indelWindowLength", default=50, type=int, required=False)
# Parametres
args = parser.parse_args()
parameters = defaultdict(bool)
parameters['input'] = args.input
parameters['minFreq'] = args.minFreq
parameters['minDepth'] = args.minDepth
parameters['strandBias'] = args.strandBias
parameters['fisherStrand'] = args.fisherBias
parameters['outputTrash'] = args.outputTrash
parameters['minReadsAlt'] = args.minReadsAlt
parameters['minQuality'] = args.quality
parameters['output'] = args.output
parameters['statistics'] = args.statistics
parameters['hetTreshold'] = args.hetTreshold
parameters['homTreshold'] = args.homTreshold
parameters['localStrandBias'] = args.localStrandBias
parameters['minBinSize'] = int(args.minBinSize)
parameters['maxBinSize'] = int(args.maxBinSize)
parameters['outputBins'] = args.outputBins
parameters['depthProfile'] = args.depthProfile
parameters['regions'] = args.regions
parameters['sName'] = args.input.rstrip('.mpileup')
parameters['onlyBed'] = True
parameters['ll'] = args.ll
parameters['toleranceForBigIndels'] = args.tolerance
parameters['refGenome'] = args.refGenome
parameters['minDeltaDepth'] = args.minDeltaDepth
parameters['maxDeltaDelta'] = int(args.maxDeltaDelta)
parameters['bigIndelsFile'] = args.bigIndelsFile
parameters['prefix'] = args.prefix
parameters['toleranceFalsePositiveEnd'] = args.toleranceFalsePositiveEnd
parameters['IndelWindowLength'] = args.indelWindowLength

if args.nofilter is True:
    parameters['minFreq'] = 0
    parameters['minDepth'] = 0
    parameters['strandBias'] = 0
    parameters['fisherStrand'] = 10000000
    parameters['minReadsAlt'] = 0
    parameters['minQuality'] = 0
    parameters['localStrandBias'] = 0

if args.clinics is True:
    parameters['minFreq'] = 15
    parameters['minDepth'] = 20
    parameters['strandBias'] = 5
    parameters['fisherStrand'] = 250
    parameters['minReadsAlt'] = 20
    parameters['minQuality'] = 50
    parameters['localStrandBias'] = 90
    parameters['hetTreshold'] = 20
    parameters['homTreshold'] = 70
    parameters['toleranceForBigIndels'] = 0.20
    parameters['toleranceFalsePositiveEnd'] = 3

if parameters['prefix'] is not None:
    parameters['output'] = parameters['prefix'] + "." + parameters['output'] + ".vcf"
    parameters['statistics'] = parameters['prefix'] + "." + parameters['statistics']
    parameters['outputBins'] = parameters['prefix'] + "." + parameters['outputBins']
    parameters['depthProfile'] = parameters['prefix'] + "." + parameters['depthProfile']
    parameters['bigIndelsFile'] = parameters['prefix'] + "." + parameters['bigIndelsFile']
    parameters['outputTrash'] = parameters['prefix'] + "." + parameters['outputTrash'] + ".vcf"
    parameters['regionsw'] = parameters['prefix'] + ".bed"

#global io

io = IO(parameters['ll'])
parameters['io'] = io
io.register("mpileupCount", parameters['input'], "r")
getLinesH = io.getIO("mpileupCount")

lineNb = getLineCount(getLinesH)

parameters['lineNb'] = lineNb
print parameters

io.giveupIO("mpileupCount")

io.register("mpileup", parameters['input'], "r")
io.register("binningW", parameters['outputBins'], "w")
io.register("depthProfile", parameters['depthProfile'], "w")
#io.register("bigIndels", parameters['bigIndelsFile'], "w")
io.register("statistics", parameters['statistics'], "w")

if (((parameters['output'] == '-') or (parameters['output'] is None)) and parameters['prefix'] is not None):
    output = sys.stdout
else:
    output = parameters['output']

io.register("GoodVariants", output, "w")
if parameters['outputTrash'] is not None:
    io.register("TrashedVariants", parameters['outputTrash'], 'w')
    #trashVariants = TrashCollection()
else:
    trashVariants = None
if parameters['regions'] is not None:
    io.register('regions', parameters['regions'], 'r')
    io.register("bedW", parameters['regionsw'], "w")

    # Instanciation de Bed et parsing du fichier.
    bed = Bed()
    bedIn = io.getIO('regions')
    bed.parse(bedIn)
    io.giveupIO('regions')

# Instanciation de MPileup et parsing du fichier.
pileuph = io.getIO('mpileup')
#goodVariants = GoodCollection()
depthProfile = DepthProfile(io.getIO('depthProfile'), io)
pileup = MPileup(pileuph, depthProfile, io, parameters, bed=bed)

unclassified_variants = pileup.parse()
good, trashed = unclassified_variants.classify_variants().getCollections()
AllCollections = AllClassifiedVariantCollections(parameters)
AllCollections.add('Good', good)
AllCollections.add('Trash', trashed)


io.giveupIO('mpileup')

io.register("depthProfileR1", parameters['depthProfile'], "r")
io.register("depthProfileR2", parameters['depthProfile'], "r")
indels = scanForIndels(parameters)
allIndels = InDelDetector(indels)
allIndels = allIndels.extract().genotype(parameters).getIndels()
#io.register("weirdRegionW", parameters['bigIndelsFile'], "w")
readDepthProfile = io.getIO('depthProfileR1')
binningW = io.getIO('binningW')
#weirdRegionH = io.getIO('weirdRegionW')
binning = Binning(binningW, parameters)
binning.bin(readDepthProfile)

#pileup.addIndelsToVCF()
io.giveupIO('binningW')
io.giveupIO('depthProfileR1')
io.giveupIO('depthProfileR2')

# Output des résultats.
# Binning
#io.register('binningR', parameters['outputBins'], "r")
#binningR = io.getIO('binningR')
# # Regular: show all features
#pileup.printBinning(binningR)
# # Bed: show each region in plot.
#pileup.alternateBinningPrinting(binningR)

goodVCF = VCF(good, parameters, indelsToAdd=allIndels)
goodVCF.printVCF(io.getIO('GoodVariants'), printHeader=True)
#pileup.printVCF(goodVariants)
io.giveupIO('GoodVariants')
if parameters['outputTrash'] is not None:
    trashVCF = VCF(trashed, parameters, indelsToAdd=allIndels)
    trashVCF.printVCF(io.getIO('TrashedVariants'), printHeader=True)
    io.giveupIO('TrashedVariants')

if parameters['statistics']:
    stath = io.getIO('statistics')
    AllCollections.do_stats(stath)
    io.giveupIO('statistics')

bedW = io.getIO('bedW')
bed.write(bedW)
io.giveupIO('bedW')
io.endlog()
sys.exit(io.unregister("all"))
