#!/usr/bin/env python2
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

# This is the wrapper for parallel processing.
# Rationale: I cannot use threads since the GIL is out there, I cannot use multiprocessing module since I have a bunch of unpickleable stuff,
#            so the third easiest way is to wrap my program and launch it several time in a map/reduce way
#            This wrapper will have almost the same interface as p2v and will become p2v in the end.
#            It will take as input bam files or pileup files, or pileup.gz files and do what it takes to index them / tabix index them using htslib.
#            Then, using a queue and the BED file, we will create datastreams and launch forks (n limited by command line arg) of p2v that will only see
#            the aforementioned genomic data. 
#            Maybe I will use os.Popen() to spawn the new forks (must see if I can make it non blocking and track the spawned forks), using pipes to feed genomic data
#            Maybe I will use os.fork() to spawn new forks in which case I might have to use temporary files. 
#
#            In all cases, a temporary directory will be created, in which we will store all the children results. Once there is more than one children result,
#            I will begin merging the files. The question is: should I manipulate the final files (easier I guess) or make p2v communicate with this process using some kind of IPC
#            system and send all results to it. This would be cleaner I guess, but there might be issues regarding pickling (well, we'll see when we're there).


from operator import itemgetter, attrgetter
from argparse import ArgumentParser
from collections import defaultdict, OrderedDict
import gzip
import sys
from decimal import *
#from scipy.stats import fisher_exact as fisher
from array import array
import numpy as np
from time import clock
#from matplotlib import pyplot as plt
#from Bio.Statistics.lowess import lowess

from GenoPy.Genomic import *
from GenoPy.Variants import *
from General.IO import *
from GenoPy.MPileup import *
from GenoPy.InDels import *
from GenoPy.DepthProfile import DepthProfile
from GenoPy.Bins import *
from GenoPy.VCF import *

from subprocess import Popen, PIPE
from multiprocessing.dummy import Pool as ThreadPool
from tempfile import mkstemp, mkdtemp
import os
from collections import deque

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
parser.add_argument("--maxDeltaDeltaDepthToCallTrueIndel", "-F", dest="maxDeltaDelta", default=400, required=False, help="This is the maximum delta between two delta depth not to filter the event as a false positive [300]")
parser.add_argument("--bigIndelsFile", dest="bigIndelsFile", default="bigIndels", required=False, help="BigIndel tabulated file filename")
parser.add_argument("--outputPrefix", dest="prefix", default=None, required=False, help="This is a prefix that will be prepended to all default filenames. When using this option, you will get every output by default")
parser.add_argument("--toleranceToDetectFalseStopInBigIndels", "-T", dest="tolerance", default=0.15, type=float, required=False)
parser.add_argument("--toleranceFalsePositiveEnd", "-C", dest="toleranceFalsePositiveEnd", default=2, type=int, required=False)
parser.add_argument("--indelWindowLength", dest="indelWindowLength", default=50, type=int, required=False)
parser.add_argument('--disableQualityScoreCheck', dest="disableQualityScoreCheck", default=False, action="store_true", help="Disables the quality score check")
parser.add_argument('--ncpu', dest='ncpu', default=2, help="How much CPUs should we use ?")
parser.add_argument('--tmpDir', dest='tmpDir', default='/tmp/', help="Specify here the temporary directory p2v should use")
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
parameters['regions'] = os.path.abspath(args.regions)
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
parameters['DisableQualityScoreCheck'] = args.disableQualityScoreCheck
parameters['tmpDir'] = os.path.abspath(args.tmpDir)

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
    parameters['minDeltaDepth'] = 110
    #parameters['minDeltaDepth'] = 1000

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
io.log("Starting program. Step 1: Map")
io.register("mpileupCount", parameters['input'], "r")
getLinesH = io.getIO("mpileupCount")

lineNb = getLineCount(getLinesH)

parameters['lineNb'] = lineNb
print parameters

io.giveupIO("mpileupCount")

if parameters['regions'] is not None:
    io.register('regions', parameters['regions'], 'r')

    # Instanciation de Bed et parsing du fichier.
    bed = Bed()
    bedIn = io.getIO('regions')
    bed.parse(bedIn)
    io.giveupIO('regions')
    
    

samtools = which('samtools')
bgzip = which('bgzip')
tabix = which('tabix')

if samtools is None or bgzip is None or tabix is None:
    raise

def validate_input(input):
    '''This should take care of any file, would it be BAM, PILEUP or PILEUP.GZ'''
    # Do some stuff to ensure we end with mpileup.gz(.tbi) files
    # What kind of file do we have as input ?
    if input.endswith('.bam'):
        # Need to check for index, then pileup, then gzip, then tabix
        # In fact, we should also sort the input ...
        index = Popen([samtools, 'index', input])
        index.communicate()
        pileupFile = input.rstrip('.bam') + '.pileup'
        with open(pileupFile, "w") as ph:
            pileup = Popen([samtools, 'mpileup', '-f', parameters['refGenome'], '-A', '-s', '-O', '-B', input], stdout=ph)
            pileup.communicate()
        input = pileupFile
    if input.endswith('.mpileup') or input.endswith('.pileup'):
        # Need to gzip then tabix
        bzip = Popen([bgzip, '-f', input])
        bzip.communicate()
        input = input + ".gz"
    # Need to tabix
    tabx = Popen([tabix, '-s', '1', '-b', '2', '-e', '2', input])
    tabx.communicate()
    return input

def prepare_data_pool(input, bed):
    io.log('Generating data pool', loglevel=2)
    pool = []
    for reg in bed.getRegions():
        tabixQuery = '{}:{}-{}'.format(reg.chr, reg.start, reg.stop)
        # Make a temporary pileup file restricted to that region, and return a list of temporary files
        tmpDir = mkdtemp("." + tabixQuery, 'p2v', dir=parameters['tmpDir'])
        tmpFile = mkstemp('.pileup', 'p2v.', dir=tmpDir, text=True)
        tmpFileH = os.fdopen(tmpFile[0])
        #print ' '.join([tabix, input, tabixQuery, '>', tmpFile[1]])
        tabxQ = Popen([tabix, input, tabixQuery], stdout=tmpFileH)
        tabxQ.communicate()
        pool.append({"region": tabixQuery, "pileup": tmpFile[1], "tmpDir": tmpDir + "/"})
        tmpFileH.close()
    io.log('Finished initializing data pool with {} elements'.format(len(pool)), loglevel=2)
    return pool
    
def launch_legacy_p2v(data):
    # Wraps launching p2v into a function
    # 1) Create a temporary subdir for results, which name will show what we're working on
    tmpDir = data['tmpDir']
    region = data['region']
    pileup = data['pileup']
    
    oPrefix = tmpDir + region.split(':')[-1]
    
    cmdLine = ['python2', 'p2v'] + sys.argv[1:] + ['--overrideInputPP2V', pileup, '--overrideOutputPP2V', oPrefix, '--silent']
    p2v = Popen(cmdLine)
    p2v.communicate()
    return {'oPrefix': oPrefix, 'region': region, 'tmpDir': tmpDir, 'pileup': pileup}

if not os.path.isdir(parameters['tmpDir']):
    os.makedirs(parameters['tmpDir'])
if not parameters['tmpDir'].endswith('/'): parameters['tmpDir'] += "/"

input = validate_input(parameters['input'])
data_pool = prepare_data_pool(input, bed)

# Thread Pool initialization so we can use map
# I use a thread pool here because I do not think I need any process complications ... Yeah I really have been shocked by this awful pickling stuff ...
io.log('Initializing Thread/Process pool with {} cpu(s)'.format(args.ncpu))
pool = ThreadPool(int(args.ncpu))
results = pool.map(launch_legacy_p2v, data_pool)
pool.close()
pool.join()

io.log("Finished step 1 (Map)")
io.log("Starting step 2 (Reduce)")

writeQueue = deque()
class Reduce(object):
    def __init__(self, output):
        self.writeQueue = deque()
        self.output = output
    def reduce(self, file, writeHeader=False):
        with open(file, "r") as handle:
            for line in handle:
                if writeHeader and line.startswith("#"):
                    self.writeQueue.append(line)
                if not line.startswith('#'):
                    self.writeQueue.append(line)
    def write(self):
        with open(self.output, "w") as handle:
            for line in self.writeQueue:
                handle.write(line)
                
writeHeader = True

GoodVariantReduce = Reduce(parameters['output'])
TrashReduce = Reduce(parameters['outputTrash'])

for el in results:
    oPrefix = el['oPrefix']
    output = oPrefix + ".GoodVariants.vcf"
    outputTrash = oPrefix + ".Trash.vcf"
    GoodVariantReduce.reduce(output, writeHeader=writeHeader)
    TrashReduce.reduce(outputTrash, writeHeader=writeHeader)
    writeHeader = False

GoodVariantReduce.write()
TrashReduce.write()

io.log("Finished step 2 (Reduce)")