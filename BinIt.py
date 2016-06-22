import sys
import pyfaidx
from collections import OrderedDict, deque, defaultdict
from argparse import ArgumentParser
from GenoPy.DepthProfile import DepthProfile
from General.IO import IO
from GenoPy.Genomic import *
from GenoPy.Bins import *

parser = ArgumentParser(description='Just bin it, bin it ... Please pipe in samtools mpileup')
parser.add_argument('-b', required=True, dest='inputBed', help="Please give an input bed (with no header)")
parser.add_argument('--min-bin-size', required=False, default=10, dest='min', help="Minimal Bin Size")
parser.add_argument('--max-bin-size', required=False, default=50, dest='max', help="Maximal Bin Size")
parser.add_argument('-o', dest="output", default=None, help='If you want results written to a file, provide a file name')
args = parser.parse_args()
if args.output is None:
    output = sys.stdout
else: output = open(args.output, "w")

class MPileup(object):
    def __init__(self, io, parameters, bed):
        self.io = parameters['io']
        self.parameters = parameters
        self.bed = bed
	self.log = self.io.log

    def parse(self):
        alreadyPrinted = False
        #lines = sys.stdin.read()
        lineCounter = 0
	avancement = []
        self.dph = DepthProfile(sys.stdout, self.io)
        for line in sys.stdin:
            #print line
            lineCounter +=1 
            elem = line.split('\t')
            chr = elem[0].strip().rstrip('\n').rstrip('\t')
            pos = int(elem[1].strip().rstrip('\n').rstrip('\t'))
            refNuc = elem[2].strip().rstrip('\n').rstrip('\t')
            depth = int(elem[3].strip().rstrip('\n').rstrip('\t'))
            if not chr in avancement:
                self.dph.walk(Chromosome=chr, keepInMemory=True, binDynamically=parameters['bin'])
                avancement.append(chr)
                before = elem[1]
                self.log("Parsing and calling on {0}".format(chr))
            elif not lineCounter % 200000:
                self.dph.walk(Chromosome=chr, keepInMemory=True, binDynamically=parameters['bin'])
            if (lineCounter % 10000 == 0):
                self.log("{0} line read".format(lineCounter), loglevel=3)
            if (int(depth) > 0):
                self.dph.add(chr, pos, Position(chr, pos, '!', '!', depth, refNuc, pos_before=before, bed=self.bed))
                before = pos
        self.dph.writeDepthProfile(keepInMemory=True, binDynamically=parameters['bin'])

io = IO(3)
parameters = defaultdict(bool)
parameters['minBinSize'] = int(args.min)
parameters['maxBinSize'] = int(args.max)
parameters['io'] = io
parameters['bed'] = args.inputBed
bedF = Bed()
bedF.parse(open(args.inputBed, "r"))

dph = DepthProfile(sys.stderr, io)
bins = Binning(output, parameters)

parameters['dph'] = dph
parameters['bin'] = bins

pileup = MPileup(io, parameters, bedF)
pileup.parse()



