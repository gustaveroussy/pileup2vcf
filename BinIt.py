#!/usr/bin/env python2
# encoding: utf-8

#usage : cat My.mpileup | python BinIt.py -b My.bed --min-bin-size 0 --max-bin-size 50 -o My.ouput
#usage : samtools mpileup -A -s -O -B -d 50000 -f <reference genome> <your bam file> | python BinIt.py -b My.bed --min-bin-size 0 --max-bin-size 50 -o My.ouput

### BED file must be sort and without overlapping regions (see below to remove overlapping regions
###cat /data/bioinfo/hg19/captureKitDesigns/Sureselect_v5.padded.bed | awk '{if(last_chr==$1){ if(last_end<$2){print last_chr"\t"last_beg"\t"last_end; last_chr=$1; last_beg=$2; last_end=$3} else{last_end=$3}}else{if(last_chr!=""){print last_chr"\t"last_beg"\t"last_end;} last_chr=$1; last_beg=$2; last_end=$3}} END{print last_chr"\t"last_beg"\t"last_end;}' > Sureselect_v5.padded_modif.bed
###


import sys
import pyfaidx
from collections import OrderedDict, deque, defaultdict
from argparse import ArgumentParser
from GenoPy.DepthProfile import DepthProfile
from General.IO import IO
from GenoPy.Genomic import *
from GenoPy.Bins import *
import numpy

parser = ArgumentParser(description='Just bin it, bin it ... Please pipe in samtools mpileup')
parser.add_argument('-b', required=True, dest='inputBed', help="Please give an input bed (with no header)")
parser.add_argument('--min-bin-size', required=False, default=10, dest='min', help="Minimal Bin Size")
parser.add_argument('--max-bin-size', required=False, default=50, dest='max', help="Maximal Bin Size")
parser.add_argument('-o', dest="output", default="output_binning", help='If you want results written to a file, provide a file name (Default output_binning)')
args = parser.parse_args()
if args.output is None:
    output = sys.stdout
else: output = open(args.output, "w")


class Bed(object):
    '''This object is intented to store the genomic coordinates of a Bed file. This will be used, 
    for instance, in Binning and, later in Indel calling'''
    def __init__(self):
        self.capture = defaultdict(list)
        self.cache = defaultdict(list)
        self.cached = None
        self.order = []
        
    def parse(self, fh):
        '''Parses BED file (without structural verification, this will be fixed later) and
        stores Genomic objects in a dictionnary-like structure (in order to narrow the search).
        '''
        for line in fh:
            if not line.startswith('#') and not line.startswith('track'):
                sLine = line.split('\t')
                chr, start, stop = sLine[0], int(sLine[1]), int(sLine[2])
                self.capture[chr].append(Genomic(chr, start, stop, commentary='\t'.join(sLine[3:]).rstrip('\n')))
                self.cache[chr].append(Genomic(chr, start, stop, commentary=None))
                self.order.append(Genomic(chr, start, stop))
                
        # Now we sort the positions in bed, using "start" values as indexes
        # For now, we do not care about overlapping regions. Each time we
        # return a value, it will be cached to be returned immediatly if request
        # Each interval will have an expiration position (stop position)
        # If we get to that position, we elect a new interval (the next one) and check that
        # the current position is inside it. If it is not, we will check that the
        # stop nucleotide hasn't already expired. If it has, we elect a new interval until we get
        # one for which it hasn't expired. If it hasn't, this becomes the new cached interval.

        for k in self.cache.keys():
            self.cache[k] = sorted(self.cache[k], key=attrgetter('start'), reverse=True)
    
    def getRegions(self):
        for ch in self.cache.keys():
            for reg in self.cache[ch]:
                yield reg.chr, reg.start, reg.stop
    def getRegionsList(self):
        regs = []
        for ch in self.cache.keys():
            for reg in self.cache[ch]:
                regs.append((reg.chr, reg.start, reg.stop))
        return regs
    def getOrderRegionsList(self):
        return self.order

    def getBedAppartenance(self, o):
        '''Checks whether a Genomic-like object is related to a region described in the bed objects.
        Input: <Genomic-like>
        Output: <str> id and <bool>True
        Output: <None> and <bool>False

        NB: Maintaining two different output shoudln't be really necessary, but I am waiting for the next release'''
        chr, start, stop = str(o.chr), o.start, o.stop
        #print self.cached if self.cache is not None
        if self.cached is None:
            if self.cache[chr] != []:
                self.cached = self.cache[chr].pop()
                #print self.cache[chr]
            else:
                self.cached = None
                return None, False
        #print self.cached
        if o in self.cached:
            return self.cached.id, True
        else:
            if (self.cache[chr] is False): return None, False
            while start > self.cached.stop or chr != self.cached.chr:
                # Elect new start
                if (len(self.cache[chr]) == 1):
                    self.cached = self.cache[chr].pop()
                    self.cache[chr] = False
                    break
                if (len(self.cache[chr]) == 0):
                    #print 'Why did you request any update ? Were still not in cache'
                    return None, False
                else:
                    self.cached = self.cache[chr].pop()
            if o in self.cached:
                return self.cached.id, True
            else: return None, False


class MPileup(object):
    def __init__(self, io, parameters, bed):
        self.io = parameters['io']
        self.parameters = parameters
        self.bed = bed
        self.log = self.io.log
        self.dph = ""

    def parseAndBin(self, outBins):
        '''This takes as input the mpileup file, parses it and write bins with coverage.
        '''
        print >>outBins, "#chr\tbeg\tend\tnumberOfCoveragePoints\tmeanCoverage"        
        
        #get back regions for current chr
        chrBedRegion = []
        bedRegions = self.bed.getOrderRegionsList()
        #for genomic in bedRegions :
            #print genomic.chr, genomic.start, genomic.stop

        #print "Current"
        index_bedRegions = 0
        #print bedRegions[index_bedRegions].chr, bedRegions[index_bedRegions].start, bedRegions[index_bedRegions].stop

        #Coverage list, mean and size
        covList= []
        covListSize = 0
        covListMean = 0
        
        #to store how many bins we add on the current unfinished region
        delta = 0

        lineCounter = 0
        for line in sys.stdin:
            lineCounter += 1
            try:
                #with line.split('\t') as elem:
                elem = line.split('\t')
                chr = elem[0].strip().rstrip('\n').rstrip('\t')
                pos = int(elem[1].strip().rstrip('\n').rstrip('\t'))
                refNuc = elem[2].strip().rstrip('\n').rstrip('\t').upper()
                depth = int(elem[3].strip().rstrip('\n').rstrip('\t'))
                quality = "A"
                mapQual = "A"
                before = 1
                #sequence = elem[4].strip().rstrip('\n').rstrip('\t')
                #quality = elem[5].strip().rstrip('\n').rstrip('\t')
                #mapQual = elem[6].strip().rstrip('\n').rstrip('\t')
                #posInRead = elem[7].strip().rstrip('\n').rstrip('\t').split(',')
                #print "MPileup locations", chr, pos, depth
            except:
                print 'CrashDump:'
                print line
                print [len(x) for x in elem]
                raise MpileupFormatError("Cannot initialize at line {}.\n Please check that your Mpileup file has 8 columns and has been generated with mapping quality and base position in read using samtools 0.1.18 (further versions seems to suffer from bugs).\n MPileup with good format might be generated using: samtools mpileup -A -s -O -B -d 50000 -f ../Pileup2VCF/hg19/hg19.fa 208204422-ADN-2_S2_L001.bam".format(lineCounter))

            curPos = Position(chr, pos, quality, mapQual, depth, refNuc, pos_before=before, bed=self.bed)
            if bedRegions[index_bedRegions].chr == chr :
                #If before the current bed regions, nothing to care so continue
                if bedRegions[index_bedRegions].start > pos :
                    continue
                else:
                    #If we pass the current region
                    while bedRegions[index_bedRegions].chr == chr and bedRegions[index_bedRegions].start <= pos and bedRegions[index_bedRegions].stop < pos :
                        for i in range(bedRegions[index_bedRegions].start+delta, bedRegions[index_bedRegions].stop-parameters['maxBinSize']+1, parameters['maxBinSize']):
                            if i+parameters['maxBinSize']-1 < pos:
                                if len(covList) > 0 :
                                    covListSize = len(covList)
                                    covListMean = numpy.mean(covList)
                                else:
                                    covListSize = 0
                                    covListMean = 0
                                #print bedRegions[index_bedRegions].chr, i, i+parameters['maxBinSize']-1, 0, covListSize, covListMean, numpy.mean(covList), covList, "PASS 1"
                                print >>outBins, "{}\t{}\t{}\t{}\t{}\tPASS1".format(bedRegions[index_bedRegions].chr, i, i+parameters['maxBinSize']-1, covListSize, covListMean)
                                covList = []
                                #delta += parameters['maxBinSize']
                                boolPass1 = 1
                                boolPass2 = 0
                                boolPass3 = 0
                                
                        if parameters['minBinSize'] < bedRegions[index_bedRegions].stop-i+parameters['maxBinSize']-1:
                            if len(covList) > 0 :
                                covListSize = len(covList)
                                covListMean = numpy.mean(covList)
                            else:
                                covListSize = 0
                                covListMean = 0
                            if boolPass1 == 0:
                                #print bedRegions[index_bedRegions].chr, i, bedRegions[index_bedRegions].stop, 0, covListSize, covListMean, numpy.mean(covList), covList, "PASS 2"
                                print >>outBins, "{}\t{}\t{}\t{}\t{}\tPASS2.0".format(bedRegions[index_bedRegions].chr, i, bedRegions[index_bedRegions].stop, covListSize, covListMean)
                            else:
                                #print bedRegions[index_bedRegions].chr, i+parameters['maxBinSize']-1, bedRegions[index_bedRegions].stop, 0, covListSize, covListMean, numpy.mean(covList), covList, "PASS 2"
                                print >>outBins, "{}\t{}\t{}\t{}\t{}\tPASS2.1".format(bedRegions[index_bedRegions].chr, i+parameters['maxBinSize'], bedRegions[index_bedRegions].stop, covListSize, covListMean)
                            boolPass1 = 0
                            boolPass2 = 1
                            boolPass3 = 0
                        index_bedRegions+=1
                        delta = 0
                        if index_bedRegions == len(bedRegions):
                            sys.exit(io.unregister("all"))
                            sys.exit(0)
                        covList = []
                        if bedRegions[index_bedRegions].chr == chr and bedRegions[index_bedRegions].start <= pos and bedRegions[index_bedRegions].stop >= pos:
                            if curPos.inBed:
                                if bedRegions[index_bedRegions].stop >= pos and bedRegions[index_bedRegions].start+delta >= pos:
                                    covList.append(depth)
                    
                    #If we are in the current region
                    if bedRegions[index_bedRegions].stop >= pos :
                        #Perform bins before variant location
                        if bedRegions[index_bedRegions].start+delta == pos:
                            if curPos.inBed:
                                covList.append(depth) 
                                #print >>outBins, "PASS"
                        for i in range(bedRegions[index_bedRegions].start+delta, pos, parameters['maxBinSize']):
                            if i+parameters['maxBinSize']-1 < pos:
                                if len(covList) :
                                    covListSize = len(covList)
                                    covListMean = numpy.mean(covList)
                                else:
                                    covListSize = 0
                                    covListMean = 0
                                #print bedRegions[index_bedRegions].chr, i, i+parameters['maxBinSize']-1, 0, covListSize, covListMean, numpy.mean(covList), covList, "PASS 3"
                                print >>outBins, "{}\t{}\t{}\t{}\t{}\tPASS3".format(bedRegions[index_bedRegions].chr, i, i+parameters['maxBinSize']-1, covListSize, covListMean)
                                boolPass1 = 0
                                boolPass2 = 0
                                boolPass3 = 1
                                delta += parameters['maxBinSize']
                                covList = []
                                if curPos.inBed:
                                    if bedRegions[index_bedRegions].stop >= pos and bedRegions[index_bedRegions].start+delta >= pos:
                                        covList.append(depth)
                                        #print >>outBins, "PASS"
                            else:
                                if curPos.inBed:
                                    covList.append(depth)
                                    #print "add"

#get back IO                              
io = IO(3)

#get back parameters
parameters = defaultdict(bool)
parameters['minBinSize'] = int(args.min)
parameters['maxBinSize'] = int(args.max)
parameters['io'] = io
parameters['bed'] = args.inputBed
parameters['output'] = args.output

#bed parsing and storage
bedF = Bed()
bedF.parse(open(args.inputBed, "r"))

#binning
io.register("bins", parameters['output'], "w")
outBins = io.getIO('bins')
pileup = MPileup(io, parameters, bedF)
pileup.parseAndBin(outBins)

sys.exit(io.unregister("all"))

