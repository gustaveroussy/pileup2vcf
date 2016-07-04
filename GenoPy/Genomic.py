#!/usr/bin/env python
# encoding: utf-8

'''
Created on 3 août 2015
Last Update: 27 avril 2016

@author: Yannick Boursin
@contact: yannick.boursin@gustaveroussy.fr
@license: GNU GPLv3
@organization: Gustave Roussy
@version: 1.3
@todo: integrate binIt program in a __main__ closure
'''

# Default library
from collections import defaultdict, deque
from operator import attrgetter

# Installed library
import scipy


###############################################################################
######################## Exceptions ###########################################
###############################################################################

class GenomicException(Exception):
    '''Generic exception thrown by Genomic-like objects'''
    pass
class BadQualityLength(GenomicException):
    '''Exception thrown if len(quality) != len(nucleotides)'''
    pass
class BadQualityEncoding(GenomicException):
    '''Exception thrown if quality isn't encoded in Illumina 1.6 (Phred+33)'''
    pass
class ObjectHasNoId(GenomicException):
    '''Exception thrown when a Genomic-like object has no ID'''
    pass
class WrongTypeOfData(Exception):
    '''Exception thrown when a function doesn't receive the good kind of data'''
    pass
class MpileupFormatError(Exception):
    '''Exception thrown upon detection of a mistake in MPileup input files'''
    pass
class GenotypingError(Exception):
    '''Exception thrown when an error occured in the genotype calling process'''
    pass
class VCFWriterError(Exception):
    '''Exception thrown by the VCF writer class'''
    pass

###############################################################################
############################  Memoization  ####################################
## Found on activestate.com :
## http://code.activestate.com/recipes/578231-probably-the-fastest-memoization-decorator-in-the-/
###############################################################################

def memoize(f):
    '''Memoization function that was found on activestate.
    It takes a function as argument and can act as a decorator.
    Point is to cache a function's results, making the function faster on the next call
    * Input: x
    * Output: f(x)'''
    class memodict(dict):
        __slots__ = ()
        def __missing__(self, key):
            self[key] = ret = f(key)
            return ret
    return memodict().__getitem__

###############################################################################
#############################  tempToObjects  #################################
###############################################################################

def tempToObjects(a):
    '''Get Bin objects back from file
    * Input: file handler
    * Output: Object generator'''
    first_line = a.readline()
    if (first_line.strip() == "##Pileup2VCF - Bins - v1.1"):
        for line in a:
            if not line.startswith("#"):
                this = line.rstrip('\n').split('\t')
                yield Bin(*this)
    elif (first_line.strip() == "##Pileup2VCF - DepthProfile - v1.1"):
        for line in a:
            if not line.startswith("#"):
                this = line.rstrip('\n').split('\t')
                yield Position(*this)
    else:
        print "Unknown file type"


###############################################################################
######################## Sanger Fastq quality  ################################             
###############################################################################
@memoize
def fastqQualityCache(el):
    '''Memoized function that returns a quality integer based on a character
    * Input: [!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI]
    * Output: int'''
    t = ord(el) - 33
    return t

def fastqQuality(sequence, assertion=False):
    '''This function takes a sequence of quality characters as input and outputs the quality numbers in the same order.
    Notes : It only supports quality for Illumina 1.8+.
    * Input: sequence = [!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI]*seqlen
    * Output qualArray = [0123456789]+ * seqlen'''
    # Small check: input sequence should be same length as output sequence
    inputLength = len(sequence)
    
    # We shift quality by 33 'cause it's Illumina. But we should never have quality higher than 41 or negative quality.
    qualArray = [fastqQualityCache(x) for x in sequence]
    #print zip(sequence, qualArray)
    #print qualArray
    # Checking that everything went fine
    if assertion:
        for el in qualArray:
            if el < 0 or el > 41: raise BadQualityEncoding("Quality score {} is not in the standard ranges for Illumina 1.8+ encoding (0-41).".format(el))
        if not inputLength == len(qualArray):
            raise BadQualityLength("Starting quality sequence had {} length. Converted array length is {}".format(inputLength, len(qualArray)))
    return qualArray

###############################################################################
###########################  Genomic Class  ###################################
###############################################################################


def getId(o):
    '''This function returns the unique region ID of any Genomic object. The id is formed using chromosome, start and stop coordinates
    * Input: <Genomic Object>
    * Output: str, something like that: chr1:1-2'''
    if isinstance(o, Genomic):
        if o.stop is not None:
            return '{0}:{1}-{2}'.format(o.chr, o.start, o.stop)
        else:
            return '{0}:{1}'.format(o.chr, o.start)
    else: raise ObjectHasNoId('Only Genomic objects can have an ID in this module')

class Genomic(object):
    '''This is a generic class that breeds many other classes'''
    def __init__(self, chr, start, stop=None, commentary=[]):
        self.chr = chr
        self.start = int(start)
        if stop is None: self.stop = int(start)
        else: self.stop = int(stop)

        #Give the genomic region a Unique ID
        self.id = getId(self)

        #Eventual commentaries ?
        self.commentary = commentary

    def __contains__(self, childObject):
        chr = childObject.chr
        start = childObject.start
        stop = childObject.stop
        if not isinstance(start, int) or not isinstance(stop, int):
            raise WrongTypeOfData("During genomic object comparison, some coordinates were not of integer types")
        return str(self.chr) == str(chr) and start >= self.start and stop <= self.stop

    def __str__(self):
        return '{0}\t{1}\t{2}\t{3}'.format(self.chr, self.start, self.stop, self.id)


###############################################################################
###########################  Position Class  ##################################
###############################################################################

# Cette fonction prends en entrée un objet position et retourne une séquence flanquante !

def getContextFromPosition(chr, pos, howMuchNucleotides=4):
        '''Returns flanking nucleotides for coordinates
        * Input: <str> chr, <int> pos, <int> howMuchNucleotides
        * Output: <str> sequence'''
        pos = int(pos)
        boundaries = (pos - howMuchNucleotides, pos + howMuchNucleotides)
        #if (boundaries[0] not in self.fasta[chr] or boundaries[1] not in self.fasta[chr]):
        #    return None
        return self.fasta[chr][boundaries[0]:boundaries[1]]


class Position(object):
    '''Abstract Factory for Position objects
    * Output: either NewPosition or TempPosition objects depending of the input
    NB: Both objects are fine in the program'''
    def __new__(cls, *args, **kwargs):
        if len(args) < 10:
            return NewPosition(*args, **kwargs)
        else:
            return TempPosition(*args)



class NewPosition(Genomic):
    """Internal class created by Position class allowing to define a new Position object"""
    def __init__(self, chr, pos_start, quality, mapqual, depth, refNuc, pos_before=None, bed=None, assertion=False):
        super(NewPosition, self).__init__(chr, pos_start)
        self.depth = int(depth)
        self.refNuc = refNuc.upper()
        self.last_position = pos_before
        self.inBed = False
        self.nbAltPassFilter = 0
        self.nbAltFiltered = 0
        self.nbTotalAlt = 0
        if (self.depth != 0):
            self.quality = [int(x) for x in fastqQuality(quality, assertion=assertion)]
            self.mapQual = [int(x) for x in fastqQuality(mapqual)]
            self.avgQuality = sum(self.quality) / len(self.quality)
            self.avgMapQuality = sum(self.mapQual) / len(self.mapQual)
        else:
            self.quality = None
            self.mapQual = None
            self.avgMapQuality = None
            self.avgQuality = None
        if bed:
            self.id, self.inBed = bed.getBedAppartenance(self)

    def __str__(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}'.format(self.chr, self.start, self.depth, self.refNuc, self.last_position, self.inBed, self.quality, self.mapQual, self.avgQuality, self.avgMapQuality, self.nbAltPassFilter, self.nbAltFiltered, self.nbTotalAlt, self.id)


class TempPosition(Genomic):
    """Internal class created by Position class allowing to import Position objects"""
    def __init__(self, chr, start, depth, refNuc, last_position, inBed, quality, mapQual, avgQuality, avgMapQuality, nbAltPassFilter, nbAltFiltered, nbTotalAlt, id):
        super(TempPosition, self).__init__(chr, int(start))
        self.depth = int(depth)
        self.refNuc = refNuc.upper()
        self.last_position = int(last_position)
        if (inBed == "True"):
            self.inBed = True
        else:
            self.inBed = False
        self.quality = list(quality)
        self.mapQual = list(mapQual)
        self.avgQuality = avgQuality
        self.avgMapQuality = avgMapQuality
        self.nbAltPassFilter = nbAltPassFilter
        self.nbAltFiltered = nbAltFiltered
        self.nbTotalAlt = nbTotalAlt
        self.id = id

    def __str__(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}'.format(self.chr, self.start, self.depth, self.refNuc, self.last_position, self.inBed, self.quality, self.mapQual, self.avgQuality, self.avgMapQuality, self.nbAltPassFilter, self.nbAltFiltered, self.nbTotalAlt, self.id)

###############################################################################
################################  Bed Class  ##################################
###############################################################################


class Bed(object):
    '''This object is intented to store the genomic coordinates of a Bed file. This will be used, 
    for instance, in Binning and, later in Indel calling'''
    def __init__(self):
        self.capture = defaultdict(list)
        self.cache = defaultdict(list)
        self.sorted_cache = defaultdict(list)
        self.cached = None
        # self.sorted_cache is kinda useless now, but I keep it until my code is stable
        
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
                
        # Now we sort the positions in bed, using "start" values as indexes
        # For now, we do not care about overlapping regions. Each time we
        # return a value, it will be cached to be returned immediatly if request
        # Each interval will have an expiration position (stop position)
        # If we get to that position, we elect a new interval (the next one) and check that
        # the current position is inside it. If it is not, we will check that the
        # stop nucleotide hasn't already expired. If it has, we elect a new interval until we get
        # one for which it hasn't expired. If it hasn't, this becomes the new cached interval.
        
        for k,v in self.cache.iteritems():
            sorted_list = sorted(v, key=attrgetter('start'), reverse=True)
            self.sorted_cache[k] = sorted_list
            self.cache[k] = self.sorted_cache[k]
            del self.sorted_cache[k]
            

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
                    print "Electing new start since {}:{}-{} > {}:{}-{}".format(chr, start, stop, self.cached.chr, self.cached.start, self.cached.stop)
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
    def write(self, fh):
        '''Intended as a dumping function for BED objects.
        * Input: <file>file_handler
        * Output: prints data to file'''
        for v in self.capture.itervalues():
            for el in v:
                print >>fh, '\t'.join([el.chr, str(el.start), str(el.stop), el.commentary, str(el.id)])

###############################################################################
################################  Bin Class  ##################################
###############################################################################


class Bin(object):
    '''Factory for Bin objects.
    * Output: returns either a NewBin or a TempBin object depending on input'''
    def __new__(self, *args):
        if (len(args) == 3):
            a = NewBin(*args)
            return a
        elif (len(args) == 9):
            a = TempBin(*args)
            return a
        else:
            return None


class TempBin(Genomic):
    '''This allows to create Bin objects from temporary files'''
    def __init__(self, chr, start, stop, binSize, binAvgDepth, binMedianDepth, binVarDepth, binMadDepth, binAvgQuality, binAvgMapQuality, gcContent, filtAlt, totalAlt, id, sequence):
        super(TempBin, self).__init__(chr, start, stop)
        self.binSize = int(binSize)
        self.binAvgDepth = float(binAvgDepth)
        self.binAvgQuality = float(binAvgQuality)
        self.binAvgMapQuality = float(binAvgMapQuality)
        self.gcContent = float(gcContent)
        self.id = id
        self.filtAlt = filtAlt
        self.totalAlt = totalAlt
        self.sequence = sequence
        self.binVarDepth = float(binVarDepth)
        self.binMadDepth = float(binMadDepth)
        self.binMedianDepth = float(binMedianDepth)

    def __str__(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{12}\t{13}\t{14}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format(self.chr, self.start, self.stop, self.binSize, self.binAvgDepth, self.binAvgQuality, self.binAvgMapQuality, self.gcContent, self.filtAlt, self.totalAlt, self.id, self.sequence, self.binMedianDepth, self.binVarDepth, self.binMadDepth)


class NewBin(Genomic):
    '''This allows to create Bin objects from real data.
    Given a Position list, each of them will be parsed and the bin object will be created accordingly.
    Continuity should be checked for any bin, meaning that in the position list, Position must be linked but this won't be supported there in case I need otherwise'''
    def __init__(self, binSize, Positions, id):
        self.binSize = int(binSize)
        super(NewBin, self).__init__(Positions[0].chr, Positions[0].start)
        qualArray = []
        mapQualArray = []
        depthArray = []
        gc = 0
        filtAltArray = []
        totAltArray = []
        seqArray = []
        
        for el in Positions:
            # Here, coordinates of the bin are dynamically changed.
            if int(el.start) >= int(self.stop):
                self.stop = el.start
            if int(el.start) <= int(self.start):
                self.start = el.start
            qualArray.append(float(el.avgQuality))
            mapQualArray.append(float(el.avgMapQuality))
            depthArray.append(int(el.depth))
            seqArray.append(el.refNuc)
            filtAltArray.append(int(el.nbAltPassFilter))
            totAltArray.append(int(el.nbTotalAlt))
            if (el.refNuc == "C" or el.refNuc == "G"):
                gc += 1 # Possible bias ? Remove N nucleotides ?
        self.binAvgQuality = sum(qualArray) / len(qualArray)
        self.binAvgMapQuality = sum(mapQualArray) / len(mapQualArray)
        self.binAvgDepth = sum(depthArray) / len(depthArray)
        self.binMedianDepth = scipy.median(depthArray)
        self.binVarDepth = scipy.var(depthArray)
        self.binMadDepth = scipy.median([abs(x - self.binMedianDepth) for x in depthArray])
        self.gcContent = float(gc) / len(Positions)
        self.filtAlt = sum(filtAltArray) if filtAltArray != [] else '.'
        self.totalAlt = sum(totAltArray) if totAltArray != [] else '.'
        #self.sequence = "".join(seqArray)
        self.sequence = ''.join(seqArray)
        self.id = id

    def __str__(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{12}\t{13}\t{14}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format(self.chr, self.start, self.stop, self.binSize, self.binAvgDepth, self.binAvgQuality, self.binAvgMapQuality, self.gcContent, self.filtAlt, self.totalAlt, self.id, self.sequence, self.binMedianDepth, self.binVarDepth, self.binMadDepth)

class MinimalBin(object):
    '''WriteMe'''
    def __init__(self, chr, start, stop, gc):
        self.chr = chr
        self.start = start
        self.stop = stop
        self.binSize = stop - start
        self.gc = float(gc) * 100 / float(self.binSize)
    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}'.format(self.chr, self.start, self.stop, self.binSize, self.gc)

class BinIt(object):
    '''This class allows to bin a bedfile with no further data. It cuts the genomic intervals in Bin-like
    objects with respect for the minBinSize and binSizeMax parameters.
    * Input: <str> bedPath, <str> refGenomePath'''
    def __init__(self, bed, refGenome):
        self.bed = bed
        self.fasta = pyfaidx.Fasta(refGenome)
        
    def importBed(self):
        '''This method will get BED file in memory, in an orderedDict structure
        * Input: 
        * Output: <self>
        NB: This method needs no input and can be seen as a complementary to __init__. It is a
        prerequisite for further usage of the object but will not be done upon object creation and
        initialization.
        '''
        bh = open(self.bed, "r")
        reg = bh.readlines()
        bh.close()
        reg = [x.rstrip('\n').split('\t') for x in reg]
        self.regions = OrderedDict()
        oldChr = None
        for el in bed:
            if oldChr != el[0]:
                oldChr = el[0]
                current = self.regions[el[0]] = []
            current.append((int(el[1]), int(el[2])))
            oldChr = el[0]
        return self
    
    def makeBins(self, minBinSize, binSizeMax):
        '''This method will create Bin objects (MinimalBin) using the bed file and the fasta file.
        * Input: <int> minBinSize, <int> binSizeMax
        * Output: <self>
        NB: this method returns self. Data is stored in <OrderedDict> self.binDic'''
        self.binDic = OrderedDict()
        for chr, tup in self.regions.iteritems():
            self.binDic[chr] = []
            for start, stop in tup:
                # First, we get the bin Sequence
                seq = self.fasta[chr][start:stop]
                length = len(seq)
                seq = deque(str(seq))
                gcCounter = 0
                curPos = start
                binSizeCounter = 0
                while curPos != stop:
                    curNuc = seq.popleft()
                    # This is the first nuc
                    if curNuc in ['G', 'C']:
                        gcCounter += 1
                    binSizeCounter+=1
                    curPos += 1
                    if binSizeCounter == binSizeMax:
                        binStart = curPos - binSizeCounter
                        self.binDic[chr].append(MinimalBin(chr, curPos - binSizeCounter, curPos, gcCounter))
                        binSizeCounter = 0
                        gcCounter = 0
                    if curPos == stop:
                        if binSizeCounter > minBinSize:
                            self.binDic[chr].append(MinimalBin(chr, curPos - binSizeCounter, curPos, gcCounter))
                print 'Processed region {}-{}'.format(start, stop)
            print 'Processed {}'.format(chr)
        return self
    def printIt(self, output):
        '''This method will output bins created in this class'''
        for k,v in self.binDic:
            for bin in v:
                print >>output, str(bin)
