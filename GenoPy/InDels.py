# encoding: utf-8

from GenoPy.Genomic import *
from collections import deque
from GenoPy.Bins import Binning

class Indel(TempBin):
    '''Indel class
    This class is responsible for all Indels that are called by IndelDetector. '''
    def __init__(self, binObject, deltaDepthInit, deltaDepthFinal, deltaOfDelta, justAnotherDelta):
        listBinObject = str(binObject).split('\t')
        super(Indel, self).__init__(*listBinObject)
  
	#Inherited attributes
	# 
	#self.chr
	#self.start
	#self.stop
	#self.binSize = int(binSize)
        #self.binAvgDepth = float(binAvgDepth)
        #self.binAvgQuality = float(binAvgQuality)
        #self.binAvgMapQuality = float(binAvgMapQuality)
        #self.gcContent = float(gcContent)
        #self.id = id
        #self.filtAlt = filtAlt
        #self.totalAlt = totalAlt
        #self.sequence = sequence
        #self.binVarDepth = binVarDepth
        #self.binMadDepth = binMadDepth
        #self.binMedianDepth = binMedianDepth

        self.deltaDepthInit = float(deltaDepthInit)
        self.deltaDepthFinal = float(deltaDepthFinal)
        self.deltaOfDelta = float(deltaOfDelta)
        self.justAnotherDelta = float(justAnotherDelta)
        self.filterState = None
        self.initDepth = self.binAvgDepth - self.deltaDepthInit
        self.ad = '{},{}'.format(self.initDepth, self.binAvgDepth)
        self.dp = self.binAvgDepth

        if self.deltaDepthInit < 0:
            self.ref = self.sequence
            self.alt = self.sequence[0]
            self.type = "DEL"
            self.freq = 100 - (100 * float(self.binAvgDepth) / float(self.initDepth))
        else:
            self.ref = self.sequence[0]
            self.alt = self.sequence
            self.type = "INS"
            self.freq = 100 - (100 * float(self.initDepth) / float(self.binAvgDepth))

    def genotype(self, parameters):
        '''Sets a proper genotype for InDels found by IndelDetector'''
        self.gt = '0/0'
        if self.freq > float(parameters['hetTreshold']):
            self.gt = '0/1'
            if self.freq > float(parameters['homTreshold']):
                self.gt = '1/1'
        self.filterIt(parameters)
        return self

    def filterIt(self, parameters):
        '''Applies various filters on Indel object'''
        #Apply various filters on indels
        # TODO
        passFilter = "PASS"
        self.failedLocFilters = []
        
        maxDeltaDelta = int(parameters['maxDeltaDelta'])
        minFreq = self.parameters['minFreq']
        
        if abs(self.deltaOfDelta) > maxDeltaDelta:
            self.failedLocFilters.append("Local: DeltaOfDeltaOverflow > {}".format(maxDeltaDelta))
        if self.freq < minFreq:
            self.failedLocFilters.append("Local: Freq < {0}".format(minFreq))
        
        if (len(self.failedLocFilters) != 0):
            passFilter = '; '.join(self.failedLocFilters)
        return passFilter
        

    def __str__(self):
        try:
            return "{0}\t{1}\t.\t{2}\t{3}\tNone\t{9}\t{4}\tGT:FREQ:AD:DP\t{5}:{6:.2f}:{7}:{8}".format(self.chr,
                                                                                                    self.start,
                                                                                                    self.ref,
                                                                                                    self.alt,
                                                                                                    "TYPE={0};DELTADELTA={1:.0f};DEPTHMAD={2:.0f};DEPTHMEDIAN={3:.0f};DEPTHVARIANCE={4:.0f};BIGINDELDETECTOR".format(self.type,
                                                                                                                                                                                            self.deltaOfDelta,
                                                                                                                                                                                            self.binMadDepth,
                                                                                                                                                                                            self.binMedianDepth,
                                                                                                                                                                                            self.binVarDepth),
                                                                                                    self.gt,
                                                                                                    self.freq,
                                                                                                    self.ad,
                                                                                                    self.dp,
                                                                                                    self.filterState)
        except:
            try:
                self.gt
            except:
                raise GenotypingError('Tried to output a variant before genotyping.')
            raise #VCFWriterError('Cannot write InDel for unknown reason')


class InDelDetector:
    '''Holder class containing all Indels detected by scanForIndels function'''
    def __init__(self, bins):
        self.allIndels = []
        self.bins = bins

    def extract(self):
        '''WriteMe'''
        self.allIndels = self.bins
        return self

    def getIndels(self):
        '''WriteMe'''
        return self.allIndels
    
    def genotype(self, parameters):
        for el in self.allIndels:
            el.genotype(parameters)
        return self

class IndelDeque():
    '''Holds variables used by scanForIndels.
    Scan for indels uses a deque and must be resettable. This object implements a reset method 
    that makes it able to get some initial state without resetting everything (like self.goodCandidatesInDels), which
    is necessary maintained to its state throughout scanning.'''
    def __init__(self):
        self.lastNucleotides = deque([])
        self.nonBed = False
        self.nonBedId = None
        self.weirdStart = False
        self.weirdStop = False
        self.markIns, self.markDel = False, False
        self.captureStack = False
        self.captured = []
        self.weirdRegions = []
        self.goodCandidatesInDels = []
        self.weirdEventCounter = 0
        self.weirdoCounter = 0
        self.capturedBins = []

    def reset(self):
        self.nonBed = False
        self.nonBedId = None
        self.weirdStart = False
        self.weirdStop = False
        self.markIns, self.markDel = False, False
        self.captureStack = False
        self.lastNucleotides = deque([])
        self.captured = []
        self.weirdoCounter = 0

    def __iter__(self):
        for x in self.capturedBins:
            yield x


def scanForIndels(parameters):
    '''This function scans for InDels. 
    This is inspired by Birama's functionnality. '''
    # Uses deque and uses a scan window of 50 (default, but changeable, parameters[IndelWindowLength]) nucleotides
    capturedBins = []
    indels = IndelDeque()
    io = parameters['io']
    depthProfileH = io.getIO("depthProfileR2")
    old_id = None
    # We gather back the depthProfile 
    for o in tempToObjects(depthProfileH):
        # This maintains the dequeu to a minimum length of 50
        if (len(indels.lastNucleotides) != parameters['IndelWindowLength']):
            indels.lastNucleotides.append(o)
        else:
            indels.lastNucleotides.popleft()
            indels.lastNucleotides.append(o)

        # Making sure we're acting on a specific and continuous regions
        if indels.nonBed is True:
            id = indels.nonBedId
            if Binning.positionDrift(o):
                id = getId(o)
        else:
            id = o.id
        if (old_id is None):
            old_id = id
        if indels.captureStack:
            indels.captured.append(o)

        if Binning.checkBedAppartenance(o):
            if old_id != id:
                indels.nonBed = False
                io.log("Resetting DeQue. Old_Id = {} ; New_Id = {}".format(old_id, id), loglevel=4)
                indels.lastNucleotides = deque([])
                if (indels.weirdStart and not indels.weirdStop):
                    io.log("Resetting InDel detector, since we cannot call InDels between different regions", loglevel=4)
                    indels.reset()
            
            if len(indels.lastNucleotides) == int(parameters['IndelWindowLength']):
                beforeDepth = indels.lastNucleotides[-2].depth # lastFifty[-1] is the current position
                beforePosition = indels.lastNucleotides[-2]
                beforeDelta = int(indels.lastNucleotides[-2].depth) - int(indels.lastNucleotides[-3].depth)
                curDepth = o.depth
                deltaDepth = int(curDepth) - int(beforeDepth)
                absDelta = abs(deltaDepth)
                #print deltaDepth
                beforeDepth = o.depth
                if (absDelta >= int(parameters['minDeltaDepth'])):
                    #print str(o)
                    io.log("Depth disengagement: {} > {}.".format(absDelta, parameters['minDeltaDepth']), loglevel=4)
                    if not indels.weirdStart:
                        io.log("Begin of event", loglevel=4)
                        deltaofthedelta = deltaDepth - beforeDelta
                        indels.captureStack = True
                        initDelta = deltaDepth
                        if (initDelta < 0):
                            markDel = True
                            markIns = False
                        else :
                            markIns = True
                            markDel = False
                        indels.weirdStart = beforePosition
                        indels.captured = [beforePosition, o]
                    else:
                        finalDelta = deltaDepth
                        remaining = finalDelta + initDelta
                        tolerance = float(parameters['toleranceForBigIndels']) * abs(initDelta)
                        if (abs(initDelta) - tolerance < abs(finalDelta)) and (abs(initDelta) + tolerance > abs(finalDelta)) and (indels.weirdoCounter <= parameters['toleranceFalsePositiveEnd']):
                            #Added abs in tolerance > abs(finalDelta)
                            indels.weirdStop = beforePosition
                            weirdRegion = Genomic(indels.weirdStart.chr, indels.weirdStart.start, (int(indels.weirdStop.start) - 1))
                            remaining = finalDelta + initDelta
                            indels.captured.pop()
                            depth_list = []
                            #if (abs(remaining) >= parameters['maxDeltaDelta']):
                            #    weirdBin = str(Bin(len(indels.captured), indels.captured, getId(weirdRegion))) + "\t" + str(initDelta) + "\t" + str(finalDelta) + "\t" + str(remaining) + "\t" + "{1}\tDeltaDelta>{0}".format(str(parameters['maxDeltaDelta']), str(deltaofthedelta))
                            #else:
                            weirdBin = Indel(Bin(len(indels.captured), indels.captured, getId(weirdRegion)), initDelta, finalDelta, remaining, deltaofthedelta)
                            indels.weirdRegions.append((weirdRegion, deltaDepth))
                            indels.capturedBins.append(weirdBin)
                            io.log("End of event. Event size: {}, false positive ends: {}".format(len(indels.captured), indels.weirdoCounter))
                            io.log("Resetting InDel detector after calling", loglevel=4)
                            indels.reset()
                            indels.weirdoCounter = 0
                            indels.weirdEventCounter += 1
                        else:
                            indels.weirdoCounter += 1
                            io.log("This is a false positive event end. Value: {}".format(indels.weirdoCounter))
                else:
                    pass
        else:
            indels.reset()
        old_id = id
    #print '\n'.join([str(x) for x in indels.capturedBins])
    return indels
