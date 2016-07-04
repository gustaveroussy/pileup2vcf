# encoding: utf-8

"""
TODO: Finaliser la détection des segments valides
=> Maximum d'évenements entre deux evenements valide pour considérer un segment valide
=> Explosion combinatoire ? Que faire quand on a 400 evenements sur une région ? Peut-être ajuster automatiquement la valeur et rescanner ? Mais c'est impossible de rescanner sans
reparser. Quoi que ... Toujours possible de devenir plus sévère sur le delta minimum
=> Métrique de qualité des InDels ? Basée sur le nombre d'evenements ?
=> Sortir un plot de la profondeur avec une notation des segments retenus et l'emplacement des points de cassure en pointillés
"""

from GenoPy.Genomic import *
from collections import deque
from GenoPy.Bins import Binning
from itertools import combinations

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
        minFreq = int(parameters['minFreq'])
        
        if abs(self.deltaOfDelta) > maxDeltaDelta:
            self.failedLocFilters.append("Local: DeltaOfDeltaOverflow > {}".format(maxDeltaDelta))
        if self.freq < minFreq:
            self.failedLocFilters.append("Local: Freq < {0}".format(minFreq))
        
        if (len(self.failedLocFilters) != 0):
            passFilter = '; '.join(self.failedLocFilters)
        self.filterState = passFilter
        

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
        self.bins = bins.weirdRegions

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
        self.captureStart = None
        self.weirdRegions = []
        self.goodCandidatesInDels = []
        self.weirdEventCounter = 0
        self.weirdoCounter = 0
        self.capturedBins = []
        self.dequeStack = []

    def reset(self):
        self.nonBed = False
        self.nonBedId = None
        self.weirdStart = False
        self.weirdStop = False
        self.captureStart = None
        self.markIns, self.markDel = False, False
        self.captureStack = False
        self.lastNucleotides = deque([])
        self.captured = []
        self.weirdoCounter = 0

    def __iter__(self):
        for x in self.capturedBins:
            yield x

class AllDisengagements(object):
    def __init__(self):
        self.dic = defaultdict(list)
        self.weirdRegions = []
        
    def add(self, wd):
        bedId = wd.getId()
        self.dic[bedId].append(wd)
        
    def computeSegments(self, bedId, parameters):
        print "in computeSegments"
        #print self.dic
        v = self.dic[bedId]
        stack = self.dic['stack_{}'.format(bedId)]
        #v should be a list
        if len(v) <= 1:
            # Skipping this region
            return
        else:
            # 1) Get combinations of elements 2 by 2
            combs = combinations(enumerate(v), 2)
            #print len(combs)
            # We get an iterator which goes like that: ((1, p1), (2, p2)). Enumerate is used in order to get a relative distance between events
            # so that we can distinguish closely detected events (1, 2 by instance) from more broadly detected events (1, 4 for instance, which means
            # there are two other events between those.
            for ((r1, p1), (r2, p2)) in combs:
                # Get associated stack
                #print ((r1, p1), (r2, p2))
                associatedStack = stack[p1.stackStart - 1:p2.stackStart]
                # r = rank, p = position
                d1 = p1.delta
                d2 = p2.delta
                p1p, p2p = p1.pos.start, p2.pos.start
                remaining = abs(d1) - abs(d2)
                tolerance = float(parameters['toleranceForBigIndels']) * abs(d1) if abs(d1) > abs(d2) else float(parameters['toleranceForBigIndels']) * abs(d2)
                if d1 + d2 > tolerance:
                    print "Discarded {}, {}, {}, {}: d1({}) + d2({}) > tol ({})".format(r1, p1p, r2, p2p, d1, d2, tolerance)
                    continue
                
                if (d1 > 0 and d2 < 0):
                    type = "INS"
                    
                elif (d1 < 0 and d2 > 0):
                    type = "DEL"
                else:
                    # The events have the same delta signs: invalid combination
                    print "Discarded {}, {}, {}, {}: d1({}) and d2({}) have the same sign".format(r1, p1p, r2, p2p, d1, d2)
                    continue
                print "Accepting {}, {}, {}, {} pair: d1({}) and d2({}) < tol({}). DeltaOfDelta: {}".format(r1, p1p, r2, p2p, d1, d2, tolerance, remaining)
                weirdRegion = Genomic(associatedStack[0].chr, associatedStack[0].start, int(associatedStack[-1].start) - 1)
                #weirdBin = Indel(Bin(len(associatedStack), associatedStack, getId(weirdRegion)), initDelta, finalDelta, remaining, deltaofthedelta)
                weirdBin = Indel(Bin(len(associatedStack), associatedStack, getId(weirdRegion)), d1, d2, remaining, 0)
                self.weirdRegions.append(weirdBin)
                        
    def getInfos(self, indels, bedId):
        # 1) Extract coordinates of list indels.captureStack in order to get easy indexing on the regions
        cs = indels.captured
        #print cs, bedId
        csStart = cs[0]
        csStop = cs[-1]
        csStartCoor = csStart.start
        csStopCoor = csStop.start
        
        for wd in self.dic[bedId]:
            wdStart = wd.pos.start
            # Get starting index
            # We know that we started capturing at the first event. So if wd is the first event, then this is 0
            # If wd is a later event, then this is the coordinate where we should slice captureStack.
            # In order not to consume too much memory, I will not give every wd its own stack as I thought at first, but I will
            # record the full stack inside self.dic in the following way: "stack_{}".format(bedId)
            # This way, when we'll wan't to compute the segments, we'll have to 
            wd.stackStart = wdStart - csStartCoor
        self.dic['stack_{}'.format(bedId)] = cs
        return self
            

class WeirdDisengagement(object):
    def __init__(self, pos, beforePos):
        self.pos = pos
        self.beforePos = beforePos
        self.delta = self.pos.depth - self.beforePos.depth
        self.id = pos.id
        self.stackStart = None
        
    def getId(self):
        return self.id

def scanForIndels(parameters):
    '''This function scans for InDels. 
    This is inspired by Birama's functionnality. '''
    # Uses deque and uses a scan window of 50 (default, but changeable, parameters[IndelWindowLength]) nucleotides
    capturedBins = []
    indels = IndelDeque()
    dis = AllDisengagements()
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
                io.log("Resetting DeQue. Old_Id = {} ; New_Id = {}".format(old_id, id), loglevel=4)
                if indels.captureStack is True:
                    dis.getInfos(indels, old_id).computeSegments(old_id, parameters)
                indels.reset()
            
            if len(indels.lastNucleotides) == int(parameters['IndelWindowLength']):
                beforePosition = indels.lastNucleotides[-2]
                beforeDepth = beforePosition.depth
#                 
#                 beforeDelta = int(indels.lastNucleotides[-2].depth) - int(indels.lastNucleotides[-3].depth)
                curDepth = o.depth
                deltaDepth = int(curDepth) - int(beforeDepth)
                absDelta = abs(deltaDepth)
                #print deltaDepth
                #beforeDepth = o.depth
                if (absDelta >= int(parameters['minDeltaDepth'])):
                    io.log("Depth disengagement: {} > {}.".format(absDelta, parameters['minDeltaDepth']), loglevel=1)
                    weird = WeirdDisengagement(o, indels.lastNucleotides[-2])
                    dis.add(weird)
                    indels.captureStack = True
                    if indels.captured == []:
                        indels.captured = [beforePosition, o]
                        #indels.captured = [o]
                    
#                     if not indels.weirdStart:
#                         io.log("Begin of event", loglevel=4)
#                         deltaofthedelta = deltaDepth - beforeDelta
#                         indels.captureStack = True
#                         initDelta = deltaDepth
#                         indels.weirdStart = beforePosition
#                         indels.captured = [beforePosition, o]
#                     else:
#                         finalDelta = deltaDepth
#                         remaining = finalDelta + initDelta
#                         tolerance = float(parameters['toleranceForBigIndels']) * abs(initDelta)
#                         if (abs(remaining) > abs(tolerance)) and (indels.weirdoCounter <= parameters['toleranceFalsePositiveEnd']):
#                             #Added abs in tolerance > abs(finalDelta)
#                             indels.weirdStop = beforePosition
#                             weirdRegion = Genomic(indels.weirdStart.chr, indels.weirdStart.start, (int(indels.weirdStop.start) - 1))
#                             remaining = finalDelta + initDelta
#                             indels.captured.pop()
#                             depth_list = []
#                             weirdBin = Indel(Bin(len(indels.captured), indels.captured, getId(weirdRegion)), initDelta, finalDelta, remaining, deltaofthedelta)
#                             indels.weirdRegions.append((weirdRegion, deltaDepth))
#                             indels.capturedBins.append(weirdBin)
#                             io.log("End of event. Event size: {}, false positive ends: {}".format(len(indels.captured), indels.weirdoCounter))
#                             io.log("Resetting InDel detector after calling", loglevel=4)
#                             indels.reset()
#                             indels.weirdoCounter = 0
#                             indels.weirdEventCounter += 1
#                         else:
#                             indels.weirdoCounter += 1
#                             io.log("This is a false positive event end. Value: {}".format(indels.weirdoCounter))
                else:
                    pass
        else:
            if indels.captureStack:
                dis.getInfos(indels, old_id).computeSegments(old_id, parameters)
            indels.reset()
        old_id = id
    #print '\n'.join([str(x) for x in indels.capturedBins])
    return dis
