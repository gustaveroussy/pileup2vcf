#encoding: utf-8
#Bins
from GenoPy.Genomic import tempToObjects, Bin, getId

class BinTemporaryObject(object):
    '''This object holds information required during Bin calling.'''
    def __init__(self, binSizeMin, binSizeMax, binh, io):
        self.binSizeCounter = 1
        self.stack = []
        self.nonBed = False
        self.nonBedId = False
        self.old_id = None
        self.segments = []
        self.beforeDepth = 0
        self.beforePosition = None
        self.mean = []
        self.capturedBins = []
        self.binSizeMin = binSizeMin
        self.binSizeMax = binSizeMax
        self.id = None
        self.binh = binh
        self.o = None
        self.log = io.log

    def addBin(self):
        '''Prints a bin into the bin file'''
        if (self.binSizeCounter > self.binSizeMin):
            if Binning.checkStackBedAppartenance(self.stack):
                print >>self.binh, str(Bin(self.binSizeCounter, self.stack, self.id))
            else: self.log("Bin is not in BED", loglevel=3)
                #self.log(str(Bin(self.binSizeCounter, self.stack, self.id)), loglevel=5)
        self.reset()
            
    def reset(self):
        '''Resets the BinTemporaryObject stack'''
        self.stack = [self.o]
        self.old_id = self.o.id
        self.binSizeCounter = 1



class Binning(object):
    '''Binning class
    This class takes charge of binning.'''
    def __init__(self, binh, parameters):
        self.binh = binh
        self.parameters = parameters
        self.io = parameters['io']
        self.header = True
        
    def bin(self, depthProfileH, binDynamically=False):
        '''Bins the depthProfile'''
        self.depthProfileH = depthProfileH
        # What we want:
        # Découper self.depthProfile(chr) par binSize, et découper par décrochage de depth (adaptation aux captures)
        # En pratique: si une position a 20 bp de différence avec sa plus proche voisine, caller la fin du bin
        # En pratique: on découpe une première fois le tableau en scannant ces décrochages
        # En pratique: puis on découpe les sous tableaux avec binSize.
        # Proposer un binning adaptatif (avec un facteur de découpe, genre 5, qui signifierait découper chaque sous tableau en 5)
        # Proposer la configuration de la différence de position pour caller la fin du bin.
        
        binSizeMin, binSizeMax = self.parameters['minBinSize'], self.parameters['maxBinSize']
        if self.header:
            print >>self.binh, "##Pileup2VCF - Bins - v1.1"
            print >>self.binh, "#Chromosome\tbinStart\tbinStop\tbinSize\tAvgDepth\tmedianDepth\tvarDepth\tMADDepth\tAvgQuality\tAvgMapQuality\tAvgGCPercent\tFilteredAltEvents\tTotalAltEvents\tBedRecord"
            self.header = False
            
        t = BinTemporaryObject(binSizeMin, binSizeMax, self.binh, self.io)
        
        def binFeed():
            '''Feeds data into bin method'''
            for k, v in self.depthProfileH.items():
                for v2 in v.values():
                    yield v2
        # Monkey-patching
        if not binDynamically:
            binFeed = lambda: tempToObjects(self.depthProfileH)
        ## WE SHOULD CHECK ABOUT SAME CHROMOSOMENESS AGAIN !!!!!!!!!!!!!!!!
        for o in binFeed(): #We can use this since position is an Ordered Dict.
            t.o = o
            t.mean.append(o.depth)
            if t.nonBed is True:
                t.id = t.nonBedId
                if self.positionDrift(o):
                    t.id = getId(o)
            else:
                t.id = o.id
            if t.old_id is None:
                t.old_id = t.id
            # Arrive quand on a deux zones d'un BED qui overlappent.
            if self.checkBedAppartenance(o):
                if t.old_id != t.id:
                    self.io.log('Switching region in bed.', loglevel=4)
                    t.nonBed = False
                    t.addBin()
                    continue
                if t.nonBed:
                    self.io.log("We're back in bed.", loglevel=4)
                    t.nonBed = False
                    t.addBin()
                    continue
            else:
                if not t.nonBed:
                    #Switching to non-bed regions
                    self.io.log('Not in bed anymore', loglevel=4)
                    t.nonBedId = getId(o)
                    t.nonBed = True
                    t.addBin()
                    continue

            if (self.positionDrift(o)): #If the bin driffted x nucleotides ... call the current bin if we can and start a new one, or discard current bin and start a new one.
                self.io.log("We're in position drift", loglevel=4)
                t.addBin()
                continue

            if (t.binSizeCounter == binSizeMax): #If the bin is at the maximum size, call it and start a new one.
                self.io.log("We're at maximum binSize", loglevel=4)
                t.addBin()
                continue

            t.stack.append(o)
            t.binSizeCounter += 1
            t.old_id = t.id

        self.binh.flush()

    @staticmethod
    def positionDrift(position, maxDrift=50):
        '''Static method that detects position drifting in pileup files'''
        old_position = int(position.last_position)
        current_position = int(position.start)
        if (current_position - old_position > maxDrift):
            return True
        return False

    @staticmethod
    def checkBedAppartenance(position):
        '''returns a boolean that states if a position is inside bed file'''
        return position.inBed

    @staticmethod
    def checkStackBedAppartenance(stack):
        '''Checks if the current stack is inside bed as a whole (all positions are in bed)'''
        toReturn = True
        #print len(stack)
        if len(stack) == 0:
            return False
        for el in stack:
            toReturn = toReturn and el.inBed
            #print el
        #print toReturn
        return toReturn
        
    def printBinning(self, binningR):
        '''Not used anymore. Used to perform segmentation and normalization 'CGH-like' on previously
        called bins'''
        depth_array, gc_array = [], []
        for el in tempToObjects(binningR):
            #self.log(el, loglevel=2)
            depth_array.append(float(el.binAvgDepth))
            gc_array.append(float(el.gcContent))
        logarized = [np.log10(x) for x in depth_array]
        #LOESS = [x[0] for x in lowess(gc_array, logarized, r eturn_sorted=True)]
        LOESS = lowess(np.array(gc_array, np.float), np.array(logarized,np.float))
        #print LOESS
        #print np.arange(len(self.bins))[1]
        toPlot = LOESS
        #varNormalizeBy = np.var(depth_array)
        #meanNormalizeBy = np.mean(depth_array)
        plt.subplot(3, 1, 1)
        fig = plt.plot(np.arange(len(depth_array)), toPlot)
        plt.subplot(3, 1, 2)
        fig2 = plt.plot(np.arange(len(depth_array)), logarized) #> kinda nice with binSize 100 
        plt.subplot(3, 1, 3)
        fig3 = plt.plot(np.arange(len(depth_array)), depth_array)
        #print meanNormalizeBy
        #print varNormalizeBy
        plt.show()

    def alternateBinningPrinting(self, binningR):
        '''Not used anymore. Used to display graphical representation of copy number profile'''
        a = defaultdict(list)
        for el in tempToObjects(binningR):
            a[el.id].append(el)
        #dico_len = len(a)
        i = 0
        for k, v in a.iteritems():
            i += 1
            tada = []
            depth_array = []
            lines = []
            id = None
            for el in v:
                depth_array.append(float(el.binAvgDepth))
                tada.append(int(el.start))
                id = el.id
            for s in self.segments:
                if (s.id == id):
                    lines.append([(int(s.start), s.binAvgDepth), (int(s.stop), int(s.binAvgDepth))])
            print lines
            #plt.subplot(3,1,1)
            #plt.plot(tada, depth_array2)
            #plt.subplot(3,1,2)
            #plt.plot(tada, depth_array)
            #plt.show()
            #print tada
            #print depth_array
            lc = mc.LineCollection(lines, linewidths=2)
            fig = plt.figure()
            ax = fig.add_subplot(111, axisbg=(1.0,1.0,1.0,0))
            ax.set_title('Report for {0}'.format(k), y=1.05)
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            ax2.add_collection(lc)
            ax1.margins(0.1)
            ax2.margins(0.1)
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
            ax.tick_params(labelcolor=(1.0,1.0,1.0,0), top='off', bottom='off', left='off', right='off')
            ax1.plot(tada, depth_array)
            start, end = ax1.get_xlim()
            #ax1.xaxis.set_ticks(np.arange(start, end))
            #ax2.xaxis.set_ticks(np.arange(start, end))
            ax1.xaxis.set_ticklabels(tada, rotation=30)
            ax2.xaxis.set_ticklabels(tada, rotation=30)
            ax.set_xlabel('Position')
            ax.set_ylabel('Depth variations')
            ax1.set_title('Raw depth')
            ax2.set_title('Segmented profile')
            #fig.show()
            #plt.show()
            plt.tight_layout()
            fig.savefig("Fig{0}.png".format(i))
            plt.close()            
