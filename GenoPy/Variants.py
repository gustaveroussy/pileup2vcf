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
@todo: Allele.computeQuality method, biggerContext and context cohabitation is useless
'''
# Default library
from operator import itemgetter, attrgetter
from itertools import islice
from re import compile as recompile

# Installed Libraries
from fisher import pvalue as fisher
import numpy as np

# My modules
from GenoPy.Genomic import *
from GenoPy.Statistics import Statistics


class UnclassifiedVariantCollection(object):
    '''UnclassifiedVariantCollection object: contains various types of variants
    
    Variants contained may be splitted in two categories: non-filtered (referred as "good variants"),
    filtered (referred as "trashed variants"). Filters may be different whether the variant comes from
    classical variant detection or depth-based InDel detection.
    This class allows for applying filters on the entire collection
    
    * Input: <list> arrayPileup, <defaultdict> parameters'''
    def __init__(self, arrayPileup, parameters):
        self.col = arrayPileup
        self.call = []
        self.trash, self.good = TrashCollection(parameters), GoodCollection(parameters)
        self.GlobalStatistics = Statistics()
        self.counter = 0
        self.parameters = parameters
        self.io = parameters['io']

    def classify_variants(self):
        '''Here, we apply filtering, genotyping and all on all variants. Two VariantCollection are outputed. VCF object takes a VariantCollection and
	a file handler as input to output VariantCollections to files'''
        passFilter = [cur.getGlobalFilterState() for cur in self.col] # Apply global filters on variants
        self.io.log("There are {} variant positions with no filtering or calling".format(len(self.col)), loglevel=2)
        [self._getCall(x) if y is True else None for x,y in zip(self.col, passFilter)] # Get call for variants that passed global filters
        return self
    
    def _getCall(self, cur):
        '''Classify alleles between variants and outputs a genotype
        
        This method calls several Variant-like object methods: call() and genotypeVariant(). This is where classifying really
        occurs as call() methods returns three lists containing Allele objects:
        1) RefAllele: always present (gives information regarding the reference allele
        2) goodAlleles: list containing only informations about the alleles that DID pass filters
        3) trashedAlleles: list containing only informations about the alleles that DID NOT pass the filters
        
        goodAlleles and trashedAlleles: either <None> or [Allele, ... Allele]
        genotypeVariant method shouldn't be there I think, i will link it inside Variant'''
        self.counter += 1
        refAllele, goodAlleles, trashedAlleles = cur.call()
        if goodAlleles is None and trashedAlleles is None:
            return None
        elif goodAlleles is None and trashedAlleles is not None:
            self.trash.add(cur)
        elif goodAlleles is not None and trashedAlleles is not None:
            self.good.add(cur)
        elif goodAlleles is not None and trashedAlleles is None:
            self.good.add(cur)
        else:
            self.io.log('I do not understand', loglevel=1)
        return True
    
    def getCollections(self):
        '''Getter function for good and trash collection'''
        return self.good, self.trash

class AllClassifiedVariantCollections(object):
    '''This class gathers all other ClassifiedVariantCollection classes'''
    def __init__(self, parameters):
        self.all = {}
        self.statistics = Statistics()
        self.parameters = parameters
        self.statistics.init_positions(parameters)
        
    def add(self, type, vc):
        self.all[type] = vc
        
    def do_stats(self, stath):
        for k,v in self.all.iteritems():
            self.statistics.init_variants(v, type=k)
        for k,v in self.all.iteritems():
            self.statistics.printStatistics(k, stath)
        

class ClassifiedVariantCollection(object):
    '''Class intended to contain the final step of processed variant
    
    This class is only a container for final results of SNP detection.
    All variants are stored in <list> self.vl
    
    * Input: <defaultdict> parameters'''
    def __init__(self, parameters):
        self.variantList = self.vl = []
        self.counter = 0
        self.parameters = parameters
        self.io = parameters['io']
        self.iterCounter = 0
        
    def add(self, cur):
        '''Method intended to add a single variant to a classified collection
        
        * Input: <Variant> cur'''
        #self.statistics.add(cur)
        self.counter += 1
        self.vl.append(cur)
        
    def getAll(self):
        '''Returns a list of all classified variants'''
        return self.vl


class TrashCollection(ClassifiedVariantCollection):
    '''Specialized class for ClassifiedVariantCollection: Trash'''
    def __init__(self, parameters):
        super(TrashCollection, self).__init__(parameters)
        self.isTrash = True

class GoodCollection(ClassifiedVariantCollection):
    '''Specialized class for ClassifiedVariantCollection: Good'''
    def __init__(self, parameters):
        super(GoodCollection, self).__init__(parameters)
        self.isTrash = False

class Fisher(object):
    '''Class intended as a fast Fisher-test worker with memoization'''
    def __init__(self):
        self.cache = defaultdict(bool)
        
    def computeFisher(self, ref1, ref2, alt1, alt2):
        '''Fisher computing function with memoization (t-test)
        
        * Input: <int> ref1, <int> ref2, <int> alt1, <int> alt2
        * Output: <float> phred_pval
        '''
        hash = '{}-{}-{}-{}'.format(ref1, ref2, alt1, alt2)
        if not self.cache[hash]:
            pvalue = fisher(ref1, ref2, alt1, alt2)
            phred_pval = -10 * np.log10(pvalue.two_tail)
            self.cache[hash] = phred_pval
            return phred_pval
        else:
            return self.cache[hash]

fi = Fisher()

###############################################################################
############################  INDEL Class  ####################################
###############################################################################
class Indel(object):
    '''Regular base holder class  for Indel variants'''
    def __init__(self):
        self.sequence = ''
        self.depth = 0
        self.size = 0
        self.freq = 0
        self.plus = 0
        self.minus = 0
        self.type = ''

    def __str__(self):
        return 'Seq: {0}, DP: {1}, Size: {2}, Freq: {3}, +: {4}, -:{5}, type: {6}'.format(self.sequence, self.depth, self.size, self.freq, self.plus, self.minus, self.type)

###############################################################################
############################  Allele Class  ####################################
###############################################################################

class Allele(object):
    '''Allele class
    
    Each instance of this class represents one Allele. We expect variants to have at least
    two alleles. Local filters are applied on Allele objects.
    
    * Input: <defaultdict> parameters'''
    def __init__(self, parameters):
        self.Allele = ''
        self.rank = "NA"
        self.Depth = ''
        self.insertion = False
        self.deletion = False
        self.Quality = 0
        self.Frequency = 0
        self.stackedQuality = 0

        self.parameters = parameters

        #Utile pour calculer le SB
        self.forward = 0
        self.reverse = 0
        self.referenceForward = 0
        self.referenceReverse = 0
        self.baseQualArray = []
        self.readMapQualArray = []

    def __str__(self):
        a = []
        a.append("Allelic description")
        a.append("Allele: {0}; Depth: {1}; Frequency: {2}".format(self.Allele, self.Depth, self.Frequency))
        try:
            a.append("Local Filters: {0}".format(self.localFilter))
        except:
            pass
        try:
            a.append("Rank: {0}".format(self.rank))
        except:
            pass
        return '\n'.join(a)

    def computeQuality(self):
        '''Dummy quality function'''
        #print self.stackedQuality
        #print np.log10(float(self.stackedQuality) + float(0.00001))
        #print "Computed Quality for {0}".format(self.Allele)
        #print 2 * (float(self.Frequency) / np.log10(float(self.stackedQuality) + float(0.00001)))
        #return 2 * (float(self.Frequency) / np.log10(float(self.stackedQuality) + float(0.00001)))
        return 50

    def getLocalFilterState(self):
        '''Method intended to apply filters on Allele.
        
        For now, there are several filters applied:
            - Relative strand bias (LocalStrandBias)
            - Allelic frequency (Freq)
            - Number of reads supporting (AltReads)
            - Fisher computed strand bias (LocalFisherStrand)
        
        * Output: <list> filters'''
        passFilter = True
        self.failedLocFilters = []
        
        localsb = float(self.forward) * 100 / (float(self.reverse) + float(self.forward))
        self.localsb = localsb
        reverseLocalStrandBias = abs(100 - self.parameters['localStrandBias'])
        paramLocalStrandBias = [self.parameters['localStrandBias'], reverseLocalStrandBias]
        lowerParamLocalStrandBias = min(paramLocalStrandBias)
        greaterParamLocalStrandBias = max(paramLocalStrandBias)
        
        if (float(localsb) > float(greaterParamLocalStrandBias)):
            self.failedLocFilters.append("Local: LocalStrandBias > {0}".format(greaterParamLocalStrandBias))
        if (float(localsb) < float(lowerParamLocalStrandBias)):
            self.failedLocFilters.append("Local: LocalStrandBias < {0}".format(lowerParamLocalStrandBias))
        
        if (float(self.Frequency) < float(self.parameters['minFreq'])):
            self.failedLocFilters.append("Local: Freq < {0}".format(self.parameters['minFreq']))
        if (int(self.Depth) < int(self.parameters['minReadsAlt'])):
            self.failedLocFilters.append("Local: AltReads < {0}".format(self.parameters['minReadsAlt']))
        
        self.lfs = self.computeFisherStrandBias()
        
        if (float(self.lfs) > float(self.parameters['fisherStrand'])):
            self.failedLocFilters.append("Local: LocalFisherStrand > {0}".format(self.parameters['fisherStrand']))
        
        if (len(self.failedLocFilters) != 0):
            passFilter = '; '.join(self.failedLocFilters)
        self.passFilter = passFilter
        return passFilter

    # Cette fonction effectue le test de Fisher exact entre la référence et les alternatifs. Elle donne une mesure de strand bias.
    def computeFisherStrandBias(self):
        '''Method to compute Fisher strand bias on Allele object'''
        
        # On créé la table de contingence
        forwardRef = float(self.referenceForward)
        reverseRef = float(self.referenceReverse)
        forwardAlt = float(self.forward)
        reverseAlt = float(self.reverse)
        
        #On lance le test proprement dit
        return fi.computeFisher(forwardRef, reverseRef, forwardAlt, reverseAlt)

###############################################################################
###########################  Variant Class  ###################################
###############################################################################


class Variant(NewPosition):
    '''Variant class
    
    This class holds Variant objects.
    * Input: <str> chr, <int> pos, <str> refNuc, <int> depth, <str> sequence,
            <list of int> qual, <list of int> mapQual, <defaultdict> parameters,
            <str> context, <str> biggerContext
    context and biggerContext are flanking sequences for the variant. context is to be removed.
    Associated methods:
        - getInfos
        - computeNucCount
        - computeFreq
        - call
        - _genotypeVariant
        - getLocalSequenceContext
        - window
        - suffixWindow
        - getLocalFilterState
        - computeIndelFreq
        - computeStrandBias
        - titv
        - statsAllAllele
        - computeFisherStrandBias
        - isVariant
        - getGlobalFilterState
        - getBaseInformations
        - getFrequencies
        - getInDels
    '''
    def __init__(self, chr, pos, refNuc, depth, sequence, qual, mapQual, parameters, context, biggerContext):
        #Initialisation des variables
        super(Variant, self).__init__(chr, pos, qual, mapQual, depth, refNuc)

        # Reads supportés par le brin sens
        self.forward = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "ins": 0, "del": 0}
        # Reads supportés par le brin anti-sens
        self.reverse = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "ins": 0, "del": 0}
        # Fréquence de chaque évenement
        self.frequencies = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "ins": 0, "del": 0}
        # Comptage de chaque évenement
        self.nucCount = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "ins": 0, "del": 0}

        # Informations sur la taille et la séquence des indels. NB: seule la dernière séquence pour une taille donnée est donnée à titre d'information.
        self.insertionsPlus = defaultdict(int)
        self.deletionsPlus = defaultdict(int)
        self.insertionsMinus = defaultdict(int)
        self.deletionsMinus = defaultdict(int)
        self.insertionsInfo = defaultdict(list)
        self.deletionsInfo = defaultdict(list)
        self.parameters = parameters
        self.allInsertions = []
        self.allDeletions = []

        # Contient des champs du mpileup
        self.sequence = sequence
        self.refCount = 0
        self.filtered = False
        self.context = context
        self.biggerContext = biggerContext

    def __str__(self):
        a = ['\nContext']
        a.append('Chromosome: {0}, position: {1}, reference: {2}'.format(self.chr, self.start, self.refNuc))
        a.append('Reads')
        a.append('Forward: ' + str(self.forward))
        a.append('Reverse: ' + str(self.reverse))
        a.append('Frequencies: ' + str(self.frequencies))
        a.append('NucCount: ' + str(self.nucCount))
        a.append('Indels')
        a.append('Insertions+: ' + str(self.insertionsPlus))
        a.append('Insertions-: ' + str(self.insertionsMinus))
        a.append('Insertion Info: ' + str(self.insertionsInfo))
        a.append('Deletion+: ' + str(self.deletionsPlus))
        a.append('Deletion-: ' + str(self.deletionsMinus))
        a.append('Deletion Info: ' + str(self.deletionsInfo))
        a.append('Pileup context')
        a.append('Sequence: ' + str(self.sequence))
        a.append('Quality: ' + str(self.quality))
        a.append('Depth: ' + str(self.depth))
        a.append('Other informations')
        a.append('Refcount: {0}, refAllele: {2}, Filtered: {1}'.format(self.refCount, self.filtered, self.refNuc))
        try:
            a.append('Called Allele: {0}'.format(self.calledAllele))
        except:
            a.append('Called Allele: None')
        try:
            a.append('Allele Depth | Frequency: {0} | {1}'.format(self.calledDepth, self.calledFrequency))
        except:
            a.append('Allele Depth | Frequency: None | None')
        try:
            a.append('Called Secondary Allele ?: {0}'.format(self.calledSecondaryAllele))
        except:
            a.append('Called Secondary Allele ?: None')
        try:
            a.append('Filters: {0}'.format(self.localFilter))
        except:
            a.append('No Filters')
        return '\n'.join(a)

    def getInfos(self):
        '''WriteMe'''
        # Variables de travail
        deletion, insertion, self.totCount, ignoreNext, non_starting_deletion, size = 0, 0, 0, 0, 0, 0
        numeric = recompile('[0-9]')
        nucleotides = recompile('[ATGCNatgcn]')
        #print self.sequence, self.quality, self.mapQual
        for item, itemQual, itemMapQual in zip(self.sequence, self.quality, self.mapQual): # 4ATTTTTT
            # Lorsque le charactère '^' est rencontré, on ignore le prochain charactère car il désigne la qualité du read.
            if ignoreNext == 1:
                ignoreNext = 0
                continue
            # Si une insertion ou une délétion est détectée, on marque sa taille dans un dictionnaire et on initialise sa séquence.
            if deletion == 1 or insertion == 1:
                size = int(item)
                next_el = ""
                if (insertion == 1):
                    insertion = 2
                else:
                    deletion = 2
                continue
            # Size a été défini lorsqu'on a vu un InDel. Tant qu'on est dans l'InDel, on ne compte pas les nucléotides. On ajoute également, lorsqu'on a fini
            # de parser l'InDel, la séquence de l'InDel à un dictionnaire dédié.
            if size != 0:
                if numeric.match(item):
                    #tmp_size = str(size)
                    size = '{0}{1}'.format(size, item)
                    size = int(size)
                    continue
                # On détecte si l'InDel est sur le brin sens ou anti-sens
                if deletion == 2 or insertion == 2:
                    if (deletion == 2):
                        if item.isupper():
                            self.forward["del"] += 1
                            self.deletionsPlus[size] += 1
                        else:
                            self.reverse["del"] += 1
                            self.deletionsMinus[size] += 1
                        last = ("del", size)
                    else:
                        if item.isupper():
                            self.forward["ins"] += 1
                            self.insertionsPlus[size] += 1
                        else:
                            self.reverse["ins"] += 1
                            self.insertionsMinus[size] += 1
                        last = ("ins", size)
                    IDSequence = ''
                    deletion = 0
                    insertion = 0
                # Les séquences associées sont stockées séparément et individuellement
                if nucleotides.match(item):
                    IDSequence += item
                    size -= 1
                    if (size == 0):
                        if (last[0] == "del"):
                            self.deletionsInfo[last[1]] += [IDSequence]
                        elif (last[0] == "ins"):
                            self.insertionsInfo[last[1]] += [IDSequence]
                    continue

            # Les charactères '.' ou ',' signifie "référence". On compte ainsi +1 dans la référence.
            if item == "." or item == ",":
                # Sur le brin sens
                if item == ".":
                    self.forward[self.refNuc] += 1
                    self.refCount += 1
                    self.totCount += 1
                    next_el = item
                # Ou le brin anti-sens
                else:
                    self.reverse[self.refNuc] += 1
                    self.refCount += 1
                    self.totCount += 1
                    next_el = item
            # On compte les variants. Une lettre non "." ou "," signifie un nucléotide différent de la référence.
            elif nucleotides.match(item):
                if item.isupper():
                    # Sur le brin sens
                    self.forward[item] += 1
                    next_el = item
                else:
                    # et sur le brin anti-sens si le charactère est minuscule
                    self.reverse[item.upper()] += 1
                    next_el = item
                self.totCount += 1
            # On ignore les informations quant au mapping quality
            elif item == "^":
                ignoreNext = 1
                continue
            # Symbole marquant une délétion.
            elif item == '-':
                deletion = 1
            # Symbole marquant une insertion. Les insertions sont comptées dans le "total depth"
            elif item == '+':
                # On ne doit pas compter le "." et le "," dans le total depth.
                if (next_el == "."):
                    self.forward[self.refNuc] -= 1
                    self.refCount -= 1
                    self.totCount -= 1
                elif (next_el == ","):
                    self.reverse[self.refNuc] -= 1
                    self.refCount -= 1
                    self.totCount -= 1
                else:
                    pass
                insertion = 1
                self.totCount += 1
            # Les non_starting_deletions sont des traces d'une délétion ayant commencée avant le read actuel. Celles-ci ne sont pas comptées ici.
            elif item == '*':
                non_starting_deletion += 1

    # Cette fonction permet de compter les évenements sur le brin sens et anti-sens, ainsi que les InDels
    def computeNucCount(self):
        '''WriteMe'''
        for el in ["A", "T", "C", "G", "N"]:
            nucCount = 0
            nucCount += self.reverse[el]
            nucCount += self.forward[el]
            self.nucCount[el] = nucCount
        # Comptage global sans prendre en compte le caractère sens ou reverse de chaque évenement.
        for i in self.insertionsPlus.itervalues():
            self.nucCount["ins"] += i
        for i in self.deletionsPlus.itervalues():
            self.nucCount["del"] += i
        for i in self.insertionsMinus.itervalues():
            self.nucCount["ins"] += i
        for i in self.deletionsMinus.itervalues():
            self.nucCount["del"] += i

    # Cette fonction permet de calculer la fréquence de chaque évenement.
    def computeFreq(self):
        '''WriteMe'''
        self.computeNucCount()
        try:
            self.refFreq = 100 * float(self.refCount) / float(self.totCount)
        except:
            self.refFreq = 0
        self.globalVarFreq = 100 - self.refFreq
        for el in ["A", "T", "C", "G", "N", "ins", "del"]:
            count = self.nucCount[el]
            try:
                freq = 100 * float(count) / float(self.totCount)
            except:
                freq = 0
            self.frequencies[el] = freq

    # Cette fonction permet de choisir un variant ou de filtrer le variant.
    def call(self):
        '''WriteMe'''
        # Le variant choisi aura la plus haute fréquence.
        #highestFreq = '',0
        #bestFreqDel = '', 0
        #bestFreqIns = '', 0
        #qual = 0
        #arrayFreq = []
        #calledFreq = '', 0
        #sensitivity = 0
        # Fonctionnalité à ajouter.
        self.callFrequencies = {}
        self.variationsList = []

        # Etape 0: preprocesser les indels à cette position.
        #insertions0 = []
        self.insertions = defaultdict(int)
        for e, i in self.insertionsInfo.iteritems():
            i = [x.upper() for x in i]
            for elem in i:
                self.insertions[elem] += 1

        #deletions0 = []
        self.deletions = defaultdict(int)
        for e, i in self.deletionsInfo.iteritems():
            i = [x.upper() for x in i]
            for elem in i:
                self.deletions[elem] += 1

        #1ère étape: créer une liste de tuples récapitulant l'ensemble des variations à la position ainsi que leurs fréquences (sans méta (ins / del))
        for e, i in self.frequencies.iteritems():
            if (e in ["ins", "del"]):
                continue
            self.variationsList.append((e, i))

        bestFreqIns, insList = self.computeIndelFreq('ins')
        bestFreqDel, delList = self.computeIndelFreq('del')

        for el in insList:
            self.variationsList.append(el)
        for el in delList:
            self.variationsList.append(el)

        #2ème étape: sorter self.variationsList par fréquences
        # Format de self.variationsList: [ (A, 44), (T, 10), (+AG, 50) ... ]
        self.sortedVariationsList = [x for x in self.variationsList if ((x[1] != float(0)))]
        self.sortedVariationsList = sorted(self.sortedVariationsList, key=itemgetter(1), reverse=True)
        # Format de self.sortedVariationsList: [ (+AG, 50), (A, 44), (T, 10) ... ]
        # On se retrouve avec un tableau de fréquences n'incluant pas la référence.
        # Tous les allèles ne matchant pas la référence devraient être décrits dans le champ ALT. Leur fréquence et leur profondeur devrait également être
        # Décrite. Néanmoins, il faudra tester chaque allèle pour vérifier qu'elle est utilisable en l'état (filterstate)
        #calledList = []
        stackedFrequencies = []
        tempFrequencies = [x[1] for x in self.sortedVariationsList]
        i = 0
        while i != len(self.sortedVariationsList):
            stackedFrequencies.append(sum(tempFrequencies[i:]))
            i += 1

        stackedAlleles = []
        correctedVariationsList = []
        for el in zip(self.sortedVariationsList, stackedFrequencies):
            correctedVariationsList.append((el[0][0], el[0][1],el[1]))
        for e, i, s in correctedVariationsList:
            a = Allele(self.parameters)
            a.Allele = e
            a.Frequency = i
            a.stackedQuality = s
            if (a.Allele.startswith('-')):
                a.deletion = True
                type = "del"
            elif (a.Allele.startswith('+')):
                a.insertion = True
                type = "ins"
            else:
                type = e

            a.forward = self.forward[type]
            a.reverse = self.reverse[type]
            a.referenceForward = self.forward[self.refNuc]
            a.referenceReverse = self.reverse[self.refNuc]

            a.Depth = float(a.Frequency) * 0.01 * float(self.totCount)
            if (int(a.Depth) < a.Depth):
                a.Depth += 0.000001
            elif (int(a.Depth) > a.Depth):
                a.Depth -= 0.000001
            else:
                pass
            a.Depth = int(a.Depth)

            passing = a.getLocalFilterState()
            a.localFilter = passing
            if (e != self.refNuc):
                stackedAlleles.append(a)
            else:
                self.refAllele = a
        try:
            self.refAllele
        except:
            self.refAllele = Allele(self.parameters)
            self.refAllele.Allele = self.refNuc
            self.refAllele.Frequency = 0
            self.refAllele.Depth = 0

        self.reference = self.refAllele.Allele
        self.stackedAlleles = stackedAlleles
        self.goodAlleles = []
        self.trashedAlleles = []

        for el in self.stackedAlleles:
            if el.localFilter is not True:
                self.trashedAlleles.append(el)
            else:
                self.goodAlleles.append(el)

        if (len(self.goodAlleles) == 0):
            if (len(self.trashedAlleles) == 0): 
                self.discard = True
                self.calledAlleles = (self.refAllele, None, None)
            else:
                self.discard = True
                self.calledAlleles = (self.refAllele, None, self.trashedAlleles)
                self._genotypeVariant(self.refAllele, self.goodAlleles, self.trashedAlleles, True)
            return self.calledAlleles
        else:
            if (len(self.trashedAlleles) == 0):
                self.calledAlleles = (self.refAllele, self.goodAlleles, None)
                self.discard = False
                self._genotypeVariant(self.refAllele, self.goodAlleles, self.trashedAlleles, False)
            else:
                self.calledAlleles = (self.refAllele, self.goodAlleles, self.trashedAlleles)
                self.discard = False
                self._genotypeVariant(self.refAllele, self.goodAlleles, self.trashedAlleles, False)
            return self.calledAlleles

    def _genotypeVariant(self, ref, good, filtered, switch):
        '''WriteMe'''
        def alleleChooser(ref, good):
            '''WriteMe'''
            temp_array = []
            temp_array.append(ref)
            hetTreshold = self.parameters['hetTreshold']
            homTreshold = self.parameters['homTreshold']
            #print ref
            for el in good:
                temp_array.append(el)
                #print el
            sorted_temp = sorted(temp_array, key=attrgetter('Frequency'))
            isHet = [True if x.Frequency > hetTreshold else False for x in sorted_temp]
            isHom = [True if x.Frequency > homTreshold else False for x in sorted_temp]
            markHet = False
            for el, het, hom in zip(sorted_temp, isHet, isHom):
                #print el
                #print het, hom
                if (hom and not markHet):
                    return '{0}/{0}'.format(el.rank)
                elif (het and not markHet):
                    markHet = [str(el.rank)]
                elif (het and markHet):
                    markHet.append(str(el.rank))
                elif (hom and markHet):
                    markHet += [str(el.rank), str(el.rank)]
            if not (markHet):
                return '0/0'
            else:
                markHet = sorted(markHet)
            return '/'.join(markHet)

        def postProcessRefAlt(ref, alts):
            '''WriteMe'''
            # This snippet was provided to me by Jeremie Pagnac - IGR - 2015
            # Yannick Boursin corrected a bug in InDel display
            alts = alts.split(',')
            biggestDel = ''
            for alt in alts:
                if '-' in alt:
                    if len(alt) > len(biggestDel):
                        biggestDel = alt
            new_ref = ref + biggestDel[1:]
            new_alts = []
            for alt in alts:
                if '-' in alt:
                    new_alt = new_ref[:-len(alt) + 1]
                elif '+' in alt:
                    new_alt = new_ref[0] + alt[1:] + new_ref[1:]
                else:
                    new_alt = alt + new_ref[1:]
                new_alts.append(new_alt)
            new_alts = ','.join(new_alts)
            return [new_ref, new_alts]

        if (switch is True):
            if (filtered is not None):
                sortedGoodAlleles = sorted(filtered, key=attrgetter('Frequency'), reverse=True)
            else:
                self.gtFreq = '0/0'
                self.gtDepth = ref.Depth
                self.alt = ref.Allele
                self.calledQuality = 100
                self.calledLocalsb = 0
                self.calledLfs = 0
                return
        else:
            sortedGoodAlleles = sorted(good, key=attrgetter('Frequency'), reverse=True)

        #resultsToReturn = []
        #called0 = ref
        #called0Freq = ref.Frequency
        self.alt = []
        self.gtFreq = []
        self.gtDepth = []
        self.LSB = []
        self.LFS = []

        self.majorAllele = sortedGoodAlleles[0]
        i = 0
        ref.rank = 0
        for el in sortedGoodAlleles:
            el.rank = i + 1
            i += 1
        #On choisi le génotype (ex: 0/1 ou 0/1/2 ou 0/1/1 ...)
        self.calledGenotype = alleleChooser(ref, sortedGoodAlleles)
        #print self.calledGenotype
        AllAlleles = []
        AllAlleles.append(ref)
        for el in sortedGoodAlleles:
            AllAlleles.append(el)

        #sortedFreqs = [x.Frequency for x in sortedGoodAlleles]

        # Maintenant, on expose les variables afin de ne pas créer de problème pour les fonctions d'affichage.
        # Si possible, on rends l'affichage plus simple en préconcaténant l'information.
        self.reference = ref.Allele
        for el in AllAlleles:
            self.gtFreq.append(el.Frequency)
            self.gtDepth.append(el.Depth)
        for el in sortedGoodAlleles:
            self.LSB.append(el.localsb)
            self.LFS.append(el.lfs)
            self.alt.append(el.Allele)

        self.gtFreq = ','.join(['{0:.2f}'.format(x) for x in self.gtFreq])
        self.gtDepth = ','.join([str(int(x)) for x in self.gtDepth])
        self.alt = ','.join(self.alt)
        self.calledQuality = 50
        self.calledLocalsb = max(self.LSB)
        self.calledLfs = max(self.LFS)
        self.oldalt = self.alt
        self.reference, self.alt = postProcessRefAlt(self.reference, self.alt)

        good_depth = 0
        filtered_depth = 0

        All_alleles_even_bad = [ref]
        if (good is not None):
            for el in good:
                All_alleles_even_bad.append(el)
                good_depth += el.Depth
        if (filtered is not None):
            for el in filtered:
                All_alleles_even_bad.append(el)
                filtered_depth += el.Depth

        score_good = 0
        score_filtered = 0
        if good is not None:
            for el in good:
                #print el
                if (str(el.rank) not in self.calledGenotype.split('/') and str(el.Allele) != self.reference):
                    if (str(el.Allele) in self.oldalt.split(',')[:2]):
                        pass
                    elif (str(el.Allele) in self.oldalt.split(',')[2:3]):
                        score_good += 250 * (float(el.Depth) / float(good_depth))
                    else:
                        score_good += 750 * (float(el.Depth) / float(good_depth))
        if filtered is not None:
            for el in filtered:
                score_filtered += 100 * (float(el.Depth) / (float(filtered_depth) + float(good_depth)))

        #other_mutation_rate = frequence de mutation des autres allèles en rapport avec la fréquence de mutation totale
        # Un allèle qui passe un filtre vaut un poids de 5
        # Un allèle filtré vaut un poids de 1
        # La profondeur est à multiclasser: profondeur des bons candidats, profondeur des mauvais candidats. Si le rapport est inférieur à 1, retirer 1.
        # Score de qualité = 100 - Somme(Poids_allèle * Depth_allele/Depth_category)
        self.calledQuality = 100 - score_good - score_filtered
        if (self.calledQuality < 0):
            self.calledQuality = 0

    def getLocalSequenceContext(self):
        """This function intends to return the genomic context of the variant."""
        ## 1) Adéquation des fréquences nucléotidiques.
        context = self.context
        biggerContext = self.biggerContext
        frequencies = {"A": float(context.count("A"))/len(context), "T": float(context.count("T"))/len(context), "G": float(context.count("G"))/len(context), "C": float(context.count("C"))/len(context)}
        biggerFrequencies = {"A": float(biggerContext.count("A"))/len(biggerContext), "T": float(biggerContext.count("T"))/len(biggerContext), "G": float(biggerContext.count("G"))/len(biggerContext), "C": float(biggerContext.count("C"))/len(biggerContext)}

        def getDeltaFrequencies(frequencies, expFrequencies):
            '''WriteMe'''
            return sum([(f - ef)**2 for f, ef in zip(frequencies, expFrequencies)])

        ## 2) Perfect Sequentiality score = LRS simplification to chars

        ## 3) LRS (Longest Repeated Substring) => Implies to make a specific data structure

        def getsubs(loc, s):
            '''WriteMe'''
            substr = s[loc:]
            i = -1
            while(substr):
                yield substr
                substr = s[loc:i]
                i -= 1

        def longestRepetitiveSubstring(r, minocc=3):
            '''WriteMe'''
            occ = defaultdict(int)
            # tally all occurrences of all substrings
            for i in range(len(r)):
                for sub in getsubs(i, r):
                    occ[sub] += 1

            # filter out all substrings with fewer than minocc occurrences
            occ_minocc = [k for k,v in occ.items() if v >= minocc]

            if occ_minocc:
                maxkey = max(occ_minocc, key=len)
                return maxkey, occ[maxkey]
            else:
                raise ValueError("no repetitions of any substring of '%s' with %d or more occurrences" % (r,minocc))

        deltaFreq = getDeltaFrequencies(frequencies)
        l, occ = longestRepetitiveSubstring(context, minocc=1)


    @staticmethod
    def window(seq, k=4):
        '''WriteMe'''
        it = iter(seq)
        result = tuple(islice(it, k))
        if (len(result) == k):
            i = 0
            yield result, i
        for elem in it:
            result = result[1:] + (elem,)
            i += 1
            yield result, i

    @staticmethod
    def suffixWindow(seq, k=4):
        '''WriteMe'''
        it = iter(seq)
        result = tuple(islice(it, k))
        if (len(result) == k):
            i = 0
            yield result, i
        for elem in it:
            pass ## Here should be continued

    def getLocalFilterState(self):
        '''WriteMe'''
        passFilter = True
        self.failedLocFilters = []

        # Defines parameters
        # There is some identical code in Allele Class ... maybe we should dedup ?
        localStrandBias = self.parameters['localStrandBias']
        minFreq = self.parameters['minFreq']
        minReadsAlt = self.parameters['minReadsAlt']
        fisherStrand = self.parameters['fisherStrand']
        minQuality = self.parameters['minQuality']
        # End define parameters

        localsb = float(self.forward[self.type]) * 100 / (float(self.reverse[self.type]) + float(self.forward[self.type]))
        self.localsb = localsb
        reverseLocalStrandBias = abs(100 - localStrandBias)
        paramLocalStrandBias = [localStrandBias, reverseLocalStrandBias]
        lowerParamLocalStrandBias = min(paramLocalStrandBias)
        greaterParamLocalStrandBias = max(paramLocalStrandBias)

        if (float(localsb) > float(greaterParamLocalStrandBias)):
            self.failedLocFilters.append("Local: LocalStrandBias > {0}".format(greaterParamLocalStrandBias))
        if (float(localsb) < float(lowerParamLocalStrandBias)):
            self.failedLocFilters.append("Local: LocalStrandBias < {0}".format(lowerParamLocalStrandBias))

        if (float(self.calledFreq) < float(minFreq)):
            self.failedLocFilters.append("Local: Freq < {0}".format(minFreq))
        if (int(self.calledDepth) < int(minReadsAlt)):
            self.failedLocFilters.append("Local: AltReads < {0}".format(minReadsAlt))
        if (self.markIns):
            allele = "ins"
        elif (self.markDel):
            allele = 'del'
        else:
            allele = self.calledAllele
        self.lfs = self.computeFisherStrandBias(allele)
        if (float(self.lfs) > float(fisherStrand)):
            self.failedLocFilters.append("Local: LocalFisherStrand > {0}".format(fisherStrand))
        if (float(self.calledQuality) < float(minQuality)):
            self.failedLocFilters.append("Quality < {0}".format(minQuality))
        if (len(self.failedLocFilters) != 0):
            passFilter = '; '.join(self.failedLocFilters)
        return passFilter

    def computeIndelFreq(self, indel):
        '''WriteMe'''
        bestFreq = '', 0
        indelList = []
        if (indel == 'ins'):
            for e, i in self.insertions.iteritems():
                freq = 100 * float(i) / float(self.totCount) + float(0.000000001)
                indelList.append(('+{0}'.format(str(e)), freq))
                if (freq > bestFreq[1]):
                    bestFreq = str(e), float(freq)
            bestFreq = '+{0}'.format(bestFreq[0]), bestFreq[1] 
        elif (indel == 'del'):
            for e, i in self.deletions.iteritems():
                freq = 100 * float(i) / float(self.totCount) + float(0.000000001)
                indelList.append(('-{0}'.format(str(e)), freq))
                if (freq > bestFreq[1]):
                    bestFreq = str(e), float(freq)
            bestFreq = '-{0}'.format(bestFreq[0]), bestFreq[1]
        else:
            return False
        return bestFreq, indelList

    # Cette fonction calcule l'équilibre entre le brin sens et non sens. Elle donne une mesure de strand bias.
    def computeStrandBias(self):
        '''WriteMe'''
        self.fw = 0
        self.rv = 0
        for el in ["A", "T", "C", "G", "N"]:
            self.fw += float(self.forward[el])
            self.rv += float(self.reverse[el])
        # Déséquilibre en faveur du brin sens ?
        self.fw += self.forward["ins"]
        self.fw += self.forward["del"]
        self.rv += self.reverse["ins"]
        self.rv += self.reverse["del"]
        try:
            ratio1 = (100 * self.fw) / (self.rv + self.fw)
        except:
            ratio1 = 0
        return ratio1

    def titv(self, alt=None):
        '''WriteMe'''
        ref = self.reference
        if alt is None:
            alt = self.majorAllele.Allele
        #print ref, alt, len(ref), len(alt)
        transition = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]
        if (ref, alt) in transition:
            return 'ti'
        elif (len(ref.lstrip('-').lstrip('+')) > 1):
            return 'del'
        elif (len(alt.lstrip('-').lstrip('+')) > 1):
            return 'ins'
        else:
            return 'tv'

    def statsAllAllele(self, type="good"):
        '''WriteMe'''
        titvAll = []
        if type == "Good":
            for allele in self.goodAlleles:
                titvAll.append(self.titv(allele.Allele))
        elif type == "all":
            for allele in self.goodAlleles:
                titvAll.append(self.titv(allele.Allele))
            for allele in self.trashedAlleles:
                titvAll.append(self.titv(allele.Allele))
        elif type == "Trash":
            for allele in self.trashedAlleles:
                titvAll.append(self.titv(allele.Allele))
        
        return titvAll

    # Cette fonction effectue le test de Fisher exact entre la référence et les alternatifs. Elle donne une mesure de strand bias.
    def computeFisherStrandBias(self, allele=""):
        '''WriteMe'''
        # On créé la table de contingence
        #print allele
        if (allele != ""):
            localref = self.rawrefnuc
        else:
            localref = self.refNuc
        forwardRef = self.forward[localref]
        reverseRef = self.reverse[localref]
        items = ["A", "T", "G", "C", "N", "del", "ins"]
        items.remove(localref)
        forwardAlt = 0
        reverseAlt = 0
        if (allele == ""):
            for el in items:
                forwardAlt += float(self.forward[el])
                reverseAlt += float(self.reverse[el])
        else:
            forwardAlt = float(self.forward[allele])
            reverseAlt = float(self.reverse[allele])
        # On lance le test proprement dit
        #print "fwr: {0}, rvr: {1}, fwa: {2}, rva: {3}".format(forwardRef, reverseRef, forwardAlt, reverseAlt)
        return fi.computeFisher(forwardRef, reverseRef, forwardAlt, reverseAlt)

    def isVariant(self):
        '''WriteMe'''
        self.computeFreq()
        if self.globalVarFreq == 0 or int(self.totCount == 0):
            if float(self.frequencies["del"]) == 0:
                return False
        return True

    # Cette fonction permet de définir l'état filtré ou présent du variant. Tout nouveau filtre doit y être enregistré.
    def getGlobalFilterState(self):
        '''WriteMe'''
        # Define Parameters
        minFreq = self.parameters['minFreq']
        minDepth = self.parameters['minDepth']
        strandBias = self.parameters['strandBias']
        fisherStrand = self.parameters['fisherStrand']

        self.computeFreq()
        passFilter = True
        self.failedGlobFilters = []
        if ((self.globalVarFreq == 0) or (int(self.totCount) == 0)):
            if (float(self.frequencies["del"]) == 0):
                self.failedGlobFilters.append("NoVariant")
        if (float(self.globalVarFreq) < float(minFreq)):
            if (float(self.frequencies["del"]) <= float(minFreq)):
                self.failedGlobFilters.append("Global: VarFreq < {0}".format(minFreq))
            else: pass
        if (int(self.totCount) < int(minDepth)):
            self.failedGlobFilters.append("Global: DP < {0}".format(minDepth))

        self.sb = self.computeStrandBias()

        reverseStrandBias = abs(100 - float(strandBias))
        paramStrandBias = [float(strandBias), float(reverseStrandBias)]
        lowerStrandBias = min(paramStrandBias)
        greaterStrandBias = max(paramStrandBias)

        if (float(self.sb) > float(greaterStrandBias)):
            self.failedGlobFilters.append("Global: SB > {0}".format(greaterStrandBias)) 
        if (float(self.sb) < float(lowerStrandBias)):
            self.failedGlobFilters.append("Global: SB < {0}".format(lowerStrandBias))

        self.gfs = self.computeFisherStrandBias()
        if (float(self.gfs) > float(fisherStrand)):
            self.failedGlobFilters.append("Global: GlobalFisherStrand > {0}".format(fisherStrand))
        if (len(self.failedGlobFilters) != 0):
            passFilter = '; '.join(self.failedGlobFilters)
        self.globalFilter = passFilter
        return passFilter

    # Les fonctions suivantes servent à afficher des informations.
    # Fonction d'affichage des informations de base.
    def getBaseInformations(self):
        '''WriteMe'''
        return "{0}\t{1}\t{2}\t{3}".format(self.chr, self.start, self.refNuc, self.totCount)

    # Fonction d'affichage des fréquences de chaque nucléotide.
    def getFrequencies(self):
        '''WriteMe'''
        return "{0:.2f}%[{1}]\t{2:.2f}%[{3}]\t{4:.2f}%[{5}]\t{6:.2f}%[{7}]\t{8:.2f}%[{9}]".format(self.frequencies["A"],
                                                                                            self.nucCount["A"],
                                                                                            self.frequencies["T"],
                                                                                            self.nucCount["T"],
                                                                                            self.frequencies["G"],
                                                                                            self.nucCount["G"],
                                                                                            self.frequencies["C"],
                                                                                            self.nucCount["C"],
                                                                                            self.frequencies["N"],
                                                                                            self.nucCount["N"])

    # Fonction d'affichage des InDels
    def getInDels(self):
        '''WriteMe'''
        toReturn = ""
        deletions0 = []
        insertions0 = []
        self.deletions = defaultdict(int)
        self.insertions = defaultdict(int)
        # Pour chaque délétion de taille x, on créé un tuple contenant la taille, le nombre de reads supportant la délétion, ainsi que la séquence de la délétion.
        for e, i in self.deletionsInfo.iteritems():
            i = [x.upper() for x in i]
            for elem in i:
                self.deletions[elem] += 1

        for e, i in self.insertionsInfo.iteritems():
            i = [x.upper() for x in i]
            for elem in i:
                self.insertions[elem] += 1

        for e, i in self.deletions.iteritems():
            deletions0.append(("size: " + str(len(e)), "number: " + str(i), "sequence: " + str(e)))

        # Idem
        for e, i in self.insertions.iteritems():
            insertions0.append(("size: " + str(len(e)), "number: " + str(i), "sequence: " + str(e)))

        # On retourne None ou l'information selon la présence ou l'absence d'InDel
        if (len(deletions0) != 0) and (len(insertions0) != 0):
            toReturn = toReturn + "DELETIONS: {2:.2f}% | {0}\tINSERTIONS: {3:.2f}% | {1}".format(deletions0, insertions0, self.frequencies["del"], self.frequencies["ins"])
        elif (len(deletions0) != 0):
            toReturn = toReturn + "DELETIONS: {1:.2f}% | {0}\tNone".format(str(deletions0), self.frequencies["del"])
        elif (len(insertions0) != 0):
            toReturn = toReturn + "None\tINSERTIONS: {1:.2f}% | {0}".format(str(insertions0), self.frequencies["ins"])
        else:
            toReturn = toReturn + "None\tNone"
        return toReturn
    
