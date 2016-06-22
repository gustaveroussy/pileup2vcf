#!/usr/bin/env python
# encoding: utf-8

'''
Created on 3 août 2015

@author: boursin
'''

from matplotlib import pyplot as plt
from matplotlib import collections as mc
from Bio.Statistics.lowess import lowess
from collections import OrderedDict, deque
from GenoPy.Genomic import *
from GenoPy.Variants import *
#from GenoPy.Segmentation import pelt
import pyfaidx
import sys


###############################################################################
###########################  MPILEUP Class  ###################################
###############################################################################

class MPileup(object):
    '''WriteMe'''
    def __init__(self, file, dph, io, parameters, bed=None):
        self.log = io.log
        self.parameters = parameters
        self.bed = bed
        self.file = file # io.getIO('mpileup')
        self.dph = dph
        self.arrayPileup = []
        
        self.wroteHeader = False
        self.fasta = pyfaidx.Fasta(parameters['refGenome'])
        self.indelsToAddToVCF = []

    def getContextFromPosition(self, chr, pos, howMuchNucleotides=4):
        '''Returns flanking nucleotides for coordinates'''
        pos = int(pos)
        boundaries = (pos - howMuchNucleotides, pos + howMuchNucleotides)
        #if (boundaries[0] not in self.fasta[chr] or boundaries[1] not in self.fasta[chr]):
        #    return None
        return self.fasta[chr][boundaries[0]:boundaries[1]]

    # Fonction permettant de parser le fichier mpileup et de définir pour chaque ligne un objet variant. Cet objet est ammené à être callé ou non lors de l'étape
    # getFilterState. Si passFilter = cur.getFilterState() = True, alors on a affaire à un variant.
    def parse(self):
        '''This takes as input the mpileup file and parses it. It has disparate functions, but this is the data provider function.
        * Returns: self.arrayPileup
        '''
        avancement = []
        alreadyPrinted = False
        lineCounter = 0
        for line in self.file:
            lineCounter += 1
            try:
                elem = line.split('\t')
                chr = elem[0].strip().rstrip('\n').rstrip('\t')
                pos = int(elem[1].strip().rstrip('\n').rstrip('\t'))
                refNuc = elem[2].strip().rstrip('\n').rstrip('\t').upper()
                depth = int(elem[3].strip().rstrip('\n').rstrip('\t'))
                sequence = elem[4].strip().rstrip('\n').rstrip('\t')
                quality = elem[5].strip().rstrip('\n').rstrip('\t')
                mapQual = elem[6].strip().rstrip('\n').rstrip('\t')
                posInRead = elem[7].strip().rstrip('\n').rstrip('\t').split(',')
                #print len(posInRead), len(mapQual), len(quality)
                # Assertions
                assert (len(mapQual) == len(quality)), "Weird formating: mapping quality track length is not the same as nucleotide quality track length"
                assert (len(quality) == len(posInRead)), "Weird formating: nucleotide quality track length is not the same as position in read track length"
                assert (len(refNuc) == 1), "Weird formating: there is {} nucleotide at this position !".format(len(refNuc))
            
            except:
                elem = line.split('\t')
                print elem
                sequence = elem[4].strip().rstrip('\n').rstrip('\t')
                quality = elem[5].strip().rstrip('\n').rstrip('\t')
                print quality, len(quality)
                mapQual = elem[6].strip().rstrip('\n').rstrip('\t')
                print mapQual, len(mapQual)
                posInRead = elem[7].strip().rstrip('\n').rstrip('\t').split(',')
                print posInRead, len(posInRead)
                print len(quality), len(mapQual), len(posInRead), len(sequence)
                raise MpileupFormatError("Cannot initialize at line {}.\n Please check that your Mpileup file has 8 columns and has been generated with mapping quality and base position in read using samtools 0.1.18 (further versions seems to suffer from bugs).\n MPileup with good format might be generated using: samtools mpileup -A -s -O -B -o test.mpileup -f ../Pileup2VCF/hg19/hg19.fa 208204422-ADN-2_S2_L001.bam".format(lineCounter))
            if depth == 0:
                #print "Unsequenced region at {}-{}".format(chr, pos)
                continue
            if not chr in avancement:
                self.dph.walk(Chromosome=chr)
                avancement.append(chr)
                before = elem[1]
                self.log("Parsing and calling on {0}".format(chr))
            if (lineCounter % 500 == 0):
                percents = float(lineCounter) * 100 / float(self.parameters['lineNb'])
                self.log("Parsed around {0:.2f}% of the file ({1} line read / {2} total)".format(percents, lineCounter, self.parameters['lineNb']), loglevel=3)
            context = self.getContextFromPosition(chr, pos)
            biggerContext = self.getContextFromPosition(chr, pos, howMuchNucleotides=100)
            if (int(depth) > 0):
                self.dph.add(chr, pos, Position(chr, pos, quality, mapQual, depth, refNuc, pos_before=before, bed=self.bed))
                before = pos
            # Création d'un object contenant la ligne et pouvant accueillir un éventuel variant
            cur = Variant(chr, pos, refNuc, depth, sequence, quality, mapQual, self.parameters, context, biggerContext)
            #print cur
            cur.getInfos()
            if cur.isVariant() is True:
                #print cur
                self.arrayPileup.append(cur)
            # A ce stade, on enregistre tous les objets variants en se basant uniquement sur la présence éventuelle d'un nucléotide non ref. Le problème de cette méthode est qu'elle
            # peut entrainer beaucoup de variants selon la qualité du séquençage...
            # On verra bien ce que cela donne. L'étape suivante est
        self.dph.writeDepthProfile() # We write the last positions in depthProfile before returning
        return UnclassifiedVariantCollection(self.arrayPileup, self.parameters)
            # # Le variant éventuel passe-t-il les filtres ?
            # #passFilter = cur.getGlobalFilterState()
            # cur.globalFilter = passFilter
            # #passFilter = "NoVariant;"
            # if passFilter is True:
            #     ####### CORRIGER A PARTIR DE LA ... #####
            #     call = cur.call()
            #     #Il y a 3 cas globaux de figure:
            #     # -Soit on a l'allèle de référence, suivi de None et None. En ce cas, il n'y a aucun variant à la position.
            #     #      Ce cas ne devrait pas se présenter vu qu'on a déjà passé le filtre global. On sait donc qu'il n'est qu'utopique et appartient à la gestion d'erreur
            #     # - Soit on a l'allèle de référence, suivi d'un tableau d'allèles alternatifs, suivi de None
            #     #      Ce cas arrive dans le cas où de bons alternatifs ont été trouvés, et aucun allèle n'a été rejeté par les filtres.
            #     #      Il faut donc encore effectuer l'étape de genotyping
            #     # - Soit on a l'allèle de référence, suivi d'un tableau d'allèles alternatifs, suivi d'un tableau d'allèles rejetés
            #     #      Dans ce cas, il faut également effectuer l'étape de genotyping sur les bons alternatifs, et voir si on reporte ou pas les cas où les filtres n'ont pas été passés.
            #     if (call[1] is None and call[2] is None):
            #         self.log("There's something fishy around here ! Code: A7B9", loglevel=2)
            #     elif (call[1] is not None and call[2] is None):
            #         alreadyPrinted = False
            #         cur.genotypeVariant(call[0], call[1], call[2], False)
            #         cur.discarded = "None"
            #         self.called += 1
            #         self.total += 1
            #         self.arrayPileup.append(cur)
            #         if (cur.titv() == "ti"):
            #             self.ti += 1
            #         elif (cur.titv() == "tv"):
            #             self.tv += 1
            #         elif (cur.titv() == "ins"):
            #             self.insertionNb += 1
            #         elif (cur.titv() == "del"):
            #             self.deletionNb += 1
            #         else:
            #             self.log("There's something fishy around here ! Code: G4A2", loglevel=2)
            #         #Now statistics on all good alleles:
            #         allAllelesStats = cur.statsAllAllele()
            #         self.calledAlleleNumber += len(allAllelesStats)
            #         for el in allAllelesStats:
            #             if el == "ti":
            #                 self.allti += 1
            #             elif el == "tv":
            #                 self.alltv += 1
            #             elif el == "ins":
            #                 self.allInsNb += 1
            #             elif el == "del":
            #                 self.allDelNb += 1
            #             else:
            #                 self.log("There's something fishy around here ! Code: G4A2", loglevel=2)
            #     elif (call[1] is not None and call[2] is not None):
            #         alreadyPrinted = False
            #         cur.genotypeVariant(call[0], call[1], call[2], False)
            #         cur.discarded = True #cur.getDiscardedString(call[2])
            #         self.called += 1
            #         self.total += 1
            #         self.arrayPileup.append(cur)
            #         if (cur.titv() == "ti"):
            #             self.ti += 1
            #         elif (cur.titv() == "tv"):
            #             self.tv += 1
            #         elif (cur.titv() == "ins"):
            #             self.insertionNb += 1
            #         elif (cur.titv() == "del"):
            #             self.deletionNb += 1
            #         else:
            #             self.log("There's something fishy around here ! Code: T7E2", loglevel=2)
            #         #Now statistics on all good alleles:
            #         allAllelesStats = cur.statsAllAllele()
            #         self.calledAlleleNumber += len(allAllelesStats)
            #         for el in allAllelesStats:
            #             if el == "ti":
            #                 self.allti += 1
            #             elif el == "tv":
            #                 self.alltv += 1
            #             elif el == "ins":
            #                 self.allInsNb += 1
            #             elif el == "del":
            #                 self.allDelNb +=1
            #             else:
            #                 self.log("There's something fishy around here ! Code: G4A2", loglevel=2)
            #     elif (call[1] is None and call[2] is not None):
            #         # Cas où toutes les variations ont été rejetées
            #         alreadyPrinted = False
            #         self.trashed += 1
            #         self.total += 1
            #         cur.genotypeVariant(call[0], call[1], call[2], True)
            #         self.discardedVariants.append(cur)
            #     self.depthProfile[chr][pos].nbAltPassFilter = len(call[1]) if call[1] is not None else 0
            #     self.depthProfile[chr][pos].nbAltFiltered = len(call[2]) if call[2] is not None else 0
            #     self.depthProfile[chr][pos].nbTotalAlt = self.depthProfile[chr][pos].nbAltFiltered + self.depthProfile[chr][pos].nbAltPassFilter
            # else:
            #     if "NoVariant" in passFilter.split(';'):
            #         pass
            #     else:
            #         alreadyPrinted = False
            #         self.trashed += 1
            #         self.total += 1
            #         call = cur.call()
            #         cur.genotypeVariant(call[0], call[1], call[2], True)
            #         self.discardedVariants.append(cur)
            # line = self.file.readline()
            # self.totalPositions += 1
            # #self.log(self.totalPosition)
            # if (((self.total % 1000) == 0) and (self.total != 0) and not (alreadyPrinted)):
            #     alreadyPrinted = True
            #     if self.total == 1000:
            #         self.printVCF(self.goodVariants, printHeader=True)
            #         self.log("Saving " + str(len(self.arrayPileup)) + " variants ! Code: Z1D2", loglevel=2)
            #         self.writeDepthProfile()
            #         del self.depthProfile
            #         self.depthProfile = OrderedDict()
            #         self.depthProfile[chr] = OrderedDict()
            #         self.arrayPileup = []
            #         if (self.parameters['outputTrash'] != None):
            #             self.outputTrash(self.trashVariants, printHeader=True)
            #             self.discardedVariants = []
            #     else:
            #         self.printVCF(self.goodVariants)
            #         self.log("Saving " + str(len(self.arrayPileup)) + " variants ! Code: Q4A1", loglevel=2)
            #         self.writeDepthProfile()
            #         del self.depthProfile
            #         self.depthProfile = OrderedDict()
            #         self.depthProfile[chr] = OrderedDict()
            #         self.arrayPileup = []
            #         if (self.parameters['outputTrash'] != None):
            #             self.outputTrash(self.trashVariants)
            #             self.discardedVariants = []


    



    

    
