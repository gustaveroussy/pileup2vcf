#!/usr/bin/env python
# encoding: utf-8

'''
Created on 23 november 2016
Last Update: 23 november 2016

@author: LEBEURRIER Manuel
@contact: manuel.lebeurrier@gustaveroussy.fr
@license: GNU GPLv3
@organization: Gustave Roussy
@version: 1.1
@todo: Apply on real data
'''

from operator import itemgetter, attrgetter
from argparse import ArgumentParser, SUPPRESS, RawDescriptionHelpFormatter
from collections import defaultdict, OrderedDict
import gzip
import sys
from decimal import *
from array import array
import numpy as np
from time import clock

from GenoPy.Genomic import *
from GenoPy.Variants import *
from General.IO import *
from GenoPy.MPileup import *
from GenoPy.InDels import *
from GenoPy.DepthProfile import DepthProfile
from GenoPy.Bins import *
from GenoPy.VCF import *
import gc
from statsmodels.stats.multitest import multipletests
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
stats = importr('stats')


def somaticMod(parameters, args, defaults):
    '''This takes as input the parameters, the arguments and the default parameters to compute the somatic mod.
    * Compute: for the tumoral mpileup
    * Compute: for the normal mpileup
    * Compare: variants between tumoral and normal
    '''
    tumoralProceed(parameters, args, defaults)
    NormalProceed(parameters, args, defaults)
    comparing(parameters, args, defaults)


def tumoralProceed(parameters, args, defaults):
    '''This takes as input the parameters, the arguments and the default parameters to compute the tumoral variant analysis.
    '''
    parameters['prefix'] = args.prefix
    if parameters['prefix'] is not None:
        parameters['output'] = parameters['prefix'] + ".T.GoodVariants.vcf"
        parameters['statistics'] = parameters['prefix'] + ".T.statistics"
        parameters['outputBins'] = parameters['prefix'] + ".T.bins"
        parameters['depthProfile'] = parameters['prefix'] + ".T.depthprofile"
        parameters['outputTrash'] = parameters['prefix'] + ".T.Trash.vcf"
        parameters['regionsw'] = parameters['prefix'] + ".T.bed"

    parameters['inputTumoral'] = args.inputTumoral
    parameters['minBinSize'] = int(args.minBinSize) if args.minBinSize is not None else int(defaults['minBinSize'])
    parameters['maxBinSize'] = int(args.maxBinSize) if args.maxBinSize is not None else int(defaults['maxBinSize'])
    parameters['regions'] = args.regions
    parameters['sTumoralName'] = args.inputTumoral.rstrip('.mpileup') if args.inputTumoral.endswith('.mpileup') else args.inputTumoral.rstrip('.pileup')
    parameters['onlyBed'] = True
    parameters['ll'] = int(args.ll) if args.ll is not None else int(defaults['ll'])
    parameters['refGenome'] = args.refGenome
    parameters['maxDeltaDelta'] = int(args.maxDeltaDelta) if args.maxDeltaDelta is not None else int(defaults['maxDeltaDelta'])

    parameters['IndelWindowLength'] = args.indelWindowLength if args.indelWindowLength is not None else defaults['indelWindowLength']
    parameters['DisableQualityScoreCheck'] = args.disableQualityScoreCheck if args.disableQualityScoreCheck is not None else defaults['disableQualityScoreCheck']

    if args.overrideInput: parameters['inputTumoral'] = args.overrideInput
    if args.overrideOutput: parameters['prefix'] = args.overrideOutput
    if args.silent: parameters['ll'] = -1

    io = IO(parameters['ll'])
    parameters['io'] = io
    io.log("Somatic mod : Tumoral", loglevel=1)
    io.log("Command line: {}".format(' '.join(sys.argv)), loglevel=1)
    io.register("mpileupCount", parameters['inputTumoral'], "r")
    getLinesH = io.getIO("mpileupCount")

    lineNb = getLineCount(getLinesH)

    parameters['lineNb'] = lineNb
    io.log(parameters, loglevel=2)

    io.giveupIO("mpileupCount")

    io.register("mpileup", parameters['inputTumoral'], "r")
    io.register("binningW", parameters['outputBins'], "w")
    io.register("depthProfile", parameters['depthProfile'] + ".withoutMut", "w")
    io.register("statistics", parameters['statistics'], "w")
    io.register("GoodVariants", parameters['output'], "w")
    io.register("TrashedVariants", parameters['outputTrash'], 'w')
    io.register('regions', parameters['regions'], 'r')
    io.register("bedW", parameters['regionsw'], "w")

    # Instanciation de Bed et parsing du fichier.
    bed = Bed()
    bedIn = io.getIO('regions')
    bed.parse(bedIn)
    io.giveupIO('regions')

    # Instanciation de MPileup et parsing du fichier.
    pileuph = io.getIO('mpileup')
    depthProfileTumoral = DepthProfile(io.getIO('depthProfile'), io)
    pileup = MPileup(pileuph, depthProfileTumoral, io, parameters, bed=bed)

    if(args.windowKAll):
        io.register("WindowKAllCov", parameters['prefix'] + ".T.All.WindowK.cov.bed", "w") #Fichier pour la couverture
        io.register("WindowKAllGroup", parameters['prefix'] + ".T.All.WindowK.group", "w") #Fichier pour l'équivalence des groupes 
        unclassified_variants = pileup.parseAndWindowK(True, False, False, io.getIO('WindowKAllCov'), io.getIO('WindowKAllGroup'), parameters['sizeK'], parameters['clusterN'], True)
        io.giveupIO('WindowKAllCov')
        io.giveupIO('WindowKAllGroup')
    else:
        unclassified_variants = pileup.parse(True,False,False)
    good, trashed = unclassified_variants.classify_variants().getCollections()

    AllCollections = AllClassifiedVariantCollections(parameters)
    AllCollections.add('Good', good)
    AllCollections.add('Trash', trashed)

    # We will parse depthProfileTumoral and write it again to integrate variants (sadly, I have no other choice)
    io.register('depthProfileAddMutW', parameters['depthProfile'], "w")
    io.register('depthProfileAddMutR', parameters['depthProfile'] + ".withoutMut", "r")

    depthProfileTumoral.addVariants(unclassified_variants)

    io.giveupIO('mpileup')
    io.unregister('depthProfileAddMutR')
    io.unregister('depthProfileAddMutW')

    io.register("depthProfileR1", parameters['depthProfile'], "r")
    io.register("depthProfileR2", parameters['depthProfile'], "r")
    indels = scanForIndels(parameters)
    allIndels = InDelDetector(indels)
    allIndels = allIndels.extract().genotype(parameters).getIndels()
    readDepthProfile = io.getIO('depthProfileR1')
    binningW = io.getIO('binningW')
    binning = Binning(binningW, parameters)
    binning.bin(readDepthProfile)

    #pileup.addIndelsToVCF()
    io.giveupIO('binningW')
    io.giveupIO('depthProfileR1')
    io.giveupIO('depthProfileR2')

    goodVCFTumoral = VCF(good, parameters, indelsToAdd=allIndels)
    if(args.windowKGood):
        io.register("WindowKGoodCov", parameters['prefix'] + ".T.GoodVariants.WindowK.cov.bed", "w")
        io.register("WindowKGoodGroup", parameters['prefix'] + ".T.GoodVariants.WindowK.group", "w")
        goodVCFTumoral.printWindowK(io.getIO('WindowKGoodCov'), io.getIO('WindowKGoodGroup'), parameters['sizeK'], True)
        io.giveupIO('WindowKGoodCov')
        io.giveupIO('WindowKGoodGroup')

    goodVCFTumoral.printVCF(io.getIO('GoodVariants'), printHeader=True)
    #pileup.printVCF(goodVariants)
    io.giveupIO('GoodVariants')
    trashVCFTumoral = VCF(trashed, parameters, indelsToAdd=allIndels)
    trashVCFTumoral.printVCF(io.getIO('TrashedVariants'), printHeader=True)
    io.giveupIO('TrashedVariants')

    stath = io.getIO('statistics')
    AllCollections.do_stats(stath)
    io.giveupIO('statistics')

    bedW = io.getIO('bedW')
    bed.write(bedW)
    io.giveupIO('bedW')
    io.endlog()

    #Memory clean up
    del bed
    del bedIn
    del stath
    del indels
    del allIndels
    del depthProfileTumoral
    del readDepthProfile
    del binningW
    del binning
    del pileuph
    del unclassified_variants
    del trashed
    del AllCollections
    del trashVCFTumoral
    gc.collect()

##########

def NormalProceed(parameters, args, defaults):
    '''This takes as input the parameters, the arguments and the default parameters to compute the normal variant analysis.
    '''
    parameters['prefix'] = args.prefix
    if parameters['prefix'] is not None:
        parameters['output'] = parameters['prefix'] + ".N.GoodVariants.vcf"
        parameters['statistics'] = parameters['prefix'] + ".N.statistics"
        parameters['outputBins'] = parameters['prefix'] + ".N.bins"
        parameters['depthProfile'] = parameters['prefix'] + ".N.depthprofile"
        parameters['outputTrash'] = parameters['prefix'] + ".N.Trash.vcf"
        parameters['regionsw'] = parameters['prefix'] + ".N.bed"

    parameters['inputNormal'] = args.inputNormal
    parameters['minBinSize'] = int(args.minBinSize) if args.minBinSize is not None else int(defaults['minBinSize'])
    parameters['maxBinSize'] = int(args.maxBinSize) if args.maxBinSize is not None else int(defaults['maxBinSize'])
    parameters['regions'] = args.regions
    parameters['sNormalName'] = args.inputNormal.rstrip('.mpileup') if args.inputNormal.endswith('.mpileup') else args.inputNormal.rstrip('.pileup')
    parameters['onlyBed'] = True
    parameters['ll'] = int(args.ll) if args.ll is not None else int(defaults['ll'])
    parameters['refGenome'] = args.refGenome
    parameters['maxDeltaDelta'] = int(args.maxDeltaDelta) if args.maxDeltaDelta is not None else int(defaults['maxDeltaDelta'])

    parameters['IndelWindowLength'] = args.indelWindowLength if args.indelWindowLength is not None else defaults['indelWindowLength']
    parameters['DisableQualityScoreCheck'] = args.disableQualityScoreCheck if args.disableQualityScoreCheck is not None else defaults['disableQualityScoreCheck']

    if args.overrideInput: parameters['inputNormal'] = args.overrideInput
    if args.overrideOutput: parameters['prefix'] = args.overrideOutput
    if args.silent: parameters['ll'] = -1

    io = IO(parameters['ll'])
    parameters['io'] = io
    io.log("Somatic mod : Normal", loglevel=1)
    io.log("Command line: {}".format(' '.join(sys.argv)), loglevel=1)
    io.register("mpileupCount", parameters['inputNormal'], "r")
    getLinesH = io.getIO("mpileupCount")

    lineNb = getLineCount(getLinesH)

    parameters['lineNb'] = lineNb
    io.log(parameters, loglevel=2)

    io.giveupIO("mpileupCount")

    io.register("mpileup", parameters['inputNormal'], "r")
    io.register("binningW", parameters['outputBins'], "w")
    io.register("depthProfile", parameters['depthProfile'] + ".withoutMut", "w")
    io.register("statistics", parameters['statistics'], "w")
    io.register("GoodVariants", parameters['output'], "w")
    io.register("TrashedVariants", parameters['outputTrash'], 'w')
    io.register('regions', parameters['regions'], 'r')
    io.register("bedW", parameters['regionsw'], "w")

    # Instanciation de Bed et parsing du fichier.
    bed = Bed()
    bedIn = io.getIO('regions')
    bed.parse(bedIn)
    io.giveupIO('regions')

    # Instanciation de MPileup et parsing du fichier.
    pileuph = io.getIO('mpileup')
    depthProfileNormal = DepthProfile(io.getIO('depthProfile'), io)
    pileup = MPileup(pileuph, depthProfileNormal, io, parameters, bed=bed)

    if(args.windowKAll):
        io.register("WindowKAllCov", parameters['prefix'] + ".N.All.WindowK.cov.bed", "w") #Fichier pour la couverture
        io.register("WindowKAllGroup", parameters['prefix'] + ".N.All.WindowK.group", "w") #Fichier pour l'équivalence des groupes 
        unclassified_variants = pileup.parseAndWindowK(True, True, False, io.getIO('WindowKAllCov'), io.getIO('WindowKAllGroup'), parameters['sizeK'], parameters['clusterN'], True)
        io.giveupIO('WindowKAllCov')
        io.giveupIO('WindowKAllGroup')
    else:
        unclassified_variants = pileup.parse(True,True,False)
    good, trashed = unclassified_variants.classify_variants().getCollections()
    
    AllCollections = AllClassifiedVariantCollections(parameters)
    AllCollections.add('Good', good)
    AllCollections.add('Trash', trashed)
    # We will parse depthProfileNormal and write it again to integrate variants (sadly, I have no other choice)

    io.register('depthProfileAddMutW', parameters['depthProfile'], "w")
    io.register('depthProfileAddMutR', parameters['depthProfile'] + ".withoutMut", "r")

    depthProfileNormal.addVariants(unclassified_variants)

    io.giveupIO('mpileup')
    io.unregister('depthProfileAddMutR')
    io.unregister('depthProfileAddMutW')

    io.register("depthProfileR1", parameters['depthProfile'], "r")
    io.register("depthProfileR2", parameters['depthProfile'], "r")
    indels = scanForIndels(parameters)
    allIndels = InDelDetector(indels)
    allIndels = allIndels.extract().genotype(parameters).getIndels()
    readDepthProfile = io.getIO('depthProfileR1')
    binningW = io.getIO('binningW')
    binning = Binning(binningW, parameters)
    binning.bin(readDepthProfile)

    #pileup.addIndelsToVCF()
    io.giveupIO('binningW')
    io.giveupIO('depthProfileR1')
    io.giveupIO('depthProfileR2')

    goodVCFNormal = VCF(good, parameters, indelsToAdd=allIndels)
    if(args.windowKGood):
        io.register("WindowKGoodCov", parameters['prefix'] + ".N.GoodVariants.WindowK.cov.bed", "w")
        io.register("WindowKGoodGroup", parameters['prefix'] + ".N.GoodVariants.WindowK.group", "w")
        goodVCFNormal.printWindowK(io.getIO('WindowKGoodCov'), io.getIO('WindowKGoodGroup'), parameters['sizeK'], True)
        io.giveupIO('WindowKGoodCov')
        io.giveupIO('WindowKGoodGroup')

    goodVCFNormal.printVCF(io.getIO('GoodVariants'), printHeader=True)
    #pileup.printVCF(goodVariants)
    io.giveupIO('GoodVariants')
    trashVCFNormal = VCF(trashed, parameters, indelsToAdd=allIndels)
    trashVCFNormal.printVCF(io.getIO('TrashedVariants'), printHeader=True)
    io.giveupIO('TrashedVariants')

    stath = io.getIO('statistics')
    AllCollections.do_stats(stath)
    io.giveupIO('statistics')

    bedW = io.getIO('bedW')
    bed.write(bedW)
    io.giveupIO('bedW')
    io.endlog()

    #Memory clean up
    del trashVCFNormal
    del bed
    del bedIn
    del stath
    del indels
    del allIndels
    del readDepthProfile
    del binningW
    del binning
    del pileuph
    del unclassified_variants
    del trashed
    del AllCollections
    gc.collect()


def comparing(parameters, args, defaults):
    '''This takes as input the parameters, the arguments and the default parameters to compare results from the tumoral and the normal variant analyses.
    This fonction gets back variant depthProfiles for both tumoral and normal analyses.
    These counts are compared using a fisher test on the corresponding contingency table.
    The pvalues are stored in order to perform a mutiple test correction.
    * Write: prints somatic in *.NT.somatic (even if below the alphaFDR) (all tumoral variants having different variant in normal and all tumoral variants having no variant in normal)
    * Write: prints germline in *.NT.germline (all tumoral variants having the same variants in normal)
    '''
    io = IO(parameters['ll'])
    io.log("Somatic analysis", loglevel=1)

    #Lists for string split
    goodNormalFreq = []
    goodTumoralFreq = []

    #Lists for pvalues and types (=somatic/N.variant!=T.variant) to compute multiple test correction
    pvals = []
    pvalsTypes = []

    #Save the data corresponding to pvalues (in order to print their lines in outputing VCF)
    pvalsCorrespondingChr = []
    pvalsCorrespondingStart = []
    pvalsCorrespondingRef = []
    pvalsCorrespondingAltNormal = []
    pvalsCorrespondingAltTumoral = []

    #Separate Germline in one side and somatic/N.variant!=T.variant in the other side in order to reduce the noise due to germline mutation
    germlineCorrespondingChr = []
    germlineCorrespondingStart = []
    germlineCorrespondingRef = []
    germlineCorrespondingAltNormal = []
    germlineCorrespondingAltTumoral = []

    for goodTumoral in goodVCFTumoral.collection.getAll():
        booleanHaveHomolog = False
        for goodNormal in goodVCFNormal.collection.getAll():
            if (goodNormal.chr == goodTumoral.chr and goodNormal.start == goodTumoral.start and goodNormal.reference != goodTumoral.reference):
                print("{} {} {} {} is a case of reference deletion in only one of the 2 datasets".format(goodNormal.chr,goodNormal.start,goodNormal.reference,goodTumoral.reference))
            if (goodNormal.chr == goodTumoral.chr and goodNormal.start == goodTumoral.start and goodNormal.reference == goodTumoral.reference):
                booleanHaveHomolog = True #Position is a variant in both Normal and tumoral datasets
                goodTumoralAlt = goodTumoral.alt.split(",")
                goodNormalAlt = goodNormal.alt.split(",")
                goodTumoralAltSet = set(goodTumoralAlt)
                goodNormalAltSet = set(goodNormalAlt)
                if(len(goodTumoralAltSet.union(goodNormalAltSet)) == len(goodTumoralAltSet.intersection(goodNormalAltSet))): #If the variants are the same in Normal and Tumoral datasets
                    germlineCorrespondingChr.append(goodTumoral.chr)
                    germlineCorrespondingStart.append(goodTumoral.start)
                    germlineCorrespondingRef.append(goodTumoral.reference)
                    germlineCorrespondingAltNormal.append(goodNormal.alt)
                    germlineCorrespondingAltTumoral.append(goodTumoral.alt)
                else: #At least one variant is different
                    pvalsCorrespondingChr.append(goodTumoral.chr)
                    pvalsCorrespondingStart.append(goodTumoral.start)
                    pvalsCorrespondingRef.append(goodTumoral.reference)
                    pvalsCorrespondingAltNormal.append(goodNormal.alt)
                    pvalsCorrespondingAltTumoral.append(goodTumoral.alt)
                    goodTumoralFreq = goodTumoral.gtFreq.split(",")
                    goodNormalFreq = goodNormal.gtFreq.split(",")
                    goodNormalCountRef = int(round(float(goodNormalFreq[0]) * int(goodNormal.totCount)) / 100)
                    goodNormalCountAlt = int(round(float(goodNormalFreq[1]) * int(goodNormal.totCount)) / 100)
                    goodTumoralCountRef = int(round(float(goodTumoralFreq[0]) * int(goodTumoral.totCount)) / 100)
                    goodTumoralCountAlt = int(round(float(goodTumoralFreq[1]) * int(goodTumoral.totCount)) / 100)
                    pvalsTypes.append("N.variant!=T.variant")


                    if (len(goodNormalAlt) == 1 and len(goodTumoralAlt) == 1): #If 1 variant in normal and 1 variant in tumoral, fisher test on contingency table 2*2
                        pval = fisher(goodNormalCountRef, goodTumoralCountRef, goodNormalCountAlt, goodTumoralCountAlt)
                        pvals.append(pval.two_tail)
                    else:
                        if (parameters['mutipleTestSum'] and len(goodNormalAlt) >= 1 and len(goodTumoralAlt) >= 1 and not(len(goodNormalAlt) == 1 and len(goodTumoralAlt) == 1)): #Perform the sum of alternative frequencies in order to be in the case of an 2*2 contingency table -> fisher test
                            sumNormalFreq = 0
                            sumTumoralFreq = 0
                            for i in range(2, len(goodNormalAlt)):
                                goodNormalCountAlt += int(round(float(goodNormalFreq[i]) * int(goodNormal.totCount)) / 100)
                                sumNormalFreq += goodNormalFreq[i]
                            for i in range(2, len(goodTumoralAlt)):
                                goodTumoralCountAlt += int(round(float(goodTumoralFreq[i]) * int(goodTumoral.totCount)) / 100)
                                sumTumoralFreq
                            pval = fisher(goodNormalCountRef, goodTumoralCountRef, goodNormalCountAlt, goodTumoralCountAlt)
                            pvals.append(pval.two_tail)
#                            print "{} {}".format(sumNormalFreq, sumTumoralFreq)
                        else: 
                            if (not(parameters['mutipleTestSum']) and len(goodNormalAlt) >= 1 and len(goodTumoralAlt) >= 1 and not(len(goodNormalAlt) == 1 and len(goodTumoralAlt) == 1)): # Compute fisher test for a table greater than 2*2 (need R)
                                #Compare variants in Normal and tumoral to deduce the size of the n * m contingency table
                                normalNucList = list(goodNormal.reference) + goodNormalAlt
                                tumoralNucList = list(goodTumoral.reference) + goodTumoralAlt
                                unionNucList = set(normalNucList) | set(tumoralNucList)
                                normalNucFreqList = []
                                tumoralNucFreqList = []
                                i = 0
                                for unionNuc in unionNucList:
                                    normalNucFreqList.append(0)
                                    tumoralNucFreqList.append(0)
                                    for j in range(0, len(goodNormalFreq)):
                                        if(unionNuc == normalNucList[j]):
                                            normalNucFreqList[i] = int(round(float(goodNormalFreq[j]) * int(goodNormal.totCount) / 100))

                                    for j in range(0, len(goodTumoralFreq)):
                                        if(unionNuc == tumoralNucList[j]):
                                            tumoralNucFreqList[i] = int(round(float(goodTumoralFreq[j]) * int(goodTumoral.totCount) / 100))
                                    i += 1

                                unionNucFreqList = []
                                for i in range(0, len(unionNucList)):
                                    unionNucFreqList.append(normalNucFreqList[i])
                                    unionNucFreqList.append(tumoralNucFreqList[i])
                                if(sum(normalNucFreqList+tumoralNucFreqList) == 0):
                                    #print("Pos {} Tref {} Talt {}".format(goodTumoral.start,goodTumoral.reference,goodTumoral.alt))
                                    pvals.append(1)
                                else:
                                    #Use R because fiher exact test does not exist in python for contingency tables greater than 2*2
                                    v = robjects.IntVector(normalNucFreqList+tumoralNucFreqList)
                                    m = robjects.r['matrix'](v,nrow=len(unionNucList))
                                    res = stats.fisher_test(m)
                                    pvals.append(res[0][0])

        if(not(booleanHaveHomolog)): #The tumoral variant has no homolog in the normal dataset -> which mean that the normal allele is the same as reference # or maybe in Trash
            pvalsCorrespondingChr.append(goodTumoral.chr)
            pvalsCorrespondingStart.append(goodTumoral.start)
            pvalsCorrespondingRef.append(goodTumoral.reference)
            pvalsCorrespondingAltTumoral.append(goodTumoral.alt)
            goodTumoralAlt = goodTumoral.alt.split(",")
            goodTumoralFreq = goodTumoral.gtFreq.split(",")
            goodTumoralCountRef = int(round(float(goodTumoralFreq[0]) * int(goodTumoral.totCount)) / 100)
            goodTumoralCountAlt = int(round(float(goodTumoralFreq[1]) * int(goodTumoral.totCount)) / 100)

            #Get back the normal coverage which is depthProfileNormal (not del)
            depthProfileNormalFields = str(depthProfileNormal.depthProfile[goodTumoral.chr][goodTumoral.start]).split('\t')
            goodNormalCountRef = int(depthProfileNormalFields[2])
        #    pvalsCorrespondingAltNormal.append(depthProfileNormalFields[3])

            if(len(goodTumoral.reference) > len(goodTumoral.alt)): # To avoid case with ref = GA and altT = G, if no mutation in normal then altN = G whereas it must be altN = GA
                pvalsCorrespondingAltNormal.append(goodTumoral.reference)
            else:
                pvalsCorrespondingAltNormal.append(depthProfileNormalFields[3])

            goodNormalCountAlt = 0
            goodTumoralCountRef = int(round(float(goodTumoralFreq[0]) * int(goodTumoral.totCount)) / 100)
            goodTumoralCountAlt = int(round(float(goodTumoralFreq[1]) * int(goodTumoral.totCount)) / 100)

            pvalsTypes.append("somatic")
            if(len(goodTumoralAlt) == 1): #If reference in normal and 1 variant in tumoral, fisher test on contingency table 2*2
                pval = fisher(goodNormalCountRef, goodTumoralCountRef, goodNormalCountAlt, goodTumoralCountAlt)
                pvals.append(pval.two_tail)
            else:
                if (parameters['mutipleTestSum'] and len(goodTumoralAlt) > 1): #Perform the sum of alternative frequencies in order to be in the case of an 2*2 contingency table -> fisher test
                    for i in range(2, len(goodTumoralAlt)):
                        goodTumoralCountAlt += int(round(float(goodTumoralFreq[i]) * int(goodTumoral.totCount)) / 100)
                    pval = fisher(goodNormalCountRef, goodTumoralCountRef, goodNormalCountAlt, goodTumoralCountAlt)
                    pvals.append(pval.two_tail)
                else: # Compute fisher test for a table greater than 2*2 (need R)
                    #Compare variants in Normal and tumoral to deduce the size of the n * m contingency table
                    if (not(parameters['mutipleTestSum']) and len(goodTumoralAlt) > 1):
                        normalNucList = list(goodTumoral.reference)
                        tumoralNucList = list(goodTumoral.reference) + goodTumoralAlt
                        unionNucList = set(normalNucList) | set(tumoralNucList)
                        normalNucFreqList = []
                        tumoralNucFreqList = []
                        i = 0
                        for unionNuc in unionNucList:
                            normalNucFreqList.append(0)
                            tumoralNucFreqList.append(0)

                        for j in range(0, len(goodTumoralFreq)):
                            if(unionNuc == tumoralNucList[j]):
                                tumoralNucFreqList[i] = int(round(float(goodTumoralFreq[j]) * int(goodTumoral.totCount) / 100))
                        i += 1
                        unionNucFreqList = []
                        for i in range(0, len(unionNucList)):
                            unionNucFreqList.append(normalNucFreqList[i])
                            unionNucFreqList.append(tumoralNucFreqList[i])
                        if(sum(normalNucFreqList+tumoralNucFreqList) == 0):
                            #pvalsTypes.pop()
                            #pvalsTypes.append("N.Filtered")
                        #    print("Pos {} Tref {} Talt {}".format(goodTumoral.start,goodTumoral.reference,goodTumoral.alt))
                            pvals.append(1)
                        else:
                            #Use R because fiher exact test does not exist in python for contingency tables greater than 2*2
                            v = robjects.IntVector(normalNucFreqList+tumoralNucFreqList)
                            m = robjects.r['matrix'](v,nrow=len(unionNucList))
                            res = stats.fisher_test(m)
                            pvals.append(res[0][0])

    #Multiple test correction according the given parameters
    if(parameters['methodFDR'] == "bonferroni"):
        fdrs = multipletests(pvals, alpha=parameters['alphaFDR'], method='bonferroni')
    if(parameters['methodFDR'] == "fdr_bh"):
        fdrs = multipletests(pvals, alpha=parameters['alphaFDR'], method='fdr_bh')
    if(parameters['methodFDR'] == "fdr_tsbh"):
        fdrs = multipletests(pvals, alpha=parameters['alphaFDR'], method='fdr_tsbh')

    #Output germline
    parameters['germline'] = parameters['prefix'] + ".NT.germline"
    io.register("germline", parameters['germline'], "w")
    germline = io.getIO('germline')
    print >>germline, "#status = 'germline' if the normal_alternative_bases set is the same as the tumoral_alternative_bases one || 'somatic' if reference_bases set the same as the normal_alternative_bases || 'N.variant!=T.variant' all other cases (at least one base differente between normal_alternative_bases and tumoral_alternative_bases one or between reference_bases and normal_alternative_bases"
    print >>germline, "#chromosomes\tlocations\treference_bases\tnormal_alternative_bases\ttumoral_alternative_bases\tstatus"
    for i in range(0, len(germlineCorrespondingChr)):
        print >>germline, "{}\t{}\t{}\t{}\t{}\tgermline".format(germlineCorrespondingChr[i],germlineCorrespondingStart[i],germlineCorrespondingRef[i],germlineCorrespondingAltNormal[i],germlineCorrespondingAltTumoral[i])

    #Output somatic
    parameters['somatic'] = parameters['prefix'] + ".NT.somatic"
    io.register("somatic", parameters['somatic'], "w")
    somatic = io.getIO('somatic')
    print >>somatic, "#status = 'germline' if the normal_alternative_bases set is the same as the tumoral_alternative_bases one || 'somatic' if reference_bases set the same as the normal_alternative_bases || 'N.variant!=T.variant' all other cases (at least one base differente between normal_alternative_bases and tumoral_alternative_bases one or between reference_bases and normal_alternative_bases"
    print >>somatic, "#test_results = 'True' if corrected_pvalue < alpha (so significative difference of allele frequencies between normal and tumoral data) || 'False' if corrected_pvalue > alpha (so non-significative difference of allele frequencies between normal and tumoral data)"
    print >>somatic, "#chromosomes\tlocations\treference_bases\tnormal_alternative_bases\ttumoral_alternative_bases\tstatus\tcorrected_pvalues\ttest_results"
    for i in range(0, len(fdrs[0])):
        if(fdrs[0][i]):
            print >>somatic, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(pvalsCorrespondingChr[i],pvalsCorrespondingStart[i],pvalsCorrespondingRef[i],pvalsCorrespondingAltNormal[i],pvalsCorrespondingAltTumoral[i],pvalsTypes[i],fdrs[1][i],fdrs[0][i])
        else:
            print >>somatic, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(pvalsCorrespondingChr[i],pvalsCorrespondingStart[i],pvalsCorrespondingRef[i],pvalsCorrespondingAltNormal[i],pvalsCorrespondingAltTumoral[i],pvalsTypes[i],fdrs[1][i],fdrs[0][i])

    io.endlog()
    sys.exit(io.unregister("all"))

####Varscan strategy (mostly the same except the Call reference and heterozygous parts) 

#If tumor match normal
##If tumor and normal do not match the reference
### Call germline
##Else
###Call reference

#else (tumor do not match normal)
##Calculate significance of allele frequeny difference by Fisher's Exact Test
###If difference is significant (pvalue<threshold)
####If normal matches reference
#####Call somatic
####If normal is heterozygous
#####Call LOH
####If normal and tumor are variant but different
######Call IndelFilter or Unknown
###Else (non significant)
####Call germline 

############
