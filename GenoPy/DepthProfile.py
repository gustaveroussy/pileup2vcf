#!/usr/bin/env python
# encoding: utf-8

'''
Created on 3 août 2015
Last Update: 24 novembre 2016

@author: Yannick Boursin
@contact: yannick.boursin@gustaveroussy.fr
@license: GNU GPLv3
@organization: Gustave Roussy
@version: 1.4
@todo: Inspect keepInMemory and binDynamically parameters
'''

# Default library
from collections import OrderedDict
from GenoPy.Genomic import *

###############################################################################
########################    DepthProfile Class    #############################
###############################################################################

class DepthProfile(object):
    '''This class will handle Position objects which are is a very thorough way to represent data
    exposed in MPileup files.
    
    Position objects are pileup-like. They contain the following informations:
    - Chromosome (str). Ex: "chr1", "1"
    - Start (int). Ex: 8721852
    - Depth (int). How many reads cover this position ? Ex: 50
    - refNuc (str). Ex: A
    - LastPosition (int). Equals: Start - 1
    - inBed (bool). If position is inside provided BED-file: True, else False. Ex: False
    - quality (list). Ex: [12, 11, 24, 40, ..., 32]
    - mapQual (list). Ex: [12, 11, 24, 40, ..., 32]
    - avgQuality (float). Equals: scipy.mean(quality)
    - avgMapQuality (float). Equals: scipy.mean(mapQual)
    - nbAltPassFilter (int). How many validated alternative alleles for this position ?. Ex : 1
    - nbAltFiltered (int). How many filtered alternative alleles for this position ?. Ex: 1
    - nbTotalAlt (int). How many alternative alleles for this position ? Equals: nbAltPassFilter + nbAltFiltered
    - id (str). Equals: see getId()
    
    DepthProfile class is intended to centralize them and their methods:
    * writeDepthProfile => Supposed to help processing those Position objects by dumping them
    * add => add a Position object to DepthProfile
    * walk => allows for calling writeDepthProfile and switching orderedDicts if chromosome change'''
    def __init__(self, handler, io):
        self.dph = handler
        self.depthProfile = OrderedDict() # Key: chr. Inside: key=Position. Value=Depth
        self.wroteHeader = False
        self.io = io

    def writeDepthProfile(self, printDepthProfile, keepInMemory, binDynamically):
        '''Method that print and/or dumps depthProfile data
        
        This method is intended to print and/or dump DepthProfile data. It can do so in different ways:
            - if binDynamically is True: it will call the binDynamically.bin bound method
            - if binDynamically is False: it will print Position objects to a filehandler (self.dph)
        # TODO: inspect arguments of this function (keepInMemory is weird)
        
        * Input: <bool> printDepthProfile, <bool> keepInMemory, <bool> binDynamically'''
        fh = self.dph
        header = self.wroteHeader
        if not header and printDepthProfile:
            print >>fh, "##Pileup2VCF - DepthProfile - v1.1"
            print >>fh, "#Chromosome\tStart\tDepth\trefNuc\tlast_position\tinBed\tquality\tmapQual\tavgQuality\tavgMapQuality\tnbAltPassFilter\tnbAltFiltered\tnbTotalAlt\tid"
            self.wroteHeader = True
        if binDynamically is not False:
            print 'Binning'
            binDynamically.bin(self.depthProfile, binDynamically=True)
        if printDepthProfile:
            for k, v in self.depthProfile.items():
                for v2 in v.values():
                    print >>fh, str(v2)
            fh.flush()
            self.io.log("Reseting DepthProfile", loglevel=3)
            if not keepInMemory:
                del self.depthProfile
                self.depthProfile = OrderedDict()
        

    def add(self, chromosome, position, PositionObject):
        '''This method adds a Position object to DepthProfile
        
        * Input: <str> chromosome, <int> position, <Position-like> PositionObject'''
        self.depthProfile[chromosome][position] = PositionObject

    def walk(self, Chromosome=None, printDepthProfile=False, keepInMemory=False, binDynamically=False):
        '''This method allows for calling data dump method and creating new OrderedDict upon chromosome change
        
        * Input: <str> Chromsome, <bool> keepInMemory, <bool> binDynamically'''
        self.writeDepthProfile(printDepthProfile=printDepthProfile, keepInMemory=keepInMemory, binDynamically=binDynamically)
        if Chromosome is not None:
            self.depthProfile[Chromosome] = OrderedDict()

    def addVariants(self, uncl):
        dphW = self.io.getIO('depthProfileAddMutW')
        dphR = self.io.getIO('depthProfileAddMutR')
        print >>dphW, "##Pileup2VCF - DepthProfile - v1.1"
        print >>dphW, "#Chromosome\tStart\tDepth\trefNuc\tlast_position\tinBed\tquality\tmapQual\tavgQuality\tavgMapQuality\tnbAltPassFilter\tnbAltFiltered\tnbTotalAlt\tid"
            
        for pos in tempToObjects(dphR):
            posHash = '{}-{}'.format(pos.chr, pos.start)
            al = uncl.varDic[posHash]
            if al is not False:
                pos.nbAltPassFilter = al[0]
                pos.nbAltFiltered = al[1]
                pos.nbTotalAlt = pos.nbAltPassFilter + pos.nbAltFiltered
            else: pass
            print >>dphW, str(pos)
            
