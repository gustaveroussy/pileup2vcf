# encoding: utf-8

'''
Created on ???
Last Update: 24 novembre 2016

@author: Yannick Boursin
@contact: yannick.boursin@gustaveroussy.fr
@license: GNU GPLv3
@organization: Gustave Roussy
@version: 1.2
@todo:
'''

class VCF(object):
    '''VCF class
    This class is intended to take in charge VCF output format.'''
    def __init__(self, collection, parameters, indelsToAdd=None):
        self.collection = collection
        self.parameters = parameters
        self.indelsToAdd = indelsToAdd
        
    #Fonction permettant de switcher à un affichage VCF
    def printVCF(self, fh,  printHeader=False):
        '''Prints a VCF file using variant in collections. This distinguishes variants from two origins:
        - trash
        - good
        It also adds indels detected by IndelDetector at the end (again distinguishing between good and trash)'''
        if self.collection.isTrash:
            if (printHeader):
                print >>fh, '##fileformat=VCFv4.2'
                print >>fh, '##INFO=<ID=A,Number=1,Type=Integer,Description="Number of seen A nucleotides">'
                print >>fh, '##INFO=<ID=G,Number=1,Type=Integer,Description="Number of seen G nucleotides">'
                print >>fh, '##INFO=<ID=T,Number=1,Type=Integer,Description="Number of seen T nucleotides">'
                print >>fh, '##INFO=<ID=C,Number=1,Type=Integer,Description="Number of seen C nucleotides">'
                print >>fh, '##INFO=<ID=INS,Number=1,Type=Integer,Description="Number of seen inserted nucleotides">'
                print >>fh, '##INFO=<ID=DEL,Number=1,Type=Integer,Description="Number of seen deleted nucleotides">'
                print >>fh, '##INFO=<ID=LFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[alt+,alt-]]">'
                print >>fh, '##INFO=<ID=GFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[global alt+, global alt-]]">'
                print >>fh, '##INFO=<ID=GSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand. 50 means great, extremes means bad">'
                print >>fh, '##INFO=<ID=LSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand for the called allele. 50 means great, extremes means bad">'
                print >>fh, '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Total depth at this position">'
                print >>fh, '##INFO=<ID=CLG,Number=1,Type=Integer,Description="Clustering group number">'
                print >>fh, '##INFO=<ID=MiRQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a read quality upper than minReadQuality">'
                print >>fh, '##INFO=<ID=MiMQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a mapping quality upper than minMappingQuality">'
                print >>fh, '##INFO=<ID=MeRQ,Number=1,Type=Integer,Description="Supporting read quality mean (including read quality below minReadQuality)">'
                print >>fh, '##INFO=<ID=MeMQ,Number=1,Type=Integer,Description="Supporting mapping quality mean (including mapping quality below minMappingQuality)">'
                print >>fh, '##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Frequency of called variant">'
                print >>fh, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of called variant">'
                print >>fh, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype. 0/1 under 75% of allelic frequency otherwise 1/1">'
                print >>fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}".format(self.parameters['sName'])
            for el in self.collection.getAll():
                if (el.globalFilter == "True"):
                    try:
                        if (el.localFilter != True): filter = el.localFilter
                        else: print "WTF ?!?!"
                    except:
                        filter = []
                        for el2 in el.trashedAlleles:
                            filter.append(el2.localFilter)
                        filter = "; ".join(filter)
                else: filter = el.globalFilter
    
                try: el.calledLfs = '{0:.2f}'.format(float(el.calledLfs))
                except: el.calledLfs = 'NA'
                try: el.gfs = '{0:.2f}'.format(float(el.gfs))
                except: el.gfs = 'NA'
                try: el.sb = '{0:.2f}'.format(float(el.sb))
                except: el.sb = 'NA'
                try: el.calledLocalsb = '{0:.2f}'.format(float(el.calledLocalsb))
                except: el.calledLocalsb = 'NA'
                try: el.gtFreq = el.gtFreq
                except: el.gtFreq = 'NA'
                try: el.calledQuality = el.calledQuality
                except: el.calledQuality = 'NA'
                try: el.alt = el.alt
                except: el.alt    = 'NA'
                try: el.gtDepth = el.gtDepth
                except: el.gtDepth = 'NA'
                try: el.calledGenotype = str(el.calledGenotype)
                except: el.calledGenotype='NA'
                try: el.reference = str(el.reference)
                except: el.reference = "NA"

                #Round problem du to multiple rounds 
                comma = ",";
                gtFreqList = [float(i) for i in el.gtFreq.split(",")]
                if(sum(gtFreqList)!=100):
                    gtFreqList[0] += 100 - float(sum(gtFreqList))
                gtFreqList = [str(i) for i in gtFreqList]
                el.gtFreq = comma.join(gtFreqList)
                print >>fh, "{0}\t{1}\t.\t{2}\t{3}\t{4}\t{19}\tA={7};C={8};T={9};G={10};LFS={17};GFS={14};LSB={16};GSB={15};INS={12};DEL={13};TDP={11};CLG={20};MiRQ={21};;MiMQ={22};MeRQ={23};MeMQ={24}\tGT:FREQ:DP\t{18}:{5}:{6}".format(el.chr,
                                                                                            el.start,
                                                                                            el.reference,
                                                                                            el.alt,
                                                                                            el.calledQuality,
                                                                                            el.gtFreq,
                                                                                            el.gtDepth,
                                                                                            el.nucCount['A'],
                                                                                            el.nucCount['C'],
                                                                                            el.nucCount['T'],
                                                                                            el.nucCount['G'],
                                                                                            el.totCount,
                                                                                            el.nucCount['ins'],
                                                                                            el.nucCount['del'],
                                                                                            el.gfs,
                                                                                            el.sb,
                                                                                            el.calledLocalsb,
                                                                                            el.calledLfs,
                                                                                            el.calledGenotype,
                                                                                            filter,el.clusterGP,el.minReadQuality,el.minMappingQuality,el.meanReadQuality,el.meanMappingQuality)
            fh.flush()
            if self.indelsToAdd is not None:
                for el in self.indelsToAdd:
                    #print el.filterState
                    if el.filterState != "PASS":
                        print >>fh, str(el)
        else:
            if (printHeader):
                print >>fh, '##fileformat=VCFv4.2'
                print >>fh, '##INFO=<ID=A,Number=1,Type=Integer,Description="Number of seen A nucleotides">'
                print >>fh, '##INFO=<ID=G,Number=1,Type=Integer,Description="Number of seen G nucleotides">'
                print >>fh, '##INFO=<ID=T,Number=1,Type=Integer,Description="Number of seen T nucleotides">'
                print >>fh, '##INFO=<ID=C,Number=1,Type=Integer,Description="Number of seen C nucleotides">'
                print >>fh, '##INFO=<ID=INS,Number=1,Type=Integer,Description="Number of seen inserted nucleotides">'
                print >>fh, '##INFO=<ID=DEL,Number=1,Type=Integer,Description="Number of seen deleted nucleotides">'
                print >>fh, '##INFO=<ID=NUMMUTALL,Number=1,Type=Integer,Description="Number of seen allele at one position">'
                print >>fh, '##INFO=<ID=NUMMUTFILT,Number=1,Type=Integer,Description="Number of seen considered good allele at one position">'
                print >>fh, '##INFO=<ID=LFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[alt+,alt-]]">'
                print >>fh, '##INFO=<ID=GFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[global alt+, global alt-]]">'
                print >>fh, '##INFO=<ID=GSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand. 50 means great, extremes means bad">'
                print >>fh, '##INFO=<ID=LSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand for the called allele. 50 means great, extremes means bad">'
                print >>fh, '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Total depth at this position">'
                print >>fh, '##INFO=<ID=CLG,Number=1,Type=Integer,Description="Clustering group number">'
                print >>fh, '##INFO=<ID=MiRQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a read quality upper than minReadQuality">'
                print >>fh, '##INFO=<ID=MiMQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a mapping quality upper than minMappingQuality">'
                print >>fh, '##INFO=<ID=MeRQ,Number=1,Type=Integer,Description="Supporting read quality mean (including read quality below minReadQuality)">'
                print >>fh, '##INFO=<ID=MeMQ,Number=1,Type=Integer,Description="Supporting mapping quality mean (including mapping quality below minMappingQuality)">'
                print >>fh, '##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Frequency of called variant">'
                print >>fh, '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of called alleles">'
                print >>fh, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth at this Position">'
                print >>fh, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype. 0/1 under 75% of allelic frequency otherwise 1/1">'
                print >>fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}".format(self.parameters['sName'])
            for el in self.collection.getAll():
                print >>fh, "{0}\t{1}\t.\t{2}\t{3}\t{4:.2f}\tPASS\tNUMMUTALL={19};NUMMUTFILT={20};A={7};C={8};T={9};G={10};LFS={17:.2f};GFS={14:.2f};LSB={16:.2f};GSB={15:.2f};INS={12};DEL={13};TDP={11};CLG={22};MiRQ={23};MiMQ={24};MeRQ={25};MeMQ={26}\tGT:FREQ:AD:DP\t{18}:{5}:{6}:{21}".format(el.chr,
                                                                                            el.start,
                                                                                            el.reference,
                                                                                            el.alt,
                                                                                            float(el.calledQuality),
                                                                                            el.gtFreq,
                                                                                            el.gtDepth,
                                                                                            el.nucCount['A'],
                                                                                            el.nucCount['C'],
                                                                                            el.nucCount['T'],
                                                                                            el.nucCount['G'],
                                                                                            el.totCount,
                                                                                            el.nucCount['ins'],
                                                                                            el.nucCount['del'],
                                                                                            el.gfs,
                                                                                            float(el.sb),
                                                                                            float(el.calledLocalsb),
                                                                                            float(el.calledLfs),
                                                                                            el.calledGenotype,
                                                                                            len(el.stackedAlleles),
                                                                                            len(el.goodAlleles),
                                                                                            el.totCount,el.clusterGP,el.minReadQuality,el.minMappingQuality,el.meanReadQuality,el.meanMappingQuality)
            if self.indelsToAdd is not None:
                for el in self.indelsToAdd:
                    if el.filterState == "PASS":
                        #print str(el)
                        print >>fh, str(el)
            fh.flush()


    #Fonction de calculer le nombre de SNP inclus dans une fenetre de taille k (paramétrable)
    def printWindowK(self, fh_cov, fh_group, k, printHeader=False):
        '''Prints 2 WindowK file using variant in collections. Use only with variants from good.
        * Set: cur.globalFilter, cur.clusterGP
        * Write: *.cov.bed (bed of segments(=replace variant position by segment of size k) according the coverage)
        * Write: *.group (equivalence between clusterID and variant positions included in this cluster)
        '''
        if (printHeader):
            print >>fh_cov, "#Chromosome\tStartWindow\tEndWindow\tClusterID\tCoverage\tLocationVariantID"
            print >>fh_group, "#ClusterID\tLocationVariantID"
        lastChr = ""
        lastStarts = []
        lastStops = []
        listID = []
        clusterElmt = []
        clusterIDList = []
        clusterNum = 1
        cov = 1
        for el in self.collection.getAll(): #Pour tous les éléments de la collection 
            if(len(lastStops)+len(lastStarts) == 0): # Si c'est le premier élément à ajouter
                lastChr = el.chr
                lastStarts.append(int(el.start-k/2))
                lastStops.append(int(el.start+k/2))
                listID.append(str(el.start))
                clusterElmt.append(el)
                clusterIDList.append(str(el.start))
                el.clusterGP = clusterNum
            else: #Si on a déjà stocké des éléments
                if(lastChr == el.chr and lastStops[len(lastStops)-1] > el.start-k/2): #Si on est sur le même chomosome et qu'il y a chevauchement
                    lastStarts.append(int(el.start-k/2))
                    lastStops.append(int(el.start+k/2))
                    listID.append(str(el.start))
                    clusterElmt.append(el)
                    clusterIDList.append(str(el.start))
                    el.clusterGP = clusterNum
                else: #Sinon on imprime
                    while len(lastStops)+len(lastStarts) > 0: #Tant qu'on a pas tout imprimé
                        if(len(lastStarts)>0):
                            if(lastStops[0]>lastStarts[0]):
                                if(len(lastStarts)>1): #On imprime les segments correspondant aux débuts et on conserve le dernier
                                    print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStarts[0],lastStarts[0+1]-1,clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                                    lastStarts.remove(lastStarts[0])
                                    cov += 1
                                if(len(lastStarts)==1): #On imprime le segment correspondant au dernier début et à la première fin, et on enlève le dernier début
                                    print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStarts[0],lastStops[0],clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                                    lastStarts.remove(lastStarts[0])
                                    listID.remove(listID[0])
                                    cov -= 1
                            else:
                                if(len(lastStops)>1): #On imprime les segments correspondant aux fins
                                    print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStops[0]+1,lastStops[0+1],clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                                    lastStops.remove(lastStops[0])
                                    listID.remove(listID[0])
                                    cov -= 1
                                if(len(lastStops)==1): #On enlève la dernière fin pour sortir de la boucle while
                                    lastStops.remove(lastStops[0])
                        else:
                            if(len(lastStops)>1): #On imprime les segments correspondant aux fins
                                print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStops[0]+1,lastStops[0+1],clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                                lastStops.remove(lastStops[0])
                                listID.remove(listID[0])
                                cov -= 1
                            if(len(lastStops)==1): #On enlève la dernière fin pour sortir de la boucle while
                                lastStops.remove(lastStops[0])

                    print >>fh_group, "{}\t{}".format(clusterNum,",".join(clusterIDList)) #Affiche l'équivalence des groupes
                    lastChr = el.chr
                    lastStarts.append(int(el.start-k/2))
                    lastStops.append(int(el.start+k/2))
                    listID.append(str(el.start))
                    cov = 1
                    clusterNum += 1
                    clusterElmt = []
                    clusterElmt.append(str(el.start))
                    clusterIDList = []
                    clusterIDList.append(str(el.start))
                    el.clusterGP = clusterNum

#Pour le dernier cluster
        while len(lastStops)+len(lastStarts) > 0: #Tant qu'on a pas tout imprimé
            if(len(lastStarts)>0):
                if(lastStops[0]>lastStarts[0]):
                    if(len(lastStarts)>1): #On imprime les segments correspondant aux débuts et on conserve le dernier
                        print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStarts[0],lastStarts[0+1]-1,clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                        lastStarts.remove(lastStarts[0])
                        cov += 1
                    if(len(lastStarts)==1): #On imprime le segment correspondant au dernier début et à la première fin, et on enlève le dernier début
                        print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStarts[0],lastStops[0],clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                        lastStarts.remove(lastStarts[0])
                        listID.remove(listID[0])
                        cov -= 1
                else:
                    if(len(lastStops)>1): #On imprime les segments correspondant aux fins
                        print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStops[0]+1,lastStops[0+1],clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                        lastStops.remove(lastStops[0])
                        listID.remove(listID[0])
                        cov -= 1
                    if(len(lastStops)==1): #On enlève la dernière fin pour sortir de la boucle while
                        lastStops.remove(lastStops[0])
            else:
                if(len(lastStops)>1): #On imprime les segments correspondant aux fins
                    print >>fh_cov, "{}\t{}\t{}\t{}\t{}\t{}".format(lastChr,lastStops[0]+1,lastStops[0+1],clusterNum,cov,",".join(listID[:cov])) #Afficher le dernier segment sauvegardé
                    lastStops.remove(lastStops[0])
                    listID.remove(listID[0])
                    cov -= 1
                if(len(lastStops)==1): #On enlève la dernière fin pour sortir de la boucle while
                    lastStops.remove(lastStops[0])

        print >>fh_group, "{}\t{}".format(clusterNum,",".join(clusterIDList))#Affiche l'équivalence des groupes

    #Fonction de filtrer les GoodVariants qui tombent dans des SSR
    def printVCFUsingSSR(self,SSRlist,GoodVariantsOutOfSSR,TrashVariantsInSSR,printHeader):
        '''Prints a VCF file using variant in collections. This distinguishes variants from two origins:
        - trash
        It also adds indels detected by IndelDetector at the end (again distinguishing between good and trash)'''
        f = open(SSRlist, "r")
        line = f.readline().rstrip("\n")
        elem = line.split('\t')
        bool_same_chromosome = 0
        bool_EOF = 0
        chr = elem[0]
        beg = int(elem[1])
        end = int(elem[2])
        lastchr = elem[0]
        lastbeg = int(elem[1])
        lastend = int(elem[2])
        if not self.collection.isTrash:
            if (printHeader):
                print >>GoodVariantsOutOfSSR, '##fileformat=VCFv4.2'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=A,Number=1,Type=Integer,Description="Number of seen A nucleotides">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=G,Number=1,Type=Integer,Description="Number of seen G nucleotides">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=T,Number=1,Type=Integer,Description="Number of seen T nucleotides">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=C,Number=1,Type=Integer,Description="Number of seen C nucleotides">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=INS,Number=1,Type=Integer,Description="Number of seen inserted nucleotides">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=DEL,Number=1,Type=Integer,Description="Number of seen deleted nucleotides">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=NUMMUTALL,Number=1,Type=Integer,Description="Number of seen allele at one position">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=NUMMUTFILT,Number=1,Type=Integer,Description="Number of seen considered good allele at one position">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=LFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[alt+,alt-]]">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=GFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[global alt+, global alt-]]">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=GSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand. 50 means great, extremes means bad">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=LSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand for the called allele. 50 means great, extremes means bad">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Total depth at this position">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=CLG,Number=1,Type=Integer,Description="Clustering group number">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=MiRQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a read quality upper than minReadQuality">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=MiMQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a mapping quality upper than minMappingQuality">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=MeRQ,Number=1,Type=Integer,Description="Supporting read quality mean (including read quality below minReadQuality)">'
                print >>GoodVariantsOutOfSSR, '##INFO=<ID=MeMQ,Number=1,Type=Integer,Description="Supporting mapping quality mean (including mapping quality below minMappingQuality)">'
                print >>GoodVariantsOutOfSSR, '##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Frequency of called variant">'
                print >>GoodVariantsOutOfSSR, '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of called alleles">'
                print >>GoodVariantsOutOfSSR, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth at this Position">'
                print >>GoodVariantsOutOfSSR, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype. 0/1 under 75% of allelic frequency otherwise 1/1">'
                print >>GoodVariantsOutOfSSR, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}".format(self.parameters['sName'])
                
                print >>TrashVariantsInSSR, '##fileformat=VCFv4.2'
                print >>TrashVariantsInSSR, '##INFO=<ID=A,Number=1,Type=Integer,Description="Number of seen A nucleotides">'
                print >>TrashVariantsInSSR, '##INFO=<ID=G,Number=1,Type=Integer,Description="Number of seen G nucleotides">'
                print >>TrashVariantsInSSR, '##INFO=<ID=T,Number=1,Type=Integer,Description="Number of seen T nucleotides">'
                print >>TrashVariantsInSSR, '##INFO=<ID=C,Number=1,Type=Integer,Description="Number of seen C nucleotides">'
                print >>TrashVariantsInSSR, '##INFO=<ID=INS,Number=1,Type=Integer,Description="Number of seen inserted nucleotides">'
                print >>TrashVariantsInSSR, '##INFO=<ID=DEL,Number=1,Type=Integer,Description="Number of seen deleted nucleotides">'
                print >>TrashVariantsInSSR, '##INFO=<ID=NUMMUTALL,Number=1,Type=Integer,Description="Number of seen allele at one position">'
                print >>TrashVariantsInSSR, '##INFO=<ID=NUMMUTFILT,Number=1,Type=Integer,Description="Number of seen considered good allele at one position">'
                print >>TrashVariantsInSSR, '##INFO=<ID=LFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[alt+,alt-]]">'
                print >>TrashVariantsInSSR, '##INFO=<ID=GFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[global alt+, global alt-]]">'
                print >>TrashVariantsInSSR, '##INFO=<ID=GSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand. 50 means great, extremes means bad">'
                print >>TrashVariantsInSSR, '##INFO=<ID=LSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand for the called allele. 50 means great, extremes means bad">'
                print >>TrashVariantsInSSR, '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Total depth at this position">'
                print >>TrashVariantsInSSR, '##INFO=<ID=CLG,Number=1,T13ype=Integer,Description="Clustering group number">'
                print >>TrashVariantsInSSR, '##INFO=<ID=MiRQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a read quality upper than minReadQuality">'
                print >>TrashVariantsInSSR, '##INFO=<ID=MiMQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a mapping quality upper than minMappingQuality">'
                print >>TrashVariantsInSSR, '##INFO=<ID=MeRQ,Number=1,Type=Integer,Description="Supporting read quality mean (including read quality below minReadQuality)">'
                print >>TrashVariantsInSSR, '##INFO=<ID=MeMQ,Number=1,Type=Integer,Description="Supporting mapping quality mean (including mapping quality below minMappingQuality)">'
                print >>TrashVariantsInSSR, '##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Frequency of called variant">'
                print >>TrashVariantsInSSR, '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of called alleles">'
                print >>TrashVariantsInSSR, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth at this Position">'
                print >>TrashVariantsInSSR, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype. 0/1 under 75% of allelic frequency otherwise 1/1">'
                print >>TrashVariantsInSSR, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}".format(self.parameters['sName'])
            for el in self.collection.getAll():
                if(elem[0] != el.chr):
                    bool_same_chromosome = 0
                else:
                    bool_same_chromosome = 1
                print("### {} {}".format(elem[0], el.chr))
                if(elem[0] != el.chr): #different chromosome than last SSR
                    print("different chromosome")
                    while(elem[0] != el.chr and bool_EOF==0 and end < el.start): #while we are on different chromosome
                        #print("while")
                        line = f.readline().rstrip("\n")
                        lastbeg = elem[1]
                        lastend = elem[2]
                        lastchr = elem[0]
                        if(line==""):
                            bool_EOF = 1
                        elem = line.split('\t')
                        chr = elem[0]
                        beg = int(elem[1])
                        end = int(elem[2])
                        if(elem[0] != el.chr):
                            bool_same_chromosome = 0
                            #print("Bool same chr")
                        else:
                            bool_same_chromosome = 1
                            print("different chr")
                    print("Passing 1 {} {} {}".format(elem[0], elem[1], elem[2]))
                    #break
                    
                if(lastchr == el.chr): #same chromosome than last SSR
                    #print("same chromosome")
                    #print("current {} {} {}".format(elem[0], elem[1], elem[2]))
                    while(end < el.start): #while the SSR end is lower to the variants begining
                        #print("last {} {} {}".format(elem[0], elem[1], elem[2]))
                        #print("elmt {} {}".format(el.chr, el.start))
                        lastchr = elem[0]
                        lastbeg = elem[1]
                        lastend = elem[2]
                        #print("while location")
                        line = f.readline().rstrip("\n")
                        elem = line.split('\t')
                        chr = elem[0]
                        beg = int(elem[1])
                        end = int(elem[2])
                        if(chr != el.chr):
                            break

                    if(lastchr == el.chr and el.start > lastbeg and el.start < lastend): #variant in a SSR
                        print("print OUT")
                        print >>TrashVariantsInSSR, "{0}\t{1}\t.\t{2}\t{3}\t{4:.2f}\tPASS\tNUMMUTALL={19};NUMMUTFILT={20};A={7};C={8};T={9};G={10};LFS={17:.2f};GFS={14:.2f};LSB={16:.2f};GSB={15:.2f};INS={12};DEL={13};TDP={11};CLG={22};MiRQ={23};MiMQ={24};MeRQ={25};MeMQ={26}\tGT:FREQ:AD:DP\t{18}:{5}:{6}:{21}".format(el.chr,
                                                                                            el.start,
                                                                                            el.reference,
                                                                                            el.alt,
                                                                                            float(el.calledQuality),
                                                                                            el.gtFreq,
                                                                                            el.gtDepth,
                                                                                            el.nucCount['A'],
                                                                                            el.nucCount['C'],
                                                                                            el.nucCount['T'],
                                                                                            el.nucCount['G'],
                                                                                            el.totCount,
                                                                                            el.nucCount['ins'],
                                                                                            el.nucCount['del'],
                                                                                            el.gfs,
                                                                                            float(el.sb),
                                                                                            float(el.calledLocalsb),
                                                                                            float(el.calledLfs),
                                                                                            el.calledGenotype,
                                                                                            len(el.stackedAlleles),
                                                                                            len(el.goodAlleles),
                                                                                            el.totCount,el.clusterGP,el.minReadQuality,el.minMappingQuality,el.meanReadQuality,el.meanMappingQuality)
                    else:
                        print("print IN")
                        print >>GoodVariantsOutOfSSR, "{0}\t{1}\t.\t{2}\t{3}\t{4:.2f}\tPASS\tNUMMUTALL={19};NUMMUTFILT={20};A={7};C={8};T={9};G={10};LFS={17:.2f};GFS={14:.2f};LSB={16:.2f};GSB={15:.2f};INS={12};DEL={13};TDP={11};CLG={22};MiRQ={23};MiMQ={24};MeRQ={25};MeMQ={26}\tGT:FREQ:AD:DP\t{18}:{5}:{6}:{21}".format(el.chr,
                                                                                            el.start,
                                                                                            el.reference,
                                                                                            el.alt,
                                                                                            float(el.calledQuality),
                                                                                            el.gtFreq,
                                                                                            el.gtDepth,
                                                                                            el.nucCount['A'],
                                                                                            el.nucCount['C'],
                                                                                            el.nucCount['T'],
                                                                                            el.nucCount['G'],
                                                                                            el.totCount,
                                                                                            el.nucCount['ins'],
                                                                                            el.nucCount['del'],
                                                                                            el.gfs,
                                                                                            float(el.sb),
                                                                                            float(el.calledLocalsb),
                                                                                            float(el.calledLfs),
                                                                                            el.calledGenotype,
                                                                                            len(el.stackedAlleles),
                                                                                            len(el.goodAlleles),
                                                                                            el.totCount,el.clusterGP,el.minReadQuality,el.minMappingQuality,el.meanReadQuality,el.meanMappingQuality)

            GoodVariantsOutOfSSR.flush()
            TrashVariantsInSSR.flush()
            print("END")
    
    
def printTrashWhileReading(el, fh, printHeader):
    '''Prints a variant in a VCF file. The variants given to this function must be mpileup line which not pass the filters'''
    if (printHeader):
        print >>fh, '##fileformat=VCFv4.2'
        print >>fh, '##INFO=<ID=A,Number=1,Type=Integer,Description="Number of seen A nucleotides">'
        print >>fh, '##INFO=<ID=G,Number=1,Type=Integer,Description="Number of seen G nucleotides">'
        print >>fh, '##INFO=<ID=T,Number=1,Type=Integer,Description="Number of seen T nucleotides">'
        print >>fh, '##INFO=<ID=C,Number=1,Type=Integer,Description="Number of seen C nucleotides">'
        print >>fh, '##INFO=<ID=INS,Number=1,Type=Integer,Description="Number of seen inserted nucleotides">'
        print >>fh, '##INFO=<ID=DEL,Number=1,Type=Integer,Description="Number of seen deleted nucleotides">'
        print >>fh, '##INFO=<ID=LFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[alt+,alt-]]">'
        print >>fh, '##INFO=<ID=GFS,Number=1,Type=Integer,Description="Phred-scaled p-value of exact fisher test on the following contingency table: [[ref+, ref-],[global alt+, global alt-]]">'
        print >>fh, '##INFO=<ID=GSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand. 50 means great, extremes means bad">'
        print >>fh, '##INFO=<ID=LSB,Number=1,Type=Integer,Description="Percentaged ratio of reads on positive and negative strand for the called allele. 50 means great, extremes means bad">'
        print >>fh, '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Total depth at this position">'
        print >>fh, '##INFO=<ID=CLG,Number=1,Type=Integer,Description="Clustering group number">'
        print >>fh, '##INFO=<ID=MiRQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a read quality upper than minReadQuality">'
        print >>fh, '##INFO=<ID=MiMQ,Number=1,Type=Integer,Description="Number of read supporting the variant with a mapping quality upper than minMappingQuality">'
        print >>fh, '##INFO=<ID=MeRQ,Number=1,Type=Integer,Description="Supporting read quality mean (including read quality below minReadQuality)">'
        print >>fh, '##INFO=<ID=MeMQ,Number=1,Type=Integer,Description="Supporting mapping quality mean (including mapping quality below minMappingQuality)">'
        print >>fh, '##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Frequency of called variant">'
        print >>fh, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of called variant">'
        print >>fh, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype. 0/1 under 75% of allelic frequency otherwise 1/1">'
        #print >>fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}".format(self.parameters['sName'])
        print >>fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNNNN"
    else:
        if (el.globalFilter == True):
            try:
                if (el.localFilter != True): filter = el.localFilter
                else: print "WTF ?!?!"
            except:
                filter = []
                for el2 in el.trashedAlleles:
                    filter.append(el2.localFilter)
                filter = "; ".join(filter)
        else: filter = el.globalFilter
    
        try: el.calledLfs = '{0:.2f}'.format(float(el.calledLfs))
        except: el.calledLfs = 'NA'
        try: el.gfs = '{0:.2f}'.format(float(el.gfs))
        except: el.gfs = 'NA'
        try: el.sb = '{0:.2f}'.format(float(el.sb))
        except: el.sb = 'NA'
        try: el.calledLocalsb = '{0:.2f}'.format(float(el.calledLocalsb))
        except: el.calledLocalsb = 'NA'
        try: el.gtFreq = el.gtFreq
        except: el.gtFreq = 'NA'
        try: el.calledQuality = el.calledQuality
        except: el.calledQuality = 'NA'
        try: el.alt = el.alt
        except: el.alt    = 'NA'
        try: el.gtDepth = el.gtDepth
        except: el.gtDepth = 'NA'
        try: el.calledGenotype = str(el.calledGenotype)
        except: el.calledGenotype='NA'
        try: el.reference = str(el.reference)
        except: el.reference = "NA"

        #Round problem du to multiple rounds 
#        comma = ",";
#        gtFreqList = [float(i) for i in el.gtFreq.split(",")]
#        if(sum(gtFreqList)!=100):
#            gtFreqList[0] += 100 - float(sum(gtFreqList))
#        gtFreqList = [str(i) for i in gtFreqList]
#        el.gtFreq = comma.join(gtFreqList)
        print >>fh, "{0}\t{1}\t.\t{2}\t{3}\t{4}\t{19}\tA={7};C={8};T={9};G={10};LFS={17};GFS={14};LSB={16};GSB={15};INS={12};DEL={13};TDP={11};CLG={20};MiRQ={21};MiMQ={22};MeRQ={23};MeMQ={24}\tGT:FREQ:DP\t{18}:{5}:{6}".format(el.chr,
                                                                                            el.start,
                                                                                            el.reference,
                                                                                            el.alt,
                                                                                            el.calledQuality,
                                                                                            el.gtFreq,
                                                                                            el.gtDepth,
                                                                                            el.nucCount['A'],
                                                                                            el.nucCount['C'],
                                                                                            el.nucCount['T'],
                                                                                            el.nucCount['G'],
                                                                                            el.totCount,
                                                                                            el.nucCount['ins'],
                                                                                            el.nucCount['del'],
                                                                                            el.gfs,
                                                                                            el.sb,
                                                                                            el.calledLocalsb,
                                                                                            el.calledLfs,
                                                                                            el.calledGenotype,
                                                                                            filter,el.clusterGP,el.minReadQuality,el.minMappingQuality,el.meanReadQuality,el.meanMappingQuality)
        fh.flush()

