# encoding: utf-8

from collections import Counter, defaultdict

### Got it from: http://stackoverflow.com/questions/13414366/python-multilevel-nested-lists
def flatten_all(iterable):
    for elem in iterable:
        if not isinstance(elem, list):
            yield elem
        else:
            for x in flatten_all(elem):
                yield x

class Statistics(object):
    '''Class responsible for various statistics collection
    
    This class takes as input a GenotypedVariantCollection (still to be created) 
    and is initialized after variant caller as occured. Any variant collection
    should be OK as long as it contains Variant objects or Variant-like objects 
    (displaying the right methods and attributes), such as titv() or statsAllAllele()
    Still working in it (this shouldn't work as for now)'''
    def __init__(self):
        self.dico_stats = {}
        self.dico_stats['all'] = defaultdict(int)
    
    def init_positions(self, parameters):
        self.dico_stats['positions'] = parameters['lineNb']
        
    def init_variants(self, VariantCollection, type=None):
        '''This function will allow for delaying initialization of Statistics object
        
        Computation of statistics may be delayed. We might add a method "add()" that would compute on-the-fly
        statistics for each variant later, but as for now, this method must be called after variant calling.
        It computes titv ratios, counts the called variants, counts total variants ...'''
        if type is not None:
            self.dico_stats[type] = {}
            titv = [variant.titv() for variant in VariantCollection.vl]
            allTiTv = list(flatten_all([variant.statsAllAllele(type=type) for variant in VariantCollection.vl]))
            self.dico_stats[type]['majorAlleleTiTv'] = Counter(titv)
            self.dico_stats[type]["allAlleleTiTv"] = Counter(allTiTv)
            self.dico_stats[type]["called"] = len(VariantCollection.getAll())
            t = float(sum(self.dico_stats[type]["allAlleleTiTv"].values()))
            #print self.dico_stats[type]["allAlleleTiTv"]
            self.dico_stats[type]['calledAlleleNumber'] = t if t != 0 else 1
            # Maintain statistics for all collections added
            if not isinstance(self.dico_stats['all']['majorAlleleTiTv'], Counter):
                self.dico_stats['all']['majorAlleleTiTv'] = Counter(titv)
                self.dico_stats['all']["allAlleleTiTv"] = Counter(allTiTv)
            else:
                self.dico_stats['all']['majorAlleleTiTv'] = self.dico_stats['all']['majorAlleleTiTv'] + Counter(titv)
                self.dico_stats['all']["allAlleleTiTv"] = self.dico_stats['all']["allAlleleTiTv"] + Counter(allTiTv)
            self.dico_stats['all']['called'] += self.dico_stats[type]["called"]
            self.dico_stats['all']['calledAlleleNumber'] += self.dico_stats[type]['calledAlleleNumber']
        
#     def add(self, var):
#         '''WriteMe'''
#         pass
        
    def printStatistics(self, type, stath):
        '''Prints statistics'''
        # Attention à bien changer la référence, ici self désigne l'object pileup, il faut donc qu'un object statistics
        # soit créé dans l'object pileup, et que celui-ci soit utilisé.

        # Nous ferons attention à donner l'occasion à l'objet statistique de compter tout ce qui est utilisable en matière de rapport
        # sur le variant calling.

        # Nous lui offrirons également une méthode pour calculer de lui même certaines choses. (genre ti/tv)
        globalStat = True
        if type in self.dico_stats.keys():
            called = self.dico_stats[type]["called"]
            if called == 0: called = 1 # POTENTIAL BUG ?
            calledAlleleNumber = self.dico_stats[type]['calledAlleleNumber']
            total = self.dico_stats['positions']
            if (self.dico_stats[type]['majorAlleleTiTv']['tv'] == 0):
                dividedTiTv = "** Ti/Tv:\t+INF\t+INF%"
            else:
                dividedTiTv = "** Ti/Tv:\t{0:.2f}\t{1:.2f}%".format( float( float(self.dico_stats[type]['majorAlleleTiTv']['ti']) / float(self.dico_stats[type]['majorAlleleTiTv']['tv']) ), float( float(self.dico_stats[type]['majorAlleleTiTv']['ti']) * 100 / float(self.dico_stats[type]['majorAlleleTiTv']['tv']) ) )
            if (self.dico_stats[type]["allAlleleTiTv"]['tv'] == 0):
                dividedAllTiTv = "** Ti/Tv:\t+INF\t+INF%"
            else:
                dividedAllTiTv = "** Ti/Tv:\t{0:.2f}\t{1:.2f}%".format( float( float(self.dico_stats[type]["allAlleleTiTv"]['ti']) / float(self.dico_stats[type]["allAlleleTiTv"]['tv']) ), float( float(self.dico_stats[type]["allAlleleTiTv"]['ti']) * 100 / float(self.dico_stats[type]["allAlleleTiTv"]['tv']) ) )
            if globalStat:
                print >>stath, "* Global statistics"
                #print >>stath, "** Called variants:\t{0:.0f}\t{1:.2f}%".format(called, float(called * 100 / self.dico_stats['all']['called']) )
                #print >>stath, "Trashed variants:\t{0:.0f}\t{1:.2f}%".format(trashed, float(trashed * 100 / total) )
                print >>stath, "** Total variants:\t{0:.0f}".format(self.dico_stats['all']['called'])
                print >>stath, "** Total called alleles:\t{0:.0f}\t{1:.2f}%".format(self.dico_stats['all']['calledAlleleNumber'], float(self.dico_stats['all']['calledAlleleNumber'] * 100 / self.dico_stats['all']['called']))
                print >>stath, "** Total region size: {0} bases".format(self.dico_stats['positions'])
                globalStat = False
            print >>stath, "Statistics for collection {}:".format(type)
            print >>stath, "* In this collection:"
            print >>stath, "* --------------"
            print >>stath, "* Called variants"
            print >>stath, "** Insertions:\t{0}\t{1:.2f}%".format(self.dico_stats[type]['majorAlleleTiTv']['ins'], (float( float(self.dico_stats[type]['majorAlleleTiTv']['ins']) * 100 / called )) )
            print >>stath, "** Deletions:\t{0}\t{1:.2f}%".format(self.dico_stats[type]['majorAlleleTiTv']['del'], (float( float(self.dico_stats[type]['majorAlleleTiTv']['del']) * 100 / called )) )
            print >>stath, "** Transitions:\t{0}\t{1:.2f}%".format(self.dico_stats[type]['majorAlleleTiTv']['ti'], (float( float(self.dico_stats[type]['majorAlleleTiTv']['ti']) * 100 / called )) )
            print >>stath, "** Transversion:\t{0}\t{1:.2f}%".format(self.dico_stats[type]['majorAlleleTiTv']['tv'], (float( float(self.dico_stats[type]['majorAlleleTiTv']['tv']) * 100 / called )) )
            print >>stath, dividedTiTv
            print >>stath, "** Mutation rate: {0:.5f} events/kb".format( float( called ) * 1000 / float( self.dico_stats['positions'] ) )
            print >>stath, "** Total called variants in this collection: {}".format(called)
            print >>stath, "* --------------"
            print >>stath, "* On all alleles (counting multiple sites at one positions)"
            print >>stath, "** Insertions:\t{0}\t{1:.2f}%".format(self.dico_stats[type]["allAlleleTiTv"]['ins'], (float( float(self.dico_stats[type]["allAlleleTiTv"]['ins']) * 100 / calledAlleleNumber )) )
            print >>stath, "** Deletions:\t{0}\t{1:.2}%".format(self.dico_stats[type]["allAlleleTiTv"]['del'], (float( float(self.dico_stats[type]["allAlleleTiTv"]['del']) * 100 / calledAlleleNumber )) )
            print >>stath, "** Transitions:\t{0}\t{1:.2f}%".format(self.dico_stats[type]["allAlleleTiTv"]['ti'], (float( float(self.dico_stats[type]["allAlleleTiTv"]['ti']) * 100 / calledAlleleNumber )) )
            print >>stath, "** Transversion:\t{0}\t{1:.2f}%".format(self.dico_stats[type]["allAlleleTiTv"]['tv'], (float( float(self.dico_stats[type]["allAlleleTiTv"]['tv']) * 100 / calledAlleleNumber )) )
            print >>stath, dividedAllTiTv
            print >>stath, "** Mutation rate: {0:.5f} events/kb".format( float( calledAlleleNumber ) * 1000 / float( self.dico_stats['positions'] ) )
            print >>stath, "** Total alleles in this collection: {}".format(calledAlleleNumber)
        else: print "No such type: {}".format(type)