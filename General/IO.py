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
'''

# Default library
import sys
from time import clock

# import sqlite3 # Deprecated and never used

def which(program):
    '''Got that on stackoverflow: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python'''
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

###############################################################################
###########################  Colors  ##########################################
# Got this snippet on stackoverflow: http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
###############################################################################

class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = "\033[1m"


###############################################################################
###########################  Time functions  ##################################
### Got this snippet on stackoverflow: http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution ###
###############################################################################

def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
            [(t*1000,),1000,60,60])

def now():
    return secondsToStr(clock())

start = clock()

###############################################################################
#############################    I/O Class    #################################
###############################################################################

def getLineCount(fh):
    '''Gets number of lines in a file
    * Input: <File> FileHandler
    * Output: <int> NumberOfLines'''
    num_lines = sum(1 for line in fh)
    return num_lines

class IO(object):
    """Designed to handle all input/output (e.g. file reading and writing, or writing to standard output)"""
    def __init__(self, ll):
        self.handles = {}
        self.ll = int(ll)
        self.dbs = {}
        self.register("log", sys.stderr, "w")
        self.log("Start Program", header=True, loglevel=0)

    def register(self, type, object, mode='r'):
        '''Gets a File-like object registered in IO object.
        
        This is useful for keeping all our I/O in the same place.
        The 'busy' tag stands here if I ever want to introduce multiprocessing. It will act as a
        semaphore: IO won't give a file handler for which "busy" tag is True using method getIO(type).
        * Input: <str> type (NameOfHandler), <str> PathToFile or <sys.stderr> or <sys.stdout>, <str> mode (see open function), <bool> Busy
        * Output:'''
        if (object in [sys.stdout, sys.stderr]):
            self.handles[type] = {'io': object, 'mode': 'w', 'busy': False }
        else:
            self.handles[type] = {'io': open(object, mode), 'mode': mode, 'busy': False } 

    def unregister(self, type):
        '''Removes a File-like object in I/O, flush it and close it
        
        This method will provide user with the possibility to immediatly close and flush all File-like objects registered.
        It will also provide methods for properly closing any handler IO would have registered.
        
        * Input: <str> type (NameOfHandler)'''
        
        if (type == "all"):
            for ioDict in self.handles.itervalues():
                if ('w' in ioDict['mode'] or ioDict['io'] in [sys.stderr, sys.stdout]):
                    ioDict['io'].flush()
                if (ioDict['io'] not in [sys.stdout, sys.stderr]):
                    ioDict['io'].close()
            self.handles.clear()
        else:
            ioDict = self.handles[type]
            ioDict['io'].flush()
            if (ioDict['io'] not in [sys.stdout, sys.stderr]):
                ioDict['io'].close()
            del self.handles[type]

    def isRegistered(self, type):
        '''Returns True if registered, else False'''
        
        if (type in self.handles.iterkeys()):
            return True
        return False

    def getIO(self, type):
        '''Returns the File-like object for a registered entry or None if it is already used
        
        * Input: <str> type (NameOfHandler)
        * Output: - if not busy: <File> object
        * Output: - if busy: <None>
        '''
        ioDict = self.handles[type]
        if ioDict['busy'] == True:
            return None
        else: 
            ioDict['busy'] = True
            return ioDict['io']

    def giveupIO(self, type):
        '''Allows to switch busy tag from True to False, disabling File-like objects retention
        
        * Input: <str> type (NameOfHandler)'''
        ioDict = self.handles[type]
        if ioDict['busy'] == True:
            ioDict['busy'] = False
            return True
        else:
            return False

    def writeIO(self, type, toWrite):
        '''Used as a print-like function with thread-safety using the "busy" flag
        
        * Input: <str> type (NameOfHandler), <str> toWrite
        * Output: True upon success'''
        ioDict = self.handles[type]
        if 'w' in ioDict['mode'] and not ioDict['busy']:
            handle = ioDict['io']
            ioDict['busy'] = True
            for el in toWrite:
                handle.write('{0}: {1}\n'.format(type, el))
            ioDict['busy'] = False
            return True
        else: return False

    def log(self, s, elapsed=None, header=False, loglevel=0):
        '''Simple logging method that provides log levels.
        
        Log levels are an easy way to level up verbosity. This log method will add a timestamp
        to everything it will report. It is thread-safe thanks to the "busy" flag
        
        * Input: <str> s (toLog), <str> elapsed, <bool> header, <int> loglevel
        '''
        if (self.ll >= loglevel):
            logh = self.getIO('log')
            if (header):
                print >>logh, "="*40
            print >>logh, secondsToStr(clock()), '-', s
            if (elapsed):
                print >>logh, "Elapsed time:", elapsed
            if (header):
                print >>logh, "="*40
                print >>logh, ""
            self.giveupIO('log')
        else:
            pass

    def endlog(self):
        '''This is a peculiar form of the log function called at the end of programs.
        
        This method will display the elapsed time and inform users of good execution'''
        end = clock()
        elapsed = end-start
        self.log("End Program", secondsToStr(elapsed), header=True)

#     def createDB(self, dbName):
#         '''Deprecated and never used'''
#         conn = sqlite3.connect(dbName)
#         self.dbs[dbName] = conn
#         return self.dbs[dbName]
