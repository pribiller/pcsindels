"""  Data structure used to parse file from UCSC. 
""" 
from collections import namedtuple

Chain = namedtuple('Chain', 'score chainId tChrom tSize tStrand tStart tEnd qChrom qSize qStrand qStart qEnd')
Block = namedtuple('Block', 'igs tGap qGap')
Pcs   = namedtuple('Pcs', 'size tChrom tStrand tPosBeg qChrom qStrand qPosBeg')

Time  = namedtuple('Time', 'desc t')
