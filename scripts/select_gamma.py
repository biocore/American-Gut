#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:13:48 2013

@author: amnon

### 80 char max please

Look at all the gammaproteobacteria and select candidate contamination sequence
    OTUs
output: a list of sorted gammaproteobacteria (or other) otuids, according to 
    mean frequency
"""

import sys
import argparse
import numpy as np
# to load a BIOM table
from biom.parse import parse_biom_table
from biom.util import biom_open

def TestAll(biomfile, outputfile, taxonomyclass, taxonomyname,level):
    """doc string here, a one liner 

    ...and then more detail
    """
    odat=[]
    t = parse_biom_table(biom_open(biomfile,'U'))
    
    t2 = t.normObservationBySample()

    # to iterate over the table by observation, doing something based on the
    # taxonomy:
    class_idx = taxonomyclass
    for values, ids, metadata in t2.iterObservations():
        tname=metadata['taxonomy'][class_idx].lstrip()
        if tname == taxonomyname:
            mv = np.mean(values)
            odat.append((ids,mv))
    
    # odat.sort(key=lambda tup: tup[1], reverse=True)
    odat.sort(key=lambda tup: tup[1])

    csum=[(odat[0][0],odat[0][1],odat[0][1])]
    for cval in odat[1:]:
        csum.append((cval[0],cval[1],csum[-1][2]+cval[1]))
    
    # no get it from big to small
    csum.reverse()
    
    # and write everything above the threshold (to filter)
    snames=open(outputfile,'w')
    for cval in csum:
        if cval[2]>=level:
            snames.write(cval[0]+"\t"+str(cval[1])+"\t"+str(cval[2])+'\n')
    snames.close()

def main(argv):
    parser=argparse.ArgumentParser(description='Select Gammaproteobacteria (or other group) contamination candidates')
    parser.add_argument('-i','--biom',help='biom file of the experiment')
    parser.add_argument('-o','--output',help='output file name')
    parser.add_argument('-c','--classpos',help='class of taxonomy name (0-kingdom,1-phylum etc.',default=2)
    parser.add_argument('-t','--taxonomy',help='taxonomy name (including c__ or equivalent)',default='c__Gammaproteobacteria')
    parser.add_argument('-l','--level',help='minimal cumulative level for OTUs to filter (use 0 to get all of them)',default='0.03')
    
    args=parser.parse_args(argv)
    TestAll(args.biom,args.output,int(args.classpos),args.taxonomy,float(args.level))

if __name__ == "__main__":
    main(sys.argv[1:])          
