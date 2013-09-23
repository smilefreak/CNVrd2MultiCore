#!/usr/bin/python3
# Genotypes CSV to Nexus file format
#
#
# url = http://paup.csit.fsu.edu/nfiles.html
#
# 

import os
import sys

from  optparse import OptionParser


def write_nexus(in):
    with open(in + "nex",'w') as input:
        input.write("#

def main():
    parser= OptionParser()
    parser.add_option("-i",'--input',dest="input_file")
    (options,args) = parser.parser_args()
    write_nexus(options.input_file)       


if __name__=="__main__":main()
