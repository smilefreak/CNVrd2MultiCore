from optparse import OptionParser
from boto.s3.connection import S3Connection
import os
import logging
import time
import progressbar
#
#  Creates coordinates files for each file in the 1000 genomes data set
#
#
#
#

pbar = None
logging.basicConfig(filename="sample.log", level=logging.INFO)

import pysam
import pickle

# Class representing all the reads and positions for all the reads in the file

class sizeAndPosition(Object):
    def __init__(self,positions,size,chromosomes):
        self._positions=positions
        self._sizes=size
        self._chromosomes=chromosomes
        
    @property
    def sizes(self):
        return self._sizes
    @property
    def positions(self):
        return self._positions
    @property
    def chromosomes(self):
        return self._chromosomes

def read_bams(bam_file):
    return pysam.SamFile(bam_file,'r')


# dumb the coordinates out to file.
def write_coordinates_and_length(bam_file,output_file):
    sizes=[]
    chromosomes=[]
    positions[]
    i = 0 
    for bam in bam_file:
        sizes[i]=int(bam.rlen)
        positions[i]=long(bam.position)
        chromosomes[i]=1
        i = i + 1
    dump=sizeAndPosition(position,sizes,chromosomes)
    with open(bam_file,'wb') as out
        pickle.dump(dump,out)


def open_1kg_connection(access_key,secret_key):
    conn=S3Connection(access_key,secret_key)
    buck=conn.get_bucket('1000genomes')
    return buck

def progress_callback(current,total):
    try:
        pbar.update(current)
    except AssertionError, e:
        logging.error(e)

def download_file(bam_file,key):
    # get_file_size
    global pbar
    size  = key.size
    print(key)
    print(key.size)
    logging.info("File Size " + str(size))
    widgets = [
        unicode(bam_file, errors='ignore').encode('utf-8'), ' ',
        progressbar.FileTransferSpeed(),
        ' <<<', progressbar.Bar(), '>>> ',
        progressbar.Percentage(), ' ', progressbar.ETA()
    ]
    pbar = progressbar.ProgressBar(widgets=widgets, maxval=size)
    pbar.start()
    try:
        key.get_contents_to_filename(bam_file,cb=progress_callback,num_cb=10000)
    except Exception, e:
        logging.error("Failed to Download " + bam_file)
        logging.error(e)
    pbar.finish()

def main():
    parser = OptionParser()
    parser.add_option('-a',dest="access_key")
    parser.add_option('-s',dest="secret_key")
    parser.add_option('-i',dest="sample_list")
    (options,args) = parser.parse_args()
    working_dir = "Temp"
    try:
        os.mkdir(working_dir)
    except:
        print("Directory Already Exists")
    access_key = options.access_key
    secret_key = options.secret_key
    sample_list = options.sample_list
    buck = open_1kg_connection(access_key,secret_key)
    with open(sample_list, 'r') as f:
        i = 0
        line= f.readline() 
        for line in f:
            if( i ==0 ):
                line = line.strip()
                key = buck.get_key(line.strip())
                bam_file = os.path.join(working_dir,os.path.basename(line))
                download_file(bam_file,key)
                # Read bam into R get readlength + position as a 64 bit interger.# 
                read_bams = read_bams(bam_file)

                # This is where it gets a bit obscure 
                # first byte = chromosome number in ASCII
                # second byte = read length
                # third byte = position start of the read mapped to 
                output = os.path.basename(bam_file) + ".cnvb"
                write_coordinates_and_length(bam_file,output)
                i = i + 1
            else:
                break

if __name__=="__main__":main()
