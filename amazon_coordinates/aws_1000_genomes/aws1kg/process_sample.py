from optparse import OptionParser
from boto.s3.connection import S3Connection
from boto.s3.key import Key
import os
import xmlrpclib
import logging
import subprocess
import gzip
import time
import progressbar
import xmlrpclib
#
#  Creates coordinates files for each file in the 1000 genomes data set
#
#
#
#

pbar = None
logging.basicConfig(filename='process_bam.log',level=logging.INFO)

import pysam
import pickle

# Class representing all the reads and positions for all the reads in the file

class sizeAndPosition(object):
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
    print(bam_file)
    return pysam.Samfile(bam_file)


# dumb the coordinates out to file.
def write_coordinates_and_length(bam_file,output_file):
    #sizes=[]
    #chromosomes=[]
    #positions=[]
    #i = 0 
    with open(output_file,'wb') as out:
        for bam in bam_file:
            # Prety sure all the reads are exactly the same size
            # Tell hoang that is why his window approximation works
            #size = bam.rlen
            positions = bam.pos
            chromosomes = bam_file.getrname(bam.rname)
            out.write(str(chromosomes) + '\t' + str(positions)+'\n')
        #sizes.append(int(bam.rlen))
        #positions.append(long(bam.pos))
        #chromosomes.append(bam_file.getrname(bam.rname))
            #out.write(str(chromosomes)  + '\t'+  str(size) + '\t' + str(positions)+'\n')
            #i = i + 1
    #dump=sizeAndPosition(position,sizes,chromosomes)
    #with open(output_file,'wb') as out:
    #    pickle.dump(dump,out)


def open_1kg_connection():

    conn=S3Connection()
    buck=conn.get_bucket('1000genomes')
    return (conn,buck)

def progress_callback(current,total):
    try:
        pbar.update(current)
    except AssertionError, e:
        logging.error(e)

def download_file(bam_file,key):
    # get_file_size
    global pbar
    size  = key.size
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
def gzip_file(output,gzip_out):
    subprocess.check_call(['gzip',output])
def connect_to_xml(server_ip):
    return xmlrpclib.ServerProxy(server_ip)

def main():
    parser = OptionParser()
    #parser.add_option('-a',dest="access_key")
    #parser.add_option('-s',dest="secret_key")
    parser.add_option('-o',dest="output_folder")
    parser.add_option('-x',dest="server_ip_and_port")
    (options,args) = parser.parse_args()
    working_dir = options.output_folder 
    #access_key = options.access_key
    #secret_key = options.secret_key
    server = connect_to_xml(options.server_ip_and_port)
    (conn,buck) = open_1kg_connection()
    print(server.system.listMethods())
    bam=server.get_bam_file()
    i = 0
    while( bam!= False):
        filestart =time.time()
        key = buck.get_key(bam)
        bam_file = os.path.join(working_dir,os.path.basename(bam))
        start_time=time.time()
        download_file(bam_file,key)
        end_time=time.time()
        logging.info("Elapsed time to download bam : {0} = {1}".format(bam_file,str(end_time-start_time)))
        # Read bam into R get readlength + position as a 64 bit interger.# 
        bam=read_bams(bam_file)

        # This is where it gets a bit obscure 
        # first byte = chromosome number in ASCII
        # second byte = read length
        # third byte = position start of the read mapped to 
        output = (bam_file) + ".cnv"
        gzip_out = (output) + ".gz"
        start_time=time.time()
        write_coordinates_and_length(bam,output)
        gzip_file(output,gzip_out)
        end_time=time.time()
        logging.info("Elapsed time to get positions of reads from bam file : {0} = {1}".format(bam_file,str(end_time-start_time)))
        os.remove(bam_file)        
        os.remove(output)
        b =conn.get_bucket('1kg_cnvrd2')
        k=Key(b)
        k.key = os.path.basename(gzip_out)
        k.set_contents_from_filename(gzip_out)
        os.remove(gzip_out)
        bam=server.get_bam_file()
        fileend=time.time()
        logging.info("Full elapsed time for file : {0}  = {1}".format(bam_file,str((fileend-filestart)/60)))
        i = i + 1

if __name__=="__main__":main()
