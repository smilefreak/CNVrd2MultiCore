from optparse import OptionParser
from boto.s3.connection import S3Connection
import os
#
#  Creates coordinates files for each file in the 1000 genomes data set
#
#
#
#

def open_1kg_connection(access_key,secret_key):
    conn=S3Connection(access_key,secret_key)
    buck=conn.get_bucket('1000genomes')
    return buck

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
                key = buck.get_key(line.strip())
                bam_file = os.path.join(working_dir,os.path.basename(line)))
                key.get_contents_to_filename(bam_file)
                
                i = i + 1
            else:
                break

if __name__=="__main__":main()
