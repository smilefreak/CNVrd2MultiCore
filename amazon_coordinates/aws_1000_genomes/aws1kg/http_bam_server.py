from optparse import OptionParser
import time
import xmlrpclib
from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler
from Queue import Queue

#Simple xml server for hosting and running the server.

import logging
logging.basicConfig(level=logging.INFO)

class RequestHandler(SimpleXMLRPCRequestHandler):
    rpc_paths =('/RPC2',)


def main():
    parser = OptionParser()
    parser.add_option('-l','--bam-list',dest="bam_list")
    (options,args) = parser.parse_args()
    bams = Queue()
    worker_threads = 1
    with open(options.bam_list,'r') as input_file:
        print(options.bam_list)
        for line in input_file:
            bams.put(line.strip())
    print("Queue Length {0} ".format(bams.qsize()))
    def get_bam_file():
        if(bams.qsize() == 0):
            logging.info("Finished processing the bam list")
            return False
        bam = bams.get()
        logging.info("Bam file : {0}".format(bam))
        bams.task_done()
        return bam
    server = SimpleXMLRPCServer(("0.0.0.0",8000),
                                requestHandler=RequestHandler)
    server.register_introspection_functions()
    server.register_function(get_bam_file)
    logging.info("Starting Server")
    server.serve_forever()
    
if __name__=="__main__":main() 
     
