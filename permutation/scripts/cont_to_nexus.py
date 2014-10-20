import os
import sys
from optparse import OptionParser

def write_nexus(input_f):
    with open(input_f,'r') as input_file:
        with open(input_f + ".nex",'w') as output:
            output.write("#NEXUS\n\n")
            output.write("begin data; \n")
            continuous_list=[]
            varieties=[]
            for line in input_file:
                    line=line.split()
                    continuous_list.append(line[1:])
                    varieties.append(line[0].strip('"'))
            output.write("Dimensions ntax={0} nchar={1};\n".format(len(varieties),len(continuous_list[1])))
            output.write("""format datatype=Continuous interleave=no missing=?;\n""")        
            output.write('Matrix\n')
            for i, variety in enumerate(varieties):
                output.write(variety + '\t')
                output.write(' '.join(continuous_list[i]))
                output.write('\n')
            output.write(';\n')
            output.write('end;\n')
            
 


def main():
    parser = OptionParser()
    parser.add_option('-i',dest="input_file")
    (options,args) = parser.parse_args()
    write_nexus(options.input_file) 



if __name__=="__main__":main()
