#!/usr/bin/env python

import sys


if __name__=="__main__":
    table=open(sys.argv[1],"r")
    clusters=open(sys.argv[2], "r")
    otp=open(sys.argv[3],"wa")

    final_dic={}
    hd=table.readline()
    hd=hd.rstrip()+"\tcluster_num\tcluster_name\n"
    for line in table:
        line=line.rsplit()
        final_dic[line[3]]=line

    cl_dic={}
    for line in clusters:
        #print line
        line=line.rsplit()
        cl_dic[line[1]]=[line[3],line[4]]

    otp.write(hd)
    for elem in final_dic.keys():
        final_dic[elem].extend(cl_dic.get(elem,["NC","NC"]))
        otpline="\t".join(final_dic[elem])
        otp.write(otpline+"\n")
