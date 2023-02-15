#!/usr/bin/env python2.7
import sys

def fill_TFdic(peaks):
    pkdic={}
    for pk in peaks:
        if pk[0]=="#":
            continue
        else:
            line=pk.rsplit()
            id="-".join(line[0:3])
            pkdic[id]=pk
    pkset=set(pkdic.keys())
    return pkdic, pkset

if __name__=='__main__':
    allPeaks=sys.argv[1]
    alignPeaks=sys.argv[2]
    otpDir=sys.argv[3]

    # define tf and sp from filename
    sp=allPeaks.split("/")[-1].split("_")[1]
    tf=allPeaks.split("/")[-1].split("_")[0]

    # open input files
    allPeaks=open(allPeaks,"r")
    #allPeaks.readline() # skip header line

    alignPeaks=open(alignPeaks,"r")

    # fill set and dictionaries
    allpkdic,allpkset=fill_TFdic(allPeaks)
    alignpkdic,alignpkset=fill_TFdic(alignPeaks)

    # now take the peaks which were not alignable by set subtraction
    notAlignset =allpkset - alignpkset

    # open output file
    otp=otpDir+ "/"+ tf + "_" + sp + "_Liver_AlignmentAndBindingConservation.tab"
    otp=open(otp,"wa")

    # now compose the final line for alignabe and not-alignable pks
    # and store to file
    for pk in allpkset:
        if pk in alignpkset:
            finalLine=allpkdic[pk].rstrip()+"\t"+sp+"\t"+"\t".join(alignpkdic[pk].rsplit()[10:17])
            otp.write(finalLine+"\n")
        elif pk in notAlignset:
        #for pk in notAlignset:
            line=allpkdic[pk].rsplit()
            alignsp=sp
            alignnum="1"
            bindingsp=sp
            bindingnum="1"
            cat="SpSpecificGains"
            label=sp+"SpecificGains"
            binding="bound"
            spfrom=sp
            finalLine=line + [spfrom, alignsp, alignnum, bindingsp, bindingnum, cat, label, binding]
            finalLine="\t".join(finalLine)
            otp.write(finalLine+"\n")
