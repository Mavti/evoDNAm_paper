#!/usr/bin/env python2.7
import sys

def fill_dics(inp,tosp,orisp,tf,dic,outDir,lp):
    index_align={"Human":10, "Macaque":11,"Mouse":12,"Rat":13, "Dog":14}
    index_conserved={"Human":15, "Macaque":16,"Mouse":17,"Rat":18, "Dog":19}

    for line in inp:
        line=line.rsplit()
        ####
        tosp_chr=line[0]
        tosp_start=line[1]
        tosp_end=line[2]
        ori_chr=line[3]
        ori_start=line[4]
        ori_end=line[5]
        summit=line[6]
        foldE=line[7]
        pval=line[8]
        nCpGtot=line[9]
        nCpGfilt=line[10]
        avgMethtot=line[11]
        avgMethfilt=line[12]

        id=tf+"_"+orisp+"_"+ori_chr+"_"+ori_start+"_"+ori_end
        details=[ori_chr,ori_start,ori_end, summit,foldE,pval,nCpGtot,nCpGfilt,avgMethtot,avgMethfilt,"0","0","0","0","0","0","0","0","0","0"]
        dic[id]=dic.get(id,details)

        if len(line)==23:
            tosppk_chr=line[13]
            tosppk_start=line[14]
            tosppk_end=line[15]
            tosppk_summit=line[16]
            tosp_foldE=line[17]
            tosp_pval=line[18]
            tosp_ncpgtot=line[19]
            tosp_ncpgfilt=line[20]
            tosp_avgmethtot=line[21]
            tosp_avgmethfilt=line[22]

            dic[id][index_align[orisp]]="1"
            dic[id][index_align[tosp]]="1"
            dic[id][index_conserved[orisp]]="1"
            dic[id][index_conserved[tosp]]="1"

        else:
            dic[id][index_align[orisp]]="1"
            dic[id][index_align[tosp]]="1"
            dic[id][index_conserved[orisp]]="1"
            lp[id]=lp.get(id,[])
            finalLine=[tosp_chr,tosp_start,tosp_end,tf,tosp, orisp]
            finalLine="\t".join(finalLine)
            lp[id].append(finalLine)

    return dic, lp

def make_table(line, otp, orisp,tf):
    index={1:"Human",2:"Macaque",3:"Mouse",4:"Rat",5:"Dog"}

    chr=line[0]
    start=line[1]
    end=line[2]
    summit=line[3]
    foldE=line[4]
    pval=line[5]
    ncpgtot=line[6]
    ncpgfilt=line[7]
    avgmethtot=line[8]
    avgmethfilt=line[9]

    # summarise align and cons
    align=""
    cons=""
    tmp=0
    for ind in range(10,15):
        tmp+=1
        if line[ind]=="1":
            align=align+index[tmp]+","
    align=align[0:-1]
    alignNum=len(align.split(","))

    tmp=0
    for ind in range(15,20):
        tmp+=1
        if line[ind]=="1":
            cons=cons+index[tmp]+","
        else:
            loss=index[tmp]
    cons=cons[0:-1]
    consNum=len(cons.split(","))
    ######

    # get the evo cat
    ## from 10-15 alignability
    ## 15-20 binding conservation
    if tf=="CTCF" or tf=="CEBPA" or tf=="HNF4a":
        if int(line[15])+int(line[16])==2 and int(line[17])+int(line[18])+int(line[19])==0:
            cat="LineageGains"
            evo="PrimateGains"
        elif int(line[15])+int(line[16])+int(line[19])==0 and int(line[17])+int(line[18])==2:
            cat="LineageGains"
            evo="RodentGains"
        elif int(line[15])+int(line[16])+int(line[17])+int(line[18])+int(line[19])==5:
            cat="UltraCons"
            evo="UltraCons"
        elif int(line[15])+int(line[16])+int(line[19])==3 and int(line[17])+int(line[18])==0:
            # and int(line[12])+int(line[13])==0
            cat="LineageLoss"
            evo="RodentLoss"
        elif int(line[17])+int(line[18])+int(line[19])==3 and int(line[15])+int(line[16])==0:
            # and int(line[10])+int(line[11])==0
            cat="LineageLoss"
            evo="PrimateLoss"
        elif int(line[15])+int(line[16])+int(line[17])+int(line[18])+int(line[19])==1:
            # and int(line[10])+int(line[11])+int(line[12])+int(line[13])+int(line[14])==5
            cat="SpSpecificGains"
            evo=cons+"SpecificGains"
        elif int(line[15])+int(line[16])+int(line[17])+int(line[18])+int(line[19])==4:
            cat="SpSpecificLoss"
            evo=loss+"SpecificLoss"
        else:
            cat="evoDynamic"
            evo="evoDynamic"
    elif tf=="FoxA1": # no Macaque
        if int(line[15])+int(line[16])==2 and int(line[17])+int(line[18])+int(line[19])==0: # never going to happen
            cat="LineageGains"
            evo="PrimateGains"
        elif int(line[15])+int(line[16])+int(line[19])==0 and int(line[17])+int(line[18])==2: # ok, Macaque always 0
            cat="LineageGains"
            evo="RodentGains"
        elif int(line[15])+int(line[17])+int(line[18])+int(line[19])==4: #mod
            cat="UltraCons"
            evo="UltraCons"
        elif int(line[15])+int(line[19])==2 and int(line[17])+int(line[18])==0: #mod
            # and int(line[12])+int(line[13])==0
            cat="LineageLoss"
            evo="RodentLoss"
        elif int(line[15])+int(line[17])+int(line[18])+int(line[19])==1:
            # and int(line[10])+int(line[11])+int(line[12])+int(line[13])+int(line[14])==5
            cat="SpSpecificGains"
            evo=cons+"SpecificGains"
        elif int(line[15])+int(line[17])+int(line[18])+int(line[19])==3: #mod
            cat="SpSpecificLoss"
            evo=loss+"SpecificLoss"
        else:
            cat="evoDynamic"
            evo="evoDynamic"
    elif tf=="HNF6": # no Dog
        if int(line[15])+int(line[16])==2 and int(line[17])+int(line[18])==0:
            cat="LineageLoss"
            evo="RodentLoss"
        elif int(line[15])+int(line[16])==0 and int(line[17])+int(line[18])==2:
            cat="LineageLoss"
            evo="PrimateLoss"
        elif int(line[15])+int(line[16])+int(line[17])+int(line[18])==4: # mod
            cat="UltraCons"
            evo="UltraCons"
        elif int(line[15])+int(line[16])+int(line[17])+int(line[18])==1: #mod
            # and int(line[10])+int(line[11])+int(line[12])+int(line[13])+int(line[14])==5
            cat="SpSpecificGains"
            evo=cons+"SpecificGains"
        elif int(line[15])+int(line[16])+int(line[17])+int(line[18])==3:
            cat="SpSpecificLoss"
            evo=loss+"SpecificLoss"
        else:
            cat="evoDynamic"
            evo="evoDynamic"
        #######

    finalLine=line[0:10]+[align, str(alignNum), cons, str(consNum), cat, evo, "bound", tf, orisp, orisp]
    finalLine="\t".join(finalLine)

    otp=open(otp, "a")
    otp.write(finalLine+"\n")
    return finalLine


if __name__=='__main__':
    filelistpath=sys.argv[1]
    outDir=sys.argv[2]

    # dic is a dictionary whose keys are TFBRs and values are details of TFBRs with
    ## info about alignability and binding conservation across all species.
    dic={}
    # lp is a dictionary whose keys are TFBRs and values are details about all
    ## lostPeaks associated to that TFBRs
    lp={}

    filelist=open(filelistpath, "r")
    path=filelistpath.split("/")
    for fl in filelist:
        fl=fl.rstrip()
        fl="ProjectionsIntersectionWithTFBRs/"+fl
        path[-1]=fl

        name=fl.split("/")[-1]
        tosp=name.split("_")[3].split(".")[0]
        tf=name.split("_")[0]
        orisp=name.split("_")[1]

        inp="/".join(path)
        inp=open(inp, "r")
        dic,lp=fill_dics(inp,tosp,orisp,tf,dic,outDir,lp)
    ###

    # once dic and lp have been filled with info from all species,
    ## define evo categories and compose final lines
    for pk in dic:
        tf=pk.split("_")[0]
        orisp=pk.split("_")[1]
        name=outDir+"/Tables/"+tf+"_"+orisp+"_Liver_AlignmentAndBindingConservation_alignedPeaks.txt"

        # make_table summarises info about alignability and binding conservation,
        ## plus defines evolutionary categories
        tableLine=make_table(dic[pk], name,orisp,tf)
        tableLine=tableLine.split("\t")
        if pk in lp:
            # now takes care of all the lost peaks associated to TFBRs
            for lost in lp[pk]:
                sp=lost.split("\t")[4]
                otpLine=lost + "\t" + "\t".join(tableLine[10:-4])+"\t"+ "unbound"
                name=outDir+"/LostPeak/"+tf+"_"+sp+"_lostPeak_from_"+orisp+".tmp1"
                otp=open(name,"a")
                otp.write(otpLine+"\n")
