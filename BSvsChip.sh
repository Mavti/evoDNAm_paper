#!/usr/bin/env bash

### HELP

usageHelp="Usage: ${0##*/} [-c config file] [-s stage]"
configHelp="* [-c config file]: provide config file as input"
configPolish="* [-p] : remove temporary files and polish folder"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$configHelp"
    echo "$helpHelp"
    exit 1
}

unset polish
while getopts "c:s:p" opt; do
    case $opt in
	c) config=$OPTARG ;;
    p) polish=1 ;;
	h) printHelpAndExit 0;;
    ?]) printHelpAndExit 1;;
    esac
done


if [ ! -z "$config" ]
then
    source $config
else
    echo "ERROR: config file not provied"
    echo "$usageHelp"
    echo "$configHelp"
    exit 1
fi
############# VARIABLES ###################
mkdir -p ${OutDir}/Logs
rm -f ${OutDir}/Logs/*
dt=`date '+%H:%M:%S'`
mkdir -p ${OutDir}/BindingRegions
rm -f ${OutDir}/BindingRegions/*tmp
script=`dirname $(realpath $0)`
################################

if [[ -z "$polish" ]]; then
    for elem in ${comb[@]}; do
        # INTERSECT PEAKS WITH CpGS
        species=$(echo ${elem} | cut -d'_' -f2)


        ###############################
        ## PEAKS THAT OVERLAP WITH CpGs
        ###############################
        ### use ufiltered methylome file
        #echo -e "#chr\tstart\tend\tnCpG\tsummit\tfoldE\tpval\tavgMeth" > ${OutDir}/BindingRegions/${elem}_+CpG_notFilt.bed

        bsub -M 20000 -R "rusage[mem=20000]" -n2 -g /prova -J "${elem}_wCpG_notFilt_${dt}" \
            -oo ${OutDir}/Logs/${elem}_+CpG_notFilt.out -e ${OutDir}/Logs/${elem}_+CpG_notFilt.err \
            "bedtools intersect -a ${PeaksDir}/${elem}.bed -b ${MethDir}/${species}_${tissue}.bed -wa -wb |
            sort -u | sort -k1,1 -k2,2n  > ${OutDir}/BindingRegions/${elem}_+CpG_notFilt.tmp"

        ## use filtered methylome file
        #echo -e "#chr\tstart\tend\tnCpG\tsummit\tfoldE\tpval\tavgMeth" > ${OutDir}/BindingRegions/${elem}_+CpG_filt.bed
        bsub -M 20000 -R "rusage[mem=20000]" -n2 -g /prova -J "${elem}_wCpG_filt_${dt}" \
            -oo ${OutDir}/Logs/${elem}_+CpG_filt.out -e ${OutDir}/Logs/${elem}_+CpG_filt.err \
            "bedtools intersect -a ${PeaksDir}/${elem}.bed -b ${MethDir}/${species}_${tissue}_filtered.bed -wa -wb |
            sort -u | sort -k1,1 -k2,2n  > ${OutDir}/BindingRegions/${elem}_+CpG_filt.tmp"


        #######################################
        ## PEAKS THAT DO NOT OVERLAP WITH CpGs
        #######################################
        ## unfiltered methylome
        #echo -e "#chr\tstart\tend\tnCpG\tsummit\tfoldE\tpval\tavgMeth" > ${OutDir}/BindingRegions/${elem}_-CpG_notFilt.bed

        bsub -M 20000 -R "rusage[mem=20000]" -n2 -g /prova -J "${elem}_woCpG_notFilt_${dt}" \
            -oo ${OutDir}/Logs/${elem}_-CpG_notFilt.out -e ${OutDir}/Logs/${elem}_-CpG_notFilt.err \
            "bedtools intersect -a ${PeaksDir}/${elem}.bed -b ${MethDir}/${species}_${tissue}.bed -v |
                sort -u | sort -k1,1 -k2,2n |
                awk '{ print \$1,\$2,\$3,\$4,\$5,\$6,0,0}' |
                sed 's/ /\t/g' | sort -u > ${OutDir}/BindingRegions/${elem}_-CpG_notFilt.tmp"

        ## filtered methylome
        #echo -e "#chr\tstart\tend\tnCpG\tsummit\tfoldE\tpval\tavgMeth" > ${OutDir}/BindingRegions/${elem}_-CpG_filt.bed

        bsub -M 20000 -R "rusage[mem=20000]" -n2 -g /prova -J "${elem}_woCpG_filt_${dt}" \
            -oo ${OutDir}/Logs/${elem}_-CpG_filt.out -e ${OutDir}/Logs/${elem}_-CpG_filt.err \
            "bedtools intersect -a ${PeaksDir}/${elem}.bed -b ${MethDir}/${species}_${tissue}_filtered.bed -v |
                sort -u | sort -k1,1 -k2,2n |
                awk '{ print \$1,\$2,\$3,\$4,\$5,\$6,0,0}' |
                sed 's/ /\t/g' | sort -u > ${OutDir}/BindingRegions/${elem}_-CpG_filt.tmp"

        ##################
        #### METH AVG
        ##################
        # filtered methylome
        bsub -M 5000 -R "rusage[mem=5000]"  -g /prova -J "${elem}_filt_methavg_${dt}" \
        -oo ${OutDir}/Logs/${elem}_filt_methavg.out -e ${OutDir}/Logs/${elem}_filt_methavg.err \
        -w "done(${elem}_wCpG_filt_${dt})" \
            "python2.7 ${script}/BSvsChip_finalOtp.py  \
                ${elem} \
                ${OutDir}/BindingRegions \
                filt "

        # unfiltered methylome
        bsub -M 5000 -R "rusage[mem=5000]"  -g /prova -J "${elem}_notFilt_methavg_${dt}" \
        -oo ${OutDir}/Logs/${elem}_notFilt_methavg.out -e ${OutDir}/Logs/${elem}_notFilt_methavg.err \
        -w "done(${elem}_wCpG_notFilt_${dt})" \
            "python2.7 ${script}/BSvsChip_finalOtp.py  \
                ${elem} \
                ${OutDir}/BindingRegions \
                notFilt "

        ############################
        ### CONCATENATE ALL PEAKS
        ###########################
        ## filtered
        #echo -e "#chr\tstart\tend\tnCpG\tsummit\tfoldE\tpval\tavgMeth" > ${OutDir}/BindingRegions/${elem}_CpG_filt.bed

        bsub -M 5000 -R "rusage[mem=5000]" -n2 -g /prova -J "${elem}_union_filt_${dt}" \
        -oo ${OutDir}/Logs/${elem}_union_filt.out -e ${OutDir}/Logs/${elem}_union_filt.err \
        -w "done(${elem}_filt_methavg_${dt}) && done(${elem}_woCpG_filt_${dt})" \
        "cat   ${OutDir}/BindingRegions/${elem}_+CpG_filt_avgmeth.tmp ${OutDir}/BindingRegions/${elem}_-CpG_filt.tmp |
            sort -u   > ${OutDir}/BindingRegions/${elem}_CpG_filt.bed"

        bsub -M 5000 -R "rusage[mem=5000]" -n2 -g /prova -J "${elem}_union_notFilt_${dt}" \
        -oo ${OutDir}/Logs/${elem}_union_notFilt.out -e ${OutDir}/Logs/${elem}_union_notFilt.err \
        -w "done(${elem}_notFilt_methavg_${dt}) && done(${elem}_woCpG_notFilt_${dt})" \
        "cat   ${OutDir}/BindingRegions/${elem}_+CpG_notFilt_avgmeth.tmp ${OutDir}/BindingRegions/${elem}_-CpG_notFilt.tmp |
            sort -u   > ${OutDir}/BindingRegions/${elem}_CpG_notFilt.bed"


        ####################################
        ### MERGE FILTERED AND NOT FILTERED
        ####################################
        bsub -M 5000 -R "rusage[mem=5000]" -n2 -g /prova -J "${elem}_final_${dt}" \
        -oo ${OutDir}/Logs/${elem}_final.out -e ${OutDir}/Logs/${elem}_final.err \
        -w "done(${elem}_union_notFilt_${dt}) && done(${elem}_union_filt_${dt})" \
            "rm -f ${OutDir}/BindingRegions/${elem}.bed
            python2.7 ${script}/BSvsChip_finalOtp.py  \
                ${elem} \
                ${OutDir}/BindingRegions \
                final "



        ########################
        ## CpG METH AT BS
        ########################
        #echo -e "#chr\tstart\tend\tstrand\tcov\tM\tavgMeth" > CpGin_${pk}.bed
        bsub -J "CpGs_at_${elem}_${dt}" -M 30000 -R "rusage[mem=30000]" -n2 -g /prova   \
            -oo ${OutDir}/Logs/${elem}_onlyCPG.out -e ${OutDir}/Logs/${elem}_onlyCPG.err \
                "bedtools intersect -b ${OutDir}/BindingRegions/${elem}.bed \
                        -a ${MethDir}/${species}_${tissue}_filtered.bed -wa |
                        sort -u | sort -k1,1 -k2,2n -u > ${OutDir}/BindingRegions/CpGin_${elem}.bed"

    done
else
    echo "polishing folder.."
    bsub -J "tmpRemove"  -g /prova \
        -oo ${OutDir}/Logs/rmTMP.out -e ${OutDir}/Logs/rmTMP.err \
        "rm ${OutDir}/BindingRegions/*tmp"

    bsub -J "iltRemove"  -g /prova \
        -oo ${OutDir}/Logs/rmILT.out -e ${OutDir}/Logs/rmILT.err \
        "rm ${OutDir}/BindingRegions/*ilt.bed"
fi
