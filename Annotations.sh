#!/usr/bin/env bash

### HELP

usageHelp="Usage: ${0##*/} [-c config file ]"
configHelp="* [-c config file]: provide config file"
stageHelp="* [-s stage]: if specified, it must be either \"CGI\", \"MethFrag\", \"cluster\" "
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$configHelp"
    echo -e "$stageHep"
    echo "$helpHelp"
    exit $1
}

while getopts "c:s:" opt; do
    case $opt in
	c) config=$OPTARG ;;
    s) stage=$OPTARG ;;
	h) printHelpAndExit 0;;
    ?]) printHelpAndExit 1;;
    esac
done

if [ ! -z "$config" ]
then
    source ${config}
    source ${condaActivate}
    conda activate ${condEnv}
else
    echo "ERROR: config file not provied"
    echo "$usageHelp"
    echo "$configHelp"
    echo "$stageHelp"
    exit 1
fi

if [ ! -z "$stage" ]
then
    case $stage in
        CGI) CGI=1 ;;
        MethFrag) MethFrag=1 ;;
        cluster) cluster=1 ;;
        *)  echo "$usageHelp"
	        echo "$stageHelp"
	        exit 1
    esac
else
    all=1
fi

############# VARIABLES ###################
mkdir -p ${OutDir}/Logs
rm -f ${OutDir}/Logs/*
dt=`date '+%H:%M:%S'`
script=`dirname $(realpath $0)`

################################

if [[ ! -z "$all" ]] || [[ ! -z "$CGI" ]]; then
    mkdir -p ${OutDir}/CGI
    for sp in ${species[@]}; do
        bsub -J"${sp}_cgi_${dt}" -g /cpg  -M 30000 -R "rusage[mem=30000]" -n 3 \
        -oo  ${OutDir}/Logs/${conv[${sp}]}_cgi.out  -e ${OutDir}/Logs/${conv[${sp}]}_cgi.err \
        "cpgplot -sequence ${Genomes}/${conv[${sp}]}/*filtered.fa \
            -window ${window} \
            -minlen ${minlen} \
            -minoe ${minoe} \
            -minpc ${minpc} \
            -outfile ${OutDir}/CGI/${conv[${sp}]}.cpg  \
            -graph ${graph} \
            -outfeat ${OutDir}/CGI/${conv[${sp}]}_feature.txt";
    done
fi

if [[ ! -z "$all" ]] ; then
    mkdir -p ${OutDir}/MethFragmentation
    for sp in in ${species[@]}; do
        for tiss in ${tissues[@]}; do
            bsub -J"${sp}_${tiss}_MF_${dt}" -g /mf -M 20000 -R "rusage[mem=20000]" -n 2 \
                -oo ${OutDir}/Logs/${conv[${sp}]}_${tiss}_MF.out -e ${OutDir}/Logs/${conv[${sp}]}_${tiss}_MF.err \
                -w "done(${sp}_cgi_${dt})" \
                "Rscript ${script}/Ann_methFragmentation.R \
                    ${conv[${sp}]} \
                    ${OutDir}/Meth/${conv[${sp}]}_${tiss}_meth.bed \
                    ${OutDir}/CGI/${conv[${sp}]}_feature.txt"
        done
    done

elif [[ ! -z "${MethFrag}" ]]; then
    mkdir -p ${OutDir}/MethFragmentation
    for sp in ${species[@]}; do
        for tiss in ${tissues[@]}; do
            bsub -J"${sp}_${tiss}_MF_${dt}" -g /mf -M 20000 -R "rusage[mem=20000]" -n 2 \
                -oo ${OutDir}/Logs/${conv[${sp}]}_${tiss}_MF.out -e ${OutDir}/Logs/${conv[${sp}]}_${tiss}_MF.err \
                "Rscript ${script}/Ann_methFragmentation.R \
                    ${conv[${sp}]} \
                    ${OutDir}/Meth/${conv[${sp}]}_${tiss}_meth.bed \
                    ${OutDir}/CGI/${conv[${sp}]}_feature.txt \
                    ${OutDir}/MethFragmentation"
        done
    done
fi
