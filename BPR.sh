#!/usr/bin/env bash

### HELP

usageHelp="Usage: ${0##*/} [-c config file ]"
configHelp="* [-c config file]: provide config file"
stageHelp="* [-s stage]: if specified, it must be either \"input\", \"parameters\", \"clusters\" or \"features\" "
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
        input) input=1 ;;
        parameters) parameters=1 ;;
        clusters) clusters=1 ;;
        features) features=1 ;;
        predictions) predictions=1 ;;
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

mkdir -p ${OutDir}/Meth2
mkdir -p ${OutDir}/Peaks2
mkdir -p ${OutDir}/ClustersVB/ParamSelection3
################################


if [[ ! -z "$all" ]] || [[ ! -z "$input" ]]; then
    for elem in ${comb[@]}; do
        # REFORMAT METHYLOME
        sptis=$( echo ${elem} | cut -d"_" -f2-3 )
        bsub -M 5000  -R "rusage[mem=5000]" -g /bpr -J"${sptis}_meth_${dt}" \
            -oo ${OutDir}/Logs/${sptis}_meth.out -e ${OutDir}/Logs/${sptis}_meth.err \
            " awk 'BEGIN{OFS=\"\t\"};{
                print \$1, \$2, \$5, \$7
            }' ${MethDir}/${sptis}.bed > ${OutDir}/Meth2/${sptis}_meth.bed"

        # REFORMAT PEAKS
        tfsptis=$( echo ${elem} | cut -d"_" -f1-3 )
        bsub -M 5000  -R "rusage[mem=5000]" -g /bpr -J"${elem}_peaks_${dt}" \
            -oo ${OutDir}/Logs/${tfsptis}_pk.out -e ${OutDir}/Logs/${tfsptis}_pk.err \
            "awk 'BEGIN{OFS=\"\t\"};{
            if(\$1!=\"#chr\")
                print \$1, \$2, \$3, \"*\", \$4
         }'  ${PeaksDir}/${tfsptis}.bed |
            sort -k1,1 -k2,2n > ${OutDir}/Peaks2/${elem}.bed"
    done
fi


if [[ ! -z "$all" ]] ; then
    for elem in ${comb[@]}; do
        sptis=$( echo ${elem} | cut -d"_" -f2-3 )
        bsub -M 15000 -R "rusage[mem=15000]" -n 7 -g /bpr -q "long" -J"${elem}_VB_KandM_${dt}" \
            -w "done(${sptis}_meth_${dt})" \
            "Rscript ${script}/BPR_clustersVB_ParamSelection.R \
                ClustersVB/ParamSelection3 \
                ${OutDir}/Meth2/${sptis}_meth.bed \
                ${OutDir}/Peaks2/${elem}.bed \
                ${elem} \
                ${elem}_${n}_MandKevolution_VB.RData \
                ${n}"
    done
elif [[ ! -z "$parameters" ]] ; then
    for elem in ${comb[@]}; do
        sptis=$( echo ${elem} | cut -d"_" -f2-3 )
        bsub -M 15000 -R "rusage[mem=15000]" -n 7 -g /bpr -q "long" -J"${elem}_VB_KandM_${dt}" \
            "Rscript ${script}/BPR_clustersVB_ParamSelection.R \
                ${OutDir}/ClustersVB/ParamSelection3 \
                ${OutDir}/Meth2/${sptis}_meth.bed \
                ${OutDir}/Peaks2/${elem}.bed \
                ${elem} \
                ${elem}_${n}_MandKevolution_VB.RData \
                ${n}"
    done
fi

if [[ ! -z "$clusters" ]]; then
    mkdir -p ${OutDir}/ClustersVB/Clusters
    for elem in ${!param[@]}; do
        M=$(echo ${param[$elem]} | cut -d" " -f1)
        K=$(echo ${param[$elem]} | cut -d" " -f2)
        bsub -J"${elem}_clusters_${dt}" -M 5000 -R "rusage[mem=5000]" -g /bpr \
            -oo ${OutDir}/Logs/${elem}_clusters.out -e ${OutDir}/Logs/${elem}_clusters.err \
            "Rscript ${script}/BPR_extractClusters.R \
                ${OutDir}/ClustersVB/ParamSelection3/${elem}_${n}_MandKevolution_VB.RData \
                ${M} \
                ${K} \
                ${OutDir}/ClustersVB/Clusters "

        bsub -J"${elem}_rename_${dt}" -M 5000 -R "rusage[mem=5000]" -g /bpr \
            -oo ${OutDir}/Logs/${elem}_rename.out -e ${OutDir}/Logs/${elem}_rename.err \
            -w "done(${elem}_clusters_${dt})" \
            "awk -v name1=${clust[${elem}_1]} \
                -v name2=${clust[${elem}_2]} \
                -v name3=${clust[${elem}_3]} \
                -v name4=${clust[${elem}_4]} \
                'BEGIN{OFS=\"\t\"; getline};
                {if(\$4==1) \$5=name1;
                else if(\$4==2) \$5=name2;
                else if(\$4==3) \$5=name3;
                else if(\$4==4) \$5=name4;
                print \$0}' ${OutDir}/ClustersVB/Clusters/${elem}_clusters_tmp.txt > ${OutDir}/ClustersVB/Clusters/${elem}_clusters.txt"

        bsub -J"${elem}_table_${dt}" -M 5000 -R "rusage[mem=5000]" -g /bpr \
            -oo ${OutDir}/Logs/${elem}_table.out -e ${OutDir}/Logs/${elem}_table.err \
            -w "done(${elem}_rename_${dt})" \
            "cd ${OutDir}
            python2.7 ${script}/BPR_appendClusters.py \
                ../BSvsChip/BindingRegions/${elem}.bed \
                ${OutDir}/ClustersVB/Clusters/${elem}_clusters.txt \
                ${OutDir}/ClustersVB/Final_Tables/${elem}.bed"

    done

fi
