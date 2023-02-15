#!/usr/bin/env bash



### HELP

usageHelp="Usage: ${0##*/} [-c configFile] [-s stage]"
configHelp="* [-c config]: provide config file as input"
stageHelp="* [-s stage]: must be one of \"EPO\", \"intersect\", \"CrossSpecies\", \"LostPeak\" or \"finalTables\""
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$configHelp"
    echo -e "$stageHelp"
    echo "$helpHelp"
    exit "$1"
}


while getopts "c:s:h" opt; do
    case $opt in
    c) config=$OPTARG ;;
	s) stage=$OPTARG ;;
	h) printHelpAndExit 0;;
    ?]) printHelpAndExit 1;;
    esac
done

if [ ! -z "${config}" ]
then
    source ${config}
    source ${condaActivate}
    conda activate ${condEnv}
else
    echo "ERROR: config file not provied"
    printHelpAndExit 1
fi
###
if [ ! -z "$stage" ]
then
    case $stage in
        EPO)       epo=1 ;;
        intersect)    intersect=1 ;;
        CrossSpecies) CrossSpecies=1 ;;
        LostPeak)  LostPeak=1 ;;
        finalTables)    finalTables=1 ;;
        CpGenrich) CpGenrich=1 ;;
        motifsUnbound) motifsUnbound=1 ;;
        RSscore) RSscore=1 ;;
        *)  echo "$usageHelp"
	        echo "$stageHelp"
	        exit 1
    esac
else
    all=1
fi


#### Variables and Paths
dt=`date '+%d/%m/%Y_%H:%M:%S'`
Scripts=`dirname $(realpath $0)`
mkdir -p ${otpDir}/Logs

####################################

# MAP ALL THE PEAKS ONTO ALL SPECIES GENOMES
if [[ ! -z "$all" ]] || [[ ! -z "$epo" ]]
then
    rm -f ${otpDir}/Logs/EPO_*
    conda deactivate
    mkdir -p ${otpDir}/TFalignment
    rm -f ${otpDir}/TFalignment/*
    for sp in ${!TFalign[@]}; do
        for tf in ${TFs[@]}; do
            for tosp in ${TFalign[${sp}]}; do
                echo ${tf}_${sp}_${tosp}
                bsub -J "EPO_${tf}_${sp}_${tosp}_${dt}" -g /tfdiv \
                    -oo ${otpDir}/Logs/epo_${tf}_${sp}_to_${tosp}.out -e ${otpDir}/Logs/epo_${tf}_${sp}_to_${tosp}.err \
                    "perl ${Scripts}/TFbindingPairwiseConservation.pl \
                        -quer_name ${sp} \
                        -quer_bed ${BRs}/${tf}_${nick[$sp]}_Liver.bed  \
                        -ref_name ${tosp} \
                        -out  ${otpDir}/TFalignment/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt"
            done
        done
    done

fi


## INTERSECT ALIGNED PEAKS
# Here I intersect the TFBR projections with the TFBR of the genome we have projocted them onto.
# example: we intersect Mouse's CEBPA binding regions projected onto Macaque with Macaque's CEBPA binding regions.
# at the end we want a table which tells us how many species the TFBRs are alignable to and conserved.
# from the same table we will also extract the "lost peak": region where the TFBR is alignable to but not conserved.
if [[ ! -z "$all" ]]
then
    echo ciao1
    conda activate ${condEnv}
    rm -f ${otpDir}/Logs/cs_*
    mkdir -p ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs
    rm -f ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs/*
    rm -f ${otpDir}/CrossSpecies/fileList.txt

    for sp in ${!TFalign[@]}; do
        for tf in ${TFs[@]}; do
            for tosp in ${TFalign[${sp}]}; do
                echo ${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt >> ${otpDir}/CrossSpecies/fileList.txt
                ## by sorting -u on columns 4 and 5, if the TFBRs was split into more blocks, we randomly pick one of them.
                ## by sorting -u on columns 1 and 2, if the TFBRs overlap with multiple bound regions you pick just one of them.
                bsub -J "CrossSpecies_${tf}_${sp}_${tosp}_${dt}" -g /tfdiv -M 2000 -R "rusage[mem=2000]"  \
                    -oo ${otpDir}/Logs/cs_${tf}_${sp}_to_${tosp}.out -e ${otpDir}/Logs/cs_${tf}_${sp}_to_${tosp}.err \
                    -w "done(EPO_*_${dt})" \
                "set -euo pipefail
                sort -k4,4 -k5,5n -u ${otpDir}/TFalignment/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt | sort -k1,1 -k2,2n |
                bedtools intersect -a stdin -b ${BRs}/${tf}_${nick[$tosp]}_Liver.bed -wa -wb |
                sort -u | sort -k1,1 -k2,2n -u > ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt

                sort -k4,4 -k5,5n -u ${otpDir}/TFalignment/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt | sort -k1,1 -k2,2n |
                bedtools intersect -a stdin -b ${BRs}/${tf}_${nick[$tosp]}_Liver.bed -v |
                sort -u | sort -k1,1 -k2,2n |
                awk '{OFS=\"\t\"; print \$0, \"none\"}' >> ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt"
            done
        done
    done


elif [[ ! -z "$intersect" ]]
then
    rm -f ${otpDir}/Logs/cs_*
    mkdir -p ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs
    rm -f ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs/*
    rm -f ${otpDir}/CrossSpecies/fileList.txt
    for sp in ${!TFalign[@]}; do
        for tf in ${TFs[@]}; do
            for tosp in ${TFalign[${sp}]}; do
                echo ${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt >> ${otpDir}/CrossSpecies/fileList.txt
                ## by sorting -u on columns 4 and 5, if the TFBRs was split into more blocks, we randomly pick one of them.
                ## by sorting -u on columns 1 and 2, if the TFBRs overlap with multiple bound regions you pick just one of them.
                bsub -J "CrossSpecies_${tf}_${sp}_${tosp}_${dt}" -g /tfdiv -M 2000 -R "rusage[mem=2000]"  \
                    -oo ${otpDir}/Logs/cs_${tf}_${sp}_to_${tosp}.out -e ${otpDir}/Logs/cs_${tf}_${sp}_to_${tosp}.err \
                "set -euo pipefail
                sort -k4,4 -k5,5n -u ${otpDir}/TFalignment/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt | sort -k1,1 -k2,2n |
                bedtools intersect -a stdin -b ${BRs}/${tf}_${nick[$tosp]}_Liver.bed -wa -wb |
                sort -u | sort -k1,1 -k2,2n -u > ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt

                sort -k4,4 -k5,5n -u ${otpDir}/TFalignment/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt | sort -k1,1 -k2,2n |
                bedtools intersect -a stdin -b ${BRs}/${tf}_${nick[$tosp]}_Liver.bed -v |
                sort -u | sort -k1,1 -k2,2n |
                awk '{OFS=\"\t\"; print \$0, \"none\"}' >> ${otpDir}/CrossSpecies/ProjectionsIntersectionWithTFBRs/${tf}_${nick[${sp}]}_to_${nick[${tosp}]}.txt"
            done
        done
    done

fi

# MAKE FINAL TABLE AND INTERSECT LOSTPEAK WITH METHYLOME

# DEFINE EVOLUTIONARY CATEGORIES OF BINDING CONSERVATION
# and extract LOST PEAKs
if [[ ! -z "$CrossSpecies" ]]
then
    rm -f ${otpDir}/Logs/cs_table*
    mkdir -p ${otpDir}/CrossSpecies/LostPeak
    rm -f ${otpDir}/CrossSpecies/LostPeak/*tmp1
    mkdir -p ${otpDir}/CrossSpecies/Tables
    rm -f ${otpDir}/CrossSpecies/Tables/*AlignmentAndBindingConservation_alignedPeaks.txt

    # NB: for losses, check if alignment is there when binding is not
    # lineage conserved: binding in lineage and not in other species
    # lineage loss: binding in others and not in lineage + alignment in lineage
    # ultraconserved: binding (and therefore alignment) conserved in all species
    # spSpecific: alignment in all, binding in one
    # spSpecificLoss: binding in 4 (no need to have alignment in 5)
    bsub -J "CrossSpecies_table_${dt}" -g /tfdiv -M 5000 -R "rusage[mem=5000]" -n2  \
        -oo ${otpDir}/Logs/cs_table.out -e ${otpDir}/Logs/cs_table.err \
        "python2.7 ${Scripts}/TFbindingDivergence_CrossSpecies.py ${otpDir}/CrossSpecies/fileList.txt ${otpDir}/CrossSpecies"
fi

# INTERSECT LOST PEAK WITH METHYLOME
if [[ ! -z "$LostPeak" ]]
then
    rm -f ${otpDir}/Logs/lp_*
    rm -f ${otpDir}/CrossSpecies/LostPeak/*LostPeak_avgMeth.txt

    for lp in `ls ${otpDir}/CrossSpecies/LostPeak/*tmp1`; do
        sp=$(echo ${lp} |  rev | cut -d"/" -f1 | rev | cut -d"_" -f2 )
        tf=$(echo ${lp} |  rev | cut -d"/" -f1 | rev | cut -d"_" -f1 )
        from=$( basename $(echo ${lp} |  rev | cut -d"/" -f1 | rev | cut -d"_" -f5) .tmp1 )
        bsub -J "lp_${tf}_${sp}_${from}_${dt}" -g /tfdiv -M 20000 -R "rusage[mem=20000]" \
            -oo ${otpDir}/Logs/lp_${tf}_${sp}_${from}.out -e ${otpDir}/Logs/lp_${tf}_${sp}_${from}.err \
        "set -o pipefail
        sort -k1,1 -k2,2n ${lp} |
        bedtools intersect -a stdin -b ${Meth}/${sp}_Liver.bed -wa -wb | sort -u |
            awk  'BEGIN{OFS=\"\t\"};{
                id=\$1 OFS \$2 OFS \$3;
                tf=\$4;
                sp=\$5;
                from=\$6;
                details[id]=\$7 OFS \$8 OFS \$9 OFS \$10 OFS \$11 OFS \$12 OFS \$13;
                tot[id]++;
                sum[id]=sum[id]+\$17
            };END{
                for (i in tot) print i, \"NA\",\"NA\",\"NA\", tot[i], sum[i]/tot[i], details[i], tf, sp, from}' > ${lp}.tmp2
        [[ $? -eq 0 ]] || exit 1

        sort -k1,1 -k2,2n ${lp} |
        bedtools intersect -a stdin -b ${Meth}/${sp}_Liver.bed -v | sort -u |
            awk 'BEGIN{OFS=\"\t\"};{
                tf=\$4;
                sp=\$5;
                from=\$6;
                print \$1, \$2, \$3, \"NA\",\"NA\",\"NA\", 0, 0, \$7 , \$8 , \$9 , \$10 , \$11 , \$12, \$13, tf, sp, from}' >> ${lp}.tmp2
        [[ $? -eq 0 ]] || exit 1

        sort -k1,1 -k2,2n ${lp}.tmp2 |
        bedtools intersect -a stdin -b ${Meth}/${sp}_Liver_filtered.bed -wa -wb | sort -u |
            awk 'BEGIN{OFS=\"\t\"};{
                id=\$1 OFS \$2 OFS \$3;
                details[id]=\$4 OFS \$5 OFS \$6 OFS \$7 OFS \$8 OFS \$9 OFS \$10 OFS \$11 OFS \$12 OFS \$13 OFS \$14 OFS \$15 OFS \$16 OFS \$17 OFS \$18;
                tot[id]++;
                sum[id]=sum[id]+\$22
            };END{
                for (i in tot) print i, details[i], tot[i], sum[i]/tot[i]}' |
        awk 'BEGIN{OFS=\"\t\"};{
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$19,\$8,\$20, \$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17, \$18
        }' > ${otpDir}/CrossSpecies/LostPeak/${tf}_${sp}_from_${from}_LostPeak_avgMeth.txt
        [[ $? -eq 0 ]] || exit 1

        sort -k1,1 -k2,2n ${lp}.tmp2 |
        bedtools intersect -a stdin -b ${Meth}/${sp}_Liver.bed -v | sort -u |
            awk 'BEGIN{OFS=\"\t\"};{
                print \$0, 0, 0}' |
        awk 'BEGIN{OFS=\"\t\"};{
            print  \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$19,\$8,\$20, \$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17, \$18
        }' >> ${otpDir}/CrossSpecies/LostPeak/${tf}_${sp}_from_${from}_LostPeak_avgMeth.txt
        [[ $? -eq 0 ]] || exit 1"



    done

fi

# This chunk extracts the binding regions that were not found alignable
# and compiles a line for them to be added to the final AlignmentAndBindingConservation table.
# Then, this is stored final tables folder and it's combined with final table cluster files.
if [[ ! -z "${finalTables}" ]]; then
    rm -f ${otpDir}/Logs/ft_*
    rm -f ${Tab}/*AlignmentAndBindingConservation.tab

    for elem in ${comb[@]}; do
        bsub -J "ft_${elem}_${dt}" -g /tfdiv -M 5000 -R "rusage[mem=5000]" \
            -oo ${otpDir}/Logs/ft_${elem}.out -e ${otpDir}/Logs/ft_${elem}.err \
        "python2.7 ${Scripts}/TFbindingDivergence_finalTables.py \
            ${Tab}/${elem}_Liver.tab \
            ${otpDir}/CrossSpecies/Tables/${elem}_Liver_AlignmentAndBindingConservation_alignedPeaks.txt \
            ${Tab}"
    done

    bsub -J "cat_${dt}" -g /tfdiv \
        -oo ${otpDir}/Logs/ft_cleanup.out -e ${otpDir}/Logs/ft_cleanup.err \
        -w "done(ft_*_${dt})" \
        "cat ${Tab}/*AlignmentAndBindingConservation.tab > ${Tab}/AllTFs_allSpecies_AlignmentAndBindingConservation.tab "

fi
