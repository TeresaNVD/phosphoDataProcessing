#!/bin/bash
clear
#RUN this script at the root of the project folder (level 1)
printf "##########################################################\n"
printf "###Execution of Phospho(STY)Sites.txt analysis pipeline###\n"
printf "##########################################################\n"
printf "\n"

printf "1. Structure and files under the project's folder:\n"
#root for absolute path building in the scripts (parameter 1); PF: project folder.
#relative path can be achieved by making PF="."
PF=$PWD
tree $PF
printf "\n"

printf "2. Which Python version is used by default:\n"
python -c 'import sys; print(sys.version)'
printf "\n"

printf "3. Which scripts are available:\n"
ls -l

printf "4. Summarizing information with exp_summary.py:\n"
python ./scripts/exp_summary.py $PF
printf "Done!\n"
printf "Quick summary of results:\n"
for i in $PF/results/*.tsv; do
            echo item: $i
            head -3 $i
            echo "------"
        done

printf "5. Extracting phosphosite information.py:\n"
printf "Analysis will be performed for three main modes: psp, iptmnet, classical\n"
site_mode=("psp" "iptmnet" "networkin" "classical")
for i in ${site_mode[@]}; do
            echo item: $i
            python ./scripts/phosphosite.py $PF $i
            printf "Quick summary of results:\n"
            head -3 $PF/results/sites_$i.tsv
            echo "End!"
        done
printf "Done!\n"

printf "6. Extracting motifs from sequences:\n"
printf "Analysis will be performed for two main modes of sequence window:\n"
printf "length 13 (sequence6) and 15 (sequence7)\n"
site_mode=("sequence6" "sequence7")
for i in ${site_mode[@]}; do
            echo item: $i
            python ./scripts/motif.py $PF $i
            printf "Quick summary of results:\n"
            head -3 $PF/results/kinase_motif_scan_$i.tsv
            echo "End!"
        done
printf "If interested in CK2 consensus motif the number of sequences containing\n"
printf "CK2 motif in the vicinity of the identified phosphosites is:\n"
for i in ${site_mode[@]}; do
            echo item: $i
            cat "$PF/results/ck2_hit_$i.tsv" | wc -l
            echo "End!"
        done
printf "Done!\n"

printf "7. Extracting motifs from sequences:\n"
printf "Analysis will be performed for two main modes of sequence window:\n"
printf "length 13 (sequence6) and 15 (sequence7)\n"
site_mode=("sequence6" "sequence7")
for i in ${site_mode[@]}; do
            echo item: $i
            python ./scripts/motif.py $PF $i
            printf "Quick summary of results:\n"
            head -3 $PF/results/kinase_motif_scan_$i.tsv
            echo "End!"
        done
printf "If interested in the CK2 consensus motif the number of sequences containing\n"
printf "one or more CK2 motif in the vicinity of the identified phosphosites is:\n"
for i in ${site_mode[@]}; do
            echo item: $i
            cat "$PF/results/ck2_hit_$i.tsv" | wc -l
            echo "End!"
        done
printf "Done!\n"

printf "8. Retrieving sequence and functional information:\n"
printf "Analysis will be performed on all proteins found in the experiment:\n"
while read p; do
  python ./scripts/function_seq.py $PF $p
done < $PF/results/proteins.tsv
printf "Done!\n"
printf "Quick summary of functional information:\n"
head -7 $PF/results/swp_info.tsv

printf "End of analysis pipeline\n"
printf "Have a good day!\n"
