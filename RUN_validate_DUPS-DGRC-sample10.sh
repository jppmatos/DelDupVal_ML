#!/bin/bash

input_TP_4=./DGRC_sampled/DGRC_DELS_raw_5_09_Sample10.csv

mkdir ./DGRC_sampled/output/
mkdir ./DGRC_sampled/output/deletions/
output_path_TP_4=./DGRC_sampled/output/deletions/ 
# ________________________________________________________________________________________________
echo "Processing:"

# look scripts:
python look_repRegions2.py $input_TP_4 $output_path_TP_4 &
python look_segDuplications2.py $input_TP_4 $output_path_TP_4 &

# Segmental Duplcations and Repititve regions covereage:
python segDuplication_coverage.py $input_TP_4 $output_path_TP_4 &

python repRegions_coverage.py $input_TP_4 $output_path_TP_4 &
python segDUPS_find_pair2.py $input_TP_4 $output_path_TP_4 &

# Centromeres and Telomeres:
python Centro_Telo_match2.py $input_TP_4 $output_path_TP_4 &

# GCcontent:
python GCcontent2.py $input_TP_4 $output_path_TP_4 &

# Sequence Complexity:
python SeqComplexity.py $input_TP_4 $output_path_TP_4 &

# Replication Time
#python Replication_Time2.py $input_TP_4 $output_path_TP_4 &

# Recombination Rate: #(WARRING - RISK OF MISSING DATA!)
python Recombination_Rate2.py $input_TP_4 $output_path_TP_4 &

#Introns and Exons coverage:
python Genes_V3.py $input_TP_4 $output_path_TP_4 &

# Lamina Associated domains:
python Lamina_Associ_Dom.py $input_TP_4 $output_path_TP_4 &

# CpG island distance:
python CpG_island2.py $input_TP_4 $output_path_TP_4 &

# TAD Boundaries:
python TADs2.py $input_TP_4 $output_path_TP_4 &


wait
echo "Dups data ready!!!"  
echo "Assembling dataset!"
python makedataset-deletions.py $output_path_TP_4 $output_path_TP_4 #"./output/duplications/" "./output/duplication/duplication_"
echo "Start predicting!"
python predictClassify.py $output_path_TP_4 DEL $output_path_TP_4
echo "Done!"