#!/bin/bash

PATH_REMAP="/scratch/bballester/thaliana/remap-thaliana"
PATH_COLOR_FILE="scratch/bballester/thaliana/remap-plants/analyses/A.thaliana_RGB.tsv"
PATH_LIST_EXP="/scratch/bballester/thaliana/remap-thaliana/list_exp.txt"


while read EXPERIEMENT; do

    #path pf narrowPeak
	EXPERIMENT_PATH=$(echo $PATH_REMAP"/6.peakcalling/"$EXPERIEMENT"/macs2/"$EXPERIEMENT"_peaks.narrowPeak")
    # Get TF name
	TF=$(echo $EXPERIEMENT | cut -d"." -f2)

    #echo $TF

    # Getting right color from color file (tf_name  R   G   B)
	while read COLORS; do
  		CURRENT_TF=$(echo "$COLORS" | awk -F"\t" '{print $1}')

  		if [ "$TF" == "$CURRENT_TF" ]; then
  			#COLOR_TF=$(echo "$COLORS" | awk -F"\t" '{print $2","$3","$4}')
        COLOR_TF=$(echo "$COLORS" | awk -F"\t" '{print $2}')
  		fi

	done <$PATH_COLOR_FILE

	#echo $COLOR_TF

    #Printing output with peak position
    awk -v color=$COLOR_TF -F"\t" '{sum=$10+$2;  if(sum<=$3) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t"sum"\t"sum+1"\t"color; else print $1"\t"$2"\t"$3"\t"$4"\t0\t.\t"sum-1"\t"sum"\t"color}' $EXPERIMENT_PATH


done <$PATH_LIST_EXP

