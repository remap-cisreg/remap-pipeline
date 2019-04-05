#!/bin/sh 
# Count the nuber of files in each directory
# then when a directory does contain only 1 file
# it find a *fq.gz file
# print it, and remove it

find . -maxdepth 1 -type d -print0 | while read -d '' -r dir; do 

num=$(find $dir -name *.fq.gz | wc -l); 
#printf "%5d files in directory %s\n" "$num" "$dir";
if [ $num -eq 1 ]
then 
printf "%5d files in dir: %s\n" "$num" "$dir"; 
find $dir -name *.fq.gz -print
#find $dir -name *.fq.gz -exec rm {} \; 
echo
fi

done
