import argparse


parser = argparse.ArgumentParser()



parser.add_argument('-b', '--bed', action='store', help='Path to the bed file')




args = parser.parse_args()

Error_input = argparse.ArgumentTypeError('input file necessary')



if args.bed is None:

    bed = "remap2022_all_macs2_hg38_v1_0.bed"

else:

    bed = args.bed








dict_GSE = {}
with open(bed, 'r') as f:

    for line in f:
        GSE = line.rstrip().split('\t')[3]

        if GSE in dict_GSE.keys():
            dict_GSE[GSE].write(line.rstrip()+"\n")
            

        else:

            dict_GSE[GSE] = open("9.bed/datasets/"+GSE+".bed", "w")
            dict_GSE[GSE].write(line.rstrip()+"\n")
