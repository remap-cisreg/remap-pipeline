import argparse


parser = argparse.ArgumentParser()



parser.add_argument('-b', '--bed', action='store', help='Path to the bed file')

parser.add_argument('-o', '--output', action='store', help='Path to the output files')

parser.add_argument('-nb_peaks', '--nb_peaks', action='store', help='Path to the output files')

parser.add_argument('-len_peaks', '--len_peaks', action='store', help='Path to the output files')

args = parser.parse_args()

Error_input = argparse.ArgumentTypeError('input file necessary')

if args.nb_peaks is None:

    nb_peaks = 80 000

else:

    bed = int(args.nb_peaks)

if args.len_peaks is None:

    nb_peaks = 1500

else:

    bed = int(args.len_peaks)

if args.bed is None:

    bed = "remap2022.bed"

else:

    bed = args.bed





if args.output is None:

    out = "remap2022"

else:

    out = args.output


out1 = open(out+'_len_peak.bed', 'w')

out2 = open(out+'_filter_len_peak.bed', 'w')

dict_GSE = {}

with open(bed, 'r') as f:

    for line in f:

        pos1 = int(line.rstrip().split('\t')[1])

        pos2 = int(line.rstrip().split('\t')[2])

        if pos2 - pos1 >= lim_len_peaks:

            out1.write(line)

        else:

            out2.write(line)

        GSE = line.rstrip().split('\t')[3].split('_')[0]

        if GSE in dict_GSE.keys():

            dict_GSE[GSE] += 1

        else:

            dict_GSE[GSE] = 1



out1.close()

out2.close()



list_GSE = [line.rstrip().split('\t')[0] for line in open("GSE_peaks.txt","r")]
list_nb = [line.rstrip().split('\t')[1] for line in open("GSE_peaks.txt","r")]
dict_GSE = dict(zip(list_GSE, list_nb))

out_GSE = open('GSE_peaks.txt', 'w')
dict_GSE = {}
with open(bed, 'r') as f:

    for line in f:
        GSE = line.rstrip().split('\t')[3].split('_')[0]

        if GSE in dict_GSE.keys():

            dict_GSE[GSE] += 1

        else:

            dict_GSE[GSE] = 1

for GSE in dict_GSE.keys():
    out_GSE.write(GSE + "\t" + str(dict_GSE[GSE]) + "\n")
out_GSE.close()

out3 = open(out+'_nb_peak.bed', 'w')

out4 = open(out+'_filter_nb_peak.bed', 'w')

f = open(out+'_filter_len_peak.bed','r')
i = 0
for line in f:
    
    GSE = line.rstrip().split('\t')[3]
    i+=1
    print(i)
    nb_peak = int(dict_GSE[GSE])

    if nb_peak >= lim_nb_peaks:

       out3.write(line)

    else:

        out4.write(line)



out3.close()

out4.close()
