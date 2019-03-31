#!
import sys


#Read each record storing the first, if current record is missing a header than find it in the previous quality and split that line.

prev_header = ''
prev_seq = ''
prev_score_header = ''



fastq = sys.argv[1]
line_count = 1


for line in open(fastq):
    line = line.rstrip()
    if line_count/1 == 1:
        prev_header = line

    if line_count/2 == 1:
        prev_seq = line

    if line_count/3 == 1:
        prev_score_header = line

    if line_count/4 == 1:
        print (prev_header)
        print (prev_seq)
        print (prev_score_header)
        len_score = len(line)
        if len_score != len(prev_seq):
            line_count = 1
            header = '@' + line.rstrip().split('@')[-1]
            fixed_scores = '@'.join(line.rstrip().split('@')[:-1])
            prev_header = header
            print (fixed_scores)

        else:
            print (line)

    line_count += 1
    if line_count == 5:
        line_count = 1