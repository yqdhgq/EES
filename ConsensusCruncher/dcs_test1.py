import pysam  # Need to install
import collections
import re
import array
from random import randint
from argparse import ArgumentParser
import math

from consensus_helper import *


###############################
#        Main Function        #
###############################

######################
#       SETUP        #
######################
# start_time = time.time()
# ===== Initialize input and output bam files =====
sscs_bam = pysam.AlignmentFile("/mnt/nfs/chenhm/project/all/HD701/fs_2_bayes/consensus/SRR7903524.fastq.sorted/sscs/SRR7903524.fastq.sorted.sscs.sorted.bam", "rb")

# ===== Initialize dictionaries and counters=====
read_dict = collections.OrderedDict()
tag_dict = collections.defaultdict(int)
pair_dict = collections.defaultdict(list)
csn_pair_dict = collections.defaultdict(list)

unmapped = 0
unmapped_mate = 0
multiple_mapping = 0  # Secondary/supplementary reads
counter = 0
sscs_singletons = 0  # Single strand consensus sequences without a complementary strand
multiple_mappings = 0
base_number = 0
base_N = 0

duplex_count = 0
duplex_dict = collections.defaultdict(int)

#######################
#   SPLIT BY REGION   #
#######################
# ===== Determine data division coordinates =====
# division by bed file if provided
division_coor = bed_separator("/mnt/nfs/chenhm/project/ConsensusCruncher-master/ConsensusCruncher/hg19_cytoBand.txt")


# ===== Process data in chunks =====
for x in division_coor:
    if division_coor == [1]:
        read_chr = None
        read_start = None
        read_end = None
    else:
        read_chr = x.split('_', 1)[0]
        read_start = division_coor[x][0]
        read_end = division_coor[x][1]

    chr_data = read_bam(sscs_bam,
                        pair_dict=pair_dict,
                        read_dict=read_dict,
                        csn_pair_dict=csn_pair_dict,
                        tag_dict=tag_dict,
                        badRead_bam=None,
                        duplex=True,
                        read_chr=read_chr,
                        read_start=read_start,
                        read_end=read_end
                        )

    read_dict = chr_data[0]
    tag_dict = chr_data[1]
    pair_dict = chr_data[2]
    csn_pair_dict = chr_data[3]

    counter += chr_data[4]
    unmapped += chr_data[5]
    multiple_mapping += chr_data[6]

for i in list(read_dict.keys()):
    print("read_dict:")
    print(i)
    print(read_dict[i])
    break
for i in list(tag_dict.keys()):
    print("tag_dict:")
    print(i)
    print(tag_dict[i])
    break
for i in list(pair_dict.keys()):
    print("pair_dict:")
    print(i)
    print(pair_dict[i])
    break
for i in list(csn_pair_dict.keys()):
    print("csn_pair_dict:")
    print(i)
    print(csn_pair_dict[i])
    break