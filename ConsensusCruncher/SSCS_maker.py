#!/usr/bin/env python3

###############################################################
#
#      Single Stranded Consensus Sequence (SSCS) Generator
#
# Author: Nina Wang
# Date Created: Mar 24, 2016
###############################################################
# Function:
# To generate single strand consensus sequences for strand based error suppression.
# - Consensus sequence from most common base with quality score >= Q30 and greater than <cutoff> representation
# - Consensus quality score from addition of quality scores (i.e. product of error probabilities)
#
# Written for Python 3.5.1
#
# Usage:
# python3 SSCS_maker.py [--cutoff CUTOFF] [--infile INFILE] [--outfile OUTFILE] [--bedfile BEDFILE]
#
# Arguments:
# --cutoff CUTOFF     Proportion of nucleotides at a given position in a sequence required to be identical to form a
#                     consensus (Recommendation: 0.7 - based on previous literature Kennedy et al.)
#                        Example (--cutoff = 0.7):
#                           Four reads (readlength = 10) are as follows:
#                              Read 1: ACTGATACTT
#                              Read 2: ACTGAAACCT
#                              Read 3: ACTGATACCT
#                              Read 4: ACTGATACTT
#                           The resulting SSCS is: ACTGATACNT
# --infile INFILE     Input BAM file
# --outfile OUTFILE   Output BAM file
# --bedfile BEDFILE   Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt -
#                     See bed_separator.R for making your own bed file based on a target panel / specific coordinates)
#
# Inputs:
# 1. A position-sorted BAM file containing paired-end reads with duplex barcode in the header
# 2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing
#
# Outputs:
# 1. A SSCS BAM file containing paired single stranded consensus sequences - "sscs.bam"
# 2. A singleton BAM file containing single reads - "singleton.bam"
# 3. A bad read BAM file containing unpaired, unmapped, and multiple mapping reads - "badReads.bam"
# 4. A text file containing summary statistics (Total reads, Unmmaped reads, Secondary/Supplementary reads, SSCS reads,
#    and singletons) - "stats.txt"
# 5. A tag family size distribution plot (x-axis: family size, y-axis: number of reads) - "tag_fam_size.png"
# 6. A text file tracking the time to complete each genomic region (based on bed file) - "time_tracker.txt"
#
# Concepts:
#    - Read family: reads that share the same molecular barcode, genome
#                   coordinates for Read1 and Read2, cigar string, strand, flag, and read number
#    - Singleton: a read family containing only one member (a single read)
#
###############################################################

##############################
#        Load Modules        #
##############################
import pysam  # Need to install
import collections
import re
import array
from random import *
from itertools import chain
import argparse
import matplotlib.pyplot as plt
import math
import time

from consensus_helper import *


###############################
#       Helper Functions      #
###############################
def consensus_maker(readList, cutoff):
    """(list, int) -> str, list, list
    Return consensus sequence and quality score.

    Arguments:
        - readList: list of reads sharing the same unique molecular identifier
        - cutoff: Proportion of nucleotides at a given position in a sequence required to be identical to form a consensus

    Concept:
        Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position.
        - At each position, reads supporting each nucleotide is recorded along with the quality score corresponding to
          each nucleotide
        - Bases below the Phred quality cutoff (Q30) are excluded from consensus making
        - The most frequent base is added to the consensus sequence, given that the proportion of reads supporting this
          base is greater than the cutoff
        - A molecular phred quality score (consensus quality score) is determined by taking the product of errors of the
          most frequent base
        - If a majority can't be determined (i.e. a tie with 2 maximums), N will be assigned as these bases won't pass
          the proportion cut-off
    """
    # Initialize counters
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = []
    a = 0
    b = 0

    # Determine consensus for every position across read length
    readLength = readList[0].infer_query_length()     #看提取出来的这个片段的长度infer query length from CIGAR alignment. This method deduces the query length from the CIGAR alignment but does not include hard-clipped bases.

    for i in range(readLength):  #这组序列有多长
        # Positions in the following lists corresponds to A, C, G, T, N
        nuc_count = [0, 0, 0, 0, 0]     #统计保留的碱基中A, C, G, T, N出现的次数
        failed_nuc_count = [0, 0, 0, 0, 0]   #分别对应A, C, G, T, N，统计质量分数不到30的碱基
        quality_score = [[], [], [], []]
        phred_fail = 0

        # Count bases and quality scores for position i across all reads in list
        for j in range(len(readList)):    #这一组有多少条序列
            # Filter bases < phred quality 30 into separate list     #去掉质量分数小于30的
            if readList[j].query_qualities[i] < 30:   #read sequence base qualities, including soft clipped bases (None if not present).
                nuc = readList[j].query_sequence[i]   #质量分数小于30对应的碱基
                nuc_index = nuc_lst.index(nuc)    #index() 方法返回指定值首次出现的位置
                failed_nuc_count[nuc_index] += 1    #统计出去掉的碱基的个数
                phred_fail += 1
            else:
                nuc = readList[j].query_sequence[i]      #得到该位置的碱基是什么
                nuc_index = nuc_lst.index(nuc)
                nuc_count[nuc_index] += 1     #统计碱基的出现次数
                
                quality_score[nuc_index].append(readList[j].query_qualities[i])   #把碱基对应的质量值添加进对应的列表

        # Find most frequent nucleotide base and quality score (don't worry about ties (2 maxes) as it won't pass the
        # proportion cut-off and N will be assigned)
        max_nuc_index = nuc_count.index(max(nuc_count))   #选取A, C, G, T, N出现次数最多的，返回出现的位置，比如A对应0
        max_nuc_quality = quality_score[max_nuc_index]    #出现次数最多的碱基的所有质量分数的列表

        # Determine consensus phred quality through addition of quality scores (i.e. product of error probabilities)
        base_fail = False     #表示这个位置的碱基不存在，或者说为N （从下面判断出来）
        if max_nuc_quality is not []:
            mol_qual = sum(max_nuc_quality)    #计算所有的质量分数的和
            # Set to max quality score if sum of qualities is greater than the threshold (Q60) imposed by genomic tools
            if mol_qual > 60:
                mol_qual = 60
        else:
            mol_qual = 0
            base_fail = True

        # Consensus only made if proportion of most common base is > cutoff (e.g. 70%)
        phred_pass_reads = len(readList) - phred_fail  # Remove number of failed bases from total count 计算每行去掉质量分数不到30的碱基后还有多少个
        if phred_pass_reads != 0:
            prop_score = nuc_count[max_nuc_index]/phred_pass_reads   #出现次数最多的碱基个数除以总个数，计算该碱基出现的概率
            if prop_score >= cutoff:
                consensus_read += nuc_lst[max_nuc_index]   #如果出现次数最多的碱基概率超过0.7，就认为是一致性碱基
                quality_consensus.append(mol_qual)         #把对应的质量分数作为该碱基的质量分数
                a += 1
            else:
                base_fail = True
        else:
            base_fail = True

        # Set base to N if no consensus could be made
        if base_fail:
            consensus_read += 'N'
            quality_consensus.append(mol_qual)    #如果该位置的碱基是N，则对应的质量分数为0
            b =+ 1

    return consensus_read, quality_consensus, a, b


# Improve readability of argument help documentation
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


###############################
#        Main Function        #
###############################
def main():
    # Command-line parameters #import argparse 这是这个 argparse 模块中导入的，用于创建一个解析对象，这里是输入-h时的帮助行
    parser = ArgumentParser(formatter_class=SmartFormatter)
    parser.add_argument("--cutoff", action="store", dest="cutoff", type=float,
                        help="R|Proportion of nucleotides at a given position in a\nsequence required to be identical"
                        " to form a consensus\n(Recommendation: 0.7 - based on previous literature\nKennedy et al.)\n"
                        "   Example (--cutoff = 0.7):\n"
                        "       Four reads (readlength = 10) are as follows:\n"
                        "       Read 1: ACTGATACTT\n"
                        "       Read 2: ACTGAAACCT\n"
                        "       Read 3: ACTGATACCT\n"
                        "       Read 4: ACTGATACTT\n"
                        "   The resulting SSCS is: ACTGATACNT",
                        required=True)
    parser.add_argument("--infile", action="store", dest="infile", help="Input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="Output SSCS BAM file", required=True)
    parser.add_argument("--bdelim", action="store", dest="bdelim", default="|",
                        help="Delimiter to differentiate barcodes from read name, default: '|'")
    parser.add_argument("--bedfile", action="store", dest="bedfile",
                        help="Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt - \
                        See bed_separator.R for making your own bed file based on a target panel/specific coordinates)",
                        required=False)
    args = parser.parse_args()

    ######################
    #       SETUP        #
    ######################
    start_time = time.time()   #获取当前的时间
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")      #读取bam文件
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)    #输出一个bam文件
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    # set up time tracker
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')

    # ===== Initialize dictionaries =====
    read_dict = collections.OrderedDict()          #定义一个字典，该字典使用OrderedDict会根据放入元素的先后顺序进行排序。所以输出的值是排好序的
    tag_dict = collections.defaultdict(int)        #defaultdict可以被用来计数，统计每个字母出现的次数，具体可以参考流程文件夹中的文档函数注释
    pair_dict = collections.defaultdict(list)      #统计相同字符串出现的位置
    csn_pair_dict = collections.defaultdict(list)

    # ===== Initialize counters =====
    unmapped = 0
    multiple_mapping = 0  # secondary/supplementary reads
    counter = 0
    singletons = 0
    SSCS_reads = 0
    bad_spacer = 0
    base_number = 0
    base_N = 0
    del_sscs = 0

    #######################
    #   SPLIT BY REGION   #
    #######################
    # ===== Determine data division coordinates =====
    # division by bed file if provided
    if args.bedfile is not None:     #函数bed_separator在consensus_helper.py中
        division_coor = bed_separator(args.bedfile)     #这里返回的是一个字典，字典的键是hg19_cytoBand.txt文件中的第一列和第四列，格式为{}_{}，即chr1_p36.33
    else:                                               #对应的值是一个数组，包含起始位置和终点位置，即{0， 2300000}
        division_coor = [1]

    # ===== Process data in chunks =====
    region = 0
    for x in division_coor:
        if division_coor == [1]:
            read_chr = None
            read_start = None
            read_end = None
        else:
            read_chr = x.split('_', 1)[0]      #染色体的名称  eg：chr1
            read_start = division_coor[x][0]   #染色体特征起始位置 eg：0
            read_end = division_coor[x][1]     #染色体特征终止位置 eg：2300000

        # === Construct dictionaries for consensus making ===   函数read_bam在consensus_helper.py中
        chr_data = read_bam(bamfile,  #读取的bam文件
                            read_dict=read_dict,    #定义的一个有顺序的字典
                            tag_dict=tag_dict,      #定义的一个用于计数的字典
                            pair_dict=pair_dict,
                            csn_pair_dict=csn_pair_dict,
                            badRead_bam=badRead_bam,    #创建的一个bam文件
                            duplex=None,  # this indicates bamfile is not for making DCS (thus headers are diff)
                            read_chr=read_chr,    #染色体的名称  eg：chr1
                            read_start=read_start,   #染色体特征起始位置 eg：0
                            read_end=read_end,       #染色体特征终止位置 eg：2300000
                            barcode_delim=args.bdelim)
        #返回的内容：read_dict, tag_dict, pair_dict, csn_pair_dict, counter, unmapped_mate, multiple_mapping, bad_spacer

        # Set dicts and update counters
        read_dict = chr_data[0]  #OrderedDict([('GT.TC_22_4292_22_4336_123M_123M_fwd_R2', [<pysam.libcalignedsegment.AlignedSegment object at 0x7f1d732d99a0>]), ('GT.TC_22_4336_22_4292_123M_123M_rev_R1', [<pysam.libcalignedsegment.AlignedSegment object at 0x7f1d732d98e0>])])
        tag_dict = chr_data[1]  #defaultdict(<class 'int'>, {'GT.TC_22_4292_22_4336_123M_123M_fwd_R2': 1, 'GT.TC_22_4336_22_4292_123M_123M_rev_R1': 1})
        pair_dict = chr_data[2]
        csn_pair_dict = chr_data[3] #defaultdict(<class 'list'>, {'GT.TC_22_4292_22_4336_123M_123M_neg_167': ['GT.TC_22_4292_22_4336_123M_123M_fwd_R2', 'GT.TC_22_4336_22_4292_123M_123M_rev_R1']})

        counter += chr_data[4]
        unmapped += chr_data[5]
        multiple_mapping += chr_data[6]
        bad_spacer += chr_data[7]

        ######################
        #     CONSENSUS      #
        ######################
        # ===== Create consensus sequences for paired reads =====
        for readPair in list(csn_pair_dict.keys()):      #键值的列表，键表示consensus_tag 例如: GT.TC_22_4292_22_4336_123M_123M_neg_167，值表示两条序列的方向等信息
            if len(csn_pair_dict[readPair]) == 2:        #指的是一条链的R1和R2？
                for tag in csn_pair_dict[readPair]:
                    # Check for singletons
                    if tag_dict[tag] == 1:
                        singletons += 1
                        # Assign singletons our unique query name
                        read_dict[tag][0].query_name = readPair + ':' + str(tag_dict[tag])
                        singleton_bam.write(read_dict[tag][0])   #这一部分if的意思就是把正负链只含有一条的单例更改一个新的序列名称，并存储到singleton_bam文件中
                    else:
                        # Create collapsed SSCSs    这里的read_dict[tag]就是fetch提取出来的内容
                        SSCS = consensus_maker(read_dict[tag], float(args.cutoff))     #调用上面的函数
                        base_number = SSCS[2] + SSCS[3] + base_number
                        base_N += SSCS[3]

                        num_N = 0
                        for b in SSCS[0]:
                            if b == 'N':
                                num_N += 1
                        if num_N / len(SSCS[0]) <= 0.3:
                            query_name = readPair + ':' + str(tag_dict[tag])
                            SSCS_read = create_aligned_segment(read_dict[tag], SSCS[0], SSCS[1], query_name)   #create_aligned_segment在consensus_helper.py中
                            # Write consensus bam
                            SSCS_bam.write(SSCS_read)
                            SSCS_reads += 1
                        else:
                            del_sscs += 1

                    # Remove read from dictionary after writing
                    del read_dict[tag]

                # Remove key from dictionary after writing
                del csn_pair_dict[readPair]

        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')
        except:
            # When no genomic coordinates (x) provided for data division
            continue

    ######################
    #       SUMMARY      #
    ######################
    # === STATS ===
    # Note: total reads = unmapped + secondary + SSCS uncollapsed + singletons
    summary_stats = '''# === SSCS ===
Uncollapsed - Total reads: {}
Uncollapsed - Unmapped reads: {}
Uncollapsed - Secondary/Supplementary reads: {}
SSCS reads: {}
Singletons: {}
base number: {}
base N: {}
Bad spacers: {}
del sscs num: {}\n'''.format(counter, unmapped, multiple_mapping, SSCS_reads, singletons, base_number, base_N,  bad_spacer, del_sscs)

    stats.write(summary_stats)
    print(summary_stats)

    # === QC to see if there's remaining reads ===
    print('# QC: Total uncollapsed reads should be equivalent to mapped reads in bam file.')
    print('Total uncollapsed reads: {}'.format(counter))
    print('Total mapped reads in bam file: {}'.format(bamfile.mapped))

    print("QC: check dictionaries to see if there are any remaining reads")
    print('=== pair_dict remaining ===')
    if bool(pair_dict):
        for i in pair_dict:
            try:
                print(i)
                print('read remaining:')
                print(pair_dict[i][0])
                print('mate:')
                print(bamfile.mate(pair_dict[i][0]))
            except ValueError:
                print("Mate not found")
    print('=== read_dict remaining ===')
    if bool(read_dict):
        for i in read_dict:
            try:
                print(i)
                print('read remaining:')
                print(read_dict[i][0])
                print('mate:')
                print(bamfile.mate(read_dict[i][0]))
            except ValueError:
                print("Mate not found")
    print('=== csn_pair_dict remaining ===')
    if bool(csn_pair_dict):
        for i in csn_pair_dict:
            try:
                print(i)
                print(csn_pair_dict[i])
            except ValueError:
                print("Mate not found")

    # ===== write tag family size dictionary to file =====
    tags_per_fam = collections.Counter([i for i in tag_dict.values()])  # count of tags within each family size collections.Counter统计出现次数
    lst_tags_per_fam = list(tags_per_fam.items())  # convert to list [(fam, numTags)]
    with open(args.outfile.split('.sscs')[0] + '.read_families.txt', "w") as stat_file:
        stat_file.write('family_size\tfrequency\n')
        stat_file.write('\n'.join('%s\t%s' % x for x in lst_tags_per_fam))

    # ===== Create tag family size plot =====
    total_reads = sum(tag_dict.values())
    # Read fraction = family size * frequency of family / total reads
    read_fraction = [(i*j)/total_reads for i, j in lst_tags_per_fam]

    plt.bar(list(tags_per_fam), read_fraction)
    # Determine read family size range to standardize plot axis
    plt.xlim([0, math.ceil(lst_tags_per_fam[-1][0]/10) * 10])
    plt.savefig(args.outfile.split('.sscs')[0]+'_tag_fam_size.png')

    # ===== Close files =====
    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    badRead_bam.close()


###############################
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)

