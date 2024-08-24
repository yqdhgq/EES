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
def consensus_maker(readList):
    nuc_list = ['A', 'T', 'C', 'G', 'N']  # 碱基列表
    base_list = ['A', 'T', 'C', 'G']
    consensus_read = ''
    quality_consensus = []
    x = []
    a = 0
    b = 0

    readLength = readList[0].infer_query_length()

    # print (len(readList))
    # print(readLength)

    for i in range(readLength):  # 这条序列的长度
        nuc_count = [0, 0, 0, 0, 0]  # 统计同一位置上A,T,C,G,N 出现的次数
        quality_score = [[], [], [], []]  # 统计同一位置上A,T,C,G,N 出现对应的质量分数
        failed_nuc_count = [0, 0, 0, 0, 0]
        probability_list = [0, 0, 0, 0, 0]  # 分别对应真实碱基为A、T、C、G 时的概率
        true_pro_list = []

        phred_fail = 0

        for j in range(len(readList)):  # 这组序列有多少条，不同序列的同一位置
            nuc = readList[j].query_sequence[i]  # 碱基是什么
            qul = readList[j].query_qualities[i]  # 碱基对应的质量值
            if qul < 30:
                nuc_index = nuc_list.index(nuc)
                failed_nuc_count[nuc_index] += 1  # 统计出去掉的碱基的个数
                phred_fail += 1
            else:
                nuc_index = nuc_list.index(nuc)
                nuc_count[nuc_index] += 1  # 统计碱基的出现次数
                quality_score[nuc_index].append(readList[j].query_qualities[i])  # 把碱基对应的质量值添加进对应的列表

        # print(nuc_count)
        # print(quality_score)
        if phred_fail / len(readList) >= 0.7:
            quality_consensus.append(10)
            consensus_read += nuc_list[4]
            #x.append(1)
            a = + 1
        else:
            for base in base_list:  # 当真实碱基为base时
                probability = 1
                # print ("while base is : {}".format(base))
                for j in range(len(readList)):  # 这组序列有多少条，不同序列的同一位置
                    nuc = readList[j].query_sequence[i]  # 碱基是什么
                    qul = readList[j].query_qualities[i]  # 碱基对应的质量值
                    if nuc == base:
                        prior_probability = 1 - pow(10, -0.1 * qul)
                    else:
                        prior_probability = pow(10, -0.1 * qul) / 3.0

                    probability = probability * prior_probability
                nuc_index = nuc_list.index(base)
                probability_list[nuc_index] = probability  # 当前假设的碱基产生的概率，贝叶斯公式的分子部分，添加到probability_list对应的位置中
            # print('probability_list: {}'.format(probability_list) )

            sum_probability = sum(probability_list)  # 计算贝叶斯公式的分母部分
            # print('sum_probability: {}'.format(sum_probability))
            if sum_probability == 0:
                quality_consensus.append(10)
                consensus_read += nuc_list[4]
            else:
                for num in range(4):
                    true_probability = probability_list[num] / sum_probability  # 分别计算A,T,G,C为真实碱基的概率
                    true_pro_list.append(true_probability)
                # sprint(true_pro_list)

                true_index = true_pro_list.index(max(true_pro_list))

                if max(true_pro_list) == 1:
                    quality_consensus.append(60)
                    consensus_read += base_list[true_index]
                    #x.append(2)
                    b = + 1
                else:
                    true_quality = round(-10 * math.log(1 - max(true_pro_list), 10))
                    if true_quality > 60:
                        quality_consensus.append(60)
                    else:
                        quality_consensus.append(true_quality)
                    consensus_read += base_list[true_index]
                    #x.append(3)
                    b = + 1
                    # print(base_list[true_index])

    # print(x)
    #print('consensus_read : {}'.format(consensus_read))
    #print('quality_consensus : {}'.format(quality_consensus))

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
    start_time = time.time()  # 获取当前的时间
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")  # 读取bam文件
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template=bamfile)  # 输出一个bam文件
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb",
                                        template=bamfile)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template=bamfile)

    # set up time tracker
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')

    # ===== Initialize dictionaries =====
    read_dict = collections.OrderedDict()  # 定义一个字典，该字典使用OrderedDict会根据放入元素的先后顺序进行排序。所以输出的值是排好序的
    tag_dict = collections.defaultdict(int)  # defaultdict可以被用来计数，统计每个字母出现的次数，具体可以参考流程文件夹中的文档函数注释
    pair_dict = collections.defaultdict(list)  # 统计相同字符串出现的位置
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
    if args.bedfile is not None:  # 函数bed_separator在consensus_helper.py中
        division_coor = bed_separator(
            args.bedfile)  # 这里返回的是一个字典，字典的键是hg19_cytoBand.txt文件中的第一列和第四列，格式为{}_{}，即chr1_p36.33
    else:  # 对应的值是一个数组，包含起始位置和终点位置，即{0， 2300000}
        division_coor = [1]

    # ===== Process data in chunks =====
    region = 0
    for x in division_coor:
        if division_coor == [1]:
            read_chr = None
            read_start = None
            read_end = None
        else:
            read_chr = x.split('_', 1)[0]  # 染色体的名称  eg：chr1
            read_start = division_coor[x][0]  # 染色体特征起始位置 eg：0
            read_end = division_coor[x][1]  # 染色体特征终止位置 eg：2300000

        # === Construct dictionaries for consensus making ===   函数read_bam在consensus_helper.py中
        chr_data = read_bam(bamfile,  # 读取的bam文件
                            read_dict=read_dict,  # 定义的一个有顺序的字典
                            tag_dict=tag_dict,  # 定义的一个用于计数的字典
                            pair_dict=pair_dict,
                            csn_pair_dict=csn_pair_dict,
                            badRead_bam=badRead_bam,  # 创建的一个bam文件
                            duplex=None,  # this indicates bamfile is not for making DCS (thus headers are diff)
                            read_chr=read_chr,  # 染色体的名称  eg：chr1
                            read_start=read_start,  # 染色体特征起始位置 eg：0
                            read_end=read_end,  # 染色体特征终止位置 eg：2300000
                            barcode_delim=args.bdelim)
        # 返回的内容：read_dict, tag_dict, pair_dict, csn_pair_dict, counter, unmapped_mate, multiple_mapping, bad_spacer

        # Set dicts and update counters
        read_dict = chr_data[
            0]  # OrderedDict([('GT.TC_22_4292_22_4336_123M_123M_fwd_R2', [<pysam.libcalignedsegment.AlignedSegment object at 0x7f1d732d99a0>]), ('GT.TC_22_4336_22_4292_123M_123M_rev_R1', [<pysam.libcalignedsegment.AlignedSegment object at 0x7f1d732d98e0>])])
        tag_dict = chr_data[
            1]  # defaultdict(<class 'int'>, {'GT.TC_22_4292_22_4336_123M_123M_fwd_R2': 1, 'GT.TC_22_4336_22_4292_123M_123M_rev_R1': 1})
        pair_dict = chr_data[2]
        csn_pair_dict = chr_data[
            3]  # defaultdict(<class 'list'>, {'GT.TC_22_4292_22_4336_123M_123M_neg_167': ['GT.TC_22_4292_22_4336_123M_123M_fwd_R2', 'GT.TC_22_4336_22_4292_123M_123M_rev_R1']})

        counter += chr_data[4]
        unmapped += chr_data[5]
        multiple_mapping += chr_data[6]
        bad_spacer += chr_data[7]

        ######################
        #     CONSENSUS      #
        ######################
        # ===== Create consensus sequences for paired reads =====
        for readPair in list(csn_pair_dict.keys()):  # 键值的列表，键表示consensus_tag 例如: GT.TC_22_4292_22_4336_123M_123M_neg_167，值表示两条序列的方向等信息
            if len(csn_pair_dict[readPair]) == 2:  # 比如在正链pos中同时包含R1和R2
                for tag in csn_pair_dict[readPair]:
                    # Check for singletons
                    if tag_dict[tag] == 1:
                        singletons += 1
                        # Assign singletons our unique query name
                        read_dict[tag][0].query_name = readPair + ':' + str(tag_dict[tag])
                        singleton_bam.write(read_dict[tag][0])  # 这一部分if的意思就是把正负链只含有一条的单例更改一个新的序列名称，并存储到singleton_bam文件中
                    else:
                        # Create collapsed SSCSs    这里的read_dict[tag]就是fetch提取出来的内容
                        #print('tag : {}'.format(tag))
                        #print('read_dict[tag] : {}'.format(read_dict[tag]))

                        SSCS = consensus_maker(read_dict[tag])  # 调用上面的函数
                        base_number = base_number + SSCS[2] + SSCS[3]
                        base_N += SSCS[2]

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
            time_tracker.write(str((time.time() - start_time) / 60) + '\n')
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
base number: {}
base num N :{}
SSCS reads: {}
Singletons: {}
Bad spacers: {}
del sscs num: {}\n'''.format(counter, unmapped, multiple_mapping, base_number, base_N, SSCS_reads, singletons, bad_spacer, del_sscs)

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
    tags_per_fam = collections.Counter(
        [i for i in tag_dict.values()])  # count of tags within each family size collections.Counter统计出现次数
    lst_tags_per_fam = list(tags_per_fam.items())  # convert to list [(fam, numTags)]
    with open(args.outfile.split('.sscs')[0] + '.read_families.txt', "w") as stat_file:
        stat_file.write('family_size\tfrequency\n')
        stat_file.write('\n'.join('%s\t%s' % x for x in lst_tags_per_fam))

    # ===== Create tag family size plot =====
    total_reads = sum(tag_dict.values())
    # Read fraction = family size * frequency of family / total reads
    read_fraction = [(i * j) / total_reads for i, j in lst_tags_per_fam]

    plt.bar(list(tags_per_fam), read_fraction)
    # Determine read family size range to standardize plot axis
    plt.xlim([0, math.ceil(lst_tags_per_fam[-1][0] / 10) * 10])
    plt.savefig(args.outfile.split('.sscs')[0] + '_tag_fam_size.png')

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
    print((time.time() - start_time) / 60)