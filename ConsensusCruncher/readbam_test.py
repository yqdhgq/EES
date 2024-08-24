import pysam
import collections
import math
from random import randint




def bed_separator(bedfile):
    """(str) -> dict
    Return dictionary of coordinates based on bed file.
    """
    coor = collections.OrderedDict()

    with open(bedfile) as f:
        for line in f:
            chr_arm = line.split('\t')
            chr_key = '{}_{}'.format(chr_arm[0], chr_arm[3])
            start = int(chr_arm[1])
            end = int(chr_arm[2])
            chr_val = (start, end)

            coor[chr_key] = chr_val

    return coor

def which_strand(read):   #确定链的方形
    """

    Note: Strand is needed
    for the common identifier to replace read orientation and number (see sscs_qname for example)
    flag_pairings = {
        # paired and mapped
        99: 147, 147: 99, 83: 163, 163: 83,
        # mapped within insert size, but wrong orientation (++, --)
        67: 131, 131: 67, 115: 179, 179: 115,
        # mapped uniquely, but wrong insert size
        81: 161, 161: 81, 97: 145, 145: 97,
        # wrong insert size and wrong orientation
        65: 129, 129: 65, 113: 177, 177: 113
    }
    """
    # Flags indicating strand direction
    pos = [99, 147, 67, 131]    #正链
    neg = [83, 163, 115, 179]    #负链
    no_ori = [65, 129, 113, 177, 81, 161, 97, 145]  # direction not defined

    if read.flag in pos:
        strand = 'pos'
    elif read.flag in neg:
        strand = 'neg'
    elif read.flag in no_ori:
        # Determine orientation of flags with no defined direction using order of chr coor
        if (read.reference_id < read.next_reference_id and which_read(read.flag) == 'R1') or \
                (read.reference_id > read.next_reference_id and which_read(read.flag) == 'R2') or \
                (read.reference_id == read.next_reference_id and which_read(read.flag) == 'R1' and
                 read.reference_start < read.next_reference_start) or \
                (read.reference_id == read.next_reference_id and which_read(read.flag) == 'R2' and
                 read.reference_start > read.next_reference_start):
            strand = 'pos'
        else:
            strand = 'neg'
    else:
        # Only uniquely mapped reads (with flags indicated above) should be retained, as 'bad reads' were filtered out
        # in a previous step
        print('STRAND ERROR')
        print(read.flag)
        strand = None

    return strand

def which_read(flag):   #确定是R1还是R2
    """(int) -> str
    Returns read number based on flag.

    Test cases:
    >>> which_read(83)
    'R1'
    >>> which_read(131)
    'R2'
    >>> which_read(177)
    'R2'
    """
    read1 = [99, 83, 67, 115, 81, 97, 65, 113]
    read2 = [147, 163, 131, 179, 161, 145, 129, 177]

    if flag in read1:
        read = 'R1'
    elif flag in read2:
        read = 'R2'
    else:
        print('UNMAPPED READ ERROR')
        print(flag)
        read = None

    return read


def cigar_order(read, mate):
    ori_strand = which_strand(read)  #Return DNA strand of origin based on flags
    #print(ori_strand)
    read_num = which_read(read.flag)   #Returns read number based on flag.
    #print(read_num)

    if (ori_strand == 'pos' and read_num == 'R1') or (ori_strand == 'neg' and read_num == 'R2'):
        cigar = '{}_{}'.format(read.cigarstring,
                               mate.cigarstring)
    else:
        cigar = '{}_{}'.format(mate.cigarstring,
                               read.cigarstring)

    return cigar

def sscs_qname(read, mate, barcode, cigar):  #Return new query name for consensus sequences
    read_chr = read.reference_id
    #print(read_chr)
    mate_chr = mate.reference_id
    #print(mate_chr)
    read_coor = read.reference_start
   # print(read_coor)
    mate_coor = mate.reference_start
   #print(mate_coor )

    if (read_chr == mate_chr and int(read_coor) > int(mate_coor)) or \
            (int(read_chr) > int(mate_chr)):
        read_chr = mate.reference_id
        mate_chr = read.reference_id
        read_coor = mate.reference_start
        mate_coor = read.reference_start

    strand = which_strand(read)
    query_tag = '{}_{}_{}_{}_{}_{}_{}_{}'.format(barcode,
                                                 read_chr,
                                                 read_coor,
                                                 mate_chr,
                                                 mate_coor,
                                                 cigar,
                                                 strand,
                                                 abs(read.template_length))

    return query_tag


def unique_tag(read, barcode, cigar):
    """
    Return unique identifier tag for one read of a strand of a molecule.

    Tag uses following characteristics to group reads belonging to the same strand of an individual molecule (PCR dupes):
    [Barcode]_[Read Chr]_[Read Start]_[Mate Chr]_[Mate Start]_[Cigar String]_[Orientation]_[ReadNum]
    e.g. TTTG_24_58847416_24_58847448_137M10S_147M_fwd_R1
    """
    orientation = 'fwd'
    if read.is_reverse:
        orientation = 'rev'

    readNum = which_read(read.flag)

    # Unique identifier for strand of individual molecules
    tag = '{}_{}_{}_{}_{}_{}_{}_{}'.format(barcode,  # mol barcode
                                           read.reference_id,  # chr
                                           read.reference_start,  # start (0-based)
                                           read.next_reference_id,  # mate chr
                                           read.next_reference_start,  # mate start
                                           cigar,
                                           orientation,  # strand direction
                                           readNum
                                           )
    return tag

bamfile = pysam.AlignmentFile('LargeMid_58_L005.sorted.bam', "rb")
SSCS_bam = pysam.AlignmentFile('LargeMid_58_L005.sorted_test.bam', "wb", template = bamfile)
badRead_bam = pysam.AlignmentFile('LargeMid_58_L005.sorted.badReads.bam', "wb", template = bamfile)

division_coor = bed_separator('hg19_cytoBand.txt')

read_dict = collections.OrderedDict()          #定义一个字典，该字典使用OrderedDict会根据放入元素的先后顺序进行排序。所以输出的值是排好序的
tag_dict = collections.defaultdict(int)        #defaultdict可以被用来计数，统计每个字母出现的次数，具体可以参考流程文件夹中的文档函数注释
pair_dict = collections.defaultdict(list)      #统计相同字符串出现的位置
csn_pair_dict = collections.defaultdict(list)

i = 0

# ===== Initialize counters =====
unmapped = 0
multiple_mapping = 0  # secondary/supplementary reads
counter = 0
singletons = 0
SSCS_reads = 0
bad_spacer = 0

for x in division_coor:
    #print(x)
    if division_coor == [1]:
        read_chr = None
        read_start = None
        read_end = None
    else:
        read_chr = x.split('_', 1)[0]
        #print(read_chr)
        read_start = division_coor[x][0]
        #print(read_start)
        read_end = division_coor[x][1]
        #print(read_end)

    #a = bamfile.count(read_chr, read_start, read_end)
    #print(a)

    unmapped = 0
    unmapped_mate = 0
    multiple_mapping = 0  # secondary/supplementary reads
    counter = 0
    bad_spacer = 0
    duplex = None

    bamLines = bamfile.fetch(read_chr, read_start, read_end)
    for line in bamLines:
        #print(line)
        if read_chr is not None:
            #print(read_start)
            #print(read_end)
            #print(line.reference_start)
            if line.reference_start < read_start or line.reference_start > read_end:
                continue
            #else:
                #print('yes')

        counter += 1
        #print(counter)

        mate_unmapped = [73, 89, 121, 153, 185, 137]
        badRead = True
        barcode_delim = '|'
        """
        print(line.qname)
        print(line.is_unmapped)
        print(line.flag)
        print(line.is_secondary)
        print(line.is_supplementary)
        print(i)
        """

        # Check if delimiter is found in read
        if barcode_delim is not None and barcode_delim not in line.qname:   #barcode_delim = args.bdelim 在ConsensusCruncher.py中大多时候 = ’|‘
            bad_spacer += 1
        elif line.is_unmapped:
            unmapped += 1
            counter -= 1
        elif line.flag in mate_unmapped:
            unmapped_mate += 1
        elif line.is_secondary:
            multiple_mapping += 1
        elif line.is_supplementary:
            multiple_mapping += 1
        else:
            badRead = False

        # Write bad reads to file
        if badRead and badRead_bam is not None:
            badRead_bam.write(line)
        else:
            pair_dict[line.qname].append(line)
            """
            print(pair_dict)
            print(pair_dict[line.qname])
            print(len(pair_dict[line.qname]))
"""
            # === 2) ASSIGN UNIQUE IDENTIFIER TO READ PAIRS ===
            if len(pair_dict[line.qname]) == 2:
                #print(pair_dict[line.qname])
                read = pair_dict[line.qname][0]
                #print(read)
                mate = pair_dict[line.qname][1]
                #print(mate)
                # === Create consensus identifier ===
                # Extract molecular barcode, barcodes in diff position for SSCS vs DCS generation
                if duplex == None or duplex == False:  #this indicates bamfile is not for making DCS (thus headers are diff)
                    if barcode_delim is None:   #指定分割符，这里默认为’|‘
                        # SSCS query name: H1080:278:C8RE3ACXX:6:1308:18882:18072|CACT
                        barcode = read.qname.split("|")[1]
                        #print(barcode)
                    else:
                        barcode = read.qname.split(barcode_delim)[1]
                        #print(barcode)
                else:    #这里暂时不考虑，这是生成DCS
                    # DCS query name: CCTG_12_25398000_12_25398118_neg:5
                    barcode = read.qname.split("_")[0]

                # Consensus_tag cigar (ordered by strand and read)
                cigar = cigar_order(read, mate)
                #print(cigar)
                # Assign consensus tag as new query name for paired consensus reads
                consensus_tag = sscs_qname(read, mate, barcode, cigar)
                #print('consensus_tag : ' + consensus_tag)

                for i in range(2):
                    read_i = pair_dict[line.qname][i]
                    #print(read_i)
                    # Molecular identifier for grouping reads belonging to the same read of a strand of a molecule
                    tag = unique_tag(read_i, barcode, cigar)
                    #print(tag)

                    ######################
                    #   Assign to Dict   #
                    ######################
                    # === 3) ADD READ PAIRS TO DICTIONARIES ===
                    if tag not in read_dict and tag not in tag_dict:
                        read_dict[tag] = [read_i]
                        #print(read_dict)
                        tag_dict[tag] += 1
                        #print(tag_dict)

                        # Group paired unique tags using consensus tag
                        if consensus_tag not in csn_pair_dict:
                            csn_pair_dict[consensus_tag] = [tag]
                            #print(csn_pair_dict)

                        elif len(csn_pair_dict[consensus_tag]) == 2:
                            # Honestly this shouldn't happen anymore with these identifiers
                            print("Consensus tag NOT UNIQUE -> multiple tags (4) share same consensus tag [due to poor "
                                  "strand differentiation as a result of identifiers lacking complexity]")
                            """
                            print(consensus_tag)
                            print(tag)
                            print(read_i)
                            print(csn_pair_dict[consensus_tag])
                            print(read_dict[csn_pair_dict[consensus_tag][0]][0])
                            print(read_dict[csn_pair_dict[consensus_tag][1]][0])
                            """
                            # Manual inspection should be done on these reads
                        else:
                            csn_pair_dict[consensus_tag].append(tag)
                            #print(csn_pair_dict)

                    elif tag in tag_dict and read not in read_dict[tag]:
                        # Append reads sharing the same unique tag together (PCR dupes)
                        read_dict[tag].append(read_i)
                        tag_dict[tag] += 1
                    else:
                        # Data fetch error - line read twice (if its found in tag_dict and read_dict)
                        print('Pair already written: line read twice - check to see if read overlapping / near cytoband'
                              ' region (point of data division)')

                # remove read pair qname from pair_dict once reads added to read_dict
                pair_dict.pop(line.qname)


    """
    i = i + 1
    if i > 7:
        break
    """

def consensus_maker(readList, cutoff):
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = []

    # Determine consensus for every position across read length
    #print(readList[0])  #这里就是fetch提取出来的相关内容，readList表示的是两条及以上条序列内容
    readLength = readList[0].infer_query_length()  #看提取出来的这个片段的长度
    #print(readLength)
    for i in range(readLength):  #片段长度
        # Positions in the following lists corresponds to A, C, G, T, N
        nuc_count = [0, 0, 0, 0, 0]  #统计保留的碱基中A, C, G, T, N出现的次数
        failed_nuc_count = [0, 0, 0, 0, 0]  #分别对应A, C, G, T, N，统计质量分数不到30的碱基
        quality_score = [[], [], [], []]   #分别对应A，T，C，G
        phred_fail = 0
        #print(readList)
        #print(len(readList))

        # Count bases and quality scores for position i across all reads in list
        for j in range(len(readList)):  #表示有几条序列，不同序列的同一个位置
            # Filter bases < phred quality 30 into separate list     #去掉质量分数小于30的
            if readList[j].query_qualities[i] < 30:  # read sequence base qualities, including soft clipped bases (None if not present).
                #print(readList[j])
                #print(i)
                #print(readList[j].query_qualities[i])
                nuc = readList[j].query_sequence[i]  #质量小于30的碱基
                #print(nuc)
                nuc_index = nuc_lst.index(nuc)  #在nuc_list中的位置，比如如果是T，则对应3
                #print(nuc_index)
                failed_nuc_count[nuc_index] += 1  #在failed_nuc_count中对应的位置上加1，作为统计，比如说是一个T，则变为[0,0,0,1,0]
                #print(failed_nuc_count)
                phred_fail += 1
                #print(phred_fail)

            else:
                # print(readList[j].query_qualities[i])
                nuc = readList[j].query_sequence[i]
                # print(nuc)
                nuc_index = nuc_lst.index(nuc)
                # print(nuc_index)
                nuc_count[nuc_index] += 1
                # print(nuc_count[nuc_index])
                quality_score[nuc_index].append(readList[j].query_qualities[i])   #把保留的这个碱基所对应的质量值保存进quality_score中对应的列表里

        # Find most frequent nucleotide base and quality score (don't worry about ties (2 maxes) as it won't pass the
        # proportion cut-off and N will be assigned)
        max_nuc_index = nuc_count.index(max(nuc_count))   #选取A, C, G, T, N出现次数最多的，返回出现的位置，比如A对应0
        max_nuc_quality = quality_score[max_nuc_index]    #出现次数最多的碱基的所有质量分数的列表,如果是N则返回一个空列表
        #print(max_nuc_quality)

        # Determine consensus phred quality through addition of quality scores (i.e. product of error probabilities)
        base_fail = False     #表示这个位置的碱基不存在，或者说为N
        if max_nuc_quality is not []:
            mol_qual = sum(max_nuc_quality)  #计算所有的质量分数的和
            #print(mol_qual)
            # Set to max quality score if sum of qualities is greater than the threshold (Q60) imposed by genomic tools
            if mol_qual > 60:
                mol_qual = 60   #这里不知道为啥要等于60
        else:
            mol_qual = 0
            base_fail = True

        # Consensus only made if proportion of most common base is > cutoff (e.g. 70%)
        phred_pass_reads = len(readList) - phred_fail  # Remove number of failed bases from total count 计算每行去掉质量分数不到30的碱基后还有多少个
        if phred_pass_reads != 0:
            prop_score = nuc_count[max_nuc_index]/phred_pass_reads   #出现次数最多的碱基个数除以总个数，计算该碱基出现的概率
            if prop_score >= cutoff:
                consensus_read += nuc_lst[max_nuc_index]   #如果出现次数最多的碱基概率超过0.7，就认为是一致性碱基
                quality_consensus.append(mol_qual)         #把对应的质量分数作为该碱基的质量分数，一般为60
            else:
                base_fail = True  #表示这个位置为N
        else:
            base_fail = True

        # Set base to N if no consensus could be made
        if base_fail:
            consensus_read += 'N'
            quality_consensus.append(mol_qual)    #如果该位置的碱基是N，则对应的质量分数为0

    #print(consensus_read)
    #print(quality_consensus)

    return consensus_read, quality_consensus


def read_mode(field, bam_reads):
    field = 'i.{}'.format(field)
    # Rank by number of occurrences
    field_lst = collections.Counter(eval(field) for i in bam_reads).most_common()
    # Take max occurrences
    common_field_lst = [i for i, j in field_lst if j == field_lst[0][1]]
    # Randomly select max if there's multiple
    common_field = common_field_lst[randint(0, len(common_field_lst)-1)]

    return common_field

def consensus_flag(bam_reads):
    # Rank flags by number of occurrences
    count_flags = collections.Counter(i.flag for i in bam_reads).most_common()  # [(97, 1), (99, 1)]
    # List all flags with max count (will show multiple if there's a tie for the max count)
    max_flag = [i for i, j in count_flags if j == count_flags[0][1]]

    if len(max_flag) != 1:
        if 99 in max_flag:
            flag = 99
        elif 83 in max_flag:
            flag = 83
        elif 147 in max_flag:
            flag = 147
        elif 163 in max_flag:
            flag = 163
        else:
            flag = max_flag[randint(0, len(max_flag)-1)]  # If flag not properly paired/mapped, randomly select from max
    else:
        flag = max_flag[0]

    return flag


def create_aligned_segment(bam_reads, sscs, sscs_qual, query_name):
    # Use first read in list as template (all reads should share same cigar, template length, and coor)
    template_read = bam_reads[0]

    # Create consensus read based on template read
    SSCS_read = pysam.AlignedSegment()
    SSCS_read.query_name = query_name
    SSCS_read.query_sequence = sscs
    SSCS_read.reference_id = template_read.reference_id
    SSCS_read.reference_start = template_read.reference_start
    SSCS_read.mapping_quality = read_mode('mapping_quality', bam_reads)  # Most common mapping quality
    SSCS_read.cigar = template_read.cigar
    SSCS_read.next_reference_id = template_read.next_reference_id
    SSCS_read.next_reference_start = template_read.next_reference_start
    SSCS_read.template_length = read_mode('template_length', bam_reads)
    SSCS_read.query_qualities = sscs_qual

    # Most common flag used unless there's a tie, then flags are ranked if its 99/83/147/163, otherwise randomly picked
    SSCS_read.flag = consensus_flag(bam_reads)

    # Optional fields
    try:
        SSCS_read.set_tag('RG', read_mode("get_tag('RG')", bam_reads))
    except:
        pass

    return SSCS_read



for readPair in list(csn_pair_dict.keys()):
    # 键值的列表，键表示consensus_tag 例如: GT.TC_22_4292_22_4336_123M_123M_neg_167，值表示两条序列的方向等信息
    #print(readPair)
    #print(csn_pair_dict[readPair])
    if len(csn_pair_dict[readPair]) == 2:
        for tag in csn_pair_dict[readPair]:
            # Check for singletons
            #print(tag)
            if tag_dict[tag] == 1:
                #print(tag_dict[tag])
                singletons += 1
                # Assign singletons our unique query name
                #print(read_dict[tag])
                #print(read_dict[tag][0])
                #print(read_dict[tag][0].query_name)
                read_dict[tag][0].query_name = readPair + ':' + str(tag_dict[tag])
                #print(read_dict[tag][0].query_name)
                #singleton_bam.write(read_dict[tag][0])
            else:
                # Create collapsed SSCSs
                #print(tag)
                #print(read_dict[tag])
                SSCS = consensus_maker(read_dict[tag], 0.7)  # 调用上面的函数

                query_name = readPair + ':' + str(tag_dict[tag])
                #print(query_name)

                SSCS_read = create_aligned_segment(read_dict[tag], SSCS[0], SSCS[1], query_name)  # 调用上面函数

                # Write consensus bam
                SSCS_bam.write(SSCS_read)
                SSCS_reads += 1




    

"""
print('read_dict:')
print(read_dict)
print('tag_dict:')
print(tag_dict)
print('pair_dict:')
print(pair_dict)
print('csn_pair_dict:')
print(csn_pair_dict)
print('counter:')
print(counter)
print('unmapped_mate:')
print(unmapped_mate)
print('multiple_mapping:')
print(multiple_mapping)
print('bad_spacer:')
print(bad_spacer)
"""