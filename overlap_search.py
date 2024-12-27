# coding=utf-8
# pzw
# 20241227

import pandas as pd

# 读取BED文件
def read_bed(file_path):
    return pd.read_csv(file_path, sep='\t', header=0)

# 检查两个区间是否有交集
def has_overlap(row1, row2):
    if row1['end'] >= row2['start']:
        return True
    return False

# 寻找存在交集的行
def find_overlapping_rows(bed_df):
    overlapping_pairs = []
    for i in range(len(bed_df)):
        j = i + 1
        if j < len(bed_df):
            if bed_df.iloc[i]['#chrom'] == bed_df.iloc[j]['#chrom'] and has_overlap(bed_df.iloc[i], bed_df.iloc[j]):
                overlapping_pairs.append((i, j))
    return overlapping_pairs

# bed文件修正，如果chrom、start、end完全一致，此时保留refseq转录本编号值靠前的
def deduplicate_bed(df):
    """
    对BED文件进行去重，按照#chrom, start, end去重，保留refseq列中转录本编号最小的行。
    并返回去重后的DataFrame以及被去除的转录本列表。

    参数:
    df (pd.DataFrame): 输入的BED文件DataFrame，包含#chrom, start, end, refseq等列。

    返回:
    pd.DataFrame: 去重后的DataFrame。
    list: 被去除的转录本列表。
    """
    # 定义一个函数来提取转录本编号的数值部分
    def extract_transcript_number(refseq):
        if refseq.startswith('NM_'):
            return float(refseq[3:])
        return float('inf')  # 如果不是NM_开头的，返回无穷大，确保它们排在后面

    # 添加临时列，用于排序
    df['transcript_number'] = df['refseq'].apply(extract_transcript_number)

    # 添加一列，记录在原始文件中的行号
    df["origin_num"] = df.index + 1

    # 按照#chrom, start, end进行分组，并在每组中保留转录本编号最小的行
    # 同时记录被去除的行的refseq
    deduplicated_df = df.sort_values(by=['#chrom', 'start', 'end', 'transcript_number'])\
                        .drop_duplicates(subset=['#chrom', 'start', 'end'], keep='first')
    
    # 找出被去除的行
    removed_df = df[~df.index.isin(deduplicated_df.index)]
    removed_transcripts = removed_df['refseq'].tolist()

    # 删除临时列
    deduplicated_df = deduplicated_df.drop(columns=['transcript_number'])

    # 从整个DataFrame中去除被去除的转录本
    final_df = df[~df['refseq'].isin(removed_transcripts)]
    final_df = final_df.drop(columns=['transcript_number'])
    final_df.reset_index(drop=True, inplace=True)
    return final_df

# cds校正，当cds bed中存在相同的外显子编号，用cds的记录代替exon记录
def correct_start_end(df, cds_df):
    """
    根据cds_df校正df中的start和end列。
    当symbol列和location列相同时，使用cds_df的start和end替换df中的start和end。

    参数:
    df (pd.DataFrame): 原始的BED文件DataFrame。
    cds_df (pd.DataFrame): 用于校正的CDS文件DataFrame。

    返回:
    pd.DataFrame: 校正后的DataFrame。
    """
    # 创建一个字典，用于存储cds_df中的gene和location对应的start和end
    correction_dict = cds_df.set_index(['symbol', 'location'])[['start', 'end']].to_dict(orient='index')

    # 定义一个函数来校正start和end
    def correct_row(row):
        key = (row['symbol'], row['location'])
        if key in correction_dict:
            row['start'] = correction_dict[key]['start']
            row['end'] = correction_dict[key]['end']
        return row

    # 应用校正函数
    corrected_df = df.apply(correct_row, axis=1)
    
    # 自定义排序规则
    def chromosome_key(chrom):
        if chrom == "chrX":
            return 23
        elif chrom == "chrY":
            return 24
        else:
            return int(chrom[3:])  # 提取 chr 后面的数字

    # 按染色体和起始位置排序
    corrected_df["chromosome_order"] = corrected_df["#chrom"].apply(chromosome_key)
    corrected_df = corrected_df.sort_values(by=["chromosome_order", "start"])

    # 删除临时列
    corrected_df = corrected_df.drop(columns=["chromosome_order"])
    return corrected_df

# 最后的UTR修正
# 采用一刀切方案，在上面的位点的end坐标改变为下面start坐标-1
# 如果修正后，小于等于自己的start，就删除这个记录
def hard_filter(df, overlap):
    """
    根据overlap列表调整df的end列，并删除无效的行。

    参数:
    df (pd.DataFrame): 输入的DataFrame。
    overlap (list): 包含索引对的列表，例如 [(i1, j1), (i2, j2), ...]。

    返回:
    pd.DataFrame: 调整后的DataFrame。
    """
    # 记录需要删除的行索引
    to_drop = set()

    for pair in overlap:
        i, j = pair
        # 调整end列
        df.loc[i, 'end'] = df.loc[j, 'start'] - 1
        # 如果end <= start，标记为需要删除
        if df.loc[i, 'end'] <= df.loc[i, 'start']:
            to_drop.add(i)

    # 删除标记的行
    df.drop(index=to_drop, inplace=True)
    # 重置索引
    df.reset_index(drop=True, inplace=True)
    return df

# 主流程
def main(exon_bed, cds_bed, output_bed):
    input_df = read_bed(exon_bed)
    input_df_dedup = deduplicate_bed(input_df)
    cds_df = read_bed(cds_bed)
    cds_df_dedup = deduplicate_bed(cds_df)
    correct_df = correct_start_end(input_df_dedup, cds_df_dedup)

    n = 1
    while True:
        overlapping_pairs = find_overlapping_rows(correct_df)
        if len(overlapping_pairs) == 0:
            break
        print("[Hard Filter Run]", n)
        correct_df = hard_filter(correct_df, overlapping_pairs)
        n += 1

    # 结果输出
    correct_df = correct_df.drop(columns=["origin_num"])
    correct_df.to_csv(output_bed, sep="\t", header=True, index=False)

main('GRCh37/Gencode.GRCh37.exon.bed', "GRCh37/Gencode.GRCh37.cds.bed", "GRCh37/Gencode.GRCh37.exon.cor.bed")
main('GRCh38/Gencode.GRCh38.exon.bed', "GRCh38/Gencode.GRCh38.cds.bed", "GRCh38/Gencode.GRCh38.exon.cor.bed")

