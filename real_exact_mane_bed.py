# coding=utf-8
# pzw
# 20241227
# 更精确的mane select转录本bed文件制作
# 风险点：只选择mane select转录本，会遗漏其他转录本的信息

import os
import requests
from urllib import request
from bs4 import BeautifulSoup
import gzip
from collections import defaultdict
import csv
import pandas as pd

# 检索网站获取最新的转录本
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v47.annotation.gff3.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/gencode.v47lift37.annotation.gff3.gz
def get_gencode_version():
    url = "https://www.gencodegenes.org/human/"
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        title = soup.title.string  # 获取标题内容
        # 从标题中提取版本号
        version = title.split(" - Human Release ")[-1]  # 假设版本号在最后
        return version.strip()
    return None

# ncbi版本
def get_ncbi_version():
    url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        # 找到MANE.GRCh38.v{version}.refseq_genomic.gff.gz文件的链接
        version_link = soup.find('a', href=lambda x: x and 'MANE.GRCh38.v' in x and 'refseq_genomic.gff.gz' in x)
        if version_link:
            # 从链接中提取版本号
            version = version_link['href'].split('MANE.GRCh38.v')[1].split('.refseq_genomic')[0]
            return version
    return None

# 下载器
def download_gencode_file(version, ncbi_version):
    path = 'database'
    release_file = 'release.txt'
    version_change = False
    db_version = "None"
    with open(release_file, 'r') as f:
        db_version = f.readline().strip()
    
    if db_version != version:
        # 下载文件
        hg19_url = f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/gencode.v{version}lift37.annotation.gff3.gz'
        hg38_url = f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v{version}.annotation.gff3.gz'
        transcript_match_file = f'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v{ncbi_version}.refseq_genomic.gff.gz'

        if not os.path.exists(path):
            os.makedirs(path)

        # 使用 urllib.request 下载文件
        request.urlretrieve(hg19_url, os.path.join(path, 'gencode.GRCh37.annotation.gff3.gz'))
        request.urlretrieve(hg38_url, os.path.join(path, 'gencode.GRCh38.annotation.gff3.gz'))
        request.urlretrieve(transcript_match_file, os.path.join(path, 'ncbi.GRCh38.gff.gz'))

        # 下载完成后，改写release文件
        with open(release_file, 'w') as f:
            f.write(version)
        version_change = True
    return version_change

# 从NCBI获得转录本对照表
def mane_select_match_refseq_ensembl(ncbi_gtf, output_txt):
    transcript_dict = {}
    
    with gzip.open(ncbi_gtf, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            lines = line.rstrip().split("\t")
            if lines[2] == "mRNA" and "tag=MANE Select" in line:
                infos = lines[8].split(";")
                info_dict = {k: v for k, v in (info.split("=") for info in infos if "=" in info)}
                
                ensembl = next((ens.replace("Ensembl:", "") for ens in info_dict.get('Dbxref', 'Ensembl:.').split(",") if ens.startswith("Ensembl:")), ".")
                refseq = info_dict.get("Name", ".")
                gene_name = info_dict.get("gene", ".")
                
                if gene_name != "." and ensembl != "." and refseq != ".":
                    transcript_dict[gene_name] = {"ensembl": ensembl, "refseq": refseq}

    # 输出
    with open(output_txt, "w", encoding="utf-8") as o:
        o.write("#gene\trefseq\tensembl\n")
        for gene, info in transcript_dict.items():
            o.write(f"{gene}\t{info['refseq']}\t{info['ensembl']}\n")

# 结果bed文件排序
def sort_bed_and_filter(input_bed, output_bed):
    # 读取 BED 文件
    df = pd.read_csv(input_bed, sep="\t", header=0, low_memory=False)

    # 过滤染色体，只保留 chr1-chr22, chrX, chrY
    valid_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    df = df[df["#chrom"].isin(valid_chromosomes)]

    # 自定义排序规则
    def chromosome_key(chrom):
        if chrom == "chrX":
            return 23
        elif chrom == "chrY":
            return 24
        else:
            return int(chrom[3:])  # 提取 chr 后面的数字

    # 按染色体和起始位置排序
    df["chromosome_order"] = df["#chrom"].apply(chromosome_key)
    df = df.sort_values(by=["chromosome_order", "start"])

    # 删除临时列
    df = df.drop(columns=["chromosome_order"])

    # 保存结果
    df.to_csv(output_bed, sep="\t", index=False, header=True)

# 解析gencode gtf
def gencode_gtf_extract(gencode_gtf, transcript, exon_output, cds_output, utr5_output, utr3_output):
    # 读取转录本信息
    transcript_dict = {}
    with open(transcript, "r") as t:
        for line in t:
            if not line.startswith("#"):
                lines = line.rstrip("\n").split("\t")
                # 不要求版本精确匹配
                transcript_dict[lines[2].split(".")[0]] = lines[1]
    
    # 打开输出文件
    with open(exon_output, "w", encoding="utf-8") as exon_output_open, \
         open(cds_output, "w", encoding="utf-8") as cds_output_open, \
         open(utr5_output, "w", encoding="utf-8") as utr5_output_open, \
         open(utr3_output, "w", encoding="utf-8") as utr3_output_open:

        # 标题抬头
        header = "#chrom\tstart\tend\tlocation\tsymbol\trefseq\tensembl\tstrand\n"
        exon_output_open.write(header)
        cds_output_open.write(header)
        utr5_output_open.write(header)
        utr3_output_open.write(header)

        # 遍历 GTF 文件
        with gzip.open(gencode_gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                lines = line.rstrip("\n").split("\t")
                chrom = lines[0]
                start = lines[3]
                end = lines[4]
                strand = lines[6]
                info_dict = {k: v for k, v in (info.split("=") for info in lines[8].split(";") if "=" in info)}
                ensembl = info_dict.get("Parent", ".")
                refseq = transcript_dict.get(ensembl.split(".")[0], ".")
                gene_name = info_dict.get("gene_name", ".")
                exon_location = info_dict.get("exon_number", ".")

                # 跳过找不到的信息
                if gene_name == "." or ensembl == "." or refseq == ".":
                    continue

                # 根据特征类型写入相应的文件
                output_line = "\t".join([chrom, start, end, exon_location, gene_name, refseq, ensembl, strand]) + "\n"
                if lines[2] == "exon":
                    exon_output_open.write(output_line)
                elif lines[2] == "CDS":
                    cds_output_open.write(output_line)
                elif lines[2] == "five_prime_UTR":
                    utr5_output_open.write(output_line)
                elif lines[2] == "three_prime_UTR":
                    utr3_output_open.write(output_line)

# 内含子bed生成
def create_intron_bed(exon_bed_file, intron_bed_file):
    # 读取外显子信息并按染色体和基因分组
    exon_dict = defaultdict(list)

    with open(exon_bed_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # 跳过标题行
        for row in reader:
            chrom, start, end, location, symbol, refseq, ensembl, strand = row
            exon_dict[(chrom, ensembl, strand)].append((int(start), int(end), int(location), symbol, refseq))

    # 构造内含子
    intron_list = []
    for (chrom, ensembl, strand), exons in exon_dict.items():
        # 跳过仅有一个外显子的基因
        if len(exons) < 2:
            continue

        # 按起始位置升序排列外显子
        exons_sorted = sorted(exons, key=lambda x: x[0])

        # 找到所有内含子
        introns = []
        for i in range(1, len(exons_sorted)):
            prev_exon_end = exons_sorted[i - 1][1]
            curr_exon_start = exons_sorted[i][0]
            # 确保外显子之间有空隙（即存在内含子）
            if curr_exon_start > prev_exon_end:
                intron_start = prev_exon_end + 1  # 前一个外显子的 end + 1
                intron_end = curr_exon_start - 1  # 后一个外显子的 start - 1
                introns.append((intron_start, intron_end))

        # 根据链的方向为内含子编号
        if strand == "+":
            for idx, (intron_start, intron_end) in enumerate(introns, start=1):
                intron_location = idx  # 小编号到大编号
                intron_list.append((chrom, intron_start, intron_end, intron_location, exons_sorted[0][3], exons_sorted[0][4], ensembl, strand))
        elif strand == "-":
            for idx, (intron_start, intron_end) in enumerate(introns, start=1):
                intron_location = len(introns) - idx + 1  # 大编号到小编号
                intron_list.append((chrom, intron_start, intron_end, intron_location, exons_sorted[0][3], exons_sorted[0][4], ensembl, strand))

    # 输出内含子 BED 文件
    with open(intron_bed_file, "w", encoding="utf-8") as f:
        f.write("#chrom\tstart\tend\tlocation\tsymbol\trefseq\tensembl\tstrand\n")
        for intron in intron_list:
            f.write("\t".join(map(str, intron)) + "\n")

# 所有bed
def all_bed_create(input_gtf, transcript, output_dir, basename):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    gencode_gtf_extract(input_gtf, transcript,
                        os.path.join(output_dir, "Gencode." + basename + ".exon.bed"),
                        os.path.join(output_dir, "Gencode." + basename + ".cds.bed"),
                        os.path.join(output_dir, "Gencode." + basename + ".utr5.bed"),
                        os.path.join(output_dir, "Gencode." + basename + ".utr3.bed"))
    for bed in os.listdir(output_dir):
        unsorted_bed = os.path.join(output_dir, bed)
        sort_bed_and_filter(unsorted_bed, unsorted_bed)
    create_intron_bed(os.path.join(output_dir, "Gencode." + basename + ".exon.bed"), os.path.join(output_dir, "Gencode." + basename + ".intron.bed"))


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

# 校正脚本
def correct_pipe(exon_bed, cds_bed, output_bed):
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

# 全流程
def main():
    # 检查当前版本
    gencode_version = get_gencode_version()
    ncbi_version = get_ncbi_version()

    # 确保都成功获取
    version_update = False
    if gencode_version and ncbi_version:
        version_update = download_gencode_file(gencode_version, ncbi_version)

    # 开始制作新的bed
    if version_update:
        # 制作新的转录本对照表
        mane_select_match_refseq_ensembl("database/ncbi.GRCh38.gff.gz", "transcript.txt")

        # GRCh37
        all_bed_create("database/gencode.GRCh37.annotation.gff3.gz", "transcript.txt", "GRCh37", "GRCh37")
        correct_pipe("GRCh37/Gencode.GRCh37.exon.bed", "GRCh37/Gencode.GRCh37.cds.bed", "GRCh37/Gencode.GRCh37.exon.cor.bed")
        
        # GRCh38
        all_bed_create("database/gencode.GRCh38.annotation.gff3.gz", "transcript.txt", "GRCh38", "GRCh38")
        correct_pipe("GRCh38/Gencode.GRCh38.exon.bed", "GRCh38/Gencode.GRCh38.cds.bed", "GRCh38/Gencode.GRCh38.exon.cor.bed")

if __name__ == "__main__":
    main()
