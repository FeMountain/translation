#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生物序列翻译工具 - 命令行版本
支持从NCBI下载序列并翻译为蛋白质
"""

import argparse
import sys
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import re
import pandas as pd

def clean_sequence(sequence):
    """清理序列，去除空格和非法字符"""
    cleaned = re.sub(r'\s+', '', sequence)
    cleaned = re.sub(r'[^ATCGUatcgu]', '', cleaned)
    return cleaned.upper()

def download_sequence(accession_id, email):
    """从NCBI下载序列"""
    Entrez.email = email
    try:
        print(f"正在从NCBI下载序列: {accession_id}")
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        print(f"✅ 成功下载: {record.description}")
        return str(record.seq), record.description
    except Exception as e:
        print(f"❌ 下载失败: {str(e)}")
        return None, None

def translate_sequence(sequence, start_position=0, reading_frame=1):
    """翻译DNA/RNA序列为蛋白质"""
    try:
        seq_obj = Seq(sequence)
        
        # 根据阅读框调整起始位置
        if reading_frame == 2:
            start_position += 1
        elif reading_frame == 3:
            start_position += 2
        
        # 从指定位置开始翻译
        if start_position > 0:
            seq_obj = seq_obj[start_position:]
        
        protein_seq = seq_obj.translate(to_stop=False)
        return str(protein_seq)
    except Exception as e:
        print(f"❌ 翻译失败: {str(e)}")
        return None

def get_codon_table(sequence, start_position=0, reading_frame=1):
    """生成密码子到氨基酸的对照表"""
    try:
        # 调整起始位置
        if reading_frame == 2:
            start_position += 1
        elif reading_frame == 3:
            start_position += 2
        
        if start_position > 0:
            sequence = sequence[start_position:]
        
        # 确保序列长度是3的倍数
        if len(sequence) % 3 != 0:
            sequence = sequence[:-(len(sequence) % 3)]
        
        codons = []
        amino_acids = []
        positions = []
        
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                temp_seq = Seq(codon)
                aa = str(temp_seq.translate(to_stop=False))
                
                codons.append(codon)
                amino_acids.append(aa)
                positions.append(f"{start_position + i + 1}-{start_position + i + 3}")
        
        return pd.DataFrame({
            '位置': positions,
            '密码子': codons,
            '氨基酸': amino_acids
        })
    except Exception as e:
        print(f"❌ 生成密码子表失败: {str(e)}")
        return None

def save_results(protein_seq, codon_table, output_prefix):
    """保存结果到文件"""
    # 保存蛋白质序列
    if protein_seq:
        with open(f"{output_prefix}_protein.fasta", "w") as f:
            f.write(f">translated_protein\n{protein_seq}\n")
        print(f"✅ 蛋白质序列已保存到: {output_prefix}_protein.fasta")
    
    # 保存密码子表
    if codon_table is not None:
        codon_table.to_csv(f"{output_prefix}_codon_table.csv", index=False, encoding='utf-8-sig')
        print(f"✅ 密码子表已保存到: {output_prefix}_codon_table.csv")

def main():
    parser = argparse.ArgumentParser(description="生物序列翻译工具")
    parser.add_argument("--accession", "-a", help="NCBI Accession ID")
    parser.add_argument("--sequence", "-s", help="直接输入DNA/RNA序列")
    parser.add_argument("--email", "-e", required=True, help="NCBI邮箱地址（必需）")
    parser.add_argument("--start", "-p", type=int, default=0, help="起始位置（默认: 0）")
    parser.add_argument("--frame", "-f", type=int, choices=[1, 2, 3], default=1, help="阅读框（默认: 1）")
    parser.add_argument("--output", "-o", default="result", help="输出文件前缀（默认: result）")
    
    args = parser.parse_args()
    
    if not args.accession and not args.sequence:
        print("❌ 请提供Accession ID或序列")
        parser.print_help()
        sys.exit(1)
    
    # 获取序列
    sequence = None
    description = None
    
    if args.accession:
        sequence, description = download_sequence(args.accession, args.email)
    else:
        sequence = clean_sequence(args.sequence)
        description = "用户输入序列"
    
    if not sequence:
        print("❌ 无法获取序列")
        sys.exit(1)
    
    print(f"\n📊 序列信息:")
    print(f"描述: {description}")
    print(f"长度: {len(sequence)} bp")
    print(f"GC含量: {(sequence.count('G') + sequence.count('C')) / len(sequence) * 100:.1f}%")
    print(f"起始位置: {args.start}")
    print(f"阅读框: {args.frame}")
    
    # 翻译序列
    print(f"\n🔄 翻译结果:")
    protein_seq = translate_sequence(sequence, args.start, args.frame)
    
    if protein_seq:
        print(f"蛋白质长度: {len(protein_seq)} aa")
        print(f"蛋白质序列: {protein_seq}")
        
        # 统计氨基酸组成
        aa_counts = {}
        for aa in protein_seq:
            if aa != '*':
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        if aa_counts:
            print(f"\n氨基酸组成:")
            for aa, count in sorted(aa_counts.items()):
                print(f"  {aa}: {count}")
        
        # 生成密码子表
        print(f"\n📋 密码子对照表:")
        codon_table = get_codon_table(sequence, args.start, args.frame)
        
        if codon_table is not None:
            print(codon_table.to_string(index=False))
            
            # 保存结果
            save_results(protein_seq, codon_table, args.output)
    else:
        print("❌ 翻译失败")
        sys.exit(1)

if __name__ == "__main__":
    main() 