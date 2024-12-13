#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
convert_msa_to_indels.py
Author: Yan Zhong
Date: 2024-11-06

这个脚本是开发用来处理多序列比对（multiple sequences alignment, MSA）文件的。默认把多序列比对文件的第一条序列作为参考序列，统计其他序列的indels的数量和每个indel的长度。支持批量处理多个MSA文件。已考虑参考序列和当前序列都是deletion的情况。
脚本是在chatgpt4的帮助下完成的。
用法: python convert_msa_to_indels.py /path/to/msa.fasta --min_indel_length 5 --reference_index 1
参数：
- <msa.fasta>：必需参数，多序列比对文件，支持多个MSA文件。
- --min_indel_length 0：可选参数，设置被统计的indel的最短长度，默认是0。
- --reference_index 0：可选参数，通过index设置参考序列，默认是0，代表多序列比对文件的第一条序列为参考序列。
- --output indel_statistics.tsv：指定结果的输出文件。不指定则输出到indel_statistics.tsv文件。
- -h：帮助列表。

This script counts insertions and deletions (indels) in multiple sequence alignment (MSA) files
compared to a reference sequence and outputs the statistics to a TSV file.

usage: convert_msa_to_indels.py [-h] [--min_indel_length MIN_INDEL_LENGTH]
                           [--reference_index REFERENCE_INDEX]
                           [--output OUTPUT]
                           file_paths [file_paths ...]

Count indels in multiple MSA files compared to a reference sequence.

positional arguments:
  file_paths            List of paths to the MSA files in FASTA format.

optional arguments:
  -h, --help            show this help message and exit
  --min_indel_length MIN_INDEL_LENGTH
                        Minimum indel length to consider.
  --reference_index REFERENCE_INDEX
                        Index of the reference sequence in the alignment.
  --output OUTPUT       Output CSV file name.
"""


import sys
import os
import argparse
import csv
from Bio import AlignIO

def count_indels(msa_file, min_indel_length=0, reference_index=0):
    alignment = AlignIO.read(msa_file, "fasta")
    try:
        reference_seq = alignment[reference_index].seq
    except IndexError:
        print("Error: Reference index out of range.")
        return []

    indel_stats = []

    for i, record in enumerate(alignment):
        if i == reference_index:
            continue

        current_seq = record.seq

        all_del_lengths = []
        all_ins_lengths = []

        del_length = 0
        ins_length = 0
        prev_state = None

        for ref, curr in zip(reference_seq, current_seq):
            if ref == '-' and curr == '-':
                continue
            # Both are gaps, continue with current indel state

            if ref != '-' and curr == '-':
                # Deletion
                if prev_state == 'insertion':
                    if ins_length > 0:
                        all_ins_lengths.append(ins_length)
                    ins_length = 0
                del_length += 1
                prev_state = 'deletion'

            elif ref == '-' and curr != '-':
                # Insertion
                if prev_state == 'deletion':
                    if del_length > 0:
                        all_del_lengths.append(del_length)
                    del_length = 0
                ins_length += 1
                prev_state = 'insertion'

            else:
                # Encounter non-gap in both sequences, finalize any ongoing indel
                if prev_state == 'deletion' and del_length > 0:
                    all_del_lengths.append(del_length)
                    del_length = 0
                elif prev_state == 'insertion' and ins_length > 0:
                    all_ins_lengths.append(ins_length)
                    ins_length = 0

                prev_state = None
        
        # Finalize indel lengths if still ongoing at the end of sequences
        if prev_state == 'deletion' and del_length > 0:
            all_del_lengths.append(del_length)
        elif prev_state == 'insertion' and ins_length > 0:
            all_ins_lengths.append(ins_length)

        # Filter indels by min_indel_length
        del_lengths = [length for length in all_del_lengths if length >= min_indel_length]
        ins_lengths = [length for length in all_ins_lengths if length >= min_indel_length]

        indel_stats.append({
            'file_name': os.path.basename(msa_file),
            'sequence_id': record.id,
            'deletions_count': len(del_lengths),
            'insertions_count': len(ins_lengths),
            'deletions_lengths': del_lengths,
            'insertions_lengths': ins_lengths
        })

    return indel_stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count indels in multiple MSA files compared to a reference sequence.")
    parser.add_argument("file_paths", type=str, nargs='+', help="List of paths to the MSA files in FASTA format.")
    parser.add_argument("--min_indel_length", type=int, default=0, help="Minimum indel length to consider.")
    parser.add_argument("--reference_index", type=int, default=0, help="Index of the reference sequence in the alignment.")
    parser.add_argument("--output", type=str, default="indel_statistics.tsv", help="Output CSV file name.")

    args = parser.parse_args()

    all_stats = []
    for file_path in args.file_paths:
        indel_stats = count_indels(file_path, args.min_indel_length, args.reference_index)
        all_stats.extend(indel_stats)

    # Write results to a TSV file
    with open(args.output, mode='w', newline='') as tsvfile:
        fieldnames = ['file_name', 'sequence_id', 'deletions_count', 'insertions_count', 'deletions_lengths', 'insertions_lengths']
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for stat in all_stats:
            writer.writerow(stat)

    print(f"Results have been written to {args.output}.")