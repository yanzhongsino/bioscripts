#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
msa2indels.py
Author: Yan Zhong
Date: 2024-11-06

This script was developed to handle multiple sequence alignment (MSA) file. By default, the first sequence of the multi-sequence alignment file is used as the reference sequence, and the number of indels (insertion or deletion) of other sequences and the length of each indel are counted.The case where both the reference sequence and the counted sequence are deletion has been considered.
The script was done with the help of chatGPT4.
这个脚本是开发用来处理多序列比对文件的。默认把多序列比对文件的第一条序列作为参考序列，统计其他序列的indels的数量和每个indel的长度。已考虑参考序列和当前序列都是deletion的情况。
脚本是在chatgpt4的帮助下完成的。
Usage: python msa2indels.py /path/to/msa.fasta --min_indel_length 5 --reference_index 1
parameters：
- <msa.fasta>：必需参数，脚本后第一个参数为多序列比对文件。
- --min_indel_length 5：可选参数，脚本后第二个参数，设置被统计的indel的最短长度，默认是0。
- --reference_index 0：可选参数，脚本后第三个参数，设置参考序列的顺序，默认是0，代表多序列比对文件的第一条序列。
- -h：show this help message and exit.
"""


import sys
import os
import argparse
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
            'sequence_id': record.id,
            'deletions_count': len(del_lengths),
            'insertions_count': len(ins_lengths),
            'deletions_lengths': del_lengths,
            'insertions_lengths': ins_lengths
        })

    return indel_stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count indels in an MSA file compared to a reference sequence.")
    parser.add_argument("file_path", type=str, help="Path to the MSA file in FASTA format.")
    parser.add_argument("--min_indel_length", type=int, default=0, help="Minimum indel length to consider.")
    parser.add_argument("--reference_index", type=int, default=0, help="Index of the reference sequence in the alignment.")

    args = parser.parse_args()

    indel_stats = count_indels(args.file_path, args.min_indel_length, args.reference_index)
    for stat in indel_stats:
        print(stat)