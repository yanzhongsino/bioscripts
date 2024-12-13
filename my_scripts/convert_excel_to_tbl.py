#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
脚本：convert_excel_to_tbl.py
Author: Yan Zhong
Date: 2024-12-11

这个脚本是开发用来把excel格式的注释文件转化为NCBI GenBank的tbl格式。支持批量处理多个excel文件。已考虑单个excel文件中有多条染色体的情况。
脚本是在chatgpt4的帮助下完成的。

Script: convert_excel_to_tbl.py

This script converts one or more Excel files containing genomic feature annotations
into the GenBank TBL format. It's designed to handle multiple chromosomes within a single
Excel file and process multiple Excel files in a single run.

Usage:
    python convert_excel_to_tbl.py -i input1.xlsx [input2.xlsx ...] -o output1.tbl [output2.tbl ...]

Parameters:
    -i, --input     : One or more input Excel files. The files should contain columns with 
                      headers: 'chrs', 'start', 'end', 'gene', 'orientation', 'region', and 'note'.
                      The 'notes' column is optional.

    -o, --output    : One or more output TBL files.

Execution Examples:
    1. Convert multiple Excel files, specifying output filenames:
       python convert_excel_to_tbl.py -i input1.xlsx input2.xlsx -o output1.tbl output2.tbl

Input File Example:
    The Excel file should have a structure similar to the following:

    |  chrs  | start |  end  |  gene   | orientation |  region     |   note    |
    |--------|-------|-------|---------|-------------|-------------|-----------|
    | Chr01  |  313  | 5531  |  nad4   |      -      |   gene      |    nad4   |
    | Chr01  |  313  |  401  | nad4    |      -      |   CDS       | nad4_e2   |
    | Chr01  |  2891 | 5531  | nad4    |      -      |   CDS       | nad4_i1   |
    | Chr02  |  1621 | 1695  | nad4    |      +      |   gene      | trnD      |
    | Chr02  |  1621 | 1695  | nad4    |      +      |   tRNA      | trnD      |
    | Chr03  |  5645 | 5731  | atp1*   |      -      |gene fragment| atp1*     |
    | ...    |  ...  |  ...  |  ...    |     ...     |   ...       |   ...     |
    

    Each row corresponds to an annotated genomic feature. Ensure all necessary columns exist.

Requirements:
    - Python 3.x
    - pandas library
    - openpyxl library for reading Excel files

Notes:
    - Ensure that the number of output files (-o) matches the number of input files (-i) if specified.
    - The script uses the 'orientation' column to correctly orient features in the TBL file.
    - Chromosomes are indicated with the '>' directive in the TBL file.
    - The script writes the TBL file in standard GenBank format suitable for submission.

This script is particularly useful for genomic data analysis and submission to databases,
allowing for an automated and consistent process of converting Excel-based annotations to
a standard GenBank submission format.
"""

import pandas as pd
import argparse
import os

def convert_excel_to_tbl(input_file, output_file):
    # Read the Excel file
    df = pd.read_excel(input_file)

    with open(output_file, 'w') as tbl_file:
        current_chromosome = None
        current_gene = None

        for index, row in df.iterrows():
            chromosome = row['chrs']
            start = row['start']
            end = row['end']
            orientation = row['orientation']
            region = row['region']
            note = row['note'] if pd.notna(row['note']) else ''
            gene = row['gene']

            # Handle new chromosome
            if current_chromosome != chromosome:
                current_chromosome = chromosome
                current_gene = None
                tbl_file.write(f">Feature {chromosome}\n")

            # Adjust positions based on orientation
            if orientation == '-':
                start, end = end, start  # Reverse start and end for negative strand

            # Handle gene and gene fragment entries
            if region in ['gene', 'gene fragment']:
                if current_gene != gene:
                    current_gene = gene
                    tbl_file.write(f"{start}\t{end}\t{region}\n")
                    if region == 'gene':
                        tbl_file.write(f"                     gene            {gene}\n")
                    else:  # For 'gene fragment'
                        tbl_file.write(f"                     gene fragment   {gene}\n")
                    
                    if note:
                        tbl_file.write(f"                     note            {note}\n")  # Include note if present

            # Handle other regions
            else:
                tbl_file.write(f"{start}\t{end}\t{region}\n")
                if region == 'CDS':
                    tbl_file.write(f"                     product         {note}\n")  # Use note as product name
                elif region in ['exon', 'intron']:
                    if note:
                        tbl_file.write(f"                     note            {note}\n")  # Include note

def main(input_files, output_files):
    if len(input_files) != len(output_files):
        raise ValueError("The number of input files must match the number of output files.")
    
    for input_file, output_file in zip(input_files, output_files):
        convert_excel_to_tbl(input_file, output_file)
        print(f"Converted {input_file} to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert Excel files to GenBank TBL format.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='List of input Excel files')
    parser.add_argument('-o', '--output', nargs='+', required=True, help='List of output TBL files')

    args = parser.parse_args()

    main(args.input, args.output)