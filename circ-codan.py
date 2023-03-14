#!/usr/bin/env python3

# circCodAn was developed from CodAn tool on March, 2023.
# circCodAn predicts CDS in circRNA sequences
# Author: Denilson Fagundes Barbosa (denilsonfbar@gmail.com)

import os
import pandas as pd
import datetime as dt
import multiprocessing
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# This function considers zero-indexed sequences
def calcule_circ_orf_details(circrna_seq_str, orf_start_position):

    stop_codons = ['TAA', 'TAG', 'TGA']
    circrna_len = len(circrna_seq_str)

    orf_end = -1
    translation_cycles = 0
    orf_length = 0

    codon_test_start = orf_start_position
    codon_test_end = -1
    while True:
        
        codon_test_start += 3
        codon_test = ''
        orf_length += 3
        
        if codon_test_start >= circrna_len:

            translation_cycles += 1

            if codon_test_start == circrna_len: codon_test_start = 0
            elif codon_test_start == circrna_len + 1: codon_test_start = 1
            elif codon_test_start == circrna_len + 2: codon_test_start = 2

        codon_test_end = codon_test_start + 2
        if codon_test_end >= circrna_len:

            if codon_test_end == circrna_len: 
                codon_test_end = 0
                codon_test = circrna_seq_str[-2:] + circrna_seq_str[0]
            elif codon_test_end == circrna_len + 1: 
                codon_test_end = 1
                codon_test = circrna_seq_str[-1] + circrna_seq_str[0:1]
        else:
            codon_test = circrna_seq_str[codon_test_start:codon_test_end+1]

        if codon_test in stop_codons:
            orf_end = codon_test_end
            orf_length += 3
            break

        if translation_cycles == 4:  # stop_codon not located
            break

    if translation_cycles >= 4:  # stop_codon not located
        
        # Same treatment as Transcirc for infinite ORFs:
        orf_end = orf_start_position - 1
        if orf_end == -1: orf_end = circrna_len - 2
        translation_cycles = 3
        orf_length = circrna_len * translation_cycles
    
    elif orf_end == 0:
        translation_cycles += 1

    seq_translations = ""
    if translation_cycles > 1:
        seq_translations = circrna_seq_str * (translation_cycles-1)

    orf_seq_str = circrna_seq_str[orf_start_position : ] + seq_translations + circrna_seq_str[ : orf_end+1]

    return orf_seq_str, orf_end, translation_cycles, orf_length


def create_cds_circrna_predicted_gtf_file(gtf_file_CodAn_prediction, new_gtf_file_name, circrna_fasta_file, new_cds_fasta_file):

    df_gtf_CodAn = pd.read_csv(gtf_file_CodAn_prediction, sep='\t', header=None)
    df_gtf_CodAn.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    # Selecting only the start_codon predictions
    df_gtf_CodAn = df_gtf_CodAn[df_gtf_CodAn['feature'] == 'start_codon']

    orf_seq_record_list = []

    source_gtf = "circCodAn v1.0"
    score_gtf = "."
    attribute_gtf = "."

    predictions_counter = 0
    if os.path.exists(new_gtf_file_name):
        os.remove(new_gtf_file_name)
    
    circrna_seqs_file_idx = SeqIO.index(circrna_fasta_file, "fasta")

    for i,row in df_gtf_CodAn.iterrows():

        predictions_counter += 1

        circrna_id = row["seqname"]
        
        seqname_gtf = circrna_id
        strand_gtf = row["strand"]
        frame_gtf = row["frame"]

        circrna_seq_str = circrna_seqs_file_idx[circrna_id].seq

        orf_start_position = row['start'] -1  # index adapter

        feature_gtf = "start_codon"
        start_gtf = orf_start_position + 1  # index correction

        end_gtf = start_gtf + 2
        if end_gtf == len(circrna_seq_str) + 1: end_gtf = 1
        elif end_gtf == len(circrna_seq_str) + 2: end_gtf = 2

        with open(new_gtf_file_name, "a") as gtf_file:  # start_codon
            gtf_file.write(f"{seqname_gtf}\t{source_gtf}\t{feature_gtf}\t{start_gtf}\t{end_gtf}\t{score_gtf}\t{strand_gtf}\t{frame_gtf}\t{attribute_gtf}\n")

        orf_seq_str, orf_end_position, orf_trans_cycles, orf_length = calcule_circ_orf_details(circrna_seq_str, orf_start_position)

        orf_seq = Seq(str(orf_seq_str))
        orf_seq_record = SeqRecord(orf_seq, id=circrna_id, description='')
        orf_seq_record_list.append(orf_seq_record)

        feature_gtf = "CDS"
        end_gtf = orf_end_position + 1  # index correction
        attribute_gtf = "translation_cycles=" + str(orf_trans_cycles) + "; length=" + str(orf_length)
        with open(new_gtf_file_name, "a") as gtf_file:  # CDS
            gtf_file.write(f"{seqname_gtf}\t{source_gtf}\t{feature_gtf}\t{start_gtf}\t{end_gtf}\t{score_gtf}\t{strand_gtf}\t{frame_gtf}\t{attribute_gtf}\n")

        start_gtf = end_gtf - 2
        if start_gtf == 0: start_gtf = len(circrna_seq_str)
        elif start_gtf == -1: start_gtf = len(circrna_seq_str) - 1

        attribute_gtf = "."
        if orf_length == len(circrna_seq_str) * orf_trans_cycles:  # infinite ORF
            feature_gtf = "terminal_codon"
        else:
            feature_gtf = "stop_codon"
        with open(new_gtf_file_name, "a") as gtf_file:  # stop_codon OR terminal_codon
            gtf_file.write(f"{seqname_gtf}\t{source_gtf}\t{feature_gtf}\t{start_gtf}\t{end_gtf}\t{score_gtf}\t{strand_gtf}\t{frame_gtf}\t{attribute_gtf}\n")

    circrna_seqs_file_idx.close()

    SeqIO.write(orf_seq_record_list, new_cds_fasta_file, 'fasta')

    return predictions_counter


def circCodAn_output_changes(output_dir, circrna_fasta_file):

    gtf_file_CodAn_prediction = output_dir + "ORFs.gtf"
    new_gtf_file_name = output_dir + "CDS_prediction.gtf"

    fasta_file_CodAn_prediction = output_dir + "ORF_sequences.fasta"
    new_cds_fasta_file = output_dir + "CDS_predicted_seqs.fa"

    n_predictions = create_cds_circrna_predicted_gtf_file(gtf_file_CodAn_prediction, new_gtf_file_name, circrna_fasta_file, new_cds_fasta_file)

    if os.path.exists(fasta_file_CodAn_prediction):
        os.remove(fasta_file_CodAn_prediction)

    if os.path.exists(gtf_file_CodAn_prediction):
        os.remove(gtf_file_CodAn_prediction)

    if os.path.exists(output_dir + "3utr_sequences.fasta"):
        os.remove(output_dir + "3utr_sequences.fasta")

    if os.path.exists(output_dir + "5utr_sequences.fasta"):
        os.remove(output_dir + "5utr_sequences.fasta")


def __main__():

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file", help="Mandatory - input circRNAs file (FASTA format), /path/to/circRNA_seqs.fa", metavar="file", default=None)
    parser.add_option("-o", "--output", dest="output_folder", help="Optional - path to output folder, /path/to/output/folder/\nif not declared, it will be created at the circRNAs input folder [default=\"circCodAn_output\"]", metavar="folder", default=None)
    (options, args) = parser.parse_args()
    
    if options.file == None:
        print("""
circCodAn v1.0

Use -h for help. 
Basic usage to find CDS in circRNA sequences:

circ-codan.py -f circRNA_seqs.fa

        """)
        quit()

    if options.file != None and os.path.isfile(options.file) == False:
        print("The circRNAs file indicated is not a valid file.")
        print('Please, indicate a valid FASTA file to the \"-f\" option.')
        quit()

    if options.output_folder == None:
        options.output_folder = ""
        folder_l = options.file.split("/")
        for i in range(0, len(folder_l)-1):
            options.output_folder += str(folder_l[i])+"/"
        options.output_folder += "circCodAn_output/"
    elif options.output_folder.endswith("/") == False:
        options.output_folder += "/"
    if os.path.isdir(options.output_folder) == False:
        os.mkdir(options.output_folder)
    
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" -> started circCodAn v1.0")
    
    script_path = os.path.abspath(__file__)
    CodAn_path = script_path.replace("circ-codan.py","CodAn/bin")
    model_path = script_path.replace("circ-codan.py","models/VERT_circ")
    n_threads = multiprocessing.cpu_count()

    call_os = 'export PATH=' + CodAn_path + ':$PATH' + ' && ' \
            + 'python3 ' + CodAn_path +'/codan.py -t ' + options.file \
            + ' -o ' + options.output_folder + ' -m ' + model_path + ' -s plus -c ' + str(n_threads)
    os.system(call_os)

    circCodAn_output_changes(options.output_folder, options.file)

    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" -> prediction finished.")
    print("\tPredicted CDS sequences FASTA file  -> "+options.output_folder+"CDS_predicted_seqs.fa")
    print("\tGTF file with prediction annotation -> "+options.output_folder+"CDS_prediction.gtf")


__main__()
