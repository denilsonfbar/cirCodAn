#!/usr/bin/env python3

# cirCodAn was developed from CodAn tool on August, 2023.
# cirCodAn predicts CDS in circRNA sequences
# Author: Denilson Fagundes Barbosa (denilsonfbar@gmail.com)

import os
import pandas as pd
import datetime as dt
import multiprocessing
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def count_sequences(fasta_file):
    count = 0
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            count += 1
    return count


def translate_orfs_in_file(input_file, output_file):
    
    file_translated = open(output_file,"w")

    for record in SeqIO.parse(input_file, "fasta"):
        id_seq = str(record.id)
        protein_seq = str(record.seq.translate())
        file_translated.write(">"+id_seq+"\n"+protein_seq+"\n")

    file_translated.close()


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

    source_gtf = "cirCodAn v1.0"
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

        circrna_seq_str = circrna_seqs_file_idx[circrna_id].seq

        orf_start_position = row['start'] -1  # index adapter

        feature_gtf = "start_codon"
        start_gtf = orf_start_position + 1  # index correction
        frame_gtf = (start_gtf-1) % 3

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
        frame_gtf = (start_gtf-1) % 3

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


def cirCodAn_output_changes(output_dir, circrna_fasta_file):

    gtf_file_CodAn_prediction = output_dir + "ORFs.gtf"
    new_gtf_file_name = output_dir + "CDS_predicted.gtf"

    fasta_file_CodAn_prediction = output_dir + "ORF_sequences.fasta"
    new_cds_fasta_file = output_dir + "CDS_predicted_seqs.fa"

    n_predictions = create_cds_circrna_predicted_gtf_file(gtf_file_CodAn_prediction, new_gtf_file_name, circrna_fasta_file, new_cds_fasta_file)

    translate_orfs_in_file(new_cds_fasta_file, output_dir + "CDS_predicted_seqs_aa.fa")

    if os.path.exists(fasta_file_CodAn_prediction):
        os.remove(fasta_file_CodAn_prediction)

    if os.path.exists(gtf_file_CodAn_prediction):
        os.remove(gtf_file_CodAn_prediction)

    if os.path.exists(output_dir + "3utr_sequences.fasta"):
        os.remove(output_dir + "3utr_sequences.fasta")

    if os.path.exists(output_dir + "5utr_sequences.fasta"):
        os.remove(output_dir + "5utr_sequences.fasta")
    
    return n_predictions


def __main__():

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file", help="Mandatory - input circRNAs file (FASTA format), /path/to/circRNA_seqs.fa", metavar="file", default=None)
    parser.add_option("-o", "--output", dest="output_folder", help="Optional - path to output folder, /path/to/output/folder/\nif not declared, it will be created at the circRNAs input folder [default=\"cirCodAn_output\"]", metavar="folder", default=None)
    parser.add_option("-m", "--model", dest="model_folder", help="Optional - path to model folder [default=\"models/VERT_circ\"]", metavar="folder", default=None)
    (options, args) = parser.parse_args()
    
    if options.file == None:
        print("""
cirCodAn v1.0

Use -h for help. 
Basic example to find CDS in circRNA sequences:

python3 circodan.py -f example/circRNA_seqs.fa

        """)
        quit()

    if options.file != None and os.path.isfile(options.file) == False:
        print("The circRNAs file indicated is not a valid file.")
        print('Please, indicate a valid FASTA file to the \"-f\" option.')
        quit()

    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" -> started cirCodAn v1.0")
    
    if options.output_folder == None:
        options.output_folder = ""
        folder_l = options.file.split("/")
        for i in range(0, len(folder_l)-1):
            options.output_folder += str(folder_l[i])+"/"
        options.output_folder += "cirCodAn_output/"
    elif options.output_folder.endswith("/") == False:
        options.output_folder += "/"
    if os.path.isdir(options.output_folder) == False:
        os.mkdir(options.output_folder)
        
    script_path = os.path.abspath(__file__)
    CodAn_path = script_path.replace("circodan.py","CodAn/bin")

    if options.model_folder == None:
        model_path = script_path.replace("circodan.py","models/VERT_circ")
    else:
        model_path = options.model_folder

    n_threads = multiprocessing.cpu_count()
    n_samples = count_sequences(options.file)

    call_os = 'export PATH=' + CodAn_path + ':$PATH' + ' && ' \
            + 'python3 ' + CodAn_path +'/codan.py -t ' + options.file \
            + ' -o ' + options.output_folder + ' -m ' + model_path + ' -s plus -c ' + str(n_threads) \
            + ' > /dev/null 2>&1'
    os.system(call_os)

    n_predictions = cirCodAn_output_changes(options.output_folder, options.file)

    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" -> prediction finished")
    print("Number of input sequences -> " + str(n_samples))
    print("Number of predicted CDSs  -> " + str(n_predictions))
    print("GTF file with prediction annotation -> " + options.output_folder + "CDS_predicted.gtf")
    print("Predicted CDS seqs FASTA file -> " + options.output_folder + "CDS_predicted_seqs.fa")
    print("Predicted peptides FASTA file -> " + options.output_folder + "CDS_predicted_seqs_aa.fa")


__main__()
