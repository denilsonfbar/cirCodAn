import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns


CIRCODAN_PATH = '/home/denilson/rep/cirCodAn/'


## FASTA files
#

def count_sequences_in_file(fasta_file_path):

    count = 0
    with open(fasta_file_path, 'r') as file:
        for _ in SeqIO.parse(file, 'fasta'):
            count += 1

    return count

def read_fasta_in_dataframe(fasta_file_path):

    records = []

    for record in SeqIO.parse(fasta_file_path, 'fasta'):
        records.append({
            'seqname': record.id,
            'description': record.description,
            'seq_na': str(record.seq)
        })

    df_records = pd.DataFrame(records)

    return df_records

def remove_repeated_sequences_in_file(fasta_file_path):

    unique_seqs = set()
    records_to_keep = []

    for record in SeqIO.parse(fasta_file_path, 'fasta'):
        
        if str(record.seq) not in unique_seqs:
            unique_seqs.add(str(record.seq))
            records_to_keep.append(record)

    SeqIO.write(records_to_keep, fasta_file_path, 'fasta')

    return len(records_to_keep)

def replicate_sequences_in_file(n_replications, fasta_file_path, new_fasta_file_path):

    replicated_records = []

    for record in SeqIO.parse(fasta_file_path, 'fasta'):

        replicated_sequence = record.seq * (n_replications + 1)
        replicated_record = record
        replicated_record.seq = replicated_sequence
        replicated_records.append(replicated_record)

    SeqIO.write(replicated_records, new_fasta_file_path, 'fasta')


## ORF sequences
#

# Return -1 if start codon not located 
def get_orf_start_position_zi(orf_na, seq_na):

    max_subseq_align = 50
    seq_len = len(seq_na)

    if seq_len < max_subseq_align:
        length_subseq_align = seq_len
    else:
        length_subseq_align = max_subseq_align

    start_codon_position_zi = -1
    for i in range(length_subseq_align):

        orf_subseq = orf_na[i : length_subseq_align + i]
        start_align = seq_na.find(orf_subseq)

        if start_align == -1:
            continue
        else:
            start_codon_position_zi = start_align - i
            if start_codon_position_zi < 0:
                start_codon_position_zi = seq_len + start_codon_position_zi
            break

    return start_codon_position_zi

def get_corf_end_position_zi(orf_start_position_zi, seq_na):

    stop_codons = ['TAA', 'TAG', 'TGA']
    circrna_len = len(seq_na)
    translation_cycles = 0

    codon_test_start = orf_start_position_zi
    codon_test_end = -1
    while True:
        
        codon_test_start += 3
        codon_test = ''
        
        if codon_test_start >= circrna_len:

            translation_cycles += 1

            if codon_test_start == circrna_len: codon_test_start = 0
            elif codon_test_start == circrna_len + 1: codon_test_start = 1
            elif codon_test_start == circrna_len + 2: codon_test_start = 2

        codon_test_end = codon_test_start + 2
        if codon_test_end >= circrna_len:

            if codon_test_end == circrna_len: 
                codon_test_end = 0
                codon_test = seq_na[-2:] + seq_na[0]
            elif codon_test_end == circrna_len + 1: 
                codon_test_end = 1
                codon_test = seq_na[-1] + seq_na[:2]
        else:
            codon_test = seq_na[codon_test_start:codon_test_end+1]

        if codon_test in stop_codons:
            corf_end_position_zi = codon_test_end
            break

        if translation_cycles == 4:  # stop_codon not located
            corf_end_position_zi = orf_start_position_zi - 1
            if corf_end_position_zi == -1: corf_end_position_zi = circrna_len - 2
            break

    return corf_end_position_zi


## GTF files
#

def get_gtf_columns():

    return ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

def extract_gtf_attribute(attributes, key):

    attributes = attributes.split(";")
    
    for attribute in attributes:
    
        key_value = attribute.split("=")
        if len(key_value) == 2 and key_value[0] == key:
            return key_value[1]
    
    return None

def create_corf_gtf_annotation(circrna_id, circrna_na, orf_na, orf_start_position_zi, source_gtf):
 
    seqname_gtf = circrna_id
    feature_gtf = 'cORF'
    score_gtf = '.'
    strand_gtf = '+'

    start_gtf = orf_start_position_zi + 1  # index correction
    frame_gtf = (start_gtf-1) % 3

    corf_end_position_zi = get_corf_end_position_zi(orf_start_position_zi, circrna_na)
    end_gtf = corf_end_position_zi + 1  # index correction

    aux_orf_len = len(orf_na) + orf_start_position_zi - 1
    translation_cycles = aux_orf_len // len(circrna_na)

    attribute_gtf = 'translation_cycles=' + str(translation_cycles) + ';length=' + str(len(orf_na))

    return [seqname_gtf, source_gtf, feature_gtf, start_gtf, end_gtf, score_gtf, strand_gtf, frame_gtf, attribute_gtf]


## Visualization
#

def plot_line_chart(df, xscale, x_axis_col, y_axis_cols, width, height):

    dims = (width, height)
    plt.subplots(figsize=dims)
    plt.grid()

    for col in y_axis_cols:
        sns.lineplot(data=df, x=x_axis_col, y=col, legend='brief', label=col, marker='o')
        plt.xscale(xscale)


## TransCirc
#

def create_fasta_file_TransCirc(ids, source_fasta_file, new_fasta_file):

    circrna_seqs_file = SeqIO.index(source_fasta_file, 'fasta')
    select_seqs = []

    for id in ids:
        select_seqs.append(circrna_seqs_file[id])

    SeqIO.write(select_seqs, new_fasta_file, 'fasta')

    circrna_seqs_file.close()

def create_gtf_file_TransCirc(circrnas_fasta_file, orfs_fasta_file, gtf_file):

    ls_gtf = []
    source_gtf = 'TransCirc v1.0'

    orf_id_sufixes = ['_ORF1','_ORF2','_ORF3']
    orfs_fasta_idx = SeqIO.index(orfs_fasta_file, 'fasta')

    circrna_recs = list(SeqIO.parse(circrnas_fasta_file, 'fasta'))
    for circrna_rec in circrna_recs:

        for orf_sufix in orf_id_sufixes:
            orf_id = circrna_rec.id + orf_sufix
            try:
                orf_na = str(orfs_fasta_idx[orf_id].seq)
            except KeyError:
                if orf_sufix == '_ORF1':
                    print("circRNA with no ORF annotation: ", circrna_rec.id)
                break

            circrna_na = str(circrna_rec.seq)

            start_codon_position_zi = get_orf_start_position_zi(orf_na, circrna_na)
            
            if start_codon_position_zi == -1:
                print('Start codon not located: ', circrna_rec.id, orf_id)
                continue

            ls_gtf.append(create_corf_gtf_annotation(circrna_rec.id, circrna_na, orf_na, start_codon_position_zi, source_gtf))
    
    orfs_fasta_idx.close()

    df_gtf = pd.DataFrame(ls_gtf, columns=get_gtf_columns())
    df_gtf['start'] = df_gtf['start'].astype(int)
    df_gtf['end'] = df_gtf['end'].astype(int)

    df_gtf.to_csv(gtf_file, sep='\t', index=False, header=False)


## riboCirc
#

def create_fasta_file_riboCirc(ids, source_fasta_file, new_fasta_file):

    select_seqs = []
    
    for id in ids:
    
        for record in SeqIO.parse(source_fasta_file, "fasta"):

            if id in str(record.id):
        
                sel_seq = Seq(str(record.seq))
                sel_seq_record = SeqRecord(sel_seq, id=id, description='')
                select_seqs.append(sel_seq_record)
                
                break
    
    SeqIO.write(select_seqs, new_fasta_file, 'fasta')

def create_gtf_file_riboCirc(circrna_fasta_file, orf_fasta_file, gtf_file):

    ls_gtf = []
    source_gtf = 'riboCirc v1.0'

    for circrna_rec in SeqIO.parse(circrna_fasta_file, "fasta"):

        for orf_rec in SeqIO.parse(orf_fasta_file, "fasta"):

            if circrna_rec.id in str(orf_rec.id):

                circrna_na = str(circrna_rec.seq)
                orf_na = str(orf_rec.seq)                

                start_codon_position_zi = get_orf_start_position_zi(orf_na, circrna_na)
            
                if start_codon_position_zi == -1:
                    print('Start codon not located: ', circrna_rec.id, orf_rec.id)
                    continue

                ls_gtf.append(create_corf_gtf_annotation(circrna_rec.id, circrna_na, orf_na, start_codon_position_zi, source_gtf))

    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df_gtf = pd.DataFrame(ls_gtf, columns=gtf_columns)
    df_gtf['start'] = df_gtf['start'].astype(int)
    df_gtf['end'] = df_gtf['end'].astype(int)

    df_gtf = df_gtf.drop_duplicates(subset=['seqname', 'start'])

    df_gtf.to_csv(gtf_file, sep='\t', index=False, header=False)


## MStoCIRC
#

def create_fasta_file_MStoCIRC(ids, sequences, fasta_file):

    seq_records = []
    
    for i in range(len(ids)):
    
        seq = Seq(sequences[i])
        seq_record = SeqRecord(seq, id=ids[i], description='')
        seq_records.append(seq_record)
    
    SeqIO.write(seq_records, fasta_file, 'fasta')

def create_gtf_file_MStoCIRC(df_data, gtf_file_name):

    with open(gtf_file_name, 'w') as file:

        for _, row in df_data.iterrows():

            seqname = row['circRNAname']
            source = 'MStoCIRC'
            feature = 'cORF'
            start = row['CORF_start']
            end = row['CORF_end_corrected']
            score = '.'
            strand = '+'
            frame = '.'
            attributes = f'translation_cycles={row["translation_cycles"]}'

            line = f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n"
            
            file.write(line)
