from lib import libdata
import os
import re
import time
import pandas as pd
from Bio import SeqIO


CIRCODAN_PATH = '/home/denilson/rep/cirCodAn/'
CODAN_PATH = '/home/denilson/rep/CodAn/'
ORFFINDER_PATH = '/home/denilson/rep/ORFfinder/'
FRAGGENESCAN_PATH = '/home/denilson/rep/FragGeneScan/'
CPAT_PATH = '/home/denilson/rep/CPAT/cpat-3.0.2/'
CPC2_PATH = '/home/denilson/rep/CPC2/CPC2_standalone-1.0.1/'


## Input files
#
def get_all_input_files_name():

    input_files_name = []

    input_files_name.append('01_TransCirc_no_evidence_hsa.fa')
    input_files_name.append('02_TransCirc_MS_evidence_hsa.fa')
    input_files_name.append('03_TransCirc_RP_PP_evidence_hsa.fa')
    input_files_name.append('04_riboCirc_hsa.fa')
    input_files_name.append('05_riboCirc_mmu.fa')
    input_files_name.append('06_MStoCIRC_hsa.fa')
    input_files_name.append('07_MStoCIRC_ath.fa')

    return input_files_name


## Tools
#
def get_all_tools_name():

    tools_name = []

    tools_name.append('cirCodAn')
    tools_name.append('CodAn')
    tools_name.append('CodAn')
    tools_name.append('ORFfinder')
    tools_name.append('FragGeneScan')
    tools_name.append('CPAT')
    tools_name.append('CPC2')
    
    return tools_name


## CDS prediction
#
def correct_in_replication(n_reps, input_file_path, df_pred_cdss):

    circrna_seqs_file = SeqIO.index(input_file_path, 'fasta')

    for i,row in df_pred_cdss.iterrows():

        circrna_len = int(len(circrna_seqs_file[row['seqname']].seq) / (n_reps + 1))

        while df_pred_cdss.at[i, 'end'] > circrna_len:
            df_pred_cdss.at[i, 'end'] = df_pred_cdss.at[i, 'end'] - circrna_len

    circrna_seqs_file.close()

    return df_pred_cdss

def predict_cdss(parameters, input_file_path, output_path):

    pred_start_time = time.time()

    if parameters['tool'] == 'cirCodAn':
        df_pred_cdss = run_cirCodAn(parameters, input_file_path, output_path)
    elif parameters['tool'] == 'CodAn':
        df_pred_cdss = run_CodAn(parameters, input_file_path, output_path)
    elif parameters['tool'] == 'ORFfinder':
        df_pred_cdss = run_ORFfinder(input_file_path, output_path)
    elif parameters['tool'] == 'FragGeneScan':
        df_pred_cdss = run_FragGeneScan(input_file_path, output_path)
    elif parameters['tool'] == 'CPAT':
        df_pred_cdss = run_CPAT(parameters, input_file_path, output_path)
    elif parameters['tool'] == 'CPC2':
        df_pred_cdss = run_CPC2(input_file_path, output_path)

    pred_end_time = time.time()
    pred_secs = round(pred_end_time - pred_start_time, 1)

    df_input_seqs = libdata.read_fasta_in_dataframe(input_file_path)

    df_real_cdss = pd.read_csv(input_file_path.replace('.fa','.gtf'), sep='\t', header=None)
    df_real_cdss.columns = libdata.get_gtf_columns()
    df_real_cdss = df_real_cdss[df_real_cdss['feature'] == 'cORF']

    if parameters['tool'] == 'CPAT':
        for i,row_pred in df_pred_cdss.iterrows():
            for _,row_seq in df_input_seqs.iterrows():
                if row_pred['seqname'] == row_seq['seqname'].upper():
                    df_pred_cdss.at[i,'seqname'] = row_seq['seqname']
                    break

    if 'n_reps' in parameters and parameters['n_reps'] > 0:
        df_pred_cdss = correct_in_replication(parameters['n_reps'], input_file_path, df_pred_cdss)

    df_pred_cdss.to_csv(output_path + '.tsv', sep='\t', index=False)

    metrics = evaluate_by_number_of_input_sequences(df_input_seqs, df_real_cdss, df_pred_cdss)
    metrics_cds = evaluate_by_number_of_cdss_in_input_sequences(df_real_cdss, df_pred_cdss)
    metrics = {**metrics, **metrics_cds}
    metrics['pred_secs'] = pred_secs

    return metrics


## cirCodAn
#
def set_cirCodAn_trans_prob_nocod_start(model_path, prob):

    model_main_file_path = model_path + '/ghmm/model/ghmm_intronless.model'
    model_main_file = open(model_main_file_path, 'r+')
    content = model_main_file.read()

    pattern_trans_prob_nocod_nocod = r'"nocod" \| "nocod": \d+\.\d+'
    pattern_trans_prob_nocod_start = r'"start" \| "nocod": \d+\.\d+'

    trans_prob_nocod_nocod = format(1.0-prob, '.10f')
    trans_prob_nocod_start = format(prob, '.10f')

    new_content = re.sub(pattern_trans_prob_nocod_nocod, f'"nocod" | "nocod": {trans_prob_nocod_nocod}', content)
    new_content = re.sub(pattern_trans_prob_nocod_start, f'"start" | "nocod": {trans_prob_nocod_start}', new_content)

    model_main_file.seek(0)
    model_main_file.write(new_content)
    model_main_file.truncate()
    model_main_file.close()

def run_cirCodAn(parameters, input_file_path, output_path):

    model_path = CIRCODAN_PATH + 'models/' + parameters['model']

    if 'prob' in parameters:
        set_cirCodAn_trans_prob_nocod_start(model_path, parameters['prob'])

    command = 'python3 ' + CIRCODAN_PATH + 'circodan.py' \
                + ' -f ' + input_file_path \
                + ' -o ' + output_path \
                + ' -m ' + model_path
    os.system(command)

    df_pred_cdss = pd.read_csv(output_path + '/CDS_predicted.gtf', sep='\t', header=None)
    df_pred_cdss.columns = libdata.get_gtf_columns()
    df_pred_cdss = df_pred_cdss[df_pred_cdss['feature'] == 'CDS']

    return df_pred_cdss  #  ['seqname', 'start', 'end', ...]


## CodAn
#
def get_predictions_from_CodAn_output(output_path):

    output_file_path =  output_path + '/ORFs.gtf'

    df_pred_cdss = pd.read_csv(output_file_path, sep='\t', header=None)
    df_pred_cdss.columns = libdata.get_gtf_columns()
    df_pred_cdss = df_pred_cdss[df_pred_cdss['feature'] == 'CDS']
    df_pred_cdss['end'] = df_pred_cdss['end'] + 3

    return df_pred_cdss

def run_CodAn(parameters, input_file_path, output_path):

    model_path = CODAN_PATH + 'models/' + parameters['model']

    command = 'export PATH=' + CODAN_PATH + '/bin:$PATH' + ' &&' \
              ' python3 ' + CODAN_PATH + 'bin/codan.py' \
            + ' -t ' + input_file_path \
            + ' -o ' + output_path \
            + ' -m ' + model_path\
            + ' -s plus'
    os.system(command)

    df_pred_cdss = get_predictions_from_CodAn_output(output_path)

    return df_pred_cdss  #  ['seqname', 'start', 'end', ...]


## ORFfinder
#
def get_predictions_from_ORFfinder_output(output_path):

    output_file_path =  output_path + '/ORFfinder_prediction.fa'

    ls_pred_cdss = []
    for orf_record in SeqIO.parse(output_file_path, 'fasta'):
        
        seqname_pattern = r'\|(.*?)\:'
        seqname = re.search(seqname_pattern, orf_record.id).group(1)

        start_pattern = r':(\d+)-'
        start = int(re.search(start_pattern, orf_record.id).group(1))

        end_pattern = r'-(\d+)'
        end = int(re.search(end_pattern, orf_record.id).group(1))

        ls_pred_cdss.append([seqname, start, end])

    column_names = ['seqname', 'start', 'end']
    df_pred_cdss = pd.DataFrame(ls_pred_cdss, columns=column_names)

    return df_pred_cdss

def run_ORFfinder(input_file_path, output_path):

    os.mkdir(output_path)

    command = ORFFINDER_PATH + 'ORFfinder' \
                + ' -in ' + input_file_path \
                + ' -out ' + output_path + '/ORFfinder_prediction.fa'\
                + ' -outfmt 1'\
                + ' -strand plus'
    os.system(command)

    df_pred_cdss = get_predictions_from_ORFfinder_output(output_path)

    return df_pred_cdss  #  ['seqname', 'start', 'end', ...]


## FragGeneScan
#
def get_predictions_from_FragGeneScan_output(output_path):

    output_file_path = output_path + '/FragGeneScan_prediction.gff'

    df_pred_cdss = pd.read_csv(output_file_path, sep='\t', header=None, comment='#')
    df_pred_cdss.columns = libdata.get_gtf_columns()
    df_pred_cdss = df_pred_cdss[df_pred_cdss['feature'] == 'CDS']

    return df_pred_cdss

def run_FragGeneScan(input_file_path, output_path):

    os.mkdir(output_path)
    output_file_path = output_path + '/FragGeneScan_prediction'

    command = FRAGGENESCAN_PATH + 'run_FragGeneScan.pl' \
                + ' -genome=' + input_file_path \
                + ' -out=' + output_file_path \
                + ' -complete=0'\
                + ' -train=complete'
    os.system(command)

    df_pred_cdss = get_predictions_from_FragGeneScan_output(output_path)

    return df_pred_cdss  #  ['seqname', 'start', 'end', ...]


## CPAT
#
def get_predictions_from_CPAT_output(output_path):

    output_file_path = output_path + '/CPAT_prediction.ORF_prob.best.tsv'

    df_pred_cdss = pd.read_csv(output_file_path, sep='\t', header=0)
    column_names = {'seq_ID':'seqname', 'ORF_start':'start', 'ORF_end':'end'}
    df_pred_cdss.rename(columns=column_names, inplace=True)

    return df_pred_cdss

def run_CPAT(parameters, input_file_path, output_path):

    os.mkdir(output_path)
    output_file_path = output_path + '/CPAT_prediction'

    hexamer_model_path = CPAT_PATH + 'dat/Human_Hexamer.tsv'
    logit_model_path = CPAT_PATH + 'dat/Human_logitModel.RData'
    if 'model' in parameters and parameters['model'] == 'mouse':
        hexamer_model_path = CPAT_PATH + 'dat/Mouse_Hexamer.tsv'
        logit_model_path = CPAT_PATH + 'dat/Mouse_logitModel.RData'

    command = CPAT_PATH + 'bin/cpat.py' \
                + ' -x ' + hexamer_model_path \
                + ' -d ' + logit_model_path \
                + ' -g ' + input_file_path \
                + ' -o ' + output_file_path
    os.system(command)

    df_pred_cdss = get_predictions_from_CPAT_output(output_path)

    return df_pred_cdss  #  ['seqname', 'start', 'end', ...]


## CPC2
#
def get_predictions_from_CPC2_output(output_path):

    output_file_path = output_path + '/CPC_prediction.txt'

    ls_pred_cdss = []
    df_output_file = pd.read_csv(output_file_path, sep='\t', header=0)

    for _,row in df_output_file.iterrows():

        if row['label'] == 'coding':

            start_cds = row['ORF_Start']
            end_cds = 0  # CPC does not return ORF final position

            ls_pred_cdss.append([row['#ID'], start_cds, end_cds])

    column_names = ['seqname', 'start', 'end']
    df_pred_cdss = pd.DataFrame(ls_pred_cdss, columns=column_names)

    return df_pred_cdss

def run_CPC2(input_file_path, output_path):

    os.mkdir(output_path)
    output_file_path = output_path + '/CPC_prediction'

    command = CPC2_PATH + 'bin/CPC2.py' \
                + ' -i ' + input_file_path \
                + ' -o ' + output_file_path \
                + ' --ORF'
    os.system(command)

    df_pred_cdss = get_predictions_from_CPC2_output(output_path)

    return df_pred_cdss  #  ['seqname', 'start', 'end', ...]


## Evaluation
#
''' 
Considering:
- input sequences file without repeated IDs
- one or more CDSs annotated in each positive input sequence
- input sequences file may have negative sequences (no CDS annotated)
- zero, one or more CDSs predicted for each input sequence
'''

'''
TP = positive sequence with at least one correctly predicted CDS
FP = positive or negative sequence only with incorrect predicted CDSs
FN = positive sequence without any predict CDS
TN = negative sequence without any predict CDS
'''
def evaluate_by_number_of_input_sequences(df_input_seqs, df_real_cdss, df_pred_cdss):
    
    metrics = {}

    # Input sequences (without repeated IDs)
    metrics['n_seqs'] = len(df_input_seqs)
    df_ids_pos_seqs = df_real_cdss.drop_duplicates(subset=['seqname'])[['seqname']]
    metrics['n_pos_seqs'] = len(df_ids_pos_seqs)
    metrics['n_neg_seqs'] = metrics['n_seqs'] - metrics['n_pos_seqs']

    # TP = positive sequence with at least one correctly predicted CDS (start and end positions)
    TP_start = len(pd.merge(df_real_cdss, df_pred_cdss, on=['seqname','start']).drop_duplicates(subset=['seqname']))
    TP_full = len(pd.merge(df_real_cdss, df_pred_cdss, on=['seqname','start','end']).drop_duplicates(subset=['seqname']))

    # FP = positive or negative sequence only with incorrect predicted CDSs (start and/or end positions)
    FP_start = len(df_pred_cdss.drop_duplicates(subset=['seqname'])) - TP_start
    FP_full = len(df_pred_cdss.drop_duplicates(subset=['seqname'])) - TP_full
    
    # FN = positive sequence without any predict CDS
    df_pos_without_predictions = pd.merge(df_ids_pos_seqs, df_pred_cdss, on='seqname', how='left', indicator=True)
    df_pos_without_predictions = df_pos_without_predictions[df_pos_without_predictions['_merge'] == 'left_only'].drop_duplicates(subset=['seqname'])
    FN = len(df_pos_without_predictions)

    # Calculating metrics
    try: precision_start =    TP_start / (TP_start + FP_start) * 100
    except: precision_start = 0.0
    try: recall_start =       TP_start / (TP_start + FN) * 100
    except: recall_start =    0.0
    try: f1_score_start =     2 * (precision_start * recall_start) / (precision_start + recall_start)
    except: f1_score_start =  0.0 

    try: precision_full =    TP_full / (TP_full + FP_full) * 100
    except: precision_full = 0.0
    try: recall_full =       TP_full / (TP_full + FN) * 100
    except: recall_full =    0.0
    try: f1_score_full =     2 * (precision_full * recall_full) / (precision_full + recall_full)
    except: f1_score_full =  0.0 
    

    TN = 0

    specificity_start = 0.0
    accuracy_start =    0.0

    specificity_full = 0.0
    accuracy_full =    0.0

    if metrics['n_neg_seqs'] > 0:

        # TN = negative sequence without any predict CDS
        df_neg_seqs = pd.merge(df_input_seqs, df_real_cdss, on='seqname', how='left', indicator=True)
        df_neg_seqs = df_neg_seqs[df_neg_seqs['_merge'] == 'left_only']
        df_neg_seqs = df_neg_seqs.drop(columns=['_merge'])
        df_tn = pd.merge(df_neg_seqs, df_pred_cdss, on='seqname', how='left', indicator=True)
        df_tn = df_tn[df_tn['_merge'] == 'left_only']
        TN = len(df_tn)

        try: specificity_start =    TN / (TN + FP_start) * 100
        except: specificity_start = 0.0
        try: accuracy_start =      (TP_start + TN) / (TP_start + TN + FP_start + FN) * 100
        except: accuracy_start =    0.0

        try: specificity_full =    TN / (TN + FP_full) * 100
        except: specificity_full = 0.0
        try: accuracy_full =      (TP_full + TN) / (TP_full + TN + FP_full + FN) * 100
        except: accuracy_full =    0.0

    metrics['FN_seq'] = FN
    metrics['TN_seq'] = TN

    metrics['TP_start_seq'] = TP_start
    metrics['FP_start_seq'] = FP_start
    metrics['Prec_start_seq'] = round(precision_start, 2)
    metrics['Rec_start_seq'] = round(recall_start, 2)
    metrics['F1-s_start_seq'] = round(f1_score_start, 2)
    metrics['Spec_start_seq'] = round(specificity_start, 2)
    metrics['Acc_start_seq'] = round(accuracy_start, 2)

    metrics['TP_full_seq'] = TP_full
    metrics['FP_full_seq'] = FP_full
    metrics['Prec_full_seq'] = round(precision_full, 2)
    metrics['Rec_full_seq'] = round(recall_full, 2)
    metrics['F1-s_full_seq'] = round(f1_score_full, 2)
    metrics['Spec_full_seq'] = round(specificity_full, 2)
    metrics['Acc_full_seq'] = round(accuracy_full, 2)

    return metrics

'''
TP = CDS correctly predicted (start and end positions)
FP = CDS incorrect predicted (start and/or end positions)
FN = CDS no predicted
TN = not applicable
'''
def evaluate_by_number_of_cdss_in_input_sequences(df_real_cdss, df_pred_cdss):
    
    metrics = {}

    # Real CDSs (zero, one or more for each input seq)
    metrics['n_real_cdss'] = len(df_real_cdss)

    # Predicted CDSs (zero, one or more for each input seq)
    metrics['n_pred_cdss'] = len(df_pred_cdss)

    # TP = CDS correctly predicted (start and end positions)
    TP_start = len(pd.merge(df_real_cdss, df_pred_cdss, on=['seqname','start']))
    TP_full = len(pd.merge(df_real_cdss, df_pred_cdss, on=['seqname','start','end']))

    # FP = CDS incorrect predicted (start and/or end positions)
    FP_start = metrics['n_pred_cdss'] - TP_start
    FP_full = metrics['n_pred_cdss'] - TP_full
    
    # FN = CDS no predicted
    FN_start = metrics['n_real_cdss'] - TP_start
    FN_full = metrics['n_real_cdss'] - TP_full

    # Calculating metrics
    try: precision_start =    TP_start / (TP_start + FP_start) * 100
    except: precision_start = 0.0
    try: recall_start =       TP_start / (TP_start + FN_start) * 100
    except: recall_start =    0.0
    try: f1_score_start =     2 * (precision_start * recall_start) / (precision_start + recall_start)
    except: f1_score_start =  0.0 

    try: precision_full =    TP_full / (TP_full + FP_full) * 100
    except: precision_full = 0.0
    try: recall_full =       TP_full / (TP_full + FN_full) * 100
    except: recall_full =    0.0
    try: f1_score_full =     2 * (precision_full * recall_full) / (precision_full + recall_full)
    except: f1_score_full =  0.0 

    metrics['TP_start_cds'] = TP_start
    metrics['FP_start_cds'] = FP_start
    metrics['FN_start_cds'] = FN_start
    metrics['Prec_start_cds'] = round(precision_start, 2)
    metrics['Rec_start_cds'] = round(recall_start, 2)
    metrics['F1-s_start_cds'] = round(f1_score_start, 2)

    metrics['TP_full_cds'] = TP_full
    metrics['FP_full_cds'] = FP_full
    metrics['FN_full_cds'] = FN_full
    metrics['Prec_full_cds'] = round(precision_full, 2)
    metrics['Rec_full_cds'] = round(recall_full, 2)
    metrics['F1-s_full_cds'] = round(f1_score_full, 2)

    return metrics

def add_prediction_cdss_result(results_file_path, parameters, metrics):

    log_prediction = {**parameters, **metrics}
    df_results = pd.DataFrame(log_prediction, index=[0])

    if os.path.exists(results_file_path):
        df_results_old = pd.read_csv(results_file_path, sep='\t', header=0)
        df_results = pd.concat([df_results_old, df_results], axis=0)

    df_results.to_csv(results_file_path, sep='\t', index=False)
