from lib import libcc_dev as libcc
import os
import pandas as pd

datasets_path = 'data/datasets/TransCirc_RP_PP_splitted/'
results_file_path = 'results/03.tsv'
outputs_folder = 'outputs/03/'
os.makedirs(outputs_folder)

parameters = {}
parameters['input_file'] = 'test_TransCirc_RP_PP_evidence_hsa.fa'
parameters['tool'] = 'cirCodAn'
parameters['model'] = 'VERT_circ'

parameters['prob'] = 0.5

for i in range(5):
    parameters['split'] = i
    print('\nSplit: ', i)
    split_path = datasets_path + str(i) + '/'

    input_file_path = split_path + parameters['input_file']

    outputs_path = outputs_folder + 'split-' + str(parameters['split']) + '-prob-' + str(parameters['prob'])
    metrics = libcc.predict_cdss(parameters, input_file_path, outputs_path)
    libcc.add_prediction_cdss_result(results_file_path, parameters, metrics)
