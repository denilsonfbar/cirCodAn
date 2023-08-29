from lib import libcc_dev as libcc
import os

datasets_path = 'data/datasets/'
results_file_path = 'results/08.tsv'
outputs_folder = 'outputs/08/'
os.makedirs(outputs_folder)

input_files_name = libcc.get_all_input_files_name()
input_files_name.pop(0)
input_files_name.pop(1)
tools_name = ['cirCodAn']
parameters = {}

for input_file_name in input_files_name:
    parameters['input_file'] = input_file_name
    input_file_path = datasets_path + parameters['input_file']

    for tool in tools_name:

        parameters['tool'] = tool
        parameters['model'] = None
        parameters['prob'] = None

        if parameters['tool'] == 'cirCodAn':
            parameters['prob'] = 0.5
            parameters['model'] = 'VERT_circ'
            if parameters['input_file'] == '07_MStoCIRC_ath.fa':
                parameters['model'] = 'PLANTS_circ'

        outputs_path = outputs_folder + parameters['input_file']  + '_' + parameters['tool']
        metrics = libcc.predict_cdss(parameters, input_file_path, outputs_path)

        libcc.add_prediction_cdss_result(results_file_path, parameters, metrics)
