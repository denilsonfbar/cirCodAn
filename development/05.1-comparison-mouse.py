from lib import libcc_dev as libcc
import os

datasets_path = 'data/datasets/'
results_file_path = 'results/05.tsv'
outputs_folder = 'outputs/05/'
os.makedirs(outputs_folder)

input_files_name = ['05_riboCirc_mmu.fa']
tools_name = libcc.get_all_tools_name()
parameters = {}

for input_file_name in input_files_name:
    parameters['input_file'] = input_file_name
    input_file_path = datasets_path + parameters['input_file']

    Codan_full = True
    for tool in tools_name:

        parameters['tool'] = tool
        parameters['model'] = None
        parameters['prob'] = None

        if parameters['tool'] == 'cirCodAn':
            parameters['model'] = 'VERT_circ'
            parameters['prob'] = 0.5

        elif parameters['tool'] == 'CodAn':
            if Codan_full:
                parameters['model'] = 'VERT_full'
                Codan_full = False
            else:
                parameters['model'] = 'VERT_partial'
        
        elif parameters['tool'] == 'CPAT':
            parameters['model'] = 'mouse'
    
        outputs_path = outputs_folder + parameters['input_file']  + '_' + parameters['tool']
        metrics = libcc.predict_cdss(parameters, input_file_path, outputs_path)

        libcc.add_prediction_cdss_result(results_file_path, parameters, metrics)
