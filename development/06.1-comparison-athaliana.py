from lib import libcc_dev as libcc
import os

datasets_path = 'data/datasets/'
results_file_path = 'results/06.tsv'
outputs_folder = 'outputs/06/'
os.makedirs(outputs_folder)

input_files_name = ['07_MStoCIRC_ath.fa']
tools_name = libcc.get_all_tools_name()
tools_name.insert(1, 'cirCodAn')
parameters = {}

for input_file_name in input_files_name:
    parameters['input_file'] = input_file_name
    input_file_path = datasets_path + parameters['input_file']

    cirCodan_plants = True
    Codan_full = True
    for tool in tools_name:

        parameters['tool'] = tool
        parameters['model'] = None
        parameters['prob'] = None

        if parameters['tool'] == 'cirCodAn':
            parameters['prob'] = 0.5
            if cirCodan_plants:
                parameters['model'] = 'PLANTS_circ'
                cirCodan_plants = False
            else:
                parameters['model'] = 'VERT_circ'

        elif parameters['tool'] == 'CodAn':
            if Codan_full:
                parameters['model'] = 'PLANTS_full'
                Codan_full = False
            else:
                parameters['model'] = 'PLANTS_partial'
        
        outputs_path = outputs_folder + parameters['input_file']  + '_' + parameters['tool']
        metrics = libcc.predict_cdss(parameters, input_file_path, outputs_path)

        libcc.add_prediction_cdss_result(results_file_path, parameters, metrics)
