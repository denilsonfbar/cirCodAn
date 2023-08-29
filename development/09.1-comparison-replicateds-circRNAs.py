from lib import libcc_dev as libcc
import os

replicateds_path = libcc.CIRCODAN_PATH + 'datasets/replicateds/'
results_file_path = 'results/09.tsv'
outputs_folder = 'outputs/09/'
os.makedirs(outputs_folder)

input_files_name = libcc.get_all_input_files_name()
input_files_name.pop(0)
input_files_name.pop(1)
tools_name = libcc.get_all_tools_name()
tools_name.pop(6)
parameters = {}

for n_replications in range(1, 4):

    parameters['n_reps'] = n_replications
    datasets_path = replicateds_path + 'n_reps_' + str(n_replications) + '/'

    for input_file_name in input_files_name:
        parameters['input_file'] = input_file_name
        input_file_path = datasets_path + parameters['input_file']

        Codan_full = True
        for tool in tools_name:

            parameters['tool'] = tool
            parameters['model'] = None
            parameters['prob'] = None

            if parameters['tool'] == 'cirCodAn':
                parameters['prob'] = 0.5
                parameters['model'] = 'VERT_circ'
                if parameters['input_file'] == '07_MStoCIRC_ath.fa':
                    parameters['model'] = 'PLANTS_circ'

            elif parameters['tool'] == 'CodAn':
                if Codan_full:
                    parameters['model'] = 'VERT_full'
                    Codan_full = False
                    if parameters['input_file'] == '07_MStoCIRC_ath.fa':
                        parameters['model'] = 'PLANTS_full'
                else:
                    parameters['model'] = 'VERT_partial'
                    if parameters['input_file'] == '07_MStoCIRC_ath.fa':
                        parameters['model'] = 'PLANTS_partial'
            
            elif parameters['tool'] == 'CPAT':
                if parameters['input_file'] == '05_riboCirc_mmu.fa':
                    parameters['model'] = 'mouse'
            
            outputs_path = outputs_folder + str(parameters['n_reps']) + '_' + parameters['input_file']  + '_' + parameters['tool']
            metrics = libcc.predict_cdss(parameters, input_file_path, outputs_path)

            libcc.add_prediction_cdss_result(results_file_path, parameters, metrics)
