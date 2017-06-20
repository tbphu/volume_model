import numpy as np
import tellurium as te
import pandas as pd
import pickle


mothercells_file_path = 'data/mothercells.p'
daughtercells_data_path = 'data/daughtercells.p'
time_file_path = "data/time.p"

def load_data():
    mothercells_data = pickle.load(open(mothercells_file_path, "rb"))
    daughtercells_data = pickle.load(open(daughtercells_data_path, "rb"))
    time_data = pickle.load(open(time_file_path, "rb"))
    return mothercells_data, daughtercells_data, time_data

def experimental_data_to_dict(time_vec, data_matrix):
    data_dict = dict()
    data_dict['time'] = time_vec
    for pos in range(len(data_matrix)):
        data_dict[pos] = data_matrix[pos]
    return data_dict

def plot_experimental_data(mothercells_data, daughtercells_data, time_data, subplot=True, show=True):
    mothercells_dict = model_data.experimental_data_to_dict(time_data, mothercells_data)
    daughtercells_dict = model_data.experimental_data_to_dict(time_data, daughtercells_data)
    data_tuple = (mothercells_dict, daughtercells_dict)
    plot(data_tuple, subplot, show)

def truncate_data(data):
    """ truncate data in such a way that there are no leading or trailing nan values left"""
    pos_min = np.inf
    pos_max = 0
    for variable in data:
        if variable == 'time':
            continue
        data_vec_wo_nan = np.where(~np.isnan(data[variable]))
        start_pos = data_vec_wo_nan[0][0]
        end_pos = data_vec_wo_nan[0][-1]
        if start_pos < pos_min:
            pos_min = start_pos
        if end_pos > pos_max:
            pos_max = end_pos
    for variable in data:
        data[variable] = data[variable][pos_min:pos_max]
    return data

def get_initial_values_from_data(data):
    initial_values = {}
    for variable in data:
        if variable == 'time':
            continue
        initial_values[variable] = data[variable][0]
    return initial_values

def get_model_parameters_as_dict(model):
    param_to_value = dict()
    for key in model.iterkeys():
        if key.startswith('init') or key.startswith('eigen'):
            continue
        try:
            param_to_value[key] = model[key]
        except RuntimeError:
            pass
    return param_to_value

def get_model_parameters_as_dataframe(model):
    param_to_value = get_model_parameters_as_dict(model)
    df = pd.DataFrame([param_to_value])
    return df

def set_model_parameters(model, params):
    for param_id in params:
        try:
            model[param_id] = params[param_id]
        except RuntimeError:
            print('could not set parameter : {0}'.format(param_id))
                
    return model

def set_model_parameters_from_dataframe(model, df, param_list=[],row=0):
    param_dict = dict(df.iloc[row])
    if param_list != []:
        param_dict = {key: param_dict[key] for key in param_dict if key in param_list}
    return set_model_parameters(model, param_dict)




