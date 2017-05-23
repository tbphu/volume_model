#!/usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
from scipy.optimize import basinhopping


mothercells_file_path = 'data/mothercells.p'
daughtercells_data_path = 'data/daughtercells.p'
time_file_path = "data/time.p"

def load_data():
    mothercells_data = pickle.load(open(mothercells_file_path, "rb"))
    daughtercells_data = pickle.load(open(daughtercells_data_path, "rb"))
    time_data = pickle.load(open(time_file_path, "rb"))
    return mothercells_data, daughtercells_data, time_data

def plot(data_tuple, subplot=True, show=True, legend=True):
    time = data_tuple[0]['time']
    var_names = data_tuple[0].keys()
    var_names.remove('time')
    rows_and_cols = np.ceil(np.sqrt(len(var_names)))
    for pos, variable in enumerate(var_names):
        if subplot:
            plt.subplot(rows_and_cols, rows_and_cols, pos + 1)
        for data_number, data_dict in enumerate(data_tuple):
            plt.plot(time, data_dict[variable], label=str(variable) + '_%s' % data_number)
    if legend:
        plt.legend()
    if show:
        plt.show()

def experimental_data_to_dict(time_vec, data_matrix):
    data_dict = dict()
    data_dict['time'] = time_vec
    for pos in range(len(data_matrix)):
        data_dict[pos] = data_matrix[pos]
    return data_dict

def simulation_to_dict(simulation_result):
    return {colname: simulation_result[colname] for colname in simulation_result.colnames}

def plot_experimental_data(mothercells_data, daughtercells_data, time_data, subplot=True, show=True):
    mothercells_dict = experimental_data_to_dict(time_data, mothercells_data)
    daughtercells_dict = experimental_data_to_dict(time_data, daughtercells_data)
    data_tuple = (mothercells_dict, daughtercells_dict)
    plot(data_tuple, subplot, show)

def load_model(model_path_antimony):
    model = te.loada(model_path_antimony)
    return model

def set_model_parameters(model, params):
    for param_id in params:
        model[param_id] = params[param_id]
    return model

def select_model_timecourses(model, time_course_selections):
    if 'time' in time_course_selections:
        model.timeCourseSelections = time_course_selections
    else:
        model.timeCourseSelections = ['time'] + time_course_selections
    return model

def simulate_model(model, end_time, steps=100):
    simulation_result = model.simulate(0, end_time, steps)
    return simulation_result

def compute_sqd_distance(simulation_result, data):
    dist = 0.
    for col_name in simulation_result.colnames:
        if col_name == 'time':
            continue
        dist += np.nansum((data[col_name] - simulation_result[col_name])**2)
    return dist

def truncate_data(data):

    """ truncate data in such a way that there are no leading or trailing nan values left"""
    pos_min = 0
    pos_max = np.inf
    for variable in data:
        data_vec_wo_nan = np.where(~np.isnan(data[variable]))
        start_pos = data_vec_wo_nan[0][0]
        end_pos = data_vec_wo_nan[0][-1]
        if start_pos > pos_min:
            pos_min = start_pos
        if end_pos < pos_max:
            pos_max = end_pos
    for variable in data:
        data[variable] = data[variable][pos_min:pos_max]
    return data

def time_vector_to_steps_and_stop(time_vector):
    assert time_vector[0] == 0
    stop = time_vector[-1]
    steps = len(time_vector)
    return steps, stop

def simulate_model_for_parameter_values(parameter_values, model, parameter_ids, time_vector):
    param_dict = dict(zip(parameter_ids, parameter_values))
    model.reset()
    model = set_model_parameters(model, param_dict)
    steps, end_time = time_vector_to_steps_and_stop(time_vector)
    simulation_result = simulate_model(model, end_time, steps)
    return simulation_result

def compute_objective_function(parameter_values, model, parameter_ids, data):
    simulation_result = simulate_model_for_parameter_values(parameter_values, model, parameter_ids, data['time'])
    return compute_sqd_distance(simulation_result, data)

def fit_model_to_data(model, data, parameters_to_fit, bounds={}):
    reference_params = {p_id: model[p_id] for p_id in parameters_to_fit}
    initial_params = [reference_params[p_id] for p_id in parameters_to_fit]
    model = select_model_timecourses(model, data.keys())
    additional_arguments = (model, parameters_to_fit, data)
    if bounds != {}:
        bounds_list = [bounds[p_id] for p_id in parameters_to_fit]
    else:
        bounds_list = [(0, 100*p_init) for p_init in initial_params]
    minimizer_kwargs ={'args': additional_arguments,
                        'bounds': bounds_list,
                        'method': 'L-BFGS-B' }
    return basinhopping(compute_objective_function, 
                        initial_params, 
                        minimizer_kwargs=minimizer_kwargs)


def plot_fitting_result_and_data(model, fitting_result, data, parameter_ids, subplot=True):
    parameter_values = fitting_result.x
    simulation_result = simulate_model_for_parameter_values(parameter_values, model, parameter_ids, data['time'])
    simulation_result_dict = simulation_to_dict(simulation_result)

    plot((simulation_result_dict, data))
    


if __name__ == '__main__':
    mothercells_data, daughtercells_data, time_data = load_data()
    #plot_single_cell(time_data, mothercells_data[0, :], daughtercells_data[0, :])
    #plot_data(mothercells_data, daughtercells_data, time_data)

    model = load_model('volume_reference.txt')
    #model = set_model_parameters(model, {'k_nutrient': 0})
    #model = select_model_timecourses(['V_tot_fl'])
    simulation_result = simulate_model(model, end_time=7200)


    time_starting_zero = np.array(time_data) + abs(min(time_data))
    data = {'time': time_starting_zero, 'V_tot_fl': mothercells_data[0, :]}
    parameter_ids = ['k_nutrient']
    fitting_result = fit_model_to_data(model, data, parameter_ids)
    plot_fitting_result_and_data(model, fitting_result, data, parameter_ids, subplot=True)
