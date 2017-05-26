#!/usr/bin/env python
from tools import OptimizationBounds
import cma
import pickle
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
from scipy.optimize import basinhopping
plt.style.use('ggplot')

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
    return simulation_to_dict(simulation_result)

def compute_sqd_distance(simulation_result_dict, data):
    dist = 0.
    for variable in simulation_result_dict:
        if variable == 'time':
            continue
        dist += np.nansum((data[variable] - simulation_result_dict[variable])**2)
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

def get_initial_values_from_data(data):
    initial_values = {}
    for variable in data:
        if variable == 'time':
            continue
        initial_values[variable] = data[variable][0]
    return initial_values

def get_initial_volume_osmotic_from_data(model, data):
    initial_values = get_initial_values_from_data(data)
    assert 'V_tot_fl' in initial_values
    volume_osmotic = initial_values['V_tot_fl'] - model['V_b']
    parameters = {'V_os': volume_osmotic}
    return parameters

def simulate_model_for_parameter_values(parameter_values, model, parameter_ids, time_vector, additional_model_parameters={}):
    param_dict = dict(zip(parameter_ids, parameter_values))
    model.reset()
    model = set_model_parameters(model, additional_model_parameters)
    model = set_model_parameters(model, param_dict)
    steps, end_time = time_vector_to_steps_and_stop(time_vector)
    simulation_result_dict = simulate_model(model, end_time, steps)
    return simulation_result_dict

def compute_objective_function(parameter_values, model, parameter_ids, data, additional_model_parameters):
    #print parameter_values
    simulation_result_dict = simulate_model_for_parameter_values(parameter_values, model, parameter_ids, data['time'], additional_model_parameters)
    sqd = compute_sqd_distance(simulation_result_dict, data)
    return compute_sqd_distance(simulation_result_dict, data)

def fit_basin(initial_params, additional_arguments, bounds_min, bounds_max, basin_iterations=1e3, basin_stepsize=1e-15):
    bounds_basinhopping = OptimizationBounds(xmin=bounds_min, xmax=bounds_max)
    bounds_lbfgsb = zip(bounds_min, bounds_max)
    minimizer_kwargs ={'args': additional_arguments,
                       'method': 'L-BFGS-B',
                       'bounds': bounds_lbfgsb }
    result = basinhopping(compute_objective_function, 
                        initial_params, 
                        niter=int(basin_iterations),
                        stepsize=basin_stepsize,
                        minimizer_kwargs=minimizer_kwargs,
                        accept_test=bounds_basinhopping)
    return result.x

def fit_cmaes(initial_params, additional_arguments, bounds_min, bounds_max, sigma0=4e-16, tolx=1e-17):
    options = cma.CMAOptions()
    #options.set('tolfun', 1e-14)
    options['tolx'] = tolx
    options['bounds'] = (bounds_min, bounds_max)

    result = cma.fmin(compute_objective_function, 
                      x0=initial_params,
                      sigma0=sigma0,
                      options=options,
                      args=additional_arguments)
    return result[0]

def fit_model_to_data(model,
                      data,
                      parameters_to_fit,
                      optimizer,
                      bounds={},
                      additional_model_parameters={}):
    reference_params = {p_id: model[p_id] for p_id in parameters_to_fit}
    initial_params = [reference_params[p_id] for p_id in parameters_to_fit]
    model = select_model_timecourses(model, data.keys())
    additional_arguments = (model, parameters_to_fit, data, additional_model_parameters)
    if bounds != {}:
        xmin = [bounds[p_id][0] for p_id in parameters_to_fit]
        xmax = [bounds[p_id][1] for p_id in parameters_to_fit]
        
    else:
        xmin = [0.] * len(initial_params)
        xmax = 100. * np.array(initial_params)
    if optimizer == 'basin':
        return fit_basin(initial_params, additional_arguments, xmin, xmax)
    elif optimizer == 'cmaes':
        return fit_cmaes(initial_params, additional_arguments, xmin, xmax)
    else:
        raise Exception('unknown optimization method')

def plot_fitting_result_and_data(model, 
                                 fitted_parameters, 
                                 data, 
                                 parameter_ids, 
                                 subplot=True, 
                                 additional_model_parameters={}):
    simulation_result_dict = simulate_model_for_parameter_values(fitted_parameters, 
                                                            model, 
                                                            parameter_ids, 
                                                            data['time'], 
                                                            additional_model_parameters=additional_model_parameters)
    plot((simulation_result_dict, data))
    


if __name__ == '__main__':
    mothercells_data, daughtercells_data, time_data = load_data()
    #plot_single_cell(time_data, mothercells_data[0, :], daughtercells_data[0, :])
    #plot_data(mothercells_data, daughtercells_data, time_data)

    model = load_model('volume_reference.txt')
    #model = set_model_parameters(model, {'k_nutrient': 0})
    #model = select_model_timecourses(['V_tot_fl'])
    simulation_result = simulate_model(model, end_time=7200)


    data = {'time': np.array(time_data), 'V_tot_fl': mothercells_data[1, :]}
    data = truncate_data(data)
    data['time'] = data['time'] + abs(min(data['time']))
    parameter_ids = ['k_nutrient', 'k_deg']
    additional_model_parameters = get_initial_volume_osmotic_from_data(model, data)
    fitted_parameters = fit_model_to_data(model, data, parameter_ids, 'cmaes', additional_model_parameters=additional_model_parameters)
    print fitted_parameters
    plot_fitting_result_and_data(model, fitted_parameters, data, parameter_ids, subplot=True, additional_model_parameters=additional_model_parameters)
