#!/usr/bin/env python
import model_data
import simulate
from tools import OptimizationBounds
import cma
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
from scipy.optimize import basinhopping
plt.style.use('ggplot')

def compute_sqd_distance(simulation_result_dict, data):
    dist = 0.
    for variable in simulation_result_dict:
        if variable == 'time':
            continue
        dist += np.nansum((data[variable] - simulation_result_dict[variable])**2)
    return dist

def get_initial_volume_osmotic_from_data(model, data):
    initial_values = model_data.get_initial_values_from_data(data)
    assert 'V_tot_fl' in initial_values
    volume_osmotic = initial_values['V_tot_fl'] - model['V_b']
    parameters = {'V_os': volume_osmotic}
    return parameters

def compute_objective_function(parameter_values, model, parameter_ids, data, additional_model_parameters):
    print parameter_values
    try:
        simulation_result_dict = simulate.simulate_model_for_parameter_values(parameter_values, model, parameter_ids, data['time'], additional_model_parameters)
    except RuntimeError:
        print 'Error simulating model for parameters %s' %parameter_values
        return np.nan
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

def fit_cmaes(initial_params, additional_arguments, bounds_min, bounds_max, sigma0=4e-15, tolx=1e-25):
    options = cma.CMAOptions()
    #options.set('tolfun', 1e-14)
    options['tolx'] = tolx
    options['bounds'] = (bounds_min, bounds_max)
    param_vec = np.array(initial_params)
    p_min = max(param_vec.min(), 1e-20)
    options['scaling_of_variables'] = param_vec / p_min
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
    model.resetAll()
    reference_params = {p_id: model[p_id] for p_id in parameters_to_fit}
    initial_params = [reference_params[p_id] for p_id in parameters_to_fit]
    model = simulate.select_model_timecourses(model, data.keys())
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
                                 additional_model_parameters={},
                                 subplot=True, 
                                 show=True):
    simulation_result_dict = simulate.simulate_model_for_parameter_values(fitted_parameters, 
                                                                          model, 
                                                                          parameter_ids, 
                                                                          data['time'], 
                                                                          additional_model_parameters=additional_model_parameters)
    simulate.plot((simulation_result_dict, data), subplot=subplot, show=show)
    

def fit_model_to_all_cells(model, 
                           mothercells_data,
                           daughtercells_data,
                           time_data,
                           parameters_to_fit):
    fitting_results = []
    rows_and_cols = np.ceil(np.sqrt(len(mothercells_data)))
    for mother_pos in range(len(mothercells_data)):
        data = {'time': np.array(time_data), 'V_tot_fl': mothercells_data[mother_pos, :]}
        data = model_data.truncate_data(data)
        data['time'] = data['time'] + abs(min(data['time']))
        additional_model_parameters = get_initial_volume_osmotic_from_data(model, data)
        fitted_parameters = fit_model_to_data(model, 
                                              data, 
                                              parameters_to_fit, 
                                              'cmaes', 
                                              additional_model_parameters=additional_model_parameters)
        fitting_results.append(fitted_parameters)
        plt.subplot(rows_and_cols, rows_and_cols, mother_pos + 1)
        plot_fitting_result_and_data(model,
                                     fitted_parameters,
                                     data,
                                     parameters_to_fit,
                                     additional_model_parameters=additional_model_parameters,
                                     subplot=False,
                                     show=False)
    plt.show()
        
    return fitting_results


if __name__ == '__main__':
    mothercells_data, daughtercells_data, time_data = model_data.load_data()

    model = simulate.load_model('volume_reference_radius.txt')
    
    
    parameters_to_fit = ['k_nutrient', 'k_deg', 'r_os']
    
    data = {'time': np.array(time_data), 'V_tot_fl': mothercells_data[1, :]}
    data = model_data.truncate_data(data)
    data['time'] = data['time'] + abs(min(data['time']))

    additional_model_parameters = {'[c_i]': 319} #get_initial_volume_osmotic_from_data(model, data)
    fitted_parameters = fit_model_to_data(model, 
                                          data, 
                                          parameters_to_fit,   
                                          'cmaes', 
                                          additional_model_parameters=additional_model_parameters) 
    print fitted_parameters
    plot_fitting_result_and_data(model, fitted_parameters, data, parameters_to_fit, subplot=True, additional_model_parameters=additional_model_parameters)

    #fitted_parameters = fit_model_to_all_cells(model, mothercells_data, daughtercells_data, time_data, parameters_to_fit)
    #plot_fitting_result_and_data(model, fitted_parameters, mothercells_data, parameters_to_fit, subplot=True)
