#!/usr/bin/env python

import model_data
import fit_data
import tools
import math
import tellurium as te
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

math_functions = {'pi': math.pi}


def plot(data_tuple, subplot=True, show=True, legend=True):
    time = data_tuple[0]['time']
    var_names = data_tuple[0].keys()
    var_names.remove('time')
    rows_and_cols = np.ceil(np.sqrt(len(var_names)))
    for pos, variable in enumerate(var_names):
        if subplot:
            plt.subplot(rows_and_cols, rows_and_cols, pos + 1)
            if legend:
                plt.title(str(variable))
        for data_number, data_dict in enumerate(data_tuple):
            plt.plot(time, data_dict[variable], label=str(variable) + '_%s' % data_number)
            plt.xlim(-50, time[-1])
    if legend:
        plt.legend()
    if show:
        plt.show()

def simulation_to_dict(simulation_result):
    return {colname: simulation_result[colname] for colname in simulation_result.colnames}

def load_model(model_path_antimony):
    model = te.loada(model_path_antimony)
    return model

def select_model_timecourses(model, time_course_selections):
    if 'time' in time_course_selections:
        model.timeCourseSelections = time_course_selections
    else:
        model.timeCourseSelections = ['time'] + time_course_selections
    return model

def set_integrator_options(model):
    model.integrator.relative_tolerance = 1e-10
    #model.integrator.absolute_tolerance = 1e-16
    #model.integrator.variable_step_size = True
    return model

def simulate_model(model, end_time, steps=100):
    model = set_integrator_options(model)
    simulation_result = model.simulate(0, end_time, steps)
    return simulation_to_dict(simulation_result)

def evaluate_initial_assignments(model, initial_assignments):
    param_dict = model_data.get_model_parameters_as_dict(model)
    for variable in initial_assignments:
        model[variable] = eval(initial_assignments[variable], param_dict, math_functions)

    return model

def time_vector_to_steps_and_stop(time_vector):
    assert time_vector[0] == 0
    stop = time_vector[-1]
    steps = len(time_vector)
    return steps, stop

def simulate_model_for_parameter_values(parameter_values, model, parameter_ids, time_vector, additional_model_parameters={}, additional_concentrations={}, initial_assignments={}):
    param_dict = dict(zip(parameter_ids, parameter_values))
    model.resetAll()
    if initial_assignments == {}:
        initial_assignments = tools.get_initial_assignments_dict(model)
    model = model_data.set_model_parameters(model, param_dict)    
    model = model_data.set_model_parameters(model, additional_model_parameters)
    model = model_data.set_model_parameters(model, additional_concentrations)
    model = evaluate_initial_assignments(model, initial_assignments)
    steps, end_time = time_vector_to_steps_and_stop(time_vector)

    #print model['']#
    simulation_result_dict = simulate_model(model, end_time, steps)
    return simulation_result_dict


if __name__ == '__main__':
    import sys
    model = load_model(sys.argv[1])
    
    #p_values = [  1e13,   1e13,   2.04684012e+00]
    #parameters_to_fit = ['k_nutrient', 'k_deg', 'mother_r_os']
    p_values = []
    parameters_to_fit = []

    p_values = [2.576782e-14,  3.124156e-14,  4.498614e-04,  1.783651e-02 ]

    parameters_to_fit = ['k_nutrient', 'k_deg', 'mother_phi', 'bud_phi'] 


    #model = select_model_timecourses()
    additional_model_parameters = { 'budding_start': 126}
    additional_model_parameters = {'mother_r_os': 1.239061741911077, 'budding_start': 129, 'bud_r_os': 1.2239599669205963}
    #additional_model_parameters =  { 'budding_start': 126,
    #                                'mother_r_os': 0.97506339013385745,
    #                                'bud_r_os': 1.2142117908125938} 
    #additional_model_parameters = {'[c_i]': 325,
    #                                'r_os': 10}

    additional_concentrations = {'[mother_c_i]': 325,
                                    '[bud_c_i]': 325 }
    


    #additional_model_parameters = {}
    simulation_result = simulate_model_for_parameter_values( p_values, model, parameters_to_fit,  np.linspace(0,400, 10000), additional_model_parameters=additional_model_parameters,
    additional_concentrations=additional_concentrations )
    #simulation_result = simulate_model(model, end_time=500)
    plot((simulation_result,), legend=True)


