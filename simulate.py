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

def simulate_model(model, end_time, steps=100):
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

def simulate_model_for_parameter_values(parameter_values, model, parameter_ids, time_vector, additional_model_parameters={}, initial_assignments={}):
    param_dict = dict(zip(parameter_ids, parameter_values))
    model.resetAll()
    if initial_assignments == {}:
        initial_assignments = tools.get_initial_assignments_dict(model)
    model = model_data.set_model_parameters(model, additional_model_parameters)
    model = model_data.set_model_parameters(model, param_dict)
    model = evaluate_initial_assignments(model, initial_assignments)
    steps, end_time = time_vector_to_steps_and_stop(time_vector)
    simulation_result_dict = simulate_model(model, end_time, steps)
    return simulation_result_dict


if __name__ == '__main__':
    import sys
    model = load_model(sys.argv[1])
    #model = select_model_timecourses(model, ['bud_V_b'])
    simulation_result = simulate_model(model, end_time=7200)
    plot((simulation_result,), legend=True)


