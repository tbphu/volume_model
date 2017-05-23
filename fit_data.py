#!/usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
from scipy.optimize import basinhopping


mothercells_file_path = 'mothercells.p'
daughtercells_data_path = 'daughtercells.p'
time_file_path = "time.p"

def load_data():
	mothercells_data = pickle.load(open(mothercells_file_path, "rb"))
	daughtercells_data = pickle.load(open(daughtercells_data_path, "rb"))
	time_data = pickle.load(open(time_file_path, "rb"))
	return mothercells_data, daughtercells_data, time_data

def plot_single_cell(time_vec, mother_vec, daughter_vec=None):
		
	plt.plot(time_vec, mother_vec)
	if daughter_vec is not None:
		plt.plot(time_vec, daughter_vec)

def plot_data(mothercells_data, daughtercells_data, time_data, subplot=True):
	rows_and_cols = np.ceil(np.sqrt(len(mothercells_data)))
	for cell_id in range(len(mothercells_data)):
		if subplot:
			plt.subplot(rows_and_cols, rows_and_cols, cell_id + 1)
		plot_single_cell(time_data, mothercells_data[cell_id], daughtercells_data[cell_id])
	plt.show()

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

def plot_simulation_result(simulation_result, subplot=True, show=True):
	time = simulation_result['time']
	rows_and_cols = np.ceil(np.sqrt(len(simulation_result.colnames)))
	for pos, variable in enumerate(simulation_result.colnames):
		if subplot:
			plt.subplot(rows_and_cols, rows_and_cols, pos + 1)
			plt.title(variable)
		plt.plot(time, simulation_result[variable], label=variable)
		plt.legend()
	if show:	
		plt.show()

def compute_sqd_distance(simulation_result, data):
	dist = 0.
	for col_name in simulation_result.colnames:
		if col_name == 'time':
			continue
		dist += np.nansum((data[col_name] - simulation_result[col_name])**2)
	return dist

def simulate_model_with_changed_parameters(parameter_values, model, parameter_ids, data):
	param_dict = dict(zip(parameter_ids, parameter_values))
	model.reset()
	model = set_model_parameters(model, param_dict)
	end_time = data['time'][-1]
	steps = len(data['time'])
	simulation_result = simulate_model(model, end_time, steps)
	return simulation_result

def simulate_objective_function(parameter_values, model, parameter_ids, data):
	simulation_result = simulate_model_with_changed_parameters(parameter_values, model, parameter_ids, data)
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
						"method": "L-BFGS-B" }
	return basinhopping(simulate_objective_function, 
						initial_params, 
						minimizer_kwargs=minimizer_kwargs)


def plot_model_simulation_for_fitting_result(model, fitting_result, data, parameter_ids, subplot=True):
	parameter_values=fitting_result.x
	model_result=simulate_model_with_changed_parameters(parameter_values, model, parameter_ids, data)
	rows_and_cols = np.ceil(np.sqrt(len(data)-1))
	pos=1
	for variable in data:
		if variable == 'time':
			continue
		if subplot: 
			plt.subplot(rows_and_cols,rows_and_cols, pos)	
		plot_single_cell(data['time'], data[variable])
		plt.plot(model_result['time'], model_result[variable])
		pos+=1
	plt.show()	
	


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
	plot_model_simulation_for_fitting_result(model, fitting_result, data, parameter_ids, subplot=True)
