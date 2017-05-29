#!/usr/bin/env python

import simulate
import fit_data
import model_data
import pandas as pd
import tools

model = simulate.load_model('volume_reference.txt')


ia = tools.get_initial_assignments_dict(model)
fit_data.evaluate_initial_assignments(model, ia)
#df =  model_data.get_model_parameters_as_dataframe(model)
#df.to_csv('used_model_parameters.csv')

#set paras
#param_list = ['V_b', 'V_os', '[c_e]', '[c_i]', 	'k_deg', 'k_nutrient', 'modulus_adjustment', 'phi',	'pi_t0', 'pi_tc_0']

#df = pd.read_csv('used_model_parameters.csv')
#model = model_data.set_model_parameters_from_dataframe(model, df, param_list, row=1)



#run simu

#simulation_result = simulate.simulate_model(model, end_time=7200)
#simulate.plot((simulation_result,), legend=True)

