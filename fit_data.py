#!/usr/bin/env python
import model_data
import simulate
from tools import OptimizationBounds
import cma
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
import pandas as pd
import pickle
from scipy.optimize import basinhopping

plt.style.use('ggplot')


def volume_to_radius(volume):
    r = (3./4/np.pi*volume) ** (1./3)
    return r   

def trafo(param_vec, direction):
  assert direction in [-1,1]
  if direction == 1:
    param_vec = np.log10(param_vec)
  else:
    param_vec.astype('float128')
    param_vec = 10**(param_vec)
    param_vec.astype('double')
  return param_vec

def compute_sqd_distance(simulation_result_dict, data, normalized=False):
    dist = 0.
    for variable in simulation_result_dict:
        if variable == 'time':
            continue
        if normalized:
          dist += np.nansum( ((data[variable] - simulation_result_dict[variable])**2)/data[variable]) / sum(~np.isnan(data[variable]))
        else:
          dist += np.nansum((data[variable] - simulation_result_dict[variable])**2) / sum(~np.isnan(data[variable]))
    return dist

def get_initial_volume_osmotic_from_data(model, data):
    initial_values = model_data.get_initial_values_from_data(data)
    assert 'V_tot_fl' in initial_values
    volume_osmotic = initial_values['V_tot_fl'] - model['V_b']
    parameters = {'V_os': volume_osmotic}
    return parameters

def compute_objective_function(parameter_values, model, parameter_ids, data, additional_model_parameters, additional_concentrations,log):
    
    if log:
      parameter_values = trafo(parameter_values,-1)
            
    try:
      simulation_result_dict = simulate.simulate_model_for_parameter_values(parameter_values, 
                                                                              model, 
                                                                              parameter_ids, 
                                                                              data['time'], 
                                                                              additional_model_parameters,
                                                                              additional_concentrations)
                                                                             
    except RuntimeError:
        print 'Error simulating model for parameters %s' %parameter_values
        return np.nan
    sqd = compute_sqd_distance(simulation_result_dict, data)

    #print('sqd:{0}'.format(sqd))
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

def fit_cmaes(initial_params,
              additional_arguments,
              bounds_min,
              bounds_max,
              sigma0=1e-17,
              tolx=1e-20,
              tolfun=1e-1,
              log=True):

    options = cma.CMAOptions()
    options.set('tolfun', tolfun)
    options['tolx'] = tolx
    #options['bounds'] = (bounds_min, bounds_max)

    param_vec = np.array(initial_params)
    
    #options['CMA_stds']=np.sqrt(param_vec)*1e5
    #options['CMA_stds']=param_vec
   
    if log:
      initial_params=trafo(param_vec,1)
      options['bounds'] = (trafo(bounds_min,1), trafo(bounds_max,1))
    else:   
      options['bounds'] = (bounds_min, bounds_max)
      p_min = max(param_vec.min(), 1e-20)
      options['scaling_of_variables'] = param_vec / p_min

    options['tolx'] = tolx  
    print('\n IP: {0}'.format(initial_params))
    
    result = cma.fmin(compute_objective_function, 
                      x0=initial_params,
                      sigma0=sigma0,
                      options=options,
                      args=additional_arguments)
    
    if log:
      results = trafo(np.array(result[0]),-1)
    else:
      results =result[0]
    
    return results#result[0]

def define_bounds(parameters_to_fit, reference_params, bounds={}):
    
    xmin=[0]*len(parameters_to_fit)
    xmax=[0]*len(parameters_to_fit)
    
    for k,p_id in enumerate(parameters_to_fit):

      if p_id in bounds.keys():
        xmin[k] = bounds[p_id][0] 
        xmax[k] = bounds[p_id][1]
      else:
        xmin[k] = reference_params[p_id]*0.05 
        xmax[k] = reference_params[p_id]*50  

    print('bounds: {0}'.format(bounds))    
    print('fitted parameters: {0}'.format(parameters_to_fit))    
    print('xmin : {0}, xmax: {1}'.format(xmin,xmax))

    return np.array(xmin), np.array(xmax)

def get_initial_parameter_from_data(model,data,params):
  para_ini = {}
  assert data['time'][0] == 0, 'data does not start at time 0'
  
  if 'budding_start' in params:
    para_ini['budding_start'] = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]

  if 'mother_r_os_0' in params:
    volume_mother_0 = data['mother_V_tot_fl'][0]
    r_tot_mother_0 = volume_to_radius(volume_mother_0)
    para_ini['mother_r_os_0'] = r_tot_mother_0 - model['mother_r_b']

  if 'bud_r_os_0' in params:  
    volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
    r_tot_bud_0 = volume_to_radius(volume_bud_0)
    para_ini['bud_r_os_0'] = r_tot_bud_0 - model['bud_r_b']

  return para_ini


def fit_model_to_data(model,
                      data,
                      parameters_to_fit,
                      optimizer,
                      bounds={},
                      additional_model_parameters={},
                      additional_concentrations={},
                      params_ini={},
                      log=True):
    
    model.resetAll()
    reference_params = {p_id: model[p_id] for p_id in parameters_to_fit}
    
    if params_ini!={}:
      for p_id in params_ini:
        reference_params[p_id] = params_ini[p_id]

    initial_params = [reference_params[p_id] for p_id in parameters_to_fit]
    model_output = simulate.select_model_timecourses(model, data.keys())
    additional_arguments = (model_output,
                            parameters_to_fit, 
                            data,
                            additional_model_parameters,
                            additional_concentrations,
                            log)
     
    xmin, xmax = define_bounds(parameters_to_fit, reference_params, bounds=bounds)    
    
    if log:
      sigma0=1e-1
      #tolx = 5e-3
    else:  
      sigma0=4e-17


    if optimizer == 'basin':
        return fit_basin(initial_params, additional_arguments, xmin, xmax)
    elif optimizer == 'cmaes':
        return fit_cmaes(initial_params, additional_arguments, xmin, xmax, sigma0=sigma0, log=log)
    else:
        raise Exception('unknown optimization method')

def truncate_initial_budsize(simulation_result_dict, data): # NaN values in data result in NaN values of fit
  n = np.nan
  idx_first_notnan = np.argwhere(~np.isnan(data['bud_V_tot_fl'])>0)[0][0]
  l = np.array([n for x in range(idx_first_notnan)])
  trunc_ini_bud=np.concatenate([l,simulation_result_dict['bud_V_tot_fl'][idx_first_notnan:]])
  simulation_result_dict['bud_V_tot_fl']  = trunc_ini_bud
  
  return simulation_result_dict
  
def expand_data_with_nan(data,observables):
  l=len(data['time'])
  nan_array=np.array([np.nan]*l)
  nan_dic =  {observables[i] : nan_array for i in range(len(observables))}
  data.update(nan_dic)
  return data


def plot_fitting_result_and_data(model, 
                                 fitted_parameters, 
                                 data, 
                                 parameter_ids, 
                                 additional_model_parameters={},
                                 additional_concentrations={},
                                 subplot=True, 
                                 show=True,
                                 legend=True,
                                 observables=[]):
    
    observables.extend(data.keys())
    print(observables)
    model = simulate.select_model_timecourses(model, observables)
    simulation_result_dict = simulate.simulate_model_for_parameter_values(fitted_parameters, 
                                                                          model, 
                                                                          parameter_ids, 
                                                                          data['time'], 
                                                                          additional_model_parameters=additional_model_parameters,
                                                                          additional_concentrations=additional_concentrations)
    if 0: # delete the initial budsize in the plot 
      simulation_result_dict = truncate_initial_budsize(simulation_result_dict, data)
    simulate.plot((simulation_result_dict, data),
                                 subplot=subplot,
                                       show=show,
                                   legend=legend,
                                    time_scale='min')


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
        data = model_data.limit_time_in_data(data, max_time=300)
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

'''def get_additional_model_parameters(data):
    assert data['time'][0] == 0
    #budding_start = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]
    volume_mother_0 = data['mother_V_tot_fl'][0]
    volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
    r_tot_mother_0 = volume_to_radius(volume_mother_0)
    r_tot_bud_0 = volume_to_radius(volume_bud_0)
    r_os_mother_0 = r_tot_mother_0 - model['mother_r_b']
    r_os_bud_0 = r_tot_bud_0 - model['bud_r_b']   
    additional_model_parameters = { 'mother_r_os_0': r_os_mother_0}#,
                                    #'bud_r_os': r_os_bud,                                   
                                    #'budding_start': budding_start }  
    print(additional_model_parameters)                                   
    return additional_model_parameters
'''


if __name__ == '__main__':
  max_time=500*60
  cellID=5
  #data
  mothercells_data, daughtercells_data, time_data_min = model_data.load_data()

  # convert time to seconds
  time_data = [x*60 for x in time_data_min]
  #data dict
  data = {'time': np.array(time_data), 'mother_V_tot_fl': mothercells_data[cellID, :],'bud_V_tot_fl': daughtercells_data[cellID, :]}
  
      # 1 parameter in log scale
  if 1:
    log=True
  else:
    log=False
  
  #model = simulate.load_model('volume_reference_radius.txt')
  #parameters_to_fit = ['k_nutrient', 'k_deg', 'r_os']

  model = simulate.load_model('volume_mother_and_bud.txt')
  parameters_to_fit = ['budding_start','k_nutrient', 'k_deg_0', 'mother_phi', 'bud_phi','mother_r_os_0']

  additional_concentrations = {'init([mother_c_i])': 319.17,
                              'init([bud_c_i])': 319.17 }
  #additional_concentrations = {'init([mother_c_i_0])': 325,
  #                                'init([bud_c_i_0])': 325 }

  #data = {'time': np.array(time_data), 'V_tot_fl': mothercells_data[1, :]}
  
  additional_model_parameters={}
  data_trunc = model_data.truncate_data(data)
  data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
  data_trunc = model_data.limit_time_in_data(data_trunc, max_time=max_time)
  #additional_model_parameters = get_additional_model_parameters(data_trunc)

  para_ini_ids = ['budding_start', 'mother_r_os_0']
  params_ini=get_initial_parameter_from_data(model,data_trunc,para_ini_ids)

  #bud_start_data = get_initial_bud_start_guess(data_trunc)
  bud_start_data = params_ini['budding_start']
  start_tolerance = 120*60

  bounds={}
  bounds['budding_start'] = [bud_start_data - start_tolerance, bud_start_data + start_tolerance]
  bounds['mother_r_os_0'] = [params_ini['mother_r_os_0']-0.3,params_ini['mother_r_os_0']+0.3] 
  #bounds['mother_phi']= [5.e-6, 5.e-1]
  #bounds['bud_phi']= [1.e-5, 5.e-1 ]
  #bounds['k_nutrient_']= [9.e-16, 9.e-17]
  #bounds['k_deg_0']= [1.e-15, 1.e-17]
  
  
  additional_model_parameters={'mother_bud_water_perm':1}
  
  if 1:
    fitted_parameters = fit_model_to_data(model, 
                                          data, 
                                          parameters_to_fit,   
                                          'cmaes',
                                          bounds=bounds,
                                          additional_concentrations=additional_concentrations,
                                          params_ini=params_ini,
                                          log=log,
                                          additional_model_parameters=additional_model_parameters) 
    
    with open('fitted_parameters_cellID_'+str(cellID)+'.p','w') as f:
      pickle.dump(fitted_parameters,f)
    
  else:
    with open('fitted_parameters_cellID_'+str(cellID)+'.p','rb') as f:
      fitted_parameters=pickle.load(f)

  if 1:
    print('parameters: {0} \n values: {1}'.format(parameters_to_fit,fitted_parameters))
    observables=['mother_R_ref','mother_r_os','mother_pi_t','[mother_c_i]']
    data=expand_data_with_nan(data,observables)
    
    plot_fitting_result_and_data(model,
                            fitted_parameters,
                            data, parameters_to_fit,
                            subplot=True,
                            observables=observables)#,
                            #additional_model_parameters=additional_model_parameters)

#fitted_parameters = fit_model_to_all_cells(model, mothercells_data, daughtercells_data, time_data, parameters_to_fit)
#plot_fitting_result_and_data(model, fitted_parameters, mothercells_data, parameters_to_fit, subplot=True)


