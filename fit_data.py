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
from multiprocessing import Pool
import functools as ft 
from scipy.optimize import basinhopping
from matplotlib import cm

#plt.style.use('seaborn-colorblind')



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

def compute_sqd_distance(simulation_result_dict, data, normalized=True):
    dist = 0.
    for variable in simulation_result_dict:
        if variable == 'time':
            continue
        if normalized:
          dist += np.nansum( ((data[variable] - simulation_result_dict[variable])**2)/data[variable]) #/ sum(~np.isnan(data[variable]))
        else:
          dist += np.nansum((data[variable] - simulation_result_dict[variable])**2) #/ sum(~np.isnan(data[variable]))
    return dist

def get_initial_volume_osmotic_from_data(model, data):
    initial_values = model_data.get_initial_values_from_data(data)
    assert 'V_tot_fl' in initial_values
    volume_osmotic = initial_values['V_tot_fl'] - model['V_b']
    parameters = {'V_os': volume_osmotic}
    return parameters

def compute_objective_function(parameter_values, model_output, parameter_ids, data, additional_model_parameters, additional_concentrations,log):
    
    if log:
      parameter_values = trafo(parameter_values,-1)
            
    try:
      simulation_result_dict = simulate.simulate_model_for_parameter_values(parameter_values, 
                                                                              model_output, 
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


    #opt_data = cma.CMADataLogger()

    options = cma.CMAOptions()
    options.set('tolfun', tolfun)
    options['tolx'] = tolx
    #options['verb_plot'] = 10
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
      results = result[0]

    #print(results)  
    
    return results, result #result[0]

def define_bounds(parameters_to_fit, reference_params, bounds={}, min_max_bounds_fact=(0.02,10)):
    
    xmin=[0]*len(parameters_to_fit)
    xmax=[0]*len(parameters_to_fit)
    
    for k,p_id in enumerate(parameters_to_fit):

      if p_id in bounds.keys():
        xmin[k] = bounds[p_id][0] 
        xmax[k] = bounds[p_id][1]
      else:
        xmin[k] = reference_params[p_id] * min_max_bounds_fact[0] 
        xmax[k] = reference_params[p_id] * min_max_bounds_fact[1]

    print('bounds: {0}'.format(bounds))    
    print('fitted parameters: {0}'.format(parameters_to_fit))    
    print('xmin : {0}, xmax: {1}'.format(xmin,xmax))

    return np.array(xmin), np.array(xmax)

def get_initial_parameter_from_data(model,data,params, params_ini={}):
  
  assert data['time'][0] == 0, 'data does not start at time 0'
  
  if ('budding_start' in params) & ('budding_start' not in params_ini.keys()):
    params_ini['budding_start'] = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]

  if ('mother_r_os_0' in params) & ('mother_r_os_0' not in params_ini.keys()):
    volume_mother_0 = data['mother_V_tot_fl'][0]
    r_tot_mother_0 = volume_to_radius(volume_mother_0)
    params_ini['mother_r_os_0'] = r_tot_mother_0 - model['mother_r_b']

  if ('bud_r_os_0' in params) & ('bud_r_os_0' not in params_ini.keys()):  
    volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
    r_tot_bud_0 = volume_to_radius(volume_bud_0)
    params_ini['bud_r_os_0'] = r_tot_bud_0 - model['bud_r_b']

  return params_ini


def fit_model_to_data(model,
                      data,
                      parameters_to_fit,
                      optimizer,
                      bounds={},
                      additional_model_parameters={},
                      additional_concentrations={},
                      params_ini={},
                      log=True,
                      tolfun=0.1,
                      sigma0_log_notlog = (0.5e-1,4e-17),
                      min_max_bounds_fact = (0.02,10)):
    
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
     
    xmin, xmax = define_bounds(parameters_to_fit, reference_params, bounds=bounds, min_max_bounds_fact=min_max_bounds_fact)    
    
    if log:
      sigma0=sigma0_log_notlog[0]
      #tolx = 5e-3
    else:  
      sigma0=sigma0_log_notlog[1]


    if optimizer == 'basin':
        return fit_basin(initial_params, additional_arguments, xmin, xmax)
    elif optimizer == 'cmaes':
        return fit_cmaes(initial_params, additional_arguments, xmin, xmax, sigma0=sigma0, log=log, tolfun=tolfun)
    else:
        raise Exception('unknown optimization method')

def truncate_initial_budsize(simulation_result_dict, data): # NaN values in data result in NaN values of fit
  n = np.nan
  idx_first_notnan = np.argwhere(~np.isnan(data['bud_V_tot_fl'])>0)[0][0]
  l = np.array([n for x in range(idx_first_notnan)])
  trunc_ini_bud=np.concatenate([l,simulation_result_dict['bud_V_tot_fl'][idx_first_notnan:]])
  simulation_result_dict['bud_V_tot_fl']  = trunc_ini_bud
  
  return simulation_result_dict
  
def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z


def expand_data_with_nan(data,observables):
  l=len(data['time'])
  nan_array = np.array([np.nan]*l)
  nan_dic = {observables[i] : nan_array for i in range(len(observables))}
  new_data = merge_two_dicts(data, nan_dic) 
  return new_data


def plot_fitting_result_and_data(model, 
                                 fitted_parameters, 
                                 data, 
                                 parameter_ids, 
                                 additional_model_parameters={},
                                 additional_concentrations={},
                                 subplot=True, 
                                 show=True,
                                 legend=True,
                                 observables=[],
                                 only_data=False):
    obs = list(observables)
    obs.extend(data.keys())
    #print(obs)

    model = simulate.select_model_timecourses(model, obs)
    
    simulation_result_dict = simulate.simulate_model_for_parameter_values(fitted_parameters, 
                                                                          model, 
                                                                          parameter_ids, 
                                                                          data['time'], 
                                                                          additional_model_parameters=additional_model_parameters,
                                                                          additional_concentrations=additional_concentrations)
    if 0: # delete the initial budsize in the plot 
      simulation_result_dict = truncate_initial_budsize(simulation_result_dict, data)
    '''simulate.plot((simulation_result_dict, data),
                                             subplot=subplot,
                                                   show=show,
                                               legend=legend,
                                                time_scale='min')'''
    if only_data:
         simulate.plot((data,),
                                  subplot=subplot,
                                        show=show,
                                    legend=legend,
                                     time_scale='min',
                                     only_data=only_data)
    else:  
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
        fitted_parameters, opt_res = fit_model_to_data(model, 
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


def get_data_trunc(data,max_time):
    data_trunc = model_data.truncate_data(data)
    data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
    data_trunc = model_data.limit_time_in_data(data_trunc, max_time=max_time)
    return data_trunc

def get_bounds_from_initial_params(params_ini, budstart_tol = 120*60):
    bud_start_data = params_ini['budding_start']
    bounds={}
    bounds['budding_start'] = [bud_start_data - budstart_tol, bud_start_data + budstart_tol]
    bounds['mother_r_os_0'] = [params_ini['mother_r_os_0']-0.3,params_ini['mother_r_os_0']+0.3]
    return bounds


def fit_single_cell(cellID, plotting=False, min_max_bounds_fact=(0.02,10), just_fitting_paras=True):
  #data dict
  data = {'time': np.array(time_data), 'mother_V_tot_fl': mothercells_data[cellID, :],'bud_V_tot_fl': daughtercells_data[cellID, :]}
  data_trunc = get_data_trunc(data, max_time)

  #additional_model_parameters = get_additional_model_parameters(data_trunc)
  para_ini_ids = ['budding_start', 'mother_r_os_0']
  params_ini = get_initial_parameter_from_data(model,data_trunc,para_ini_ids)
  bounds = get_bounds_from_initial_params(params_ini, budstart_tol = 120*60)

  try:
    fitted_parameters, opt_result = fit_model_to_data(model, 
                                                      data_trunc, #data
                                                      parameters_to_fit,   
                                                      'cmaes',
                                                      bounds=bounds,
                                                      additional_concentrations=additional_concentrations,
                                                      params_ini=params_ini,
                                                      log=log,
                                                      additional_model_parameters=additional_model_parameters,
                                                      tolfun=tolfun,
                                                      sigma0_log_notlog=sigma0_log_notlog,
                                                      min_max_bounds_fact = min_max_bounds_fact) 
    #print(fitted_parameters)
    model_output = simulate.select_model_timecourses(model, data.keys())
    
    msd = compute_objective_function(fitted_parameters,
                                     model_output,
                                     parameters_to_fit,
                                     data_trunc,
                                     additional_model_parameters,
                                     additional_concentrations,
                                     False)
    
    data_trunc = expand_data_with_nan(data_trunc,observables)
    fitted_values = fitted_parameters.tolist() + [msd] + [cellID]
    
  except:
    fitted_values = [np.nan]*(len(parameters_to_fit)+2)
  
  if plotting:
    plot_fitting_result_and_data(model,
                                  fitted_parameters,
                                  data_trunc,
                                  parameters_to_fit,
                                  subplot=True,
                                  observables=observables,
                                  additional_model_parameters=additional_model_parameters,
                                  additional_concentrations=additional_concentrations) 
  colnames = list(parameters_to_fit)
  colnames.extend(['MSD', 'CellID']) 
  fitted_paras_dict = dict(zip(colnames,fitted_values))

  if just_fitting_paras == True:
     return  fitted_paras_dict
  else:
     return  opt_result[-1].load().data()    


def plot_liklihood(df_dict, para_to_test):
  
  df_test = df_dict[df_dict.keys()[0]]  
  

  exlcluded_para_names = ['CellID', 'MSD']
  para_names = [x for x in df_test.keys() if x not in exlcluded_para_names]

  k = len(para_names)

  
  for value in df_dict.keys():
    df = df_dict_2[maxF]

    for  i, para in enumerate(para_names): 
      plt.figure(1) 
      ax = plt.subplot(k/2+1,2,i+1)
      y = df['MSD'].mean()
      x = df[para].mean()
      
      ax.plot(x,y,'.')
      ax.set_ylabel(r'$\chi ^2$')
      ax.set_xlabel(para)
  
  plt.tight_layout()
  plt.show()
  plt.savefig('liklihood.png')



if __name__ == '__main__':
  import sys
  
  max_time = 450*60
  #cellID = 5
  
  try:
    cellID = sys.argv[1]
  except IndexError:  
    cellID = 3
  
  #data
  mothercells_data, daughtercells_data, time_data_min = model_data.load_data()

  # convert time to seconds
  time_data = [x*60 for x in time_data_min]
  
  # 1 parameter in log scale
  if 1:
    log = True
  else:
    log = False

  tolfun = 1e-4  
  sigma0_log_notlog = (1e-1, 4e-17)
  min_max_bounds_fact = (0.02, 10000)
  
  #model = simulate.load_model('volume_reference_radius.txt')
  #parameters_to_fit = ['k_nutrient', 'k_deg', 'r_os']

  model = simulate.load_model('volume_mother_and_bud.txt')
  parameters_to_fit = ['budding_start','k_nutrient', 'k_scaling_factor', 'mother_phi', 'extens_factor','mother_r_os_0']

  # global model parameters: 
  additional_concentrations = {'init([mother_c_i])': 319.4,
                              'init([bud_c_i])': 319.4 }
  
  additional_model_parameters = {'mother_bud_water_perm':1,
                                'withSF':1}
  #extra plotted model values
  observables = ['mother_R_ref','mother_r_os','mother_pi_t','[mother_c_i]']                              
  #observables=[]
  
  #data dict
  data = {'time': np.array(time_data), 'mother_V_tot_fl': mothercells_data[cellID, :],'bud_V_tot_fl': daughtercells_data[cellID, :]}
  data_trunc = get_data_trunc(data,max_time)

  #additional_model_parameters = get_additional_model_parameters(data_trunc)
  para_ini_ids = ['budding_start', 'mother_r_os_0']
  params_ini = get_initial_parameter_from_data(model,data_trunc,para_ini_ids)
  bounds = get_bounds_from_initial_params(params_ini, budstart_tol=120*60)

  
  
  if 0: # fit one cell with defined CellID
    fitted_parameters, opt_result = fit_model_to_data(model, 
                                          data_trunc, #data
                                          parameters_to_fit,   
                                          'cmaes',
                                          bounds=bounds,
                                          additional_concentrations=additional_concentrations,
                                          params_ini=params_ini,
                                          log=log,
                                          additional_model_parameters=additional_model_parameters,
                                          tolfun=tolfun,
                                          sigma0_log_notlog=sigma0_log_notlog,
                                          min_max_bounds_fact=min_max_bounds_fact) 
    print(fitted_parameters)
    
    with open('fitted_parameters_cellID_'+str(cellID)+'.p','w') as f:
      pickle.dump(fitted_parameters,f)
    with open('optimization_result.p','w') as f:
      pickle.dump(opt_result,f)   
    
  else:
    with open('fitted_parameters_cellID_'+str(cellID)+'.p','rb') as f:
      fitted_parameters=pickle.load(f)

  if 0:
    print('parameters: {0} \n values: {1}'.format(parameters_to_fit,fitted_parameters))
    observables=['mother_R_ref','mother_r_os','mother_pi_t','[mother_c_i]']
    #observables=[]
    data_trunc=expand_data_with_nan(data_trunc,observables)
    
    plot_fitting_result_and_data(model,
                            fitted_parameters,
                            data_trunc,
                            parameters_to_fit,
                            subplot=True,
                            observables=observables,
                            additional_model_parameters=additional_model_parameters,
                            additional_concentrations=additional_concentrations)
  


  cellIDs=range(21) # evalutated cells
  
  if 0: #parallel fit all
      
      just_fitting_paras=True #!!!!!
      #cellIDs=[3] # testcase
      p=Pool(20)
      fitted_parameters_list=[]
      fitted_opt_res_list=[]
      f = ft.partial(fit_single_cell, plotting=False, min_max_bounds_fact=min_max_bounds_fact, just_fitting_paras=just_fitting_paras)
      
      if just_fitting_paras:
        fitted_parameters_list = p.map(f, cellIDs)
        df = pd.DataFrame(fitted_parameters_list)
        df.to_csv('fitted_parameters_parallel.csv')
      else:
        fitted_opt_res_list  = p.map(f, cellIDs)
        with open('df_res_opt_parallel.p', 'w') as f:
          pickle.dump(fitted_opt_res_list,f) 
        
      df.to_csv('fitted_parameters_parallel.csv')
  else:
      df = pd.read_csv('fitted_parameters_parallel.csv', index_col=0) 


  #check_likelihood_for_para(para_to_test, values)

  if 0: #parallel fit all test parameter values and chi^2 
      para_to_test = 'extens_factor'
      df_dict={}
      opt_results_dict={}
      values = np.linspace(10, 400, 20)

      parameters_to_fit.remove(para_to_test)
      
      for value in values:
          additional_model_parameters[para_to_test] =  value
          p=Pool(20)
          fitted_parameters_list=[]
          #cellIDs=range(1) # testcase
          f = ft.partial(fit_single_cell, plotting=False, min_max_bounds_fact=min_max_bounds_fact)
          fitted_parameters_list = p.map(f, cellIDs)
          df_test = pd.DataFrame(fitted_parameters_list)
          df_dict[value] = df_test 
          p.close() 
          plt.figure(4)
          plt.plot(value, df_test.MSD.mean(),'.r')
          
      plt.xlabel(para_to_test)
      plt.ylabel(r'$\chi^2$')

      with open('df_dict_likelihood_' + para_to_test + '.p', 'w') as f:
        pickle.dump(df_dict,f) 

            #df.to_csv('fitted_parameters_parallel_maxF_'+str(maxFact)+'.csv')
  #else:
   #    with open('df_dict_likelihood_' + para_to_test + '.p', 'r') as f:
    #    df_dict = pickle.load(f)
        #df = pd.read_csv('fitted_parameters_parallel_2_norm.csv', index_col=0) 
 # if 1: #plot liklihood
  #    plot_liklihood(df_dict, para_to_test)


  if 1: # plot all
    only_data=False
    plt.style.use('seaborn-whitegrid')
    plt.style.use('seaborn-colorblind')
    #for i,style in enumerate(plt.style.available):

    plt.figure(1,figsize=(17,14))

    observables = []
    df_p = df.drop(['MSD','CellID'], axis=1)

    for cellID in cellIDs:

      fitted_parameters = df_p.iloc[cellID]
      parameters_to_fit = df_p.columns
      
      data = {'time': np.array(time_data), 'mother_V_tot_fl': mothercells_data[cellID, :],'bud_V_tot_fl': daughtercells_data[cellID, :]}
      data_trunc = get_data_trunc(data,max_time)

      rows_and_cols = ()
      rows_and_cols = (len(cellIDs)/4 + 1, 4)

      ax = plt.subplot(rows_and_cols[0] ,rows_and_cols[1], cellID + 1)
      ax.set_title('cell ID: {0}'.format(cellID))
      try:

        plot_fitting_result_and_data(model,
                                      fitted_parameters,
                                      data_trunc,
                                      parameters_to_fit,
                                      subplot=False,
                                      observables=observables,
                                      additional_model_parameters=additional_model_parameters,
                                      additional_concentrations=additional_concentrations,
                                      legend=False,
                                      only_data=only_data)

      except:
        print('fail to plot simulation of cellID' + str(cellID))
      ax.set_ylabel('volume, fl', fontsize='large')
      ax.set_xlabel('time, min', fontsize='large')
      
    plt.legend(bbox_to_anchor=(1.5, 1),
               loc=2, borderaxespad=0.,
               labels=['bud volume (fitted)', 'bud volume ', 'mother volume (fitted)', 'mother volume ']) 
    #plt.annotate(style, (0,0))


    plt.tight_layout()
    plt.savefig('plots/all_cells_fitted_data.png', dpi=300)




    
  ''' for cellID in range(1):
                    fit_paras  = fit_single_cell(cellID)
                    print(fit_paras)
              '''
      
      #plt.savefig('plots/fits_and_add_Observables_cellID'+str(cellID)+'.png')
#fitted_parameters = fit_model_to_all_cells(model, mothercells_data, daughtercells_data, time_data, parameters_to_fit)
#plot_fitting_result_and_data(model, fitted_parameters, mothercells_data, parameters_to_fit, subplot=True)


