import fit_data
import model_data
import simulate
import matplotlib.pyplot as plt
import numpy as np
import tellurium as tet
import pandas as pd
import seaborn as sns
from multiprocessing import Pool

plt.style.use('ggplot')



def volume_to_radius(volume):
    r = (3./4/np.pi*volume) ** (1./3)
    return r

def get_initial_bud_start_guess(data):
  budding_start = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]
  return budding_start 


def get_additional_model_parameters(data):
  assert data['time'][0] == 0
  #budding_start = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]
  volume_mother_0 = data['mother_V_tot_fl'][0]
  volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
  r_tot_mother = volume_to_radius(volume_mother_0)
  r_tot_bud = volume_to_radius(volume_bud_0)
  r_os_mother = r_tot_mother - model['mother_r_b']
  r_os_bud = r_tot_bud - model['bud_r_b']   
  additional_model_parameters = { 'mother_r_os': r_os_mother}#,
                                  #'bud_r_os': r_os_bud,                                   
                                  #'budding_start': budding_start }     
  return additional_model_parameters

def get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id):
  data = {'time': np.array(time_data), 
          'mother_V_tot_fl': mothercells_data[cell_id, :],
          'bud_V_tot_fl': daughtercells_data[cell_id, :]}
  return data




def fit_single_mother_and_bud(model,
                              mothercells_data,
                              daughtercells_data,
                              time_data,
                              cell_id,
                              parameters_to_fit, 
                              additional_concentrations={},
                              max_time=400,
                              params_ini={},
                              tolerance_factor=1):
  
  bounds={}
  data = get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id)
  data_trunc = model_data.truncate_data(data)
  data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
  data_trunc = model_data.limit_time_in_data(data_trunc, max_time=max_time)
  additional_model_parameters = get_additional_model_parameters(data_trunc)
  
  bud_start_data = get_initial_bud_start_guess(data_trunc)
  start_tolerance = 100

  if params_ini=={}:
    params_ini['budding_start'] = bud_start_data
    print(params_ini)
  
  bounds={}
  bounds['budding_start'] = [bud_start_data - start_tolerance, bud_start_data + start_tolerance]
  

  fitted_parameters = fit_data.fit_model_to_data(model, 
                                                 data_trunc, 
                                                 parameters_to_fit,   
                                                 'cmaes', 
                                                 additional_model_parameters=additional_model_parameters,
                                                 additional_concentrations=additional_concentrations,
                                                 params_ini=params_ini,
                                                 tolerance_factor=tolerance_factor,
                                                 bounds=bounds)
  msd = fit_data.compute_objective_function(fitted_parameters, 
                                            model, 
                                            parameters_to_fit, 
                                            data_trunc, 
                                            additional_model_parameters, 
                                            additional_concentrations)
  return fitted_parameters, msd

def f(cell_id): # function for parallel fitting 
    print 'fitting cell no %d out of %d' %(cell_id, len(mothercells_data))

    fitted_parameters, msd = fit_single_mother_and_bud(model, 
                                                       mothercells_data, 
                                                       daughtercells_data, 
                                                       time_data, 
                                                       cell_id, 
                                                       parameters_to_fit, 
                                                       additional_concentrations,
                                                       max_time=max_time,
                                                       params_ini={},#params_ini,
                                                       tolerance_factor=1)#tolerance_factor)
    
    return fitted_parameters.tolist() + [msd]
  

def fit_all_mother_bud(model,
                       mothercells_data,
                       daughtercells_data,
                       time_data,
                       parameters_to_fit,
                       additional_concentrations={},
                       max_time=400,
                       params_ini={},
                       tolerance_factor=1):


  
  assert len(mothercells_data) == len(daughtercells_data)
  
  p = Pool(40) 
  cell_ids = range(len(mothercells_data))
  fitted_params_list = p.map(f, cell_ids)
  df_params = pd.DataFrame(fitted_params_list)
  df_params.columns = parameters_to_fit + ['MSD']

  return df_params      

def plot_fitting_for_all(model,
                         mothercells_data,
                         daughtercells_data,
                         time_data,
                         df_params,
                         additional_concentrations={},
                         max_time=400):
  assert len(df_params) == len(mothercells_data)
  df_params = df_params.copy()
  df_params = df_params.drop('MSD', axis=1)
  no_cols = 4
  no_rows = np.ceil(float(len(mothercells_data))/no_cols)
  fig = plt.figure(1, figsize=(10,12))
  ax = plt.subplot(no_rows, no_cols, 1)
  for cell_id in range(len(mothercells_data)):
    data = get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id)
    data_trunc = model_data.truncate_data(data)
    data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
    data_trunc = model_data.limit_time_in_data(data_trunc, max_time=max_time)
    additional_model_parameters = get_additional_model_parameters(data_trunc)
    print cell_id
    print additional_model_parameters
    print 
    fitted_parameters = df_params.iloc[cell_id]

    ax = plt.subplot(no_rows, no_cols, cell_id + 1, sharey=ax)
    fit_data.plot_fitting_result_and_data(model, 
                                          fitted_parameters, 
                                          data_trunc, 
                                          parameters_to_fit, 
                                          subplot=False, 
                                          additional_model_parameters=additional_model_parameters,
                                          additional_concentrations=additional_concentrations,
                                          legend=False)
    plt.ylabel('volume, fl')
    plt.xlabel('time, min')
    plt.title('Cell %s' %cell_id)
  plt.legend(bbox_to_anchor=(1.5, 1),
             loc=2, borderaxespad=0.,
             labels=['bud volume (fitted)',
             'bud volume ',
             'mother volume (fitted)',
             'mother volume '])

  plt.tight_layout()




def plot_parameter_distribution(df_params):
  #df_params = df_params.copy()
  #df_params = df_params.drop('MSD', axis=1)
  df_stacked = df_params[['k_nutrient', 'k_deg']].stack()
  df_stacked = df_stacked.reset_index()
  df_stacked.columns = ['cell_id', 'parameter', 'value']
  
  plt.figure(2)
  sns.swarmplot(x="parameter", y="value", data=df_stacked)
  plt.ylim(0, 1e-13)

  #plt.figure(3)
  ax = df_params.plot(x='mother_phi', y='bud_phi', kind='scatter', color=df_params['MSD'], colormap='winter', s=100, fontsize='x-large')
  plt.xlabel('extensility (mother),  $Pa^-1 s^-1$',fontsize='x-large')
  plt.ylabel('extensility (bud),  $Pa^-1 s^-1$',fontsize='x-large')
  #cb = plt.colorbar(ax)
  #cb.set_label('mean square displacement', rotation=270)



if __name__ == '__main__':

    #budding_start=100  
    max_time=400
    cells_to_test=4
    mothercells_data, daughtercells_data, time_data = model_data.load_data()
    
    #mothercells_data = reduce_time_in_data(mothercells_data, time_data, max_time=400)
    #daughtercells_data = reduce_time_in_data(daughtercells_data, time_data, max_time=400)
    #time_data = reduce_time_in_data(time_data, time_data, max_time=400)
    
    # 1 for test on cells_to_test
    if 0:
      mothercells_data = mothercells_data[:cells_to_test, :]
      daughtercells_data = daughtercells_data[:cells_to_test, :]

    
    time_data = time_data * 60 # convert time to seconds
    model = simulate.load_model('volume_mother_and_bud.txt')
    
    parameters_to_fit = ['budding_start','k_nutrient', 'k_deg', 'mother_phi', 'bud_phi']
    

    additional_concentrations = {'[mother_c_i]': 325,
                                 '[bud_c_i]': 325 }
    #params_ini = {'budding_start': budding_start}   # not mentioned initial parameters are taken from the model                          
    
    # 1 for fitting                             
    if 1:
      df_params = fit_all_mother_bud(model, 
                                     mothercells_data, 
                                     daughtercells_data, 
                                     time_data, 
                                     parameters_to_fit, 
                                     additional_concentrations,
                                     max_time=max_time,
                                     tolerance_factor = 1)
  
      df_params.to_csv('fitted_parameters_parallel.csv')
    else:
      df_params = pd.read_csv('fitted_parameters_parallel.csv', index_col=0)

    if 0:
      plot_fitting_for_all(model,
                           mothercells_data,
                           daughtercells_data,
                           time_data,
                           df_params,
                           additional_concentrations,
                           max_time=max_time)
      if 0:
        plt.savefig('plots/fits_parallel.png')
    else:
      df_stacked = plot_parameter_distribution(df_params)




