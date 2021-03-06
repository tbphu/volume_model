import fit_data
import model_data
import simulate
import matplotlib.pyplot as plt
import numpy as np
import tellurium as tet
import pandas as pd
import seaborn as sns
from multiprocessing import Pool
#from functools import partial
import fit_kn_vs_kdeg

plt.style.use('ggplot')



def volume_to_radius(volume):
    r = (3./4/np.pi*volume) ** (1./3)
    return r

def get_initial_bud_start_guess(data):
  budding_start = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]
  return budding_start 


'''def get_initial_parameter_from_data(data,params):
  para_ini = {}
  assert data['time'][0] == 0, 'data does not start at time 0'
  
  if 'budding_start' in params:
    para_ini['budding_start'] = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]

  if 'mother_r_os_0' in paras:
    volume_mother_0 = data['mother_V_tot_fl'][0]
    r_tot_mother_0 = volume_to_radius(volume_mother_0)
    para_ini['mother_r_os_0'] = r_tot_mother_0 - model['mother_r_b']

  if 'bud_r_os_0' in paras:  
    volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
    r_tot_bud_0 = volume_to_radius(volume_bud_0)
    para_ini['bud_r_os_0'] = r_tot_bud_0 - model['bud_r_b']

  return para_ini'''


'''def get_additional_model_parameters(data):
  assert data['time'][0] == 0
  #budding_start = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]
  volume_mother_0 = data['mother_V_tot_fl'][0]
  volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
  r_tot_mother_0 = volume_to_radius(volume_mother_0)
  r_tot_bud_0 = volume_to_radius(volume_bud_0)
  r_os_mother_0 = r_tot_mother_0 - model['mother_r_b']
  r_os_bud_0 = r_tot_bud_0 - model['bud_r_b']   
  additional_model_parameters = { 'mother_r_os': r_os_mother_0}#,
                                  #'bud_r_os': r_os_bud,                                   
                                  #'budding_start': budding_start } 
  #print(r_os_mother_0)                                    
  return additional_model_parameters'''

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
                              additional_model_parameters={},
                              max_time=400,
                              params_ini={},
                              log=True,
                              bounds={},
                              tolfun=0.1,
                              sigma0_log_notlog = (1e-1,4e-17)):
  additional_model_parameters={}
  data = get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id)
  data_trunc = model_data.truncate_data(data)
  data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
  data_trunc = model_data.limit_time_in_data(data_trunc, max_time=max_time)
  #additional_model_parameters =  get_additional_model_parameters(data_trunc)
  para_ini_ids = ['budding_start', 'mother_r_os_0']
  params_ini=fit_data.get_initial_parameter_from_data(model,data_trunc,para_ini_ids)
  #bud_start_data = get_initial_bud_start_guess(data_trunc)
    
  bud_start_data = params_ini['budding_start']
  start_tolerance = 120*60
  
  bounds['budding_start'] = [bud_start_data - start_tolerance, bud_start_data + start_tolerance]
  bounds['mother_r_os_0'] = [params_ini['mother_r_os_0']-0.3,params_ini['mother_r_os_0']+0.3] 
  #bounds['mother_phi']= [5.e-6, 5.e-1]
  #bounds['bud_phi']= [1.e-5, 5.e-1 ]
  #bounds['k_nutrient_']= [9.e-16, 9.e-17]
  #bounds['k_deg_0']= [1.e-15, 1.e-17]
  #bounds['mother_r_os'] = [0.1,2]
  
  fitted_parameters = fit_data.fit_model_to_data(model, 
                                                 data_trunc, 
                                                 parameters_to_fit,   
                                                 'cmaes', 
                                                 additional_model_parameters=additional_model_parameters,
                                                 additional_concentrations=additional_concentrations,
                                                 params_ini=params_ini,
                                                 bounds=bounds,
                                                 log=log,
                                                 tolfun=tolfun,
                                                 sigma0_log_notlog = sigma0_log_notlog)
  print(fitted_parameters)
  log=False
  msd = fit_data.compute_objective_function(fitted_parameters, 
                                            model, 
                                            parameters_to_fit, 
                                            data_trunc, 
                                            additional_model_parameters, 
                                            additional_concentrations,log)
  return fitted_parameters, msd



def f(cell_id): ##importent): 

      # function for parallel fitting 
    print 'fitting cell no %d out of %d' %(cell_id, len(mothercells_data))
    print(sigma0_log_notlog)

    fitted_parameters, msd = fit_single_mother_and_bud(model, 
                                                       mothercells_data, 
                                                       daughtercells_data, 
                                                       time_data, 
                                                       cell_id, 
                                                       parameters_to_fit, 
                                                       additional_concentrations=additional_concentrations,
                                                       additional_model_parameters=additional_model_parameters,
                                                       max_time=max_time,
                                                       params_ini=params_ini,
                                                       log=log,
                                                       bounds=bounds,
                                                       tolfun=tolfun,
                                                       sigma0_log_notlog=sigma0_log_notlog)
    
    return fitted_parameters.tolist() + [msd]
  

def fit_all_mother_bud(mothercells_data,
                       daughtercells_data,
                       parameters_to_fit):


  
  assert len(mothercells_data) == len(daughtercells_data)
  p = Pool(30) 
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
                         additional_model_parameters={},
                         max_time=400,
                         cols=4,
                         max_cell_ID=0):
  
  if max_cell_ID == 0:
    Cell_IDs = range(len(mothercells_data))
  else:
    Cell_IDs = range(max_cell_ID) 




  assert len(df_params) == len(mothercells_data)
  df_params = df_params.copy()
  df_params.drop('MSD', inplace=True, axis=1)
  no_cols = cols
  no_rows = np.ceil(float(len(Cell_IDs))/no_cols)
  fig = plt.figure(1, figsize=(8,8))
  ax = plt.subplot(no_rows, no_cols, 1)
    

  for cell_id in Cell_IDs:
    data = get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id)
    data_trunc = model_data.truncate_data(data)
    data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
    data_trunc = model_data.limit_time_in_data(data_trunc, max_time=max_time)
    #additional_model_parameters = get_additional_model_parameters(data_trunc)
    print cell_id
    #print additional_model_parameters
    fitted_parameters = df_params.iloc[cell_id]

    ax = plt.subplot(no_rows, no_cols, cell_id + 1, sharey=ax)

    try:
      fit_data.plot_fitting_result_and_data(model, 
                                            fitted_parameters, 
                                            data_trunc, 
                                            parameters_to_fit, 
                                            subplot=False, 
                                            additional_model_parameters=additional_model_parameters,
                                            additional_concentrations=additional_concentrations,
                                            legend=False,
                                            observables=[])
    except:
      print('could not simulate for cell Nr. {0}'.format(cell_id))

    budstart = fitted_parameters['budding_start']/60
    nr_of_points = 10
    y = np.arange(0,100,nr_of_points)
    x = [budstart]*nr_of_points
    plt.plot(x,y)       



    plt.ylabel('volume, fl', fontsize='x-large')
    plt.xlabel('time, min', fontsize='x-large')
    plt.title('Cell %s' %cell_id)
  plt.legend(bbox_to_anchor=(1.5, 1),
             loc=2, borderaxespad=0.,
             labels=['bud volume (fitted)',
             'bud volume ',
             'mother volume (fitted)',
             'mother volume '])

  plt.tight_layout()



def set_colorbar_title(text,cbID=2):
    f1 = plt.gcf()
    cax=f1.get_axes()[cbID]
    cax.set_ylabel(text, fontsize='x-large')


def plot_parameter_distribution(df_params):
  #df_params = df_params.copy()
  #df_params = df_params.drop('MSD', axis=1)
  
  # linear fit through 0,0

  max_k = 3.e-14
  k_nut = np.arange(0,max_k + max_k / 10,max_k / 10)
  
  pop,pcov = fit_kn_vs_kdeg.fit_linear_k(df_params)
  fitted_kdeg = fit_kn_vs_kdeg.linF(k_nut, pop[0]) 
  
  print('f(x) = m*x  \nm:{0}'.format(pop[0]))
  print('pcov: {0}'.format(pcov[0][0]))
  print('standard deviation: {0}'.format(np.sqrt(np.diag(pcov))))



  df_stacked = df_params[['k_nutrient', 'k_deg']].stack()
  df_stacked = df_stacked.reset_index()
  df_stacked.columns = ['cell_id', 'parameter', 'value']
  df_stacked['value']= df_stacked['value']*1e+15
  
  
  plt.figure(2,figsize=(4,4))
  sns.swarmplot(x="parameter", y="value", data=df_stacked, size=6)
  plt.savefig('plots/fitted_ks_parallel.png')

  #plt.figure(3)
  fig, ax = plt.subplots(2,1)
  fig.set_size_inches(7, 12)

  lim0 = (1.e-7,5.e-1)
  
  
  df_params.plot(x='mother_phi',
                 y='bud_phi',
                 kind='scatter', 
                 color=df_params['MSD'],
                 colormap='winter',
                 s=75,
                 fontsize='xx-large',
                 ax=ax[0])
  
  set_colorbar_title('MSD',cbID=2)



  ax[0].plot(np.arange(1.e-7, 5.e-1, 5e-6),
            np.arange(1.e-7, 5.e-1, 5e-6),
            'r--',
            alpha = 0.2,
            linewidth=2,
            color = 'gray')

  ax[0].set_xscale("log")#, nonposx='clip')
  ax[0].set_xlim(lim0)
  ax[0].set_yscale("log")#, nonposy='clip')
  ax[0].set_ylim(lim0)
  ax[0].set_xlabel('extensibility $\phi$ (mother),  $Pa^{-1} s^{-1}$',fontsize='xx-large')
  ax[0].set_ylabel('extensibility $\phi$ (bud),  $Pa^{-1} s^{-1}$',fontsize='xx-large')
  
  df_params.plot(x='k_nutrient',
                 y='k_deg',
                 kind='scatter',
                 color=df_params['MSD'],
                 colormap='winter',
                 s=75,
                 fontsize='xx-large',
                 ax=ax[1])
  
  set_colorbar_title('MSD', cbID=3)



  ax[1].plot(k_nut, fitted_kdeg,'g--')
  ax[1].plot(np.arange(0, max_k + max_k / 10, max_k / 10),
            np.arange(0, max_k+ max_k / 10, max_k / 10),
            '--',
            alpha = 0.2,
            linewidth=2,
            color='gray')

  ax[1].set_xlabel('$k_{nutrient}$, mM $s^{-1}$ $um^{-2}$', fontsize='xx-large')
  ax[1].set_ylabel('$k_{deg}$, mM $s^{-1}$ $um^{-3}$', fontsize='xx-large')
  lim1 = (1.e-15, max_k)
  ax[1].set_xlim(lim1)
  ax[1].set_ylim(lim1)
  





  plt.tight_layout()
  plt.savefig('plots/fitted__MB_phi__k_parallel.png')

  

  plt.figure(5)
  df_stacked_budstart = df_params[['budding_start']].stack()
  df_stacked_budstart = df_stacked_budstart.reset_index()
  df_stacked_budstart.columns = ['cell_id', 'parameter', 'value']
  #df_stacked['value']= df_stacked['value']*1e+15

  sns.swarmplot(x="parameter", y="value", data=df_stacked_budstart, size=6)
  plt.savefig('plots/fitted_budding_start.png')




  #cb = plt.colorbar(ax)
  #cb.set_label('mean square displacement', rotation=270)
  return df_stacked  




if __name__ == '__main__':

    max_time=500*60
    cells_to_test=4
    mothercells_data, daughtercells_data, time_data_min = model_data.load_data()

    # convert time to seconds
    time_data = [x*60 for x in time_data_min]

    # 1 parameter in log scale
    if 1:
      log=True
    else:
      log=False

    sigma0_log_notlog = (1e-1,1e-17)  

    # 1 tolerance of the optimizer at the minimum
    tolfun=1e-6  
    
    # 1 for test on cells_to_test
    if 1:
      mothercells_data = mothercells_data[:cells_to_test, :]
      daughtercells_data = daughtercells_data[:cells_to_test, :]
    
    model = simulate.load_model('volume_mother_and_bud.txt')
    
    parameters_to_fit = ['budding_start','k_nutrient', 'k_scaling_factor', 'mother_phi', 'bud_phi','mother_r_os_0']
#    parameters_to_fit = ['budding_start','k_nutrient', 'k_deg_0', 'mother_phi', 'bud_phi','mother_r_os_0']
     

    additional_concentrations = {'init([mother_c_i])': 319.4,
                                 'init([bud_c_i])': 319.4 }
    additional_model_parameters = {'withSF': 1,
                                   'mother_bud_water_perm':1,}                             
    #params_ini = {'budding_start': budding_start}   # not mentioned initial parameters are taken from the model                          
    params_ini = {'k_scaling_factor': 1.2}  
    
    bounds={}
    #bounds['mother_phi']= [5.e-6, 5.e-1]
    #bounds['bud_phi']= [1.e-5, 5.e+2 ]
    #bounds['k_nutrient_']= [9.e-16, 9.e-17]
    #bounds['k_deg_0']= [1.e-15, 1.e-17]
    #bounds['mother_r_os'] = [0.1,2]



    # 1 for fitting                             
    if 1:
      df_params = fit_all_mother_bud(mothercells_data, 
                                     daughtercells_data, 
                                     parameters_to_fit)
     
      df_params.to_csv('fitted_parameters_parallel_norm.csv')
    else:
      df_params = pd.read_csv('fitted_parameters_parallel_norm.csv', index_col=0)

    if 1:
      
      plot_fitting_for_all(model,
                           mothercells_data,
                           daughtercells_data,
                           time_data,
                           df_params,
                           additional_concentrations=additional_concentrations,
                           additional_model_parameters=additional_model_parameters,
                           max_time=max_time,
                           cols=4)#,
                           #max_cell_ID=3)

      if 1:
        plt.savefig('plots/fits_parallel.png')
       
      #df_stacked = plot_parameter_distribution(df_params)  
    #else:
     # df_stacked = plot_parameter_distribution(df_params)
 





