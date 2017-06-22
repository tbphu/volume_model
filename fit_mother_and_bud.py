import fit_data
import model_data
import simulate
import matplotlib.pyplot as plt
import numpy as np
import tellurium as tet
import pandas as pd
plt.style.use('ggplot')



def volume_to_radius(volume):
    r = (3./4/np.pi*volume) ** (1./3)
    return r

def get_additional_model_parameters(data):
  assert data['time'][0] == 0
  budding_start = data['time'][~np.isnan(data['bud_V_tot_fl'])][0]
  volume_mother_0 = data['mother_V_tot_fl'][0]
  volume_bud_0 = data['bud_V_tot_fl'][~np.isnan(data['bud_V_tot_fl'])][0]
  r_tot_mother = volume_to_radius(volume_mother_0)
  r_tot_bud = volume_to_radius(volume_bud_0)
  r_os_mother = r_tot_mother - model['mother_r_b']
  r_os_bud = r_tot_bud - model['bud_r_b']   
  additional_model_parameters = { 'mother_r_os': r_os_mother,
                                  'bud_r_os': r_os_bud,                                   
                                  'budding_start': budding_start }     
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
                              additional_concentrations={}):
  data = get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id)
  data_trunc = model_data.truncate_data(data)
  data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
  additional_model_parameters = get_additional_model_parameters(data_trunc)
  fitted_parameters = fit_data.fit_model_to_data(model, 
                                                 data_trunc, 
                                                 parameters_to_fit,   
                                                 'cmaes', 
                                                 additional_model_parameters=additional_model_parameters,
                                                 additional_concentrations=additional_concentrations)
  return fitted_parameters

def fit_all_mother_bud(model,
                       mothercells_data,
                       daughtercells_data,
                       time_data,
                       parameters_to_fit,
                       additional_concentrations={}):
  assert len(mothercells_data) == len(daughtercells_data)
  df_params = pd.DataFrame(columns=parameters_to_fit)
  for cell_id in range(len(mothercells_data)):
    print 'fitting cell no %d out of %d' %(cell_id, len(mothercells_data))
    fitted_parameters = fit_single_mother_and_bud(model, 
                                                  mothercells_data, 
                                                  daughtercells_data, 
                                                  time_data, 
                                                  cell_id, 
                                                  parameters_to_fit, 
                                                  additional_concentrations)
    df_params.loc[len(df_params)] = fitted_parameters
  return df_params

def plot_fitting_for_all(model,
                         mothercells_data,
                         daughtercells_data,
                         time_data,
                         df_params,
                         additional_concentrations={}):
  assert len(df_params) == len(mothercells_data)
  for cell_id in range(len(mothercells_data)):
    data = get_data_dict_for_cell(mothercells_data, daughtercells_data, time_data, cell_id)
    data_trunc = model_data.truncate_data(data)
    data_trunc['time'] = data_trunc['time'] + abs(min(data_trunc['time']))
    additional_model_parameters = get_additional_model_parameters(data_trunc)
    fitted_parameters = df_params.iloc[cell_id]

    fit_data.plot_fitting_result_and_data(model, 
                                          fitted_parameters, 
                                          data_trunc, 
                                          parameters_to_fit, 
                                          subplot=True, 
                                          additional_model_parameters=additional_model_parameters,
                                          additional_concentrations=additional_concentrations)





    


if __name__ == '__main__':
    mothercells_data, daughtercells_data, time_data = model_data.load_data()

    mothercells_data = mothercells_data[:3, :]
    daughtercells_data = daughtercells_data[:3, :]

    
    time_data = time_data * 60 # convert time to seconds
    model = simulate.load_model('volume_mother_and_bud.txt')
    
    parameters_to_fit = ['k_nutrient', 'k_deg','mother_phi','bud_phi']
    

    additional_concentrations = {'[mother_c_i]': 325,
                                 '[bud_c_i]': 325 }

    if 0:
      df_params = fit_all_mother_bud(model, 
                                     mothercells_data, 
                                     daughtercells_data, 
                                     time_data, 
                                     parameters_to_fit, 
                                     additional_concentrations)
  
      df_params.to_csv('fitted_parameters.csv')
    else:
      df_params = pd.read_csv('fitted_parameters.csv', index_col=0)

    plot_fitting_for_all(model,
                         mothercells_data,
                         daughtercells_data,
                         time_data,
                         df_params,
                         additional_concentrations)




df=pd.DataFrame(fit_para_list) 
df.columns=parameters_to_fit   
plt.figure(2)
df.mean().plot(kind='bar')
    #fitted_parameters = fit_model_to_all_cells(model, mothercells_data, daughtercells_data, time_data, parameters_to_fit)
    #plot_fitting_result_and_data(model, fitted_parameters, mothercells_data, parameters_to_fit, subplot=True)
