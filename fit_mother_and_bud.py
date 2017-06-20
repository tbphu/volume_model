import fit_data
import model_data
import simulate
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
plt.style.use('ggplot')



def volume_to_radius(volume):
    r = (3./4/np.pi*volume) ** (1./3)
    return r


if __name__ == '__main__':
    mothercells_data, daughtercells_data, time_data = model_data.load_data()
    
    time_data = time_data * 60
    
    #plot_single_cell(time_data, mothercells_data[0, :], daughtercells_data[0, :])
    #plot_data(mothercells_data, daughtercells_data, time_data)

    model = simulate.load_model('volume_mother_and_bud.txt')
    #model = model_data.set_model_parameters(model, {'k_nutrient': 0})
    #model = model_data.select_model_timecourses(['V_tot_fl'])
    #simulation_result = simulate.simulate_model(model, end_time=7200)
    #simulate.plot((simulation_result,), legend=True)

    
    parameters_to_fit = ['k_nutrient', 'k_deg']
    

    cell_id = 5
    data = {'time': np.array(time_data), 
            'mother_V_tot_fl': mothercells_data[cell_id, :],
            'bud_V_tot_fl': daughtercells_data[cell_id, :]}

    data = model_data.truncate_data(data)


    #budding_start = abs(min(data['time']))
    data['time'] = data['time'] + abs(min(data['time']))
    budding_start = data['time'][ ~np.isnan(data['bud_V_tot_fl']) ][0]

    #data['time'] = data['time'] + budding_start
    

    volume_mother_0 = data['mother_V_tot_fl'][0]
    volume_bud_0 = data['bud_V_tot_fl'][ ~np.isnan(data['bud_V_tot_fl']) ][0]

    r_tot_mother = volume_to_radius(volume_mother_0)
    r_tot_bud = volume_to_radius(volume_bud_0)
    r_os_mother = r_tot_mother - model['mother_r_b']
    r_os_bud = r_tot_bud - model['bud_r_b']



    additional_model_parameters = { 'mother_r_os': r_os_mother,
                                   'bud_r_os': r_os_bud,                                   
                                   'budding_start': budding_start }  

    additional_concentrations = {'[mother_c_i]': 325,
                                   '[bud_c_i]': 325 }

    #print additional_model_parameters
    fitted_parameters = fit_data.fit_model_to_data(model, 
                                          data, 
                                          parameters_to_fit,   
                                          'cmaes', 
                                          additional_model_parameters=additional_model_parameters,
                                          additional_concentrations=additional_concentrations) 

    print 'fit params%s'  %fitted_parameters
    fit_data.plot_fitting_result_and_data(model, 
                                          fitted_parameters, 
                                          data, 
                                          parameters_to_fit, 
                                          subplot=True, 
                                          additional_model_parameters=additional_model_parameters,
                                          additional_concentrations=additional_concentrations)

    #fitted_parameters = fit_model_to_all_cells(model, mothercells_data, daughtercells_data, time_data, parameters_to_fit)
    #plot_fitting_result_and_data(model, fitted_parameters, mothercells_data, parameters_to_fit, subplot=True)
