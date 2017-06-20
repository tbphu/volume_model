import fit_data
import model_data
import simulate
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
plt.style.use('ggplot')


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




    
    
    parameters_to_fit = ['k_nutrient', 'k_deg', 'mother_r_os', 'bud_r_os']
    
    data = {'time': np.array(time_data), 
            'mother_V_tot_fl': mothercells_data[1, :], 
            'bud_V_tot_fl': daughtercells_data[1, :]}

    data = model_data.truncate_data(data)


    budding_start = abs(min(data['time']))
    data['time'] = data['time'] + budding_start
    additional_model_parameters = {'[mother_c_i]': 319,
                                   '[bud_c_i]': 319,
                                   'budding_start': budding_start }  
    #print additional_model_parameters
    fitted_parameters = fit_data.fit_model_to_data(model, 
                                          data, 
                                          parameters_to_fit,   
                                          'cmaes', additional_model_parameters=additional_model_parameters) 
                                          #additional_model_parameters=additional_model_parameters)
    print fitted_parameters
    fit_data.plot_fitting_result_and_data(model, fitted_parameters, data, parameters_to_fit, subplot=True, additional_model_parameters=additional_model_parameters)

    #fitted_parameters = fit_model_to_all_cells(model, mothercells_data, daughtercells_data, time_data, parameters_to_fit)
    #plot_fitting_result_and_data(model, fitted_parameters, mothercells_data, parameters_to_fit, subplot=True)
