import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cma
from scipy.optimize import curve_fit


plt.style.use('ggplot')


def calculate_msd(data,func):
	
	msd = 0
	msd = np.mean(((data-func)**2)/len(data))
	
	return msd


def linear_f(x, paras):
	return x * paras[0] + paras[1]

def linF(x,m):
	return x*m	


def objective_function(paras, data_x, data_y):
	
	y = linear_f(data_x, paras)
	msd = calculate_msd(data_y,y)
	return msd

def fit_cmaes(paras_0, data, sigma0=4e-16, tolx=1e-20): #bounds_min, bounds_max,
    
    options = cma.CMAOptions()
    #options.set('tolfun', 1e-14)
    #options['tolx'] = tolx
    #options['bounds'] = (bounds_min, bounds_max)
    #param_vec = np.array(initial_params)
    #p_min = max(param_vec.min(), 1e-20)
    #options['scaling_of_variables'] = param_vec / p_min
    result = cma.fmin(objective_function, 
                      x0=paras_0,
                      sigma0=sigma0,
                      options=options,
                      args=data)
    return result



def plot(x,y,f_y):


	plt.figure(1, figsize=(5,5))
	plt.scatter(x,y)#, color='b')
	plt.plot(x,f_y, ls='--', dashes=(3,8))
	plt.xlabel('$k_{nutient}$',fontsize='large')
	plt.ylabel('$k_{deg}$',fontsize='large')
	plt.xlim(0,3.e-14)
	plt.ylim(0,3.e-14)
	plt.show()



def fit_k(fitted_parameters, paras_0=[1.,0.],sigma0=3):
	
	fitted_parameters['k_nutrient'] = fitted_parameters['k_nutrient']*1.e+14
	fitted_parameters['k_deg'] = fitted_parameters['k_deg']*1.e+14
	

	options = cma.CMAOptions()
	options['tolx'] = 1e-20

	f_paras =  cma.fmin(objective_function,
						 paras_0,
						 sigma0,
						 args=(fitted_parameters['k_nutrient'],fitted_parameters['k_deg']))#, options=options)
    #fit_cmaes(paras_0, data)#, [0.,-1],[5,1])
	print('m={0}, n={1}'.format(f_paras[0][0],f_paras[0][1]))
	msd = calculate_msd(fitted_parameters['k_deg'],linear_f(fitted_parameters['k_nutrient'],f_paras[0]))
	print('MSD:{0}'.format(msd)	)
	return f_paras




if __name__ == '__main__':


	fitted_parameters = pd.read_csv('fitted_parameters_parallel.csv', index_col=0)
	fitted_parameters.sort(['k_nutrient'],inplace=True)
	#fitted_parameters['k_nutrient'] = fitted_parameters['k_nutrient']*1.e+14
	#fitted_parameters['k_deg'] = fitted_parameters['k_deg']*1.e+14
	#fitted_parameters.plot(['k_nutrient'],['k_deg'], kind='scatter')


	
	#data=np.array([fitted_parameters['k_nutrient'],fitted_parameters['k_deg']])
	#x=data[0]
	#x=fitted_parameters['k_nutrient']#*1.e+14
	#y=fitted_parameters['k_deg']#*1.e+14
	#f_paras = fit_k(fitted_parameters)

	
 	#calculated_k_nut = linear_f(fitted_parameters['k_nutrient'],f_paras[0])
 	#plot(fitted_parameters['k_nutrient'],fitted_parameters['k_deg'],calculated_k_nut)

 	bounds=([-np.inf,-1e-14],[np.inf,1e-14])
 	pop, pcov = curve_fit(linF,fitted_parameters['k_nutrient'],fitted_parameters['k_deg'])#,bounds=bounds)
 	
 	plot (fitted_parameters['k_nutrient'],
 			fitted_parameters['k_deg'],
 			linF(fitted_parameters['k_nutrient'], pop[0]))
 	plt.savefig('plots/k_nutr_VS_k_deg.png', dpi=300)
 	print('f(x) = m*x  \nm:{0}'.format(pop[0]))
 	print('pcov: {0}'.format(pcov[0][0]))
 	print('standard deviation: {0}'.format(np.sqrt(np.diag(pcov))))

	
		


	

