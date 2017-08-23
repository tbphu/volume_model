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


	plt.figure(1, figsize=(8,8))
	plt.scatter(x,y, color='gray', s=14)#, color='b')
	plt.plot(x,f_y, ls='--', dashes=(8,8), color='red')
	plt.xlabel('$k_{nutrient}$, mM $s^{-1}$ $um^{-2}$',fontsize='large')
	plt.ylabel('$k_{deg}$, mM $s^{-1}$ $um^{-3}$',fontsize='large')
	plt.xlim(0,5.e-16)
	plt.ylim(0,5.e-16)
	plt.show()
	



def fit_k(fitted_parameters, paras_0=[1.,0.],sigma0=3):
	
	fitted_parameters['k_nutrient'] = fitted_parameters['k_nutrient']*1.e+14
	fitted_parameters['k_deg_0'] = fitted_parameters['k_deg_0']*1.e+14
	

	options = cma.CMAOptions()
	options['tolx'] = 1e-20

	f_paras =  cma.fmin(objective_function,
						 paras_0,
						 sigma0,
						 args=(fitted_parameters['k_nutrient'],fitted_parameters['k_deg_0']))#, options=options)
    #fit_cmaes(paras_0, data)#, [0.,-1],[5,1])
	print('m={0}, n={1}'.format(f_paras[0][0],f_paras[0][1]))
	msd = calculate_msd(fitted_parameters['k_deg_0'],linear_f(fitted_parameters['k_nutrient'],f_paras[0]))
	print('MSD:{0}'.format(msd)	)
	return f_paras

def fit_linear_k(fitted_parameters):
	fitted_parameters.sort_values('k_nutrient',inplace=True)
	pop, pcov = curve_fit(linF,fitted_parameters['k_nutrient'],fitted_parameters['k_deg_0'])
	return pop, pcov

if __name__ == '__main__':


	fitted_parameters = pd.read_csv('fitted_parameters_parallel.csv', index_col=0)
	fitted_parameters.sort_values(['k_nutrient'],inplace=True)
	#fitted_parameters['k_nutrient'] = fitted_parameters['k_nutrient']*1.e+14
	#fitted_parameters['k_deg_0'] = fitted_parameters['k_deg_0']*1.e+14
	#fitted_parameters.plot(['k_nutrient'],['k_deg_0'], kind='scatter')


	
	#data=np.array([fitted_parameters['k_nutrient'],fitted_parameters['k_deg_0']])
	#x=data[0]
	#x=fitted_parameters['k_nutrient']#*1.e+14run 
	#y=fitted_parameters['k_deg_0']#*1.e+14
	#f_paras = fit_k(fitted_parameters)

	
 	#calculated_k_nut = linear_f(fitted_parameters['k_nutrient'],f_paras[0])
 	#plot(fitted_parameters['k_nutrient'],fitted_parameters['k_deg_0'],calculated_k_nut)

 	max_deg_0=1e-14
 	x=fitted_parameters[fitted_parameters['k_deg_0']<max_deg_0]['k_nutrient']
 	y=fitted_parameters[fitted_parameters['k_deg_0']<max_deg_0]['k_deg_0']
 	#pop, pcov  = curve_fit(linF,fitted_parameters['k_nutrient'],fitted_parameters['k_deg_0'])#,bounds=bounds)
 	

 	pop, pcov  = curve_fit(linF,x,y)#,bounds=bounds)
 	
 	plot (fitted_parameters['k_nutrient'],
 			fitted_parameters['k_deg_0'],
 			linF(fitted_parameters['k_nutrient'], pop[0]))

 	plt.savefig('plots/k_nutr_VS_k_deg_0.png', dpi=300)
 	print('f(x) = m*x  \nm:{0}'.format(pop[0]))
 	print('pcov: {0}'.format(pcov[0][0]))
 	print('standard deviation: {0}'.format(np.sqrt(np.diag(pcov))))



 	plt.figure(2)
 	plt.scatter(fitted_parameters['k_nutrient'],
 			fitted_parameters['k_deg_0'],
 			c=fitted_parameters['MSD'])
 	plt.colorbar()
 	plt.xlim(0,5.e-16)
	plt.ylim(0,5.e-16)
	plt.show()

	plt.figure(3)
 	plt.scatter(fitted_parameters['mother_phi'],
 			fitted_parameters['bud_phi'],
 			c=fitted_parameters['MSD'])
 	plt.colorbar()
 	plt.xlim(-0.1e-5,1.e-5)
	plt.ylim(-0.1,0.6)
	plt.show()
 	
	
		


	

