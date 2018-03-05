import tellurium as te
te.setDefaultPlottingEngine('matplotlib')
#%matplotlib inline 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import simulate
import numpy as np
import seaborn as sns
sns.set_style("ticks")
import matplotlib.gridspec as gridspec
import FigureTools as FT




def apply_param_dict_to_model(model,param_dict):
    for param in param_dict.keys():
        if param_dict[param] is not None: # do not change parameters marked with value==None.
            #print param_dict[param]
            model[param] = param_dict[param]
    return model

volume_and_hog_model = te.loadAntimonyModel("volume_and_hog.txt")
volume_model = te.loadAntimonyModel("volume_reference_radius.txt")

volume_and_hog_model.integrator.relative_tolerance = 1e-10
volume_model.integrator.relative_tolerance = 1e-10

volume_and_hog_model_species=['time','vol_V_tot_fl', '[hog_Glyc_in]','[vol_c_i]', 'hog_Glyc_in','vol_c_i', 'hog_totalHog1PP','hog_n_totalHog1', 'hog_Hog1PPn', 'hog_Hog1PPc','hog_Hog1n', 'hog_Hog1c', '[hog_Hog1PPn]', '[hog_Hog1PPc]','[hog_Hog1n]', '[hog_Hog1c]', 'vol_pi_t', 'vol_pi_i']
volume_model_species=['time','V_tot_fl','V_ref','r_os','r_b', '[c_e]','[c_i]','c_i', 'pi_t', 'pi_i']

end_time = 25000

n = 10

def param_plot(i,what):
	global volume_model
	
	volume_model.resetAll()
	# Initiate with a small ~5-10 fL, for HOG model - Discussion! :
	#volume_and_hog_model.resetAll()
	no_shock = volume_model['shock_set_c_e']
	no_hyposhock = volume_model['hyposhock_set_c_e']
	if what == 'event_time':
		shock_event_time = (i*30 + 40)*60 #volume_model['shock_event_time'] 
		hyposhock_event_time = (i*30 + 40+30)*60 #volume_model['hyposhock_event_time']
	else:
		shock_event_time = (1*30 + 40)*60 #volume_model['shock_event_time'] 
		hyposhock_event_time = (1*30 + 40+30)*60 #volume_model['hyposhock_event_time']

	if what == 'shock_strength':
		shock_add = (i*30)+110 # raise external osmolarity by +260mM to 500mM overall.
	else:
		shock_add = 260 # raise external osmolarity by +260mM to 500mM overall.
		#shock_add = 500
	shock = no_shock + shock_add
	hyposhock = no_hyposhock + shock_add 
	# volume_hog_param_dict = {'vol_k_nutrient': None,
	#                          'vol_k_deg': None,
	#                          'vol_shock_event_time': shock_event_time,
	#                          'vol_hyposhock_event_time': hyposhock_event_time,
	#                          'vol_shock_set_c_e': shock,
	#                          'vol_hyposhock_set_c_e': hyposhock
	#                     }

	volume_param_dict = {'k_nutrient': None,
	                     'k_deg': None,
	                     'shock_event_time': shock_event_time,
	                     'hyposhock_event_time': hyposhock_event_time,
	                     'shock_set_c_e': shock,
	                     'hyposhock_set_c_e': hyposhock
	                    }


	#volume_and_hog_model = apply_param_dict_to_model(volume_and_hog_model,volume_hog_param_dict)
	volume_model = apply_param_dict_to_model(volume_model,volume_param_dict)

	#volume_and_hog_sim=volume_and_hog_model.simulate(1, end_time, end_time, selections = volume_and_hog_model_species)
	volume_sim=volume_model.simulate(1, end_time, end_time, selections = volume_model_species)

	#sns.set_palette(sns.cubehelix_palette(n,start=3,reverse=False))
	sns.set_palette(sns.cubehelix_palette(n,start=2,rot=0,reverse=False))
		
	if what == 'event_time':
		if i==2 or i==5:
			plt.subplot(gs[0,1])
			plt.gca().add_artist(FT.figure_numbering('B'))

			plt.plot(volume_sim['time']/60,volume_sim['V_tot_fl'],color = sns.color_palette()[i],label=str(shock_event_time/60))
			plt.annotate('Total volume at event', xy=(193, 51), xycoords="data",xytext=(2,62),arrowprops=dict(arrowstyle="->"))
			#plt.annotate(r'Shrinkage: $\Delta V$', xy=(200, 41), xycoords="data",xytext=(370,40),rotation=90,arrowprops=dict(arrowstyle="-[,widthB=2.3",connectionstyle='angle,angleA=90,angleB=0,rad=0.0'))
			plt.annotate('Adapted total volume', xy=(189, 30), xycoords="data",xytext=(200,15),arrowprops=dict(arrowstyle="->"))
			plt.ylim(-2,70)
			plt.xlabel('Time [minutes]')
			plt.ylabel('Total volume [fL]')

			plt.subplot(gs[1,1])
			plt.gca().add_artist(FT.figure_numbering('C'))

			plt.plot(volume_sim['time']/60,volume_sim['pi_t']/1e6,color = sns.color_palette()[i],label=str(shock_event_time/60))
			plt.xlabel('Time [minutes]')
			plt.ylabel('Turgor pressure [MPa]')

		#diff = volume_sim['V_tot_fl'] - volume_sim['V_ref']
		#y  = volume_sim['V_tot_fl'][shock_event_time-1]-np.min(volume_sim['V_tot_fl'][shock_event_time-1:])
		y  = np.min(volume_sim['V_tot_fl'][shock_event_time-1:])
		x = volume_sim['V_tot_fl'][shock_event_time-1]

		plt.subplot(gs[:,0])
		plt.gca().add_artist(FT.figure_numbering('A'))
		
		if i==2 or i == 5:
			plt.plot(x,y,'x', color = sns.color_palette()[i],label=str(shock_event_time/60) + 'min.',marker = 'o')
		else:
			plt.plot(x,y,'x', color = sns.color_palette()[i],label=str(shock_event_time/60) + 'min.')
		


		#plt.ylabel(r'Shrinkage: $\Delta V$ [fL]')
		plt.ylabel('Adapted total volume [fL]')
		plt.xlabel('Total volume at event [fL]')
		plt.xlim(-10,60)
		plt.legend(ncol=1)

	elif what == 'shock_strength':
		if i==5:
			plt.subplot(gs[0,1])
			plt.gca().add_artist(FT.figure_numbering('B'))
			plt.plot(volume_sim['time']/60,volume_sim['V_tot_fl'], color = 'orange',label = 'total volume')#, color = sns.color_palette()[i] )
			plt.plot(volume_sim['time']/60,volume_sim['V_ref'], color = 'black',label = 'reference volume')
			plt.axhline(11.5,linestyle=':',color='grey')
			plt.axhline(13,linestyle=':',color='grey')
			plt.annotate(r'$\Delta V$', xy=(300,14), xycoords="data",xytext=(300,14))
			#plt.annotate('At event: Total volume', xy=(310, 57), xycoords="data",xytext=(2,62),arrowprops=dict(arrowstyle="->"))
			#plt.annotate(r'Shrinkage: $\Delta V$', xy=(330, 45), xycoords="data",xytext=(370,40),rotation=90,arrowprops=dict(arrowstyle="-[,widthB=1.8",connectionstyle='angle,angleA=90,angleB=0,rad=0.0'))
			plt.ylim(-2,70)
			plt.xlabel('Time [minutes]')
			plt.ylabel('Volume [fL]')
			plt.legend()

			
			plt.subplot(gs[1,1])
			plt.gca().add_artist(FT.figure_numbering('C'))
			plt.plot(volume_sim['time']/60,volume_sim['pi_t']/1e6,color='black')
			plt.annotate('Adapted turgor', xy=(70-2, 0.058), xycoords="data",xytext=(190,0.1),arrowprops=dict(arrowstyle="->"))
			plt.xlabel('Time [minutes]')
			plt.ylabel('Turgor pressure [MPa]')



		diff = volume_sim['V_tot_fl'] - volume_sim['V_ref']
		y  = np.min(diff[shock_event_time:])
		x = np.min(volume_sim['pi_t'])/1e6 # 1e6 Pa -> MPa

		plt.subplot(gs[:,0])
		plt.gca().add_artist(FT.figure_numbering('A'))
		if i==5:
			plt.plot(x,y,'x', color = 'orange',marker = 'o',label=str(int(shock)) + 'mM')
		else:
			plt.plot(x,y,'x', color = sns.color_palette()[i],label=str(int(shock)) + 'mM')
		plt.ylabel(r'Adapted: $\Delta V = V-V_{ref}$ [fL]')
		plt.xlabel('Adapted tugor pressure [MPa]')
		plt.legend(ncol=1)

plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(2, 2)

for i in range(n):
	param_plot(i , 'shock_strength')
	#param_plot(i , 'event_time')

plt.tight_layout()
plt.savefig('figures/check_model_shock.png',dpi=300)
plt.savefig('figures/check_model_shock.eps')

plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(2, 2)

for i in range(n):
	#param_plot(i , 'shock_strength')
	param_plot(i , 'event_time')

plt.tight_layout()
plt.savefig('figures/check_model_event_time.png',dpi=300)
plt.savefig('figures/check_model_event_time.eps')


plt.show()


