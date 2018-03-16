from scipy.special import lambertw
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import FigureTools as FT
sns.set_style("ticks")

import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(3,5)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# constants
R = 8.314; # J/mol/K
T = 303; # K 30C

# nutrient uptake and consumptions
k_nutrient = 2.e-1 ## nmol/s/m^2 

k_scaling_factor = 1.27 # 1/um 
k_deg = k_scaling_factor*k_nutrient

c_e = 240.
pi_t = 200000.
pi_e = c_e * R * T

# water flow over membrane
Lp = 1.19e-6; # um/s/Pa #Klipp 2005

a = - ((Lp * (pi_t + pi_e))/2)

def set_b(Lp,pi_t,pi_e,k_nutrient,R,T):
    return(Lp**2 * (pi_t + pi_e)**2)/ 4 + (k_nutrient * Lp * R  * T )

def set_c(k_deg,Lp,R ,T):
    return 1./3. * (k_deg * Lp * R  * T )

def set_F0(a,b,c,r0):
    return (2. * (a * np.log(np.sqrt(b - c * r0) + a) - np.sqrt(b - c * r0) ))/c

r0 = 0.1
t0 = 0.


b = set_b(Lp,pi_t,pi_e,k_nutrient,R,T)
c = set_c(k_deg,Lp,R ,T)
F0 = set_F0(a,b,c,r0)


print np.sqrt(b - c * r0) + a


def inner_exp(t,F0,a,b,c,t0):
    return -(np.exp((-1 + (c * (F0 + t - t0))/(2 * a)))/a)


#def radius(t,F0,a,b,c,t0):
def radius(t,a,b,c,t0,F0):
    return (-a**2 + b - 2 * a**2 * lambertw(inner_exp(t,F0,a,b,c,t0)) - a**2 * lambertw(inner_exp(t,F0,a,b,c,t0))**2)/c

def radius_taylor_series(t,a,b,c,t0,F0):
    return (-a**2 + b - 2 * a**2 * inner_exp(t,F0,a,b,c,t0) - a**2 * inner_exp(t,F0,a,b,c,t0)**2)/c

def radius_taylor_series_simple(t,a,b,c,t0,F0):
    return ((-a**2 + b)/c) - (((-a**2 + b)/c) - r0) * np.exp(c/(2*a) * t)

def radius_polynom(t,a,b,c,t0,F0):
    return (a + np.sqrt(a**2 + k_nutrient * Lp * R * T))*(t - t0) + r0

def volume_exponential(t,k,V0):
    return V0 * np.exp(k*t)

def volume_linear(t,k,V0):
    return k * t + V0

time = np.linspace(0,24000,24000)
r_t = map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time)

#r_t = map(lambda foo: radius(foo,[F0,a,b,c,t0]),time)
b = set_b(Lp,pi_t,pi_e,k_nutrient,R,T)
c = set_c(k_deg,Lp,R ,T)
F0 = set_F0(a,b,c,r0)
r_t = np.array(map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
r_taylor_series_t = np.array(map(lambda foo: radius_taylor_series(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
r_taylor_series_simple_t = np.array(map(lambda foo: radius_taylor_series_simple(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
r_polynom_t = np.array(map(lambda foo: radius_polynom(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))

V0_exp = 0.3
k_exp  = 0.001
V0_lin = 0.3
k_lin  = 0.003
V_exp_t = np.array(map(lambda foo: volume_exponential(foo,**{'V0': V0_exp,'k':k_exp}),time))
V_lin_t = np.array(map(lambda foo: volume_linear(foo,**{'V0': V0_lin,'k': k_lin}),time))

fig = plt.figure(figsize=(8,8))

ax1 = plt.subplot(2,1,1)
plt.gca().add_artist(FT.figure_numbering('A'))
plt.plot(time/60,4./3. * np.pi * (r_t)**3,'-',label='analytical',color='grey')
plt.plot(time/60,4./3. * np.pi * (r_taylor_series_simple_t)**3,':',linewidth=3.,label='approx. analytical',color='red')
plt.axhline(0,color = 'grey',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
#plt.ylim(-5,110)
plt.legend()

ax1 = plt.subplot(2,1,2)
plt.gca().add_artist(FT.figure_numbering('B'))
plt.plot(time/60,4./3. * np.pi * (r_t)**3 - (4./3. * np.pi * (r_taylor_series_simple_t)**3),label='diff',color='red')
plt.axhline(0,color = 'grey',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
#plt.ylim(-5,110)
plt.legend()

#plt.savefig('figures/analytical_solution_supplement.png',dpi=300)
#plt.savefig('figures/analytical_solution_supplement.eps')


#################################################
################## 


fig = plt.figure(figsize=(7,6))




plt.subplot(gs[0:2,0:5])

#plt.gca().add_artist(FT.figure_numbering('A'))
plt.text(0,0,'A',size=18)
from PIL import Image
im = Image.open('figures/sketch_volume_model.svg.png')
plt.imshow(im)# ,extent=[-.2, .8, 0, 1]
plt.axis('off')


ax1 = plt.subplot(gs[2,0:2])
plt.gca().add_artist(FT.figure_numbering('B'))



#plt.plot(time/60,4./3. * np.pi * (r_t)**3,label='analytical',color='grey')
#plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3,':',color = sns.color_palette()[2],linewidth=4.,label='b')
plt.plot(time/60,4./3. * np.pi * (r_taylor_series_simple_t)**3,label='approx. analytical')
#plt.plot(time/60,4./3. * np.pi * r_polynom_t**3,label='c')
plt.plot(time/60,V_exp_t,':',label='exponential: %s exp(%s t)' %(V0_exp,k_exp))
plt.plot(time/60,V_lin_t,':',label='linear: %s t + %s' %(k_lin,V0_lin))


#plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3 - 4./3. * np.pi * r_t**3,color = 'red',label='d')

plt.axhline(0,color = 'grey',linestyle=':')
r_inf = (-a**2 + b)/c
plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylim(-5,110)
#plt.legend(loc='upper right',ncol=2)
plt.legend(loc='upper right',fontsize=8,frameon=True)


#vec_of_scaling_both_k = [0.5,1.,1.5]
vec_of_scaling_both_k = [0.8,1.,1.2]
vec_of_k_scaling_factor = [1./0.75,1./0.5]
#vec_of_k_scaling_factor = [1.,1./0.75]

colors = sns.color_palette()

plt.subplot(gs[2,2:4])
plt.gca().add_artist(FT.figure_numbering('C'))

handles=[]
labels=[]
for i,some_scaling_factor_for_both in enumerate(reversed(vec_of_scaling_both_k)):
    new_k_deg = some_scaling_factor_for_both * k_nutrient
    new_k_nutrient = some_scaling_factor_for_both * k_nutrient  
    b = set_b(Lp,pi_t,pi_e,new_k_nutrient,R,T)
    c = set_c(new_k_deg,Lp,R ,T)
    F0 = set_F0(a,b,c,r0)
    r_t = np.array(map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
    r_polynom_t = map(lambda foo: radius_polynom(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time)
    

    #label=r'$\frac{k_{upt.}}{k_{cons.}} = %s \mu m$; '%str(1.) + r'$k_{scaling} = %s $'%(str(some_scaling_factor_for_both))
    label=r'$ %s; \quad %s $'%(str(1.),str(new_k_nutrient))
    labels.append(label)
    #if some_scaling_factor_for_both != 1.0:
    #    label=r'$%s \cdot k_{upt.}$; $%s \cdot k_{cons.}$'%(some_scaling_factor_for_both,some_scaling_factor_for_both)
    #else:
    #    label=r'$k_{upt.}$ = %s; $k_{cons.}$ = %s'%(k_nutrient,k_deg)
    handle, = plt.plot(time/60,4./3. * np.pi * r_t**3, alpha = 1-(i*.3), color = 'green',label=label)
    handles.append(handle)
    #plt.plot(time/60,r_polynom_t)
    r_inf = (-a**2 + b)/c
    plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
#plt.ylim(ymax=80)
#plt.legend()


plt.subplot(gs[2,2:4])
plt.gca().add_artist(FT.figure_numbering('C'))
for i,some_k in enumerate(vec_of_k_scaling_factor):
    new_k_scaling_factor = some_k
    new_k_deg = new_k_scaling_factor * k_nutrient
    b = set_b(Lp,pi_t,pi_e,k_nutrient,R,T)
    c = set_c(new_k_deg,Lp,R ,T)
    F0 = set_F0(a,b,c,r0)
    r_t = np.array(map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
    r_polynom_t = map(lambda foo: radius_polynom(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time)
    label =r'$ %s; \quad %s $'%(str(1./(new_k_scaling_factor)),'')
    labels.append(label)
    handle, = plt.plot(time/60,4./3. * np.pi * r_t**3, color = colors[i],label=label)
    handles.append(handle)
    #plt.plot(time/60,r_polynom_t)
    r_inf = (-a**2 + b)/c
    plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
plt.ylim(-5,130)

#plt.legend(loc=(0.5,0.45))
#plt.legend(loc=(0.01,0.6),fontsize=7)

leg = plt.legend(handles, labels, 
             loc=(.85,-.1),
             #title=r'$    k_{u/c}    k_{upt.}$' + '\n' + r'$\left[ \mu m \right] \quad \left[ mM/s/um^2 \right] $',  
             ncol=1, numpoints=1, handletextpad=0.6, handlelength = 1.2, fontsize = 10)
title=r'$ \qquad \qquad \qquad k_{u/c}  \quad \quad  k_{upt.}$' + '\n' + r'$ \qquad \qquad \qquad \left[ \mu m \right] \quad \left[ nmol/s/m^2 \right] $'
leg.set_title(title,prop={'size':8})

#leg = plt.legend(loc=(0.01,0.6),title=r'$k_{u/c}\ \left[ \mu m \right] \quad k_{upt.} \ \left[ mM/s/um^2 \right] $')
#leg = plt.legend(loc=(1.1,-.1),fontsize=10,title=r'$ \qquad k_{u/c} \quad k_{upt.}$' + '\n' + r'$\left[ \mu m \right] \quad \left[ mM/s/um^2 \right] $')
#leg._legend_box.align = "right"

############################################################################################
# ## Fancy table legend:
# from matplotlib.patches import Rectangle
# extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

# #Create organized list containing all handles for table. Extra represent empty space
# legend_handle = handles #[extra, extra, extra, extra, extra, im1, im2, im3, extra, im4, im5, im6, extra, im7, im8, im9]

# #Define the labels
# label_row_1 = [r"$f_{i,j}$", r"$i = 1$", r"$i = 2$", r"$i = 3$"]
# label_j_1 = [r"$j = 1$"]
# label_j_2 = [r"$j = 2$"]
# label_j_3 = [r"$j = 3$"]
# label_empty = [""]

# #organize labels for table construction
# legend_labels = np.concatenate([label_row_1, label_j_1, label_empty * 3, label_j_2, label_empty * 3, label_j_3, label_empty * 3])

# #Create legend
# plt.legend(legend_handle, legend_labels, 
#           loc = 9, ncol = 4, shadow = True, handletextpad = -2)
############################################################################################
######################################

plt.tight_layout()
plt.savefig('figures/analytical_solution.png',dpi=300)
plt.savefig('figures/analytical_solution.eps')


plt.show()

# # 'Doubling time' is depending on r0:
# def doubling_time(r0,a,b,c):
#     return (set_F0(a,b,c,2.*r0)-set_F0(a,b,c,r0))



# plt.plot(np.linspace(0.1,1,100),np.array(map(lambda foo: doubling_time(foo,**{'a':a,'b':b,'c':c}),np.linspace(.1,1.,100)))/60.)
# plt.ylabel("Doubling time [minutes]")
# plt.xlabel(r"Initial radius $[\mu m]$")




# plt.plot(np.linspace(0,20000,20000),np.array(map(lambda foo: inner_exp(foo,**{'a':a,'b':b,'c':c,'t0':t0,'F0':F0}),np.linspace(0,20000,20000))))
# plt.show()
