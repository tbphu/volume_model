from scipy.special import lambertw
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import FigureTools as FT
sns.set_style("ticks")

import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(4, 4)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# constants
R = 8.314; # J/mol/K
T = 303; # K 30C

# nutrient uptake and consumptions
k_nutrient = 2.e-1 #3.5e-15; # mM/s/um^2 
#k_nutrient = 3.5e-15; # mM/s/um^2 
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

V_exp_t = np.array(map(lambda foo: volume_exponential(foo,**{'V0': 0.3,'k':0.001}),time))
V_lin_t = np.array(map(lambda foo: volume_linear(foo,**{'V0': .3,'k':0.003}),time))

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



fig = plt.figure(figsize=(8,8))

ax1 = plt.subplot(gs[2:4,0:2])
plt.gca().add_artist(FT.figure_numbering('C'))



#plt.plot(time/60,4./3. * np.pi * (r_t)**3,label='analytical',color='grey')
#plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3,':',color = sns.color_palette()[2],linewidth=4.,label='b')
plt.plot(time/60,4./3. * np.pi * (r_taylor_series_simple_t)**3,label='approx. analytical')
#plt.plot(time/60,4./3. * np.pi * r_polynom_t**3,label='c')
plt.plot(time/60,V_exp_t,':',label='exponential')
plt.plot(time/60,V_lin_t,':',label='linear')


#plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3 - 4./3. * np.pi * r_t**3,color = 'red',label='d')

plt.axhline(0,color = 'grey',linestyle=':')
r_inf = (-a**2 + b)/c
plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylim(-5,110)
#plt.legend(loc='upper right',ncol=2)
plt.legend()

#plt.savefig('figures/analytical_solution_supplement.png',dpi=300)
#plt.savefig('figures/analytical_solution_supplement.eps')

#print a,b,c,np.exp(-1./a * np.sqrt(b-c*r0)-1),(-1./a * np.sqrt(b-c*r0))

#plt.figure(figsize=(8,8))

#plt.plot(time/60,4./3. * np.pi * r_taylor_series_simple_t**3 - 4./3. * np.pi * r_t**3,label='error')

#plt.show()


# ax1 = plt.subplot(gs[2:3,0:2])
# plt.gca().add_artist(FT.figure_numbering('C'))

# #k_nutrient = 0.1
# #k_deg = 0.1
# b = set_b(Lp,pi_t,pi_e,k_nutrient,R,T)
# c = set_c(k_deg,Lp,R ,T)
# F0 = set_F0(a,b,c,r0)
# r_t = np.array(map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
# r_taylor_series_t = np.array(map(lambda foo: radius_taylor_series(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
# r_taylor_series_simple_t = np.array(map(lambda foo: radius_taylor_series_simple(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
# r_polynom_t = np.array(map(lambda foo: radius_polynom(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))


# plt.plot(time/60,4./3. * np.pi * r_t**3,color = sns.color_palette()[2],label='a')
# plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3,':',color = sns.color_palette()[2],linewidth=4.,label='b')
# #plt.plot(time/60,4./3. * np.pi * r_taylor_series_simple_t**3,label='simple')

# plt.plot(time/60,4./3. * np.pi * r_polynom_t**3,label='c')

# plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3 - 4./3. * np.pi * r_t**3,color = 'red',label='d')

# plt.axhline(0,color = 'grey',linestyle=':')
# r_inf = (-a**2 + b)/c
# plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
# plt.ylabel("Total volume [fL]")
# plt.xlabel("Time [minutes]")
# #plt.setp(ax1.get_xticklabels(), visible=False)
# plt.ylim(-5,120)
# plt.legend(loc='upper right',ncol=2)

# ax2 = plt.subplot(gs[3:4,0:2])
# #ax2 = plt.subplot(gs[3:4,0:2],sharex=ax1)
# plt.gca().add_artist(FT.figure_numbering('D'))
# #ax2 = inset_axes(ax, width=2, height=1, loc = 3, bbox_to_anchor=(0.2, 0.1),bbox_transform=ax.figure.transFigure)
# plt.plot(time/60,4./3. * np.pi * r_taylor_series_t**3 - 4./3. * np.pi * r_t**3,color = 'red',label='Difference w/o Taylor series')
# plt.axhline(0,color = 'grey',linestyle=':')
# plt.ylim(-0.1,0.1)
# plt.ylabel("Difference [fL]")
# plt.xlabel("Time [minutes]")
# plt.legend()


vec_of_scaling_both_k = [0.5,1.,1.5]
vec_of_k_scaling_factor = [1.,1./0.75,1./0.5]

plt.subplot(gs[0:2,2:])
plt.gca().add_artist(FT.figure_numbering('B'))
for some_scaling_factor_for_both in reversed(vec_of_scaling_both_k):
    new_k_deg = some_scaling_factor_for_both * k_scaling_factor * k_nutrient
    new_k_nutrient = some_scaling_factor_for_both * k_nutrient  
    b = set_b(Lp,pi_t,pi_e,new_k_nutrient,R,T)
    c = set_c(new_k_deg,Lp,R ,T)
    F0 = set_F0(a,b,c,r0)
    r_t = np.array(map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
    r_polynom_t = map(lambda foo: radius_polynom(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time)
    
    if some_scaling_factor_for_both != 1.0:
        label=r'$%s \cdot k_{upt.}$; $%s \cdot k_{cons.}$'%(some_scaling_factor_for_both,some_scaling_factor_for_both)
    else:
        label=r'$k_{upt.}$ = %s; $k_{cons.}$ = %s'%(k_nutrient,k_deg)
    plt.plot(time/60,4./3. * np.pi * r_t**3,label=label)
    #plt.plot(time/60,r_polynom_t)
    r_inf = (-a**2 + b)/c
    plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")
#plt.ylim(ymax=80)
plt.legend()



plt.subplot(gs[2:,2:4])
plt.gca().add_artist(FT.figure_numbering('D'))
for some_k in vec_of_k_scaling_factor:
    new_k_scaling_factor = some_k
    new_k_deg = new_k_scaling_factor * k_nutrient
    b = set_b(Lp,pi_t,pi_e,k_nutrient,R,T)
    c = set_c(new_k_deg,Lp,R ,T)
    F0 = set_F0(a,b,c,r0)
    r_t = np.array(map(lambda foo: radius(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time))
    r_polynom_t = map(lambda foo: radius_polynom(foo,**{'F0':F0,'a':a,'b':b,'c':c,'t0':t0}),time)
    plt.plot(time/60,4./3. * np.pi * r_t**3,label=r'$\frac{k_{upt.}}{k_{cons.}} = %s$'%str(1./new_k_scaling_factor))
    #plt.plot(time/60,r_polynom_t)
    r_inf = (-a**2 + b)/c
    plt.axhline((4./3. * np.pi * r_inf**3),color = 'black',linestyle=':')
plt.ylabel("Total volume [fL]")
plt.xlabel("Time [minutes]")

plt.legend(loc=(0.6,0.45))




plt.subplot(gs[0:2,0:2])
plt.gca().add_artist(FT.figure_numbering('A'))
from PIL import Image
im = Image.open('figures/sketch_volume_model.svg.png')
plt.imshow(im)
plt.axis('off')

plt.tight_layout()
#plt.savefig('figures/analytical_solution.png',dpi=300)
#plt.savefig('figures/analytical_solution.eps')


plt.show()

# # 'Doubling time' is depending on r0:
# def doubling_time(r0,a,b,c):
#     return (set_F0(a,b,c,2.*r0)-set_F0(a,b,c,r0))



# plt.plot(np.linspace(0.1,1,100),np.array(map(lambda foo: doubling_time(foo,**{'a':a,'b':b,'c':c}),np.linspace(.1,1.,100)))/60.)
# plt.ylabel("Doubling time [minutes]")
# plt.xlabel(r"Initial radius $[\mu m]$")




# plt.plot(np.linspace(0,20000,20000),np.array(map(lambda foo: inner_exp(foo,**{'a':a,'b':b,'c':c,'t0':t0,'F0':F0}),np.linspace(0,20000,20000))))
# plt.show()
