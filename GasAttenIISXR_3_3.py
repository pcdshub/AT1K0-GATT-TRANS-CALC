# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 21:59:18 2012

@author: schaferd
"""
name='GasAttenIISXR_3_3'
#3_3 is adding another pumping stage ~5/28/15
#3_2 is adding Yiping's steady state solution-notebook ~9/25/14
#3_1 is adding thermal diffusion calc (Milchberg)-notebook ~8/8/14
#3_0 is splitting input file out, moving speed function S to pump_data_3_0_0 and gas choice to gas_data_3_0_0
#2_2 is putting pump speed and pressure and beamline speed and pressure calcs in
#2_1 is changing pump effective speed (notebook 11/12/13)
#GasAttenIIHXR_2_0 is Super Conducting Linac on HXR 
#GasAttenIISXR_2_0 is Super Conducting Linac on SXR 
#2_0-modify for LCLSII replan
#1_6-add orifice calculations for choked flow (Euler)
#1.5 based on 1.4 - roughing pump and foreline design.
#1.4 is based on 1.3a, adding the beam stay clear in Michael Rowen email 1/29/13
#stay clear(w/o steering error)=z[m]*25.4*photon energy[eV]^-(5/6)
#stay clear(w/steering error)=z[m]*(25.4*photon energy[eV]^-(5/6)+0.01)+1
#add min_energy variable as lowest beam energy 
#GasAttenII_1.4 - 2/4/13
#predecessor in SAGE is GasAtten5-LCLS-I-8.10.12


import matplotlib.pyplot as plt
import numpy
import math
import sympy as sp
import csv
from scipy.integrate import simps, trapz
from matplotlib.ticker import FormatStrFormatter

#import gas_data_3_0_0 as gd
from gas_data_3_0_0 import *
from pump_data_3_0_0 import *
from PfeifferEVR116 import *

#input file: 
from GasAttenIISXR_3_3_0_in import *

#DEFINE FUNCTIONS
    

#check if number is positive, if negative return 0
def is_pos(value):
    if value < 0:
        value = 0
    return value

#D_beamsize    
def D_beamsize(a_keV):
    D_beamsize=(1000*F_opr*B_FWHM(a_keV))*(z_adjaper-(z_e+z_0(a_keV)))  #XTES PRD (LCLSII-3.5-PR-0051-R1) eqn (1)-notebook 1/20/15-using B_FWHM from Beam source parameters for LCLS II Krzywinski_ypf.xlsx     
    #D_beamsize=(1000*F_opr*B_FWHM(a_keV)*(a_keV/E_0)**(n_exp))*(z_adjaper-(z_e+z_0(a_keV)))  #XTES PRD (LCLSII-3.5-PR-0051-R1) eqn (1)
    #D_beamsize=(1000*2*38e-6*(a_keV*1000./250.)**(-1))*(z_adjaper-source)  #XTES PRD (LCLSII-3.5-PR-0051-R1) eqn (1)
    #D_beamsize=2*1000*(z_adjaper-source)*math.tan(11/(a_keV*2e6))    #FWHM, see notebook 11/6/13
    #D_beamsize=(2*63.1e-6*(a_keV*1000./250.)**(-5./6.))*(z_adjaper-source)*1000   #PRD R002, eqn (1)
    return D_beamsize

#D_aperture
def D_aperture(a_keV):
    D_aperture=F_diff*(1000*F_opr*B_FWHM(a_keV))*(z_adjaper-(z_e+z_0(a_keV)))  #XTES PRD eqn (2)-notebook 1/20/15-using B_FWHM from Beam source parameters for LCLS II Krzywinski_ypf.xlsx    
    #D_aperture=F_diff*(1000*F_opr*B_FWHM(a_keV)*(a_keV/E_0)**(n_exp))*(z_adjaper-(z_e+z_0(a_keV)))  #XTES PRD eqn (2)
    #D_aperture=2*(1000*2*38e-6*(a_keV*1000./250.)**(-1))*(z_adjaper-source)  #XTES PRD eqn (2)    
    #D_aperture=4*2*1000*(z_adjaper-source)*math.tan(11/(a_keV*2e6))    #4*FWHM, see notebook 11/6/13    
    #D_aperture=2.*(2*63.1e-6*(a_keV*1000./250.)**(-5./6.))*(z_adjaper-source)*1000   #PRD R002, eqn (2)  
    #D_aperture=z_loc*25.4*(a_keV*1000)**(-5./6)    #Rowen email 1/29/13-stay clear w/o steering
    #D_aperture=5.5/a_keV             #value based on email 1/14/13, new stay clear    
    #D_aperture=3.60604/a_keV        #value based on last 3 unds not installed in baseline(notebook 11/6/12), also at end of device, not center    
    #D_aperture=2.8469/a_keV          #value based on released PRD
    #D_aperture=3.3922/a_keV         #value based on old PRD
    return D_aperture

    
#D_stayclear
def D_stayclear(a_keV):
    D_stayclear=max(D_min,F_diff*(1000*F_opr*B_FWHM(a_keV))*(z_adjaper-(z_e+z_0(a_keV)))+1000*B_steer*1e-6*(z_adjaper-z_e)+D_offset)    #XTES PRD eqn (3)-notebook 1/20/15-using B_FWHM from Beam source parameters for LCLS II Krzywinski_ypf.xlsx    
    #D_stayclear=max(D_min,F_diff*(1000*F_opr*B_FWHM(a_keV)*(a_keV/E_0)**(n_exp))*(z_adjaper-(z_e+z_0(a_keV)))+1000*B_steer*1e-6*(z_adjaper-z_e)+D_offset)    #XTES PRD eqn (3)
    #D_stayclear=max(5,2*(1000*2*38e-6*(a_keV*1000./250.)**(-1))*(z_adjaper-source)+1000*10e-6*(z_adjaper-und_end)+1)    #XTES PRD eqn (3)
    #D_stayclear=D_aperture(a_keV)+10e-6*(z_end-und_end)*1000+1    
    #D_stayclear=D_aperture(a_keV)+10e-6*(z_adjaper-und_end)*1000+1              #PRD R002 eqn (3)
    return D_stayclear

def closest(target, collection):
    return min((abs(target - i), i) for i in collection)[1]
            


#find closest aperture dia from choices, if it is >= beam dia then use, if < use the next larger
def ap_pick(m,apXdia):   
    if m<=min_energy:    
        beam_dia=D_aperture(min_energy)
        #print m,"m<=min_energy",min_energy
    if m>min_energy:
        beam_dia=D_aperture(m)
        #print m,"m>min_energy",min_energy
    #print "beam_dia", beam_dia
    ap_pick=closest(beam_dia,apXdia)
    #if only one diameter, use that dia    
    if len(apXdia)==1:
        ap_pick=apXdia[0]
        #print "ap_pick1", ap_pick
        
    if ap_pick>=beam_dia:
        bogus=1
        #print "if ap_pick",ap_pick,">=beam_dia",beam_dia
 
    if len(apXdia)!=1 and ap_pick<beam_dia:
        ai=apXdia.index(ap_pick)
        ap_pick=apXdia[ai+1]
        #print "ai", ai,"beam dia",beam_dia,"ap_pick",ap_pick,apXdia[ai]        
        #print "else ap_pick=", ap_pick,"[mm]"
           
    return ap_pick

# free-mol. throughput (eq. 4) see notebook 1/7/13 - output is pressure[units]*[l/s]
def qfm(pressure,aperture): 
    qfm=1.e-3*((3.1415926)**.5*T0/(8.*T))*pressure*vm*aperture**2 
    return qfm 
                
# rarefaction parameter (eq. 3), delta(pressure[Torr],aperture[mm])
def delta(pressure,aperture):
    delta=pressure*133.32*(aperture/2.e3)/(amu*vm)
    return delta

# reduce flow rate
def W(de,ap_type):    
    
    #plate apertures
    if ap_type=="plate":        
        #from 'Final Report', Dr. Sharipov, FORTRAN Code, curve fit of P1/P0=0 in Fig 2. 
        #'Numerical simulation of rarefied gas flow through a thin orifice'        
        W=1.+de*(0.000555*de**2+0.0375*de+0.1468)/(1.92*0.000555*de**3+0.06675*de**2+0.3722*de+1.)      

    #tube apertures, with L/D=5
    if ap_type=="tube":        
        W=numpy.interp(de,L_D5_de,L_D5_W)
        
    return W
    
#integrate using Simpsons rule, h is x offset step    
def integrate(y_vals, h):
    i=1
    total=y_vals[0]+y_vals[-1]
    for y in y_vals[1:-1]:
        if i%2 == 0:
            total+=2*y
        else:
            total+=4*y
        i+=1
    return total*(h/3.0)
    
############################
#make logarithmic count list L4_1 so pumps modeled with equations can be plotted.

L4_1=[]
L5_1=[]
L6_1=[]
s_mult=[]

L1=[-8,-7,-6,-5,-4,-3,-2,-1,0]
L2=[10**p for p in L1]
L3=[p for p in xrange(1,10,1)]
L6=[.01*h for h in xrange(20,250,1)]


for i in L2:
    for j in L3:
        L4=[i*j]
        L4_1=L4_1+L4
#print L4_1       
#run L4_1 logarithmic number list through the pump equation and make y list. 
def get_y(function,list_used):
    L5_1=[]    
    for i in list_used:
        L5=[function(i)]
        L5_1=L5_1+L5
    return L5_1

#need list to plot plate reduction in flow rate (W) from equation
L7=[10**i for i in range(-6,5,1)]
L8=[W(i,"plate") for i in L7]

#multiply pump speed by the number of pumps
def s_mult(mult,pump_choice):
    s_mult=[]    
    for i in L4_1:
        s_mult1=[mult*S(i,pump_choice,gas)[0]]
        s_mult=s_mult+s_mult1
    return s_mult
    
##############
#chooses x and y values from main output lists
def x_list(list_4):
    x_list=[list_4[i:i+1][0] for i in xrange(0,len(list_4),3)]
    return x_list

def y_list(list_4):
    y_list=[list_4[i+1:i+2][0] for i in xrange(0,len(list_4),3)]
    return y_list
            
###########
#Transmission pump pipe probabilty (alpha) data in "pump_data.py", from Table 5.1 "Modern Vacuum Physics"
def trans_prob(pxL_D):
    pumpL_D=pxL_D[0]/pxL_D[1]
    trans_prob=numpy.interp(pumpL_D,L_D,alpha_c)
    return trans_prob         

###########
#Conductance of an aperture, used for short pipes, see eq 5.38 "Mod Vac Phys"
def C_a(dia):
    C_a=(math.pi/4)*1000*(dia**2.0)*(Rg*T/(2*math.pi*am))**.5
    return C_a
    
########### calculate pressure, iterate eq. 9 to get correct pump speed and pressure
#"p0=",p0,"W0",W0,"aperture0",aperture0,"[mm]","# pumps",stage_1[0],"pump",stage_1[1],"pump1 L/D",p1LD
def pcalc(p_c,W_c,ap_c,N_c,S_c,pxL_D,qfm_ap_c):
    
    i=0
    pconverge=[]
    #init pump speed
    #x=at beam
    P=trans_prob(pxL_D)    
    C_a1=C_a(pxL_D[1])*P    #eq. 5.41 Mod Vac Phys - for molecular flow...  
    S_eff=N_c*S(p_c/100,S_c,gas)[0]*C_a1/(S(p_c/100,S_c,gas)[0]+C_a1)
    #print "pxLD",pxL_D,"trans_prob",P,"S_eff",S_eff
    #init pcalc, eq.9
    #pcalc_x=(((p_c*(math.pi)**.5*T0*W_c*(ap_c)**2*vm)/(8*T*1000))+q_outgas)/(S_eff)
    pcalc_x=(qfm_ap_c*W_c+q_outgas)/(S_eff)    
    while i<25:
        S_x_a=N_c*S(pcalc_x,S_c,gas)[0]*C_a1/(S(pcalc_x,S_c,gas)[0]+C_a1)
        #(eq. 9)
        #pcalc_x=(((p_c*(math.pi)**.5*T0*W_c*(ap_c)**2*vm)/(8*T*1000))+q_outgas)/(S_eff)
        pcalc_x=(qfm_ap_c*W_c+q_outgas)/(S_eff)        
        S_x_b=N_c*S(pcalc_x,S_c,gas)[0]*C_a1/(S(pcalc_x,S_c,gas)[0]+C_a1)
        S_eff=(S_x_a+S_x_b)/2.
        pconverge.append(pcalc_x)
        if i>=1:        
            #print i,"conv",abs((pconverge[-1]-pconverge[-2])/pconverge[-1]),pconverge[-1],pconverge[-2]         
            if abs((pconverge[-1]-pconverge[-2])/pconverge[-1])<.001:
                #print pconverge[-1],pconverge[-2]                
                i=25
        i=i+1
    if pcalc_x < S(pcalc_x,S_c,gas)[1]:
        pcalc_x = S(pcalc_x,S_c,gas)[1]
    pcalc=[pcalc_x,S_eff]
    #print "pcalc",pcalc
    #print "pconverge",pconverge
    return pcalc

########### calculate pressure to get correct pump speed and pressure
#"p0=",p0,"W0",W0,"aperture0",aperture0,"[mm]","# pumps",stage_1[0],"pump",stage_1[1],"pump1 L/D",p1LD
def pzcalc(p_c,ap_c,ap_t,N_c,S_c,pxL_D):
    
    i=0
    pconverge=[]
    #init pump speed
    #x=at beam
    P=trans_prob(pxL_D)         #alpha
    C_a1=P*C_a(pxL_D[1])      #conductance of pump tube, eq (5.41) Mod Vac Phys
    if ap_t == 'plate':    
        C_a2=C_a(ap_c/1000)     #conductance of aperture, conv ap_c[mm] to [m] 
    if ap_t == 'plate_sq':
        C_a2=(4/math.pi)*C_a(ap_c/1000)     #conductance of aperture, conv ap_c[mm] to [m], 4/pi is difference in area between circle and square 
    if ap_t == 'tube':
        C_a2=trans_prob([5.0,1.0])*C_a(ap_c/1000)     #Mod Vac Phys eq (5.42), tube L/D=5.0
    if ap_t == 'tube_sq':
        C_a2=trans_prob([5.0,1.0])*(4/math.pi)*C_a(ap_c/1000)     #Mod Vac Phys eq (5.42), tube L/D=5.0    
    S_eff_tot=1/(1/(N_c*S(p_c/100,S_c,gas)[0])+1/(N_c*C_a1)+1/C_a2)    #eff speed from aper to pump
    S_eff_beam=1/(1/(N_c*S(p_c/100,S_c,gas)[0])+1/(N_c*C_a1))           #eff speed from beam to pump   
    S_pump=N_c*S(p_c/100,S_c,gas)[0]   
    #print "pxLD",pxL_D,"trans_prob",P,"C_a1",C_a1,"C_a2",C_a2,"S_eff_tot",S_eff_tot,"S_eff_pump",S_eff_pump
    #init pcalc, eq.9
    q_aper=S_eff_tot*p_c+q_outgas
    pzcalc_pump=q_aper/S_pump
    pzcalc_beam=q_aper/S_eff_beam
    while i<25:
        #S_x_a=N_c*S(pzcalc_beam,S_c)*C_a1*P/(S(pzcalc_beam,S_c)+P*C_a1)
        S_eff_tot_a=1/(1/S(pzcalc_beam,S_c,gas)[0]+1/(N_c*C_a1)+1/C_a2)    #eff speed from aper to pump
        S_eff_beam_a=1/(1/S(pzcalc_beam,S_c,gas)[0]+1/(N_c*C_a1))           #eff speed from beam to pump        
        S_pump_a=N_c*S(pzcalc_pump,S_c,gas)[0]        
        q_aper=S_eff_tot*p_c+q_outgas
        pzcalc_pump=q_aper/S_pump
        pzcalc_beam=q_aper/S_eff_beam
        #S_x_b=N_c*S(pzcalc_beam,S_c)*C_a1*P/(S(pzcalc_beam,S_c)+P*C_a1)
        S_eff_tot_b=1/(1/S(pzcalc_beam,S_c,gas)[0]+1/(N_c*C_a1)+1/C_a2)    #eff speed from aper to pump
        S_eff_beam_b=1/(1/S(pzcalc_beam,S_c,gas)[0]+1/(N_c*C_a1))           #eff speed from beam to pump         
        S_pump_b=N_c*S(pzcalc_pump,S_c,gas)[0]
        #S_eff=(S_x_a+S_x_b)/2.
        S_eff_tot=(S_eff_tot_a+S_eff_tot_b)/2
        S_eff_beam=(S_eff_beam_a+S_eff_beam_b)/2
        S_pump=(S_pump_a+S_pump_b)/2
        pconverge.append(pzcalc_beam)
        if i>=1:        
            #print i,"pzconv",abs((pconverge[-1]-pconverge[-2])/pconverge[-1]),pconverge[-1],pconverge[-2]         
            if abs((pconverge[-1]-pconverge[-2])/pconverge[-1])<.001:
                #print pconverge[-1],pconverge[-2]                
                i=25
        i=i+1
    if pzcalc_pump < S(pzcalc_beam,S_c,gas)[1]:
        pzcalc_beam = S(pzcalc_beam,S_c,gas)[1]    
    pzcalc=[pzcalc_beam,pzcalc_pump,q_aper,C_a1,S_eff_tot,S_eff_beam,S_pump]
    #print "pzcalc",pzcalc,p_c,ap_c,N_c,S_c,pxL_D
    #print "pconverge",pconverge
    return pzcalc

################# calc rough line pressure
def rcalc(p_c,Q_x,N_c,S_c,trans):
    i=0
    pconverge=[]
    sconv=[]
    #init pump speed
    S_x=N_c*S(.05,S_c,gas)[0]
    #init rcalc, P=Q/S
    rcalc_x=Q_x/(S_x*trans)
    if rcalc_x<.01:
        rcalc_x=.01
    #print "rcalc init",rcalc_x,"q",Q_x,"S_x",S_x,N_c,S_c,trans
    while i<100:
        S_x_a=N_c*S(rcalc_x,S_c,gas)[0]
        rcalc_x=Q_x/(S_x*trans)
        #base pressure of rough pump        
        if i>1 and rcalc_x<.01:        #[Torr]
            rcalc_x=.01                 
        S_x_b=N_c*S(rcalc_x,S_c,gas)[0]
        S_x=(S_x_a+S_x_b)/2.        
        pconverge.append(rcalc_x)
        sconv.append(S_x)
        S_x=math.fsum(sconv)/len(sconv)        
        if i>=1:        
            #print i,"sxa",S_x_a,"sxb",S_x_b,"sx",S_x,"pconv",pconverge[-1],pconverge[-2]         
            if abs((pconverge[-1]-pconverge[-2])/pconverge[-1])<.001:
                #print pconverge[-1],pconverge[-2]                
                i=100           
        i=i+1
    #rcalc=[rcalc_x,S_x*trans]
    rcalc=[rcalc_x,pconverge]
    #print "rcalc",rcalc
    return rcalc                
#################
def foreline_max(p):    
    foreline_max=[i*1e6 for i in p]     #multiplier is compression ratio or turbo
    #print foreline_max      
    i=0    
    while i<len(foreline_max):
        if foreline_max[i]>1:
        #print foreline_max(i)
            foreline_max[i]=1
        i=i+1
    return foreline_max

################
def main_calc(e_list,pressure,aperture,ap_type,num_pump,pump,pumpL_D):    
    apdia=[ap_pick(self,aperture) for self in e_list]
    #print "apdia", apdia, "pressure",pressure    
    de=[self1*133.32*(self2/2.e3)/(amu*vm) for self1,self2 in zip(pressure,apdia)]
    
    #plate apertures-round hole
    if ap_type=="plate":        
        W_ap=[1.+self*(0.000555*self**2+0.0375*self+0.1468)/(1.92*0.000555*self**3+0.06675*self**2+0.3722*self+1.) for self in de]
        qfm_ap=[1.e-3*((math.pi)**.5*T0/(8.*T))*self1*vm*self2**2 for self1,self2 in zip(pressure,apdia)]

    #plate aperture-square hole
    if ap_type=="plate_sq":        
        W_ap=[1.+self*(0.000555*self**2+0.0375*self+0.1468)/(1.92*0.000555*self**3+0.06675*self**2+0.3722*self+1.) for self in de]
        qfm_ap=[1.e-3*(T0/(2.*T*(math.pi)**.5))*self1*vm*self2**2 for self1,self2 in zip(pressure,apdia)]     #see notebook 5/27/15

    #tube apertures, with L/D=5
    if ap_type=="tube":        
        W_ap=[numpy.interp(self,L_D5_de,L_D5_W) for self in de]
        qfm_ap=[1.e-3*((math.pi)**.5*T0/(8.*T))*self1*vm*self2**2 for self1,self2 in zip(pressure,apdia)]
        
    #square tube apertures, with L/D=5
    if ap_type=="tube_sq":        
        W_ap=[numpy.interp(self,L_D5_de,L_D5_W) for self in de]
        qfm_ap=[1.e-3*(T0/(2.*T*(math.pi)**.5))*self1*vm*self2**2 for self1,self2 in zip(pressure,apdia)]     #see notebook 5/27/15        
        
    #qfm_ap=[1.e-3*((3.1415926)**.5*T0/(8.*T))*self1*vm*self2**2 for self1,self2 in zip(pressure,apdia)] 
    q_ap=[self1*self2 for self1,self2 in zip(W_ap,qfm_ap)] 
    if num_pump==0:
        p_ap=[0,0]
    if num_pump!=0:
        p_ap=[pcalc(self0,self1,self2,self3,self4,self5,self6) 
        for self0,self1,self2,self3,self4,self5,self6 
        in zip(pressure,W_ap,apdia,len(e_list)*[num_pump],[pump]*len(e_list),[pumpL_D]*len(e_list),qfm_ap)]
    main_calc=[apdia,de,W_ap,qfm_ap,q_ap,p_ap] 
    #print main_calc
    return main_calc


###########################################
#begin main

[amu,am,gas_keV,gas_MA,gas_cost,gas_mass,gas_g,gas_c_p,gas_c_v,gas_k,gas_T,gas_Cp2,gas_K]=gas_test(gas,T)

#calculate pump cost
t_cost=2*(stage_1[0]*stage_1[2]+stage_2[0]*stage_2[2]+stage_3[0]*stage_3[2]+
          stage_4[0]*stage_4[2]+stage_5[0]*stage_5[2]+r_2[0]*r_2[2]+
          r_3[0]*r_3[2]+r_4[0]*r_4[2]+r_5[0]*r_5[2])
#print "total pump cost", t_cost
print("total pump cost", t_cost)

vm=(2.*Rg*T/am)**.5                            #probable molecular speed (eq. 3)
    
#generate energy list
e_list=[keV+e_index*self for self in xrange(0,int((keV_end-keV)/e_index)+1,1)]
e_count=range(0,len(e_list),1)

#test beam source location and divergence values
divergence=[B_FWHM(self) for self in e_list]
source_location=[z_0(self) for self in e_list]
dist_to_source=[(z_adjaper-(z_e+z_0(self))) for self in e_list]

#p0=[3.4347*self-.5869 for self in e_list]       #delete - test 6/19/14
p0=[((-(math.log(atten))*Rg*T)/(((numpy.interp(self,gas_keV,gas_MA))/10)*cell_l*am))/133.32 for self in e_list]
p0_Pa=[self*133.32 for self in p0]

#main_calc(e_list,pressure,aperture,ap_type,num_pump,pump,pumpL_D): 
    #return [apdia,de,W_ap,qfm_ap,q_ap,p_ap] 
    
[ap0,de0,W0,qfm0,q0,ps1]=main_calc(e_list,p0,ap0dia,ap0type,stage_1[0],stage_1[1],p1L_D)
p1=[ps1[self][0] for self in e_count]
s1=[ps1[self][1] for self in e_count]

[ap1,de1,W1,qfm1,q1,ps2]=main_calc(e_list,p1,ap1dia,ap1type,stage_2[0],stage_2[1],p2L_D)  
p2=[ps2[self][0] for self in e_count]
s2=[ps2[self][1] for self in e_count]

[ap2,de2,W2,qfm2,q2,ps3]=main_calc(e_list,p2,ap2dia,ap2type,stage_3[0],stage_3[1],p3L_D)
p3=[ps3[self][0] for self in e_count]
s3=[ps3[self][1] for self in e_count]

[ap3,de3,W3,qfm3,q3,ps4]=main_calc(e_list,p3,ap3dia,ap3type,stage_4[0],stage_4[1],p4L_D)
p4=[ps4[self][0] for self in e_count]
s4=[ps4[self][1] for self in e_count]

[ap4,de4,W4,qfm4,q4,ps5]=main_calc(e_list,p4,ap4dia,ap4type,stage_5[0],stage_5[1],p5L_D)
[p5,s5]=[ps5[self][0] for self in e_count],[ps5[self][1] for self in e_count]

#[ap5,de5,W5,qfm5,q5,ps6]=main_calc(e_list,p5,ap5dia,ap5type,0,0,0)
[ap5,de5,W5,qfm5,q5,ps6]=main_calc(e_list,p5,ap5dia,ap5type,stage_6[0],stage_6[1],p6L_D)
[p6,s6]=[ps6[self][0] for self in e_count],[ps6[self][1] for self in e_count]

[ap6,de6,W6,qfm6,q6,ps7]=main_calc(e_list,p6,ap6dia,ap6type,0,0,0)

#p6=[ps6[self][0] for self in range(0,len(e_list),1)]
#s6=[ps6[self][1] for self in range(0,len(e_list),1)]


#use p0 and gauge pressure error to calculate attenuation error. 
p0_plus=[h+p0_error/200*h for h in p0]
p0_minus=[h-p0_error/200*h for h in p0] 
p_diff_1=[x-y for x, y in zip(p0_plus,p0_minus)]
#print p_diff_1
#p_in [Torr], /133.32 to convert [Pa] to [Torr]
#print zip(p0_plus,p0x)
att_diff=[math.exp(((-r*133.3*(numpy.interp(k1,gas_keV,gas_MA))/10)*cell_l*am)/(Rg*T)) for r,k1 in zip(p0_minus,e_list)]
att_error=[w/atten for w in att_diff]
#print att_diff


#calculate total flow rate t_q[mole/hr] and gas cost, see notebook 2/12/13 and 1/29/14
t_q=[2*3600*.1333*self/(Rg*T) for self in q0] 
gcost=[self*gas_cost for self in t_q]

############################################################
#try to figure out losses in first pump connection
#mass flow [kg/s] for q0
Mdot0=[am*.1333*self/Rg*T for self in 2*q0]

#p1range=[x/100. for x in range(1,300)]
qfm1a=[qfm(self,p1p_conn[0]) for self in p1]
de1a=[delta(self,p1p_conn[0]) for self in p1]
W1a=[numpy.interp(self,L_D5_de,L_D5_W) for self in de1a]
q1a=[self*self1 for self,self1 in zip(qfm1a,W1a)]

#"Rarefied gas flow through a long tube at any pressure ratio", Sharipov, J. Vac. Sci. Tech. A 12(5), Sep/Oct 1994 
de_P1,P_P1=sp.symbols('de_P1,P_P1')
de_P1=P_P1*p1p_conn[0]/(2*amu*vm)
de_P1_avg=[(self+de_P1)/2. for self in de1a]
#print 'inputs, de_P1_avg[0],p1[0]
print('inputs',de_P1_avg[0],p1[0])
#eq (11),(12), and (16) to find mass flow rate Mdot:
MdotP1=[((self1/4.)+1.0162+(.5489/self1)-(.6801/self1**2))/((-1/(p1p_conn[0]/2))*
         (p1p_conn[1]/(P_P1-self2))*(2*k*T/gas_mass)**.5) for self1,self2 in zip(de_P1_avg,p1)]

#print len(MdotP1)==len(Mdot0)
#print zip(Mdot0,MdotP1)
#ans=[solve(i-j,P_P1) for i,j in zip(Mdot0,MdotP1)]
print ("MdotP1[0]-Mdot0[0]",MdotP1[0]-Mdot0[0])
ans=sp.solve([MdotP1[10]-Mdot0[10]],[P_P1])
#print ans

area1=integrate(t_q,e_index)

# Compute the area using the composite trapezoidal rule.
area2 = trapz(t_q, dx=e_index)

# Compute the area using the composite Simpson's rule.
area3 = simps(t_q, dx=e_index)

avg=(area1+area2+area3)/3

print ("area under gas flow curve:","area1",area1,"area2",area2,"area3",area3,"avg",avg)

#how effective are orifices?
p_diff10=[i/j for i, j in zip(p1,p0)]
p_diff21=[i/j for i, j in zip(p2,p1)]
p_diff32=[i/j for i, j in zip(p3,p2)]
p_diff43=[i/j for i, j in zip(p4,p3)]
p_diff54=[i/j for i, j in zip(p5,p4)]

S_p1d=[i/(.06*trans_prob(p1L_D)) for i in q0]
S_p1=[i/(.06) for i in q0]
q_diff=[i/j for i,j in zip(q1,q0)]
print ('trans_prob',trans_prob(p1L_D))

#estimate voltage input for Pfeiffer EVR116 gas dosing valve
#need total flow (2 sides). q0y2=q0_tot
q0y2=[i*2 for i in q0]
inj_V=[numpy.interp(i,EVR116_q_air,EVR116_V_air) for i in q0y2]

################################################
#estimate flow rate assuming choked flow using Euler ean.
#eqn (3.2) from Numerical simulation of rarefied flow through a thin orifice" Dr. Sharipov
 
CPR=(2/(gas_g+1))**(gas_g/(gas_g-1))
print ('critical pressure ratio, (p_1/p_0)*=', CPR)

vo=(2*k*T/gas_mass)**.5

Mdot_o=[(((math.pi**.5)*(a/2000)**2)*Po*133.3)/vo for a,Po in zip(ap0,p0)]      #Eq. 2.4, a is radius[m], dia is [mm]/2000=rad[m]

factor=((2*math.pi*gas_g)**.5)*(2/(gas_g+1))**((gas_g+1)/(2*(gas_g-1)))

Mdot_n=[self*factor for self in Mdot_o]     #Eq. 3.2, [kg/s]

print ('factor',factor) #,'Mdot_n[kg/s]', Mdot_n

cell_vol=(cell_l*math.pi*cell_dia**2)/4                #[m^3]
dens=[(am*self*133.3)/(Rg*T) for self in p0]   #density of gas[kg/m^3], 133.3[Pa]=1[Torr], notebook 10/2/13
num_dens=[self*avogadro/am for self in dens]    #number density of gas[atoms/m**3]
gas_mass_in_cell=[self*cell_vol for self in dens]   #[kg] gas mass in gas cell

Mdot_n2=[C*self1*self2*1000/self3 for self1,self2,self3 in zip(Mdot_n,p0,dens)]      #[Torr*l/s]

#######################################
#estimate gas heating 
#from LCLS-TN-09-5, section 3.3 (in notebook~10/7/13, 4/11/14, 4/23/14)

chi=[gas_k/(self*gas_c_p) for self in dens]                             #calculate thermal diffusivity, better approx than eq. 3.22?
sigma=[(numpy.interp(self,gas_keV,gas_MA))/10 for self in e_list]       #[m^2/kg], mass atten coeff[cm^2/g], div by 10 in prog to get [m^2/kg]
q_length=[i*j*Wi for i,j in zip(dens,sigma)]                            #[J/m]  Eq 3.12
p_elstat=[(160*self*(self-.4))/((self2*10)*Wi*1000) for self,self2 in zip(e_list,sigma)]     #eq 3.14 - notebook 4/11/14
beam_r=[(D_beamsize(self))/2000 for self in e_list]                     # beam radius[m]
F_0=[2*Wi/(math.pi*(self)**2) for self in beam_r]                       #[kg/(pulse*s**2)]   eq 3.16
Q_vol=[(eta*self*self1*self2) for self,self1,self2 in zip(dens,sigma,F_0)]        #[J/(pulse*m**3)]  eq 3.17
pwr_dens=['%.2e' %(self*rep_rate) for self in Q_vol]                    #[W/m^3] power density 
del_T=[self1/(gas_c_v*self2) for self1,self2 in zip(Q_vol,dens)]        #temp rise from one pulse[K], similar to eq 3.19
Tau_egas=[2e-5/self for self in p0]                                     #eq. 3.20 energy exchange time[s] b/w e- and gas
r_therm=[self2+(2*self*pulse_period)**.5 for self,self2 in zip(chi,beam_r)]                   #eq. 3.21, with beam_r added
del5_T=[4.2e-4*self1*self2*10*Wi*1000 for self1,self2 in zip(p0,sigma)] #eq. 3.23 [K] on axis temp increase/pulse
del2_T=[self1/(math.pi*gas_c_p*self2*self3**2.0) for self1,self2,self3 in zip(q_length,dens,r_therm)]
#del3_T=[self1/(math.pi*self2**2*gas_c_p*self3) for self1,self2,self3 in zip(q_length,beam_r,dens)]
#[K],bounding temp input, qdot=2*mdot*cp*delT (notebook 10/25/13), Mdot_n is approximation, really should use q0y2 and convert to [kg/s]
del4_T=[P_i/(2*self*gas_c_p) for self in Mdot_n]

###########################################
#gas density hole-thermal diffusion
#"The effect of long timescale gas dynamics on femtosecond filementation". Milchberg
#http://lasermatter.umd.edu/Publications%20PDF/2013_ChengOE.pdf
#also see notebook 8/1-8/8/14
t=[1e-6*self for self in [25,100,200,300,400,600,800,1000]]  #[s] time list
r2=[1e-6*self for self in range(1,30000,300)]    #[m] hole size range
r1=[-self for self in r2]
r1.reverse()
r=r1+r2
#add change in temp to backround temp
T_list=[T]*len(del2_T)
del_T1=[self+self1 for self,self1 in zip(del2_T,T_list)]#del_T#[self/4 for self in del_T]#[1000]*len(e_list)

print ('T',del_T1[0],'dens',dens[0],'K',numpy.interp(del_T1[0],gas_T,gas_K),'Cp',numpy.interp(del_T1[0],gas_T,gas_Cp2))
alpha=[(numpy.interp(self,gas_T,gas_K))/((numpy.interp(self,gas_T,gas_Cp2))*self1) for self,self1 in zip(del_T1,dens)]

delN=[[[self3*(self4/T)*(self5**2/(self5**2+4*self6*self))*math.exp(-(self2**2)/(self5**2+4*self6*self)) 
    for self2 in r] 
    for self in t] 
    for self3,self4,self5,self6 in zip(num_dens,del_T1,beam_r,alpha)]

totN=[[[is_pos(num_dens[self]-delN[self][self2][self3])     
    for self3 in range(len(r))]
    for self2 in range(len(t))]
    for self in range(len(num_dens))] 

#estimate attenuation change - see notebook 8/7/14
#need mass density change
mass_dens=[[[totN[self][self2][self3]*am/avogadro 
    for self3 in range(len(r))] 
    for self2 in range(len(t))] 
    for self in range(len(num_dens))]    #mass density of gas[kg/m**3]

#this is a lot of values and slows down simulation!
atten_dens=[[[math.exp(((-numpy.interp(self4,gas_keV,gas_MA))/10)*mass_dens[self][self2][self3]*cell_l) 
    for self3 in range(len(r))] 
    for self2 in range(len(t))] 
    for self,self4 in zip(range(len(num_dens)),e_list)]  

#get peak values from attenuation change
atten_chg_peak=[[atten_dens[self][self2][100] for self in range(len(atten_dens))] for self2 in range(len(t))]
e1_index=[0,4,5,20,21,110]
atten_times=[[atten_dens[self][self2][100] for self2 in range(len(t))] for self in e1_index] 
totN_times=[[totN[self][self2][100] for self2 in range(len(t))] for self in e1_index]

###########################################
#"Heat flux between parallel plates through a binary gaseous mixture over the whole range of the Knudsen number." Sharipov et al. Phisica A 378 (2007) 183-193
H=[cell_dia/2-self for self in beam_r]
d_hf=[(self*self1)/(self2*self3) for self,self1,self2,self3 in zip(H,p0_Pa,[amu]*len(e_list),[vo]*len(e_list))]    #eq.2, delta [unitless]
q_hf=[1.875*(1/self)*(1+((2*1.954)/self))**-1 for self in d_hf]                                                 #eq.37  heat flux[unitless]
q_xhf=[self*self1*self2*self3/T for self,self1,self2,self3 in zip(q_hf,p0_Pa,[vo]*len(e_list),del2_T)]             #eq.5   heat flux [W/m**2]
beam_area=[math.pi*2.0*self*cell_l for self in beam_r]                                                                 #surface area of beam [m**2]
q_total=[self*self1 for self,self1 in zip(q_xhf,beam_area)]                                                 #[W] assuming heat transfer from outside of beam surface area
t_hf=[self/self1 for self,self1 in zip([Wi]*len(e_list),q_total)]                                       #[J]/[J/s]=[s] heat flux time

##########################################
#Yiping's steady state solution-notebook~9/25/14

def ss_calc(T_0z_prev,I_0z_prev,r_beam,z_1,z_0,p_nom):
    I_0z=I_0z_prev*math.exp(-(z_1-z_0)*const[0]/T_0z_prev)    
    T_0z=(T_Rz**2.5+(5*T_Rz**.5*const[0]*I_0z*P_i*math.log((cell_dia/2)/r_beam))/(4*math.pi*numpy.interp(T_0z_prev,gas_T,gas_K)))**(1/2.5) 
    n_0z=(p_nom*avogadro)/(Rg*T_0z)
    
    return [I_0z,T_0z,n_0z]

l_chunk=100         #number of length steps
l_steps=[(cell_l/l_chunk)*self for self in range(l_chunk)]      #[m]step length
I_0z_0=1.0                    #initial transmission value, I(0,z)
T_Rz=T                      #[K] T(R,z)
const=[am*self*(numpy.interp(self2,gas_keV,gas_MA))/(10*Rg) for self,self2 in zip(p0_Pa,e_list)]    #[K/m] see notebook 9/22/14

T_0z_test=(T_Rz**2.5+(5*T_Rz**.5*1300.8*I_0z_0*P_i*math.log((.01)/.001))/(4*math.pi*.024))**(1/2.5)        #[K]
T_0z_0=(T_Rz**2.5+(5*T_Rz**.5*const[0]*I_0z_0*P_i*math.log((cell_dia/2)/beam_r[0]))/(4*math.pi*numpy.interp(T_Rz,gas_T,gas_K)))**(1/2.5)        #[K]

#print 'T_0z_test',1300.8,math.log((.01)/.001),.024
#print 'T_0z_0 in',const[0],math.log((cell_dia/2)/beam_r[0]),numpy.interp(T_Rz,gas_T,gas_K)
ss_calc_n=[[] for self in range(l_chunk)]
#print 'ss_calc_n',ss_calc_n
ss_calc_n[0]=ss_calc(T_0z_0,I_0z_0,beam_r[0],l_steps[1],l_steps[0],p0_Pa[0])
print ('ss_calc_n',ss_calc_n[0][1])

#for self in range(l_chunk):    
#    ss_calc_n[self]=(ss_calc(ss_calc_n[self][1],ss_calc_n[self][0],beam_r[0],l_steps[1],l_steps[0],p0_Pa[0]))

###########################################
#estimate ions
#from LCLS-TN-09-5, section 3.6 (in notebook~10/21/13)
#eqn 3.26
ion_dens_pulse=[self/(self1*1.602176565e-19*1000) for self,self1 in zip(Q_vol,e_list)]    #[ions/(pulse*m**3)]

#approximately eq 3.27
ions=[(2*self1/(math.pi*self2**2*self3*1.602176e-19*1000))/(self4/gas_mass) for self1,self2,self3,self4 in zip(q_length,beam_r,e_list,dens)]

#discussion w/ Ryan Coffee 5/7/14
ions2=[N_ph_SCSXR*self*am/(avogadro*math.pi*self2**2) for self,self2 in zip(sigma,beam_r)]      #[photons/(atom*pulse)] 

#estimate mean free path[m], see notebook 5/6/14
gas_dia=2*188e-12                                                                       #[m] van der Waals dia,2* 188e-12 Ar, 155e-12 N2
lambda_kin=[(k*T)/((2**.5)*math.pi*(gas_dia**2)*(self*133.3)) for self in p0]           #conv p0[Torr]*133.3=[Pa]
lambda_atten=[1/(self*self1) for self,self1 in zip(sigma,dens)]                         #[m] 1/e length, average dist a photon travels between collisions w/ atoms of target 

#estimate sonic time to refill beam radius
v_sonic=[(gas_g*Rg*self/am)**.5 for self in [300,500,1000]]
sonic_t_300=[self/v_sonic[0] for self in beam_r]
sonic_t_500=[self/v_sonic[1] for self in beam_r]
sonic_t_1000=[self/v_sonic[2] for self in beam_r]

sonic_f_300=[1/self for self in sonic_t_300]
sonic_f_500=[1/self for self in sonic_t_500]
sonic_f_1000=[1/self for self in sonic_t_1000]


#########################################
#Eval space charge effects: LCLS-TN-07-11, notebook 10/25/13, 1/15/15
N_x=[6.4e15*Wi/self for self in e_list]        #eq 6, [xray quanta/pulse]
N_i_p=[self1*self2*self3 for self1,self2,self3 in zip(dens,sigma,N_x)]      #eq 3 [ions/(m*pulse)]
ionpercent=[100*self1*cell_l*am/(self2*avogadro*cell_vol) for self1,self2 in zip(N_i_p,dens)]           #%[ions/total atoms*pulse]

p_ion=[7.5e-21*self*(self-.4)/(10000*self1*self2*cell_vol) for self,self1,self2 in zip(e_list,sigma,dens)]       #eq 8, this doesn't really make sense...


#######################################
#conventinal calc
#[pzcalc_beam,pzcalc_pump,q_aper,C_a1,S_eff_tot,S_eff_beam,S_pump]

p2zcalc=[pzcalc(self1,self2,ap2type,stage_2[0],stage_2[1],p2L_D) for self1,self2 in zip(p1,ap1)]
[pzb2,pzp2,qz2,Cza1,S_t2,S_b2,S_p2]=[[p2zcalc[self][self1] for self in e_count] for self1 in range(0,len(p2zcalc[0]))]

p3zcalc=[pzcalc(self1,self2,ap3type,stage_3[0],stage_3[1],p3L_D) for self1,self2 in zip(p2,ap2)]
[pzb3,pzp3,qz3,Cza2,S_t3,S_b3,S_p3]=[[p3zcalc[self][self1] for self in e_count] for self1 in range(0,len(p3zcalc[0]))]

p4zcalc=[pzcalc(self1,self2,ap4type,stage_4[0],stage_4[1],p4L_D) for self1,self2 in zip(pzb3,ap3)]
[pzb4,pzp4,qz4,Cza3,S_t4,S_b4,S_p4]=[[p4zcalc[self][self1] for self in e_count] for self1 in range(0,len(p4zcalc[0]))]

p5zcalc=[pzcalc(self1,self2,ap5type,stage_5[0],stage_5[1],p5L_D) for self1,self2 in zip(pzb4,ap4)]
[pzb5,pzp5,qz5,Cza4,S_t5,S_b5,S_p5]=[[p5zcalc[self][self1] for self in e_count] for self1 in range(0,len(p5zcalc[0]))]

p6zcalc=[pzcalc(self1,self2,ap6type,stage_6[0],stage_6[1],p6L_D) for self1,self2 in zip(pzb5,ap5)]
[pzb6,pzp6,qz6,Cza5,S_t6,S_b6,S_p6]=[[p6zcalc[self][self1] for self in e_count] for self1 in range(0,len(p6zcalc[0]))]

#rouging pump calculations
p2r=foreline_max(p2)
p3r=foreline_max(p3)
p4r=foreline_max(p4)
p5r=foreline_max(p5)
p6r=foreline_max(p6)

s2r=[i/j for i,j in zip(q1,p2r)]
s3r=[i/j for i,j in zip(q2,p3r)]
s4r=[i/j for i,j in zip(q3,p4r)]
s5r=[i/j for i,j in zip(q4,p5r)]
s6r=[i/j for i,j in zip(q5,p6r)]

#use pXr and sXr above as guesses and iterate using rcalc to find roughing pressure
#r_x=[#pumps,pump choice,pump cost,pump pipe transmission]
#rcalc(p_c,Q_x,N_c,S_c,trans)
#p2r_1=[rcalc(i,j,r_2[0],r_2[1],r_2[3]) for i,j in zip(p2r,q1)]
#p3r_1=[rcalc(i,j,r_3[0],r_3[1],r_3[3]) for i,j in zip(p3r,q2)]
#p4r_1=[rcalc(i,j,r_4[0],r_4[1],r_4[3]) for i,j in zip(p4r,q3)]
#p5r_1=[rcalc(i,j,r_5[0],r_5[1],r_5[3]) for i,j in zip(p5r,q4)]

p2f=[rcalc(i,j,r_2[0],r_2[1],r_2[3]) for i,j in zip(pzp2,q1)]
[p2fore,pconv2]=[[p2f[self][self1] for self in e_count] for self1 in range(0,len(p2f[0]))]

p3f=[rcalc(i,j,r_3[0],r_3[1],r_3[3]) for i,j in zip(pzp3,q2)]
[p3fore,pconv3]=[[p3f[self][self1] for self in e_count] for self1 in range(0,len(p3f[0]))]

p4f=[rcalc(i,j,r_4[0],r_4[1],r_4[3]) for i,j in zip(pzp4,q3)]
[p4fore,pconv4]=[[p4f[self][self1] for self in e_count] for self1 in range(0,len(p4f[0]))]

p5f=[rcalc(i,j,r_5[0],r_5[1],r_5[3]) for i,j in zip(pzp5,q4)]
[p5fore,pconv5]=[[p5f[self][self1] for self in e_count] for self1 in range(0,len(p5f[0]))]

p6f=[rcalc(i,j,r_6[0],r_6[1],r_6[3]) for i,j in zip(pzp6,q5)]
[p6fore,pconv6]=[[p6f[self][self1] for self in e_count] for self1 in range(0,len(p6f[0]))]

#init columns of data
#[idx,c209,p222,p221,c221,p231,c231,p232,b240,v240,p251,c251,p252,p261,b311,p351,c351,p352,b360,v360,p371]=[[] for self in range(0,21)]
#read data from columns and assign names
#with open('LCLSIEpicsData.csv') as ed:
#    data_in = csv.reader(ed)    
#    for col in data_in:
#        idx.append(col[0])
#        c209.append(col[2])
#        p222.append(col[3])
#        p221.append(col[4])
#        c221.append(col[5])
#        p231.append(col[6])
#        c231.append(col[7])
#        p232.append(col[8])
#        b240.append(col[9])
#        v240.append(col[10])
#        p251.append(col[11])
#        c251.append(col[12])
#        p252.append(col[13])
#        p261.append(col[14])
#        b311.append(col[15])
#        p351.append(col[16])
#        c351.append(col[17])
#        p352.append(col[18])
#        b360.append(col[19])
#        v360.append(col[20])
#        p371.append(col[21])
#             
#ins=[idx,c209,p222,p221,c221,p231,c231,p232,b240,v240,p251,c251,p252,p261,b311,p351,c351,p352,b360,v360,p371]
#ins=[self[350:4000] for self in ins]    #select interesting part of data
#ins=[map(float,self) for self in ins]   #change values from string to float
#[idx,c209,p222,p221,c221,p231,c231,p232,b240,v240,p251,c251,p252,p261,b311,p351,c351,p352,b360,v360,p371]=ins

#z=[closest(self,b311) for self in p0]   #takes a value from p0 and finds closest 
#c1=[[i for i,val in enumerate(b311) if val==z[self]] for self in range(0,len(z))]       #finds index of all entries with each value in z
#z1=[]
#for self in range(0,len(c1)):
#    z1.extend(c1[self])

#z1=[b311.index(self) for self in z]     #this will only pull one index from b311, even if multiple w/ same value (problem?).
#z1=[self+2 for self in z1]
#use z1 to make lists of other useful values.
#ins2=[[self[self2] for self in ins] for self2 in z1]   #select
#ins6=[[[ins[self2][c1[self1][self]] for self in xrange(0,len(c1[self1]))] for self1 in xrange(0,len(c1))] for self2 in xrange(0,len(ins))]  #pulls lists of values matching the indices of matching p0 values (c1)
#ins7=[[sum(ins6[self][self1])/float(len(ins6[self][self1])) for self1 in xrange(0,len(ins6[self]))] for self in xrange(0,len(ins6))]        #average of each list of matching values
#ins8=[[ins7[self][self2] for self in xrange(0,len(ins7))] for self2 in xrange(0,len(c1))]           #reformat back to rows of .csv type list (like ins2)

#[i_dx,c_209,p_222,p_221,c_221,p_231,c_231,p_232,b_240,v_240,p_251,c_251,p_252,p_261,b_311,p_351,c_351,p_352,b_360,v_360,p_371]=[[[] for self in range(0,len(p0))] for self in range(0,21)]
#for self in range(0,len(p0)):
#    [i_dx[self],c_209[self],p_222[self],p_221[self],c_221[self],p_231[self],c_231[self],p_232[self],b_240[self],v_240[self],p_251[self],c_251[self],p_252[self],p_261[self],b_311[self],p_351[self],c_351[self],p_352[self],b_360[self],v_360[self],p_371[self]]=ins8[self]


#init columns of data
[idx,stamp,c209,c221,c231,b240,v240,p251,p261,b311,v311]=[[] for self in xrange(0,11)]
#read data from columns and assign names
with open('LCLSIEpics4-4-13.csv') as ed:
    data_in = csv.reader(ed)    
    for col in data_in:
        idx.append(col[0])
        c209.append(col[2])
        c221.append(col[3])
        c231.append(col[4])
        b240.append(col[5])
        v240.append(col[6])
        p251.append(col[7])
        p261.append(col[8])
        b311.append(col[9])
        v311.append(col[10])
             
ins=[idx,c209,c221,c231,b240,v240,p251,p261,b311,v311]
#filter index for unstabilized values.
inidx=range(174,219)+range(264,345)+range(370,403)+range(480,705)+range(790,908)+range(1054,1469)+range(1547,1685)+range(1712,1799)  
ins=[[self[self2] for self2 in inidx] for self in ins]    #select interesting part of data
#print 'ins', ins
ins=[map(float,self) for self in ins]   #change values from string to float
[idx,c209,c221,c231,b240,v240,p251,p261,b311,v311]=ins

z=[closest(self,b311) for self in p0]   #takes a value from p0 and finds closest 
c1=[[i for i,val in enumerate(b311) if val==z[self]] for self in range(0,len(z))]       #finds index of all entries with each value in z
z1=[]
for self in range(0,len(c1)):
    z1.extend(c1[self])

#z1=[b311.index(self) for self in z]     #this will only pull one index from b311, even if multiple w/ same value (problem?).
#z1=[self+2 for self in z1]
#use z1 to make lists of other useful values.
ins2=[[self[self2] for self in ins] for self2 in z1]   #select
ins6=[[[ins[self2][c1[self1][self]] for self in xrange(0,len(c1[self1]))] 
        for self1 in xrange(0,len(c1))] for self2 in xrange(0,len(ins))]  #pulls lists of values matching the indices of matching p0 values (c1)
ins7=[[sum(ins6[self][self1])/float(len(ins6[self][self1])) 
       for self1 in xrange(0,len(ins6[self]))] for self in xrange(0,len(ins6))]        #average of each list of matching values
ins8=[[ins7[self][self2] for self in xrange(0,len(ins7))] for self2 in xrange(0,len(c1))]           #reformat back to rows of .csv type list (like ins2)

[i_dx,c_209,c_221,c_231,b_240,v_240,p_251,p_261,b_311,v_311]=[[[] for self in range(0,len(p0))] for self in range(0,10)]
for self in range(0,len(p0)):
    [i_dx[self],c_209[self],c_221[self],c_231[self],b_240[self],v_240[self],p_251[self],p_261[self],b_311[self],v_311[self]]=ins8[self]

v_311Tls=[self*.01267/2 for self in v_311]       #convert v_311 from [sccm] to [Torr*l/s], and div by 2 for one side

#end main()

###################################################################################  
name1='%s, %s[m] of %s, Attenuation=%.1e, Temp %.1f[K]' %(name, cell_l, gas, atten, T)

plt.figure(1)
plt.suptitle(name1) 
#fig, ax = plt.subplots(1)
plt.subplot(131)
pt_form= '{} {}'.format(stage_1[0],stage_1[1])
p0_plot=plt.plot(e_list,p0,'r-',label='$p_0$, gas cell',ms=6,mew=0.0)
p_elstat_plot=plt.plot(e_list,p_elstat,'r-.',label='p_elstat, electrostatic crit pressure',ms=6,mew=0.0)
p1_new=plt.plot(e_list,p1,'g-',label='$p_1$,({}) {}'.format(stage_1[0],stage_1[1]),ms=6,mew=0.0)
p2_new=plt.plot(e_list,p2,'b-',label='$p_2$,({}) {}'.format(stage_2[0],stage_2[1]),ms=6,mew=0.0)
p3_new=plt.plot(e_list,p3,'c-',label='$p_3$,({}) {}'.format(stage_3[0],stage_3[1]),ms=6,mew=0.0)
p4_new=plt.plot(e_list,p4,'m-',label='$p_4$,({}) {}'.format(stage_4[0],stage_4[1]),ms=6,mew=0.0)
p5_new=plt.plot(e_list,p5,'k-',label='$p_5$,({}) {}'.format(stage_5[0],stage_5[1]),ms=6,mew=0.0)
p6_new=plt.plot(e_list,p6,c='orange',ls='-',label='$p_6$,({}) {}'.format(stage_6[0],stage_6[1]),ms=6,mew=0.0)

#b_311_plot=plt.plot(e_list,b_311,'r+',label='$b311$ gauge gas cell',ms=6)
#p_261_plot=plt.plot(e_list,p_261,'g+',label='$p261$ gauge',ms=6)
#p_251_plot=plt.plot(e_list,p_251,'b+',label='$p251$ gauge',ms=6)
#p_351_plot=plt.plot(e_list,p_351,'bx',label='$p351$ gauge',ms=6)
#b_240_plot=plt.plot(e_list,b_240,'c+',label='$b240$',ms=6)
#p_231_plot=plt.plot(e_list,p_231,'c+',label='$p231$',ms=6)
#c_231_plot=plt.plot(e_list,c_231,'c+',label='$c231$',ms=6)
#c_371_plot=plt.plot(e_list,p_371,'cx',label='$p371$',ms=6)
#c_221_plot=plt.plot(e_list,c_221,'m+',label='$c221$',ms=6)
#c_209_plot=plt.plot(e_list,c_209,'k+',label='$c209$',ms=6)

#p2z_plot=plt.plot(e_list,pzb2,'b.-',label='$p_2$ std')
#p3z_plot=plt.plot(e_list,pzb3,'c.-',label='$p_3$ std')
#p4z_plot=plt.plot(e_list,pzb4,'m.-',label='$p_4$ std')
#p5z_plot=plt.plot(e_list,pzb5,'k.-',label='$p_5$ std')
#p6z_plot=plt.plot(e_list,pzb6,c='orange',ls='-',marker='.',label='$p_6$ std')

#plt.vlines(1.56,1e-8,100,color='k',linestyles='--')
plt.title('%s System Pressures @ Attenuation=%.1e' %(gas,atten))
plt.xlabel('[keV]')
plt.ylabel('[Torr]')
plt.semilogy()
plt.grid(True)
leg=plt.legend(loc=4, ncol=1,fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)
#plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)

plt.subplot(132)
plt.markersize=10
plt.markeredgewidth=0.0
s1_plot=plt.plot(p1,s1,'go',label='$s_1$',ms=6,mew=0.0)
s2_plot=plt.plot(p2,s2,'bo',label='$s_2$@beam',ms=6,mew=0.0)
s2z_plot=plt.plot(pzb2,S_b2,'b+',label='$s_2$std@beam',ms=6)
s2z_pump_plot=plt.plot(pzp2,S_p2,'b+',label='$s_2$@pump',ms=6)

s3_plot=plt.plot(p3,s3,'co',label='$s_3$@beam',ms=6,mew=0.0)
s3z_plot=plt.plot(pzb3,S_b3,'c+',label='$s_3$std@beam',ms=6)
s3z_pump_plot=plt.plot(pzp3,S_p3,'c+',label='$s_3$@pump',ms=6)

s4_plot=plt.plot(p4,s4,'mo',label='$s_4$',ms=6,mew=0.0)
s4z_plot=plt.plot(pzb4,S_b4,'m+',label='$s_4$std@beam',ms=6)
s4z_pump_plot=plt.plot(pzp4,S_p4,'m+',label='$s_4$@pump',ms=6)

s5_plot=plt.plot(p5,s5,'ko',label='$s_5$',ms=6,mew=0.0)
s5z_plot=plt.plot(pzb5,S_b5,'k+',label='$s_5$std@beam',ms=6)
s5z_pump_plot=plt.plot(pzp5,S_p5,'k+',label='$s_5$@pump',ms=6)

s6_plot=plt.plot(p6,s6,c='orange',ls='none',marker='o',label='$s_6$',ms=6,mew=0.0)
s6z_plot=plt.plot(pzb6,S_b6,c='orange',ls='none',marker='+',label='$s_6$std@beam',ms=6)
s6z_pump_plot=plt.plot(pzp6,S_p6,c='orange',ls='none',marker='+',label='$s_6$@pump',ms=6)

sa1=plt.plot(L4_1,s_mult(stage_1[0],stage_1[1]),'g-',label='({}) {}'.format(stage_1[0],stage_1[1]))
sa2=plt.plot(L4_1,s_mult(stage_2[0],stage_2[1]),'b-',label='({}) {}'.format(stage_2[0],stage_2[1]))
sa3=plt.plot(L4_1,s_mult(stage_3[0],stage_3[1]),'c-',label='({}) {}'.format(stage_3[0],stage_3[1]))
sa4=plt.plot(L4_1,s_mult(stage_4[0],stage_4[1]),'m-',label='({}) {}'.format(stage_4[0],stage_4[1]))
sa5=plt.plot(L4_1,s_mult(stage_5[0],stage_5[1]),'k-',label='({}) {}'.format(stage_5[0],stage_5[1]))
sa6=plt.plot(L4_1,s_mult(stage_6[0],stage_6[1]),c='orange',ls='-',label='({}) {}'.format(stage_6[0],stage_6[1]))

plt.axis([1e-10,1e1,1,1e4])
plt.title('Pumping Speed, S')
plt.xlabel('[Torr]')
plt.ylabel('[l/s]')
plt.loglog()
plt.grid(True)
leg=plt.legend(loc=4, ncol=2,fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)


plt.subplot(133)
q0_plot=plt.plot(p0,q0,'ro',label='ap0, {}, {}[mm]'.format(ap0type,ap0dia),ms=6,mew=0.0)
Mdot_n2_plot=plt.plot(p0,Mdot_n2,'r.-',label='Euler, choked, ap0, round hole in thin plate')
#v_311_plot=plt.plot(p0,v_311Tls,'r+',label='real flow,v311')
q1_plot=plt.plot(p1,q1,'go',label='ap1, {}, {}[mm]'.format(ap1type,ap1dia),ms=6,mew=0.0)
q2_plot=plt.plot(p2,q2,'bo',label='ap2, {}, {}[mm]'.format(ap2type,ap2dia),ms=6,mew=0.0)
q3_plot=plt.plot(p3,q3,'co',label='ap3, {}, {}[mm]'.format(ap3type,ap3dia),ms=6,mew=0.0)
q4_plot=plt.plot(p4,q4,'mo',label='ap4, {}, {}[mm]'.format(ap4type,ap4dia),ms=6,mew=0.0)
q5_plot=plt.plot(p5,q5,'ko',label='ap5, {}, {}[mm]'.format(ap5type,ap5dia),ms=6,mew=0.0)
q6_plot=plt.plot(p6,q6,c='orange',ls='none',marker='o',label='ap6, {}, {}[mm]'.format(ap6type,ap6dia),ms=6,mew=0.0)
plt.title('Beamline Orifice Throughputs, q')
plt.xlabel('[Torr]')
plt.ylabel('[Torr*l/s]')
plt.loglog()
plt.grid(True)
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

##################################################
plt.figure(2)
plt.suptitle(name1)
plt.subplot(131)
de0_plot=plt.plot(p0,de0,'ro',label='$\delta_0$',ms=6,mew=0.0)
de1_plot=plt.plot(p1,de1,'go',label='$\delta_1$',ms=6,mew=0.0)
de2_plot=plt.plot(p2,de2,'bo',label='$\delta_2$',ms=6,mew=0.0)
de3_plot=plt.plot(p3,de3,'co',label='$\delta_3$',ms=6,mew=0.0)
de4_plot=plt.plot(p4,de4,'mo',label='$\delta_4$',ms=6,mew=0.0)
de5_plot=plt.plot(p5,de5,'ko',label='$\delta_5$',ms=6,mew=0.0)
de6_plot=plt.plot(p6,de6,c='orange',ls='none',marker='o',label='$\delta_6$',ms=6,mew=0.0)
plt.hlines(1, 1e-8,10, colors='k', linestyles='solid')
plt.title('Rarefaction Parameter')
plt.xlabel('[Torr]')
plt.ylabel('Rarefaction Parameter, $\delta$')
plt.loglog()
plt.grid(True)
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.subplot(132)
gcost6_plot=plt.plot(e_list,gcost,'b-',label="Attenuation %.1e" %atten)
plt.title('Gas cost for total loss system')
plt.xlabel('[keV]')
plt.ylabel('[$/hr]')
leg=plt.legend(loc=8,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)
plt.grid(True)

plt.subplot(133)
t_q_plot=plt.plot(e_list,t_q,'r-',label=gas)
plt.markersize=6
plt.title('Gas flow for total loss system')
plt.xlabel('[kev]')
plt.ylabel('[moles/hr]')
leg=plt.legend(loc=8,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)
plt.grid(True)


################################################
plt.figure(3)
plt.suptitle(name1)
plt.subplot(121)
V301_c_plot=plt.plot(V301_c_pr,V301_c_sp,'r-',label='V301')
V701_c_plot=plt.plot(V701_c_pr,V701_c_sp,'g-',label='V701')
V1001_c_plot=plt.plot(V1001_c_pr,V1001_c_sp,'b-',label='V1001')
MU100_c1_plot=plt.plot(MU100_c1_pr,MU100_c1_sp,'c-',label='MU100_c1')
MU100_c2_plot=plt.plot(MU100_c2_pr,MU100_c2_sp,'m-',label='MU100_c2')
plt.title('Connection Conductance')
plt.xlabel('[Torr]')
plt.ylabel('[l/s]')
plt.semilogx()
plt.grid(True)
leg=plt.legend(loc=2,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.subplot(122)
D_beamsize_plot=plt.plot(e_list,get_y(D_beamsize,e_list),'k-',label='$D_{beamsize}$')
D_aperture_plot=plt.plot(e_list,get_y(D_aperture,e_list),'k-.',label='$D_{aperture}$')
D_stayclear_plot=plt.plot(e_list,get_y(D_stayclear,e_list),'k--',label='$D_{stayclear}$')
#ap_plot=plt.plot(e_list,ap0,'bo',label='aperture0 dia',ms=6,mew=0.0)
ap_plot=plt.plot(e_list,ap0,'bo',label='Aperture 0 {}:\n{},{},{},{} [mm]'.format(ap0type,ap0dia[0],ap0dia[1],ap0dia[2],ap0dia[3]),ms=6,mew=0.0)
plt.vlines(min_energy,0,25,color='r',linestyle='-', linewidth=.5)
plt.title('SC SXR Beam Diameters at {}[m] from und end'.format(z_adjaper-z_e))
plt.xlabel('[keV]')
plt.ylabel('[mm]')
plt.grid(True)
leg=plt.legend(loc=1,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(4)
plt.suptitle(name1)
plt.subplot(121)
p_diff10_plot=plt.plot(e_list,p_diff10,'ro',label='p1/p0',ms=6,mew=0.0)
p_diff21_plot=plt.plot(e_list,p_diff21,'go',label='p2/p1',ms=6,mew=0.0)
p_diff32_plot=plt.plot(e_list,p_diff32,'bo',label='p3/p2',ms=6,mew=0.0)
p_diff43_plot=plt.plot(e_list,p_diff43,'co',label='p4/p3',ms=6,mew=0.0)
p_diff54_plot=plt.plot(e_list,p_diff54,'mo',label='p5/p4',ms=6,mew=0.0)
CPR_plot=plt.hlines(CPR, 0,2.5, colors='k', linestyles='solid',label='Critical Pressure Ratio, < is choked for hydrodynamic')
plt.title('Pressure Drop Through Each Orifice, Lower # = Better')
plt.xlabel('[keV]')
plt.ylabel('pressure ratio')
plt.grid(True)
leg=plt.legend(loc=8,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)
plt.semilogy()

plt.subplot(122)
#S_p1_plot=plt.plot(q0x,S_p1,'ro',label='pump 1,speed at beam')
S_p1d_plot=plt.plot(e_list,S_p1d,'ro',label='pump 1,speed at pump',ms=6,mew=0.0)
plt.title('Pump Speed Required, Lower # = Better')
plt.xlabel('[Torr]')
plt.ylabel('[l/s]')
plt.grid(True)
leg=plt.legend(loc=8,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)
plt.semilogx()
plt.semilogy()

################################################
plt.figure(5)
plt.suptitle(name1)
plt.subplot(121)
MU100_plot=plt.plot(MU100_pr,MU100_sp,'r-',label='MU100')
MU300_plot=plt.plot(MU300_pr,MU300_sp,'b-',label='MU300')
nXDS15i_plot=plt.plot(nXDS15i_pr,nXDS15i_sp,'c--.',label='nXDS15i')
nXDS20i_plot=plt.plot(nXDS20i_pr,nXDS20i_sp,'k--.',label='nXDS20i')
XDS35i_plot=plt.plot(XDS35i_pr,XDS35i_sp,'g--.',label='XDS35i')
AgilentTS300_plot=plt.plot(AgilentTS300_pr,AgilentTS300_sp,'b--.',label='AgilentTS300')
AgilentTS600_plot=plt.plot(AgilentTS600_pr,AgilentTS600_sp,'r--.',label='AgilentTS600')
PfeifferACP15_plot=plt.plot(PfeifferACP15_pr,PfeifferACP15_sp,'g-o',label='PfeifferACP15')
PfeifferACP28_plot=plt.plot(PfeifferACP28_pr,PfeifferACP28_sp,'b-o',label='PfeifferACP28')
PfeifferACP40_plot=plt.plot(PfeifferACP40_pr,PfeifferACP40_sp,'c-o',label='PfeifferACP40')
AI_ISP250C_plot=plt.plot(AI_ISP250C_pr,AI_ISP250C_sp,'r-',label='AI_ISP250C')
AI_ISP500C_plot=plt.plot(AI_ISP500C_pr,AI_ISP500C_sp,'g-',label='AI_ISP500C')
AI_ISP1000_plot=plt.plot(AI_ISP1000_pr,AI_ISP1000_sp,'b-',label='AI_ISP1000')
NeoDry15E_plot=plt.plot(NeoDry15E_pr,NeoDry15E_sp,'y-',label='NeoDry15E')
NeoDry30E_plot=plt.plot(NeoDry30E_pr,NeoDry30E_sp,'c-',label='NeoDry30E')
EbaraEV_S20_plot=plt.plot(EbaraEV_S20_pr,EbaraEV_S20_sp,'m-',label='EbaraEV_S20')
EbaraEV_S50_plot=plt.plot(EbaraEV_S50_pr,EbaraEV_S50_sp,'k-',label='EbaraEV_S50')
s2r_plot=plt.plot(p2fore,s2r,'b1',label='rough 2')
s3r_plot=plt.plot(p3fore,s3r,'c1',label='rough 3')
s4r_plot=plt.plot(p4fore,s4r,'m1',label='rough 4')
s5r_plot=plt.plot(p5fore,s5r,'k1',label='rough 5')
s6r_plot=plt.plot(p6fore,s6r,c='orange',ls='none',marker='1',label='rough 6')
plt.title('Available Pumps-Speed')
plt.xlabel('[Torr]')
plt.ylabel('[l/s]')
plt.semilogx()
plt.semilogy()
plt.grid(True)
plt.axis([1e-4,1e4,1e-6,100])
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.subplot(122)
#p2r_1_plot=plt.plot(p2x,p2r_1,'b--.',label='p2r_1 available')
#p3r_1_plot=plt.plot(p3x,p3r_1,'c--.',label='p3r_1 available')
#p4r_1_plot=plt.plot(p4x,p4r_1,'m--.',label='p4r_1 available')
#p5r_1_plot=plt.plot(p5x,p5r_1,'k--.',label='p5r_1 available')
p2r_plot=plt.plot(e_list,p2r,'b-',label="p2r max")
p3r_plot=plt.plot(e_list,p3r,'c-',label="p3r max")
p4r_plot=plt.plot(e_list,p4r,'m-',label="p4r max")
p5r_plot=plt.plot(e_list,p5r,'k-',label="p5r max")
p6r_plot=plt.plot(e_list,p6r,c='orange',ls='-',label="p6r max")
p2fore_plot=plt.plot(e_list,p2fore,'b--',label="p2r calc")
p3fore_plot=plt.plot(e_list,p3fore,'c--',label="p3r calc")
p4fore_plot=plt.plot(e_list,p4fore,'m--',label="p4r calc")
p5fore_plot=plt.plot(e_list,p5fore,'k--',label="p5r calc")
p6fore_plot=plt.plot(e_list,p6fore,c='orange',ls='--',label="p6r calc")
#p_252_plot=plt.plot(e_list,p_252,'b+',label='$p252$ gauge',ms=6)
#p_232_plot=plt.plot(e_list,p_232,'c+',label='$p232$ gauge',ms=6)
#p_222_plot=plt.plot(e_list,p_222,'m+',label='$p222$ gauge',ms=6)
plt.title('Foreline Pressures')
plt.xlabel('[keV]')
plt.ylabel('[Torr]')
plt.semilogy()
plt.grid(True)
plt.axis([0,2.5,1e-4,2])
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(6)
plt.suptitle(name1)
plt.subplot(121)
#MAGW300iPKrXe_plot=plt.plot(MAGW300iPKrXe_pr,MAGW300iPKrXe_sp,'g-',label='MAGW300iPKrXe')
#MAGW300P_N2_plot=plt.plot(L4_1,get_y(MAGW300P_N2,L4_1),'r-.',label='MAGW300P_N2')
#MAGW1300Kr_plot=plt.plot(MAGW1300Kr_pr,MAGW1300Kr_sp,'b-',label='MAGW1300Kr')
#MAGW1300Xe_plot=plt.plot(MAGW1300Xe_pr,MAGW1300Xe_sp,'c-',label='MAGW1300Xe')
#MAGW1300N2_plot=plt.plot(MAGW1300N2_pr,MAGW1300N2_sp,'m-',label='MAGW1300N2')
MAGW1300Ar_plot=plt.plot(MAGW1300Ar_pr,MAGW1300Ar_sp,'m-',label='MAGW1300Ar')
#MAGW2000Kr_plot=plt.plot(MAGW2000Kr_pr,MAGW2000Kr_sp,'y-',label='MAGW2000Kr')
#MAGW2000Xe_plot=plt.plot(MAGW2000Xe_pr,MAGW2000Xe_sp,'k-',label='MAGW2000Xe')
#MAGW2000N2_plot=plt.plot(MAGW2000N2_pr,MAGW2000N2_sp,'r-',label='MAGW2000N2')
V1KG_plot=plt.plot(V1KG_pr,V1KG_sp,'r-.',label='Agilent V1KG_Ar')
V1001_plot=plt.plot(V1001_pr,V1001_sp,'b-.',label='Agilent V1001Ar')
V701_plot=plt.plot(V701_pr,V701_sp,'g-.',label='Agilent V701Ar')
TMP1003LMC_Ar_plot=plt.plot(TMP1003LMC_Ar_pr,TMP1003LMC_Ar_sp,'m--.',label='Shimadzu TMP1003LMC_Ar')
TMP1303LMC_Ar_plot=plt.plot(TMP1303LMC_Ar_pr,TMP1303LMC_Ar_sp,'k--.',label='Shimadzu TMP1303LMC_Ar')
#AlcatelMDP5011_N2_plot=plt.plot(AlcatelMDP5011_N2_pr,AlcatelMDP5011_N2_sp,'r:',label='AlcatelMDP5011_N2')
HiPace700_plot=plt.plot(HiPace700_pr,HiPace700_sp,'r-+',label='Pfeiffer HiPace700_Ar')
HiPace1200_plot=plt.plot(HiPace1200_pr,HiPace1200_sp,'c-+',label='Pfeiffer HiPace1200_Ar')
ATH1603M_plot=plt.plot(ATH1603M_pr,ATH1603M_sp,'k-+',label='Pfeiffer ATH1603M_Ar')
TPU1201_Ar_plot=plt.plot(TPU1201_Ar_pr,TPU1201_Ar_sp,'g-+',label='Pfeiffer TPU1201_Ar')
EbaraEMT1300_plot=plt.plot(EbaraEMT1300Ar_pr,EbaraEMT1300Ar_sp,'m-o',label='Ebara EMT1300_Ar')
#VarianTurbo600_plot=plt.plot(L4_1,get_y(VarianTurbo600,L4_1),'g--',label='VarianTurbo600')
plt.title('UHV Pumps-Speed')
plt.xlabel('[Torr]')
plt.ylabel('[l/s]')
plt.semilogx()
plt.semilogy()
plt.grid(True)
plt.axis([1e-5,1e0,1e1,3e3])
leg=plt.legend(loc=3,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.subplot(122)
MU100_plot=plt.plot(MU100_pr,MU100_sp,'r--',label='MU100')
MU300_plot=plt.plot(MU300_pr,MU300_sp,'b--',label='MU300')
nXDS15i_plot=plt.plot(nXDS15i_pr,nXDS15i_sp,'c--.',label='nXDS15i')
nXDS20i_plot=plt.plot(nXDS20i_pr,nXDS20i_sp,'k--.',label='nXDS20i')
XDS35i_plot=plt.plot(XDS35i_pr,XDS35i_sp,'g--.',label='XDS35i')
#TriScroll600Varian_plot=plt.plot(L4_1,get_y(TriScroll600Varian,L4_1),'b--.',label='TriScroll600Varian')
AgilentTS300_plot=plt.plot(AgilentTS300_pr,AgilentTS300_sp,'b--.',label='AgilentTS300')
AgilentTS600_plot=plt.plot(AgilentTS600_pr,AgilentTS600_sp,'r--.',label='AgilentTS600')
PfeifferACP15_plot=plt.plot(PfeifferACP15_pr,PfeifferACP15_sp,'g-o',label='PfeifferACP15')
PfeifferACP28_plot=plt.plot(PfeifferACP28_pr,PfeifferACP28_sp,'b-o',label='PfeifferACP28')
PfeifferACP40_plot=plt.plot(PfeifferACP40_pr,PfeifferACP40_sp,'c-o',label='PfeifferACP40')
AI_ISP250C_plot=plt.plot(AI_ISP250C_pr,AI_ISP250C_sp,'r-',label='AI_ISP250C')
AI_ISP500C_plot=plt.plot(AI_ISP500C_pr,AI_ISP500C_sp,'g-',label='AI_ISP500C')
AI_ISP1000_plot=plt.plot(AI_ISP1000_pr,AI_ISP1000_sp,'b-',label='AI_ISP1000')
NeoDry15E_plot=plt.plot(NeoDry15E_pr,NeoDry15E_sp,'y-',label='NeoDry15E')
NeoDry30E_plot=plt.plot(NeoDry30E_pr,NeoDry30E_sp,'c-',label='NeoDry30E')
EbaraEV_S20_plot=plt.plot(EbaraEV_S20_pr,EbaraEV_S20_sp,'m-',label='EbaraEV_S20')
EbaraEV_S50_plot=plt.plot(EbaraEV_S50_pr,EbaraEV_S50_sp,'k-',label='EbaraEV_S50')
plt.title('Rough Pumps-Speed')
plt.xlabel('[Torr]')
plt.ylabel('[l/s]')
plt.semilogx()
plt.semilogy()
plt.grid(True)
plt.axis([1e-3,1e3,1e-1,1e3])
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(7)
plt.suptitle(name1)
W0_plot=plt.plot(de0,W0,'rx',label='W0')
W1_plot=plt.plot(de1,W1,'gx',label='W1')
W2_plot=plt.plot(de2,W2,'bx',label='W2')
W3_plot=plt.plot(de3,W3,'cx',label='W3')
W4_plot=plt.plot(de4,W4,'mx',label='W4')
W5_plot=plt.plot(de5,W5,'kx',label='W5')
W6_plot=plt.plot(de6,W6,c='orange',ls='none',marker='x',label='W6')
W_pl_plot=plt.plot(L7,L8,'c-',label='plate')
W_tube_plot=plt.plot(L_D5_de,L_D5_W,'y-',label='tube')
plt.title('Reduction in flow rate (W)')
plt.xlabel('Rarefaction Parameter, $\delta$')
plt.ylabel('W')
plt.semilogx()
plt.grid(True)
#plt.axis([1e-8,1e4,1e-1,1e3])
leg=plt.legend(loc=2,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(8)
plt.suptitle(name1)
EVR116_plot=plt.plot(EVR116_V_air,EVR116_q_air,'b-',label='EVR116')
injection_plot=plt.plot(inj_V,q0y2,'gx',label='flow')
plt.title('Injector Range')
plt.xlabel('[Volts]')
plt.ylabel('[Torr*l/s]')
plt.semilogy()
plt.grid(True)
#plt.axis([1e-8,1e4,1e-1,1e3])
leg=plt.legend(loc=2,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(9)
plt.suptitle(name1)
plt.subplot(121)
del_T_plot=plt.plot(e_list,del_T,'r-',label='del_T, similar to eq 3.19 LCLS-TN-09-5')
del2_T_plot=plt.plot(e_list,del2_T,'g-',label='del2_T')
del4_T_plot=plt.plot(e_list,del4_T,'b-',label='del4_T')
del5_T_plot=plt.plot(e_list,del5_T,'m-',label='del5_T, eq. 3.23 LCLS-TN-09-5')
plt.title('Temp rise/pulse')
plt.xlabel('[keV]')
plt.ylabel('[K]')
plt.semilogy()
plt.grid(True)
leg=plt.legend(loc=2,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)


#fig = plt.figure()
ax = plt.subplot(122)
plt.title('Temp [K] increase on axis before next pulse')
lns1=ax.plot(p0,del5_T,'r-',label = 'Temp[K]')
ax2 = ax.twinx()
lns2=ax2.plot(p0,r_therm,'g-',label = 'thermal rad[m]')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
#ax.legend(lns, labs, loc=0)
ax.grid()
ax.set_xlabel("p0[Torr]")
ax.set_ylabel('$\Delta$T[K]',color='r')
for tl in ax.get_yticklabels():
    tl.set_color('r')
ax2.set_ylabel('radius[m]',color='g')
for tl in ax2.get_yticklabels():
    tl.set_color('g')
ax.set_ylim(0,35)
ax2.set_ylim(0,.002)
leg=plt.legend(lns, labs, loc=0,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(10)
plt.suptitle(name1)
plt.subplot(121)
sonic_f_300_plot=plt.plot(e_list,sonic_f_300,'g+')
sonic_f_500_plot=plt.plot(e_list,sonic_f_500,'g.')
sonic_f_1000_plot=plt.plot(e_list,sonic_f_1000,'g.-')
plt.title('Sonic frequency to fill beam rad')
plt.xlabel('[keV]')
plt.ylabel('[Hz]')
plt.semilogy()
plt.grid(True)

#fig = plt.figure()
ax = plt.subplot(122)
plt.title('Sonic time to fill beam rad')
lns1=ax.plot(e_list,beam_r,'r-',label = 'beam radius[m]')
ax2 = ax.twinx()
lns2=ax2.plot(e_list,sonic_t_300,'g+',label = '@300[K]')
lns3=ax2.plot(e_list,sonic_t_500,'g.',label = '@500[K]')
lns4=ax2.plot(e_list,sonic_t_1000,'g.-',label = '@1000[K]')
lns = lns1+lns2+lns3+lns4
labs = [l.get_label() for l in lns]
#ax.legend(lns, labs, loc=0)
ax.grid()
ax.set_xlabel("[keV]")
ax.set_ylabel('beam radius[m]',color='r')
for tl in ax.get_yticklabels():
    tl.set_color('r')
ax2.set_ylabel('time[s]',color='g')
for tl in ax2.get_yticklabels():
    tl.set_color('g')
ax.set_ylim(0,.005)
ax2.set_ylim(0,1e-5)
leg=plt.legend(lns, labs, loc=0,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

plt.figure(11)
plt.suptitle('Attenuation change from gas density hole \n {}'.format(name1)) 
#fig, ax = plt.subplots(1)
plt.subplot(231)
#pt_form= '{} {}'.format(stage_1[0],stage_1[1])
u0=plt.plot(r,atten_dens[e1_index[0]][0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(r,atten_dens[e1_index[0]][1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(r,atten_dens[e1_index[0]][2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(r,atten_dens[e1_index[0]][3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(r,atten_dens[e1_index[0]][4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(r,atten_dens[e1_index[0]][5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(r,atten_dens[e1_index[0]][6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(r,atten_dens[e1_index[0]][7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))

plt.title(r'energy={}[keV],beam rad={:.4f}[m],T={:.0f}[K]'.format(e_list[e1_index[0]],beam_r[e1_index[0]],del_T1[e1_index[0]]))
plt.xlabel('[m]')
plt.ylabel('Transmission')
plt.ylim(atten,1)
plt.grid(True)
plt.semilogy()

#plt.vlines(1.56,1e-8,100,color='k',linestyles='--')
plt.subplot(232)
#pt_form= '{} {}'.format(stage_1[0],stage_1[1])
u0=plt.plot(r,atten_dens[e1_index[1]][0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(r,atten_dens[e1_index[1]][1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(r,atten_dens[e1_index[1]][2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(r,atten_dens[e1_index[1]][3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(r,atten_dens[e1_index[1]][4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(r,atten_dens[e1_index[1]][5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(r,atten_dens[e1_index[1]][6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(r,atten_dens[e1_index[1]][7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))

plt.title(r'energy={}[keV],beam rad={:.4f}[m],T={:.0f}[K]'.format(e_list[e1_index[1]],beam_r[e1_index[1]],del_T1[e1_index[1]]))
plt.xlabel('[m]')
plt.ylabel('Transmission')
plt.ylim(atten,1)
plt.grid(True)
plt.semilogy()

plt.subplot(233)
#pt_form= '{} {}'.format(stage_1[0],stage_1[1])
u0=plt.plot(r,atten_dens[e1_index[2]][0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(r,atten_dens[e1_index[2]][1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(r,atten_dens[e1_index[2]][2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(r,atten_dens[e1_index[2]][3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(r,atten_dens[e1_index[2]][4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(r,atten_dens[e1_index[2]][5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(r,atten_dens[e1_index[2]][6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(r,atten_dens[e1_index[2]][7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))

plt.title(r'energy={}[keV],beam rad={:.4f}[m],T={:.0f}[K]'.format(e_list[e1_index[2]],beam_r[e1_index[2]],del_T1[e1_index[2]]))
plt.xlabel('[m]')
plt.ylabel('Transmission')
plt.ylim(atten,1)
plt.grid(True)
plt.semilogy()

plt.subplot(234)
#pt_form= '{} {}'.format(stage_1[0],stage_1[1])
u0=plt.plot(r,atten_dens[e1_index[3]][0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(r,atten_dens[e1_index[3]][1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(r,atten_dens[e1_index[3]][2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(r,atten_dens[e1_index[3]][3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(r,atten_dens[e1_index[3]][4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(r,atten_dens[e1_index[3]][5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(r,atten_dens[e1_index[3]][6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(r,atten_dens[e1_index[3]][7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))

plt.title(r'energy={}[keV],beam rad={:.4f}[m],T={:.0f}[K]'.format(e_list[e1_index[3]],beam_r[e1_index[3]],del_T1[e1_index[3]]))
plt.xlabel('[m]')
plt.ylabel('Transmission')
plt.ylim(atten,1)
plt.grid(True)
plt.semilogy()

plt.subplot(235)
#pt_form= '{} {}'.format(stage_1[0],stage_1[1])
u0=plt.plot(r,atten_dens[e1_index[4]][0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(r,atten_dens[e1_index[4]][1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(r,atten_dens[e1_index[4]][2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(r,atten_dens[e1_index[4]][3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(r,atten_dens[e1_index[4]][4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(r,atten_dens[e1_index[4]][5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(r,atten_dens[e1_index[4]][6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(r,atten_dens[e1_index[4]][7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))

plt.title(r'energy={}[keV],beam rad={:.4f}[m],T={:.0f}[K]'.format(e_list[e1_index[4]],beam_r[e1_index[4]],del_T1[e1_index[4]]))
plt.xlabel('[m]')
plt.ylabel('Transmission')
plt.ylim(atten,1)
plt.grid(True)
plt.semilogy()

plt.subplot(236)
#pt_form= '{} {}'.format(stage_1[0],stage_1[1])
u0=plt.plot(r,atten_dens[e1_index[5]][0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(r,atten_dens[e1_index[5]][1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(r,atten_dens[e1_index[5]][2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(r,atten_dens[e1_index[5]][3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(r,atten_dens[e1_index[5]][4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(r,atten_dens[e1_index[5]][5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(r,atten_dens[e1_index[5]][6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(r,atten_dens[e1_index[5]][7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))

plt.title(r'energy={}[keV],beam rad={:.4f}[m],T={:.0f}[K]'.format(e_list[e1_index[5]],beam_r[e1_index[5]],del_T1[e1_index[5]]))
plt.xlabel('[m]')
plt.ylabel('Transmission')
plt.ylim(atten,1)
plt.grid(True)
plt.semilogy()
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)


plt.figure(12)
plt.suptitle('Attenuation change from gas density hole \n {}'.format(name1))
u0=plt.plot(e_list,atten_chg_peak[0],'r-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[0],1/t[0]))
u1=plt.plot(e_list,atten_chg_peak[1],'g-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[1],1/t[1]))
u2=plt.plot(e_list,atten_chg_peak[2],'b-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[2],1/t[2]))
u3=plt.plot(e_list,atten_chg_peak[3],'y-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[3],1/t[3]))
u4=plt.plot(e_list,atten_chg_peak[4],'m-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[4],1/t[4]))
u5=plt.plot(e_list,atten_chg_peak[5],'k-',label='{:.2e}[s],{:.0f}[Hz]'.format(t[5],1/t[5]))
u6=plt.plot(e_list,atten_chg_peak[6],'r-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[6],1/t[6]))
u7=plt.plot(e_list,atten_chg_peak[7],'g-.',label='{:.2e}[s],{:.0f}[Hz]'.format(t[7],1/t[7]))
plt.semilogy()
plt.xlabel('[keV]')
plt.ylabel('Transmission')
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

fig, ax = plt.subplots()
plt.suptitle('Attenuation change from gas density hole \n {}'.format(name1))
plt.plot(t,atten_times[0],'r-',label='{}[keV]'.format(e_list[e1_index[0]]))
plt.plot(t,atten_times[1],'g-',label='{}[keV]'.format(e_list[e1_index[1]]))
plt.plot(t,atten_times[2],'b-',label='{}[keV]'.format(e_list[e1_index[2]]))
plt.plot(t,atten_times[3],'y-',label='{}[keV]'.format(e_list[e1_index[3]]))
plt.plot(t,atten_times[4],'m-',label='{}[keV]'.format(e_list[e1_index[4]]))
plt.plot(t,atten_times[5],'k-',label='{}[keV]'.format(e_list[e1_index[5]]))
plt.xlabel('time')
plt.ylabel('Transmission')
plt.semilogy()
plt.grid(True)
majorFormatter = FormatStrFormatter('%.1e')
ax.xaxis.set_major_formatter(majorFormatter)
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

fig, ax = plt.subplots()
plt.suptitle('Density change from gas density hole \n {}'.format(name1))
plt.plot(t,totN_times[0],'r-',label='{}[keV]'.format(e_list[e1_index[0]]))
plt.plot(t,totN_times[1],'g-',label='{}[keV]'.format(e_list[e1_index[1]]))
plt.plot(t,totN_times[2],'b-',label='{}[keV]'.format(e_list[e1_index[2]]))
plt.plot(t,totN_times[3],'y-',label='{}[keV]'.format(e_list[e1_index[3]]))
plt.plot(t,totN_times[4],'m-',label='{}[keV]'.format(e_list[e1_index[4]]))
plt.plot(t,totN_times[5],'k-',label='{}[keV]'.format(e_list[e1_index[5]]))
plt.axhline(num_dens[e1_index[0]], color='k', ls='-.')
plt.axhline(num_dens[e1_index[1]], color='k', ls='-.')
plt.axhline(num_dens[e1_index[2]], color='k', ls='-.')
plt.axhline(num_dens[e1_index[3]], color='k', ls='-.')
plt.axhline(num_dens[e1_index[4]], color='k', ls='-.')
plt.axhline(num_dens[e1_index[5]], color='k', ls='-.')

plt.xlabel('time[s]')
plt.ylabel('total number density[particles/m^3]')
plt.semilogy()
majorFormatter = FormatStrFormatter('%.1e')
ax.xaxis.set_major_formatter(majorFormatter)
plt.grid(True)
leg=plt.legend(loc=4,ncol=1, fancybox=True, shadow=True, numpoints=1)
leg.get_frame().set_alpha(0.5)

print ("plots ready!")
plt.show()
#
#

