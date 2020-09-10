# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 21:59:18 2012

@author: schaferd
"""
name='GasAttenIISXR_3_3_0_in'
#this is the input file for 3_3




#####################################################################
#inputs

#input for pumping stages, 
#stage_x=[#pumps,pump choice,pump cost]
#pump choice:("MAGW300P","MAGW300iP","MAGW1300","MAGW2000","MU100","VarianTurbo600","TMP1003LMC_Ar","HiPace1200") Check pump_data.py
#r_x=[#pumps,pump choice,pump cost,pump pipe transmission]
#r_x pump choice:"nXDS15i","nXDS20i","XDS35i","AgilentTS300","AgilentTS600","PfeifferACP15","PfeifferACP28","PfeifferACP40","AI_ISP250C","AI_ISP500C","AI_ISP1000","NeoDry15E","NeoDry30E","DS600"
#apXtype="plate","plate_sq", or "tube". plate=plate with round aperture, plate_sq= plate with square aperture, tube=round tube with L/D=5
#apXdia=[list of available aperture diameters, must be increasing]
#assuming "tube" applies to pipes with L/D=5. 

#choose gas: ("Kr","Xe","N2","Ar","Ne","O2")
gas="N2"

ap0type="plate_sq"
ap0dia=[5.5,8.0,10.0,13.0]
#ap0dia=[5.5]
p0_error=.12            #gauge error, % of reading (MKS)

stage_1=[1,"EbaraEV_S50",21000]
p1L_D=[.001*.0602,.0602]#[.3,.0602]   #pump connection [length, dia] [m]
p1p_conn=[60.198,1060.91]           #pump1 tube [dia,length][mm]

ap1type="plate"
ap1dia=[15.0]

stage_2=[1,"EbaraEV_S50",21000]
p2L_D=[.001*.0602,.0602]                       #pump connection [length, dia] [m]

ap2type="plate_sq"
ap2dia=[5.5,8.0,10.0,13.0]
#ap2dia=[13.0]

stage_3=[1,"EbaraEMT1300",21000]
p3L_D=[.7*.18,.197]                       #pump connection [length, dia] [m]

ap3type="tube"
ap3dia=[15.0]

stage_4=[1,"EbaraEMT1300",21000]
p4L_D=[.7*.18,.197]                       #pump connection [length, dia] [m]

ap4type="tube"
ap4dia=[15.0]

stage_5=[1,"EbaraEMT1300",21000]
p5L_D=[.7*.18,.197]                       #pump connection [length, dia] [m]

ap5type="tube"
ap5dia=[15.0]

stage_6=[1,"EbaraEMT1300",21000]
p6L_D=[.18,.197]                       #pump connection [length, dia] [m]

ap6type="tube"
ap6dia=[15.0]

#roughing pumps
r_2=[1,"AI_ISP500C",8000,.7]
r_3=[1,"AI_ISP500C",8000,.7]
r_4=[1,"AI_ISP250C",8000,.7]
r_5=[1,"AI_ISP250C",8000,.7]
r_6=[1,"AI_ISP250C",8000,.7]

#which turbos to which foreline (roughing pump)? [st1,st2,...st4]
#fore1=[st2]
#fore2=[st3,st4,st5]

q_outgas=2.0e-8 #3.0e-6       #[Torr*l/s], estimated chamber outgassing/# of pumps
keV=.2      #iteration start beam energy[keV], tables have .16-2.5[keV]
keV_end=1.3
min_energy=.4   #[keV] below this energy allowed to clip beam (Tor email 2/1/13), updated 12/14 to .4 (Yiping) 
atten=1e-10       #start attenuation for loop, requested attenuation 1<atten<1e-6, actually specifying transmission
cell_l=14.7    #cell length [m]
cell_dia=.0413  #cell dia [m]=1.625[in]
T=300.0     # temperature in [K]

##########################################
#beam source parameters-see XTES System PRD and Yiping email 12/3/14 and saved in gas atten/ref on sharepoint
linac='SXR_SC'      #choose linac, choices are 'SXR_SC', 'HXR_SC', and 'HXR_Cu'

#Beam Divergence (FWHM)[microrad] - From: LCLSII-TN-15-20
def B_FWHM(keV):
    if linac=='SXR_SC':    
        B_FWHM= 1.28531e-5*keV**(-8.37861e-1)
    if linac=='HXR_SC':
        B_FWHM= 1.37652e-5*keV**(-8.88487e-1)
    if linac=='HXR_Cu':
        B_FWHM= 1.18385e-5*keV**(-8.29988e-1)
    return B_FWHM

#Effective Source Location [m] - From: LCLSII-TN-15-20
def z_0(keV):
    if linac=='SXR_SC':
        z_0= -6.8658*keV**4 + 28.8017*keV**3 - 46.1806*keV**2 + 46.1743*keV - 51.3787
    if linac=='HXR_SC':
        z_0= 1.24313e-1*keV**4 - 5.61271e-2*keV**3 - 1.52348*keV**2 + 1.14013e1*keV - 8.21218e1
    if linac=='HXR_Cu':
        z_0= -3.695930e-6*keV**6 + 2.787051e-4*keV**5 - 7.848709e-3*keV**4 + 1.088181e-1*keV**3 - 1.057413*keV**2 + 9.881049*keV - 7.310326e1
    return z_0

if linac=='SXR_SC':    
    z_e=655.516492       #[m] end of SXR undulator
    B_steer=10.0    #[microrad]
    D_offset=1.0      #[mm]
    F_opr=2.0
    F_diff=2.0
    n_exp=-.838
    E_0=.4         #[keV]
    D_min=5.0       #[mm]
    z_loc=715.625768  #728.22    #[m]
    dev_len=19.286                          #device length [m]    
    z_adjaper=z_loc+7.35           #adjustable aperture location wrt device center, 14.7[m]/2=7.35[m]
    z_start=z_loc-dev_len/2.
    z_end=z_loc+dev_len/2.
    
    
if linac=='HXR_SC':    
    z_e=650.9399844       #[m] end of HXR undulator
    B_steer=10.0    #[microrad]
    D_offset=1.0      #[mm]
    F_opr=2.0
    F_diff=2.0
    n_exp=-.888
    E_0=1.8         #[keV]
    D_min=5.0       #[mm]
    z_loc=729.295995    #[m]
    dev_len=10.183                          #device length [m]
    z_adjaper=z_loc+2.950845            #adjustable aperture location wrt device center
    z_start=z_loc-dev_len/2.
    z_end=z_loc+dev_len/2.

    
if linac=='HXR_Cu':    
    z_e=650.9399844       #[m] end of HXR undulator
    B_steer=10.0    #[microrad]
    D_offset=1.0      #[mm]
    F_opr=2.0
    F_diff=2.0
    n_exp=-.83
    E_0=1.5         #[keV]
    D_min=5.0       #[mm]
    z_loc=729.295995    #[m]
    dev_len=10.183                          #device length [m]
    z_adjaper=z_loc+2.950845            #adjustable aperture location wrt device center
    z_start=z_loc-dev_len/2.
    z_end=z_loc+dev_len/2.


#estimate gas heating 
#from LCLS-TN-09-5, section 3.3 (in notebook~10/7/13)
eta=1.0                                  #absorbed energy %
Wi=2e-3                                 #initial pulse energy[J]
P_i=200                                 #[W], input power from beam, max
rep_rate=1000000.0                      #[Hz]
pulse_period=1/rep_rate               #[sec]
N_ph_SCSXR=46.0e12                      #[photons/pulse]
N_ph_SCHXR=9.1e12                       #[photons/pulse]
N_ph_CuHXR=18.0e12                      #[photons/pulse]

#source=664.0                        #FEL source [m] PRD z_0
#und_end=655.546                       #undulator end [m] PRD z_e
#z_loc=727.25                     #device z location in [m]


Rg=8.3144621                                     # gas constant[J/(mole*K)]
T0=273.15                                      # reference temperature
k=1.3806488e-23                                 #Boltzmann [J/K]
C=.83                               #Discharge coeff, "Numerical sim of rare gas flow thr thin orifice",Table 2, Sharipov 
#index to next energy loop, remark print statements when increment gets smaller, .05-.01 work
e_index=.01
        

# end of input data
######################################################################



