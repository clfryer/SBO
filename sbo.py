import numpy as np
import matplotlib.pyplot as plt
import math
import random
from matplotlib import pyplot;
from pylab import genfromtxt;

#Tells the code what opacity file to read in
opa_file_name = "kappa_file_C_O.txt"

time_counter = 500
dtime = 2.0
wave_counter = 1000
gam_counter = 10000
rho = 1e-13
vfac = 0.1
dde=1.015
emergt=300.0

jet_emit = np.zeros((gam_counter))
jet_gam = np.zeros((gam_counter))
jet_temp = np.zeros((gam_counter))
jet_area = np.zeros((gam_counter))
jet_ene = np.zeros((gam_counter))
igam = np.zeros((gam_counter))

p1=4.0
jetmax=2.0
efac=10.0
djet=jetmax/10000.

jetenetot=0.0
for x in range(gam_counter):
   jet_gam[x]=1.0*(djet*x+1.0+djet)
   facnr=1.0-1.0/(jet_gam[x]*jet_gam[x])
   jet_temp[x]=vfac*facnr*1900.*pow(jet_gam[x]/100.,0.5)*pow(rho/1e-10,-0.25)
#   if (x < 100): 
#      print(jet_temp[x],facnr,jet_gam[x])
   if jet_gam[x] > 0. and jet_gam[x] < jetmax:
#      jet_area[x]=1.e20*200.*pow(jet_gam[x],-1.5)
      jet_area[x]=1.e23
#*1./jetmax/efac
      if (jet_gam[x]*facnr < 10.0):
          jet_ene[x]=3e42/1.15*pow(jet_gam[x]*facnr,-p1)
#         jet_ene[x]=jet_area[x]*1.e24/efac*137.0*pow(jet_temp[x],4.0)
          ijet=x
      else:
         jet_ene[x]=jet_ene[ijet]/100.*jet_gam[ijet]/jet_gam[x]
#         jet_ene[x]=1e43/jet_gam[x]/efac
#         jet_area[x]=0.
#         jet_temp[x]=0.00001
   else:
      jet_area[x]=0.
      jet_temp[x]=0.00001
      jet_ene[x]=0.0
#   print(jet_ene[x])
   if (facnr*jet_gam[x] > 0.1):
      jetenetot=jetenetot+jet_ene[x]
   igam[x]=int(math.log(jet_gam[x])/math.log(dde))
#   print(facnr,igam[x],jet_gam[x],jet_temp[x])
#   print(x,facnr,jet_gam[x],igam[x])
#   print(jet_gam[x],igam[x],jet_temp[x]*jet_gam[x])

print(jetenetot)

#Gets the opacity data
def kappa_file(file_name):
    # Open a file
    kapfile = open(file_name, "r")

    den_counter = 10
    temp_counter = 69
    counter = 0
    tempstep = 0
    denstep = 0

    useless = 0

    temp_opa = np.zeros((temp_counter))
    den_opa = np.zeros((den_counter))
    kappa_opa = np.zeros((temp_counter,den_counter))

    kappa_file = kapfile.readlines()

    slice_temp = slice(65,75,1)
    slice_den = slice(2,12,1)
    slice_kappa = slice(14,24,1)

    for line in kappa_file:
        #Is the counter for what time step it is on
        tempstep = (counter // 12)
        #Gets rid of the last section of useless text
        if(counter >= ((den_counter + 2) * temp_counter)):
            useless = 0
        #Skips the first useless line of each unit
        elif(counter == ((2 + den_counter) * (tempstep))):
            useless = 0
        #Gets the temperature from the file
        elif(counter == (1 + (2 + den_counter) * (tempstep))):
            temp_opa[tempstep] = kappa_file[counter][slice_temp]
            denstep = -1
        #Gets the density and kappa from the file
        else:
            denstep = denstep + 1
            den_opa[denstep] = kappa_file[counter][slice_den]
            kappa_opa[tempstep, denstep] = kappa_file[counter][slice_kappa]
        #Updates the counter
        counter = counter + 1

    #Close opened file
    kapfile.close()

    return temp_opa, den_opa, kappa_opa

#Interpolates the different kappas
def interpolate(den_inter, temp_inter, kappa, den, temp):
    useless = 0
    #Finds the upper and lower known value limit
    den_n = 0
    den_n_plus = 0
    den_outbound = 1
    for dencount in range(len(den)):
        if(den[dencount] <= den_inter):
            den_n = dencount
            den_n_plus = dencount + 1
            den_outbound = 0
        else:
            useless = 0
        if(den_n_plus >= len(den)):
            den_outbound = 1
            den_n = dencount
            den_n_plus = dencount
        else:
            useless = 0
    temp_n = 0
    temp_n_plus = 0
    temp_outbound = 1
    for tempcount in range(len(temp)):
        if(temp[tempcount] <= temp_inter):
            temp_n = tempcount
            temp_n_plus = tempcount + 1
            temp_outbound = 0
        else:
            useless = 0
        if(temp_n_plus >= len(temp)):
            temp_outbound = 1
            temp_n = tempcount
            temp_n_plus = tempcount
        else:
            useless = 0
    
    if(den_outbound == 0 and temp_outbound == 0):
        #Calculates the percentage of each value the interpolate value is
        x = ((den_inter - den[den_n]) / (den[den_n_plus] - den[den_n]))
        y = ((temp_inter - temp[temp_n]) / (temp[temp_n_plus] - temp[temp_n]))
    elif(den_outbound == 1 and temp_outbound == 0):
        x = 0
        y = ((temp_inter - temp[temp_n]) / (temp[temp_n_plus] - temp[temp_n]))
    elif(den_outbound == 0 and temp_outbound == 1):
        x = ((den_inter - den[den_n]) / (den[den_n_plus] - den[den_n]))
        y = 0
    else:
        x = 0
        y = 0

    #Calculates the interpolated value of the Rossland Opacity
    kappa_inter = ((x * y * kappa[temp_n_plus][den_n_plus]) + ((1 - x) * y * kappa[temp_n_plus][den_n]) + (x * (1-y) * kappa[temp_n][den_n_plus]) + ((1-x) * (1-y) * kappa[temp_n][den_n]))
    return kappa_inter
    
    #Calculates the percentage of each value the interpolate value is
    x = ((den_inter - den[den_n]) / (den_n_plus - den_n))
    y = ((temp_inter - temp[temp_n]) / (temp_n_plus - temp_n))

    #Calculates the interpolated value of the Rossland Opacity
    kappa_inter = ((x * y * kappa[temp_n_plus][den_n_plus]) + ((1 - x) * y * kappa[temp_n_plus][den_n]) + (x * (1-y) * kappa[temp_n][den_n_plus]) + ((1-x) * (1-y) * kappa[temp_n][den_n]))

    return kappa_inter

#Reads in the opacity file
temp_opa, den_opa, kappa_opa = kappa_file(opa_file_name)

def emitk_function(kappak,enek,dek,tempk):
   if (enek/tempk<50.):
      expval=math.exp(enek/tempk)-1.0
      rat=dek*kappak*pow(enek,3)/expval
   else:
      rat=0.e0
   return rat


file=open("spec.out","w")
filel=open("lum.out","w")
files=open("specene.out","w")

rattot=0.   
enek = np.zeros((wave_counter))
dek = np.zeros((wave_counter))

for x in range(wave_counter):
   enek[x]=1.0*pow(dde,x)
   if (x == 0):
      dek[x]=enek[x]
   else:
      dek[x]=enek[x]-enek[x-1]


time = 0.0
for z in range(time_counter):
   emit = np.zeros((wave_counter))
   time = time + dtime
   print(time) 

   for y in range(len(jet_temp)):
#  Calls the Kappa Interpolate function
      kappak = interpolate(rho, (jet_temp[y]/1000.), kappa_opa, den_opa, temp_opa)
      emittot=0.0
      for x in range(wave_counter):
         emittot=emittot+emitk_function(kappak,enek[x],dek[x],jet_temp[y])
#         if (enek[x] > 2000 and jet_temp[y] > 1000 and enek[x] < jet_temp[y]*40.0):
#            print(emittot,enek[x],jet_temp[y],emitk_function(kappak,enek[x],dek[x],jet_temp[y]))

      emittot=1.0283e12*pow(jet_temp[y],4)/(emittot+1e-20)
      jet_emit[y]=1.0283e12*pow(jet_temp[y],4)*jet_area[y]*(math.exp(-emergt/time))

      for x in range(wave_counter):
         z=int(x+igam[y])
         z = min(z,wave_counter-1)
         emit[z]=emit[z]+emittot*jet_area[y]*emitk_function(kappak,enek[x],dek[x],jet_temp[y])*(math.exp(-emergt/time))
#      print(z,emit[z],emitk_function(kappak,enek[x],dek[x],jet_temp[y]))

      if (jet_ene[y] > 1e5):
         if (jet_ene[y] > dtime*jet_emit[y]):
            jfac = (jet_ene[y] - dtime*jet_emit[y])/jet_ene[y]
            jet_ene[y] = jet_ene[y] - dtime*jet_emit[y]
            jet_temp[y] = pow(jfac,0.25)*jet_temp[y]
         else:
            jfac = 0.0
            jet_ene[y] = 0.0
            jet_temp[y] = 0.00001
            jet_area[y] = 0.0
      else:
            jfac = 0.0
            jet_ene[y] = 0.0
            jet_temp[y] = 0.00001
            jet_area[y] = 0.0
#      print(y,jet_ene[y],dtime*jet_emit[y],jfac)
      files.write('%s'%time + ' ' + '%s'%jet_ene[y] + ' ' + '%s'%jet_emit[y] + '\n')

   lumgp1kev=0
   lumgp3kev=0
   lumgp5kev=0
   lumg1p5kev=0
   lumg100kev=0
   for x in range(wave_counter):
      if (enek[x] > 100):
         lumgp1kev=lumgp1kev+emit[x]
      if (enek[x] > 300):
         lumgp3kev=lumgp3kev+emit[x]
      if (enek[x] > 500):
         lumgp5kev=lumgp5kev+emit[x]
      if (enek[x] > 1500):
         lumg1p5kev=lumg1p5kev+emit[x]
      if (enek[x] > 100000):
         lumg100kev=lumg100kev+emit[x]

   file.write('%s'%time + ' ' + '%s'%lumgp1kev + '\n')
   filel.write('%s'%time + ' ' + '%s'%lumgp1kev + ' ' + '%s'%lumgp3kev + ' ' + '%s'%lumgp5kev + ' '  + '%s'%lumg1p5kev + ' '  + '%s'%lumg100kev + '\n')

   for x in range(wave_counter):
      enekp=enek[x]/1000.
      emit[x]=emit[x]/dek[x]
      file.write('%s'%enekp + ' ' + '%s'%emit[x] + '\n')

