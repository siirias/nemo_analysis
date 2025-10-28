#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 09:42:49 2020

@author: oikkonea
"""



from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
#
##### This script calculates the length of ice season and the annual maximum ice thickness in coastal stations
# Input data: time series of sea ice concentration and volume, extracted from model output with Ice_data_stations.py



########## KEMI ######################
A001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_A001.txt', delimiter = ",")
B001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_B001.txt', delimiter = ",")
D001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_D001.txt', delimiter = ",")

A002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_A002.txt', delimiter = ",")
B002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_B002.txt', delimiter = ",")
D002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_D002.txt', delimiter = ",")

A005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_A005.txt', delimiter = ",")
B005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_B005.txt', delimiter = ",")
D005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_D005.txt', delimiter = ",")

HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kemi_HC.txt', delimiter = ",")

Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Kemi_Ajos_Number_of_days_with_ice.csv', delimiter = ",")
ObsH=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Kemi_maxH.csv', delimiter = ",")


year_A001=A001[:,0]
day_A001=A001[:,1]
C_A001=A001[:,2]
Vol_A001=A001[:,3]
year_B001=B001[:,0]
day_B001=B001[:,1]
C_B001=B001[:,2]
Vol_B001=B001[:,3]
year_D001=D001[:,0]
day_D001=D001[:,1]
C_D001=D001[:,2]
Vol_D001=D001[:,3]

year_A002=A002[:,0]
day_A002=A002[:,1]
C_A002=A002[:,2]
Vol_A002=A002[:,3]
year_B002=B002[:,0]
day_B002=B002[:,1]
C_B002=B002[:,2]
Vol_B002=B002[:,3]
year_D002=D002[:,0]
day_D002=D002[:,1]
C_D002=D002[:,2]
Vol_D002=D002[:,3]

year_A005=A005[:,0]
day_A005=A005[:,1]
C_A005=A005[:,2]
Vol_A005=A005[:,3]
year_B005=B005[:,0]
day_B005=B005[:,1]
C_B005=B005[:,2]
Vol_B005=B005[:,3]
year_D005=D005[:,0]
day_D005=D005[:,1]
C_D005=D005[:,2]
Vol_D005=D005[:,3]

year_HC=HC[:,0]
day_HC=HC[:,1]
C_HC=HC[:,2]
Vol_HC=HC[:,3]

Year_A002=np.hstack((year_A001,year_A002))
Day_A002=np.hstack((day_A001,day_A002))
Conc_A002=np.hstack((C_A001,C_A002))
Volume_A002=np.hstack((Vol_A001,Vol_A002))
Thickness_A002=Volume_A002/Conc_A002

Year_B002=np.hstack((year_B001,year_B002))
Day_B002=np.hstack((day_B001,day_B002))
Conc_B002=np.hstack((C_B001,C_B002))
Volume_B002=np.hstack((Vol_B001,Vol_B002))
Thickness_B002=Volume_B002/Conc_B002

Year_D002=np.hstack((year_D001,year_D002))
Day_D002=np.hstack((day_D001,day_D002))
Conc_D002=np.hstack((C_D001,C_D002))
Volume_D002=np.hstack((Vol_D001,Vol_D002))
Thickness_D002=Volume_D002/Conc_D002

Year_A005=np.hstack((year_A001,year_A005))
Day_A005=np.hstack((day_A001,day_A005))
Conc_A005=np.hstack((C_A001,C_A005))
Volume_A005=np.hstack((Vol_A001,Vol_A005))
Thickness_A005=Volume_A005/Conc_A005

Year_B005=np.hstack((year_B001,year_B005))
Day_B005=np.hstack((day_B001,day_B005))
Conc_B005=np.hstack((C_B001,C_B005))
Volume_B005=np.hstack((Vol_B001,Vol_B005))
Thickness_B005=Volume_B005/Conc_B005

Year_D005=np.hstack((year_D001,year_D005))
Day_D005=np.hstack((day_D001,day_D005))
Conc_D005=np.hstack((C_D001,C_D005))
Volume_D005=np.hstack((Vol_D001,Vol_D005))
Thickness_D005=Volume_D005/Conc_D005

Thickness_HC=Vol_HC/C_HC

Annual_max_A002=[]
Annual_max_B002=[]
Annual_max_D002=[]
Ice_days_A002=[]
Ice_days_B002=[]
Ice_days_D002=[]

Annual_max_A005=[]
Annual_max_B005=[]
Annual_max_D005=[]
Ice_days_A005=[]
Ice_days_B005=[]
Ice_days_D005=[]

Annual_max_HC=[]
Ice_days_HC=[]


for i in range(1975,2059):

    Hi_A002_a=Thickness_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.8)]
    Hi_A002_s=Thickness_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.8)]
    Hi_A002=np.hstack((Hi_A002_a,Hi_A002_s))
    Annual_max_A002.append(np.max(Hi_A002))
    Ci_A002_a=Conc_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.2)]
    Ci_A002_s=Conc_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.2)]
    Ci_A002=np.hstack((Ci_A002_a,Ci_A002_s))
    Ice_days_A002.append(len(Ci_A002))
  
    Hi_B002_a=Thickness_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.8)]
    Hi_B002_s=Thickness_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.8)]
    Hi_B002=np.hstack((Hi_B002_a,Hi_B002_s))
    Annual_max_B002.append(np.max(Hi_B002))
    Ci_B002_a=Conc_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.2)]
    Ci_B002_s=Conc_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.2)]
    Ci_B002=np.hstack((Ci_B002_a,Ci_B002_s))
    Ice_days_B002.append(len(Ci_B002))

    Hi_D002_a=Thickness_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.8)]
    Hi_D002_s=Thickness_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.8)]
    Hi_D002=np.hstack((Hi_D002_a,Hi_D002_s))
    Annual_max_D002.append(np.max(Hi_D002))
    Ci_D002_a=Conc_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.2)]
    Ci_D002_s=Conc_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.2)]
    Ci_D002=np.hstack((Ci_D002_a,Ci_D002_s))
    Ice_days_D002.append(len(Ci_D002))
    
    Hi_A005_a=Thickness_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.8)]
    Hi_A005_s=Thickness_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.8)]
    Hi_A005=np.hstack((Hi_A005_a,Hi_A005_s))
    Annual_max_A005.append(np.max(Hi_A005))
    Ci_A005_a=Conc_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.2)]
    Ci_A005_s=Conc_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.2)]
    Ci_A005=np.hstack((Ci_A005_a,Ci_A005_s))
    Ice_days_A005.append(len(Ci_A005))
  
    Hi_B005_a=Thickness_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.8)]
    Hi_B005_s=Thickness_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.8)]
    Hi_B005=np.hstack((Hi_B005_a,Hi_B005_s))
    Annual_max_B005.append(np.max(Hi_B005))
    Ci_B005_a=Conc_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.2)]
    Ci_B005_s=Conc_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.2)]
    Ci_B005=np.hstack((Ci_B005_a,Ci_B005_s))
    Ice_days_B005.append(len(Ci_B005))

    Hi_D005_a=Thickness_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.8)]
    Hi_D005_s=Thickness_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.8)]
    Hi_D005=np.hstack((Hi_D005_a,Hi_D005_s))
    Annual_max_D005.append(np.max(Hi_D005))
    Ci_D005_a=Conc_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.2)]
    Ci_D005_s=Conc_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.2)]
    Ci_D005=np.hstack((Ci_D005_a,Ci_D005_s))
    Ice_days_D005.append(len(Ci_D005))

for i in range(1980,2008):
 
    Hi_HC_a=Thickness_HC[(year_HC==i) & (day_HC>180) & (C_HC>=0.8)]
    Hi_HC_s=Thickness_HC[(year_HC==i+1) & (day_HC<180) & (C_HC>=0.8)]
    Hi_HC=np.hstack((Hi_HC_a,Hi_HC_s))
    Annual_max_HC.append(np.max(Hi_HC))
    Ice_days_HC.append(len(Hi_HC))



x1=np.linspace(1976,2059,84)
x2=np.linspace(1981,2008,28)

data_=[x1,Ice_days_A002,Ice_days_B002,Ice_days_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kemi_ABD002.txt',data_.T)

data_=[x1,Ice_days_A005,Ice_days_B005,Ice_days_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kemi_ABD005.txt',data_.T)

data_=[x1,Annual_max_A002,Annual_max_B002,Annual_max_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kemi_ABD002.txt',data_.T)

data_=[x1,Annual_max_A005,Annual_max_B005,Annual_max_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kemi_ABD005.txt',data_.T)


data_=[x2,Ice_days_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kemi_HC.txt',data_.T)

data_=[x2,Annual_max_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kemi_HC.txt',data_.T)



H_002=[(A002[i,3]+B002[i,3]+D002[i,3])/3 for i in range(len(A002))]
Max_H_002=[(Annual_max_A002[i]+Annual_max_B002[i]+Annual_max_D002[i])/3 for i in range(len(Annual_max_A002))]
ice_days_002=[(Ice_days_A002[i]+Ice_days_B002[i]+Ice_days_D002[i])/3 for i in range(len(Annual_max_A002))]

H_005=[(A005[i,3]+B005[i,3]+D005[i,3])/3 for i in range(len(A005))]
Max_H_005=[(Annual_max_A005[i]+Annual_max_B005[i]+Annual_max_D005[i])/3 for i in range(len(Annual_max_A005))]
ice_days_005=[(Ice_days_A005[i]+Ice_days_B005[i]+Ice_days_D005[i])/3 for i in range(len(Annual_max_A005))]

z_mean_002=np.polyfit(x1, ice_days_002, 1)
p_mean_002=np.poly1d(z_mean_002)
print ("IceSeason ABD002 y=%.6fx+(%.6f)"%(z_mean_002[0],z_mean_002[1]))

z_mean_005=np.polyfit(x1, ice_days_005, 1)
p_mean_005=np.poly1d(z_mean_005)
print ("IceSeason ABD005 y=%.6fx+(%.6f)"%(z_mean_005[0],z_mean_005[1]))

z_mean_obs=np.polyfit(Obs[:,0], Obs[:,1],1)
p_mean_obs=np.poly1d(z_mean_obs)
print ("ICE season Obs y=%.6fx+(%.6f)"%(z_mean_obs[0],z_mean_obs[1]))

z_mean_002_H=np.polyfit(x1, Max_H_002, 1)
p_mean_002_H=np.poly1d(z_mean_002_H)
print ("Thickness ABD002 y=%.6fx+(%.6f)"%(z_mean_002_H[0],z_mean_002_H[1]))

z_mean_005_H=np.polyfit(x1, Max_H_005, 1)
p_mean_005_H=np.poly1d(z_mean_005_H)
print ("Thickness ABD005 y=%.6fx+(%.6f)"%(z_mean_005_H[0],z_mean_005_H[1]))

z_mean_obsH=np.polyfit(ObsH[:,0], ObsH[:,1]/100,1)
p_mean_obsH=np.poly1d(z_mean_obs)
print ("Thickness Obs y=%.6fx+(%.6f)"%(z_mean_obsH[0],z_mean_obsH[1]))




########## KALAJOKI ######################
#year_A001,day_A001,c_A001,V_A001=np.loadtxt('Kemi/Ice_season_Kemi_A001.txt', delimiter = ",")
A001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_A001.txt', delimiter = ",")
B001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_B001.txt', delimiter = ",")
D001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_D001.txt', delimiter = ",")

A002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_A002.txt', delimiter = ",")
B002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_B002.txt', delimiter = ",")
D002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_D002.txt', delimiter = ",")

A005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_A005.txt', delimiter = ",")
B005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_B005.txt', delimiter = ",")
D005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_D005.txt', delimiter = ",")

HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kalajoki_HC.txt', delimiter = ",")

Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Kalajoki_Number_of_days_with_ice.csv', delimiter = ",")
ObsH=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Kalajoki_maxH.csv', delimiter = ",")



year_A001=A001[:,0]
day_A001=A001[:,1]
C_A001=A001[:,2]
Vol_A001=A001[:,3]
year_B001=B001[:,0]
day_B001=B001[:,1]
C_B001=B001[:,2]
Vol_B001=B001[:,3]
year_D001=D001[:,0]
day_D001=D001[:,1]
C_D001=D001[:,2]
Vol_D001=D001[:,3]

year_A002=A002[:,0]
day_A002=A002[:,1]
C_A002=A002[:,2]
Vol_A002=A002[:,3]
year_B002=B002[:,0]
day_B002=B002[:,1]
C_B002=B002[:,2]
Vol_B002=B002[:,3]
year_D002=D002[:,0]
day_D002=D002[:,1]
C_D002=D002[:,2]
Vol_D002=D002[:,3]

year_A005=A005[:,0]
day_A005=A005[:,1]
C_A005=A005[:,2]
Vol_A005=A005[:,3]
year_B005=B005[:,0]
day_B005=B005[:,1]
C_B005=B005[:,2]
Vol_B005=B005[:,3]
year_D005=D005[:,0]
day_D005=D005[:,1]
C_D005=D005[:,2]
Vol_D005=D005[:,3]

year_HC=HC[:,0]
day_HC=HC[:,1]
C_HC=HC[:,2]
Vol_HC=HC[:,3]

Year_A002=np.hstack((year_A001,year_A002))
Day_A002=np.hstack((day_A001,day_A002))
Conc_A002=np.hstack((C_A001,C_A002))
Volume_A002=np.hstack((Vol_A001,Vol_A002))
Thickness_A002=Volume_A002/Conc_A002

Year_B002=np.hstack((year_B001,year_B002))
Day_B002=np.hstack((day_B001,day_B002))
Conc_B002=np.hstack((C_B001,C_B002))
Volume_B002=np.hstack((Vol_B001,Vol_B002))
Thickness_B002=Volume_B002/Conc_B002

Year_D002=np.hstack((year_D001,year_D002))
Day_D002=np.hstack((day_D001,day_D002))
Conc_D002=np.hstack((C_D001,C_D002))
Volume_D002=np.hstack((Vol_D001,Vol_D002))
Thickness_D002=Volume_D002/Conc_D002

Year_A005=np.hstack((year_A001,year_A005))
Day_A005=np.hstack((day_A001,day_A005))
Conc_A005=np.hstack((C_A001,C_A005))
Volume_A005=np.hstack((Vol_A001,Vol_A005))
Thickness_A005=Volume_A005/Conc_A005

Year_B005=np.hstack((year_B001,year_B005))
Day_B005=np.hstack((day_B001,day_B005))
Conc_B005=np.hstack((C_B001,C_B005))
Volume_B005=np.hstack((Vol_B001,Vol_B005))
Thickness_B005=Volume_B005/Conc_B005

Year_D005=np.hstack((year_D001,year_D005))
Day_D005=np.hstack((day_D001,day_D005))
Conc_D005=np.hstack((C_D001,C_D005))
Volume_D005=np.hstack((Vol_D001,Vol_D005))
Thickness_D005=Volume_D005/Conc_D005

Thickness_HC=Vol_HC/C_HC

Annual_max_A002=[]
Annual_max_B002=[]
Annual_max_D002=[]
Ice_days_A002=[]
Ice_days_B002=[]
Ice_days_D002=[]

Annual_max_A005=[]
Annual_max_B005=[]
Annual_max_D005=[]
Ice_days_A005=[]
Ice_days_B005=[]
Ice_days_D005=[]

Annual_max_HC=[]
Ice_days_HC=[]


for i in range(1975,2059):

    Hi_A002_a=Thickness_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.8)]
    Hi_A002_s=Thickness_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.8)]
    Hi_A002=np.hstack((Hi_A002_a,Hi_A002_s))
    Annual_max_A002.append(np.max(Hi_A002))
    Ci_A002_a=Conc_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.2)]
    Ci_A002_s=Conc_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.2)]
    Ci_A002=np.hstack((Ci_A002_a,Ci_A002_s))
    Ice_days_A002.append(len(Ci_A002))
  
    Hi_B002_a=Thickness_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.8)]
    Hi_B002_s=Thickness_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.8)]
    Hi_B002=np.hstack((Hi_B002_a,Hi_B002_s))
    Annual_max_B002.append(np.max(Hi_B002))
    Ci_B002_a=Conc_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.2)]
    Ci_B002_s=Conc_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.2)]
    Ci_B002=np.hstack((Ci_B002_a,Ci_B002_s))
    Ice_days_B002.append(len(Ci_B002))

    Hi_D002_a=Thickness_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.8)]
    Hi_D002_s=Thickness_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.8)]
    Hi_D002=np.hstack((Hi_D002_a,Hi_D002_s))
    Annual_max_D002.append(np.max(Hi_D002))
    Ci_D002_a=Conc_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.2)]
    Ci_D002_s=Conc_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.2)]
    Ci_D002=np.hstack((Ci_D002_a,Ci_D002_s))
    Ice_days_D002.append(len(Ci_D002))
    
    Hi_A005_a=Thickness_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.8)]
    Hi_A005_s=Thickness_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.8)]
    Hi_A005=np.hstack((Hi_A005_a,Hi_A005_s))
    Annual_max_A005.append(np.max(Hi_A005))
    Ci_A005_a=Conc_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.2)]
    Ci_A005_s=Conc_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.2)]
    Ci_A005=np.hstack((Ci_A005_a,Ci_A005_s))
    Ice_days_A005.append(len(Ci_A005))
  
    Hi_B005_a=Thickness_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.8)]
    Hi_B005_s=Thickness_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.8)]
    Hi_B005=np.hstack((Hi_B005_a,Hi_B005_s))
    Annual_max_B005.append(np.max(Hi_B005))
    Ci_B005_a=Conc_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.2)]
    Ci_B005_s=Conc_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.2)]
    Ci_B005=np.hstack((Ci_B005_a,Ci_B005_s))
    Ice_days_B005.append(len(Ci_B005))

    Hi_D005_a=Thickness_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.8)]
    Hi_D005_s=Thickness_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.8)]
    Hi_D005=np.hstack((Hi_D005_a,Hi_D005_s))
    Annual_max_D005.append(np.max(Hi_D005))
    Ci_D005_a=Conc_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.2)]
    Ci_D005_s=Conc_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.2)]
    Ci_D005=np.hstack((Ci_D005_a,Ci_D005_s))
    Ice_days_D005.append(len(Ci_D005))

for i in range(1980,2008):
 
    Hi_HC_a=Thickness_HC[(year_HC==i) & (day_HC>180) & (C_HC>=0.8)]
    Hi_HC_s=Thickness_HC[(year_HC==i+1) & (day_HC<180) & (C_HC>=0.8)]
    Hi_HC=np.hstack((Hi_HC_a,Hi_HC_s))
    Annual_max_HC.append(np.max(Hi_HC))
    Ice_days_HC.append(len(Hi_HC))



x1=np.linspace(1976,2059,84)
x2=np.linspace(1981,2008,28)

data_=[x1,Ice_days_A002,Ice_days_B002,Ice_days_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kalajoki_ABD002.txt',data_.T)

data_=[x1,Ice_days_A005,Ice_days_B005,Ice_days_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kalajoki_ABD005.txt',data_.T)

data_=[x1,Annual_max_A002,Annual_max_B002,Annual_max_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kalajoki_ABD002.txt',data_.T)

data_=[x1,Annual_max_A005,Annual_max_B005,Annual_max_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kalajoki_ABD005.txt',data_.T)


data_=[x2,Ice_days_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kalajoki_HC.txt',data_.T)

data_=[x2,Annual_max_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kalajoki_HC.txt',data_.T)



H_002=[(A002[i,3]+B002[i,3]+D002[i,3])/3 for i in range(len(A002))]
Max_H_002=[(Annual_max_A002[i]+Annual_max_B002[i]+Annual_max_D002[i])/3 for i in range(len(Annual_max_A002))]
ice_days_002=[(Ice_days_A002[i]+Ice_days_B002[i]+Ice_days_D002[i])/3 for i in range(len(Annual_max_A002))]

H_005=[(A005[i,3]+B005[i,3]+D005[i,3])/3 for i in range(len(A005))]
Max_H_005=[(Annual_max_A005[i]+Annual_max_B005[i]+Annual_max_D005[i])/3 for i in range(len(Annual_max_A005))]
ice_days_005=[(Ice_days_A005[i]+Ice_days_B005[i]+Ice_days_D005[i])/3 for i in range(len(Annual_max_A005))]

z_mean_002=np.polyfit(x1, ice_days_002, 1)
p_mean_002=np.poly1d(z_mean_002)
print ("IceSeason ABD002 y=%.6fx+(%.6f)"%(z_mean_002[0],z_mean_002[1]))

z_mean_005=np.polyfit(x1, ice_days_005, 1)
p_mean_005=np.poly1d(z_mean_005)
print ("IceSeason ABD005 y=%.6fx+(%.6f)"%(z_mean_005[0],z_mean_005[1]))

z_mean_obs=np.polyfit(Obs[:,0], Obs[:,1],1)
p_mean_obs=np.poly1d(z_mean_obs)
print ("ICE season Obs y=%.6fx+(%.6f)"%(z_mean_obs[0],z_mean_obs[1]))

z_mean_002_H=np.polyfit(x1, Max_H_002, 1)
p_mean_002_H=np.poly1d(z_mean_002_H)
print ("Thickness ABD002 y=%.6fx+(%.6f)"%(z_mean_002_H[0],z_mean_002_H[1]))

z_mean_005_H=np.polyfit(x1, Max_H_005, 1)
p_mean_005_H=np.poly1d(z_mean_005_H)
print ("Thickness ABD005 y=%.6fx+(%.6f)"%(z_mean_005_H[0],z_mean_005_H[1]))

z_mean_obsH=np.polyfit(ObsH[:,0], ObsH[:,1]/100,1)
p_mean_obsH=np.poly1d(z_mean_obs)
print ("Thickness Obs y=%.6fx+(%.6f)"%(z_mean_obsH[0],z_mean_obsH[1]))


########## SÄLGRUND ######################
#year_A001,day_A001,c_A001,V_A001=np.loadtxt('Kemi/Ice_season_Kemi_A001.txt', delimiter = ",")
A001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_A001.txt', delimiter = ",")
B001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_B001.txt', delimiter = ",")
D001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_D001.txt', delimiter = ",")

A002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_A002.txt', delimiter = ",")
B002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_B002.txt', delimiter = ",")
D002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_D002.txt', delimiter = ",")

A005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_A005.txt', delimiter = ",")
B005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_B005.txt', delimiter = ",")
D005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_D005.txt', delimiter = ",")

HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Salgrund_HC.txt', delimiter = ",")

Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Sälgrund_Number_of_days_with_ice.csv', delimiter = ",")
ObsH=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Sälgrund_maxH.csv', delimiter = ",")


year_A001=A001[:,0]
day_A001=A001[:,1]
C_A001=A001[:,2]
Vol_A001=A001[:,3]
year_B001=B001[:,0]
day_B001=B001[:,1]
C_B001=B001[:,2]
Vol_B001=B001[:,3]
year_D001=D001[:,0]
day_D001=D001[:,1]
C_D001=D001[:,2]
Vol_D001=D001[:,3]

year_A002=A002[:,0]
day_A002=A002[:,1]
C_A002=A002[:,2]
Vol_A002=A002[:,3]
year_B002=B002[:,0]
day_B002=B002[:,1]
C_B002=B002[:,2]
Vol_B002=B002[:,3]
year_D002=D002[:,0]
day_D002=D002[:,1]
C_D002=D002[:,2]
Vol_D002=D002[:,3]

year_A005=A005[:,0]
day_A005=A005[:,1]
C_A005=A005[:,2]
Vol_A005=A005[:,3]
year_B005=B005[:,0]
day_B005=B005[:,1]
C_B005=B005[:,2]
Vol_B005=B005[:,3]
year_D005=D005[:,0]
day_D005=D005[:,1]
C_D005=D005[:,2]
Vol_D005=D005[:,3]

year_HC=HC[:,0]
day_HC=HC[:,1]
C_HC=HC[:,2]
Vol_HC=HC[:,3]

Year_A002=np.hstack((year_A001,year_A002))
Day_A002=np.hstack((day_A001,day_A002))
Conc_A002=np.hstack((C_A001,C_A002))
Volume_A002=np.hstack((Vol_A001,Vol_A002))
Thickness_A002=Volume_A002/Conc_A002

Year_B002=np.hstack((year_B001,year_B002))
Day_B002=np.hstack((day_B001,day_B002))
Conc_B002=np.hstack((C_B001,C_B002))
Volume_B002=np.hstack((Vol_B001,Vol_B002))
Thickness_B002=Volume_B002/Conc_B002

Year_D002=np.hstack((year_D001,year_D002))
Day_D002=np.hstack((day_D001,day_D002))
Conc_D002=np.hstack((C_D001,C_D002))
Volume_D002=np.hstack((Vol_D001,Vol_D002))
Thickness_D002=Volume_D002/Conc_D002

Year_A005=np.hstack((year_A001,year_A005))
Day_A005=np.hstack((day_A001,day_A005))
Conc_A005=np.hstack((C_A001,C_A005))
Volume_A005=np.hstack((Vol_A001,Vol_A005))
Thickness_A005=Volume_A005/Conc_A005

Year_B005=np.hstack((year_B001,year_B005))
Day_B005=np.hstack((day_B001,day_B005))
Conc_B005=np.hstack((C_B001,C_B005))
Volume_B005=np.hstack((Vol_B001,Vol_B005))
Thickness_B005=Volume_B005/Conc_B005

Year_D005=np.hstack((year_D001,year_D005))
Day_D005=np.hstack((day_D001,day_D005))
Conc_D005=np.hstack((C_D001,C_D005))
Volume_D005=np.hstack((Vol_D001,Vol_D005))
Thickness_D005=Volume_D005/Conc_D005

Thickness_HC=Vol_HC/C_HC

Annual_max_A002=[]
Annual_max_B002=[]
Annual_max_D002=[]
Ice_days_A002=[]
Ice_days_B002=[]
Ice_days_D002=[]

Annual_max_A005=[]
Annual_max_B005=[]
Annual_max_D005=[]
Ice_days_A005=[]
Ice_days_B005=[]
Ice_days_D005=[]

Annual_max_HC=[]
Ice_days_HC=[]


for i in range(1975,2059):

    Hi_A002_a=Thickness_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.8)]
    Hi_A002_s=Thickness_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.8)]
    Hi_A002=np.hstack((Hi_A002_a,Hi_A002_s))
    if len(Hi_A002)>0:
        Annual_max_A002.append(np.max(Hi_A002))
    else:
        Annual_max_A002.append(0)
    Ci_A002_a=Conc_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.2)]
    Ci_A002_s=Conc_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.2)]
    Ci_A002=np.hstack((Ci_A002_a,Ci_A002_s))
    Ice_days_A002.append(len(Ci_A002))
  
    Hi_B002_a=Thickness_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.8)]
    Hi_B002_s=Thickness_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.8)]
    Hi_B002=np.hstack((Hi_B002_a,Hi_B002_s))
    if len(Hi_B002)>0:
        Annual_max_B002.append(np.max(Hi_B002))
    else:
        Annual_max_B002.append(0)
    Ci_B002_a=Conc_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.2)]
    Ci_B002_s=Conc_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.2)]
    Ci_B002=np.hstack((Ci_B002_a,Ci_B002_s))
    Ice_days_B002.append(len(Ci_B002))

    Hi_D002_a=Thickness_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.8)]
    Hi_D002_s=Thickness_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.8)]
    Hi_D002=np.hstack((Hi_D002_a,Hi_D002_s))
    if len(Hi_D002)>0:
        Annual_max_D002.append(np.max(Hi_D002))
    else:
        Annual_max_D002.append(0)
    Ci_D002_a=Conc_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.2)]
    Ci_D002_s=Conc_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.2)]
    Ci_D002=np.hstack((Ci_D002_a,Ci_D002_s))
    Ice_days_D002.append(len(Ci_D002))
    
    Hi_A005_a=Thickness_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.8)]
    Hi_A005_s=Thickness_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.8)]
    Hi_A005=np.hstack((Hi_A005_a,Hi_A005_s))
    if len(Hi_A005)>0:
        Annual_max_A005.append(np.max(Hi_A005))
    else:
        Annual_max_A005.append(0)
    Ci_A005_a=Conc_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.2)]
    Ci_A005_s=Conc_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.2)]
    Ci_A005=np.hstack((Ci_A005_a,Ci_A005_s))
    Ice_days_A005.append(len(Ci_A005))
  
    Hi_B005_a=Thickness_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.8)]
    Hi_B005_s=Thickness_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.8)]
    Hi_B005=np.hstack((Hi_B005_a,Hi_B005_s))
    if len(Hi_B005)>0:
        Annual_max_B005.append(np.max(Hi_B005))
    else:
        Annual_max_B005.append(0)
    Ci_B005_a=Conc_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.2)]
    Ci_B005_s=Conc_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.2)]
    Ci_B005=np.hstack((Ci_B005_a,Ci_B005_s))
    Ice_days_B005.append(len(Ci_B005))

    Hi_D005_a=Thickness_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.8)]
    Hi_D005_s=Thickness_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.8)]
    Hi_D005=np.hstack((Hi_D005_a,Hi_D005_s))
    if len(Hi_D005)>0:
        Annual_max_D005.append(np.max(Hi_D005))
    else:
        Annual_max_D005.append(0)
    Ci_D005_a=Conc_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.2)]
    Ci_D005_s=Conc_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.2)]
    Ci_D005=np.hstack((Ci_D005_a,Ci_D005_s))
    Ice_days_D005.append(len(Ci_D005))

for i in range(1980,2008):
 
    Hi_HC_a=Thickness_HC[(year_HC==i) & (day_HC>180) & (C_HC>=0.8)]
    Hi_HC_s=Thickness_HC[(year_HC==i+1) & (day_HC<180) & (C_HC>=0.8)]
    Hi_HC=np.hstack((Hi_HC_a,Hi_HC_s))
    Annual_max_HC.append(np.max(Hi_HC))
    Ice_days_HC.append(len(Hi_HC))



x1=np.linspace(1976,2059,84)
x2=np.linspace(1981,2008,28)

data_=[x1,Ice_days_A002,Ice_days_B002,Ice_days_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Salgrund_ABD002.txt',data_.T)

data_=[x1,Ice_days_A005,Ice_days_B005,Ice_days_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Salgrund_ABD005.txt',data_.T)

data_=[x1,Annual_max_A002,Annual_max_B002,Annual_max_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Salgrund_ABD002.txt',data_.T)

data_=[x1,Annual_max_A005,Annual_max_B005,Annual_max_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Salgrund_ABD005.txt',data_.T)


data_=[x2,Ice_days_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Salgrund_HC.txt',data_.T)

data_=[x2,Annual_max_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Salgrund_HC.txt',data_.T)



H_002=[(A002[i,3]+B002[i,3]+D002[i,3])/3 for i in range(len(A002))]
Max_H_002=[(Annual_max_A002[i]+Annual_max_B002[i]+Annual_max_D002[i])/3 for i in range(len(Annual_max_A002))]
ice_days_002=[(Ice_days_A002[i]+Ice_days_B002[i]+Ice_days_D002[i])/3 for i in range(len(Annual_max_A002))]

H_005=[(A005[i,3]+B005[i,3]+D005[i,3])/3 for i in range(len(A005))]
Max_H_005=[(Annual_max_A005[i]+Annual_max_B005[i]+Annual_max_D005[i])/3 for i in range(len(Annual_max_A005))]
ice_days_005=[(Ice_days_A005[i]+Ice_days_B005[i]+Ice_days_D005[i])/3 for i in range(len(Annual_max_A005))]

z_mean_002=np.polyfit(x1, ice_days_002, 1)
p_mean_002=np.poly1d(z_mean_002)
print ("IceSeason ABD002 y=%.6fx+(%.6f)"%(z_mean_002[0],z_mean_002[1]))

z_mean_005=np.polyfit(x1, ice_days_005, 1)
p_mean_005=np.poly1d(z_mean_005)
print ("IceSeason ABD005 y=%.6fx+(%.6f)"%(z_mean_005[0],z_mean_005[1]))

z_mean_obs=np.polyfit(Obs[:,0], Obs[:,1],1)
p_mean_obs=np.poly1d(z_mean_obs)
print ("ICE season Obs y=%.6fx+(%.6f)"%(z_mean_obs[0],z_mean_obs[1]))

z_mean_002_H=np.polyfit(x1, Max_H_002, 1)
p_mean_002_H=np.poly1d(z_mean_002_H)
print ("Thickness ABD002 y=%.6fx+(%.6f)"%(z_mean_002_H[0],z_mean_002_H[1]))

z_mean_005_H=np.polyfit(x1, Max_H_005, 1)
p_mean_005_H=np.poly1d(z_mean_005_H)
print ("Thickness ABD005 y=%.6fx+(%.6f)"%(z_mean_005_H[0],z_mean_005_H[1]))

z_mean_obsH=np.polyfit(ObsH[:,0], ObsH[:,1]/100,1)
p_mean_obsH=np.poly1d(z_mean_obs)
print ("Thickness Obs y=%.6fx+(%.6f)"%(z_mean_obsH[0],z_mean_obsH[1]))



########## KYLMÄPIHLAJA ######################
#year_A001,day_A001,c_A001,V_A001=np.loadtxt('Kemi/Ice_season_Kemi_A001.txt', delimiter = ",")
A001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_A001.txt', delimiter = ",")
B001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_B001.txt', delimiter = ",")
D001=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_D001.txt', delimiter = ",")

A002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_A002.txt', delimiter = ",")
B002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_B002.txt', delimiter = ",")
D002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_D002.txt', delimiter = ",")

A005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_A005.txt', delimiter = ",")
B005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_B005.txt', delimiter = ",")
D005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_D005.txt', delimiter = ",")

HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Ice_season_Kylmapihlaja_HC.txt', delimiter = ",")

Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Kylmäpihlaja_Number_of_days_with_ice.csv', delimiter = ",")
#ObsH=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/ice_seasons/stations/Kummelgrund_maxH.csv', delimiter = ",")

year_A001=A001[:,0]
day_A001=A001[:,1]
C_A001=A001[:,2]
Vol_A001=A001[:,3]
year_B001=B001[:,0]
day_B001=B001[:,1]
C_B001=B001[:,2]
Vol_B001=B001[:,3]
year_D001=D001[:,0]
day_D001=D001[:,1]
C_D001=D001[:,2]
Vol_D001=D001[:,3]

year_A002=A002[:,0]
day_A002=A002[:,1]
C_A002=A002[:,2]
Vol_A002=A002[:,3]
year_B002=B002[:,0]
day_B002=B002[:,1]
C_B002=B002[:,2]
Vol_B002=B002[:,3]
year_D002=D002[:,0]
day_D002=D002[:,1]
C_D002=D002[:,2]
Vol_D002=D002[:,3]

year_A005=A005[:,0]
day_A005=A005[:,1]
C_A005=A005[:,2]
Vol_A005=A005[:,3]
year_B005=B005[:,0]
day_B005=B005[:,1]
C_B005=B005[:,2]
Vol_B005=B005[:,3]
year_D005=D005[:,0]
day_D005=D005[:,1]
C_D005=D005[:,2]
Vol_D005=D005[:,3]

year_HC=HC[:,0]
day_HC=HC[:,1]
C_HC=HC[:,2]
Vol_HC=HC[:,3]

Year_A002=np.hstack((year_A001,year_A002))
Day_A002=np.hstack((day_A001,day_A002))
Conc_A002=np.hstack((C_A001,C_A002))
Volume_A002=np.hstack((Vol_A001,Vol_A002))
Thickness_A002=Volume_A002/Conc_A002

Year_B002=np.hstack((year_B001,year_B002))
Day_B002=np.hstack((day_B001,day_B002))
Conc_B002=np.hstack((C_B001,C_B002))
Volume_B002=np.hstack((Vol_B001,Vol_B002))
Thickness_B002=Volume_B002/Conc_B002

Year_D002=np.hstack((year_D001,year_D002))
Day_D002=np.hstack((day_D001,day_D002))
Conc_D002=np.hstack((C_D001,C_D002))
Volume_D002=np.hstack((Vol_D001,Vol_D002))
Thickness_D002=Volume_D002/Conc_D002

Year_A005=np.hstack((year_A001,year_A005))
Day_A005=np.hstack((day_A001,day_A005))
Conc_A005=np.hstack((C_A001,C_A005))
Volume_A005=np.hstack((Vol_A001,Vol_A005))
Thickness_A005=Volume_A005/Conc_A005

Year_B005=np.hstack((year_B001,year_B005))
Day_B005=np.hstack((day_B001,day_B005))
Conc_B005=np.hstack((C_B001,C_B005))
Volume_B005=np.hstack((Vol_B001,Vol_B005))
Thickness_B005=Volume_B005/Conc_B005

Year_D005=np.hstack((year_D001,year_D005))
Day_D005=np.hstack((day_D001,day_D005))
Conc_D005=np.hstack((C_D001,C_D005))
Volume_D005=np.hstack((Vol_D001,Vol_D005))
Thickness_D005=Volume_D005/Conc_D005

Thickness_HC=Vol_HC/C_HC

Annual_max_A002=[]
Annual_max_B002=[]
Annual_max_D002=[]
Ice_days_A002=[]
Ice_days_B002=[]
Ice_days_D002=[]

Annual_max_A005=[]
Annual_max_B005=[]
Annual_max_D005=[]
Ice_days_A005=[]
Ice_days_B005=[]
Ice_days_D005=[]

Annual_max_HC=[]
Ice_days_HC=[]


for i in range(1975,2059):

    Hi_A002_a=Thickness_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.8)]
    Hi_A002_s=Thickness_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.8)]
    Hi_A002=np.hstack((Hi_A002_a,Hi_A002_s))
    if len(Hi_A002)>0:
        Annual_max_A002.append(np.max(Hi_A002))
    else:
        Annual_max_A002.append(0)
    Ci_A002_a=Conc_A002[(Year_A002==i) & (Day_A002>180) & (Conc_A002>=0.2)]
    Ci_A002_s=Conc_A002[(Year_A002==i+1) & (Day_A002<180) & (Conc_A002>=0.2)]
    Ci_A002=np.hstack((Ci_A002_a,Ci_A002_s))
    Ice_days_A002.append(len(Ci_A002))
  
    Hi_B002_a=Thickness_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.8)]
    Hi_B002_s=Thickness_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.8)]
    Hi_B002=np.hstack((Hi_B002_a,Hi_B002_s))
    if len(Hi_B002)>0:
        Annual_max_B002.append(np.max(Hi_B002))
    else:
        Annual_max_B002.append(0)
    Ci_B002_a=Conc_B002[(Year_B002==i) & (Day_B002>180) & (Conc_B002>=0.2)]
    Ci_B002_s=Conc_B002[(Year_B002==i+1) & (Day_B002<180) & (Conc_B002>=0.2)]
    Ci_B002=np.hstack((Ci_B002_a,Ci_B002_s))
    Ice_days_B002.append(len(Ci_B002))

    Hi_D002_a=Thickness_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.8)]
    Hi_D002_s=Thickness_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.8)]
    Hi_D002=np.hstack((Hi_D002_a,Hi_D002_s))
    if len(Hi_D002)>0:
        Annual_max_D002.append(np.max(Hi_D002))
    else:
        Annual_max_D002.append(0)
    Ci_D002_a=Conc_D002[(Year_D002==i) & (Day_D002>180) & (Conc_D002>=0.2)]
    Ci_D002_s=Conc_D002[(Year_D002==i+1) & (Day_D002<180) & (Conc_D002>=0.2)]
    Ci_D002=np.hstack((Ci_D002_a,Ci_D002_s))
    Ice_days_D002.append(len(Ci_D002))
    
    Hi_A005_a=Thickness_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.8)]
    Hi_A005_s=Thickness_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.8)]
    Hi_A005=np.hstack((Hi_A005_a,Hi_A005_s))
    if len(Hi_A005)>0:
        Annual_max_A005.append(np.max(Hi_A005))
    else:
        Annual_max_A005.append(0)    
    Ci_A005_a=Conc_A005[(Year_A005==i) & (Day_A005>180) & (Conc_A005>=0.2)]
    Ci_A005_s=Conc_A005[(Year_A005==i+1) & (Day_A005<180) & (Conc_A005>=0.2)]
    Ci_A005=np.hstack((Ci_A005_a,Ci_A005_s))
    Ice_days_A005.append(len(Ci_A005))
  
    Hi_B005_a=Thickness_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.8)]
    Hi_B005_s=Thickness_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.8)]
    Hi_B005=np.hstack((Hi_B005_a,Hi_B005_s))
    if len(Hi_B005)>0:
        Annual_max_B005.append(np.max(Hi_B005))
    else:
        Annual_max_B005.append(0)
    Ci_B005_a=Conc_B005[(Year_B005==i) & (Day_B005>180) & (Conc_B005>=0.2)]
    Ci_B005_s=Conc_B005[(Year_B005==i+1) & (Day_B005<180) & (Conc_B005>=0.2)]
    Ci_B005=np.hstack((Ci_B005_a,Ci_B005_s))
    Ice_days_B005.append(len(Ci_B005))

    Hi_D005_a=Thickness_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.8)]
    Hi_D005_s=Thickness_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.8)]
    Hi_D005=np.hstack((Hi_D005_a,Hi_D005_s))
    if len(Hi_D005)>0:
        Annual_max_D005.append(np.max(Hi_D005))
    else:
        Annual_max_D005.append(0)
    Ci_D005_a=Conc_D005[(Year_D005==i) & (Day_D005>180) & (Conc_D005>=0.2)]
    Ci_D005_s=Conc_D005[(Year_D005==i+1) & (Day_D005<180) & (Conc_D005>=0.2)]
    Ci_D005=np.hstack((Ci_D005_a,Ci_D005_s))
    Ice_days_D005.append(len(Ci_D005))

for i in range(1980,2008):
 
    Hi_HC_a=Thickness_HC[(year_HC==i) & (day_HC>180) & (C_HC>=0.8)]
    Hi_HC_s=Thickness_HC[(year_HC==i+1) & (day_HC<180) & (C_HC>=0.8)]
    Hi_HC=np.hstack((Hi_HC_a,Hi_HC_s))
    Annual_max_HC.append(np.max(Hi_HC))
    Ice_days_HC.append(len(Hi_HC))



x1=np.linspace(1976,2059,84)
x2=np.linspace(1981,2008,28)

data_=[x1,Ice_days_A002,Ice_days_B002,Ice_days_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kylmapihlaja_ABD002.txt',data_.T)

data_=[x1,Ice_days_A005,Ice_days_B005,Ice_days_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kylmapihlaja_ABD005.txt',data_.T)

data_=[x1,Annual_max_A002,Annual_max_B002,Annual_max_D002]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kylmapihlaja_ABD002.txt',data_.T)

data_=[x1,Annual_max_A005,Annual_max_B005,Annual_max_D005]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kylmapihlaja_ABD005.txt',data_.T)


data_=[x2,Ice_days_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/Ice_season_L_Kylmapihlaja_HC.txt',data_.T)

data_=[x2,Annual_max_HC]
data_=np.array(data_)
np.savetxt('/home/oikkonea/Documents/SmartSea/scenarios/data/AnnualMax_H_Kylmapihlaja_HC.txt',data_.T)



H_002=[(A002[i,3]+B002[i,3]+D002[i,3])/3 for i in range(len(A002))]
Max_H_002=[(Annual_max_A002[i]+Annual_max_B002[i]+Annual_max_D002[i])/3 for i in range(len(Annual_max_A002))]
ice_days_002=[(Ice_days_A002[i]+Ice_days_B002[i]+Ice_days_D002[i])/3 for i in range(len(Annual_max_A002))]

H_005=[(A005[i,3]+B005[i,3]+D005[i,3])/3 for i in range(len(A005))]
Max_H_005=[(Annual_max_A005[i]+Annual_max_B005[i]+Annual_max_D005[i])/3 for i in range(len(Annual_max_A005))]
ice_days_005=[(Ice_days_A005[i]+Ice_days_B005[i]+Ice_days_D005[i])/3 for i in range(len(Annual_max_A005))]

z_mean_002=np.polyfit(x1, ice_days_002, 1)
p_mean_002=np.poly1d(z_mean_002)
print ("IceSeason ABD002 y=%.6fx+(%.6f)"%(z_mean_002[0],z_mean_002[1]))

z_mean_005=np.polyfit(x1, ice_days_005, 1)
p_mean_005=np.poly1d(z_mean_005)
print ("IceSeason ABD005 y=%.6fx+(%.6f)"%(z_mean_005[0],z_mean_005[1]))

z_mean_obs=np.polyfit(Obs[:,0], Obs[:,1],1)
p_mean_obs=np.poly1d(z_mean_obs)
print ("ICE season Obs y=%.6fx+(%.6f)"%(z_mean_obs[0],z_mean_obs[1]))

z_mean_002_H=np.polyfit(x1, Max_H_002, 1)
p_mean_002_H=np.poly1d(z_mean_002_H)
print ("Thickness ABD002 y=%.6fx+(%.6f)"%(z_mean_002_H[0],z_mean_002_H[1]))

z_mean_005_H=np.polyfit(x1, Max_H_005, 1)
p_mean_005_H=np.poly1d(z_mean_005_H)
print ("Thickness ABD005 y=%.6fx+(%.6f)"%(z_mean_005_H[0],z_mean_005_H[1]))

z_mean_obsH=np.polyfit(ObsH[:,0], ObsH[:,1]/100,1)
p_mean_obsH=np.poly1d(z_mean_obs)
print ("Thickness Obs y=%.6fx+(%.6f)"%(z_mean_obsH[0],z_mean_obsH[1]))
