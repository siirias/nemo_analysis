#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:46:17 2020

@author: oikkonea
"""
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
#import netCDF4 as nc
#import mpl_toolkits.basemap as bm
#import datetime as dt
#


#####################################
########## KEMI ######################

###### Length of ice season #########
###### Data #########################
#runs A002, B002, D002
ice_season_002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kemi_ABD002.txt')

year=ice_season_002[:,0]
A002=ice_season_002[:,1]
B002=ice_season_002[:,2]
D002=ice_season_002[:,3]
ave_002=[(A002[i]+B002[i]+D002[i])/3 for i in range(len(A002))] #average of three models, yearly time series

ice_season_005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kemi_ABD005.txt')

#runs A002, B002, D002
A005=ice_season_005[:,1]
B005=ice_season_005[:,2]
D005=ice_season_005[:,3]
ave_005=[(A005[i]+B005[i]+D005[i])/3 for i in range(len(A005))]

#hindcast
HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kemi_HC.txt')
HC_L=HC[:,1]

#observations
Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/observations/Kemi_Ajos_Number_of_days_with_ice.csv', delimiter = ",")
Obs_L=Obs[5:,1] #1981->

##################################################################
##Time series are divided into 10 year periods
# 10a periods of 002 and 005 runs: 1976-1985, 1986-1995, ..., 2046-2055 (3 first of these are identical for 002 and 005)
# For plotting, mean, std, min and max are needed for each 10 a period. For the ensembles (002 and 005 runs), these are determined as mean of three models.
#I.e. std is the mean of std of each model. Min/max is the mean of min and max values of each model.


# 002 runs:
min_002 = [[]]*8
max_002 = [[]]*8
ave_002 = [[]]*8
med_002 = [[]]*8
std_002 = [[]]*8
upper_002 = [[]]*8
lower_002 = [[]]*8
upper_0022 = [[]]*8
lower_0022 = [[]]*8
#diff_002 = [[]]*8
data_002 = [[]]*8
p25_002 = [[]]*8
p75_002 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_002[i]=np.hstack((listA,listB,listD))
    min_002[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_002[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_002[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_002[i]=np.median(data_002[i])
    std_002[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_002[i]=max_002[i]-ave_002[i]  #for plotting, the difference between mean and min/max is needed
    lower_002[i]=ave_002[i]-min_002[i]

    upper_0022[i]=max_002[i]-med_002[i]  #for plotting, the difference between median and min/max is needed
    lower_0022[i]=med_002[i]-min_002[i]

    p25_002[i]=med_002[i]-np.percentile(data_002[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_002[i]=np.percentile(data_002[i],75)-med_002[i]

extremes_002=np.vstack((lower_002,upper_002))
extremes_0022=np.vstack((lower_0022,upper_0022))
p25_75_002=np.vstack((p25_002,p75_002))


# 005 runs:    
min_005 = [[]]*8
max_005 = [[]]*8
ave_005 = [[]]*8
med_005 = [[]]*8
std_005 = [[]]*8
upper_005 = [[]]*8
lower_005 = [[]]*8
upper_0052 = [[]]*8
lower_0052 = [[]]*8
#diff_005 = [[]]*8
data_005 = [[]]*8
p25_005 = [[]]*8
p75_005 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_005[i]=np.hstack((listA,listB,listD))
    min_005[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_005[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_005[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_005[i]=np.median(data_005[i])
    std_005[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_005[i]=max_005[i]-ave_005[i]  #for plotting, the difference between mean and min/max is needed
    lower_005[i]=ave_005[i]-min_005[i]

    upper_0052[i]=max_005[i]-med_005[i]  #for plotting, the difference between median and min/max is needed
    lower_0052[i]=med_005[i]-min_005[i]

    p25_005[i]=med_005[i]-np.percentile(data_005[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_005[i]=np.percentile(data_005[i],75)-med_005[i]

extremes_005=np.vstack((lower_005,upper_005))
extremes_0052=np.vstack((lower_0052,upper_0052))
p25_75_005=np.vstack((p25_005,p75_005))

#Hindcast
# 10a periods of hindacast: 1981-1990, 1991-2000, 2001-2008 (last one only 8 years!)
data_HC = [[]]*3
ave_HC = [[]]*3
med_HC = [[]]*3
std_HC = [[]]*3
min_HC = [[]]*3
max_HC = [[]]*3
upper_HC = [[]]*3
lower_HC = [[]]*3
upper_HC2 = [[]]*3
lower_HC2 = [[]]*3
p25_HC = [[]]*3
p75_HC = [[]]*3 

for i in range(2):
    data_HC[i]=HC_L[i*10:(i+1)*10]
for i in range(2,3):
    data_HC[i]=HC_L[i*10:]     
for i in range(3):
    ave_HC[i]=np.mean(data_HC[i])
    med_HC[i]=np.median(data_HC[i])
    std_HC[i]=np.std(data_HC[i])
    min_HC[i]=np.min(data_HC[i])    
    max_HC[i]=np.max(data_HC[i])    

    print('ave,min,max',ave_HC[i],min_HC[i],max_HC[i])
    upper_HC[i]=max_HC[i]-ave_HC[i]   #for plotting, the difference between mean and min/max is needed
    lower_HC[i]=ave_HC[i]-min_HC[i]

    upper_HC2[i]=max_HC[i]-med_HC[i]   #for plotting, the difference between median and min/max is needed
    lower_HC2[i]=med_HC[i]-min_HC[i]
    
    p25_HC[i]=med_HC[i]-np.percentile(data_HC[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_HC[i]=np.percentile(data_HC[i],75)-med_HC[i]
    print('ave,25,75',ave_HC[i],med_HC[i],np.percentile(data_HC[i],25),np.percentile(data_HC[i],75))

extremes_HC=np.vstack((lower_HC,upper_HC))
extremes_HC2=np.vstack((lower_HC2,upper_HC2))
p25_75_HC=np.vstack((p25_HC,p75_HC))
       
# 10a periods of observations: 1981-1990, 1991-2000, 2001-2010,2011-2020 
data_obs = [[]]*4
ave_obs = [[]]*4
med_obs = [[]]*4
std_obs = [[]]*4
min_obs = [[]]*4
max_obs = [[]]*4
upper_obs = [[]]*4
lower_obs = [[]]*4
upper_obs2 = [[]]*4
lower_obs2 = [[]]*4
p25_obs = [[]]*4
p75_obs = [[]]*4 

for i in range(4):
    data_obs[i]=Obs_L[i*10:(i+1)*10]
    a=data_obs[i]
    if len(a[~np.isnan(a)])>8:
        ave_obs[i]=np.mean(a[~np.isnan(a)])
        med_obs[i]=np.median(a[~np.isnan(a)])
        std_obs[i]=np.std(a[~np.isnan(a)])
        min_obs[i]=np.min(a[~np.isnan(a)])    
        max_obs[i]=np.max(a[~np.isnan(a)])    
    
        upper_obs[i]=max_obs[i]-ave_obs[i]   #for plotting, the difference between mean and min/max is needed
        lower_obs[i]=ave_obs[i]-min_obs[i]

        upper_obs2[i]=max_obs[i]-med_obs[i]   #for plotting, the difference between median and min/max is needed
        lower_obs2[i]=med_obs[i]-min_obs[i]
    
        p25_obs[i]=med_obs[i]-np.percentile(a[~np.isnan(a)],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
        p75_obs[i]=np.percentile(a[~np.isnan(a)],75)-med_obs[i]
        
extremes_obs=np.vstack((lower_obs,upper_obs))
extremes_obs2=np.vstack((lower_obs2,upper_obs2))
p25_75_obs=np.vstack((p25_obs,p75_obs)) 
 
x11=np.linspace(1979.4,2049.6,8)#[1979.5, 1989.5,1999.5, 2009.5,]
x12=np.linspace(1980.4,2050.6,8)#[1980.5, 1990.5,2000.5]
x21=[1984.5, 1994.5, 2004.5]
x22=[1985.5, 1995.5, 2005.5, 2015.5]   




plt.figure(1),plt.clf()
plt.subplot(2,2,1)
#control period
plt.errorbar(x11[:3],ave_002[:3],yerr=std_002[:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')
plt.scatter(x11[:3],ave_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],ave_002[:3],yerr=extremes_002[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3,)
#002
plt.scatter(x11[3:],ave_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],ave_002[3:],yerr=std_002[3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')
plt.errorbar(x11[3:],ave_002[3:],yerr=extremes_002[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3,)
#005
plt.scatter(x12[3:],ave_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],ave_005[3:],yerr=std_005[3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')
plt.errorbar(x12[3:],ave_005[3:],yerr=extremes_005[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3,)
#refrence run
plt.scatter(x21,ave_HC,marker='D',s=60,c='k')
plt.errorbar(x21,ave_HC,yerr=std_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')
plt.errorbar(x21,ave_HC,yerr=extremes_HC,linestyle=' ',elinewidth=1.2,color='k',capsize=3,)
#observations
plt.scatter(x22,ave_obs,marker='D',s=60,c='g')
plt.errorbar(x22,ave_obs,yerr=std_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')
plt.errorbar(x22,ave_obs,yerr=extremes_obs,linestyle=' ',elinewidth=1.2,color='g',capsize=3,)
plt.grid('on')
plt.title('Kemi',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)


plt.figure(2),plt.clf() # Plotting: median, middle 50%, extremes
plt.subplot(2,2,1)
#control period
plt.scatter(x11[:3],med_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],med_002[:3],yerr=p25_75_002[:,:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')#box:50 % of values
plt.errorbar(x11[:3],med_002[:3],yerr=extremes_0022[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3)
#002
plt.scatter(x11[3:],med_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],med_002[3:],yerr=p25_75_002[:,3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')#box:50 % of values
plt.errorbar(x11[3:],med_002[3:],yerr=extremes_0022[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3)
#005
plt.scatter(x12[3:],med_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],med_005[3:],yerr=p25_75_005[:,3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')#box:50 % of values
plt.errorbar(x12[3:],med_005[3:],yerr=extremes_0052[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3)
#refrence run
plt.scatter(x21,med_HC,marker='D',s=60,c='k')
plt.errorbar(x21,med_HC,yerr=p25_75_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')#box:50 % of values
plt.errorbar(x21,med_HC,yerr=extremes_HC2,linestyle=' ',elinewidth=1.2,color='k',capsize=3)
#observations
plt.scatter(x22,med_obs,marker='D',s=60,c='g')
plt.errorbar(x22,med_obs,yerr=p25_75_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')#box:50 % of values
plt.errorbar(x22,med_obs,yerr=extremes_obs2,linestyle=' ',elinewidth=1.2,color='g',capsize=3)
plt.grid('on')
plt.title('Kemi',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)

#####################################
########## KALAJOKI ######################

###### Length of ice season #########
###### Data #########################
#runs A002, B002, D002
ice_season_002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kalajoki_ABD002.txt')

year=ice_season_002[:,0]
A002=ice_season_002[:,1]
B002=ice_season_002[:,2]
D002=ice_season_002[:,3]
ave_002=[(A002[i]+B002[i]+D002[i])/3 for i in range(len(A002))] #average of three models, yearly time series

ice_season_005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kalajoki_ABD005.txt')

#runs A002, B002, D002
A005=ice_season_005[:,1]
B005=ice_season_005[:,2]
D005=ice_season_005[:,3]
ave_005=[(A005[i]+B005[i]+D005[i])/3 for i in range(len(A005))]

#hindcast
HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kalajoki_HC.txt')
HC_L=HC[:,1]

#observations
Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/observations/Kalajoki_Number_of_days_with_ice.csv', delimiter = ",")
Obs_L=Obs[5:,1] #1981->

##################################################################
##Time series are divided into 10 year periods
# 10a periods of 002 and 005 runs: 1976-1985, 1986-1995, ..., 2046-2055 (3 first of these are identical for 002 and 005)
# For plotting, mean, std, min and max are needed for each 10 a period. For the ensembles (002 and 005 runs), these are determined as mean of three models.
#I.e. std is the mean of std of each model. Min/max is the mean of min and max values of each model.



# 002 runs:
min_002 = [[]]*8
max_002 = [[]]*8
ave_002 = [[]]*8
med_002 = [[]]*8
std_002 = [[]]*8
upper_002 = [[]]*8
lower_002 = [[]]*8
upper_0022 = [[]]*8
lower_0022 = [[]]*8
#diff_002 = [[]]*8
data_002 = [[]]*8
p25_002 = [[]]*8
p75_002 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_002[i]=np.hstack((listA,listB,listD))
    min_002[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_002[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_002[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_002[i]=np.median(data_002[i])
    std_002[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_002[i]=max_002[i]-ave_002[i]  #for plotting, the difference between mean and min/max is needed
    lower_002[i]=ave_002[i]-min_002[i]

    upper_0022[i]=max_002[i]-med_002[i]  #for plotting, the difference between median and min/max is needed
    lower_0022[i]=med_002[i]-min_002[i]

    p25_002[i]=med_002[i]-np.percentile(data_002[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_002[i]=np.percentile(data_002[i],75)-med_002[i]

extremes_002=np.vstack((lower_002,upper_002))
extremes_0022=np.vstack((lower_0022,upper_0022))
p25_75_002=np.vstack((p25_002,p75_002))


# 005 runs:    
min_005 = [[]]*8
max_005 = [[]]*8
ave_005 = [[]]*8
med_005 = [[]]*8
std_005 = [[]]*8
upper_005 = [[]]*8
lower_005 = [[]]*8
upper_0052 = [[]]*8
lower_0052 = [[]]*8
#diff_005 = [[]]*8
data_005 = [[]]*8
p25_005 = [[]]*8
p75_005 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_005[i]=np.hstack((listA,listB,listD))
    min_005[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_005[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_005[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_005[i]=np.median(data_005[i])
    std_005[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_005[i]=max_005[i]-ave_005[i]  #for plotting, the difference between mean and min/max is needed
    lower_005[i]=ave_005[i]-min_005[i]

    upper_0052[i]=max_005[i]-med_005[i]  #for plotting, the difference between median and min/max is needed
    lower_0052[i]=med_005[i]-min_005[i]

    p25_005[i]=med_005[i]-np.percentile(data_005[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_005[i]=np.percentile(data_005[i],75)-med_005[i]

extremes_005=np.vstack((lower_005,upper_005))
extremes_0052=np.vstack((lower_0052,upper_0052))
p25_75_005=np.vstack((p25_005,p75_005))

#Hindcast
# 10a periods of hindacast: 1981-1990, 1991-2000, 2001-2008 (last one only 8 years!)
data_HC = [[]]*3
ave_HC = [[]]*3
med_HC = [[]]*3
std_HC = [[]]*3
min_HC = [[]]*3
max_HC = [[]]*3
upper_HC = [[]]*3
lower_HC = [[]]*3
upper_HC2 = [[]]*3
lower_HC2 = [[]]*3
p25_HC = [[]]*3
p75_HC = [[]]*3 

for i in range(2):
    data_HC[i]=HC_L[i*10:(i+1)*10]
for i in range(2,3):
    data_HC[i]=HC_L[i*10:]     
for i in range(3):
    ave_HC[i]=np.mean(data_HC[i])
    med_HC[i]=np.median(data_HC[i])
    std_HC[i]=np.std(data_HC[i])
    min_HC[i]=np.min(data_HC[i])    
    max_HC[i]=np.max(data_HC[i])    

    print('ave,min,max',ave_HC[i],min_HC[i],max_HC[i])
    upper_HC[i]=max_HC[i]-ave_HC[i]   #for plotting, the difference between mean and min/max is needed
    lower_HC[i]=ave_HC[i]-min_HC[i]

    upper_HC2[i]=max_HC[i]-med_HC[i]   #for plotting, the difference between median and min/max is needed
    lower_HC2[i]=med_HC[i]-min_HC[i]
    
    p25_HC[i]=med_HC[i]-np.percentile(data_HC[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_HC[i]=np.percentile(data_HC[i],75)-med_HC[i]
    print('ave,25,75',ave_HC[i],med_HC[i],np.percentile(data_HC[i],25),np.percentile(data_HC[i],75))

extremes_HC=np.vstack((lower_HC,upper_HC))
extremes_HC2=np.vstack((lower_HC2,upper_HC2))
p25_75_HC=np.vstack((p25_HC,p75_HC))
       
# 10a periods of observations: 1981-1990, 1991-2000, 2001-2010,2011-2020 
data_obs = [[]]*4
ave_obs = [[]]*4
med_obs = [[]]*4
std_obs = [[]]*4
min_obs = [[]]*4
max_obs = [[]]*4
upper_obs = [[]]*4
lower_obs = [[]]*4
upper_obs2 = [[]]*4
lower_obs2 = [[]]*4
p25_obs = [[]]*4
p75_obs = [[]]*4 

for i in range(4):
    data_obs[i]=Obs_L[i*10:(i+1)*10]
    a=data_obs[i]
    if len(a[~np.isnan(a)])>8:
        ave_obs[i]=np.mean(a[~np.isnan(a)])
        med_obs[i]=np.median(a[~np.isnan(a)])
        std_obs[i]=np.std(a[~np.isnan(a)])
        min_obs[i]=np.min(a[~np.isnan(a)])    
        max_obs[i]=np.max(a[~np.isnan(a)])    
    
        upper_obs[i]=max_obs[i]-ave_obs[i]   #for plotting, the difference between mean and min/max is needed
        lower_obs[i]=ave_obs[i]-min_obs[i]

        upper_obs2[i]=max_obs[i]-med_obs[i]   #for plotting, the difference between median and min/max is needed
        lower_obs2[i]=med_obs[i]-min_obs[i]
    
        p25_obs[i]=med_obs[i]-np.percentile(a[~np.isnan(a)],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
        p75_obs[i]=np.percentile(a[~np.isnan(a)],75)-med_obs[i]
        
extremes_obs=np.vstack((lower_obs,upper_obs))
extremes_obs2=np.vstack((lower_obs2,upper_obs2))
p25_75_obs=np.vstack((p25_obs,p75_obs)) 


 
x11=np.linspace(1979.4,2049.6,8)#[1979.5, 1989.5,1999.5, 2009.5,]
x12=np.linspace(1980.4,2050.6,8)#[1980.5, 1990.5,2000.5]
x21=[1984.5, 1994.5, 2004.5]
x22=[1985.5, 1995.5, 2005.5, 2015.5]   




plt.figure(1)
plt.subplot(2,2,2)
#control period
plt.errorbar(x11[:3],ave_002[:3],yerr=std_002[:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')
plt.scatter(x11[:3],ave_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],ave_002[:3],yerr=extremes_002[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3,)
#002
plt.scatter(x11[3:],ave_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],ave_002[3:],yerr=std_002[3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')
plt.errorbar(x11[3:],ave_002[3:],yerr=extremes_002[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3,)
#005
plt.scatter(x12[3:],ave_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],ave_005[3:],yerr=std_005[3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')
plt.errorbar(x12[3:],ave_005[3:],yerr=extremes_005[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3,)
#refrence run
plt.scatter(x21,ave_HC,marker='D',s=60,c='k')
plt.errorbar(x21,ave_HC,yerr=std_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')
plt.errorbar(x21,ave_HC,yerr=extremes_HC,linestyle=' ',elinewidth=1.2,color='k',capsize=3,)
#observations
plt.scatter(x22,ave_obs,marker='D',s=60,c='g')
plt.errorbar(x22,ave_obs,yerr=std_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')
plt.errorbar(x22,ave_obs,yerr=extremes_obs,linestyle=' ',elinewidth=1.2,color='g',capsize=3,)
plt.grid('on')
plt.title('Kalajoki',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)

plt.figure(2)
plt.subplot(2,2,2)
#control period
plt.scatter(x11[:3],med_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],med_002[:3],yerr=p25_75_002[:,:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')#box:50 % of values
plt.errorbar(x11[:3],med_002[:3],yerr=extremes_0022[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3)
#002
plt.scatter(x11[3:],med_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],med_002[3:],yerr=p25_75_002[:,3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')#box:50 % of values
plt.errorbar(x11[3:],med_002[3:],yerr=extremes_0022[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3)
#005
plt.scatter(x12[3:],med_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],med_005[3:],yerr=p25_75_005[:,3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')#box:50 % of values
plt.errorbar(x12[3:],med_005[3:],yerr=extremes_0052[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3)
#refrence run
plt.scatter(x21,med_HC,marker='D',s=60,c='k')
plt.errorbar(x21,med_HC,yerr=p25_75_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')#box:50 % of values
plt.errorbar(x21,med_HC,yerr=extremes_HC2,linestyle=' ',elinewidth=1.2,color='k',capsize=3)
#observations
plt.scatter(x22,med_obs,marker='D',s=60,c='g')
plt.errorbar(x22,med_obs,yerr=p25_75_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')#box:50 % of values
plt.errorbar(x22,med_obs,yerr=extremes_obs2,linestyle=' ',elinewidth=1.2,color='g',capsize=3)
plt.grid('on')
plt.title('Kalajoki',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)
#####################################
########## SÄLGRUND ######################

###### Length of ice season #########
###### Data #########################
#runs A002, B002, D002
ice_season_002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Salgrund_ABD002.txt')

year=ice_season_002[:,0]
A002=ice_season_002[:,1]
B002=ice_season_002[:,2]
D002=ice_season_002[:,3]
ave_002=[(A002[i]+B002[i]+D002[i])/3 for i in range(len(A002))] #average of three models, yearly time series

ice_season_005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Salgrund_ABD005.txt')

#runs A002, B002, D002
A005=ice_season_005[:,1]
B005=ice_season_005[:,2]
D005=ice_season_005[:,3]
ave_005=[(A005[i]+B005[i]+D005[i])/3 for i in range(len(A005))]

#hindcast
HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Salgrund_HC.txt')
HC_L=HC[:,1]

#observations
Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/observations/Sälgrund_Number_of_days_with_ice.csv', delimiter = ",")
Obs_L=Obs[5:,1] #1981->

##################################################################
##Time series are divided into 10 year periods
# 10a periods of 002 and 005 runs: 1976-1985, 1986-1995, ..., 2046-2055 (3 first of these are identical for 002 and 005)
# For plotting, mean, std, min and max are needed for each 10 a period. For the ensembles (002 and 005 runs), these are determined as mean of three models.
#I.e. std is the mean of std of each model. Min/max is the mean of min and max values of each model.


# 002 runs:
min_002 = [[]]*8
max_002 = [[]]*8
ave_002 = [[]]*8
med_002 = [[]]*8
std_002 = [[]]*8
upper_002 = [[]]*8
lower_002 = [[]]*8
upper_0022 = [[]]*8
lower_0022 = [[]]*8
#diff_002 = [[]]*8
data_002 = [[]]*8
p25_002 = [[]]*8
p75_002 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_002[i]=np.hstack((listA,listB,listD))
    min_002[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_002[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_002[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_002[i]=np.median(data_002[i])
    std_002[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_002[i]=max_002[i]-ave_002[i]  #for plotting, the difference between mean and min/max is needed
    lower_002[i]=ave_002[i]-min_002[i]

    upper_0022[i]=max_002[i]-med_002[i]  #for plotting, the difference between median and min/max is needed
    lower_0022[i]=med_002[i]-min_002[i]

    p25_002[i]=med_002[i]-np.percentile(data_002[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_002[i]=np.percentile(data_002[i],75)-med_002[i]

extremes_002=np.vstack((lower_002,upper_002))
extremes_0022=np.vstack((lower_0022,upper_0022))
p25_75_002=np.vstack((p25_002,p75_002))


# 005 runs:    
min_005 = [[]]*8
max_005 = [[]]*8
ave_005 = [[]]*8
med_005 = [[]]*8
std_005 = [[]]*8
upper_005 = [[]]*8
lower_005 = [[]]*8
upper_0052 = [[]]*8
lower_0052 = [[]]*8
#diff_005 = [[]]*8
data_005 = [[]]*8
p25_005 = [[]]*8
p75_005 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_005[i]=np.hstack((listA,listB,listD))
    min_005[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_005[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_005[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_005[i]=np.median(data_005[i])
    std_005[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_005[i]=max_005[i]-ave_005[i]  #for plotting, the difference between mean and min/max is needed
    lower_005[i]=ave_005[i]-min_005[i]

    upper_0052[i]=max_005[i]-med_005[i]  #for plotting, the difference between median and min/max is needed
    lower_0052[i]=med_005[i]-min_005[i]

    p25_005[i]=med_005[i]-np.percentile(data_005[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_005[i]=np.percentile(data_005[i],75)-med_005[i]

extremes_005=np.vstack((lower_005,upper_005))
extremes_0052=np.vstack((lower_0052,upper_0052))
p25_75_005=np.vstack((p25_005,p75_005))

#Hindcast
# 10a periods of hindacast: 1981-1990, 1991-2000, 2001-2008 (last one only 8 years!)
data_HC = [[]]*3
ave_HC = [[]]*3
med_HC = [[]]*3
std_HC = [[]]*3
min_HC = [[]]*3
max_HC = [[]]*3
upper_HC = [[]]*3
lower_HC = [[]]*3
upper_HC2 = [[]]*3
lower_HC2 = [[]]*3
p25_HC = [[]]*3
p75_HC = [[]]*3 

for i in range(2):
    data_HC[i]=HC_L[i*10:(i+1)*10]
for i in range(2,3):
    data_HC[i]=HC_L[i*10:]     
for i in range(3):
    ave_HC[i]=np.mean(data_HC[i])
    med_HC[i]=np.median(data_HC[i])
    std_HC[i]=np.std(data_HC[i])
    min_HC[i]=np.min(data_HC[i])    
    max_HC[i]=np.max(data_HC[i])    

    print('ave,min,max',ave_HC[i],min_HC[i],max_HC[i])
    upper_HC[i]=max_HC[i]-ave_HC[i]   #for plotting, the difference between mean and min/max is needed
    lower_HC[i]=ave_HC[i]-min_HC[i]

    upper_HC2[i]=max_HC[i]-med_HC[i]   #for plotting, the difference between median and min/max is needed
    lower_HC2[i]=med_HC[i]-min_HC[i]
    
    p25_HC[i]=med_HC[i]-np.percentile(data_HC[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_HC[i]=np.percentile(data_HC[i],75)-med_HC[i]
    print('ave,25,75',ave_HC[i],med_HC[i],np.percentile(data_HC[i],25),np.percentile(data_HC[i],75))

extremes_HC=np.vstack((lower_HC,upper_HC))
extremes_HC2=np.vstack((lower_HC2,upper_HC2))
p25_75_HC=np.vstack((p25_HC,p75_HC))
       
# 10a periods of observations: 1981-1990, 1991-2000, 2001-2010,2011-2020 
data_obs = [[]]*4
ave_obs = [[]]*4
med_obs = [[]]*4
std_obs = [[]]*4
min_obs = [[]]*4
max_obs = [[]]*4
upper_obs = [[]]*4
lower_obs = [[]]*4
upper_obs2 = [[]]*4
lower_obs2 = [[]]*4
p25_obs = [[]]*4
p75_obs = [[]]*4 

for i in range(4):
    data_obs[i]=Obs_L[i*10:(i+1)*10]
    a=data_obs[i]
    if len(a[~np.isnan(a)])>8:
        ave_obs[i]=np.mean(a[~np.isnan(a)])
        med_obs[i]=np.median(a[~np.isnan(a)])
        std_obs[i]=np.std(a[~np.isnan(a)])
        min_obs[i]=np.min(a[~np.isnan(a)])    
        max_obs[i]=np.max(a[~np.isnan(a)])    
    
        upper_obs[i]=max_obs[i]-ave_obs[i]   #for plotting, the difference between mean and min/max is needed
        lower_obs[i]=ave_obs[i]-min_obs[i]

        upper_obs2[i]=max_obs[i]-med_obs[i]   #for plotting, the difference between median and min/max is needed
        lower_obs2[i]=med_obs[i]-min_obs[i]
    
        p25_obs[i]=med_obs[i]-np.percentile(a[~np.isnan(a)],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
        p75_obs[i]=np.percentile(a[~np.isnan(a)],75)-med_obs[i]
        
extremes_obs=np.vstack((lower_obs,upper_obs))
extremes_obs2=np.vstack((lower_obs2,upper_obs2))
p25_75_obs=np.vstack((p25_obs,p75_obs)) 

x11=np.linspace(1979.4,2049.6,8)#[1979.5, 1989.5,1999.5, 2009.5,]
x12=np.linspace(1980.4,2050.6,8)#[1980.5, 1990.5,2000.5]
x21=[1984.5, 1994.5, 2004.5]
x22=[1985.5, 1995.5, 2005.5, 2015.5]   




plt.figure(1)
plt.subplot(2,2,3)
#control period
plt.errorbar(x11[:3],ave_002[:3],yerr=std_002[:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')
plt.scatter(x11[:3],ave_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],ave_002[:3],yerr=extremes_002[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3,)
#002
plt.scatter(x11[3:],ave_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],ave_002[3:],yerr=std_002[3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')
plt.errorbar(x11[3:],ave_002[3:],yerr=extremes_002[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3,)
#005
plt.scatter(x12[3:],ave_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],ave_005[3:],yerr=std_005[3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')
plt.errorbar(x12[3:],ave_005[3:],yerr=extremes_005[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3,)
#refrence run
plt.scatter(x21,ave_HC,marker='D',s=60,c='k')
plt.errorbar(x21,ave_HC,yerr=std_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')
plt.errorbar(x21,ave_HC,yerr=extremes_HC,linestyle=' ',elinewidth=1.2,color='k',capsize=3,)
#observations
plt.scatter(x22,ave_obs,marker='D',s=60,c='g')
plt.errorbar(x22,ave_obs,yerr=std_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')
plt.errorbar(x22,ave_obs,yerr=extremes_obs,linestyle=' ',elinewidth=1.2,color='g',capsize=3,)
plt.grid('on')
plt.title('Sälgrund',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)

plt.figure(2)
plt.subplot(2,2,3)
#control period
plt.scatter(x11[:3],med_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],med_002[:3],yerr=p25_75_002[:,:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')#box:50 % of values
plt.errorbar(x11[:3],med_002[:3],yerr=extremes_0022[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3)
#002
plt.scatter(x11[3:],med_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],med_002[3:],yerr=p25_75_002[:,3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')#box:50 % of values
plt.errorbar(x11[3:],med_002[3:],yerr=extremes_0022[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3)
#005
plt.scatter(x12[3:],med_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],med_005[3:],yerr=p25_75_005[:,3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')#box:50 % of values
plt.errorbar(x12[3:],med_005[3:],yerr=extremes_0052[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3)
#refrence run
plt.scatter(x21,med_HC,marker='D',s=60,c='k')
plt.errorbar(x21,med_HC,yerr=p25_75_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')#box:50 % of values
plt.errorbar(x21,med_HC,yerr=extremes_HC2,linestyle=' ',elinewidth=1.2,color='k',capsize=3)
#observations
plt.scatter(x22,med_obs,marker='D',s=60,c='g')
plt.errorbar(x22,med_obs,yerr=p25_75_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')#box:50 % of values
plt.errorbar(x22,med_obs,yerr=extremes_obs2,linestyle=' ',elinewidth=1.2,color='g',capsize=3)
plt.grid('on')
plt.title('Sälgrund',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)
#####################################
########## KYLMÄPIHLAJA ######################

###### Length of ice season #########
###### Data #########################
#runs A002, B002, D002
ice_season_002=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kylmapihlaja_ABD002.txt')

year=ice_season_002[:,0]
A002=ice_season_002[:,1]
B002=ice_season_002[:,2]
D002=ice_season_002[:,3]
ave_002=[(A002[i]+B002[i]+D002[i])/3 for i in range(len(A002))] #average of three models, yearly time series

ice_season_005=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kylmapihlaja_ABD005.txt')

#runs A002, B002, D002
A005=ice_season_005[:,1]
B005=ice_season_005[:,2]
D005=ice_season_005[:,3]
ave_005=[(A005[i]+B005[i]+D005[i])/3 for i in range(len(A005))]

#hindcast
HC=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/model_results/Ice_season_L_Kylmapihlaja_HC.txt')
HC_L=HC[:,1]

#observations
Obs=np.loadtxt('/home/oikkonea/Documents/SmartSea/scenarios/data/coastal_stations/observations/Kylmäpihlaja_Number_of_days_with_ice.csv', delimiter = ",")
Obs_L=Obs[5:,1] #1981->

##################################################################
##Time series are divided into 10 year periods
# 10a periods of 002 and 005 runs: 1976-1985, 1986-1995, ..., 2046-2055 (3 first of these are identical for 002 and 005)
# For plotting, mean, std, min and max are needed for each 10 a period. For the ensembles (002 and 005 runs), these are determined as mean of three models.
#I.e. std is the mean of std of each model. Min/max is the mean of min and max values of each model.


# 002 runs:
min_002 = [[]]*8
max_002 = [[]]*8
ave_002 = [[]]*8
med_002 = [[]]*8
std_002 = [[]]*8
upper_002 = [[]]*8
lower_002 = [[]]*8
upper_0022 = [[]]*8
lower_0022 = [[]]*8
#diff_002 = [[]]*8
data_002 = [[]]*8
p25_002 = [[]]*8
p75_002 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_002[i]=np.hstack((listA,listB,listD))
    min_002[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_002[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_002[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_002[i]=np.median(data_002[i])
    std_002[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_002[i]=max_002[i]-ave_002[i]  #for plotting, the difference between mean and min/max is needed
    lower_002[i]=ave_002[i]-min_002[i]

    upper_0022[i]=max_002[i]-med_002[i]  #for plotting, the difference between median and min/max is needed
    lower_0022[i]=med_002[i]-min_002[i]

    p25_002[i]=med_002[i]-np.percentile(data_002[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_002[i]=np.percentile(data_002[i],75)-med_002[i]

extremes_002=np.vstack((lower_002,upper_002))
extremes_0022=np.vstack((lower_0022,upper_0022))
p25_75_002=np.vstack((p25_002,p75_002))


# 005 runs:    
min_005 = [[]]*8
max_005 = [[]]*8
ave_005 = [[]]*8
med_005 = [[]]*8
std_005 = [[]]*8
upper_005 = [[]]*8
lower_005 = [[]]*8
upper_0052 = [[]]*8
lower_0052 = [[]]*8
#diff_005 = [[]]*8
data_005 = [[]]*8
p25_005 = [[]]*8
p75_005 = [[]]*8 

for i in range(8):
    listA=A002[i*10:(i+1)*10]
    listB=B002[i*10:(i+1)*10]
    listD=D002[i*10:(i+1)*10]
    
    data_005[i]=np.hstack((listA,listB,listD))
    min_005[i]=(np.min(listA)+np.min(listB)+np.min(listD))/3   #mean of minimums in three models
    max_005[i]=(np.max(listA)+np.max(listB)+np.max(listD))/3   #mean of maximums in three models
    ave_005[i]=(np.mean(listA)+np.mean(listB)+np.mean(listD))/3
    med_005[i]=np.median(data_005[i])
    std_005[i]=(np.std(listA)+np.std(listB)+np.std(listD))/3   #mean of stds in three models

    upper_005[i]=max_005[i]-ave_005[i]  #for plotting, the difference between mean and min/max is needed
    lower_005[i]=ave_005[i]-min_005[i]

    upper_0052[i]=max_005[i]-med_005[i]  #for plotting, the difference between median and min/max is needed
    lower_0052[i]=med_005[i]-min_005[i]

    p25_005[i]=med_005[i]-np.percentile(data_005[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_005[i]=np.percentile(data_005[i],75)-med_005[i]

extremes_005=np.vstack((lower_005,upper_005))
extremes_0052=np.vstack((lower_0052,upper_0052))
p25_75_005=np.vstack((p25_005,p75_005))

#Hindcast
# 10a periods of hindacast: 1981-1990, 1991-2000, 2001-2008 (last one only 8 years!)
data_HC = [[]]*3
ave_HC = [[]]*3
med_HC = [[]]*3
std_HC = [[]]*3
min_HC = [[]]*3
max_HC = [[]]*3
upper_HC = [[]]*3
lower_HC = [[]]*3
upper_HC2 = [[]]*3
lower_HC2 = [[]]*3
p25_HC = [[]]*3
p75_HC = [[]]*3 

for i in range(2):
    data_HC[i]=HC_L[i*10:(i+1)*10]
for i in range(2,3):
    data_HC[i]=HC_L[i*10:]     
for i in range(3):
    ave_HC[i]=np.mean(data_HC[i])
    med_HC[i]=np.median(data_HC[i])
    std_HC[i]=np.std(data_HC[i])
    min_HC[i]=np.min(data_HC[i])    
    max_HC[i]=np.max(data_HC[i])    

    print('ave,min,max',ave_HC[i],min_HC[i],max_HC[i])
    upper_HC[i]=max_HC[i]-ave_HC[i]   #for plotting, the difference between mean and min/max is needed
    lower_HC[i]=ave_HC[i]-min_HC[i]

    upper_HC2[i]=max_HC[i]-med_HC[i]   #for plotting, the difference between median and min/max is needed
    lower_HC2[i]=med_HC[i]-min_HC[i]
    
    p25_HC[i]=med_HC[i]-np.percentile(data_HC[i],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
    p75_HC[i]=np.percentile(data_HC[i],75)-med_HC[i]
    print('ave,25,75',ave_HC[i],med_HC[i],np.percentile(data_HC[i],25),np.percentile(data_HC[i],75))

extremes_HC=np.vstack((lower_HC,upper_HC))
extremes_HC2=np.vstack((lower_HC2,upper_HC2))
p25_75_HC=np.vstack((p25_HC,p75_HC))
       
# 10a periods of observations: 1981-1990, 1991-2000, 2001-2010,2011-2020 
data_obs = [[]]*4
ave_obs = [[]]*4
med_obs = [[]]*4
std_obs = [[]]*4
min_obs = [[]]*4
max_obs = [[]]*4
upper_obs = [[]]*4
lower_obs = [[]]*4
upper_obs2 = [[]]*4
lower_obs2 = [[]]*4
p25_obs = [[]]*4
p75_obs = [[]]*4 

for i in range(4):
    data_obs[i]=Obs_L[i*10:(i+1)*10]
    a=data_obs[i]
    if len(a[~np.isnan(a)])>8:
        ave_obs[i]=np.mean(a[~np.isnan(a)])
        med_obs[i]=np.median(a[~np.isnan(a)])
        std_obs[i]=np.std(a[~np.isnan(a)])
        min_obs[i]=np.min(a[~np.isnan(a)])    
        max_obs[i]=np.max(a[~np.isnan(a)])    
    
        upper_obs[i]=max_obs[i]-ave_obs[i]   #for plotting, the difference between mean and min/max is needed
        lower_obs[i]=ave_obs[i]-min_obs[i]

        upper_obs2[i]=max_obs[i]-med_obs[i]   #for plotting, the difference between median and min/max is needed
        lower_obs2[i]=med_obs[i]-min_obs[i]
    
        p25_obs[i]=med_obs[i]-np.percentile(a[~np.isnan(a)],25)  #The range where 50 % of values fall in. For plotting, the difference between median and upper/lower limit is needed
        p75_obs[i]=np.percentile(a[~np.isnan(a)],75)-med_obs[i]
        
extremes_obs=np.vstack((lower_obs,upper_obs))
extremes_obs2=np.vstack((lower_obs2,upper_obs2))
p25_75_obs=np.vstack((p25_obs,p75_obs)) 


x11=np.linspace(1979.4,2049.6,8)#[1979.5, 1989.5,1999.5, 2009.5,]
x12=np.linspace(1980.4,2050.6,8)#[1980.5, 1990.5,2000.5]
x21=[1984.5, 1994.5, 2004.5]
x22=[1985.5, 1995.5, 2005.5, 2015.5]   




plt.figure(1)
plt.subplot(2,2,4)
#control period
plt.errorbar(x11[:3],ave_002[:3],yerr=std_002[:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')
plt.scatter(x11[:3],ave_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],ave_002[:3],yerr=extremes_002[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3,)
#002
plt.scatter(x11[3:],ave_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],ave_002[3:],yerr=std_002[3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')
plt.errorbar(x11[3:],ave_002[3:],yerr=extremes_002[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3,)
#005
plt.scatter(x12[3:],ave_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],ave_005[3:],yerr=std_005[3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')
plt.errorbar(x12[3:],ave_005[3:],yerr=extremes_005[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3,)
#refrence run
plt.scatter(x21,ave_HC,marker='D',s=60,c='k')
plt.errorbar(x21,ave_HC,yerr=std_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')
plt.errorbar(x21,ave_HC,yerr=extremes_HC,linestyle=' ',elinewidth=1.2,color='k',capsize=3,)
#observations
plt.scatter(x22,ave_obs,marker='D',s=60,c='g')
plt.errorbar(x22,ave_obs,yerr=std_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')
plt.errorbar(x22,ave_obs,yerr=extremes_obs,linestyle=' ',elinewidth=1.2,color='g',capsize=3,)
plt.grid('on')
plt.title('Kylmäpihlaja',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)


plt.figure(2)
plt.subplot(2,2,4)
#control period
plt.scatter(x11[:3],med_002[:3],marker='D',s=60,c='purple')
plt.errorbar(x11[:3],med_002[:3],yerr=p25_75_002[:,:3],linestyle=' ',elinewidth=5,color='purple',capsize=3,label='Control period')#box:50 % of values
plt.errorbar(x11[:3],med_002[:3],yerr=extremes_0022[:,:3],linestyle=' ',elinewidth=1.2,color='purple',capsize=3)
#002
plt.scatter(x11[3:],med_002[3:],marker='D',s=50,c='b')
plt.errorbar(x11[3:],med_002[3:],yerr=p25_75_002[:,3:],linestyle=' ',elinewidth=5,color='b',capsize=3,label='RCP 4.5')#box:50 % of values
plt.errorbar(x11[3:],med_002[3:],yerr=extremes_0022[:,3:],linestyle=' ',elinewidth=1.2,color='b',capsize=3)
#005
plt.scatter(x12[3:],med_005[3:],marker='D',s=60,c='r')
plt.errorbar(x12[3:],med_005[3:],yerr=p25_75_005[:,3:],linestyle=' ',elinewidth=5,color='r',capsize=3,label='RCP 8.5')#box:50 % of values
plt.errorbar(x12[3:],med_005[3:],yerr=extremes_0052[:,3:],linestyle=' ',elinewidth=1.2,color='r',capsize=3)
#refrence run
plt.scatter(x21,med_HC,marker='D',s=60,c='k')
plt.errorbar(x21,med_HC,yerr=p25_75_HC,linestyle=' ',elinewidth=5,color='k',capsize=3,label='Reference run')#box:50 % of values
plt.errorbar(x21,med_HC,yerr=extremes_HC2,linestyle=' ',elinewidth=1.2,color='k',capsize=3)
#observations
plt.scatter(x22,med_obs,marker='D',s=60,c='g')
plt.errorbar(x22,med_obs,yerr=p25_75_obs,linestyle=' ',elinewidth=5,color='g',capsize=3,label='Observations')#box:50 % of values
plt.errorbar(x22,med_obs,yerr=extremes_obs2,linestyle=' ',elinewidth=1.2,color='g',capsize=3)
plt.grid('on')
plt.title('Kylmäpihlaja',fontsize=17)
plt.legend(ncol=1,fontsize=13)
plt.ylim(0,220)
plt.ylabel('Number of days with ice cover',fontsize=15)
plt.xlabel('Year',fontsize=15)
