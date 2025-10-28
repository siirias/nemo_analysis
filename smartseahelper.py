# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:08:28 2018

@author: siirias
"""
import datetime as dt
import calendar
import numpy as np
import os
import re
import pandas as pd
class smh:
    # 1:d/h 2:startdate 3:enddate 4:grid type (T)
    #date format YYYYMMDD .strftime("%Y%m%d")
    def __init__(self):
        cnf_file='configuration.txt'  #this file on the run directory will set the directories corresponding to wherever this script is run.
        self.root_data_out='.\\'
        self.root_data_in='.\\'        
        if os.path.isfile(cnf_file):  
            params=eval("".join(open(cnf_file).readlines()))
            self.root_data_out=params['root_data_out']
            self.root_data_in=params['root_data_in']        
        self.main_data_folder=self.root_data_in+"\\OUTPUTA001\\"
        self.file_name_format="NORDIC-GOB_1{}_{}_{}_grid_{}.nc"  
        self.grid_type='T'
        self.interval='d' #'d','m','h'
        self.save_interval='year'
        self.reftime=dt.datetime(1950,1,1)
    def give_file_interval(self,unit='days'):
        multiplier=1.0
        if(unit=='months'):
            multiplier=1.0/30.0
        if(unit=='years'):
            multiplier=1.0/365.0
        if(self.save_interval == 'year'):
            return multiplier*365.
        if(self.save_interval == 'day'):
            return multiplier
        if(self.save_interval == 'month'):
            return multiplier*30.0
        return 0.0        
    def filenames_between(self,start, end):
        if self.save_interval=='month':
            filenames=[]
            current_date=start
            while(current_date<end):
                end_date=dt.datetime(current_date.year,\
                                     current_date.month,\
                                     calendar.monthrange(current_date.year,\
                                                         current_date.month)[1]
                                     )
                
                filenames.append(self.file_name_format.format(self.interval,\
                                  current_date.strftime("%Y%m%d"),\
                                  end_date.strftime("%Y%m%d"),\
                                  self.grid_type))
                if(current_date.month<12):
                    current_date=dt.datetime(current_date.year,\
                                             current_date.month+1,\
                                             current_date.day)
                else:
                    current_date=dt.datetime(current_date.year+1,\
                                             1,\
                                             current_date.day)
                
            return filenames
        if self.save_interval=='year':
            filenames=[]
            current_date=start
            while(current_date<end):
                end_date=dt.datetime(current_date.year,\
                                     12,\
                                     calendar.monthrange(current_date.year,\
                                                         12)[1]
                                     )
                
                filenames.append(self.file_name_format.format(self.interval,\
                                  current_date.strftime("%Y%m%d"),\
                                  end_date.strftime("%Y%m%d"),\
                                  self.grid_type))
                current_date=dt.datetime(current_date.year+1,\
                                             1,\
                                             current_date.day)
                
            return filenames
        
        return []
        
    def nemo_time_to_datetime(self,nemo_time):
        return self.reftime + dt.timedelta(seconds=nemo_time)
    def fix_latslons(self,lats,lons):
        #The lats lons have masked values, which causes problems for some plottings,
        #this tries to fix those.
        lats_fix=lats.data
        for i in range(len(lats[:,0])):
            lats_fix[i,:]=lats[i,:].max()
        lons_fix=lons.data
        for i in range(len(lons[0,:])):
            lons_fix[:,i]=lons[:,i].max()            
        #specific gludge, as lons seem to have few extra empty rows.
        for i in range(16):
            lons_fix[:,i]=lons_fix[:,16]
        return lats_fix,lons_fix
    def give_areas(self,lats,lons):
        dlats=np.diff(lats,1,0)
        dlats=np.concatenate((dlats,dlats[0:1,:]),0)
        dlats=dlats*1.852*60. #kilometers

        dlons=np.diff(lons,1,1) #so kilometers
        dlons=np.concatenate((dlons,dlons[:,0:1]),1)
        dlons=dlons*1.852*60.*np.cos((lats/360.0)*(2.0*np.pi)) #so kilometers
        return dlats*dlons
    def latlon_index(self,lat,lon,lats,lons):
        #gives the indices corresponding to given latitude and longitude, on the given lat/lon grid
        return (np.abs(lats[:,0]-lat).argmin(), np.abs(lons[0,:]-lon).argmin())
        
    def give_bottom_values(self, array4d):
        bottom_index=(array4d[0,:,:,:]!=0.0).sum(axis=0)-1 #first axis time, then depth, x,y
        bottom_index[bottom_index<0]=0
        values=array4d[:,0,:,:].copy()
        for i in range(values.shape[1]):
            for j in range(values.shape[2]):
                if ~bottom_index.mask[i,j]:
                    values[:,i,j]=array4d[:,bottom_index[i,j],i,j]
        return values
    def get_bottom(self, grid):
        # grid is supposed to be masked array, Time, D,Lat,Lon
        # or D,lat, lon
        # The idea in this is to shifht the mask one layer up,
        # and find the values which are masked in one (and only one) of
        # these masks. 
        full_shape = grid.shape
        dimensions = len(full_shape)
        if(dimensions == 4): # case with time, D, lat, lon
            bottom_layers = np.zeros((  full_shape[0],\
                                        full_shape[2],
                                        full_shape[3]))
            bottom_layers = np.ma.masked_array(bottom_layers,False)
            mask_roll = np.roll(grid.mask,-1,1) # move mask values one up.
            mask_roll[:,-1,:,:] = True  # And mark bottom most mask as True.
                                       # This to get bottom values if there are no mask at end
            grid.mask = ~(grid.mask ^ mask_roll)
            bottom_layers = np.sum(grid,1)
            values = np.array(np.sum(~grid.mask,1),bool)  # used to get the mask
        else: # case with D, lat, lon
            bottom_layers = np.zeros((  full_shape[1],
                                        full_shape[2]))
            bottom_layers = np.ma.masked_array(bottom_layers,False)
            mask_roll = np.roll(grid.mask,-1,0) # move mask values one up.
            mask_roll[-1,:,:] = True  # And mark bottom most mask as True.
                                       # This to get bottom values if there are no mask at end
            grid.mask = ~(grid.mask ^ mask_roll)
            bottom_layers = np.sum(grid,0)
            values = np.array(np.sum(~grid.mask,0),bool)  # used to get the mask
            
        bottom_layers.mask = ~values
        return bottom_layers

    def get_depth(self, grid, depth_axis, the_depth):
        # grid is supposed to be masked array, Time, D,Lat,Lon
        # or D,lat, lon
        # The idea in this is to shifht the mask one layer up,
        # and find the values which are masked in one (and only one) of
        # these masks. 
        full_shape = grid.shape
        dimensions = len(full_shape)
        layer = np.argmin(np.abs(depth_axis-the_depth))
        if(dimensions == 4): # case with time, D, lat, lon
            the_layer = grid[:,layer,:,:]
        else: # case with D, lat, lon
            the_layer = grid[layer,:,:]
            
        return the_layer
            
    def set_style(self, set_name,alpha=1.0):
        # returns a dictionary that cna be used to
        # define plot style. COlor, linestyles, thickness,
        # Based on the string given. Homogenizes plots when this
        # is used.
        scenario = ""
        if( '001' in set_name or 'hindcast' or 'HISTORY' in set_name):
            scenario = "history"
        if('002' in set_name or 'RCP45' in set_name):
            scenario = "rcp45"
        if('005' in set_name or 'RCP85' in set_name):
            scenario = "rcp85"
        model_type = re.search("^[A-Z]*",set_name).group()
        #execeptions for mean values:
        if(set_name.upper() == "RCP85"):
            model_type = "RCP85"
        if(set_name.upper() == "RCP45"):
            model_type = "RCP45"
        if(set_name.upper() == "CONTROL"):
            model_type = "Control"
        if(set_name.upper() == "REFERENCE"):
            model_type = "Reference"
        if(set_name.upper() == "HINDCAST"):
            model_type = "hindcast"
    
        colors = {
            'A':'b',
            'B':'#ff8c00',
            'C':'c',
            'D':'g',
            'h':'k',
            'RCP':'r',
            'HISTORY':'#66CCEE',
            'hindcast':'#66CCEE',
            'REANALYSIS':'k',
            'RCP45':'b',
            'RCP85':'#ff8c00',
            'Control':'m',
            'Reference':'m'
        }
        scen_styles = {
            'history':'--',
            'rcp45':'-',
            'rcp85':'-'
        }
        line_width = 1.0
        marker = ''
        if(scenario == 'rcp85'):
            line_width=2.0
#        print(set_name,scenario)
#        print(model_type,scenario)
        return {'color':colors[model_type],
                'linestyle':scen_styles[scenario],
                'linewidth':line_width,
                'alpha':alpha,
                'marker':marker}

    def load_boundary_data(self, period = None):
        boundary_data = {}
        change_time = dt.datetime(2006,1,1) # used to cut the forecasts before this
        if(not period): #default period
            period={'min':dt.datetime(1976,1,1), 'max':dt.datetime(2100,1,1)}
        data_dir = self.root_data_in + '\\derived_data\\boundary\\'
    #    for subset in ['boundary_mean','5meter','20meter','80meter','120meter']:
        yearly_means_bnds = {}
        for subset in ['boundary_mean','5meter','80meter']:
            files = os.listdir(data_dir)
            files = [f for f in files if subset in f]
            dat={}
            for f in files:
                set_name=re.search('_([^_]*)\.csv',f).groups()[0]
                dat[set_name]=pd.read_csv(data_dir+f,\
                                     parse_dates=[0])
                dat[set_name]=dat[set_name].set_index('time')
                if(not set_name in yearly_means_bnds.keys()):
                    yearly_means_bnds[set_name]={}
                yearly_means_bnds[set_name][subset] = \
                    dat[set_name].groupby(pd.Grouper(freq='1AS')).mean()
            #calculate the means for History, RCP4.5 and RCP8.5
            dat["Control"] = pd.concat([dat['A001'],dat['B001'],dat['D001']])
            dat["Control"].sort_index(inplace = True)
            dat["RCP45"] = pd.concat([dat['A002'],dat['B002'],dat['D002']])
            dat["RCP45"].sort_index(inplace = True)
            dat["RCP85"] = pd.concat([dat['A005'],dat['B005'],dat['D005']])
            dat["RCP85"].sort_index(inplace = True)
            for s in dat:
                d=dat[s]            
                d = d[(d.index>period['min']) & (d.index<period['max'])]
                if(s == 'hindcast'): # this to cut hindcast in similar shape than control
                    dat[s] = d[(d.index>period['min']) & (d.index<change_time)]
                dat[s] = dat[s].sort_index()
            boundary_data[subset] = dat.copy()
        return boundary_data
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
