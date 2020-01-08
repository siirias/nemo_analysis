clear all; close all;
% Name of simulation. Use for the path of results
%directory='./Run_v170912_tauc_0.2_stat'
directory='./plot_matlab' % Changed by sm_samfr

% Read new model data (=1) or use old saved .mat files (=0)
read_new_modeldata=1;

% Path to Nemo-Nordic-Scobi data

% Information about the model grid
gridres='gob';
dlat=0.0167; latsize=359; dlon=0.0278; lonsize=339;

latstart=61.5; latend=64.0; % Lat plot area (latstart=0 gives full plot area)
lonstart=17; lonend=21 % Lon plot area
lon11 = 17; lat11 = 61.5;


% Start time of period (check that the model data files  exists)
startyear=1975; startmonth=01; startday=1; % Changed by sm_samfr
% Ending time of period
%endyear=2005; endmonth=12;endday=31; % Changed by sm_samfr
endyear=1976; endmonth=12;endday=31; % Changed by sm_samfr

% Information save interval in files
% ann=Annual, files with 12 months data assumed, e.g. NORDIC-NS2_1m_20000101_20001231_ptrc_T.nc
% mon=Monthly, files with 1 months data assumed, e.g. NORDIC-GOB_1m_20000101_20000131_ptrc_T.nc
saveint='ann'; % Changed by sm_samfr

% Select if SGU or only GoB stations
station_plot='GoB'; % (Set = SGU or GoB) 

% Select the plots, the default is "plot all"
default_plot=1; % (Set =0 for manual selection below)

if(default_plot==0) % Changed by sm_samfr
% Plot all
% Good for plotting statistics or mean maps for a selected time periods 
  monthly_maps=1;          % plot monthly mean surface maps for the selected period
  annual_maps=1;           % plot annual mean surface maps for the selected period based on the monthly mean data
  time_plot=1;             % plot time series at certain stations and depths
  pools_plot=1;            % plot time series of mean concentrations calculated from pools
  plot_seasonal=1;         % plot seasonal scatter plots at certain stations and depths for the period
  plot_vertical=1;         % plot vertical mean profiles for the period
  plot_surf_to_depth=1;    % plot horizontal mean surface concentration distribution with water depth
  plot_isoplots=1;         % plot vertical isoplots in time
  plot_isoplots_surface=1; % plot vertical isoplots in time
else
% Define which plots to save (1=plot; 0=dont plot)
% Time and pool plots are good for checking performance of long time series 
 monthly_maps=1;          % plot monthly mean surface maps for the selected period % Changed by sm_samfr
 annual_maps=0;           % plot annual mean surface maps for the selected period based on the monthly mean data
 time_plot=1;             % plot time series at certain stations and depths
 pools_plot=1;            % plot time series of mean concentrations calculated from pools
 plot_seasonal=0;         % plot seasonal scatter plots at certain stations and depths for the period
 plot_vertical=0;         % plot vertical mean profiles for the period
 plot_surf_to_depth=0;    % plot horizontal mean surface concentration distribution with water depth
 plot_isoplots=0;         % plot vertical isoplots in time
 plot_isoplots_surface=0; % plot vertical isoplots in time
end

% Output data will be written into hard coded sub directories

%%%%%%%%%%%%%%%
% Unit conversions in model, see trcnam_scobi.F90
x1=16.*50./(106.*12.011);   % From chlorphyll to nitrogen
x2=1.*50./(106.*12.011);    % From chlorphyll to phosphorus
x3= 15.*50./(106.*12.011); % From chlorophyll to silicate
x4= 1./50.;                % From carbon to chlorophyll
x5= 1./(106.*12.011);      % From carbon to phosphorus
x6= 16./(106.*12.011);     % From carbon to nitrogen
x7= 15./(106.*12.011);     % From carbon to silicate

%%%%%%%%%%%%%%%
% Variables included from the present model 
% Standard list
n=1; 
% model data, file type, names and units for titles, range for c-axis of surface map, tranform unit factor of P or N , var index and       Var nr   
%
variable(n,:)='po4      ';typ(n,:)='ptrc_T';variable_name(n,:)='PO4      ';unit(n,:)='(\mu mol l^-^1) ';cax(n,:)=[0 0.3];xi(n,:)=[1;0;0;0];   n=n+1; %1
variable(n,:)='no3      ';typ(n,:)='ptrc_T';variable_name(n,:)='NO3      ';unit(n,:)='(\mu mol l^-^1) ';cax(n,:)=[0 15] ;xi(n,:)=[0;1;1;0];   n=n+1; %2
variable(n,:)='nh4      ';typ(n,:)='ptrc_T';variable_name(n,:)='NH4      ';unit(n,:)='(\mu mol l^-^1) ';cax(n,:)=[0 1]  ;xi(n,:)=[0;1;1;0];   n=n+1; %3
variable(n,:)='phy1     ';typ(n,:)='ptrc_T';variable_name(n,:)='Phy1     ';unit(n,:)='(mg Chl l^-^1)  ';cax(n,:)=[0 2]  ;xi(n,:)=[x2;0;x1;0]; n=n+1; %4
variable(n,:)='phy2     ';typ(n,:)='ptrc_T';variable_name(n,:)='Phy2     ';unit(n,:)='(mg Chl l^-^1)  ';cax(n,:)=[0 2]  ;xi(n,:)=[x2;0;x1;0]; n=n+1; %5
variable(n,:)='phy3     ';typ(n,:)='ptrc_T';variable_name(n,:)='Phy3     ';unit(n,:)='(mg Chl l^-^1)  ';cax(n,:)=[0 2]  ;xi(n,:)=[x2;0;x1;0]; n=n+1; %6
variable(n,:)='oxy      ';typ(n,:)='ptrc_T';variable_name(n,:)='Oxy      ';unit(n,:)='(ml l^-^1)      ';cax(n,:)=[0 10] ;xi(n,:)=[0;0;0;0];   n=n+1; %7
variable(n,:)='si       ';typ(n,:)='ptrc_T';variable_name(n,:)='Si       ';unit(n,:)='(\mu mol l^-^1) ';cax(n,:)=[0 20] ;xi(n,:)=[0;0;0;0];   n=n+1; %8
variable(n,:)='votemper ';typ(n,:)='grid_T';variable_name(n,:)='Temp     ';unit(n,:)='(^oC)           ';cax(n,:)=[0 10] ;xi(n,:)=[0;0;0;0];   n=n+1; %9
variable(n,:)='vosaline ';typ(n,:)='grid_T';variable_name(n,:)='Sal      ';unit(n,:)='(g kg-^1)       ';cax(n,:)=[0 8]  ;xi(n,:)=[0;0;0;0];   n=n+1; %10
variable(n,:)='zoo      ';typ(n,:)='ptrc_T';variable_name(n,:)='zoo      ';unit(n,:)='(\mug l^-^1)    ';cax(n,:)=[0 100];xi(n,:)=[x5;0;x6;0]; n=n+1; %11
variable(n,:)='det      ';typ(n,:)='ptrc_T';variable_name(n,:)='DetN     ';unit(n,:)='(\mug l^-^1)    ';cax(n,:)=[0 200];xi(n,:)=[x5;0;x6;0]; n=n+1; %12
variable(n,:)='detp     ';typ(n,:)='ptrc_T';variable_name(n,:)='DetP     ';unit(n,:)='(\mug l^-^1)    ';cax(n,:)=[0 200];xi(n,:)=[x5;0;x6;0]; n=n+1; %13
variable(n,:)='fep      ';typ(n,:)='ptrc_T';variable_name(n,:)='Fep      ';unit(n,:)='(\mu mol l^-^1) ';cax(n,:)=[0 0.3];xi(n,:)=[1;0;0;0];   n=n+1; %14
variable(n,:)='opal     ';typ(n,:)='ptrc_T';variable_name(n,:)='Opal     ';unit(n,:)='(\mu mol l^-^1) ';cax(n,:)=[0 2]  ;xi(n,:)=[0;0;0;0];   n=n+1; %15
variable(n,:)='SEDP     ';typ(n,:)='ptrc_T';variable_name(n,:)='SedP     ';unit(n,:)='(mmol m^-^2)    ';cax(n,:)=[0 200];xi(n,:)=[0;0;0;1];   n=n+1; %16
variable(n,:)='SEDN     ';typ(n,:)='ptrc_T';variable_name(n,:)='SedN     ';unit(n,:)='(mmol m^-^2)    ';cax(n,:)=[0 200];xi(n,:)=[0;0;0;0];   n=n+1; %17
variable(n,:)='SEDF     ';typ(n,:)='diad_T';variable_name(n,:)='SedFeP   ';unit(n,:)='(mmol m^-^2)    ';cax(n,:)=[0 200];xi(n,:)=[0;0;0;1];   n=n+1; %18
variable(n,:)='SEDS     ';typ(n,:)='diad_T';variable_name(n,:)='SedS     ';unit(n,:)='(mmol m^-^2)    ';cax(n,:)=[0 200];xi(n,:)=[0;0;0;0];   n=n+1; %19

% Added by sm_samfr
variable(n,:)='uos      ';typ(n,:)='grid_U';variable_name(n,:)='uos      ';unit(n,:)='(m s^-^1)       ';cax(n,:)=[-5 5];xi(n,:)=[0;0;0;0];   n=n+1; %20
variable(n,:)='vos      ';typ(n,:)='grid_V';variable_name(n,:)='vos      ';unit(n,:)='(m s^-^1)       ';cax(n,:)=[-5 5];xi(n,:)=[0;0;0;0];   n=n+1; %21
variable(n,:)='sorunoffs';typ(n,:)='grid_T';variable_name(n,:)='Runoff   ';unit(n,:)='(kg m^-^2 s^-^1)';cax(n,:)=[0 0.2]; xi(n,:)=[0;0;0;0];   n=n+1; %22

% Additional variables not in standard output
%variable(n,:)='sorunoffs';typ(n,:)='grid_T';variable_name(n,:)='Runoff   ';unit(n,:)='(kg m^-^2 s^-^1)';cax(n,:)=[0 10]; xi(n,:)=[0;0;0;0];   n=n+1; %20
%variable(n,:)='don      ';typ(n,:)='ptrc_T';variable_name(n,:)='don      ';unit(n,:)='(\mug l^-^1)    ';cax(n,:)=[0 500]; xi(n,:)=[0;0;0;0];   n=n+1; %21
%variable(n,:)='dop      ';typ(n,:)='ptrc_T';variable_name(n,:)='dop      ';unit(n,:)='(\mug l^-^1)    ';cax(n,:)=[0 200]; xi(n,:)=[0;0;0;0];   n=n+1; %22

% Caclulated sumvariables (use the factors xi to sum up N and P components)
% Note: The numbering of sumvariables follow after the length (variable) and is hard coded in read model data
% To be improved in later versions if other output is required
variable_name(n,:)='TotP     '; unit(n,:)='(mmol m^-^3)    '; cax(n,:)=[0 0.6]; xi(n,:)=[0;0;0;0];   n=n+1;   
variable_name(n,:)='DIN      '; unit(n,:)='(mmol m^-^3)    '; cax(n,:)=[0 15];  xi(n,:)=[0;0;0;0];   n=n+1;   
variable_name(n,:)='TotN     '; unit(n,:)='(mmol m^-^3)    '; cax(n,:)=[0 20];  xi(n,:)=[0;0;0;0];   n=n+1;   
variable_name(n,:)='TotSedP  '; unit(n,:)='(mmol m^-^2)    '; cax(n,:)=[0 200]; xi(n,:)=[0;0;0;0];   n=n+1;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional information used for plotting
monthname=['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];
yr1=num2str(startyear); yr2=num2str(endyear);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional information about the observations data used in func_read_obs_data.m

%%
% Shark data 
% See more information in "path:nemo_nordic_scobi_toolsother/sharkdata/readme"
% how to extract and transform the shark data 
%Data header: 1.year; 2.month; 3.day;4.Siktdjup (m);5.Provtagningsdjup (m);6.Temperatur CTD (C);
%7.Salinitet CTD (o/oo psu);8.Syrgashalt O2 vattenhämtare (ml/l);9.Fosfatfosfor PO4-P (umol/l);
%10.Nitritkväve NO2-N (umol/l);11.Nitratkväve NO3-N (umol/l);12.Ammoniumkväve NH4-N (umol/l);
%13.Total kväve Tot-N (umol/l);14.Silicate SiO3-Si (umol/l);15.Humus (ug/l);16.Löst Organiskt Kol DOC (umol/l);
%17.Klorofyll-a vattenhämtare (ug/l); 18.Totalfosfor (umol/l);

%%
% Pools data from NEST files (personal communication with Oleg Savchuk)
% Annual mean concentrations (mmol/m3)calculated from nutrient pools
% Nitrogen col 2-9 and Phosphorus col 10-17
% Year || BB BS BP0-60m BP60-bot GF GR DS KT ||  BB BS BP0-60m BP60-bot GF GR DS KT


h='Finished config'
