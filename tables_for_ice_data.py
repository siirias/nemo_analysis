# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 14:10:32 2021

@author: siirias
"""
"""
Observations, reference run and control period runs show a decreasing trend 
in both length of ice season and ice thickness, but with some differences in 
the rate of change. Observations of the ice season length give a trend 
of 
-10 days/decade in Kemi, 
-12 days/decade in Kalajoki and 
-13 days/decade in Sälgrund, but 
no change in Kylmäpihlaja. 

Reference run/control period runs show 
-2/-2 days/decade in Kemi, 
-9/-2 days/decade in Kalajoki, 
-9/0 days/decade in Sälgrund and 
-8/0 days/decade in Kylmäpihlaja.  

Scenarios (RCP4.5/RCP8.5) predict the change to continue: 
-7/-6 days/decade in Kemi, 
-10/-11 days/decade in Kalajoki, 
-9/-11 days/decade in Sälgrun and 
-8/-9 days/decade in Kylmäpihlaja. 

\huom{Oisko nämä nätimpi esittää taulukossa? tulee niin monta lukua, 
että alkaa mennä teksti sekavaksi} 
In addition to decreasing trend, scenarios predict increase in the interannual 
variability of these sea ice variables
"""

import pandas as pd

names = ['Kemi', 'Kalajoki', 'Sälgrund', 'Kylmäpihlaja']
observation_trend_length = [-10,-12,-13,0]
reference_trend_length = [-2,-9,-9,-8]
control_trend_length = [-2,-2,0,0]
rcp45_trend_length = [-7, -10, -9, -8]
rcp85_trend_length = [-6, -11, -11, -9]

column_names = ['observations', 'reference', 'control', 'RCP 4.5', 'RCP 8.5']

tmp = [observation_trend_length, reference_trend_length,
        control_trend_length, rcp45_trend_length, rcp85_trend_length]
data = {}
for i,j in zip(column_names,tmp):
    data[i] = j

dat = pd.DataFrame(data, index = names)
dat.name = "Ice season length trends in days/decade"

print(dat.to_latex(label = "tab:ice_season_length_trends", caption = dat.name))
