# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 19:31:58 2021
Huge gludge just to get some correlations out, using two earlier scripts
@author: siirias
"""
import pandas as pd
rerun_series = False
if(rerun_series):
    yearly_means={}
    yearly_means_bnds={}
    print('run model trends')
    exec(open('D:/Data/Sorsat/nemo_analysis/plot_model_trends.py').read())
    print('run forcing trends')
    exec(open('D:/Data/Sorsat/nemo_analysis/plot_forcing_trends.py').read())
    print('done')

correlations = {}
for layer in ('surface','deep'):
    correlations[layer]=[]
    print("####{}! ####".format(layer.upper()))
    for point in yearly_means.keys():
        for forcing in yearly_means_bnds.keys():
            if(layer == 'surface'):
                bnd_depth = '5meter'
                depth = list(yearly_means[point][forcing].keys())[0]
            if(layer == 'deep'):
                bnd_depth = '80meter'
                depth = list(yearly_means[point][forcing].keys())[-1]
            correlation = yearly_means_bnds[forcing][bnd_depth].corrwith(\
                             yearly_means[point][forcing][depth])['vosaline']
            print("{}:{} correlation:{}".format(forcing,point,correlation))
            correlations[layer].append([forcing,point,correlation])

for i in correlations:
    correlations[i] = pd.DataFrame(correlations[i], columns=['set','point','correlation'])
    print(correlations[i].to_latex(caption=i))