#testing

#plots for temperature:
temperature_param_sets =[\
                # basic plots from surface
             {  'var1':'votemper',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

    
                # basic plots from bottom
             {  'var1':'votemper',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

                #comparisons from surface

             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A005',\
                'other_name_marker':'A002',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

                #comparisons from bottom
             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A005',\
                'other_name_marker':'A002',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

             ]


temperature_param_sets_reduced =[\
                # basic plots from surface
             {  'var1':'votemper',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
                # basic plots from bottom
             {  'var1':'votemper',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\
             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

                #comparisons from surface

             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

                #comparisons from bottom
             {  'var1':'votemper',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/temperature/",\
                'cludge':False\
             },\

             ]
#ice parameters

ice_param_sets =[\
             {  'var1':'icecon',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/ice/",\
                'cludge':True\
             },\
             {  'var1':'icecon',\
                'name_marker':'A005',\
                'other_name_marker':'A002',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/ice/",\
                'cludge':True\
             },\
             {  'var1':'icecon',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/ice/A001/",\
                'cludge':True\
             },\
             {  'var1':'icecon',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/ice/A002/",\
                'cludge':True\
             },\
             {  'var1':'icecon',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/ice/A005/",\
                'cludge':True\
             } 
             ]


#salinity parameters

salt_param_sets =[\
                # basic plots from surface
             {  'var1':'vosaline',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\

    
                # basic plots from bottom
             {  'var1':'vosaline',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\

                #comparisons from surface

             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A005',\
                'other_name_marker':'A002',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\

                #comparisons from bottom
             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A005',\
                'other_name_marker':'A002',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\

             ]


salt_param_sets_reduced =[\
                # basic plots from surface
             {  'var1':'vosaline',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
   
                # basic plots from bottom
             {  'var1':'vosaline',\
                'name_marker':'A001',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':False,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
                #comparisons from surface

             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':0,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
                #comparisons from bottom
             {  'var1':'vosaline',\
                'name_marker':'A002',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/salinity/",\
                'cludge':False\
             },\
             ]



