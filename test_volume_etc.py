import imp
import iris
import gsw
import siri_omen.utility
import numpy as np
import siri_omen.nemo_reader as nrd

data_folder = '/arch/smartsea/analysis/new_REANALYSIS/'
#data_folder = '/lustre/tmp/siirias/o/tmp//new_REANALYSIS/'
file_name = 'SS-GOB_1m_19800101_19801231_grid_T.nc'

#data_folder = '/lustre/tmp/siirias/o/tmp//REANALYSIS_SMHI/'
#data_folder = '/arch/smartsea/analysis/REANALYSIS_SMHI/'
#file_name = 'NORDIC-GOB_1m_19850101_19850131_grid_T.nc'

sou = imp.reload(siri_omen.utility)
nrd = imp.reload(siri_omen.nemo_reader)

temperature = iris.load(data_folder+file_name, 'potential_temperature')[0]
temperature = nrd.remove_null_indices(temperature, fill_value=0.0)
nrd.fix_cube_coordinates(temperature)
temperature = temperature[0, :, :, :]

salinity = iris.load(data_folder+file_name, 'salinity')[0]
salinity = nrd.remove_null_indices(salinity, fill_value=0.0)
nrd.fix_cube_coordinates(salinity)
salinity = sou.abs_sal_from_pract_sal(salinity)
salinity = salinity[0, :, :, :]

# TODO: convert potential temperature to conservative temperature
# TODO: convert salinity to absolute salinity
#

# test the volume is right
v = sou.cube_volumes(temperature)
depth = v.coord('depth').points
thicknesses = sou.cube_cell_thicknesses(v,return_dictionary=True)
#for i in range(36):
#    print(
#        ("depth {:5.1f} m, \tcells {:5d}, \t"
#        + "volume {:4.0f} km^3\tthickness {:4.1f}")
#        .format(depth[i],
#                np.sum(~v.data.mask[i, :, :]),
#                np.sum(v.data[i, :, :]*1e-9),
#                thicknesses[depth[i]]))

total_volume = np.sum(v.data)*1e-9  # 1e-9 as w want km^3
# actual volume of Gulf of Bothnia is 6369 km^3
difference = 6369.0-total_volume

# a test separately for pressure
pressure = sou.cube_pressure(salinity).data
print("\n")
for d,p in zip(depth,pressure[:,40,40]):
    print("{:.1f} m~{:.1f} dbar | ".format(d,p), end="")


print("\n\nTotal volume is {:.0f} km^3".format(total_volume))
print("Difference to literature {:.0f} km^3, {:.1f} %"
      .format(difference, 100.0*difference/total_volume))

# test that density is reasonable.

d = sou.cube_density(salinity, temperature)
total_mass = np.sum(v.data*d.data)*1e-3  
# 1e-12 as v = m^3  d = kg/m^-3 and we want kg -> T
# actual water mass of Gulf of Bothnia is around 6433 T
print("\n\nTotal mass is {:.0f} GT".format(total_mass*1e-9))

# salinity (should still change from psu -> g/g 
total_salt = np.sum(salinity.data*v.data*d.data)*1e-6
# 1e-6 as salinity g/kg, and want T
average_salt = total_salt/total_mass
print("Total salt {:.1f} GT\tAverage salinity {:.2} g/kg".format(total_salt*1e-9,average_salt*1e3))



# should test for heat content also
heat_content = sou.cube_heat_content(salinity, temperature)
#heat_content_tot = np.sum(heat_content.data)
heat_content_tot = np.sum(heat_content.collapsed(['longitude', 'latitude'],iris.analysis.SUM).data)
print("\n")
print("Total heat content {:.3e} J".format(heat_content_tot))
print("\nSome sanity checks, dT of 1 C of 1 kg requires ~ 4 kJ ")

