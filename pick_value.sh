# Needs at least following:
# module load cdo/1.9.3
# module load netcdf

variable='SSS'
in_files='/arch/smartsea/analysis/D005/SS-GOB_1d_*_grid_T.nc'
out_dir='/arch/smartsea/analysis/tmp/'
out_file='D005JustSSS.nc'
in_files=$(ls $in_files)
for in_file in $in_files; do
    in_name=$(basename ${in_file})
    cdo -selname,$variable $in_file ${out_dir}$in_name'tmp'
done
cdo -mergetime ${out_dir}'*.nctmp' $out_file
rm ${out_dir}*nctmp
