input_dir='/scratch/project_2002540/siiriasi/smartsea_data/'
ouput_dir='/scratch/project_2002540/siiriasi/smartsea_data/derived_data/kvarken/'
for sets in B002 B005 D002 D005 A002 A005; do
#current_set='REANALYSIS'
   current_set=$sets
   files=$(ls ${input_dir}${current_set}/*1d*T.nc --color=none)
   for f in $files; do
       f="$(basename -- $f)"
       echo $f
       cdo select,name=SBS ${input_dir}${current_set}/$f ${input_dir}${current_set}/tmp_$f
   done
   cdo mergetime ${input_dir}${current_set}/tmp_* ${input_dir}${current_set}/tmp2_joined_set.nc
   rm ${input_dir}${current_set}/tmp_*
    cdo sellonlatbox,17.8,24.0,62.5,64.0 ${input_dir}${current_set}/tmp2_joined_set.nc ${ouput_dir}/kvarken_${current_set}.nc
   #cdo zonmean ${input_dir}${current_set}/tmp2_joined_set.nc ${input_dir}${current_set}/sbs_merman_${current_set}.nc
   rm ${input_dir}${current_set}/tmp2_joined_set.nc
done
