in_dir='/scratch/project_2001635/siiriasi/smartsea_data/luke_test/'
out_dir='/scratch/project_2001635/siiriasi/smartsea_data/luke_test/means/'
in_files=$(ls $in_dir/*.nc  --color=none)
for f in $in_files; do
    f_name_full=$(basename $f)
    f_extension="${f_name_full##*.}"
    f_name="${f_name_full%.*}"
    f_out_name=${f_name}_mean.${f_extension}
    cdo vertmean $f $out_dir/$f_out_name 
done

