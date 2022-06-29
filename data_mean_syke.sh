in_dir='/scratch/project_2001635/siiriasi/smartsea_data/syke_test/'
out_dir='/scratch/project_2001635/siiriasi/smartsea_data/syke_test/means/'
in_files=$(ls $in_dir  --color=none)
for f in $in_files; do
    if [[ $f == *nc ]];
    then
        f_name_full=$(basename $f)
        f_extension="${f_name_full##*.}"
        f_name="${f_name_full%.*}"
        f_out_name=${f_name}_mean.${f_extension}
        echo ${in_dir} , $f , $f_out_name
        cdo vertmean ${in_dir}$f $out_dir/$f_out_name 
    fi
done

