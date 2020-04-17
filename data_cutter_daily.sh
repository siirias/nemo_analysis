#all_series="A001 A002 A005 B001 B002 B005 D001 D002 D005"
all_series="A001 A002 A005 B001 B002 B005 D001 D002 D005"
all_variables="SST SSS"
monthly=true
for series in $all_series; do
    for variable in $all_variables; do
        in_dir='/scratch/project_2001635/siiriasi/smartsea_data/'${series}'/'
        out_dir='/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
        files=$(ls ${in_dir}*_1d_*_grid_T.nc --color=none)
        for in_file in $files; do
            year=$(echo $in_file|grep "d_[0-9][0-9][0-9][0-9]" -o| grep "[0-9]*" -o)
            echo $year $series 
            outf=${out_dir}${variable}_${series}_${year}.nc
            ncea  -O -v $variable $in_file $outf&
        done
    done
done
