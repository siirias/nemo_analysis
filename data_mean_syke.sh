in_dir='/scratch/project_2001635/siiriasi/smartsea_data/syke_test/'
out_dir='/scratch/project_2001635/siiriasi/smartsea_data/syke_test/means/'
series="D002"
variable="votemper"
year="2020"
depth="0-3"
in_files=$(ls $in_dir/*.nc  --color=none)
for f in $in_files; do
    f_name_full=$(basename $f)
    f_extension="${f_name_full##*.}"
    f_name="${f_name_full%.*}"
    f_out_name=${f_name}_mean.${f_extension}
    cdo vertmean $f $out_dir/$f_out_name 
done

#for depth in '0-3' '0-12' '12-80' '80-inf'; do
#	case $depth in
#		"0-3")
#		depthlim="deptht,0,0"
#		;;
#		"0-12")
#		depthlim="deptht,0,5"
#		;;
#		"12-80")
#		depthlim="deptht,6,19"
#		;;
#		"80-inf")
#		depthlim="deptht,19,35"
#		;;
#	esac
#	in_file='SS-GOB_1m_'${year}'0101_'${year}'1231_grid_T.nc'
#	for m in {0..11}; do
#		month=$(expr $m + 1)
#		outf=${out_dir}${variable}_${depth}_${series}_${year}_${month}.nc
#		ncea  -O -d $depthlim -d time_counter,$m,$m -v $variable $in_dir$in_file $outf
#	done
#done
