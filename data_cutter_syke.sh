all_series="A001 A002 A005 B001 B002 B005 D001 D002 D005"
all_variables="vosaline votemper"
monthly=true
for series in $all_series; do
    for variable in $all_variables; do
        year="2020"
        depth="0-3"
        in_dir='/scratch/project_2001635/siiriasi/smartsea_data/'${series}'/'
        out_dir='/scratch/project_2001635/siiriasi/smartsea_data/luke_test/'
#        for depth in '0-3' '0-12' '12-80' '80-inf'; do
        for depth in  '0-9'; do
            case $depth in
                "0-3")
                depthlim="deptht,0,0"
                ;;
                "0-12")
                depthlim="deptht,0,5"
                ;;
                "0-9")
                depthlim="deptht,0,2"
                ;;
                "12-80")
                depthlim="deptht,6,19"
                ;;
                "80-inf")
                depthlim="deptht,19,35"
                ;;
            esac
    #    	in_file='SS-GOB_1m_'${year}'0101_'${year}'1231_grid_T.nc'
            files=$(ls ${in_dir}*_1m_*_grid_T.nc --color=none)
            for in_file in $files; do
                year=$(echo $in_file|grep "m_[0-9][0-9][0-9][0-9]" -o| grep "[0-9]*" -o)
                echo $year $series 
#                for m in {0..11}; do
#                    month=$(expr $m + 1)
#                    outf=${out_dir}${variable}_${depth}_${series}_${year}_${month}.nc
#                    ncea  -O -d $depthlim -d time_counter,$m,$m -v $variable $in_file $outf&
#                done
                outf=${out_dir}${variable}_${depth}_${series}_${year}.nc
                ncea  -O -d $depthlim  -v $variable $in_file $outf&
            done
        done
    done
done
