maindir=/scratch/project_2001635/siiriasi/smartsea_data/
for letter in B D; do
    for scenario in 1 2 4 5 7; do
        theset=${letter}00${scenario}
        thedir=${maindir}/${theset}/
        echo packing $theset
        zip -v ${maindir}SmartSea_SCOBI_${theset}.zip ${thedir}/*ptrc* ${thedir}/*diad*
    done
done
