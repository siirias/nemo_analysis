for f in $@; do
    ncrename -v Maximum\ Bottom\ Stress,MaximumBottomStress -v Mean\ Bottom\ Stress,MeanBottomStress $f
done
