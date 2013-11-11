#!/bin/bash
#mode=expected
mode=observed
for i in `seq 200 25 3500`
do
  if ls dimuon_ratio_mcmc_limit_m${i}.[0-9]*_${mode}*.ascii &> /dev/null; then
     ./merge_ascii.py -f dimuon_ratio_mcmc_limit_m${i}_${mode}.ascii dimuon_ratio_mcmc_limit_m${i}.[0-9]*_${mode}*.ascii
     rm dimuon_ratio_mcmc_limit_m${i}.[0-9]*_${mode}*.ascii
  else
    echo "dimuon_ratio_mcmc_limit_m${i} do not exist"
  fi
done
./merge_ascii.py -f dimuon_ratio_${mode}_${1}.ascii dimuon_ratio_mcmc_limit_m*_${mode}.ascii
rm dimuon_ratio_mcmc_limit_m*_${mode}.ascii

