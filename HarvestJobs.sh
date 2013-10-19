#!/bin/bash
mode=expected
for i in `seq 200 25 3500`
do
   ./merge_ascii.py -f dimuon_ratio_mcmc_limit_m${i}.ascii dimuon_ratio_mcmc_limit_m${i}.[0-9]*.ascii
  rm dimuon_ratio_mcmc_limit_m${i}.[0-9]*_ascii
done
./merge_ascii.py -f dimuon_ratio_${1}.ascii dimuon_ratio_mcmc_limit_m*.ascii
rm dimuon_ratio_mcmc_limit_m*.ascii

