stand procedure on lxplus:

1. cd CMSSW_X_X_X/src/; cmsenv
2. git clone git@github.com:zhang8473/exost.git
3. cd exost
4. source setup/cmssw_setup.[c]sh
5. cd workdir
6. git clone git@github.com:zhang8473/twobody.git
7. cd twobody
8. exost -a workspace -c dimuon_ratio.cfg 
9. cp myWS.root ws_dimuon_ratio.cfg
10. root rootlogon.C
11. root [0] .L dimuon.C+
12. root [1] limit("dimuon","expected",300.0,"m300.0_test",1,10000,500,"","")
<pre>
/*=======================================================================================
         limit( std::string channel, // dimuon, dielectron, mumuee, etc
                std::string mode,    // observed, expected, mass limit (extra k-factor uncertainty) 
                Float_t peak,        // resonance mass
                std::string suffix,  // suffix for output file names
                Int_t ntoys,         // number of pseudoexperiments for expected limit
                Int_t mcmc_iter,     // number of MCMC iterations
                Int_t mcmc_burnin,   // number of MCMC burn in steps to be discarded
                std::string inputdir,// directory with workspace files
                std::string plot_name // not using
              )
=========================================================================================*/
</pre>
