### Only dimuon works currently.

Standard procedure on **_lxplus_**:
<pre>
cd CMSSW_X_X_X/src/; cmsenv
svn co svn+ssh://[your_cern_login]@svn.cern.ch/reps/exost
cd exost
source setup/cmssw_setup.[c]sh
cd workdir
git clone git@github.com:zhang8473/twobody.git
cd twobody
<comment>
PWDTMP=${PWD//\//\\\/};sed -i "s/,.*dimuon_20637invpb.root/, 'Dimuon data', ${PWDTMP}\/dimuon_20637invpb.root/g" dimuon_ratio.cfg
</comment>
exost -a workspace -c dimuon_ratio.cfg 
cp myWS.root ws_dimuon_ratio.cfg
root rootlogon.C
</pre>

Inside **_ROOT_**:
<pre>
root [0] .L dimuon.C+
root [1] limit("dimuon","expected",300.0,"m300.0_test",1,10000,500,"",-1.,600)
</pre>

Explatation of **_limit()_** function
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
                Double_t masswindow_width, // events with invmass from peak*(1-masswindow_width) to peak*(1+masswindow_width) will be considered in the profile likelihood calculation, no masswindow cut if this value is less than zero
                Int_t minEvents // the minium events to be used in the profile likelihood calculation, max(minEvents,# of events in the mass window) will be used in the profile likelihood calculation
              )
=========================================================================================*/
</pre>
