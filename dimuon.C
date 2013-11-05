#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>

#include "TFile.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"


#include "RooGlobalFunc.h"
#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "Roo1DTable.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/SequentialProposal.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

Double_t limit( std::string channel,           // dimuon, dielectron...
		std::string mode,              // obsereved, expected, mass limit (extra k-factor uncertainty)
		Float_t peak,
		std::string suffix = "",       // suffix for output file names
		Int_t ntoys = 1,               // number of pseudoexperiments for expected limit
		Int_t mcmc_iter = 100000,      // number of MCMC iterations
		Int_t mcmc_burnin = 100,       // number of MCMC burn in steps to be discarded
		std::string inputdir = "",     // directory with workspace files
		std::string plot_name = "" );

class TwoBody{
  //
  // The class combines multiple channel analyses. A workspace is created
  // with the combined model, data, model config. The class can also call
  // interval calculation routines
  //

public:

  TwoBody();
  ~TwoBody();
  void DimuonRatioLimit( Float_t peak,
			 std::string mode = "observed", // obsereved, expected, mass limit (extra k-factor uncertainty)
			 std::string suffix = "",  // suffix for output file names
			 Int_t nruns = 1,          // number of pseudoexperiments for expected limit
			 Int_t mcmc_iter = 100000, // number of MCMC iterations
			 Int_t mcmc_burnin = 100,  // number of MCMC burn in steps to be discarded
			 std::string inputdir = "",// directory with workspace files
			 std::string plot_name = "" );
  
  Int_t AddWorkspace(std::string filename,
		     std::string ws_name,
		     std::string channel_name,
		     std::string shared_vars);
  ModelConfig prepareDimuonRatioModel( std::string inputdir );
  double printMcmcUpperLimit( double peak, ModelConfig &_mc , std::string filename = "");
  MCMCInterval * GetMcmcInterval(ModelConfig mc,
				 double conf_level,
				 int n_iter,
				 int n_burn,
				 double left_side_tail_fraction,
				 int n_bins);
  LikelihoodInterval * GetPlrInterval( double conf_level, ModelConfig &_mc );

private:
  Int_t CreateDimuonToyMc();
  Double_t GetRandom( std::string pdf, std::string var );
  Int_t FixVariables( std::set<std::string> par );
  Double_t GetPoiUpperSimple(std::string channel, Double_t peak);
  Double_t GetPoiUpper(std::string channel, Double_t peak, ModelConfig &_mc);
  std::map<std::string, double> GetDataRange( RooAbsData * _data, double peak, int goal );
  RooAbsData * SetObservableRange( double peak );
  std::ofstream logfile;
  RooAbsData * data, * realdata;
  RooAbsPdf * model;
  RooArgSet * poi;
  RooArgSet * nuis;
  RooAbsPdf * prior;
  RooWorkspace * ws;
// roostats calculators results
  MCMCInterval * mcInt;
  bool bMcmcConverged;
  TRandom3 r;
};

double TwoBody::printMcmcUpperLimit( double peak, ModelConfig &_mc, std::string filename ){
  //
  // print out the upper limit on the first Parameter of Interest
  //

  std::string _legend = "[TwoBody::printMcmcUpperLimit]: ";

  char buf[1024];

  double _limit = numeric_limits<double>::max();

  if (mcInt){
    //mc.SetWorkspace(*ws);
    RooRealVar * firstPOI = (RooRealVar*) _mc.GetParametersOfInterest()->first();
    _limit = mcInt->UpperLimit(*firstPOI);
    std::cout << "\n95% upper limit on " << firstPOI->GetName() << " is : " <<
      _limit << endl;
    
    if (bMcmcConverged){
      sprintf(buf, "%7.1f   %7.6f", peak, _limit);
    }
    else{
      sprintf(buf, "# %7.1f   %7.6f  # MCMC did not converge", peak, _limit);
    }

  }
  else{
    sprintf(buf, "# MCMC did not converge");
  }

  if (filename.size()!=0){
    
    std::ofstream aFile;
    
    // append to file if exists
    aFile.open(filename.c_str(), std::ios_base::app);
    
    aFile << buf << std::endl;
    
    // close outfile here so it is safe even if subsequent iterations crash
    aFile.close();
    
  }
    
  return _limit;
}

std::map<std::string, double>
TwoBody::GetDataRange( RooAbsData * _data, 
		       double peak,
		       int goal ){
  //
  // Estimate the reduced observable range so either
  //  - ~goal events are used or
  //  - minimal range (so not to kill signal efficiency)
  // for individual channel dataset
  //

  std::string legend = "[TwoBody::GetDataRange]: ";

  int iGoal       = goal; // we want that many events in the range
  double sig_low  = peak*0.8;  // signal box that must be in the range
  double sig_high = peak*1.2; // signal box that must be in the range

  double _total = _data->sumEntries();

  // get data in a vector
  std::vector<double> v_data;
  for (int i=0; i!=_total; ++i){
    RooRealVar * _var = (RooRealVar *)(_data->get(i)->first());
    v_data.push_back(_var->getVal());
  }

  int iTotal = v_data.size();
  std::cout << legend << "data vector size: " << iTotal << "; want [" <<sig_low<<","<<sig_high<<"]"<<std::endl;

  // sort data vector
  std::sort(v_data.begin(), v_data.end());

  // find highest point in the signal box
  int iPeak = iTotal;
  while( iPeak != 0 ){
    --iPeak;
    if ( v_data[iPeak] < sig_high ) break;
  }
  int iSigHigh = iPeak; 
  if (v_data[iSigHigh]>sig_high) iSigHigh = -1;    // all data higher than the signal box
  if (v_data[iSigHigh]<sig_low) iSigHigh = iTotal; // all data lower than the signal box

  // find lowest point in the signal box
  iPeak = 0;
  while( iPeak != (iTotal-1) ){
    if ( v_data[iPeak] > sig_low ) break;
    ++iPeak;
  }
  int iSigLow = iPeak;
  std::cout << legend << "all data lower than the signal box: " << iSigLow << std::endl;
  if (v_data[iSigLow]<sig_low) iSigLow = iTotal; // all data lower than the signal box
  if (v_data[iSigLow]>sig_high) iSigLow = -1;    // all data higher than the signal box

  // now expand the range starting from signal box if necessary
  int iSignal = iSigHigh-iSigLow+1;
  if ( iSigHigh==iTotal || iSigLow==-1 ) iSignal = 0;
  while (iSignal<iGoal){
    int iLow  = std::max(iSigLow-1,0);
    int iHigh = std::min(iSigHigh+1,iTotal-1);
    
    double dLow = 0.0;
    if ( iLow >= 0 ) dLow = std::max(0.0, sig_low - v_data[iLow]);

    double dHigh = 0.0;
    if ( iHigh < iTotal ) dHigh = std::max(0.0, v_data[iHigh] - sig_high);

    if ( dLow < dHigh || iHigh <= iSigHigh ){
      iSigLow = iLow;
      sig_high += sig_low - v_data[iSigLow];
      sig_low = v_data[iSigLow];
    }
    else{
      iSigHigh = iHigh;
      sig_low -= v_data[iSigHigh] - sig_high;
      sig_high = v_data[iSigHigh];
    }
    ++iSignal;
  }

  std::cout << legend << "signal box index range: [" << iSigLow << ", " << iSigHigh << "]" << std::endl;
  std::cout << legend << "signal box range: [" << sig_low << ", " << sig_high << "]" << std::endl;
  std::cout << legend << "events in signal box: " << iSignal << std::endl;

  // events above the signal box
  int iAbove = iTotal - iSigHigh - 1;
  if ( iAbove < 0 ) iAbove = 0;
  if ( iAbove > iTotal ) iAbove= iTotal;

  // events below the signal box
  int iBelow = iSigLow;
  if ( iBelow < 0 ) iBelow = 0;

  std::cout << legend << "events below signal box: " << iBelow << std::endl;
  std::cout << legend << "events above signal box: " << iAbove << std::endl;

  
  // prepare return map
  std::map<std::string, double> _mres;
  _mres["low"] = sig_low;
  _mres["high"] = sig_high;

  return _mres;
}

RooAbsData * TwoBody::SetObservableRange( double peak ){
  //
  // Reduce the observable range so ~400 events are used
  // for the full combined dataset
  //

  std::string legend = "[TwoBody::SetObservableRange]: ";

  int iTotal;

  std::map<std::string, double> _range;

  iTotal = (int)data->sumEntries();
  _range = GetDataRange( data, peak, 600 );
  char buf[256];

  // FIXME: hardcoded POI name throughout

  // change the range
  ws->var("mass")->setRange(_range["low"], _range["high"]);
  ws->var("mass")->Print();

  // replace data
  sprintf(buf, "mass>%f && mass<%f", _range["low"], _range["high"]);
  RooAbsData * _data = data->reduce( RooFit::Cut(buf) );
  // correct the nbkg constraint accordingly
  double _nbkg = ws->var("nbkg_est_dimuon")->getVal();
  ws->var("nbkg_est_dimuon")->setVal(_nbkg*_data->sumEntries()/(double)(iTotal));
  ws->var("nbkg_est_dimuon")->Print();

  data->Print();
  _data->Print();
  delete data;
  return _data;
}

TwoBody::TwoBody(){

  std::string legend = "[TwoBody::TwoBody]: ";

  ws = new RooWorkspace("ws");
  //mc.SetName("mc");
  //mc.SetTitle("model_config");
  data = 0;
  model = 0;
  poi = 0;
  nuis = 0;
  prior = 0;

  mcInt = 0;

  bMcmcConverged = false;

  // set random seed
  r.SetSeed();
  UInt_t _seed = r.GetSeed();
  UInt_t _pid = gSystem->GetPid();
  std::cout << legend << "Random seed: " << _seed << std::endl;
  std::cout << legend << "PID: " << _pid << std::endl;
  _seed = 31*_seed+_pid;

  std::cout << legend << "New random seed (31*seed+pid): " << _seed << std::endl;
  r.SetSeed(_seed);

  // set RooFit random seed (it has a private copy)
  RooRandom::randomGenerator()->SetSeed(_seed);

  // open log file
  logfile.open("twobody.log");
}

TwoBody::~TwoBody(){
  
  std::string _legend = "[TwoBody::~TwoBody]: ";

  logfile.close();
  delete ws;
  delete data;
  delete model;
  delete poi;
  delete nuis;
  delete prior;
  delete mcInt;

}

Int_t TwoBody::FixVariables( std::set<std::string> par ){
  //
  // Set all RooRealVars except <par> to be constants 
  //

  Int_t _fixed = 0;

  RooArgSet _vars = ws->allVars();

  TIterator * iter = _vars.createIterator();

  for(TObject * _obj = iter->Next(); _obj; _obj = iter->Next() ){
    std::string _name = _obj->GetName();
    if (par.find(_name) == par.end()){
      RooRealVar * _var = (RooRealVar *)( _vars.find(_name.c_str()) );
      _var->setConstant(kTRUE);
      ++_fixed;
    }
  }
  delete iter;

  return _fixed;
}

ModelConfig TwoBody::prepareDimuonRatioModel( std::string inputdir ){
  //
  // prepare workspace and ModelConfig for the dimuon xsec ratio limit
  //
  std::string _legend = "[TwoBody::prepareDimuonRatioModel]: ";

  std::string _infile = inputdir+"ws_dimuon_ratio.root";

  AddWorkspace(_infile.c_str(),
	       "myWS",
	       "dimuon",
	       "peak,mass,ratio");

  ws->pdf("model_dimuon")->SetName("model");

  //ws->Print();
  // set all vars to const except <par>
  std::set<std::string> par;
  par.insert("mass");
  par.insert("ratio");
  par.insert("beta_nsig_dimuon");  
  // par.insert("beta_nbkg_dimuon"); 
  //par.insert("beta_mass_dimuon"); 
  FixVariables(par); 
  // POI
  RooArgSet sPoi( *(ws->var("ratio")) );

  // nuisance
  RooArgSet sNuis( *(ws->var("beta_nsig_dimuon"))
		   //*(ws->var("beta_nbkg_dimuon")),
            //*(ws->var("beta_mass_dimuon")) 
		   );

  // observables
  RooArgSet sObs( *(ws->var("mass")) );

  // prior
  
 // ModelConfig
  ModelConfig _mc("mc",ws);
  _mc.SetPdf(*(ws->pdf("model")));
  _mc.SetParametersOfInterest( sPoi );
  _mc.SetPriorPdf( *(ws->pdf("prior_dimuon")) );
  _mc.SetNuisanceParameters( sNuis );
  _mc.SetObservables( sObs );
  return _mc;
}

LikelihoodInterval *TwoBody::GetPlrInterval( double conf_level , ModelConfig &_mc){
  //
  // Profile likelihood ratio interval calculations
  //

  cerr<<"Print ModelConfig: "<<endl;
  _mc.Print();

  ProfileLikelihoodCalculator plc(*data, _mc);
  plc.SetConfidenceLevel(conf_level);
  return plc.GetInterval();
}

double TwoBody::GetPoiUpperSimple(std::string channel, Double_t peak){
  //
  // Estimate a good value for the upper boundary of the range of POI
  // using just data in a window corresponding to the signal region
  //

  std::string legend = "[TwoBody::GetPoiUpperSimple]: ";

  char buf[128];

  // special handling needed for single-channel workspaces
  bool b_single_channel;
  if (channel.size()==0) b_single_channel = true;  else b_single_channel = false;

  double _range = 1.0;

  // get data yield under the signal peak
  ws->var("peak")->setVal(peak);

  if (b_single_channel) sprintf(buf,"width");
  else sprintf(buf,"width_%s",channel.c_str());
  double _width = 0;//ws->function(buf)->getVal();

  if (b_single_channel) sprintf(buf,"sigma");
  else sprintf(buf,"sigma_%s",channel.c_str());
  double _sigma = ws->function(buf)->getVal();

  double c_min = peak - sqrt(_width*_width + _sigma*_sigma);
  double c_max = peak + sqrt(_width*_width + _sigma*_sigma);

  sprintf(buf, "mass>=%f && mass<=%f", c_min, c_max);
  double n_count = data->sumEntries( buf );

  // ad-hoc fix when there are no events in window
  if (n_count < 0.3) n_count = 0.3;
  std::cout << legend << "event yield in data in the window: "
	    << n_count << std::endl;

  // compute the corresponding POI range
  if (b_single_channel) sprintf(buf,"nsig_scale");
  else sprintf(buf,"nsig_scale_%s",channel.c_str());
  double _nsig_scale = ws->var(buf)->getVal();

  if (b_single_channel) sprintf(buf,"nz");
  else sprintf(buf,"nz_%s",channel.c_str());
  double _nz = ws->var(buf)->getVal();

  if (b_single_channel) sprintf(buf,"eff");
  else sprintf(buf,"eff_%s",channel.c_str());
  double _eff = ws->function(buf)->getVal();

  double n_excess = 3.0 * sqrt(n_count)/0.68; // let signal excess be ~ uncertainty on BG
  _range = n_excess / _nsig_scale / _nz / _eff;

  return _range;
}



Double_t TwoBody::GetPoiUpper(std::string channel, Double_t peak, ModelConfig &_mc){
  //
  // Estimate a good value for the upper boundary of the range of POI
  //

  std::string legend = "[TwoBody::GetPoiUpper]: ";

  double _range = -1.0;

  std::cout << legend << "doing a rough estimate of the POI range" << std::endl;
  // adding auto-channel option (multi) while
  // preserving backwards compatibility
  // if single channel
  _range = GetPoiUpperSimple(channel, peak);
  
  std::cout << legend << "crude estimate for poi range (x3): "
	    << 3.0*_range << std::endl;
  std::cout << legend 
	    << "this will be used if the profile likelihood ratio estimate fails"
	    << std::endl;
  std::cout << legend 
	    << "will try to estimate POI range better with profile likelihood ratio limit now"
	    << std::endl;

  Double_t result = 0.1;

  // estimate limit with profile likelihood ratio and
  // set the range to 3 times the limit
  
  // query intervals
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooRealVar * _poi = (RooRealVar *)_mc.GetParametersOfInterest()->first();
  double upper_limit = GetPlrInterval(0.95, _mc)->UpperLimit( *_poi );
  RooMsgService::instance().setGlobalKillBelow(msglevel);
  
  // safety in case upper limit == 0
  if (upper_limit<std::numeric_limits<double>::min()){
    upper_limit = _range;
  }

  result = 3.0*upper_limit;

  return result;
}

Double_t TwoBody::GetRandom( std::string pdf, std::string var ){
  //
  // generates a random number using a pdf in the workspace
  //
  
  // generate a dataset with one entry
  if (ws!=NULL) {
    RooRealVar * _par = ws->var(var.c_str());
    if (_par!=NULL) {
      RooAbsPdf * _pdf=ws->pdf(pdf.c_str());
      /*
	RooPlot* xframe = _par->frame(Title("p.d.f")) ;
      _pdf->plotOn(xframe);
      TCanvas* c = new TCanvas("test","rf101_basics",800,400) ;
      gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
      c->SaveAs("syst_nbkg.pdf");
      */
      if (_pdf!=NULL) return _pdf->generate(*_par, 1)->get(0)->getRealValue(var.c_str(),0);
      else {
	std::cerr<<"Cannot find RooPdf:"<<pdf<<std::endl;
      }
    }
    else std::cerr<<"Cannot find RooVar:"<<var<<std::endl;
  }
  else std::cerr<<"[BUG]workspace is deleted??"<<var<<std::endl;
  return 0;
}

Int_t TwoBody::AddWorkspace(std::string filename,
			    std::string ws_name,
			    std::string channel_name,
			    std::string shared_vars){
  //
  // Load a single channel model and data from a workspace
  //

  std::string _legend = "[TwoBody::AddWorkspace]: ";

   // load workspace from a file
  TFile _file(filename.c_str(), "read");
  RooWorkspace * _ws = (RooWorkspace *)_file.Get( ws_name.c_str() );

  // get the single channel model PDF
  RooAbsPdf * _model = _ws->pdf("model"),
            * _prior = _ws->pdf("prior");

  // import the channel model PDF into the combined workspace
  // add prefix channel_name to all nodes except shared_vars
  ws->import( RooArgSet(*_model,*_prior),
     	      RooFit::RenameAllNodes( channel_name.c_str() ),
    	      RooFit::RenameAllVariablesExcept(channel_name.c_str(), shared_vars.c_str()) );
  delete data; 
  data = new RooDataSet( *(RooDataSet *) _ws->data("data") );
  data->changeObservableName("vertex_m","mass");

  realdata = new RooDataSet( *(RooDataSet *) data );
  //new RooDataSet( *(RooDataSet *) ws->data("data") );
  realdata->SetName("RealData");
  //realdata->changeObservableName("vertex_m","mass");

  ws->import( *data );

  ws->Print();
  //  for (int i=0; i!=100; ++i){
  //  RooRealVar * _var = (RooRealVar *)(data->get(i)->first());
  //  cerr<<_var->getVal()<<",";
  //}
  _file.Close();

  return 0;
}

Int_t TwoBody::CreateDimuonToyMc( void ){
  //
  // generate a toy di-muon dataset with systematics
  // set mData accordingly
  //

  // generate expected number of events from its uncertainty
  //RooDataSet * _ds = ws->pdf("syst_nbkg_dimuon")->generate(*ws->var("beta_nbkg_dimuon"), 1);
  //Double_t _ntoy = ((RooRealVar *)(_ds->get(0)->first()))->getVal() * (ws->var("nbkg_est_dimuon")->getVal());
  //delete _ds;

  Double_t _beta = GetRandom("syst_nbkg_dimuon", "beta_nbkg_dimuon");
  //  Double_t _kappa = ws->var("nbkg_kappa_dimuon")->getVal();
  Double_t _nbkg_est = ws->var("nbkg_est_dimuon")->getVal();
  //Double_t _ntoy = pow(_kappa,_beta) * _nbkg_est;
  Double_t _ntoy = _beta * _nbkg_est;
 
  Int_t _n = r.Poisson(_ntoy);
  _n=Int_t(_nbkg_est);
  // all nuisance parameters:
  //   beta_nsig_dimuon, 
  //   beta_nbkg_dimuon,
  //   lumi_nuis

  // create dataset
  RooRealVar * _mass = ws->var("mass");
  RooArgSet _vars(*_mass);

  RooAbsPdf * _pdf = ws->pdf("bkgpdf_dimuon");

  RooAbsPdf::GenSpec * _spec = _pdf->prepareMultiGen(_vars,
						     Name("toys"),
						     NumEvents(_n),
						     Extended(kFALSE),
						     Verbose(kTRUE));

  //RooPlot* xframe = _mass->frame(Title("Gaussian p.d.f.")) ;
  //realdata->plotOn(xframe,LineColor(kRed),MarkerColor(kRed));

  delete data;
  data = _pdf->generate(*_spec); // class member
  delete _spec;

  //data->plotOn(xframe);
  //TCanvas* c = new TCanvas("test","rf101_basics",800,400) ;
  //gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
  //c->SaveAs("test.pdf");

  Int_t n_generated_entries = (Int_t)(data->sumEntries());

  // debug
  std::cout << "!!!!!!!!!!!!!! _beta = " << _beta << std::endl;
  //std::cout << "!!!!!!!!!!!!!! _kappa = " << _kappa << std::endl;
  std::cout << "!!!!!!!!!!!!!! _nbkg_est = " << _nbkg_est << std::endl;
  std::cout << "!!!!!!!!!!!!!! _ntoy     = " << _ntoy << std::endl;
  std::cout << "!!!!!!!!!!!!!! _n        = " << _n    << std::endl;
  std::cout << "!!!!!!!!!!!!!! n_generated_entries = " << n_generated_entries    << std::endl;
  return n_generated_entries;
}


MCMCInterval * TwoBody::GetMcmcInterval(ModelConfig _mc,
					double conf_level,
					int n_iter,
					int n_burn,
					double left_side_tail_fraction,
					int n_bins){
  // use MCMCCalculator  (takes about 1 min)
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit

  std::string legend = "[TwoBody::GetMcmcInterval]: ";

  /*
  RooFitResult* fit = ws->pdf("model")->fitTo(*data,Save());
  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  ProposalFunction* pf = ph.GetProposalFunction();
  */

  // FIXME: testing: this proposal function seems fairly robust
  SequentialProposal sp(0.5);

  MCMCCalculator mcmc( *data, _mc );
  mcmc.SetConfidenceLevel(conf_level);
  mcmc.SetNumIters(n_iter);          // Metropolis-Hastings algorithm iterations

  // FIXME: testing: different proposal function
  //mcmc.SetProposalFunction(*pf);
  mcmc.SetProposalFunction(sp);

  mcmc.SetNumBurnInSteps(n_burn); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(left_side_tail_fraction);
  mcmc.SetNumBins(n_bins);
  
  // FIXME: testing good initial values - don't seem to do anything different
  //ws->var("ratio")->setVal(0.01);
  //ws->var("beta_nsig_dielectron")->setRange(-3.0, 3.0);
  //ws->var("beta_nbkg_dielectron")->setRange(-3.0, 3.0);
  //ws->var("beta_mass_dielectron")->setRange(-3.0, 3.0);

//mcInt = mcmc.GetInterval();
  try {
    delete mcInt;
    mcInt = mcmc.GetInterval();
  } catch ( std::length_error &ex) {
    mcInt = 0;
  }
  cerr<<legend<<endl;

  // check if limit makes sense
  bMcmcConverged = false; // default
  if (mcInt){
    RooRealVar * p_first_poi = (RooRealVar*) _mc.GetParametersOfInterest()->first();
    double poi_limit = mcInt->UpperLimit(*p_first_poi);
    double u_poi_min  = p_first_poi->getMin();
    double u_poi_max  = p_first_poi->getMax();
    double u_poi_gap = (u_poi_max-poi_limit)/(u_poi_max-u_poi_min);
    std::cout << legend << "POI upper limit: " << poi_limit << std::endl;
    std::cout << legend << "POI range: [" << u_poi_min << ", " << u_poi_max << "]" << std::endl;
    std::cout << legend << "POI upper gap (fraction of range): " << u_poi_gap << std::endl;
    if (u_poi_gap<0.2){
      std::cout << legend 
		<< "POI limit too close to the upper boundary, MCMC probably failed!!!" << std::endl;
      std::cout << legend
		<< "returning interval and setting fail flag" << std::endl;
      bMcmcConverged = false;
    }
    else{
      bMcmcConverged = true;
    }
  }
  else std::cout << "No interval found!" << std::endl;
  
  return mcInt;
}

void TwoBody::DimuonRatioLimit( Float_t peak,
				std::string mode,      // obsereved, expected, mass limit (extra k-factor uncertainty)
				std::string suffix,    // suffix for output file names
				Int_t ntoys,           // number of pseudoexperiments for expected limit
				Int_t mcmc_iter,
				Int_t mcmc_burnin,
				std::string inputdir,  // directory with workspace files
				std::string plot_name ){
  //
  // limit on ratio = xsec(Z')/xsec(Z) in dimuon channel
  //

  std::string _legend = "[TwoBody::DimuonRatioLimit]: ";

  ModelConfig _mc = prepareDimuonRatioModel(inputdir);  

  // set model parameters and data
  ws->var("peak")->setVal(peak);

  int pe_counter = 0;
  std::vector<Double_t> _limits;
  while (pe_counter < ntoys){
    
    if ( mode.find("expected") != std::string::npos ){
      std::cout << _legend << std::endl;
      std::cout << _legend << "this is pseudoexperiment " << pe_counter+1 << " of " << ntoys << std::endl;
      std::cout << _legend << "for the expected limit estimate" << std::endl;
      std::cout << _legend << std::endl;
      // prepare PE data
      CreateDimuonToyMc();
      if (!pe_counter) data=SetObservableRange(peak);
    }
    else { //  "regular" observed limit
      
      std::cout << _legend << std::endl;
      std::cout << _legend << "calculating an observed limit..." << std::endl;
      std::cout << _legend << "I will do it " << ntoys << " times, so one can average. " << pe_counter+1 << " of " << ntoys << std::endl;
      std::cout << _legend << std::endl;
      if (!pe_counter) data=SetObservableRange(peak);
      //ntoys = 1;
    }

    // change POI range
    _mc.SetWorkspace(*ws);
    double poiUpper = GetPoiUpper("dimuon", peak, _mc);
    std::cout << _legend << "setting POI range to [0; " << poiUpper << "]" << std::endl;
    ws->var("ratio")->setRange(0.0, poiUpper);
    // FIXME: try to restrict the range of observable
    
    mcInt = GetMcmcInterval(_mc,
			    0.95,
			    mcmc_iter,
			    mcmc_burnin,
			    0.0,
			    100);
    
    std::string _outfile = "dimuon_ratio_mcmc_limit_" + suffix + ".ascii";
    printMcmcUpperLimit( peak, _mc, _outfile);

    ++pe_counter;
  } // end of while
 
  return;
}

Double_t limit( std::string channel, // dimuon, dielectron, mumuee, etc
		std::string mode,    // observed, expected, mass limit (extra k-factor uncertainty)
		Float_t peak,        // resonance mass
		std::string suffix,  // suffix for output file names
		Int_t ntoys,         // number of pseudoexperiments for expected limit
		Int_t mcmc_iter,     // number of MCMC iterations
		Int_t mcmc_burnin,   // number of MCMC burn in steps to be discarded
		std::string inputdir,// directory with workspace files
		std::string plot_name ){

  // time it
  TStopwatch t;
  t.Start();
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  Double_t limit = -1.0;  

  TwoBody manager;

  //dimuon single channel ratio limit
  if (channel.find("dimuon") != std::string::npos ){
    manager.DimuonRatioLimit(peak, mode, suffix,
			     ntoys, mcmc_iter, mcmc_burnin,
			     inputdir, plot_name);
  }

  t.Print();
  
  return limit;
}

void twobody(void){} // dummy
void dimuon(void){} // dummy
