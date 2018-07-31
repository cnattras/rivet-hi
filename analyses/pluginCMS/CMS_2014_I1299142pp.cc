// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2014_I1299142pp : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2014_I1299142pp);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      	
		ChargedFinalState cfs(-2.,2.,1.0*GeV);	
		declare(cfs, "CFS");

		FastJets jets(cfs, FastJets::ANTIKT, 0.3, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
		declare(jets, "jets");
		_h_pp_xe_100_300 = bookHisto1D(1,1,1);
		_h_pp_xe_100_120 = bookHisto1D(2,1,1);
		_h_pp_xe_120_150 = bookHisto1D(3,1,1);
		_h_pp_xe_150_300 = bookHisto1D(4,1,1);
		_h_pp_pt_100_300 = bookHisto1D(5,1,1);
		_h_pp_pt_100_120 = bookHisto1D(6,1,1);
		_h_pp_pt_120_150 = bookHisto1D(7,1,1);
		_h_pp_pt_150_300 = bookHisto1D(8,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

		const double weight = event.weight();
		
		const FastJets & jets = apply<FastJets>(event,"jets");
		Jets allJets = jets.jetsByPt(100.*GeV); ///// For testing 10 GeV Jets!!! For final 100+ GeV
		foreach (const Jet& jet, allJets){
			if (jet.pT() <= 300.*GeV){
				nJet_100_300 += weight;
			
				if (jet.pT() <= 120.*GeV){
					foreach (const Particle& p, jet.particles()){
					
						_h_pp_xe_100_300->fill(xe(jet,p),weight);
						_h_pp_xe_100_120->fill(xe(jet,p),weight);
						_h_pp_pt_100_120->fill(p.pT(),weight);
						_h_pp_pt_100_300->fill(p.pT(),weight);
						nJet_100_120 += weight;
					}
				
				}
				else if (jet.pT() <= 150.*GeV){
					foreach (const Particle& p, jet.particles()){
					
						_h_pp_xe_100_300->fill(xe(jet,p),weight);
						_h_pp_xe_120_150->fill(xe(jet,p),weight);
						_h_pp_pt_120_150->fill(p.pT(),weight);
						_h_pp_pt_100_300->fill(p.pT(),weight);
						nJet_120_150 += weight;
					}
				
				}
				else {
					foreach (const Particle& p, jet.particles()){
					
						_h_pp_xe_100_300->fill(xe(jet,p),weight);
						_h_pp_xe_150_300->fill(xe(jet,p),weight);
						_h_pp_pt_150_300->fill(p.pT(),weight);
						_h_pp_pt_100_300->fill(p.pT(),weight);
						nJet_150_300 += weight;
					}
				
				}
				
			}
		}

    }


    /// Normalise histograms etc., after the run
    void finalize() {
		scale(_h_pp_xe_100_300,1./nJet_100_300);
		scale(_h_pp_pt_100_300,1./nJet_100_300);
		scale(_h_pp_xe_100_120,1./nJet_100_120);
		scale(_h_pp_pt_100_120,1./nJet_100_120);
		scale(_h_pp_xe_120_150,1./nJet_120_150);
		scale(_h_pp_pt_120_150,1./nJet_120_150);
		scale(_h_pp_xe_150_300,1./nJet_150_300);
		scale(_h_pp_pt_150_300,1./nJet_150_300);
		
		_h_pp_xe_100_300->setPath("/CMS_2014_I1299142/h_pp_xe_100_300");
		_h_pp_pt_100_300->setPath("/CMS_2014_I1299142/h_pp_pt_100_300");
		_h_pp_xe_100_120->setPath("/CMS_2014_I1299142/h_pp_xe_100_120");
		_h_pp_pt_100_120->setPath("/CMS_2014_I1299142/h_pp_pt_100_120");
		_h_pp_xe_120_150->setPath("/CMS_2014_I1299142/h_pp_xe_120_150");
		_h_pp_pt_120_150->setPath("/CMS_2014_I1299142/h_pp_pt_120_150");
		_h_pp_xe_150_300->setPath("/CMS_2014_I1299142/h_pp_xe_150_300");
		_h_pp_pt_150_300->setPath("/CMS_2014_I1299142/h_pp_pt_150_300");
		

    }
private:
    //@}
	double z (const Jet& jet, const Particle& part){
			return dot(jet.p3(),part.p3()) / (jet.p3().mod()*jet.p3().mod());
	}
	double xe (const Jet& jet, const Particle& part){
		return log(1./z(jet,part));
	}

    /// @name Histograms
    //@{
    Histo1DPtr _h_pp_xe_100_300;
    Histo1DPtr _h_pp_xe_100_120;
	Histo1DPtr _h_pp_xe_120_150;
	Histo1DPtr _h_pp_xe_150_300;
	Histo1DPtr _h_pp_pt_100_300;
	Histo1DPtr _h_pp_pt_100_120;
	Histo1DPtr _h_pp_pt_120_150;
	Histo1DPtr _h_pp_pt_150_300;
	
    //@}
	double nJet_100_300 =0;
	double nJet_100_120 =0;
	double nJet_120_150 =0;
	double nJet_150_300 =0;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2014_I1299142pp);


}
