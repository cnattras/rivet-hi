// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2015_I1343112pp : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2015_I1343112pp);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
		const FinalState fs(Cuts::pT > 0.15*GeV);	  
		declare (fs, "fs");
		declare(FastJets(fs, FastJets::ANTIKT, 0.2), "Jets2");
		_histraapt5c1 = bookHisto1D(6, 1, 1);
		_histraapt5c3 = bookHisto1D(7, 1, 1);
    }
    


    /// Perform the per-event analysis
    void analyze(const Event& event) {
		
		const double weight = event.weight();
		const FastJets& jets2x = apply<FastJets>(event, "Jets2");
		Jets jets2 = jets2x.jetsByPt(.15*GeV);
		foreach (const Jet& jet, jets2)
		{
			counter2pt5++;
			Particles ps = jet.constituents();
			if (any(ps, PtGtr(5*GeV)))
			{
				_histraapt5c1->fill(jet.pT(), weight);
				_histraapt5c3->fill(jet.pT(), weight);
			}
		}

    }


    /// Normalise histograms etc., after the run
    void finalize() {
		scale(_histraapt5c1,  1/counter2pt5);
		scale(_histraapt5c3,   1/counter2pt5);
    }
	
	double counter2pt5 = 0;
	Histo1DPtr _histraapt5c1;
	Histo1DPtr _histraapt5c3;




  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2015_I1343112pp);


}
