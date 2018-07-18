// -*- C++ -*-
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "YODA/ReaderYODA.h"
#include "Rivet/ReferenceDataLoader.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/Projections/ALICEToolsHI.hh"
#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include <fstream>

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2015_I1343112 : public HeavyIonAnalysis {
  public:

    /// Constructor
    ALICE_2015_I1343112(): HeavyIonAnalysis ("ALICE_2015_I1343112"){}

    void init() 
	{
		HeavyIonAnalysis::init();
		// Initialise and register projections
		addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "ImpactParameterMethod");
		ReferenceDataLoader rdl("ALICE_2015_I1343112pp_tmp.yoda");
		_histraapt5c1pp = rdl.getPPReferenceHisto1D("/ALICE_2015_I1343112pp/d06-x01-y01");
		_histraapt5c3pp = rdl.getPPReferenceHisto1D("/ALICE_2015_I1343112pp/d07-x01-y01");

		const FinalState fs(Cuts::pT > 0.15*GeV);
		declare (fs, "Fs");
		FastJets fjr2(fs, FastJets::ANTIKT, 0.2);
		fjr2.useJetArea(new fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
		declare(fjr2, "jets");
		ALICEToolsHI athr2(fs, fjr2);
		declare(athr2, "ATHr2");
		// Book histograms
		_histjsr2c1 = bookHisto1D(2, 1, 1);
		_histjsr2c3 = bookHisto1D(3, 1, 1);
		
		_histdummypt0 = bookHisto1D("TMP/dummyhistpt0",  refData(5, 1, 1));
		_histdummypt3 = bookHisto1D("TMP/dummyhistpt3",  refData(5, 1, 2));
		
		_histdummypt7 = bookHisto1D("TMP/dummyhistpt7",  refData(5, 1, 3));
		_histdummypt10 = bookHisto1D("TMP/dummyhistpt10",  refData(5, 1, 4));
		
		_histdummypt5 = bookHisto1D("TMP/dummyhistpt5",  refData(5, 1, 1));
		_histdummypt5v2 = bookHisto1D("TMP/dummyhistpt5v2",  refData(5, 1, 2));
		_histdummypt5v3 = bookHisto1D("TMP/dummyhistpt5v3",  refData(5, 1, 3));
		_histdummypt5v4 = bookHisto1D("TMP/dummyhistpt5v4",  refData(5, 1, 4));
		
		_histrjsr2pt0by5c1 = bookScatter2D(5, 1, 1);
		_histrjsr2pt3by5c1 = bookScatter2D(5, 1, 2);
		_histrjsr2pt7by5c1 = bookScatter2D(5, 1, 3);
		_histrjsr2pt10by5c1 = bookScatter2D(5, 1, 4);
		
		_histraapt5c1 = bookScatter2D(6, 1, 1);
		_histraapt5c3 = bookScatter2D(7, 1, 1);
		
		_histdumraapt5c1 = bookHisto1D("TMP/dumhistraa1", refData(6, 1, 1));
		_histdumraapt5c3 = bookHisto1D("TMP/dumhistraa2", refData(7, 1, 1));
		
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) 
	{
		const double c = centrality(event, "ImpactParameterMethod");
		cout << c << "cent\n";
		if ((c < 0.)||(c>80.)){
			vetoEvent;
		}
		const double weight = event.weight();
		
		const FastJets& fastjetr2 = apply<FastJets>(event, "jets");
		const Jets& jets2 = fastjetr2.jetsByPt(.15*GeV);
		const auto seqr2 = fastjetr2.clusterSeqArea();
		const ALICEToolsHI &athr2 = apply<ALICEToolsHI>(event, "ATHr2");
		
		if (c >= 0 && c <= 10){
			cencounter10++;
			foreach(const Jet& jets, jets2){
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				const Particles ps = jets.constituents();
				if(any(ps, PtGtr(0*GeV))){
					_histdummypt0->fill(jetpt, weight);
				}
				if(any(ps, PtGtr(3*GeV))){
					_histdummypt3->fill(jetpt, weight);
				}
				if(any(ps, PtGtr(5*GeV))){
					_histjsr2c1->fill(jetpt, weight);
					_histdummypt5->fill(jetpt, weight);
					_histdummypt5v2->fill(jetpt, weight);
					_histdummypt5v3->fill(jetpt, weight);
					_histdummypt5v4->fill(jetpt, weight);
					_histdumraapt5c1->fill(jetpt, weight);
				}
				if(any(ps, PtGtr(7*GeV))){
					_histdummypt7->fill(jetpt, weight);
				}
				if(any(ps, PtGtr(10*GeV))){
					_histdummypt10->fill(jetpt, weight);
				}
			}
		}
		if (c >= 10 && c <= 30){
			cencounter30++;
			foreach(const Jet& jets, jets2){
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				const Particles ps = jets.constituents();
				if(any(ps, PtGtr(5*GeV))){
					_histjsr2c3->fill(jetpt, weight);
					_histdumraapt5c3->fill(jetpt, weight);
				}
			}
		}


    }


    /// Normalise histograms etc., after the run
    void finalize() 
	{
		scale(_histjsr2c1, (1e-8)/cencounter10);
		scale(_histjsr2c3, (1e-8)/cencounter30);
		scale(_histdumraapt5c1, (1e-8)/cencounter10);
		scale(_histdumraapt5c3, (1e-8)/cencounter30);
		divide(_histdummypt0, _histdummypt5, _histrjsr2pt0by5c1);
		divide(_histdummypt3, _histdummypt5v2, _histrjsr2pt3by5c1);
		divide(_histdummypt7, _histdummypt5v3, _histrjsr2pt7by5c1);
		divide(_histdummypt10, _histdummypt5v4, _histrjsr2pt10by5c1);
		divide(_histdumraapt5c1, _histraapt5c1pp, _histraapt5c1);
		divide(_histdumraapt5c3, _histraapt5c3pp, _histraapt5c3);
    }


	Histo1DPtr _histjsr2c1;
	Histo1DPtr _histjsr2c3;
	Histo1DPtr _histdummypt0;
	Histo1DPtr _histdummypt3;
	Histo1DPtr _histdummypt7;
	Histo1DPtr _histdummypt10;
	
	Histo1DPtr _histdummypt5;
	Histo1DPtr _histdummypt5v2;
	Histo1DPtr _histdummypt5v3;
	Histo1DPtr _histdummypt5v4;
	
	Scatter2DPtr _histrjsr2pt0by5c1;
	Scatter2DPtr _histrjsr2pt3by5c1;
	Scatter2DPtr _histrjsr2pt7by5c1;
	Scatter2DPtr _histrjsr2pt10by5c1;
	
	Histo1DPtr _histraapt5c1pp;
	Histo1DPtr _histraapt5c3pp;
	
	Histo1DPtr _histdumraapt5c1;
	Histo1DPtr _histdumraapt5c3;
	
	Scatter2DPtr _histraapt5c1;
	Scatter2DPtr _histraapt5c3;
	
	double cencounter10 = 0;
	double cencounter30 = 0;
	
	void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr out){
		const string path = out->path();
		*out = *h1 / *h2;
		out->setPath(path);
	}

  };
	

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2015_I1343112);


}
