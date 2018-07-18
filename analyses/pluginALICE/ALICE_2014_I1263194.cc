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
  class ALICE_2014_I1263194 : public HeavyIonAnalysis {
  public:
    /// Constructor
    ALICE_2014_I1263194(): HeavyIonAnalysis ("ALICE_2014_I1263194"){}
    void init() 
	{
	  HeavyIonAnalysis::init();
      // Initialise and register projections
	  // 0-10% centrality 
	  addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "ImpactParameterMethod");

	  //setting up the cuts 
	  //{
	  const ChargedFinalState cfs(Cuts::pT > 0.15*GeV);
	  declare (cfs, "cfs");
	  FastJets fjr2(cfs, FastJets::ANTIKT, 0.2);
	  FastJets fjr3(cfs, FastJets::ANTIKT, 0.3);
	  fjr2.useJetArea(new fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
	  fjr3.useJetArea(new fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
	  declare(fjr2, "Jetsr2");
	  declare(fjr3, "Jetsr3");
	  ALICEToolsHI athr2(cfs, fjr2);
	  ALICEToolsHI athr3(cfs, fjr3);
	  declare(athr2, "ATHr2");
	  declare(athr3, "ATHr3");
	  //}
      // Book histograms
	  //{
	  // booking histograms with pT > 1.5GeV at different centralities
      _h_CJS_R2_P15_C10 = bookHisto1D(2, 1, 1);
	  _h_CJS_R3_P15_C10 = bookHisto1D(2, 1, 2);
	  
      _h_CJS_R2_P15_C30 = bookHisto1D(3, 1, 1);
	  _h_CJS_R3_P15_C30 = bookHisto1D(3, 1, 2);
	  
	  _h_CJS_R2_P15_C50 = bookHisto1D(4, 1, 1);
	  _h_CJS_R3_P15_C50 = bookHisto1D(4, 1, 2);
	  
	  _h_CJS_R2_P15_C80 = bookHisto1D(5, 1, 1);
	  _h_CJS_R3_P15_C80 = bookHisto1D(5, 1, 2);
	  
	  // booking histograms with pT > 5GeV at different centralities
	  _h_CJS_R2_P5_C10 = bookHisto1D(6, 1, 1);
	  _h_CJS_R3_P5_C10 = bookHisto1D(6, 1, 2);
	  
      _h_CJS_R2_P5_C30 = bookHisto1D(7, 1, 1);
	  _h_CJS_R3_P5_C30 = bookHisto1D(7, 1, 2);
	  
	  _h_CJS_R2_P5_C50 = bookHisto1D(8, 1, 1);
	  _h_CJS_R3_P5_C50 = bookHisto1D(8, 1, 2);
	  
	  _h_CJS_R2_P5_C80 = bookHisto1D(9, 1, 1);
	  _h_CJS_R3_P5_C80 = bookHisto1D(9, 1, 2);
	  
	  // booking histograms with pT > 10GeV at different centralities
	  _h_CJS_R2_P10_C10 = bookHisto1D("TMP/CJSR2P10C10", refData(10, 1, 1));
	  _h_CJS_R3_P10_C10 = bookHisto1D("TMP/CJSR3P10C10", refData(10, 1, 2));
      
	  _h_CJS_R2_P10_C30 = bookHisto1D("TMP/CJSR2P10C30", refData(11, 1, 1));
	  _h_CJS_R3_P10_C30 = bookHisto1D("TMP/CJSR3P20C30", refData(11, 1, 2));
	  
	  _h_CJS_R2_P10_C50 = bookHisto1D("TMP/CJSR2P10C50", refData(12, 1, 1));
	  _h_CJS_R3_P10_C50 = bookHisto1D("TMP/CJSR3P10C50", refData(12, 1, 2));
	  
	  _h_CJS_R2_P10_C80 = bookHisto1D(13, 1, 1);
	  _h_CJS_R3_P10_C80 = bookHisto1D(13, 1, 2);
	  
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  _histRJS15by5r2c1 = bookScatter2D(14, 1, 1);
	  _histRJS15by5r2c8  = bookScatter2D(15, 1, 1);
	  _histRJS10by5r2c1 = bookScatter2D(16, 1, 1);
	  _histRJS10by5r2c8 = bookScatter2D(17, 1, 1);
	  
	  _histRJS15by5r3c1 = bookScatter2D(14, 1, 2);
	  _histRJS15by5r3c8 = bookScatter2D(15, 1, 2);
	  _histRJS10by5r3c1 = bookScatter2D(16, 1, 2);
	  _histRJS10by5r3c8 = bookScatter2D(17, 1, 2);
	  
	  _histdummy15r2c1 = bookHisto1D("TMP/dummyhist",  refData(14, 1, 1));
	  _histdummy15r3c1 = bookHisto1D("TMP/dummyhist2", refData(14, 1, 2));
	  _histdummy15r2c8 = bookHisto1D("TMP/dummyhist3", refData(15, 1, 1));
	  _histdummy15r3c8 = bookHisto1D("TMP/dummyhist4", refData(15, 1, 2));
	  
	  _histdummy10r2c1 = bookHisto1D("TMP/dummyhist5",  refData(16, 1, 1));
	  _histdummy10r3c1 = bookHisto1D("TMP/dummyhist6", refData(16, 1, 2));
	  _histdummy10r2c8 = bookHisto1D("TMP/dummyhist7", refData(17, 1, 1));
	  _histdummy10r3c8 = bookHisto1D("TMP/dummyhist8", refData(17, 1, 2));
	  
	  _histdummy5r2c1  = bookHisto1D("TMP/dummyhist9", refData(14, 1, 1));
	  _histdummy5r3c1  = bookHisto1D("TMP/dummyhist10", refData(14, 1, 2));
	  _histdummy5r2c8  = bookHisto1D("TMP/dummyhist11", refData(15, 1, 1));
	  _histdummy5r3c8  = bookHisto1D("TMP/dummyhist12", refData(15, 1, 2));
	  
	  _histdummy5r2c1v2 = bookHisto1D("TMP/dummyhist13", refData(16, 1, 1));
	  _histdummy5r3c1v2 = bookHisto1D("TMP/dummyhist14", refData(16, 1, 2));
	  _histdummy5r2c8v2 = bookHisto1D("TMP/dummyhist15", refData(17, 1, 1));
	  _histdummy5r3c8v2 = bookHisto1D("TMP/dummyhist16", refData(17, 1, 2));
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  _histRCPr2pt5c1per8 = bookScatter2D(21, 1, 1);
	  _histRCPr2pt5c3per8 = bookScatter2D(22, 1, 1);
	  _histRCPr2pt5c5per8 = bookScatter2D(23, 1, 1);
	
	  _histRCPr3pt5c1per8 = bookScatter2D(21, 1, 2);
	  _histRCPr3pt5c3per8 = bookScatter2D(22, 1, 2);
	  _histRCPr3pt5c5per8 = bookScatter2D(23, 1, 2);
	
	  _histdummyRCPr2pt5c1 = bookHisto1D("TMP/dummyhistRCPr2pt5c1", refData(21, 1, 1));
	  _histdummyRCPr2pt5c3 = bookHisto1D("TMP/dummyhistRCPr2pt5c3", refData(22, 1, 1));
	  _histdummyRCPr2pt5c5 = bookHisto1D("TMP/dummyhistRCPr2pt5c5", refData(23, 1, 1));
	
	  _histdummyRCPr3pt5c1 = bookHisto1D("TMP/dummyhistRCPr3pt5c1", refData(21, 1, 2));
	  _histdummyRCPr3pt5c3 = bookHisto1D("TMP/dummyhistRCPr3pt5c3", refData(22, 1, 2));
	  _histdummyRCPr3pt5c5 = bookHisto1D("TMP/dummyhistRCPr3pt5c5", refData(23, 1, 2));
	
	  _histdummyRCPr2pt5c8 = bookHisto1D("TMP/dummyhistRCPr2c8", refData(21, 1, 1));
	  _histdummyRCPr2pt5c8v2 = bookHisto1D("TMP/dummyhistRCPr28v2", refData(22, 1, 1));
	  _histdummyRCPr2pt5c8v3 = bookHisto1D("TMP/dummyhistRCPr28v3", refData(23, 1, 1));
	  
	  _histdummyRCPr3pt5c8 = bookHisto1D("TMP/dummyhistRCPr3c8", refData(21, 1, 2));
	  _histdummyRCPr3pt5c8v2 = bookHisto1D("TMP/dummyhistRCPr3c8v2", refData(22, 1, 2));
	  _histdummyRCPr3pt5c8v3 = bookHisto1D("TMP/dummyhistRCPr3c8v3", refData(23, 1, 2));
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  _histRCPvNpr2pt = bookScatter2D(27, 1, 1);
	  _histRCPvNpr2pt5 = bookScatter2D(28, 1, 1);
	  _histRCPvNpr2pt10 = bookScatter2D(29, 1, 1);
	  
	  _histRCPvNpr3pt = bookScatter2D(27, 1, 2);
	  _histRCPvNpr3pt5 = bookScatter2D(28, 1, 2);
	  _histRCPvNpr3pt10 = bookScatter2D(29, 1, 2);
	  //////////////////////////////////////////////////////////////////////////////////////////
	  _histRCPr2ptc1per8 = bookScatter2D(18, 1, 1);
	  _histRCPr2ptc3per8 = bookScatter2D(19, 1, 1);
	  _histRCPr2ptc5per8 = bookScatter2D(20, 1, 1);
	
	  _histRCPr3ptc1per8 = bookScatter2D(18, 1, 2);
	  _histRCPr3ptc3per8 = bookScatter2D(19, 1, 2);
	  _histRCPr3ptc5per8 = bookScatter2D(20, 1, 2);
	  
	  _histdummyRCPr2ptc1 = bookHisto1D("TMP/dummyhistRCPr2ptc1", refData(18, 1, 1));
	  _histdummyRCPr2ptc3 = bookHisto1D("TMP/dummyhistRCPr2ptc3", refData(19, 1, 1));
	  _histdummyRCPr2ptc5 = bookHisto1D("TMP/dummyhistRCPr2ptc5", refData(20, 1, 1));
	
	  _histdummyRCPr3ptc1 = bookHisto1D("TMP/dummyhistRCPr3ptc1", refData(18, 1, 2));
	  _histdummyRCPr3ptc3 = bookHisto1D("TMP/dummyhistRCPr3ptc3", refData(19, 1, 2));
	  _histdummyRCPr3ptc5 = bookHisto1D("TMP/dummyhistRCPr3ptc5", refData(20, 1, 2));
	  
	  _histdummyRCPr2ptc8 = bookHisto1D("TMP/dummyhistRCPr2c8pt", refData(18, 1, 1));
	  _histdummyRCPr2ptc8v2 = bookHisto1D("TMP/dummyhistRCPr28v2pt", refData(19, 1, 1));
	  _histdummyRCPr2ptc8v3 = bookHisto1D("TMP/dummyhistRCPr28v3pt", refData(20, 1, 1));
	  
	  _histdummyRCPr3ptc8 = bookHisto1D("TMP/dummyhistRCPr3c8pt", refData(18, 1, 2));
	  _histdummyRCPr3ptc8v2 = bookHisto1D("TMP/dummyhistRCPr3c8v2pt", refData(19, 1, 2));
	  _histdummyRCPr3ptc8v3 = bookHisto1D("TMP/dummyhistRCPr3c8v3pt", refData(20, 1, 2));
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  _histRCPr2pt10c1per8 = bookScatter2D(24, 1, 1);
	  _histRCPr2pt10c3per8 = bookScatter2D(25, 1, 1);
	  _histRCPr2pt10c5per8 = bookScatter2D(26, 1, 1);
	
	  _histRCPr3pt10c1per8 = bookScatter2D(24, 1, 2);
	  _histRCPr3pt10c3per8 = bookScatter2D(25, 1, 2);
	  _histRCPr3pt10c5per8 = bookScatter2D(26, 1, 2);
	  
	  _histdummyRCPr2pt10c1 = bookHisto1D("TMP/dummyhistRCPr2pt10c1", refData(24, 1, 1));
	  _histdummyRCPr2pt10c3 = bookHisto1D("TMP/dummyhistRCPr2pt10c3", refData(25, 1, 1));
	  _histdummyRCPr2pt10c5 = bookHisto1D("TMP/dummyhistRCPr2pt10c5", refData(26, 1, 1));
	
	  _histdummyRCPr3pt10c1 = bookHisto1D("TMP/dummyhistRCPr3pt10c1", refData(24, 1, 2));
	  _histdummyRCPr3pt10c3 = bookHisto1D("TMP/dummyhistRCPr3pt10c3", refData(25, 1, 2));
	  _histdummyRCPr3pt10c5 = bookHisto1D("TMP/dummyhistRCPr3pt10c5", refData(26, 1, 2));
	  
	  _histdummyRCPr2pt10c8 = bookHisto1D("TMP/dummyhistRCPr2c8pt10", refData(24, 1, 1));
	  _histdummyRCPr2pt10c8v2 = bookHisto1D("TMP/dummyhistRCPr28v2pt10", refData(25, 1, 1));
	  _histdummyRCPr2pt10c8v3 = bookHisto1D("TMP/dummyhistRCPr28v3pt10", refData(26, 1, 1));
	  
	  _histdummyRCPr3pt10c8 = bookHisto1D("TMP/dummyhistRCPr3c8pt10", refData(24, 1, 2));
	  _histdummyRCPr3pt10c8v2 = bookHisto1D("TMP/dummyhistRCPr3c8v2pt10", refData(25, 1, 2));
	  _histdummyRCPr3pt10c8v3 = bookHisto1D("TMP/dummyhistRCPr3c8v3pt10", refData(26, 1, 2));
	  
	  _histratiocjsr2c1 = bookHisto1D("TMP/dumratiocjs1", refData(30, 1, 1));
	  _histratiocjsr3c1 = bookHisto1D("TMP/dumratiocjs2", refData(30, 1, 1));
	  _histratiocjsr2c8 = bookHisto1D("TMP/dumratiocjs3", refData(30, 1, 2));
	  _histratiocjsr3c8 = bookHisto1D("TMP/dumratiocjs4", refData(30, 1, 2));
	  
	  _histratiocjsr2by3c1 = bookScatter2D(30, 1, 1);
	  _histratiocjsr2by3c8 = bookScatter2D(30, 1, 2);
	  //}
    }



    void analyze(const Event& event) 
	{

		const double c = centrality(event, "ImpactParameterMethod");
		cout << c << "cent\n";
		if ((c < 0.)||(c>80.)){
			vetoEvent;
		}
		const double weight = event.weight();
		//decalring things 
		//{
		const FastJets& fastjetr2 = apply<FastJets>(event, "Jetsr2");
		const FastJets& fastjetr3 = apply<FastJets>(event, "Jetsr3");
		// setting up the jets 
		const Jets& jets2 = fastjetr2.jetsByPt(.15*GeV);
		const Jets& jets3 = fastjetr3.jetsByPt(.15*GeV);
		
		const auto seqr2 = fastjetr2.clusterSeqArea();
		const auto seqr3 = fastjetr3.clusterSeqArea();
		const ALICEToolsHI &athr2 = apply<ALICEToolsHI>(event, "ATHr2");
		const ALICEToolsHI &athr3 = apply<ALICEToolsHI>(event, "ATHr3");
		//}
		// analysis and fill for histogram 
		// 0-10% centrality 

		if (c >= 0 && c <= 10) {
			cencounter10++;
			_Ncollc1 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Ncoll(): 0;
			_Nparttc1 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Npart_targ(): 0;
			_Npartpc1 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Npart_proj(): 0;
			foreach(const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				
				_h_CJS_R2_P15_C10->fill(jetpt, weight);
				_histdummy15r2c1->fill(jetpt, weight);
				_histdummyRCPr2ptc1->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R2_P5_C10->fill(jetpt, weight);
					_histdummy5r2c1->fill(jetpt, weight);
					_histdummy5r2c1v2->fill(jetpt, weight);
					_histdummyRCPr2pt5c1->fill(jetpt, weight);
					_histratiocjsr2c1->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R2_P10_C10->fill(jetpt, weight);
					_histdummy10r2c1->fill(jetpt, weight);
					_histdummyRCPr2pt10c1->fill(jetpt, weight);
				}
			}	
			foreach(const Jet& jets, jets3)
			{
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R3_P15_C10->fill(jetpt, weight);
				_histdummy15r3c1->fill(jetpt, weight);
				_histdummyRCPr3ptc1->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R3_P5_C10->fill(jetpt, weight);
					_histdummy5r3c1->fill(jetpt, weight);
					_histdummy5r3c1v2->fill(jetpt, weight);
					_histdummyRCPr3pt5c1->fill(jetpt, weight);
					_histratiocjsr3c1->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R3_P10_C10->fill(jetpt, weight);
					_histdummy10r3c1->fill(jetpt, weight);
					_histdummyRCPr3pt10c1->fill(jetpt, weight);
				}	
			}	
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 10-30% centrality 
		if (10 < c && c <= 30){
			cencounter30++;
			_Ncollc3 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Ncoll(): 0;
			_Nparttc3 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Npart_targ(): 0;
			_Npartpc3 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Npart_proj(): 0;
			foreach (const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R2_P15_C30->fill(jetpt,weight);
				_histdummyRCPr2ptc3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R2_P5_C30->fill(jetpt, weight);
					_histdummyRCPr2pt5c3->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R2_P10_C30->fill(jetpt, weight);
					_histdummyRCPr2pt10c3->fill(jetpt, weight);
				}
			}
			foreach (const Jet& jets, jets3)
			{
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R3_P15_C30->fill(jetpt,weight);
				_histdummyRCPr3ptc3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R3_P5_C30->fill(jetpt, weight);
					_histdummyRCPr3pt5c3->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R3_P10_C30->fill(jetpt, weight);
					_histdummyRCPr3pt10c3->fill(jetpt, weight);
				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 30-50% centrality 
		if ( 30 < c && c <= 50){
			cencounter50++;
			_Ncollc5 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Ncoll(): 0;
			_Nparttc5 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Npart_targ(): 0;
			_Npartpc5 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Npart_proj(): 0;
			foreach (const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R2_P15_C50->fill(jetpt,weight);
				_histdummyRCPr2ptc5->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R2_P5_C50->fill(jetpt, weight);
					_histdummyRCPr2pt5c5->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R2_P10_C50->fill(jetpt, weight);
					_histdummyRCPr2pt10c5->fill(jetpt, weight);
				}
			}
			foreach (const Jet& jets, jets3)
			{
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R3_P15_C50->fill(jetpt,weight);
				_histdummyRCPr3ptc5->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R3_P5_C50->fill(jetpt, weight);
					_histdummyRCPr3pt5c5->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R3_P10_C50->fill(jetpt, weight);
					_histdummyRCPr3pt10c5->fill(jetpt, weight);
				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// 50-80% centrality 
		if ( 50 < c && c <= 80){
			cencounter80++;
			_Ncollc8 += event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->Ncoll(): 0;
			foreach (const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R2_P15_C80->fill(jetpt,weight);
				_histdummy15r2c8->fill(jetpt, weight);
				_histdummyRCPr2ptc8->fill(jetpt, weight);
				_histdummyRCPr2ptc8v2->fill(jetpt, weight);
				_histdummyRCPr2ptc8v3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R2_P5_C80->fill(jetpt, weight);
					_histdummy5r2c8->fill(jetpt, weight);
					_histdummy5r2c8v2->fill(jetpt, weight);
					_histdummyRCPr2pt5c8->fill(jetpt, weight);
					_histdummyRCPr2pt5c8v2->fill(jetpt, weight);
					_histdummyRCPr2pt5c8v3->fill(jetpt, weight);
					_histratiocjsr2c8->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R2_P10_C80->fill(jetpt, weight);
					_histdummy10r2c8->fill(jetpt, weight);
					_histdummyRCPr2pt10c8->fill(jetpt, weight);
					_histdummyRCPr2pt10c8v2->fill(jetpt, weight);
					_histdummyRCPr2pt10c8v3->fill(jetpt, weight);
				}
			}
			foreach (const Jet& jets, jets3)
			{
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_h_CJS_R3_P15_C80->fill(jetpt,weight);
				_histdummy15r3c8->fill(jetpt, weight);
				_histdummyRCPr3ptc8->fill(jetpt, weight);
				_histdummyRCPr3ptc8v2->fill(jetpt, weight);
				_histdummyRCPr3ptc8v3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_h_CJS_R3_P5_C80->fill(jetpt, weight);
					_histdummy5r3c8->fill(jetpt, weight);
					_histdummy5r2c8v2->fill(jetpt, weight);
					_histdummyRCPr3pt5c8->fill(jetpt, weight);
					_histdummyRCPr3pt5c8v2->fill(jetpt, weight);
					_histdummyRCPr3pt5c8v3->fill(jetpt, weight);
					_histratiocjsr3c8->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_h_CJS_R3_P10_C80->fill(jetpt, weight);
					_histdummy10r3c8->fill(jetpt, weight);
					_histdummyRCPr3pt10c8->fill(jetpt, weight);
					_histdummyRCPr3pt10c8v2->fill(jetpt, weight);
					_histdummyRCPr3pt10c8v3->fill(jetpt, weight);

				}
			}
		}
    }
    /// Normalise histograms etc., after the run
    void finalize() {
		//scaling the historgram 
		// 0-10% centrality 
		//{
		double _Npartc1 = (_Nparttc1 + _Npartpc1)/cencounter10;
		double _Npartc3 = (_Nparttc3 + _Npartpc3)/cencounter30;
		double _Npartc5 = (_Nparttc5 + _Npartpc5)/cencounter50;
		cout << _Npartc1 << endl;
		cout << _Npartc3 << endl;
		cout << _Npartc5 << endl;
		if (_Ncollc1 > 0.) cencounter10*=_Ncollc1;
		if (_Ncollc3 > 0.) cencounter30*=_Ncollc3;
		if (_Ncollc5 > 0.) cencounter50*=_Ncollc5;
		if (_Ncollc8 > 0.) cencounter80*=_Ncollc8;
		scale(_h_CJS_R2_P15_C10, 1e-6/cencounter10);
		scale(_h_CJS_R3_P15_C10, 1e-6/cencounter10);
		
		scale(_h_CJS_R2_P5_C10,  1e-6/cencounter10);
		scale(_h_CJS_R3_P5_C10,  1e-6/cencounter10);
		
		scale(_h_CJS_R2_P10_C10,  1e-6/cencounter10);
		scale(_h_CJS_R3_P10_C10,  1e-6/cencounter10);
		
		scale(_histdummyRCPr2pt5c1,  64/cencounter10);
		scale(_histdummyRCPr3pt5c1,  64/cencounter10);
		
		scale(_histdummyRCPr2ptc1,  64/cencounter10);
		scale(_histdummyRCPr3ptc1,  64/cencounter10);
		
		scale(_histdummyRCPr2pt10c1,  64/cencounter10);
		scale(_histdummyRCPr3pt10c1,  64/cencounter10);
		
		scale(_histratiocjsr2c1, 1/cencounter10);
		scale(_histratiocjsr3c1, 1/cencounter10);
		//}
		//}
		
		// 10-30% centrality 
		//{
		scale(_h_CJS_R2_P15_C30, 1e-6/cencounter30);
		scale(_h_CJS_R3_P15_C30, 1e-6/cencounter30);
		
		scale(_h_CJS_R2_P5_C30,  1e-6/cencounter30);
		scale(_h_CJS_R3_P5_C30,  1e-6/cencounter30);
		
		scale(_h_CJS_R2_P10_C30,  1e-6/cencounter30);
		scale(_h_CJS_R3_P10_C30,  1e-6/cencounter30);
		
		scale(_histdummyRCPr2pt5c3,  64/cencounter30);
		scale(_histdummyRCPr3pt5c3,  64/cencounter30);
		
		scale(_histdummyRCPr2ptc3,  64/cencounter30);
		scale(_histdummyRCPr3ptc3,  64/cencounter30);
		
		scale(_histdummyRCPr2pt10c3,  64/cencounter30);
		scale(_histdummyRCPr3pt10c3,  64/cencounter30);
		//}
		
		// 30-50% centrality 
		//{
		scale(_h_CJS_R2_P15_C50, 1e-6/cencounter50);
		scale(_h_CJS_R3_P15_C50, 1e-6/cencounter50);
		
		scale(_h_CJS_R2_P5_C50,  1e-6/cencounter50);
		scale(_h_CJS_R3_P5_C50,  1e-6/cencounter50);
		
		scale(_h_CJS_R2_P10_C50,  1e-6/cencounter50);
		scale(_h_CJS_R3_P10_C50,  1e-6/cencounter50);
		
		scale(_histdummyRCPr2pt5c5,  64/cencounter50);
		scale(_histdummyRCPr3pt5c5,  64/cencounter50);
		
		scale(_histdummyRCPr2ptc5,  64/cencounter50);
		scale(_histdummyRCPr3ptc5,  64/cencounter50);
		
		scale(_histdummyRCPr2pt10c5,  64/cencounter50);
		scale(_histdummyRCPr3pt10c5,  64/cencounter50);
		//}
		
		// 50-80% centrality 
		//{
		scale(_h_CJS_R2_P15_C80, 1e-6/cencounter80);
		scale(_h_CJS_R3_P15_C80, 1e-6/cencounter80);
		
		scale(_h_CJS_R2_P5_C80,  1e-6/cencounter80);
		scale(_h_CJS_R3_P5_C80,  1e-6/cencounter80);
		
		scale(_h_CJS_R2_P10_C80,  1e-6/cencounter80);
		scale(_h_CJS_R3_P10_C80,  1e-6/cencounter80);
		
		scale(_histdummyRCPr2pt5c8,  64/cencounter80);
		scale(_histdummyRCPr2pt5c8v2,  64/cencounter80);
		scale(_histdummyRCPr2pt5c8v3,  64/cencounter80);
		
		scale(_histdummyRCPr3pt5c8,  64/cencounter80);
		scale(_histdummyRCPr3pt5c8v2,  64/cencounter80);
		scale(_histdummyRCPr3pt5c8v3,  64/cencounter80);
		
		scale(_histdummyRCPr2ptc8,  64/cencounter80);
		scale(_histdummyRCPr2ptc8v2,  64/cencounter80);
		scale(_histdummyRCPr2ptc8v3,  64/cencounter80);
		
		scale(_histdummyRCPr3ptc8,  64/cencounter80);
		scale(_histdummyRCPr3ptc8v2,  64/cencounter80);
		scale(_histdummyRCPr3ptc8v3,  64/cencounter80);
		
		scale(_histdummyRCPr2pt10c8,  64/cencounter80);
		scale(_histdummyRCPr2pt10c8v2,  64/cencounter80);
		scale(_histdummyRCPr2pt10c8v3,  64/cencounter80);
		
		scale(_histdummyRCPr3pt10c8,  64/cencounter80);
		scale(_histdummyRCPr3pt10c8v2,  64/cencounter80);
		scale(_histdummyRCPr3pt10c8v3,  64/cencounter80);
		
		scale(_histratiocjsr2c1, 1/cencounter80);
		scale(_histratiocjsr3c1, 1/cencounter80);
		//}
		//dividing histograms
		//{
		
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummy15r2c1, _histdummy5r2c1, _histRJS15by5r2c1);
		divide(_histdummy15r2c8, _histdummy5r2c8, _histRJS15by5r2c8);
		divide(_histdummy10r2c1, _histdummy5r2c1v2, _histRJS10by5r2c1);
		divide(_histdummy10r2c8, _histdummy5r2c8v2, _histRJS10by5r2c8);
		
		divide(_histdummy15r3c1, _histdummy5r3c1, _histRJS15by5r3c1);
		divide(_histdummy15r3c8, _histdummy5r3c8, _histRJS15by5r3c8);
		divide(_histdummy10r3c1, _histdummy5r3c1v2, _histRJS10by5r3c1);
		divide(_histdummy10r3c8, _histdummy5r3c8v2, _histRJS10by5r3c8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummyRCPr2pt5c1, _histdummyRCPr2pt5c8, _histRCPr2pt5c1per8);
		divide(_histdummyRCPr2pt5c3, _histdummyRCPr2pt5c8v2, _histRCPr2pt5c3per8);
		divide(_histdummyRCPr2pt5c5, _histdummyRCPr2pt5c8v3, _histRCPr2pt5c5per8);
		
		divide(_histdummyRCPr3pt5c1, _histdummyRCPr3pt5c8, _histRCPr3pt5c1per8);
		divide(_histdummyRCPr3pt5c3, _histdummyRCPr3pt5c8v2, _histRCPr3pt5c3per8);
		divide(_histdummyRCPr3pt5c5, _histdummyRCPr3pt5c8v3, _histRCPr3pt5c5per8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummyRCPr2ptc1, _histdummyRCPr2ptc8, _histRCPr2ptc1per8);
		divide(_histdummyRCPr2ptc3, _histdummyRCPr2ptc8v2, _histRCPr2ptc3per8);
		divide(_histdummyRCPr2ptc5, _histdummyRCPr2ptc8v3, _histRCPr2ptc5per8);
		
		divide(_histdummyRCPr3ptc1, _histdummyRCPr3ptc8, _histRCPr3ptc1per8);
		divide(_histdummyRCPr3ptc3, _histdummyRCPr3ptc8v2, _histRCPr3ptc3per8);
		divide(_histdummyRCPr3ptc5, _histdummyRCPr3ptc8v3, _histRCPr3ptc5per8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummyRCPr2pt10c1, _histdummyRCPr2pt10c8, _histRCPr2pt10c1per8);
		divide(_histdummyRCPr2pt10c3, _histdummyRCPr2pt10c8v2, _histRCPr2pt10c3per8);
		divide(_histdummyRCPr2pt10c5, _histdummyRCPr2pt10c8v3, _histRCPr2pt10c5per8);
		
		divide(_histdummyRCPr3pt10c1, _histdummyRCPr3pt10c8, _histRCPr3pt10c1per8);
		divide(_histdummyRCPr3pt10c3, _histdummyRCPr3pt10c8v2, _histRCPr3pt10c3per8);
		divide(_histdummyRCPr3pt10c5, _histdummyRCPr3pt10c8v3, _histRCPr3pt10c5per8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histratiocjsr2c1, _histratiocjsr3c1, _histratiocjsr2by3c1);
		divide(_histratiocjsr2c8, _histratiocjsr3c8, _histratiocjsr2by3c8);
////////////////////////////////////////////////////////////////////////////////////////////////
		//}
		//adding points to the RCP plots
		//{
////////////////////////////////////////////////////////////////////////////////////////////////
		_histRCPvNpr2pt->addPoint(_Npartc1, _histRCPr2ptc1per8->point(4).y());
		_histRCPvNpr2pt->addPoint(_Npartc3, _histRCPr2ptc3per8->point(4).y());
		_histRCPvNpr2pt->addPoint(_Npartc5, _histRCPr2ptc5per8->point(4).y());

		_histRCPvNpr3pt->addPoint(_Npartc1, _histRCPr3ptc1per8->point(3).y());
		_histRCPvNpr3pt->addPoint(_Npartc3, _histRCPr3ptc3per8->point(3).y());
		_histRCPvNpr3pt->addPoint(_Npartc5, _histRCPr3ptc5per8->point(3).y());
////////////////////////////////////////////////////////////////////////////////////////////////		
		_histRCPvNpr2pt5->addPoint(_Npartc1, _histRCPr2pt5c1per8->point(4).y());
		_histRCPvNpr2pt5->addPoint(_Npartc3, _histRCPr2pt5c3per8->point(4).y());
		_histRCPvNpr2pt5->addPoint(_Npartc5, _histRCPr2pt5c5per8->point(4).y());
		
		_histRCPvNpr3pt5->addPoint(_Npartc1, _histRCPr3pt5c1per8->point(3).y());
		_histRCPvNpr3pt5->addPoint(_Npartc3, _histRCPr3pt5c3per8->point(3).y());
		_histRCPvNpr3pt5->addPoint(_Npartc5, _histRCPr3pt5c5per8->point(3).y());
////////////////////////////////////////////////////////////////////////////////////////////////
		_histRCPvNpr2pt10->addPoint(_Npartc1, _histRCPr2pt10c1per8->point(4).y());
		_histRCPvNpr2pt10->addPoint(_Npartc3, _histRCPr2pt10c3per8->point(4).y());
		_histRCPvNpr2pt10->addPoint(_Npartc5, _histRCPr2pt10c5per8->point(4).y());
		
		_histRCPvNpr3pt10->addPoint(_Npartc1, _histRCPr3pt10c1per8->point(4).y());
		_histRCPvNpr3pt10->addPoint(_Npartc3, _histRCPr3pt10c3per8->point(4).y());
		_histRCPvNpr3pt10->addPoint(_Npartc5, _histRCPr3pt10c5per8->point(4).y());
		//}
    }
	// counters for each measurement for scalling reasons
	//{
	double cencounter10 = 0;
	double cencounter30 = 0;
	double cencounter50 = 0;
	double cencounter80 = 0;
    /// Perform the per-event analysis
//}
	//declareing histograms for this file
	//{
	//pT > 1.5 GeV/c at different centrlities 
	Histo1DPtr _h_CJS_R2_P15_C10;
	Histo1DPtr _h_CJS_R3_P15_C10;
	
	Histo1DPtr _h_CJS_R2_P15_C30;
	Histo1DPtr _h_CJS_R3_P15_C30;
	
	Histo1DPtr _h_CJS_R2_P15_C50;
	Histo1DPtr _h_CJS_R3_P15_C50;
	
	Histo1DPtr _h_CJS_R2_P15_C80;
	Histo1DPtr _h_CJS_R3_P15_C80;
	// histograms with pT > 5 GeV/c at different centrlities 
	Histo1DPtr _h_CJS_R2_P5_C10;
	Histo1DPtr _h_CJS_R3_P5_C10;
	
	Histo1DPtr _h_CJS_R2_P5_C30;
	Histo1DPtr _h_CJS_R3_P5_C30;
	
	Histo1DPtr _h_CJS_R2_P5_C50;
	Histo1DPtr _h_CJS_R3_P5_C50;
	
	Histo1DPtr _h_CJS_R2_P5_C80;
	Histo1DPtr _h_CJS_R3_P5_C80;
	// histograms with pT > 10 GeV/c at different centrlities 
	Histo1DPtr _h_CJS_R2_P10_C10;
	Histo1DPtr _h_CJS_R3_P10_C10;
	
	Histo1DPtr _h_CJS_R2_P10_C30;
	Histo1DPtr _h_CJS_R3_P10_C30;
	
	Histo1DPtr _h_CJS_R2_P10_C50;
	Histo1DPtr _h_CJS_R3_P10_C50;
	
	Histo1DPtr _h_CJS_R2_P10_C80;
	Histo1DPtr _h_CJS_R3_P10_C80;
///////////////////////////////////////////////////////////////////
	Scatter2DPtr  _histRJS15by5r2c1;
	Scatter2DPtr  _histRJS15by5r2c8;
	Scatter2DPtr  _histRJS10by5r2c1;
	Scatter2DPtr  _histRJS10by5r2c8;

	Scatter2DPtr  _histRJS15by5r3c1;
	Scatter2DPtr  _histRJS15by5r3c8;
	Scatter2DPtr  _histRJS10by5r3c1;
	Scatter2DPtr  _histRJS10by5r3c8;
	
	Histo1DPtr _histdummy15r2c1;
	Histo1DPtr _histdummy15r2c8;
	Histo1DPtr _histdummy15r3c1;
	Histo1DPtr _histdummy15r3c8;
	
	Histo1DPtr _histdummy10r2c1;
	Histo1DPtr _histdummy10r2c8;
	Histo1DPtr _histdummy10r3c1;
	Histo1DPtr _histdummy10r3c8;
	
	Histo1DPtr _histdummy5r2c1;
	Histo1DPtr _histdummy5r2c8;
	Histo1DPtr _histdummy5r3c1;
	Histo1DPtr _histdummy5r3c8;
	
	Histo1DPtr _histdummy5r2c1v2;
	Histo1DPtr _histdummy5r2c8v2;
	Histo1DPtr _histdummy5r3c1v2;
	Histo1DPtr _histdummy5r3c8v2;
///////////////////////////////////////////////////////////////////	
	Scatter2DPtr _histRCPr2pt5c1per8;
	Scatter2DPtr _histRCPr2pt5c3per8;
	Scatter2DPtr _histRCPr2pt5c5per8;
	
	Scatter2DPtr _histRCPr3pt5c1per8;
	Scatter2DPtr _histRCPr3pt5c3per8;
	Scatter2DPtr _histRCPr3pt5c5per8;
	
	Histo1DPtr _histdummyRCPr2pt5c1;
	Histo1DPtr _histdummyRCPr2pt5c3;
	Histo1DPtr _histdummyRCPr2pt5c5;
	
	Histo1DPtr _histdummyRCPr3pt5c1;
	Histo1DPtr _histdummyRCPr3pt5c3;
	Histo1DPtr _histdummyRCPr3pt5c5;
	
	Histo1DPtr _histdummyRCPr2pt5c8;
	Histo1DPtr _histdummyRCPr3pt5c8;
	Histo1DPtr _histdummyRCPr2pt5c8v2;
	Histo1DPtr _histdummyRCPr3pt5c8v2;
	Histo1DPtr _histdummyRCPr2pt5c8v3;
	Histo1DPtr _histdummyRCPr3pt5c8v3;
///////////////////////////////////////////////////////////////////
	Scatter2DPtr _histRCPvNpr2pt;
	Scatter2DPtr _histRCPvNpr2pt5;
	Scatter2DPtr _histRCPvNpr2pt10;
	
	Scatter2DPtr _histRCPvNpr3pt;
	Scatter2DPtr _histRCPvNpr3pt5;
	Scatter2DPtr _histRCPvNpr3pt10;
////////////////////////////////////////////////////////////////////	

	Scatter2DPtr _histRCPr2ptc1per8;
	Scatter2DPtr _histRCPr2ptc3per8;
	Scatter2DPtr _histRCPr2ptc5per8;
	
	Scatter2DPtr _histRCPr3ptc1per8;
	Scatter2DPtr _histRCPr3ptc3per8;
	Scatter2DPtr _histRCPr3ptc5per8;
	
	Histo1DPtr _histdummyRCPr2ptc1;
	Histo1DPtr _histdummyRCPr2ptc3;
	Histo1DPtr _histdummyRCPr2ptc5;
	
	Histo1DPtr _histdummyRCPr3ptc1;
	Histo1DPtr _histdummyRCPr3ptc3;
	Histo1DPtr _histdummyRCPr3ptc5;
	
	Histo1DPtr _histdummyRCPr2ptc8;
	Histo1DPtr _histdummyRCPr3ptc8;
	Histo1DPtr _histdummyRCPr2ptc8v2;
	Histo1DPtr _histdummyRCPr3ptc8v2;
	Histo1DPtr _histdummyRCPr2ptc8v3;
	Histo1DPtr _histdummyRCPr3ptc8v3;

/////////////////////////////////////////////////////////////////////

	Scatter2DPtr _histRCPr2pt10c1per8;
	Scatter2DPtr _histRCPr2pt10c3per8;
	Scatter2DPtr _histRCPr2pt10c5per8;
	
	Scatter2DPtr _histRCPr3pt10c1per8;
	Scatter2DPtr _histRCPr3pt10c3per8;
	Scatter2DPtr _histRCPr3pt10c5per8;
	
	Histo1DPtr _histdummyRCPr2pt10c1;
	Histo1DPtr _histdummyRCPr2pt10c3;
	Histo1DPtr _histdummyRCPr2pt10c5;
	
	Histo1DPtr _histdummyRCPr3pt10c1;
	Histo1DPtr _histdummyRCPr3pt10c3;
	Histo1DPtr _histdummyRCPr3pt10c5;
	
	Histo1DPtr _histdummyRCPr2pt10c8;
	Histo1DPtr _histdummyRCPr3pt10c8;
	Histo1DPtr _histdummyRCPr2pt10c8v2;
	Histo1DPtr _histdummyRCPr3pt10c8v2;
	Histo1DPtr _histdummyRCPr2pt10c8v3;
	Histo1DPtr _histdummyRCPr3pt10c8v3;

/////////////////////////////////////////////////////////////////////

	Histo1DPtr _histratiocjsr2c1;
	Histo1DPtr _histratiocjsr2c8;
	Histo1DPtr _histratiocjsr3c1;
	Histo1DPtr _histratiocjsr3c8;
	
	Scatter2DPtr _histratiocjsr2by3c1;
	Scatter2DPtr _histratiocjsr2by3c8;
	
	//}
	
	//{
	double _Ncollc1 = 0;
	double _Ncollc3 = 0;
	double _Ncollc5 = 0;
	double _Ncollc8 = 0;
	
	double _Nparttc1 = 0;
	double _Nparttc3 = 0;
	double _Nparttc5 = 0;
	
	double _Npartpc1 = 0;
	double _Npartpc3 = 0;
	double _Npartpc5 = 0;
	//}
	
	void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr out){
		const string path = out->path();
		*out = *h1 / *h2;
		out->setPath(path);
	}
	//}
  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2014_I1263194);
}