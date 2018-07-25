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
      //Initialise and register projections
	  //sets the number of events that will be used to create the centralty binning
	  addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "ImpactParameterMethod");
	  //setting up the cuts 
	  //{
	  //cut of jets with pt 0.15 GeV and higher
	  const ChargedFinalState cfs(Cuts::pT > 0.15*GeV);
	  declare (cfs, "cfs");
	  //cone radius parameter of 0.2 and 0.3 
	  FastJets fjr2(cfs, FastJets::ANTIKT, 0.2);
	  FastJets fjr3(cfs, FastJets::ANTIKT, 0.3);
	  //setting up to find area 
	  fjr2.useJetArea(new fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
	  fjr3.useJetArea(new fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
	  //decalring jets with R=0.2 and R=0.3
	  declare(fjr2, "Jetsr2");
	  declare(fjr3, "Jetsr3");
	  ALICEToolsHI athr2(cfs, fjr2);
	  ALICEToolsHI athr3(cfs, fjr3);
	  declare(athr2, "ATHr2");
	  declare(athr3, "ATHr3");
	  //}
      //Book histograms
	  //{
	  //booking histograms of charged jet spectra with pT > 1.5GeV at 0-10%, 10-30%, 30-50%, 50-80% centralities in R=0.2 and R=0.3
      _hist_cjs_r2_pt0_c1 = bookHisto1D(2, 1, 1);
	  _hist_cjs_r3_pt0_c1 = bookHisto1D(2, 1, 2);

      _hist_cjs_r2_pt0_c3 = bookHisto1D(3, 1, 1);
	  _hist_cjs_r3_pt0_c3 = bookHisto1D(3, 1, 2);
	  
	  _hist_cjs_r2_pt0_c5 = bookHisto1D(4, 1, 1);
	  _hist_cjs_r3_pt0_c5 = bookHisto1D(4, 1, 2);
	  
	  _hist_cjs_r2_pt0_c8 = bookHisto1D(5, 1, 1);
	  _hist_cjs_r3_pt0_c8 = bookHisto1D(5, 1, 2);
	  
	  //booking histograms of charged jet spectra with pT > 5GeV at 0-10%, 10-30%, 30-50%, 50-80% centralities in R=0.2 and R=0.3
	  _hist_cjs_r2_pt5_c1 = bookHisto1D(6, 1, 1);
	  _hist_cjs_r3_pt5_c1 = bookHisto1D(6, 1, 2);
	  
      _hist_cjs_r2_pt5_c3 = bookHisto1D(7, 1, 1);
	  _hist_cjs_r3_pt5_c3 = bookHisto1D(7, 1, 2);
	  
	  _hist_cjs_r2_pt5_c5 = bookHisto1D(8, 1, 1);
	  _hist_cjs_r3_pt5_c5 = bookHisto1D(8, 1, 2);
	  
	  _hist_cjs_r2_pt5_c8 = bookHisto1D(9, 1, 1);
	  _hist_cjs_r3_pt5_c8 = bookHisto1D(9, 1, 2);
	  
	  //booking histograms of charged jet spectra with pT > 10GeV at 0-10%, 10-30%, 30-50%, 50-80% centralities in R=0.2 and R=0.3
	  _hist_cjs_r2_pt10_c1 = bookHisto1D(10, 1, 1);
	  _hist_cjs_r3_pt10_c1 = bookHisto1D(10, 1, 2);
      
	  _hist_cjs_r2_pt10_c3 = bookHisto1D(11, 1, 1);
	  _hist_cjs_r3_pt10_c3 = bookHisto1D(11, 1, 2);
	  
	  _hist_cjs_r2_pt10_c5 = bookHisto1D(12, 1, 1);
	  _hist_cjs_r3_pt10_c5 = bookHisto1D(12, 1, 2);
	  
	  _hist_cjs_r2_pt10_c8 = bookHisto1D(13, 1, 1);
	  _hist_cjs_r3_pt10_c8 = bookHisto1D(13, 1, 2);
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  //booking scatter plots  for ratio of jet spectra at 0-10%, 10-30%, 30-50%, 50-80% centralities in R=0.2
	  _hist_rjs_pt0by5_r2_c1 = bookScatter2D(14, 1, 1);
	  _hist_rjs_pt0by5_r2_c8 = bookScatter2D(15, 1, 1);
	  _hist_rjs_pt10by5_r2_c1 = bookScatter2D(16, 1, 1);
	  _hist_rjs_pt10by5_r2_c8 = bookScatter2D(17, 1, 1);
	  //booking scatter plots  for ratio of jet spectra at 0-10%, 10-30%, 30-50%, 50-80% centralities in R=0.3
	  _hist_rjs_pt0by5_r3_c1 = bookScatter2D(14, 1, 2);
	  _hist_rjs_pt0by5_r3_c8 = bookScatter2D(15, 1, 2);
	  _hist_rjs_pt10by5_r3_c1 = bookScatter2D(16, 1, 2);
	  _hist_rjs_pt10by5_r3_c8 = bookScatter2D(17, 1, 2);
	  //booking dummy histograms for the ratio of jet spectra data 
	  //{
	  //booking dummy histograms of charged jet spectra with pT > 0.15GeV at 0-10% and 50-80% centralities in R=0.2 and R=0.3
	  _histdummypt0r2c1 = bookHisto1D("TMP/dummyhist",  refData(14, 1, 1));
	  _histdummypt0r3c1 = bookHisto1D("TMP/dummyhist2", refData(14, 1, 2));
	  _histdummypt0r2c8 = bookHisto1D("TMP/dummyhist3", refData(15, 1, 1));
	  _histdummypt0r3c8 = bookHisto1D("TMP/dummyhist4", refData(15, 1, 2));
	  //booking dummy histograms of charged jet spectra with pT > 10GeV at 0-10% and 50-80% centralities in R=0.2 and R=0.3
	  _histdummypt10r2c1 = bookHisto1D("TMP/dummyhist5",  refData(16, 1, 1));
	  _histdummypt10r3c1 = bookHisto1D("TMP/dummyhist6", refData(16, 1, 2));
	  _histdummypt10r2c8 = bookHisto1D("TMP/dummyhist7", refData(17, 1, 1));
	  _histdummypt10r3c8 = bookHisto1D("TMP/dummyhist8", refData(17, 1, 2));
	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 0-10% and 50-80% centralities in R=0.2 and R=0.3
	  _histdummypt5r2c1  = bookHisto1D("TMP/dummyhist9", refData(14, 1, 1));
	  _histdummypt5r3c1  = bookHisto1D("TMP/dummyhist10", refData(14, 1, 2));
	  _histdummypt5r2c8  = bookHisto1D("TMP/dummyhist11", refData(15, 1, 1));
	  _histdummypt5r3c8  = bookHisto1D("TMP/dummyhist12", refData(15, 1, 2));
	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 0-10% and 50-80% centralities in R=0.2 and R=0.3
	  //there is another version of the above graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummypt5r2c1v2 = bookHisto1D("TMP/dummyhist13", refData(16, 1, 1));
	  _histdummypt5r3c1v2 = bookHisto1D("TMP/dummyhist14", refData(16, 1, 2));
	  _histdummypt5r2c8v2 = bookHisto1D("TMP/dummyhist15", refData(17, 1, 1));
	  _histdummypt5r3c8v2 = bookHisto1D("TMP/dummyhist16", refData(17, 1, 2));
	  //}
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  //booking scatter plots for the nuclear mod factor(RCP) at pt>5GeV 0-10%, 10-30%, 30-50% centralities divided by 50-80% perhiperal in R=0.2
	  _hist_RCP_r2_pt5_c1per8 = bookScatter2D(21, 1, 1);
	  _hist_RCP_r2_pt5_c3per8 = bookScatter2D(22, 1, 1);
	  _hist_RCP_r2_pt5_c5per8 = bookScatter2D(23, 1, 1);
	  //booking scatter plots for the nuclear mod factor(RCP) at pt>5GeV 0-10%, 10-30%, 30-50% centralities divided by 50-80% perhiperal  in R=0.3
	  _hist_RCP_r3_pt5_c1per8 = bookScatter2D(21, 1, 2);
	  _hist_RCP_r3_pt5_c3per8 = bookScatter2D(22, 1, 2);
	  _hist_RCP_r3_pt5_c5per8 = bookScatter2D(23, 1, 2);
	  //booking dummy histgrams for RCP 
	  //{
	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 0-10%, 10-30%, and 30-50% centralities in R=0.2
	  _histdummyRCPr2pt5c1 = bookHisto1D("TMP/dummyhistRCPr2pt5c1", refData(21, 1, 1));
	  _histdummyRCPr2pt5c3 = bookHisto1D("TMP/dummyhistRCPr2pt5c3", refData(22, 1, 1));
	  _histdummyRCPr2pt5c5 = bookHisto1D("TMP/dummyhistRCPr2pt5c5", refData(23, 1, 1));
	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 0-10%, 10-30%, and 30-50% centralities in R=0.3
	  _histdummyRCPr3pt5c1 = bookHisto1D("TMP/dummyhistRCPr3pt5c1", refData(21, 1, 2));
	  _histdummyRCPr3pt5c3 = bookHisto1D("TMP/dummyhistRCPr3pt5c3", refData(22, 1, 2));
	  _histdummyRCPr3pt5c5 = bookHisto1D("TMP/dummyhistRCPr3pt5c5", refData(23, 1, 2));
	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 50-80% centralities in R=0.2
	  //there is another version of the following graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummyRCPr2pt5c8 = bookHisto1D("TMP/dummyhistRCPr2c8", refData(21, 1, 1));
	  _histdummyRCPr2pt5c8v2 = bookHisto1D("TMP/dummyhistRCPr28v2", refData(22, 1, 1));
	  _histdummyRCPr2pt5c8v3 = bookHisto1D("TMP/dummyhistRCPr28v3", refData(23, 1, 1));
	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 50-80% centralities in R=0.3
	  //there is another version of the following graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummyRCPr3pt5c8 = bookHisto1D("TMP/dummyhistRCPr3c8", refData(21, 1, 2));
	  _histdummyRCPr3pt5c8v2 = bookHisto1D("TMP/dummyhistRCPr3c8v2", refData(22, 1, 2));
	  _histdummyRCPr3pt5c8v3 = bookHisto1D("TMP/dummyhistRCPr3c8v3", refData(23, 1, 2));
	  //}
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  //booking scatter plots for Nulcear modification factor(RCP) vs Number of particle(Npart) for pt>0.15GeV, pt>5GeV, pt>10GeV in R=0.2
	  _hist_RCP_vs_Np_r2_pt0 = bookScatter2D(27, 1, 1);
	  _hist_RCP_vs_Np_r2_pt5 = bookScatter2D(28, 1, 1);
	  _hist_RCP_vs_Np_r2_pt10 = bookScatter2D(29, 1, 1);
	  //booking scatter plots for Nulcear modification factor(RCP) vs Number of particle(Npart) for pt>0.15GeV, pt>5GeV, pt>10GeV in R=0.3
	  _hist_RCP_vs_Np_r3_pt0 = bookScatter2D(27, 1, 2);
	  _hist_RCP_vs_Np_r3_pt5 = bookScatter2D(28, 1, 2);
	  _hist_RCP_vs_Np_r3_pt10 = bookScatter2D(29, 1, 2);
	  //////////////////////////////////////////////////////////////////////////////////////////
	  //booking scatter plots for the nuclear mod factor(RCP) at pt>0.15GeV 0-10%, 10-30%, 30-50% centralities divided by 50-80% perhiperal in R=0.2
	  _hist_RCP_r2_pt0_c1per8 = bookScatter2D(18, 1, 1);
	  _hist_RCP_r2_pt0_c3per8 = bookScatter2D(19, 1, 1);
	  _hist_RCP_r2_pt0_c5per8 = bookScatter2D(20, 1, 1);
	  //booking scatter plots for the nuclear mod factor(RCP) at pt>0.15GeV 0-10%, 10-30%, 30-50% centralities divided by 50-80% perhiperal in R=0.3
	  _hist_RCP_r3_pt0_c1per8 = bookScatter2D(18, 1, 2);
	  _hist_RCP_r3_pt0_c3per8 = bookScatter2D(19, 1, 2);
	  _hist_RCP_r3_pt0_c5per8 = bookScatter2D(20, 1, 2);
	  //booking dummy histgrams for RCP 
	  //{
	  //booking dummy histograms of charged jet spectra with pT > 0.15GeV at 0-10%, 10-30%, and 30-50% centralities in R=0.2
	  _histdummyRCPr2pt0c1 = bookHisto1D("TMP/dummyhistRCPr2ptc1", refData(18, 1, 1));
	  _histdummyRCPr2pt0c3 = bookHisto1D("TMP/dummyhistRCPr2ptc3", refData(19, 1, 1));
	  _histdummyRCPr2pt0c5 = bookHisto1D("TMP/dummyhistRCPr2ptc5", refData(20, 1, 1));
	  //booking dummy histograms of charged jet spectra with pT > 0.15GeV at 0-10%, 10-30%, and 30-50% centralities in R=0.3
	  _histdummyRCPr3pt0c1 = bookHisto1D("TMP/dummyhistRCPr3ptc1", refData(18, 1, 2));
	  _histdummyRCPr3pt0c3 = bookHisto1D("TMP/dummyhistRCPr3ptc3", refData(19, 1, 2));
	  _histdummyRCPr3pt0c5 = bookHisto1D("TMP/dummyhistRCPr3ptc5", refData(20, 1, 2));
	  //booking dummy histograms of charged jet spectra with pT > 0.15GeV at 50-80% centralities in R=0.2
	  //there is another version of the following graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummyRCPr2pt0c8 = bookHisto1D("TMP/dummyhistRCPr2c8pt", refData(18, 1, 1));
	  _histdummyRCPr2pt0c8v2 = bookHisto1D("TMP/dummyhistRCPr28v2pt", refData(19, 1, 1));
	  _histdummyRCPr2pt0c8v3 = bookHisto1D("TMP/dummyhistRCPr28v3pt", refData(20, 1, 1));
	  //booking dummy histograms of charged jet spectra with pT > 0.15GeV at 50-80% centralities in R=0.3
	  //there is another version of the following graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummyRCPr3pt0c8 = bookHisto1D("TMP/dummyhistRCPr3c8pt", refData(18, 1, 2));
	  _histdummyRCPr3pt0c8v2 = bookHisto1D("TMP/dummyhistRCPr3c8v2pt", refData(19, 1, 2));
	  _histdummyRCPr3pt0c8v3 = bookHisto1D("TMP/dummyhistRCPr3c8v3pt", refData(20, 1, 2));
	  //}
	  ////////////////////////////////////////////////////////////////////////////////////////////
	  //booking scatter plots for the nuclear mod factor(RCP) at pt>10GeV 0-10%, 10-30%, 30-50% centralities divided by 50-80% perhiperal in R=0.2
	  _hist_RCP_r2_pt10_c1per8 = bookScatter2D(24, 1, 1);
	  _hist_RCP_r2_pt10_c3per8 = bookScatter2D(25, 1, 1);
	  _hist_RCP_r2_pt10_c5per8 = bookScatter2D(26, 1, 1);
	  //booking scatter plots for the nuclear mod factor(RCP) at pt>10GeV 0-10%, 10-30%, 30-50% centralities divided by 50-80% perhiperal in R=0.3
	  _hist_RCP_r3_pt10_c1per8 = bookScatter2D(24, 1, 2);
	  _hist_RCP_r3_pt10_c3per8 = bookScatter2D(25, 1, 2);
	  _hist_RCP_r3_pt10_c5per8 = bookScatter2D(26, 1, 2);
	  //booking dummy histgrams for RCP 
	  //{
  	  //booking dummy histograms of charged jet spectra with pT > 10GeV at 0-10%, 10-30%, and 30-50% centralities in R=0.2
	  _histdummyRCPr2pt10c1 = bookHisto1D("TMP/dummyhistRCPr2pt10c1", refData(24, 1, 1));
	  _histdummyRCPr2pt10c3 = bookHisto1D("TMP/dummyhistRCPr2pt10c3", refData(25, 1, 1));
	  _histdummyRCPr2pt10c5 = bookHisto1D("TMP/dummyhistRCPr2pt10c5", refData(26, 1, 1));
	  // booking dummy histograms of charged jet spectra with pT > 10GeV at 0-10%, 10-30%, and 30-50% centralities in R=0.3
	  _histdummyRCPr3pt10c1 = bookHisto1D("TMP/dummyhistRCPr3pt10c1", refData(24, 1, 2));
	  _histdummyRCPr3pt10c3 = bookHisto1D("TMP/dummyhistRCPr3pt10c3", refData(25, 1, 2));
	  _histdummyRCPr3pt10c5 = bookHisto1D("TMP/dummyhistRCPr3pt10c5", refData(26, 1, 2));
	  //booking dummy histograms of charged jet spectra with pT > 10GeV at 50-80% centralities in R=0.2
	  //there is another version of the following graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummyRCPr2pt10c8 = bookHisto1D("TMP/dummyhistRCPr2c8pt10", refData(24, 1, 1));
	  _histdummyRCPr2pt10c8v2 = bookHisto1D("TMP/dummyhistRCPr28v2pt10", refData(25, 1, 1));
	  _histdummyRCPr2pt10c8v3 = bookHisto1D("TMP/dummyhistRCPr28v3pt10", refData(26, 1, 1));
	  //booking dummy histograms of charged jet spectra with pT > 10GeV at 50-80% centralities in R=0.3
	  //there is another version of the following graphs due to a binning issue if you divide 2 histograms with differnet refrence data
	  _histdummyRCPr3pt10c8 = bookHisto1D("TMP/dummyhistRCPr3c8pt10", refData(24, 1, 2));
	  _histdummyRCPr3pt10c8v2 = bookHisto1D("TMP/dummyhistRCPr3c8v2pt10", refData(25, 1, 2));
	  _histdummyRCPr3pt10c8v3 = bookHisto1D("TMP/dummyhistRCPr3c8v3pt10", refData(26, 1, 2));
	  //}
	  //booking scatter plots of ratio of R=0.2/R=0.3 charge jet spcetra in 0-10% and 50-80% centralities at pt>5GeV
	  _hist_ratio_cjs_r2by3_c1 = bookScatter2D(30, 1, 1);
	  _hist_ratio_cjs_r2by3_c8 = bookScatter2D(30, 1, 2);
	  //booking dummy histgrams for ratio of R=0.2/R=0.3 
	  //{
  	  //booking dummy histograms of charged jet spectra with pT > 5GeV at 0-10% and 50-80% centralities in R=0.2 and R=0.3
	  _histdummy_ratio_cjs_r2_c1 = bookHisto1D("TMP/dumratiocjs1", refData(30, 1, 1));
	  _histdummy_ratio_cjs_r3_c1 = bookHisto1D("TMP/dumratiocjs2", refData(30, 1, 1));
	  _histdummy_ratio_cjs_r2_c8 = bookHisto1D("TMP/dumratiocjs3", refData(30, 1, 2));
	  _histdummy_ratio_cjs_r3_c8 = bookHisto1D("TMP/dumratiocjs4", refData(30, 1, 2));
	  //}
	  //}
    }

    void analyze(const Event& event) 
	{
		//vetos any event that is not within 0-80% centralties 
		const double c = centrality(event, "ImpactParameterMethod");
		cout << c << "cent\n";
		if ((c < 0.)||(c>80.)){
			vetoEvent;
		}
		//decalring weight for events
		const double weight = 1.0;
		//decalring fast jet constants 
		//{
		const FastJets& fastjetr2 = apply<FastJets>(event, "Jetsr2");
		const FastJets& fastjetr3 = apply<FastJets>(event, "Jetsr3");
		//setting up the jets to only b looked at if pT is higher than 0.15GeV
		const Jets& jets2 = fastjetr2.jetsByPt(0.15*GeV);
		const Jets& jets3 = fastjetr3.jetsByPt(0.15*GeV);
		//decalres variabels for fining area and rho
		const auto seqr2 = fastjetr2.clusterSeqArea();
		const auto seqr3 = fastjetr3.clusterSeqArea();
		const ALICEToolsHI &athr2 = apply<ALICEToolsHI>(event, "ATHr2");
		const ALICEToolsHI &athr3 = apply<ALICEToolsHI>(event, "ATHr3");
		//}
		// analysis and fillimg histograms
		// runs through all of the events within 0-10% centrality 
		if (c >= 0 && c <= 10) {
			_eventcounter1++; //counts the number of events within this centrality 
			//runs through all of the R=0.2 jets
			foreach(const Jet& jets, jets2)
			{
				//gets area and rho0 for the background subtraction equation 
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				//fills histograms with jets with pt > 0.15GeV 
				_hist_cjs_r2_pt0_c1->fill(jetpt, weight);
				_histdummypt0r2c1->fill(jetpt, weight);
				_histdummyRCPr2pt0c1->fill(jetpt, weight);
				//creates varaible to find the constituents with a certain pt range
				const Particles ps = jets.constituents(); 
				//fills histograms with jets with pt > 5GeV
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r2_pt5_c1->fill(jetpt, weight);
					_histdummypt5r2c1->fill(jetpt, weight);
					_histdummypt5r2c1v2->fill(jetpt, weight);
					_histdummyRCPr2pt5c1->fill(jetpt, weight);
					_histdummy_ratio_cjs_r2_c1->fill(jetpt, weight);
				}
				//fills histograms with jets with pt > 10GeV
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r2_pt10_c1->fill(jetpt, weight);
					_histdummypt10r2c1->fill(jetpt, weight);
					_histdummyRCPr2pt10c1->fill(jetpt, weight);
				}
			}	
			//runs through all of the R=0.3 jets
			foreach(const Jet& jets, jets3)
			{
				//gets area and rho0 for the background subtraction equation 
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				//fills histograms with jets with pt > 0.15GeV
				_hist_cjs_r3_pt0_c1->fill(jetpt, weight);
				_histdummypt0r3c1->fill(jetpt, weight);
				_histdummyRCPr3pt0c1->fill(jetpt, weight);
				//creates varaible to find the constituents with a certain pt range
				const Particles ps = jets.constituents();
				//fills histograms with jets with pt > 5GeV
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r3_pt5_c1->fill(jetpt, weight);
					_histdummypt5r3c1->fill(jetpt, weight);
					_histdummypt5r3c1v2->fill(jetpt, weight);
					_histdummyRCPr3pt5c1->fill(jetpt, weight);
					_histdummy_ratio_cjs_r3_c1->fill(jetpt, weight);
				}
				//fills histograms with jets with pt > 10GeV
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r3_pt10_c1->fill(jetpt, weight);
					_histdummypt10r3c1->fill(jetpt, weight);
					_histdummyRCPr3pt10c1->fill(jetpt, weight);
				}	
			}	
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// runs through all of the events within 10-30% centrality 
		if (10 < c && c <= 30){
			_eventcounter3++;
			foreach (const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_hist_cjs_r2_pt0_c3->fill(jetpt,weight);
				_histdummyRCPr2pt0c3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r2_pt5_c3->fill(jetpt, weight);
					_histdummyRCPr2pt5c3->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r2_pt10_c3->fill(jetpt, weight);
					_histdummyRCPr2pt10c3->fill(jetpt, weight);
				}
			}
			foreach (const Jet& jets, jets3)
			{
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_hist_cjs_r3_pt0_c3->fill(jetpt,weight);
				_histdummyRCPr3pt0c3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r3_pt5_c3->fill(jetpt, weight);
					_histdummyRCPr3pt5c3->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r3_pt10_c3->fill(jetpt, weight);
					_histdummyRCPr3pt10c3->fill(jetpt, weight);
				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// runs through all of the evets within 30-50% centrality 
		if ( 30 < c && c <= 50){
			_eventcounter5++;
			foreach (const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_hist_cjs_r2_pt0_c5->fill(jetpt,weight);
				_histdummyRCPr2pt0c5->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r2_pt5_c5->fill(jetpt, weight);
					_histdummyRCPr2pt5c5->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r2_pt10_c5->fill(jetpt, weight);
					_histdummyRCPr2pt10c5->fill(jetpt, weight);
				}
			}
			foreach (const Jet& jets, jets3)
			{
				double area = seqr3->area(jets);
				double rho0 = athr3.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_hist_cjs_r3_pt0_c5->fill(jetpt,weight);
				_histdummyRCPr3pt0c5->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r3_pt5_c5->fill(jetpt, weight);
					_histdummyRCPr3pt5c5->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r3_pt10_c5->fill(jetpt, weight);
					_histdummyRCPr3pt10c5->fill(jetpt, weight);
				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// runs through all of the events within 50-80% centrality 
		if ( 50 < c && c <= 80){
			_eventcounter8++;
			foreach (const Jet& jets, jets2)
			{
				double area = seqr2->area(jets);
				double rho0 = athr2.RhoLocal(jets.phi());
				double jetpt = jets.pT() - area*rho0;
				_hist_cjs_r2_pt0_c8->fill(jetpt,weight);
				_histdummypt0r2c8->fill(jetpt, weight);
				_histdummyRCPr2pt0c8->fill(jetpt, weight);
				_histdummyRCPr2pt0c8v2->fill(jetpt, weight);
				_histdummyRCPr2pt0c8v3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r2_pt5_c8->fill(jetpt, weight);
					_histdummypt5r2c8->fill(jetpt, weight);
					_histdummypt5r2c8v2->fill(jetpt, weight);
					_histdummyRCPr2pt5c8->fill(jetpt, weight);
					_histdummyRCPr2pt5c8v2->fill(jetpt, weight);
					_histdummyRCPr2pt5c8v3->fill(jetpt, weight);
					_histdummy_ratio_cjs_r2_c8->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r2_pt10_c8->fill(jetpt, weight);
					_histdummypt10r2c8->fill(jetpt, weight);
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
				_hist_cjs_r3_pt0_c8->fill(jetpt,weight);
				_histdummypt0r3c8->fill(jetpt, weight);
				_histdummyRCPr3pt0c8->fill(jetpt, weight);
				_histdummyRCPr3pt0c8v2->fill(jetpt, weight);
				_histdummyRCPr3pt0c8v3->fill(jetpt, weight);
				const Particles ps = jets.constituents();
				if (any(ps, PtGtr(5*GeV))){
					_hist_cjs_r3_pt5_c8->fill(jetpt, weight);
					_histdummypt5r3c8->fill(jetpt, weight);
					_histdummypt5r3c8v2->fill(jetpt, weight);
					_histdummyRCPr3pt5c8->fill(jetpt, weight);
					_histdummyRCPr3pt5c8v2->fill(jetpt, weight);
					_histdummyRCPr3pt5c8v3->fill(jetpt, weight);
					_histdummy_ratio_cjs_r3_c8->fill(jetpt, weight);
				}
				if (any(ps, PtGtr(10*GeV))){
					_hist_cjs_r3_pt10_c8->fill(jetpt, weight);
					_histdummypt10r3c8->fill(jetpt, weight);
					_histdummyRCPr3pt10c8->fill(jetpt, weight);
					_histdummyRCPr3pt10c8v2->fill(jetpt, weight);
					_histdummyRCPr3pt10c8v3->fill(jetpt, weight);

				}
			}
		}
    }
    
	/// Normalise histograms etc., after the run
    void finalize() {
		//scaling the historgrams
		// 0-10% centrality 
		//{
        _evntandNcoll1 = _eventcounter1*_ncollc1;
        _evntandNcoll3 = _eventcounter3*_ncollc3;
        _evntandNcoll5 = _eventcounter5*_ncollc5;
        _evntandNcoll8 = _eventcounter8*_ncollc8;
		scale(_hist_cjs_r2_pt0_c1, 1/_evntandNcoll1);
		scale(_hist_cjs_r3_pt0_c1, 1/_evntandNcoll1);
		
		scale(_hist_cjs_r2_pt5_c1,  1/_evntandNcoll1);
		scale(_hist_cjs_r3_pt5_c1,  1/_evntandNcoll1);
		
		scale(_hist_cjs_r2_pt10_c1,  1/_evntandNcoll1);
		scale(_hist_cjs_r3_pt10_c1,  1/_evntandNcoll1);
		//23.5 is the TAA value for 0-10% centrality 
		scale(_histdummyRCPr2pt5c1,  1/(23.5*_eventcounter1));
		scale(_histdummyRCPr3pt5c1,  1/(23.5*_eventcounter1));
		
		scale(_histdummyRCPr2pt0c1,  1/(23.5*_eventcounter1));
		scale(_histdummyRCPr3pt0c1,  1/(23.5*_eventcounter1));
		
		scale(_histdummyRCPr2pt10c1,  1/(23.5*_eventcounter1));
		scale(_histdummyRCPr3pt10c1,  1/(23.5*_eventcounter1));
		
		scale(_histdummy_ratio_cjs_r2_c1, 1/_evntandNcoll1);
		scale(_histdummy_ratio_cjs_r3_c1, 1/_evntandNcoll1);
		
		scale(_histdummypt0r2c1, 1/_evntandNcoll1);
		scale(_histdummypt5r2c1, 1/_evntandNcoll1);
		scale(_histdummypt10r2c1, 1/_evntandNcoll1);
		scale(_histdummypt5r2c1v2, 1/_evntandNcoll1);
		scale(_histdummypt0r3c1, 1/_evntandNcoll1);
		scale(_histdummypt5r3c1, 1/_evntandNcoll1);
		scale(_histdummypt10r3c1, 1/_evntandNcoll1);
		scale(_histdummypt5r3c1v2, 1/_evntandNcoll1);
		//}
		//}
		
		// 10-30% centrality 
		//{
		scale(_hist_cjs_r2_pt0_c3, 1/_evntandNcoll3);
		scale(_hist_cjs_r3_pt0_c3, 1/_evntandNcoll3);
		
		scale(_hist_cjs_r2_pt5_c3,  1/_evntandNcoll3);
		scale(_hist_cjs_r3_pt5_c3,  1/_evntandNcoll3);
		
		scale(_hist_cjs_r2_pt10_c3,  1/_evntandNcoll3);
		scale(_hist_cjs_r3_pt10_c3,  1/_evntandNcoll3);
		//11.6 is the TAA value for 10-30% centrality
		scale(_histdummyRCPr2pt5c3,  1/(11.6*_eventcounter3));
		scale(_histdummyRCPr3pt5c3,  1/(11.6*_eventcounter3));
		
		scale(_histdummyRCPr2pt0c3,  1/(11.6*_eventcounter3));
		scale(_histdummyRCPr3pt0c3,  1/(11.6*_eventcounter3));
		
		scale(_histdummyRCPr2pt10c3,  1/(11.6*_eventcounter3));
		scale(_histdummyRCPr3pt10c3,  1/(11.6*_eventcounter3));
		//}
		
		// 30-50% centrality 
		//{
		scale(_hist_cjs_r2_pt0_c5, 1/_evntandNcoll5);
		scale(_hist_cjs_r3_pt0_c5, 1/_evntandNcoll5);
		
		scale(_hist_cjs_r2_pt5_c5,  1/_evntandNcoll5);
		scale(_hist_cjs_r3_pt5_c5,  1/_evntandNcoll5);
		
		scale(_hist_cjs_r2_pt10_c5,  1/_evntandNcoll5);
		scale(_hist_cjs_r3_pt10_c5,  1/_evntandNcoll5);
		//3.8 is the TAA value for 30-50% centrality
		scale(_histdummyRCPr2pt5c5,  1/(3.8*_eventcounter5));
		scale(_histdummyRCPr3pt5c5,  1/(3.8*_eventcounter5));
		
		scale(_histdummyRCPr2pt0c5,  1/(3.8*_eventcounter5));
		scale(_histdummyRCPr3pt0c5,  1/(3.8*_eventcounter5));
		
		scale(_histdummyRCPr2pt10c5,  1/(3.8*_eventcounter5));
		scale(_histdummyRCPr3pt10c5,  1/(3.8*_eventcounter5));
		//}
		
		// 50-80% centrality 
		//{
		scale(_hist_cjs_r2_pt0_c8, 1/_evntandNcoll8);
		scale(_hist_cjs_r3_pt0_c8, 1/_evntandNcoll8);
		
		scale(_hist_cjs_r2_pt5_c8,  1/_evntandNcoll8);
		scale(_hist_cjs_r3_pt5_c8,  1/_evntandNcoll8);
		
		scale(_hist_cjs_r2_pt10_c8,  1/_evntandNcoll8);
		scale(_hist_cjs_r3_pt10_c8,  1/_evntandNcoll8);
		//0.70 is the TAA value for 50-80% centrality
		scale(_histdummyRCPr2pt5c8,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr2pt5c8v2,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr2pt5c8v3,  1/(0.70*_eventcounter8));
		
		scale(_histdummyRCPr3pt5c8,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr3pt5c8v2, 1/(0.70*_eventcounter8));
		scale(_histdummyRCPr3pt5c8v3, 1/(0.70*_eventcounter8));
		
		scale(_histdummyRCPr2pt0c8, 1/(0.70*_eventcounter8));
		scale(_histdummyRCPr2pt0c8v2,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr2pt0c8v3,  1/(0.70*_eventcounter8));
		
		scale(_histdummyRCPr3pt0c8,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr3pt0c8v2,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr3pt0c8v3,  1/(0.70*_eventcounter8));
		
		scale(_histdummyRCPr2pt10c8,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr2pt10c8v2,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr2pt10c8v3,  1/(0.70*_eventcounter8));
		
		scale(_histdummyRCPr3pt10c8,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr3pt10c8v2,  1/(0.70*_eventcounter8));
		scale(_histdummyRCPr3pt10c8v3,  1/(0.70*_eventcounter8));
		
		scale(_histdummy_ratio_cjs_r2_c8, 1/_evntandNcoll8);
		scale(_histdummy_ratio_cjs_r3_c8, 1/_evntandNcoll8);
		scale(_histdummypt0r2c8, 1/_evntandNcoll8);
		scale(_histdummypt10r2c8, 1/_evntandNcoll8);
		scale(_histdummypt0r3c8, 1/_evntandNcoll8);
		scale(_histdummypt10r3c8, 1/_evntandNcoll8);
		scale(_histdummypt5r2c8, 1/_evntandNcoll8);
		scale(_histdummypt5r2c8v2, 1/_evntandNcoll8);
		scale(_histdummypt5r3c8, 1/_evntandNcoll8);
		scale(_histdummypt5r3c8v2, 1/_evntandNcoll8);
		//}
		//dividing histograms to fill scatter plots
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummypt0r2c1, _histdummypt5r2c1, _hist_rjs_pt0by5_r2_c1);
		divide(_histdummypt0r2c8, _histdummypt5r2c8, _hist_rjs_pt0by5_r2_c8);
		divide(_histdummypt10r2c1, _histdummypt5r2c1v2, _hist_rjs_pt10by5_r2_c1);
		divide(_histdummypt10r2c8, _histdummypt5r2c8v2, _hist_rjs_pt10by5_r2_c8);
		
		divide(_histdummypt0r3c1, _histdummypt5r3c1, _hist_rjs_pt0by5_r3_c1);
		divide(_histdummypt0r3c8, _histdummypt5r3c8, _hist_rjs_pt0by5_r3_c8);
		divide(_histdummypt10r3c1, _histdummypt5r3c1v2, _hist_rjs_pt10by5_r3_c1);
		divide(_histdummypt10r3c8, _histdummypt5r3c8v2, _hist_rjs_pt10by5_r3_c8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummyRCPr2pt5c1, _histdummyRCPr2pt5c8, _hist_RCP_r2_pt5_c1per8);
		divide(_histdummyRCPr2pt5c3, _histdummyRCPr2pt5c8v2, _hist_RCP_r2_pt5_c3per8);
		divide(_histdummyRCPr2pt5c5, _histdummyRCPr2pt5c8v3, _hist_RCP_r2_pt5_c5per8);
		
		divide(_histdummyRCPr3pt5c1, _histdummyRCPr3pt5c8, _hist_RCP_r3_pt5_c1per8);
		divide(_histdummyRCPr3pt5c3, _histdummyRCPr3pt5c8v2, _hist_RCP_r3_pt5_c3per8);
		divide(_histdummyRCPr3pt5c5, _histdummyRCPr3pt5c8v3, _hist_RCP_r3_pt5_c5per8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummyRCPr2pt0c1, _histdummyRCPr2pt0c8, _hist_RCP_r2_pt0_c1per8);
		divide(_histdummyRCPr2pt0c3, _histdummyRCPr2pt0c8v2, _hist_RCP_r2_pt0_c3per8);
		divide(_histdummyRCPr2pt0c5, _histdummyRCPr2pt0c8v3, _hist_RCP_r2_pt0_c5per8);
		
		divide(_histdummyRCPr3pt0c1, _histdummyRCPr3pt0c8, _hist_RCP_r3_pt0_c1per8);
		divide(_histdummyRCPr3pt0c3, _histdummyRCPr3pt0c8v2, _hist_RCP_r3_pt0_c3per8);
		divide(_histdummyRCPr3pt0c5, _histdummyRCPr3pt0c8v3, _hist_RCP_r3_pt0_c5per8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummyRCPr2pt10c1, _histdummyRCPr2pt10c8, _hist_RCP_r2_pt10_c1per8);
		divide(_histdummyRCPr2pt10c3, _histdummyRCPr2pt10c8v2, _hist_RCP_r2_pt10_c3per8);
		divide(_histdummyRCPr2pt10c5, _histdummyRCPr2pt10c8v3, _hist_RCP_r2_pt10_c5per8);
		
		divide(_histdummyRCPr3pt10c1, _histdummyRCPr3pt10c8, _hist_RCP_r3_pt10_c1per8);
		divide(_histdummyRCPr3pt10c3, _histdummyRCPr3pt10c8v2, _hist_RCP_r3_pt10_c3per8);
		divide(_histdummyRCPr3pt10c5, _histdummyRCPr3pt10c8v3, _hist_RCP_r3_pt10_c5per8);
////////////////////////////////////////////////////////////////////////////////////////////////
		divide(_histdummy_ratio_cjs_r2_c1, _histdummy_ratio_cjs_r3_c1, _hist_ratio_cjs_r2by3_c1);
		divide(_histdummy_ratio_cjs_r2_c8, _histdummy_ratio_cjs_r3_c8, _hist_ratio_cjs_r2by3_c8);
////////////////////////////////////////////////////////////////////////////////////////////////
		//}
		//adding points to the RCP vs Npart plots 
		//{
		//WARNING the error bars in the add point are wrong, they are merely there so that the point will be visable in the graphs 
////////////////////////////////////////////////////////////////////////////////////////////////
		_hist_RCP_vs_Np_r2_pt0->addPoint(_npartc1, _hist_RCP_r2_pt0_c1per8->point(3).y(), 20, 0);
		_hist_RCP_vs_Np_r2_pt0->addPoint(_npartc3, _hist_RCP_r2_pt0_c3per8->point(3).y(), 20, 0);
		_hist_RCP_vs_Np_r2_pt0->addPoint(_npartc5, _hist_RCP_r2_pt0_c5per8->point(4).y(), 20, 0);
		
		_hist_RCP_vs_Np_r3_pt0->addPoint(_npartc1, _hist_RCP_r3_pt0_c1per8->point(2).y(), 20, 0);
		_hist_RCP_vs_Np_r3_pt0->addPoint(_npartc3, _hist_RCP_r3_pt0_c3per8->point(2).y(), 20, 0);
		_hist_RCP_vs_Np_r3_pt0->addPoint(_npartc5, _hist_RCP_r3_pt0_c5per8->point(3).y(), 20, 0);
////////////////////////////////////////////////////////////////////////////////////////////////		
		_hist_RCP_vs_Np_r2_pt5->addPoint(_npartc1, _hist_RCP_r2_pt5_c1per8->point(4).y(), 20, 0);
		_hist_RCP_vs_Np_r2_pt5->addPoint(_npartc3, _hist_RCP_r2_pt5_c3per8->point(4).y(), 20, 0);
		_hist_RCP_vs_Np_r2_pt5->addPoint(_npartc5, _hist_RCP_r2_pt5_c5per8->point(4).y(), 20, 0);
		
		_hist_RCP_vs_Np_r3_pt5->addPoint(_npartc1, _hist_RCP_r3_pt5_c1per8->point(3).y(), 20, 0);
		_hist_RCP_vs_Np_r3_pt5->addPoint(_npartc3, _hist_RCP_r3_pt5_c3per8->point(3).y(), 20, 0);
		_hist_RCP_vs_Np_r3_pt5->addPoint(_npartc5, _hist_RCP_r3_pt5_c5per8->point(4).y(), 20, 0);
////////////////////////////////////////////////////////////////////////////////////////////////
		_hist_RCP_vs_Np_r2_pt10->addPoint(_npartc1, _hist_RCP_r2_pt10_c1per8->point(4).y(), 20, 0);
		_hist_RCP_vs_Np_r2_pt10->addPoint(_npartc3, _hist_RCP_r2_pt10_c3per8->point(4).y(), 20, 0);
		_hist_RCP_vs_Np_r2_pt10->addPoint(_npartc5, _hist_RCP_r2_pt10_c5per8->point(4).y(), 20, 0);
		
		_hist_RCP_vs_Np_r3_pt10->addPoint(_npartc1, _hist_RCP_r3_pt10_c1per8->point(4).y(), 20, 0);
		_hist_RCP_vs_Np_r3_pt10->addPoint(_npartc3, _hist_RCP_r3_pt10_c3per8->point(4).y(), 20, 0);
		_hist_RCP_vs_Np_r3_pt10->addPoint(_npartc5, _hist_RCP_r3_pt10_c5per8->point(4).y(), 20, 0);
		//}
    }
	// decalring event counters for each centrality range
	//{
	double _eventcounter1 = 0;
	double _eventcounter3 = 0;
	double _eventcounter5 = 0;
	double _eventcounter8 = 0;
    /// Perform the per-event analysis
	//}
	//declareing histograms and scatter plots 
	//{
	//pT > 1.5 GeV/c at different centrlities 
	Histo1DPtr _hist_cjs_r2_pt0_c1;
	Histo1DPtr _hist_cjs_r3_pt0_c1;
	
	Histo1DPtr _hist_cjs_r2_pt0_c3;
	Histo1DPtr _hist_cjs_r3_pt0_c3;
	
	Histo1DPtr _hist_cjs_r2_pt0_c5;
	Histo1DPtr _hist_cjs_r3_pt0_c5;
	
	Histo1DPtr _hist_cjs_r2_pt0_c8;
	Histo1DPtr _hist_cjs_r3_pt0_c8;
	// histograms with pT > 5 GeV/c at different centrlities 
	Histo1DPtr _hist_cjs_r2_pt5_c1;
	Histo1DPtr _hist_cjs_r3_pt5_c1;
	
	Histo1DPtr _hist_cjs_r2_pt5_c3;
	Histo1DPtr _hist_cjs_r3_pt5_c3;
	
	Histo1DPtr _hist_cjs_r2_pt5_c5;
	Histo1DPtr _hist_cjs_r3_pt5_c5;
	
	Histo1DPtr _hist_cjs_r2_pt5_c8;
	Histo1DPtr _hist_cjs_r3_pt5_c8;
	// histograms with pT > 10 GeV/c at different centrlities 
	Histo1DPtr _hist_cjs_r2_pt10_c1;
	Histo1DPtr _hist_cjs_r3_pt10_c1;
	
	Histo1DPtr _hist_cjs_r2_pt10_c3;
	Histo1DPtr _hist_cjs_r3_pt10_c3;
	
	Histo1DPtr _hist_cjs_r2_pt10_c5;
	Histo1DPtr _hist_cjs_r3_pt10_c5;
	
	Histo1DPtr _hist_cjs_r2_pt10_c8;
	Histo1DPtr _hist_cjs_r3_pt10_c8;
///////////////////////////////////////////////////////////////////
	Scatter2DPtr  _hist_rjs_pt0by5_r2_c1;
	Scatter2DPtr  _hist_rjs_pt0by5_r2_c8;
	Scatter2DPtr  _hist_rjs_pt10by5_r2_c1;
	Scatter2DPtr  _hist_rjs_pt10by5_r2_c8;

	Scatter2DPtr  _hist_rjs_pt0by5_r3_c1;
	Scatter2DPtr  _hist_rjs_pt0by5_r3_c8;
	Scatter2DPtr  _hist_rjs_pt10by5_r3_c1;
	Scatter2DPtr  _hist_rjs_pt10by5_r3_c8;
	
	Histo1DPtr _histdummypt0r2c1;
	Histo1DPtr _histdummypt0r2c8;
	Histo1DPtr _histdummypt0r3c1;
	Histo1DPtr _histdummypt0r3c8;
	
	Histo1DPtr _histdummypt10r2c1;
	Histo1DPtr _histdummypt10r2c8;
	Histo1DPtr _histdummypt10r3c1;
	Histo1DPtr _histdummypt10r3c8;
	
	Histo1DPtr _histdummypt5r2c1;
	Histo1DPtr _histdummypt5r2c8;
	Histo1DPtr _histdummypt5r3c1;
	Histo1DPtr _histdummypt5r3c8;
	
	Histo1DPtr _histdummypt5r2c1v2;
	Histo1DPtr _histdummypt5r2c8v2;
	Histo1DPtr _histdummypt5r3c1v2;
	Histo1DPtr _histdummypt5r3c8v2;
///////////////////////////////////////////////////////////////////	
	Scatter2DPtr _hist_RCP_r2_pt5_c1per8;
	Scatter2DPtr _hist_RCP_r2_pt5_c3per8;
	Scatter2DPtr _hist_RCP_r2_pt5_c5per8;
	
	Scatter2DPtr _hist_RCP_r3_pt5_c1per8;
	Scatter2DPtr _hist_RCP_r3_pt5_c3per8;
	Scatter2DPtr _hist_RCP_r3_pt5_c5per8;
	
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
	Scatter2DPtr _hist_RCP_vs_Np_r2_pt0;
	Scatter2DPtr _hist_RCP_vs_Np_r2_pt5;
	Scatter2DPtr _hist_RCP_vs_Np_r2_pt10;
	
	Scatter2DPtr _hist_RCP_vs_Np_r3_pt0;
	Scatter2DPtr _hist_RCP_vs_Np_r3_pt5;
	Scatter2DPtr _hist_RCP_vs_Np_r3_pt10;
////////////////////////////////////////////////////////////////////	

	Scatter2DPtr _hist_RCP_r2_pt0_c1per8;
	Scatter2DPtr _hist_RCP_r2_pt0_c3per8;
	Scatter2DPtr _hist_RCP_r2_pt0_c5per8;
	
	Scatter2DPtr _hist_RCP_r3_pt0_c1per8;
	Scatter2DPtr _hist_RCP_r3_pt0_c3per8;
	Scatter2DPtr _hist_RCP_r3_pt0_c5per8;
	
	Histo1DPtr _histdummyRCPr2pt0c1;
	Histo1DPtr _histdummyRCPr2pt0c3;
	Histo1DPtr _histdummyRCPr2pt0c5;
	
	Histo1DPtr _histdummyRCPr3pt0c1;
	Histo1DPtr _histdummyRCPr3pt0c3;
	Histo1DPtr _histdummyRCPr3pt0c5;
	
	Histo1DPtr _histdummyRCPr2pt0c8;
	Histo1DPtr _histdummyRCPr3pt0c8;
	Histo1DPtr _histdummyRCPr2pt0c8v2;
	Histo1DPtr _histdummyRCPr3pt0c8v2;
	Histo1DPtr _histdummyRCPr2pt0c8v3;
	Histo1DPtr _histdummyRCPr3pt0c8v3;

/////////////////////////////////////////////////////////////////////

	Scatter2DPtr _hist_RCP_r2_pt10_c1per8;
	Scatter2DPtr _hist_RCP_r2_pt10_c3per8;
	Scatter2DPtr _hist_RCP_r2_pt10_c5per8;
	
	Scatter2DPtr _hist_RCP_r3_pt10_c1per8;
	Scatter2DPtr _hist_RCP_r3_pt10_c3per8;
	Scatter2DPtr _hist_RCP_r3_pt10_c5per8;
	
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

	Histo1DPtr _histdummy_ratio_cjs_r2_c1;
	Histo1DPtr _histdummy_ratio_cjs_r2_c8;
	Histo1DPtr _histdummy_ratio_cjs_r3_c1;
	Histo1DPtr _histdummy_ratio_cjs_r3_c8;
	
	Scatter2DPtr _hist_ratio_cjs_r2by3_c1;
	Scatter2DPtr _hist_ratio_cjs_r2by3_c8;
	
	//}
	//decalring the variables for counting number of colliosns and target and porjectile particiants
	//{
	double _ncollc1 = 1500.5;
	double _ncollc3 = 738.8;
	double _ncollc5 = 245.6;
	double _ncollc8 = 45.9;
	
	double _npartc1 = 356.0;
	double _npartc3 = 223.0;
	double _npartc5 = 107.2;
	//}
	//divide function created to divide histograms 
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
