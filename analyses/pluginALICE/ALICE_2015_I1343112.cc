// -*- C++ -*-
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/Projections/ALICEToolsHI.hh"
#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include <fstream>

namespace Rivet
{


	/// @brief Add a short analysis description here
	class ALICE_2015_I1343112 : public HeavyIonAnalysis
	{
		public:

			/// Constructor
			ALICE_2015_I1343112(): HeavyIonAnalysis ("ALICE_2015_I1343112") {}

			void init()
			{
				HeavyIonAnalysis::init();
				
				//Initialise and register projections
				//takes 50 events and uses them for creating the centrality binning
				addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "ImpactParameterMethod");
				
				//loads the pp reference histograms for the RAA calculations
				// Tries to load from preloaded file, if unable loads pp data used
				// in paper
				try
				{
					_hist_raa_pt5_c1_pp = getHisto1D("d06-x01-y01");

				}
				catch(const exception &e)
				{
					_hist_raa_pt5_c1_pp = bookHisto1D(6, 1, 1);
					_fillRef(_hist_raa_pt5_c1_pp, 6, 1, 1);
				}
				_hist_raa_pt5_c1_pp->setPath("/ALICE_2015_I1343112/TMP/d06-x01-y01");
				try
				{
					_hist_raa_pt5_c3_pp = getHisto1D("d07-x01-y01");
				}
				catch(const exception &e)
				{
					_hist_raa_pt5_c3_pp = bookHisto1D(7, 1, 1);
					_fillRef(_hist_raa_pt5_c3_pp, 7, 1, 1);
				}
				_hist_raa_pt5_c3_pp->setPath("/ALICE_2015_I1343112/TMP/d07-x01-y01");
				//cut of jets with pt 0.15 GeV and higher
				const ChargedFinalState cfs(Cuts::pT > 0.15*GeV);
				declare (cfs, "CFS");
				//sets cone radius parameter to 0.2
				FastJets fjr2(cfs, FastJets::ANTIKT, 0.2);
				//sets up to find area
				fjr2.useJetArea(new fastjet::AreaDefinition(fastjet::active_area, fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
				declare(fjr2, "jets");
				ALICEToolsHI athr2(cfs, fjr2);
				declare(athr2, "ATHr2");
				
				// Book histograms
				//{
				//booking histograms of Jet Spectra with pT > 5GeV in R=0.2 at 0-10% and 10-30% centrality
				_hist_js_pt5_r2_c1 = bookHisto1D(3, 1, 1);
				_hist_js_pt5_r2_c3 = bookHisto1D(3, 1, 2);
				//booking scatter plots of Ratio of Jet Spectra for pT (0/5, 3/5, 7/5, 10/5)GeV at 0-10% centrality in R=0.2
				_hist_rjs_r2_c1_pt0by5 = bookScatter2D(5, 1, 1);
				_hist_rjs_r2_c1_pt3by5 = bookScatter2D(5, 1, 2);
				_hist_rjs_r2_c1_pt7by5 = bookScatter2D(5, 1, 3);
				_hist_rjs_r2_c1_pt10by5 = bookScatter2D(5, 1, 4);
				//booking dummy histograms of jet spectra with pt > (0, 3, 5, 7, 10)GeV at 0-10% centrality in R=0.2
				_histdummypt0 = bookHisto1D("TMP/dummyhistpt0",  refData(5, 1, 1));
				_histdummypt3 = bookHisto1D("TMP/dummyhistpt3",  refData(5, 1, 2));
				_histdummypt7 = bookHisto1D("TMP/dummyhistpt7",  refData(5, 1, 3));
				_histdummypt10 = bookHisto1D("TMP/dummyhistpt10",  refData(5, 1, 4));

				_histdummypt5 = bookHisto1D("TMP/dummyhistpt5",  refData(5, 1, 1));
				_histdummypt5v2 = bookHisto1D("TMP/dummyhistpt5v2",  refData(5, 1, 2));
				_histdummypt5v3 = bookHisto1D("TMP/dummyhistpt5v3",  refData(5, 1, 3));
				_histdummypt5v4 = bookHisto1D("TMP/dummyhistpt5v4",  refData(5, 1, 4));
				//there is mulitple version of the above graphs due to a binning issue if you divide 2 histograms with differnet refrence data
				//booking scatter plots of nuclear mod factor RAA with pt > 5GeV at 0-10% and 10-30% in R=0.2
				_hist_raa_pt5_c1 = bookScatter2D(6, 1, 1);
				_hist_raa_pt5_c3 = bookScatter2D(7, 1, 1);
				//booking dummy histograms of jet spectra with pt > 5GeV at 0-10% centrality in R=0.2
				_histdumraapt5c1 = bookHisto1D("TMP/dumhistraa1", refData(6, 1, 1));
				_histdumraapt5c3 = bookHisto1D("TMP/dumhistraa2", refData(7, 1, 1));
				//}
			}


			/// Perform the per-event analysis
			void analyze(const Event& event)
			{
				//vetos any event that is not within 0-30% centralties
				const double c = centrality(event, "ImpactParameterMethod");
			
				if ((c < 0.)||(c>30.))
				{
					vetoEvent;
				}
				//decalring weight for event
				const double weight = event.weight();
				//decalring fast jet constants
				const FastJets& fastjetr2 = apply<FastJets>(event, "jets");
				//setting up the jets
				const Jets& jets2 = fastjetr2.jetsByPt(.15*GeV);
				//decalres variabels for fining area and rho
				const auto seqr2 = fastjetr2.clusterSeqArea();
				const ALICEToolsHI &athr2 = apply<ALICEToolsHI>(event, "ATHr2");
				//goes through all the events in the 0-10% centraltiy
				if (c <= 10.)
				{
					_eventcounter1 += weight; //counts the number of events within this centrality
					//runs through all of the jets in R=0.2
					foreach(const Jet& jets, jets2)
					{
						//gets area and rho0 for the background subtraction equation
						double area = seqr2->area(jets);
						double rho0 = athr2.RhoLocal(jets.phi());
						double jetpt = jets.pT() - area*rho0;
						//creates varaible to find the constituents with a certain pt range
						const Particles ps = jets.constituents();
						//fills histograms with jets with pt > 0GeV
						if(any(ps, PtGtr(0*GeV)))
						{
							_histdummypt0->fill(jetpt, weight);
						}
						//fills histograms with jets with pt > 3GeV
						if(any(ps, PtGtr(3.*GeV)))
						{
							_histdummypt3->fill(jetpt, weight);
						}
						//fills histograms with jets with pt > 5GeV
						if(any(ps, PtGtr(5.*GeV)))
						{
							_hist_js_pt5_r2_c1->fill(jetpt, weight);
							_histdummypt5->fill(jetpt, weight);
							_histdummypt5v2->fill(jetpt, weight);
							_histdummypt5v3->fill(jetpt, weight);
							_histdummypt5v4->fill(jetpt, weight);
							_histdumraapt5c1->fill(jetpt, weight);
						}
						//fills histograms with jets with pt > 7GeV
						if(any(ps, PtGtr(7.*GeV)))
						{
							_histdummypt7->fill(jetpt, weight);
						}
						//fills histograms with jets with pt > 10GeV
						if(any(ps, PtGtr(10.*GeV)))
						{
							_histdummypt10->fill(jetpt, weight);
						}
					}
				}
				//goes through all the events in the 10-30% centraltiy
				else
				{
					_eventcounter3 += weight; //counts the number of events within this centrality
					//runs through all of the jets in R=0.3
					foreach(const Jet& jets, jets2)
					{
						//gets area and rho0 for the background subtraction equation
						double area = seqr2->area(jets);
						double rho0 = athr2.RhoLocal(jets.phi());
						double jetpt = jets.pT() - area*rho0;
						//creates varaible to find the constituents with a certain pt range
						const Particles ps = jets.constituents();
						//fills histograms with jets with pt > 5GeV
						if(any(ps, PtGtr(5.*GeV)))
						{
							_hist_js_pt5_r2_c3->fill(jetpt, weight);
							_histdumraapt5c3->fill(jetpt, weight);
						}
					}
				}


			}


			/// Normalise histograms etc., after the run
			void finalize()
			{
				//scaling the histograms
				scale(_hist_js_pt5_r2_c1, (1.)/(_eventcounter1*_ncollc1));
				scale(_hist_js_pt5_r2_c3, (1.)/(_eventcounter3*_ncollc3));
				scale(_histdumraapt5c1, (1.)/(_eventcounter1*_ncollc1));
				scale(_histdumraapt5c3, (1.)/(_eventcounter3*_ncollc3));
				//dividing histograms and filling the scatter plots
				divide(_histdummypt0, _histdummypt5, _hist_rjs_r2_c1_pt0by5);
				divide(_histdummypt3, _histdummypt5v2, _hist_rjs_r2_c1_pt3by5);
				divide(_histdummypt7, _histdummypt5v3, _hist_rjs_r2_c1_pt7by5);
				divide(_histdummypt10, _histdummypt5v4, _hist_rjs_r2_c1_pt10by5);
				divide(_histdumraapt5c1, _hist_raa_pt5_c1_pp, _hist_raa_pt5_c1);
				divide(_histdumraapt5c3, _hist_raa_pt5_c3_pp, _hist_raa_pt5_c3);
			}

			//decalring histograms and scatter plots
			Histo1DPtr _hist_js_pt5_r2_c1;
			Histo1DPtr _hist_js_pt5_r2_c3;
			Histo1DPtr _histdummypt0;
			Histo1DPtr _histdummypt3;
			Histo1DPtr _histdummypt7;
			Histo1DPtr _histdummypt10;

			Histo1DPtr _histdummypt5;
			Histo1DPtr _histdummypt5v2;
			Histo1DPtr _histdummypt5v3;
			Histo1DPtr _histdummypt5v4;

			Scatter2DPtr _hist_rjs_r2_c1_pt0by5;
			Scatter2DPtr _hist_rjs_r2_c1_pt3by5;
			Scatter2DPtr _hist_rjs_r2_c1_pt7by5;
			Scatter2DPtr _hist_rjs_r2_c1_pt10by5;

			Histo1DPtr _hist_raa_pt5_c1_pp;
			Histo1DPtr _hist_raa_pt5_c3_pp;

			Histo1DPtr _histdumraapt5c1;
			Histo1DPtr _histdumraapt5c3;

			Scatter2DPtr _hist_raa_pt5_c1;
			Scatter2DPtr _hist_raa_pt5_c3;
			//decalring the event counters with in each centrality
			double _eventcounter1 = 0;
			double _eventcounter3 = 0;
			//declaring the avg  number of collions counter for differnet centralities
			// numbers taken from paper
			double _ncollc1 = 1501;
			double _ncollc3 = 743;

			//declaring the divide function
			void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr out)
			{
				const string path = out->path();
				*out = *h1 / *h2;
				out->setPath(path);
			}
			
			//Fill a histogram based on reference data from paper
			void _fillRef(Histo1DPtr & h, int d, int x, int y)
			{
				const Scatter2D scat = refData(d,x,y);
				const std::vector<Point2D>&  ps = scat.points();
				foreach (const Point2D& p, ps)
				{
					h->fill(p.x(), p.y());
				}

			}

	};


	// The hook for the plugin system
	DECLARE_RIVET_PLUGIN(ALICE_2015_I1343112);


}
