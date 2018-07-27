// -*- C++ -*-
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/CMSToolsHI.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include "Rivet/Projections/HeavyIonInfo.hh"


namespace Rivet
{



/// @brief Add a short analysis description here
	class CMS_2014_I1299142 : public HeavyIonAnalysis
	{
		public:

			/// Constructor
			CMS_2014_I1299142():HeavyIonAnalysis("CMS_2014_I1299142") {}


			/// @name Analysis methods
			//@{

			/// Book histograms and initialise projections before the run
			void init()
			{
				
				
				// Load MC pp reference data from .yoda file passed with -p flag
				//   If unable to find data use pp reference
				//   data supplied by paper
				try
				{
					_h_pp_xe_100_300 = getHisto1D("h_pp_xe_100_300");
					
				}
				catch (const exception &e)
				{
					_h_pp_xe_100_300 = bookHisto1D(1,1,1);
					_fillRef(_h_pp_xe_100_300, 1,1,1);
				}
				_h_pp_xe_100_300->setPath("/CMS_2014_I1299142/TMP/h_pp_xe_100_300");
				
				try
				{
					_h_pp_xe_100_120 = getHisto1D("h_pp_xe_100_120");
				}
				catch  (const exception &e)
				{
					_h_pp_xe_100_120 = bookHisto1D(2,1,1);
					_fillRef(_h_pp_xe_100_120, 2,1,1);
				
				}
				_h_pp_xe_100_120->setPath("/CMS_2014_I1299142/TMP/h_pp_xe_100_120");
				
				try
				{
					_h_pp_xe_120_150 = getHisto1D("h_pp_xe_120_150");
					
				}
				catch  (const exception &e)
				{
					_h_pp_xe_120_150 = bookHisto1D(3,1,1);
					_fillRef(_h_pp_xe_120_150, 3,1,1);
				
				}
				_h_pp_xe_120_150->setPath("/CMS_2014_I1299142/TMP/h_pp_xe_120_150");
				
				try
				{
					_h_pp_xe_150_300 = getHisto1D("h_pp_xe_150_300");
				}
				catch  (const exception &e)
				{
					_h_pp_xe_150_300 = bookHisto1D(4,1,1);
					_fillRef(_h_pp_xe_150_300, 4,1,1);
				
				}
				_h_pp_xe_150_300->setPath("/CMS_2014_I1299142/TMP/h_pp_xe_150_300");
				
				try
				{
					_h_pp_pt_100_300 = getHisto1D("h_pp_pt_100_300");
				}
				catch  (const exception &e)
				{
					_h_pp_pt_100_300 = bookHisto1D(5,1,1);
					_fillRef(_h_pp_pt_100_300, 5,1,1);
				}
				_h_pp_pt_100_300->setPath("/CMS_2014_I1299142/TMP/h_pp_pt_100_300");
				
				try
				{
					_h_pp_pt_100_120 = getHisto1D("h_pp_pt_100_120");
				}
				catch  (const exception &e)
				{
					_h_pp_pt_100_120 = bookHisto1D(6,1,1);
					_fillRef(_h_pp_pt_100_120, 6,1,1);
				}
				_h_pp_pt_100_120->setPath("/CMS_2014_I1299142/TMP/h_pp_pt_100_120");
				
				try
				{
					_h_pp_pt_120_150 = getHisto1D("h_pp_pt_120_150");
				}
				catch  (const exception &e)
				{
					_h_pp_pt_120_150 = bookHisto1D(7,1,1);
					_fillRef(_h_pp_pt_120_150, 7,1,1);
				}
				_h_pp_pt_120_150->setPath("/CMS_2014_I1299142/TMP/h_pp_pt_120_150");
				
				try
				{
					_h_pp_pt_150_300 = getHisto1D("h_pp_pt_150_300");
				}
				catch  (const exception &e)
				{
					_h_pp_pt_150_300 = bookHisto1D(8,1,1);
					_fillRef(_h_pp_pt_150_300, 8,1,1);
				}
				_h_pp_pt_150_300->setPath("/CMS_2014_I1299142/TMP/h_pp_xe_150_300");

				

				addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod");


				ChargedFinalState cfs(Cuts::abseta < 2.1 && Cuts::pt > 1.0 * GeV);
				declare(cfs, "CFS");

				CMSToolsHI bgf(cfs,CMSToolsHI::CMS_Reflection);
				declare(bgf,"BGF");



				FastJets jets(cfs, FastJets::ANTIKT, 0.3, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
				declare(jets, "jets");


				// Book Histograms
				for (int i =0; i < 5; i++)
				{
					for (int j=0; j <4; j++) count[j][i]=0.;

					string str = "d01-x01-y03";
					str[6] = '1' + i;
					_h_f1_L[i] = bookScatter2D(str,10,0.,5.);
					str[10]--;
					_h_f1_U[i] = bookHisto1D(str,10,0.,5.);
					str[10]--;
					_h_f1_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(1,1+i,1));

					str = "d05-x01-y03";
					str[6] = '1' + i;
					_h_f5_L[i] = bookScatter2D(str);
					str[10]--;
					_h_f5_U[i] = bookHisto1D(str,refData(5,1+i,2));
					str[10]--;
					_h_f5_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(5,1+i,1));

					if (i < 4)
					{
						str = "d02-x01-y03";
						str[6] = '1' + i;
						_h_f2_L[i] = bookScatter2D(str,10,0.,5.);
						str[10]--;
						_h_f2_U[i] = bookHisto1D(str);
						str[10]--;
						_h_f2_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(2,1+i,1));

						str = "d03-x01-y03";
						str[6] = '1' + i;
						_h_f3_L[i] = bookScatter2D(str,10,0.,5.);
						str[10]--;
						_h_f3_U[i] = bookHisto1D(str);
						str[10]--;
						_h_f3_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(3,1+i,1));

						str = "d04-x01-y03";
						str[6] = '1' + i;
						_h_f4_L[i] = bookScatter2D(str,10,0.,5.);
						str[10]--;
						_h_f4_U[i] = bookHisto1D(str);
						str[10]--;
						_h_f4_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(4,1+i,1));

						str = "d06-x01-y03";
						str[6] = '1' + i;
						_h_f6_L[i] = bookScatter2D(str);
						str[10]--;
						_h_f6_U[i] = bookHisto1D(str,refData(6,1+i,2));
						str[10]--;
						_h_f6_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(6,1+i,1));

						str = "d07-x01-y03";
						str[6] = '1' + i;
						_h_f7_L[i] = bookScatter2D(str);
						str[10]--;
						_h_f7_U[i] = bookHisto1D(str,refData(7,1+i,2));
						str[10]--;
						_h_f7_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(7,1+i,1));

						str = "d08-x01-y03";
						str[6] = '1' + i;
						_h_f8_L[i] = bookScatter2D(str);
						str[10]--;
						_h_f8_U[i] = bookHisto1D(str,refData(8,1+i,2));
						str[10]--;
						_h_f8_U_BG[i]= bookHisto1D("TMP/bg"+ str,refData(8,1+i,1));

					}
				}

			}

			void analyze(const Event& event)
			{




				const double c = centrality(event, "IPMethod");

				if ((c < 0.) || (c > 100.)) vetoEvent; // ignore calibration events

				
				int bin = getBin(c);
				const double weight = event.weight();

				const FinalState & cfs = apply<ChargedFinalState>(event,"CFS");
				Particles particles = cfs.particles();
				const FastJets & jets = apply<FastJets>(event,"jets");
				Jets allJets = jets.jetsByPt(100.0*GeV); 
				const CMSToolsHI &bgf = apply<CMSToolsHI>(event,"BGF");


				foreach (const Jet& jet, allJets)
				{
					Particles pbg = bgf.CMSReflectionParticles(jet,0.3);
					double jetPt = jet.pt() - bgf.CMSReflectionFourMomentum(jet,0.3).pt();
					
					if (jetPt > 100 * GeV && abs(jet.eta())>0.3  && abs(jet.eta())<2.0 )
					{
						
						if (jetPt < 300.*GeV)
						{
							count[0][bin]+=weight;
							
							// Get ALL particles within R of 0.3 of jet centroid, fill histo
							foreach (const Particle& p, getConeParticles(particles,jet.eta(),jet.phi()))
							{
								_h_f1_U[bin]->fill(xe(jet,p,jetPt),weight);
								_h_f5_U[bin]->fill(p.pt(),weight);
							}
							// Get BG Particles within R of 0.3 of eta reflected centroid
							// CMSReflection(Jet& jetOfInterest, double R_Param)
							foreach (const Particle& p, pbg)
							{
								_h_f1_U_BG[bin]->fill(xeBG(jet,p,jetPt),weight);
								_h_f5_U_BG[bin]->fill(p.pt(),weight);
							}
						
							if (jetPt < 120.*GeV)
							{
								
								if (bin > 0) bin--;
								count[1][bin]+=weight;
								foreach (const Particle& p, getConeParticles(particles,jet.eta(),jet.phi()))
								{
									_h_f2_U[bin]->fill(xe(jet,p,jetPt),weight);
									_h_f6_U[bin]->fill(p.pt(),weight);
								}

								// Get BG Particles within R of 0.3 of eta reflected centroid
								foreach (const Particle& p, pbg)
								{
									_h_f2_U_BG[bin]->fill(xeBG(jet,p,jetPt),weight);
									_h_f6_U_BG[bin]->fill(p.pt(),weight);
								}
							}
							else if (jetPt < 150.*GeV)
							{   
								
								if (bin > 0) bin--;
								count[2][bin]+=weight;
								foreach (const Particle& p, getConeParticles(particles,jet.eta(),jet.phi()))
								{
									_h_f3_U[bin]->fill(xe(jet,p,jetPt),weight);
									_h_f7_U[bin]->fill(p.pt(),weight);
								}

								// Get BG Particles within R of 0.3 of eta reflected centroid
								foreach (const Particle& p, pbg)
								{
									_h_f3_U_BG[bin]->fill(xeBG(jet,p,jetPt),weight);
									_h_f7_U_BG[bin]->fill(p.pt(),weight);
								}
							}
							else
							{
							
								if (bin > 0) bin--;
								count[3][bin]+=weight;
								foreach (const Particle& p, getConeParticles(particles,jet.eta(),jet.phi()))
								{
									_h_f4_U[bin]->fill(xe(jet,p,jetPt),weight);
									_h_f8_U[bin]->fill(p.pt(),weight);
								}

								// Get BG Particles within R of 0.3 of eta reflected centroid
								foreach (const Particle& p, pbg)
								{
									_h_f4_U_BG[bin]->fill(xeBG(jet,p,jetPt),weight);
									_h_f8_U_BG[bin]->fill(p.pt(),weight);
								}
							}
							
						}
					}
				}
				  

			}

			void finalize()
			{

				for (int i =0; i < 4; i++)
				{

					// Remove BG contribution
					*_h_f1_U[i]-=*_h_f1_U_BG[i];
					
					*_h_f2_U[i]-=*_h_f2_U_BG[i];
					*_h_f3_U[i]-=*_h_f3_U_BG[i];
					*_h_f4_U[i]-=*_h_f4_U_BG[i];
					*_h_f5_U[i]-=*_h_f5_U_BG[i];
					*_h_f6_U[i]-=*_h_f6_U_BG[i];
					*_h_f7_U[i]-=*_h_f7_U_BG[i];
					*_h_f8_U[i]-=*_h_f8_U_BG[i];

					// Scale to 1/Njet
					scale(_h_f1_U[i],1./count[0][i]);
					scale(_h_f2_U[i],1./count[1][i]);
					scale(_h_f3_U[i],1./count[2][i]);
					scale(_h_f4_U[i],1./count[3][i]);
					scale(_h_f5_U[i],1./count[0][i]);
					scale(_h_f6_U[i],1./count[1][i]);
					scale(_h_f7_U[i],1./count[2][i]);
					scale(_h_f8_U[i],1./count[3][i]);

					// Get relative data
					divide(_h_f1_U[i],_h_pp_xe_100_300,_h_f1_L[i]);
					divide(_h_f2_U[i],_h_pp_xe_100_120,_h_f2_L[i]);
					divide(_h_f3_U[i],_h_pp_xe_120_150,_h_f3_L[i]);
					divide(_h_f4_U[i],_h_pp_xe_150_300,_h_f4_L[i]);
					divide(_h_f5_U[i],_h_pp_pt_100_300,_h_f5_L[i]);
					divide(_h_f6_U[i],_h_pp_pt_100_120,_h_f6_L[i]);
					divide(_h_f7_U[i],_h_pp_pt_120_150,_h_f7_L[i]);
					divide(_h_f8_U[i],_h_pp_pt_150_300,_h_f8_L[i]);

				}
				*_h_f1_U[4]-=*_h_f1_U_BG[4];
				*_h_f5_U[4]-=*_h_f5_U_BG[4];
				scale(_h_f1_U[4],1./count[0][4]);
				scale(_h_f5_U[4],1./count[0][4]);
				divide(_h_f1_U[4],_h_pp_xe_100_300,_h_f1_L[4]);
				divide(_h_f5_U[4],_h_pp_pt_100_300,_h_f5_L[4]);

			}



		private:
			double z (const Jet& jet, const Particle& part, double adjPt)
			{
				return dot(jet.p3(),part.p3()) / (jet.p() * adjPt );
			}
			double xe (const Jet& jet, const Particle& part, double adjPt)
			{
				return log(1./z(jet,part, adjPt));
			}

			double xeBG (const Jet& jet, const Particle& part, double adjPt)
			{
				double _z = ((-jet.pz() * part.pz()) + (jet.px() * part.px()) + (jet.py() * part.py()))/ (jet.p() * adjPt);
				return log(1./_z);
			}
			Particles getConeParticles(const Particles& ps, double eta, double phi,double R = 0.3)
			{
				Particles conePart;
				foreach (const Particle &p, ps)
				{
					if (deltaR(p,eta,phi) < R) conePart.push_back(p);
				}
				return conePart;
			}


			int getBin(double cent)
			{
				if (cent < 10.) return 4;
				else if (cent < 30.) return 3;
				else if (cent < 50.) return 2;
				else if (cent < 70.) return 1;
				else return 0;
			}
			
//			void divide(Histo1DPtr h1, Scatter2DPtr h2, Scatter2DPtr out){
//				const string path = out->path();
//				*out = *h1 / *h2;
//				out->setPath(path);
//				cout << path << " scatter div\n";
//			}
			void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr out){
				const string path = out->path();
				*out = *h1 / *h2;
				out->setPath(path);
			}

			void _fillRef(Histo1DPtr & h, int d, int x, int y){
				const Scatter2D scat = refData(d,x,y);
				const std::vector<Point2D>&  ps = scat.points();
				foreach (const Point2D& p, ps){
					h->fill(p.x(), p.y());
				}
				
			}

			// histograms;
			// "_h_" + 'f' + (figure number) + '_'
			// centrality bin in [x],
			// if 4: 0 = 50-100, 1 - 30-50, 2 - 10-30, 3 - 0-10
			// if 5: 0 = 70-100, 1 - 50-70, 2 - 30-50, 3 - 10-30, 4 - 0-10
			Histo1DPtr _h_f1_U[5]; // upper half fig 1 (1/Njet)dNtrack/dXe jet pT 100-300 GeV/c
			Scatter2DPtr _h_f1_L[5]; // lower half fig 1 (PbPb/pp) numerator


			Histo1DPtr _h_f2_U[4]; // upper half fig 2 (1/Njet)dNtrack/dXe jet pT 100-120 GeV/c
			Histo1DPtr _h_f3_U[4]; // upper half fig 3 (1/Njet)dNtrack/dXe jet pT 120-150 GeV/c
			Histo1DPtr _h_f4_U[4]; // upper half fig 4 (1/Njet)dNtrack/dXe jet pT 150-300 GeV/c
			Histo1DPtr _h_f5_U[5]; // upper half fig 5 (1/Njet)dNtrack/dpT jet pT 100-300 GeV/c
			Histo1DPtr _h_f6_U[4]; // upper half fig 6 (1/Njet)dNtrack/dpT jet pT 100-120 GeV/c
			Histo1DPtr _h_f7_U[4]; // upper half fig 7 (1/Njet)dNtrack/dpT jet pT 120-150 GeV/c
			Histo1DPtr _h_f8_U[4]; // upper half fig 8 (1/Njet)dNtrack/dpT jet pT 150-300 GeV/c
			
			Histo1DPtr _h_f1_U_out[4]; // upper half fig 2 (1/Njet)dNtrack/dXe jet pT 100-120 GeV/c
			Histo1DPtr _h_f2_U_out[4]; // upper half fig 2 (1/Njet)dNtrack/dXe jet pT 100-120 GeV/c
			Histo1DPtr _h_f3_U_out[4]; // upper half fig 3 (1/Njet)dNtrack/dXe jet pT 120-150 GeV/c
			Histo1DPtr _h_f4_U_out[4]; // upper half fig 4 (1/Njet)dNtrack/dXe jet pT 150-300 GeV/c
			Histo1DPtr _h_f5_U_out[5]; // upper half fig 5 (1/Njet)dNtrack/dpT jet pT 100-300 GeV/c
			Histo1DPtr _h_f6_U_out[4]; // upper half fig 6 (1/Njet)dNtrack/dpT jet pT 100-120 GeV/c
			Histo1DPtr _h_f7_U_out[4]; // upper half fig 7 (1/Njet)dNtrack/dpT jet pT 120-150 GeV/c
			Histo1DPtr _h_f8_U_out[4]; // upper half fig 8 (1/Njet)dNtrack/dpT jet pT 150-300 GeV/c
			
			Histo1DPtr _h_f1_U_BG[5];
			Histo1DPtr _h_f2_U_BG[5];
			Histo1DPtr _h_f3_U_BG[5];
			Histo1DPtr _h_f4_U_BG[5];
			Histo1DPtr _h_f5_U_BG[5];
			Histo1DPtr _h_f6_U_BG[5];
			Histo1DPtr _h_f7_U_BG[5];
			Histo1DPtr _h_f8_U_BG[5];

			Scatter2DPtr _h_f2_L[4]; // lower half fig 2 (PbPb/pp)
			Scatter2DPtr _h_f3_L[4]; // lower half fig 3 (PbPb/pp)
			Scatter2DPtr _h_f4_L[4]; // lower half fig 4 (PbPb/pp)
			Scatter2DPtr _h_f5_L[5]; // lower half fig 5 (PbPb/pp)
			Scatter2DPtr _h_f6_L[4]; // lower half fig 6 (PbPb/pp)
			Scatter2DPtr _h_f7_L[4]; // lower half fig 7 (PbPb/pp)
			Scatter2DPtr _h_f8_L[4]; // lower half fig 8 (PbPb/pp)


			double count[4][5];

			// pp reference data
			Histo1DPtr _h_pp_xe_100_300;
			Histo1DPtr _h_pp_xe_100_120;
			Histo1DPtr _h_pp_xe_120_150;
			Histo1DPtr _h_pp_xe_150_300;
			Histo1DPtr _h_pp_pt_100_300;
			Histo1DPtr _h_pp_pt_100_120;
			Histo1DPtr _h_pp_pt_120_150;
			Histo1DPtr _h_pp_pt_150_300;

		

	};

	// The hook for the plugin system
	DECLARE_RIVET_PLUGIN(CMS_2014_I1299142);


}
