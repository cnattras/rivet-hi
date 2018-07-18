// -*- C++ -*-
#ifndef RIVET_ATLASToolsHI_HH
#define RIVET_ATLASToolsHI_HH

#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Projections/ATLASCalorimeterTowers.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include <algorithm>
#include <cmath>

namespace Rivet
{				


	/// @brief Project Out the Background of an event
	///
	///
	/// @todo
	class ATLASToolsHI : public Projection
	{
		public:


			/// @name Constructors and destructors.
			//@{

			/// Constructor
			///
			
			ATLASToolsHI(const ATLASCalorimeterTowers &act, const Cut& c = (Cuts::eta < 2.1 && Cuts::eta > -2.1), double towerWidthRap = 0.1, double towerWidthPhi = 0.1)
			{
				setName("ATLASToolsHI");
				addProjection(act, "ACT");
				
				
				ChargedFinalState tracks(c && Cuts::pt > 0.15);
				declare(tracks,"TRACKS");
				
				double minRap = act.RapidityMin();
				double maxRap = act.RapidityMax();
				_nPhi = act.PhiTowers();
				_wPhi = act.PhiWidth();
				_nRap = act.RapTowers();
				_wRap = act.RapWidth();
				_minRap = minRap;
				_maxRap = maxRap;
				_nRap2 = _nRap /2;
				
			//	_towerParticles = std::vector<fastjet::PseudoJet>(_nPhi * _nRap);
				_towerArea = _wPhi * _wRap;
				
				
				ATLASCalorimeterTowers fcal(Cuts::abseta > 3.2 && Cuts::abseta < 4.9, 3.2, 4.9, 0.2, 0.2);
				declare(fcal,"FCAL");
				
				ParticleVn fcalVn(fcal,2,_nRap);
				declare(fcalVn,"FCALVN");
				
				
				FastJets jets(act, FastJets::ANTIKT, 0.2, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
				declare(jets, "CJETS");
				
				FastJets trackjets(tracks, FastJets::ANTIKT, 0.4, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
				declare(trackjets, "TJETS");
				//_calorimeterJets = Jets(1);
				
			}

			const Particles& calorimeterTowers() const { return _calorimeterTowers; }
			const Jets& calorimeterJets() const { cout << _calorimeterJets.size() << " return\n"; return _calorimeterJets; }

			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(ATLASToolsHI);

			//@}


			/// @name
			//@{




		protected:

			/// Apply the projection to the event.
			void project(const Event& e)
			{
				_calorimeterJets.clear();
				_calorimeterTowers = apply<ATLASCalorimeterTowers>(e, "ACT").particles();
				//const ParticleVn &midVn = apply<ParticleVn>(e, "MIDVN");
				
				// Get Calorimeter Jets
				Jets calorimeterJets(apply<FastJets>(e, "CJETS").jets());
				Jets caloJetsAdj(calorimeterJets);
				cout << "size 1 " << calorimeterJets.size() << endl;
				
				// Find seed jets 
				_tagSeedJets1(caloJetsAdj);
				
				// Get v2 from FCAL
				const ParticleVn &fcalVn = apply<ParticleVn>(e, "FCALVN");
				_v2 = fcalVn.GetVn(2);
				_psi2 = fcalVn.EventPlaneN(2);
				cout << "size 2 " << calorimeterJets.size() << endl;
				// Calculate rhoi(eta) for mid rapidities;
				_calculateRhoi1(caloJetsAdj);
				
				// Adjust Jet 4 vectors to account for initial background
				caloJetsAdj = calorimeterJets;
				foreach (const Jet& j, calorimeterJets) {
				
					if (j.Et() == 0) {
					  cout << "Jet: " << j.Et() << ' ' << j.eta() << ' ' << j.phi() << endl;
					  foreach (const Particle& p, j.particles())
								cout << "part: " << p.pid() << p.px() << ' ' << p.py() << ' ' << p.pz() << endl;
						  }
				  }
				  cout << "x1\n";
				_subtractInterimUE(caloJetsAdj);
				
				
				// Find New Seed Jets including Track Jets
				const Jets &trackJets = apply<FastJets>(e, "TJETS").jets(Cuts::pt > 10 * GeV); // per 1208.1967v2
				Jets cSeeds, tSeeds;
				_getSeedJets2(caloJetsAdj, trackJets, cSeeds, tSeeds);
				
				
				
				
				// Get new rhoi and V2
				_calculateRhoi2( _calorimeterTowers, cSeeds, tSeeds);
				_calcV2Psi2( _calorimeterTowers, cSeeds, tSeeds );
				cout << "size 3 " << calorimeterJets.size() << endl;
				_subtractUE();
				
				foreach (const Jet& j, calorimeterJets) {
				
					if (j.Et() == 0) {
					  cout << "Jet: " << j.Et() << ' ' << j.eta() << ' ' << j.phi() << endl;
					  foreach (const Particle& p, j.particles())
								cout << "part: " << p.pid() << p.px() << ' ' << p.py() << ' ' << p.pz() << endl;
						  }
				  }
				  cout << "x2\n";
				
				// Subtrack UE from towers, create final jets
				
				cout << "size 4 " << calorimeterJets.size() << endl;
				_finalizeJets(calorimeterJets);
				cout << "size 5 " << calorimeterJets.size() << endl;
				
			}

			/// Compare projections (only difference is in UFS definition)
			int compare(const Projection& p) const
			{
				return -99999;
			}

	private:
			void _tagSeedJets1(const Jets &js){
				_isSeed = std::vector<bool>(js.size(),false);
				for (int i =0; i < (signed)js.size(); ++i){
					const Particles& ps = js[i].particles();
					bool has3gevTower = false;
					int nPart = ps.size();
					double sumEt = 0.;
					double maxEt = 0.;
					foreach (const Particle& p, ps){
						if (p.Et() > 3 * GeV) has3gevTower = true;
						sumEt += p.Et();
						if (p.Et() > maxEt) maxEt = p.Et();
					}
					double avgEt = sumEt / nPart;
					if (has3gevTower && maxEt/avgEt > 4.0) _isSeed[i] = true;
				}
				
			}
			void _calculateRhoi1(const Jets& js)
			{
				_rhoi = std::vector<double>(_nRap, 0.0);
				_rhoiEntries = std::vector<int>(_nRap, 0);
				int nJet = js.size();
				
				for(int i = 0; i < nJet; ++i)
				{
					if (!_isSeed[i]){
						const Particles& ps = js[i].constituents();
						foreach (const Particle& p, ps)
						{
							
							int bin = _rapBin(p.eta());
							_rhoi[bin] += p.Et();
							++_rhoiEntries[bin];

						}
						
					}
				}
				
				for (int i = 0; i < _nRap; ++i){
					if (_rhoiEntries[i]) _rhoi[i] /= _rhoiEntries[i] * _wRap * _wPhi;
					
				}
			}
			
			void _subtractInterimUE(Jets &js){
				foreach (Jet& j, js) {
					FourMomentum mom(0,0,0,0);
					foreach (const Particle &p, j.particles()){
						mom += _adjMom(p);
					}
					j.setState(mom, j.particles());
				}
			}
			FourMomentum _adjMom(const Particle& p){
				double theta = p.theta();
				double phi = p.phi();
				double Et = p.Et();
				double Eadj = 1e-20; 
				double padj = 1e-20; // needed to 
				double EtUE = _towerArea * _rhoi[_rapBin(p.eta())] * (1 + 2 * _v2 * cos(2 * (phi - _psi2)));
				if (Et-EtUE > 0.0) Eadj = (Et - EtUE) / sin(theta);
				if (Eadj > padj) padj = Eadj;
				
				return FourMomentum(Eadj,padj*cos(phi)*sin(theta),padj*sin(phi)*sin(theta),padj*cos(theta));
				
			}
			
			void _getSeedJets2(const Jets& caloJets, const Jets& trackJets, Jets& cSeeds, Jets& tSeeds){
				//cout << "caloJet Et: ";
				foreach (const Jet& j, caloJets){
					//cout << j.Et() << '\t';
					if (j.Et() > 25 * GeV) {
						cSeeds.push_back(j);
						//cout << "\nCAL seed " << j.Et() << ' ' << j.pt() << ' '<< j.phi() << ' ' << j.eta() << endl;
					}
				}
				foreach (const Jet& j, trackJets){
					if(any(j.particles(),PtGtr(4*GeV))){
						
						tSeeds.push_back(j);
						//cout << "TRK seed " << j.Et() << ' ' << j.pt() << ' '<< j.phi() << ' ' << j.eta() << endl;
					}
				}
			}
			void _calculateRhoi2(const Particles& towers, const Jets& cSeeds, const Jets& tSeeds){
				//cout << "towers\n";
//				foreach (const Particle &p, towers){
//					if (p.phi() == 0.) cout << p.eta() << '\t' << p.phi() << endl;
//				}
				_rhoi = std::vector<double>(_nRap, 0.0);
				_rhoiEntries = std::vector<int>(_nRap, 0);
				
				//cout << "\nrhoiEntries ";
				//foreach (const int &x, _rhoiEntries) { cout << x << ' '; }
				//cout << endl;
				
				//cout << "seeds: \n";
//				foreach (const Jet& j, seeds) {
//					cout << j.eta() << '\t' << j.phi() << '\n';
//				}
//				bool ks = false;
//				if (seeds.size() > 0) ks= true;
//				 (ks) cout << endl << "towers acc: ";
				foreach (const Particle &p, towers){
					if (!any(cSeeds,DeltaRLess(p,0.4)) && !any(tSeeds,DeltaRLess(p,0.4))) {
						//if (ks) cout << "exclude "<< p.eta() << '\t' << p.phi() << endl;
						int bin = _rapBin(p.eta());
						_rhoi[bin] += p.Et();
						++_rhoiEntries[bin];
					}
					
				}
				//cout << "\nrhoiEntries ";
				//foreach (const int &x, _rhoiEntries) { cout << x << ' '; if (x==0) exit(1);}
				cout << endl;
				for (int i = 0; i < _nRap; ++i){
					if (_rhoiEntries[i]) _rhoi[i] /= _rhoiEntries[i] * _wRap * _wPhi;
					//cout << _rhoi[i] << " rhoi " << i << endl;
				}
				
				
			}
			
			void _calcV2Psi2(const Particles& towers, const Jets& cSeeds, const Jets& tSeeds){
				double Qx = 0.0;
				double Qy = 0.0;
				foreach (const Particle& p, towers){
					if (!any(cSeeds,DeltaRLess(p,0.4)) && !any(tSeeds,DeltaRLess(p,0.4))) {
						Qx += p.Et() * cos(2 * p.phi());
						Qy += p.Et() * sin(2 * p.phi());
					}
					
				}
				cout << _psi2;
				_psi2 = 0.5 * atan2(Qy,Qx);
//				cout << " old psi2 new " << _psi2 << endl;
				
				double EtSum = 0.0;
				std::vector<double> EtPhi(_nPhi, 0.0);
//				cout << "\netas for v2 ";
//				int v2Cnt = 0;
				foreach (const Particle& p, towers){
					if (!any(cSeeds,DeltaRLess(p,0.4)) && !any(tSeeds,DeltaRLess(p,0.4)))
					{
						int bin = static_cast<int>(p.phi() /_wPhi);
						EtPhi[bin] += p.Et();
						EtSum += p.Et();
//						if (bin == 0) {cout << p.eta() << '\t'; ++v2Cnt; }
					}
				}
//				cout << "\n**************" << v2Cnt << " eta bands\n";
//				cout << endl;
				int cnt = 0;
				//foreach (const double& d, EtPhi) {cout << d << '\t'; ++cnt; if (cnt ==_nPhi/2) {cout << endl; cnt =0; }}
				//cout << "\n V2 in " << _v2;
				_v2 = 0.;
				for (int i = 0; i < _nPhi; ++i){
					double phi = (0.5 + static_cast<double>(i)) * _wPhi;
					//cout << EtPhi[i] << '\t';
					_v2 += EtPhi[i] * cos(2*(phi - _psi2));
				}
				_v2 /= EtSum;
				//cout << "\t v2 out " << _v2 << endl;
			}
			
			void _subtractUE(){
				foreach (Particle &p, _calorimeterTowers)
					p.setMomentum(_adjMom(p));
			}
			
			void _finalizeJets(const Jets& js){
				cout << endl << js.size() << " sze\n";
				foreach (const Jet &j, js){
					
					
					if (j.particles().size() !=0){ 
						Particles jPart;
						
						FourMomentum mom(0,0,0,0);
						foreach (const Particle &p, j.particles()){
							mom += _calorimeterTowers[p.pid()].momentum();
							jPart.push_back(_calorimeterTowers[p.pid()]);
						}
						_calorimeterJets.push_back(Jet(mom,Particles(jPart)));
					} else {
						_calorimeterJets.push_back(Jet(j));
					}
					//j.setState(mom);
					
				}
			}

			int _rapBin(double eta){
				if (eta < 0) return static_cast<int>((eta+_maxRap)/_wRap);
					else return _nRap2 + static_cast<int>((eta-_minRap)/_wRap);
			}


			double _minRap;
			double _maxRap;
			int _nRap;
			int _nRap2;
			int _nPhi;
			double _wRap;
			double _wPhi;
			
			
			std::vector<FourMomentum> _towerMom0;
			double _towerArea;
			std::vector<double> _rhoi; // rho i as a function of eta
			std::vector<int> _rhoiEntries;
			double _v2;
			double _psi2;
			
			std::vector<bool> _isSeed;
			Jets _calorimeterJets;
			Particles _calorimeterTowers;
			
	};


}


#endif
