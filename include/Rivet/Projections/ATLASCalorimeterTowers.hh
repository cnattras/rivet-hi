// -*- C++ -*-
#ifndef RIVET_ATLASCalorimeterTowers_HH
#define RIVET_ATLASCalorimeterTowers_HH

#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


	/// @brief Project Out the Background of an event
	///
	///
	/// @todo 
	class ATLASCalorimeterTowers : public FinalState {
		public:

		/// @name Constructors and destructors.
		//@{

		/// Constructor
		///
		
		ATLASCalorimeterTowers() : ATLASCalorimeterTowers( Cuts::abseta < 2.1, 0.,  2.1, 0.1,  0.1) {}
		
		ATLASCalorimeterTowers( const Cut& c , double minAbsRap = 0, double maxAbsRap = 2.1, double towerWidthRap = 0.1, double towerWidthPhi = 0.1) {
			setName("ATLASCalorimeterTowers");
			ChargedFinalState cfs(Cuts::pT > 0.15 * GeV && c);
			addProjection(cfs, "CFS");
			_nPhi = static_cast<int>(2*pi / towerWidthPhi);
			_wPhi = 2.*pi / _nPhi;
			_nRap = 2*static_cast<int>((maxAbsRap-minAbsRap)/towerWidthRap);
			_wRap = (2*(maxAbsRap-minAbsRap))/_nRap;
			cout << _wRap << " wrap";
			_minRap = minAbsRap;
			_maxRap = maxAbsRap;
			_nRap2 = _nRap /2;
			_towerMom0 = std::vector<FourMomentum>(_nPhi*_nRap,FourMomentum(0,0,0,0));
			
			_theParticles = Particles(_nPhi * _nRap);
			
		}
		

		/// Clone on the heap.
		DEFAULT_RIVET_PROJ_CLONE(ATLASCalorimeterTowers);

		//@}

		const double PhiWidth() const { return _wPhi; }
		const double RapWidth() const { return _wRap; }
		const int PhiTowers() const { return _nPhi; }
		const int RapTowers() const { return _nRap; }
		const double RapidityMin() const { return _minRap; }
		const double RapidityMax() const { return _maxRap; }
		const double TowerArea() const { return _wPhi * _wRap; }
		
		/// @name 
		//@{
		
		//}

	
		//@}
//		void ApplyEtSubtraction(std::vector<double> Et) {
//			_theParticles.resize(_nRap * _nPhi);
//			int nPart = _theParticles.size();
//			
//			for (int i = 0; i < nPart; ++i){
//				double theta = _theParticles[i].theta();
//				double phi = _theParticles[i].phi();
//				double E = _theParticles[i].E();
//				double Eadj;
//				if ((_theParticles[i].Et() - Et[i]) > 0.0) Eadj = (_theParticles[i].Et() - Et[i])/sin(theta);
//				else Eadj = 1e-20;
//				FourMomentum mom(Eadj,Eadj*cos(phi)*sin(theta),Eadj*sin(phi)*sin(theta),Eadj*cos(theta));
//				_theParticles[i].setMomentum(mom);
//			}
//			
//		}
		
		
		protected:

		/// Apply the projection to the event.
		void project(const Event& e){
			Particles parts = apply<ChargedFinalState>(e,"CFS").particles();
		    std::vector<FourMomentum> towerMom = _towerMom0;
			
			int binRap2 =  _nRap / 2;
			
			foreach (const Particle& p, parts){
				int binPhi = static_cast<int>(p.phi()/_wPhi);
				
				int binRap;
				double pEta = p.eta();
				if (pEta < 0) binRap = static_cast<int>((p.eta()+_maxRap)/_wRap);
				else binRap = binRap2 + static_cast<int>((p.eta()-_minRap)/_wRap);
				towerMom[binPhi + (_nPhi * binRap)] += p.mom();
			}
			

			
			for (int i = 0; i < _nPhi * _nRap ; ++i){
				int binRap = i / _nPhi;
				int binPhi = i % _nPhi;
				
				double rap;
				
				if (binRap < binRap2) rap = -1 * _maxRap + ((binRap+0.5) * _wRap) ;
				else rap =  _minRap + ((binRap+0.5-binRap2) * _wRap);
			
				
				
				double phi = (binPhi+0.5) * _wPhi;
				double p = towerMom[i].p();
				double E = towerMom[i].E();
				if (E < 1e-20) E =1e-20; // for non-zero towers, make into ghosts with ~0 Et
				if (p < 1e-20) p =1e-20; // ghost 0 < .1 MeV leads to miscalculated phi
				double theta = 2*atan(exp(-rap));
//				if (p == 1e-10) {
//					cout << "ghost " << E << ' ' << p*cos(phi)*sin(theta) << ' ' << p*sin(phi)*sin(theta) << ' ' << p*cos(theta) << endl;
//					cout << atan2(p*sin(phi)*sin(theta), p*cos(phi)*sin(theta) ) << " atan\n";
//					cin.get();
//				}
				FourMomentum mom(E,p*cos(phi)*sin(theta),p*sin(phi)*sin(theta),p*cos(theta));
				//mom.setEtaPhiMPt(rap, phi, 1e-20, E);
				// use the index as the PID, for later tracking use
				_theParticles[i] = Particle(i, mom, FourVector(0,0,0,0));
//				_theParticles[i].setMomentum();

				if (_theParticles[i].phi() == 0.0) 
					cout << "BAD PHI " << _theParticles[i].eta() << ' ' << phi << ' '<< _theParticles[i].phi() << ' '<< _theParticles[i].px() << ' ' << _theParticles[i].py() << endl;//<< p << " p " << E << "E\n";
			    
			}
			
	
		}

		/// Compare projections (only difference is in UFS definition)
		int compare(const Projection& p) const {
		  return -9876;
		}
		
	private:
//		double deltaPhi2(double phi1, double phi2) {
//			return mapAngle0To2Pi(phi1 - phi2);
//		}
		
		
		std::vector<FourMomentum> _towerMom0;
		//Particles _the;
		double _wPhi;
		int _nPhi;
		double _wRap;
		int _nRap;
		int _nRap2;
		double _minRap;
		double _maxRap;
	//	Particles _theParticles;
		
		
	};


}


#endif
