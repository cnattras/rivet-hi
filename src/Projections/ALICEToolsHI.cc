// -*- C++ -*-
#include "Rivet/Projections/ALICEToolsHI.hh"

namespace Rivet {
	ALICEToolsHI::ALICEToolsHI( const ChargedFinalState& cfs, ///< Input particles
                      const FastJets& antiKtJets, ///< FastJets Object for your anti-kt Jets
                      double maxAbsEta, ///< Highest |eta| for your particles
					  std::vector<int> vn ///< Vn of interest to remove flow
					  ) : ParticleVn(cfs,vn) {
		  setName("ALICEToolsHI");
		  addProjection(antiKtJets, "antiKtJets");
		  
		  ChargedFinalState v0part(((Cuts::eta < 5.1)&&(Cuts::eta > 2.8))||((Cuts::eta < -1.7)&&(Cuts::eta > -3.7)));
		  declare(v0part,"VOParticles");
		  
		  EventPlane v0plane(v0part, std::vector<int>{2,3});
		  declare(v0plane,"EPV0");
		  
		  declare(EventPlane(ChargedFinalState((Cuts::eta < 0.9) && (Cuts::eta > 0.)),std::vector<int>{2,3}),"EPPOS");
		  declare(EventPlane(ChargedFinalState((Cuts::eta < 0.) && (Cuts::eta > -0.9)),std::vector<int>{2,3}),"EPNEG");
		  
		  FastJets fj(cfs, FastJets::KT, 0.2); // per 10.1007/JHEP03(2014)013 no significant R parameter dependence 
		  fj.useJetArea(new fastjet::AreaDefinition(fastjet::active_area,fastjet::GhostedAreaSpec(0.7, 1,1./200.)));
			              // using area spec specified in 10.1016/j.physletb.2015.12.047
          declare(fj, "KtJetsD02");
		 
		  _maxAbsEta = maxAbsEta;
		}
		
		const double ALICEToolsHI::RhoLocal(double phi) const {
			double bg = 1; // first term in decomposition
			double dPhi;
			double Rinv = 1/_jetR;
			dPhi = mapAngle0To2Pi(phi-_epN[0]);
			bg += Rinv * _vn[0] * cos(2* dPhi) * _sin2R; // second term
			dPhi = mapAngle0To2Pi(phi-_epN[1]);
			bg += (2/3) * Rinv * _vn[1] * cos(3* dPhi) * _sin3R; // third term
				
			return _backgroundMedianDensity*bg;
		} 
		
		void ALICEToolsHI::project(const Event& e){
			
			// Get the area-filtered jet inputs for computing median energy density, etc.
			vector<double> ptDensity;
			FastJets KtJets = apply<FastJets>(e, "KtJetsD02");
			Jets Jetskt = KtJets.jets();
			const auto KtCSA = KtJets.clusterSeqArea();
			FastJets AntiKtJets = apply<FastJets>(e, "antiKtJets");
			_jetR = AntiKtJets.clusterSeq()->jet_def().R();
			Jets Jetsakt = AntiKtJets.jetsByPt();
			double leadingOrderJetEta = Jetsakt[0].eta();
			
			foreach (const Jet& jet, Jetskt) {
				const double area = KtCSA->area(jet);
				if (jet.pT() >= 0.15 && area >= 0.01) 
					ptDensity.push_back(jet.pT()/area);
			}
			
			std::sort(ptDensity.begin(),ptDensity.end()); // sort ascending
			int nDensities = ptDensity.size() -2 ; // exclude 2 highest density
			_backgroundMedianDensity = ptDensity[(nDensities/2)-1]; // median 
			
			_eventPart = applyProjection<ChargedFinalState>(e, "CFS").particles((Cuts::eta > leadingOrderJetEta+ _jetR) || (Cuts::eta < leadingOrderJetEta - _jetR));
			
			double pTSum = 0.0;
			foreach (const Particle& p, _eventPart){
					pTSum += p.pT();
			}
			_backgroundRho0 = pTSum /= 4*pi*(_maxAbsEta - _jetR); //(2*_maxAbsEta - 2*_jetR)*2*pi;
			
			_epN = apply<EventPlane>(e, "EventPlane").EventPlanes();
			for (unsigned int i = 0; i < _nOfInterest.size(); i++){
				_vn[i] = _fitToCos(_eventPart,_epN[i],_nOfInterest[i]);
			}
			_sin2R = sin(2*_jetR);
			_sin3R = sin(3*_jetR);
			
			////
			// EP Resolution calcs
			//
			_nEvent++;
			
			std::vector<double> _epNV0 = apply<EventPlane>(e, "EPV0").EventPlanes();
			std::vector<double> _epNPos = apply<EventPlane>(e, "EPPOS").EventPlanes();
			std::vector<double> _epNNeg = apply<EventPlane>(e, "EPNEG").EventPlanes();
			_r2PosEta += cos(2*(_epNV0[0]-_epNPos[0]));
			_r2NegEta += cos(2*(_epNV0[0]-_epNNeg[0]));
			_r2DifEta += cos(2*(_epNPos[0]-_epNNeg[0]));
			_r3PosEta += cos(3*(_epNV0[1]-_epNPos[1]));
			_r3NegEta += cos(3*(_epNV0[1]-_epNNeg[1]));
			_r3DifEta += cos(3*(_epNPos[1]-_epNNeg[1]));
			
		}
		
		
		int ALICEToolsHI::compare(const Projection& p) const {
			const ALICEToolsHI& other = dynamic_cast<const ALICEToolsHI&>(p);
		    if (other._nOfInterest.size() != _nOfInterest.size()) return UNDEFINED;
			for (int i = 0; i < (signed)_nOfInterest.size(); ++i) {
			  if (other._nOfInterest[i] != _nOfInterest[i]) return UNDEFINED;
			}
			if (_nBins != other._nBins) return UNDEFINED;
			if (_maxAbsEta != other._maxAbsEta) return UNDEFINED;
			if (mkNamedPCmp(p,"antiKtJets") != EQUIVALENT) return UNDEFINED;
			return mkNamedPCmp(p, "CFS");
		}


}


