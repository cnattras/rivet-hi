// -*- C++ -*-
#ifndef RIVET_EventPlane_HH
#define RIVET_EventPlane_HH

#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet {


	/// @brief Project Out observed Event plane 
	///
	///
	/// @todo 
	class EventPlane : public Projection {
		public:

		
		/// @name Constructors and destructors.
		//@{

		/// Constructor
		///
		EventPlane(const ChargedFinalState& cfs, const std::vector<int>& nOfInterest, const Cut& c=(Cuts::pT < 5.0)) {
		  setName("EventPlane");
		  addProjection(cfs, "CFS");
		  _nOfInterest = nOfInterest;
		  _epN.resize(_nOfInterest.size());
		  _pCut = c;
		}
		EventPlane(const ChargedFinalState& cfs,int nOfInterest, const Cut& c=(Cuts::pT < 5.0)) {
		  setName("EventPlane");
		  addProjection(cfs, "CFS");
		  _nOfInterest.push_back(nOfInterest);
		  _epN.resize(_nOfInterest.size());
		  _pCut = c;
		}

		/// Clone on the heap.
		DEFAULT_RIVET_PROJ_CLONE(EventPlane);

		//@}
		const double GetEventPlaneOrdinal(int n = 0) const { return _epN[n]; }
		
		const double EventPlaneN(int N = 2) const { 
			for (unsigned int i = 0; i < _nOfInterest.size(); i++){
				if (_nOfInterest[i] == N) return _epN[i]; 
			}
			return -100.; 
			
		}
		const std::vector<double> EventPlanes() const { return _epN; }

		



		protected:

		/// Apply the projection to the event.
		void project(const Event& e){
			
			const Particles& parts = apply<ChargedFinalState>(e, "CFS").particles(_pCut);
			
			for (unsigned int i = 0; i < _nOfInterest.size(); i++){
				double n = static_cast<double>(_nOfInterest[i]);	
				double Qx = 0.0;
				double Qy = 0.0;
				foreach (const Particle& p, parts){
					Qx += p.pT()*cos(n*p.phi());
					Qy += p.pT()*sin(n*p.phi());
				}
				
				_epN[i] = atan2(Qy,Qx)/ n;
				
			}
			
		}

		/// Compare projections (only difference is in UFS definition)
		int compare(const Projection& p) const {
		  return mkNamedPCmp(p, "CFS");
		}
		
	private:
		
		std::vector<int> _nOfInterest;
		std::vector<double> _epN;
		Cut _pCut;
	};


}


#endif
