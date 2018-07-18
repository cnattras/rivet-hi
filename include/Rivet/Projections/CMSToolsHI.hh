// -*- C++ -*-
#ifndef RIVET_BackgroundFinder_HH
#define RIVET_BackgroundFinder_HH

//#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
//#include "Rivet/Projections/ParticleVn.hh"
//#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
//#include <algorithm>
//#include <cmath>

namespace Rivet {


	/// @brief Project Out the Background of an event
	///
	///
	/// @todo 
	class BackgroundFinder : public Projection {
		public:

		enum BackgroundMethod {
		  CMS_Reflection,
		  ATLAS,
		  Undefined
		};

		/// @name Constructors and destructors.
		//@{

		/// Constructor
		///
		BackgroundFinder(const ChargedFinalState& cfs, BackgroundMethod bgm) {
		  setName("BackgroundFinder");
		  addProjection(cfs, "CFS");
		  _bgm = bgm;
		  
		}


		/// Clone on the heap.
		DEFAULT_RIVET_PROJ_CLONE(BackgroundFinder);

		//@}


		/// @name 
		//@{


		/// 
		const Particles CMSReflectionParticles(const Jet& jetOfInterest, double R_Param) const {
			Particles particles;
			particles = getConeParticles(eventPart,-jetOfInterest.eta(),jetOfInterest.phi(), R_Param);
	    	return particles;
		}
		const FourMomentum CMSReflectionFourMomentum(const Jet& jetOfInterest, double R_Param) const {
			FourMomentum vec(0,0,0,0);
			Particles particles = CMSReflectionParticles(jetOfInterest, R_Param);
			foreach (const Particle &p, particles){
				vec += p.momentum();
			}
			return vec;
		}
		
		//@}

		const Particles getConeParticles(const Jet& jet,double R) const{
			return getConeParticles(jet.constituents(),jet.eta(), jet.phi(), R);
		}
		const Particles getConeParticles(const Particles& ps, double eta, double phi,double R) const{
			Particles conePart;
			foreach (const Particle &p, ps){
				if (deltaR(p,eta,phi) < R) conePart.push_back(p);
			}
			return conePart;
		}

		protected:

		/// Apply the projection to the event.
		void project(const Event& e){
			eventPart = applyProjection<ChargedFinalState>(e, "CFS").particles();
		}

		/// Compare projections (only difference is in UFS definition)
		int compare(const Projection& p) const {
		  return mkNamedPCmp(p, "CFS");
		}
		
	private:
		double deltaPhi2(double phi1, double phi2) {
			return mapAngle0To2Pi(phi1 - phi2);
		}
		
		
		
		Particles eventPart;
		BackgroundMethod _bgm;
		
	
	};


}


#endif
