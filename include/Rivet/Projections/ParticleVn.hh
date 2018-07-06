// -*- C++ -*-
#ifndef RIVET_ParticleVn_HH
#define RIVET_ParticleVn_HH

#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Math/MathUtils.hh"

namespace Rivet {


	/// @brief Project Out event Vn
	///
	///
	/// @todo 
	class ParticleVn : public Projection {
		public:

		
		/// @name Constructors and destructors.
		//@{

		/// Constructor
		///
		ParticleVn(const ChargedFinalState& cfs,const std::vector<int>& nOfInterest);
		
		ParticleVn(const ChargedFinalState& cfs,int nOfInterest);

		/// Clone on the heap.
		DEFAULT_RIVET_PROJ_CLONE(ParticleVn);

		/// Accessors
		const double GetVnOrdinal(int n = 0) const { return _Vn[n]; }
		const double GetVn(int n = 2) const;
		const std::vector<double> GetVnList() const { return _Vn; }

		



	protected:

		/// Apply the projection to the event.
		void project(const Event& e);

		/// Compare projections 
		int compare(const Projection& p) const {
			return mkNamedPCmp(p, "CFS");
		}
		double deltaPhi2(double phi1, double phi2) ;
		
	private:
		/// cosine fit function
		double fitToCos(const Particles& part, double EPn, int n);

		const int _n_histbins = 50;
		std::vector<int> _nOfInterest;
		std::vector<double> _Vn;
		
	};


}


#endif
