// -*- C++ -*-
#ifndef RIVET_ParticleVn_HH
#define RIVET_ParticleVn_HH

#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Math/MathUtils.hh"

namespace Rivet
{


	/// @brief Project Out Particle Flow Vn
	///
	///
	/// @todo Add Cuts
	class ParticleVn : public Projection
	{
		public:

			static const int SQRT_N = -1;

			/// @name Constructors and destructors.
			//@{


			///
			/// Constructor
			ParticleVn( const FinalState& fs, const std::vector<int>& nOfInterest, int nBins = ParticleVn::SQRT_N);
			///
			/// Constructor
			ParticleVn( const FinalState& fs, int nOfInterest, int nBins = ParticleVn::SQRT_N);

			//@}
			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(ParticleVn);

			/// Return Vn for a an entry (by order in passed vector)
			const double GetVnOrdinal( int i = 0 ///< Vn for n in ith position passed in nOfInterest
			                           ) const
			{

				if (i < (int)_vn.size() ) return _vn[i];
				else return -1.;
			}


			/// Return the nth order Vn (returns -1 if it does not exist)
			const double GetVn(int n = 2) const
			{
				for (unsigned int i = 0; i < _nOfInterest.size(); i++)
				{
					if (_nOfInterest[i] == n) return _vn[i];
				}
				return -1.;
			}

			/// Return full list of calculated Vn
			const std::vector<double> GetVnList() const
			{
				return _vn;
			}

			/// Return event plane for nth calculated entry by order passed to calculate
			const double GetEventPlaneOrdinal(int i = 0) const
			{
				return _epNvn[i];
			}

			/// Get \f$\psi_{n}\f$ 
			const double EventPlaneN(int n = 2) const
			{
				for (unsigned int i = 0; i < _nOfInterest.size(); i++)
				{
					if (_nOfInterest[i] == n) return _epNvn[i];
				}
				return -100.;

			}

			/// Return full list of calculated Vn
			const std::vector<double> EventPlanes() const
			{
				return _epNvn;
			}

		protected:

			/// Apply the projection to the event.
			void project(const Event& e);

			/// Compare projections
			int compare(const Projection& p) const;


			/// Calculate delta phi 0-2pi (not absolute value)
			double deltaPhi2(double phi1, double phi2)
			{
				return mapAngle0To2Pi(phi1 - phi2);
			}


			/// Cosine fit function
			double _fitToCos(const Particles& part, double EPn, int n);

			std::vector<int> _nOfInterest;
			std::vector<double> _vn;
			std::vector<double> _epNvn;
			int _nBins;
	};

}


#endif
