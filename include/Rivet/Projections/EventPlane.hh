// -*- C++ -*-
#ifndef RIVET_EventPlane_HH
#define RIVET_EventPlane_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet
{

	/// @brief Project Out Observed Event Plane(s)
	///
	/// Does stuff
	
	class EventPlane : public Projection
	{
		public:


			
			EventPlane( const FinalState& fs, const std::vector<int>& nOfInterest,  const Cut& c=(Cuts::pT < 5.0) );

			EventPlane( const FinalState& fs, int nOfInterest, const Cut& c=(Cuts::pT < 5.0) );

			///
			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(EventPlane);

			/////////////////////////////////////////////////////////////////////
			/// Returns an Event Plane by order of passed n of interest
			/////////////////////////////////////////////////////////////////////
			const double GetEventPlaneOrdinal(int n = 0) const
			{
				return _epN[n];
			}

			/////////////////////////////////////////////////////////////////////
			/// Returns the nth Event Plane
			/////////////////////////////////////////////////////////////////////
			const double EventPlaneN(int N = 2) const
			{
				for (unsigned int i = 0; i < _nOfInterest.size(); i++)
				{
					if (_nOfInterest[i] == N) return _epN[i];
				}
				return -100.;

			}

			/////////////////////////////////////////////////////////////////////
			/// Returns a vector of all selected event planes
			/////////////////////////////////////////////////////////////////////
			const std::vector<double> EventPlanes() const
			{
				return _epN;
			}





		protected:

			/// Apply the projection to the event.
			void project(const Event& e);

			/// Compare projections
			int compare(const Projection& p) const;


			std::vector<int> _nOfInterest;
			std::vector<double> _epN;
			Cut _pCut;
	};


	

}






#endif
