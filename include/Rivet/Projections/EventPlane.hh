// -*- C++ -*-
#ifndef RIVET_EventPlane_HH
#define RIVET_EventPlane_HH

#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"


namespace Rivet
{


	/// @brief Project Out observed Event plane
	///
	///
	/// @todo
	class EventPlane : public Projection
	{
	public:


		/// @name Constructors and destructors.
		//@{

		/// Constructor
		///
		EventPlane(const ChargedFinalState& cfs,const std::vector<int>& nOfInterest) ;

		EventPlane(const ChargedFinalState& cfs,int nOfInterest) ;

		/// Clone on the heap.
		DEFAULT_RIVET_PROJ_CLONE(EventPlane);

		//@}
		const double GetEventPlaneOrdinal(int n = 0) const
		{
			return _EPn[n];
		}

		const double EventPlaneN(int N = 2) const;

		const std::vector<double> EventPlanes() const
		{
			return _EPn;
		}





	protected:

		/// Apply the projection to the event.
		void project(const Event& e);

		/// Compare projections
		// todo: make this do something useful
		int compare(const Projection& p) const
		{
			return mkNamedPCmp(p, "CFS");
		}

	private:

		std::vector<int> _nOfInterest;
		std::vector<double> _EPn;

	};


}


#endif
