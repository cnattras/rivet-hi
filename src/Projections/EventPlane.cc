// -*- C++ -*-
#include "Rivet/Projections/EventPlane.hh"


namespace Rivet
{
		
	EventPlane::EventPlane(const ChargedFinalState& cfs,const std::vector<int>& nOfInterest)
	{
		setName("EventPlane");
		addProjection(cfs, "CFS");
		_nOfInterest = nOfInterest;
		_EPn.resize(_nOfInterest.size());
	}
	EventPlane::EventPlane(const ChargedFinalState& cfs,int nOfInterest)
	{
		setName("EventPlane");
		addProjection(cfs, "CFS");
		_nOfInterest.push_back(nOfInterest);
		_EPn.resize(_nOfInterest.size());
	}


	const double EventPlane::EventPlaneN(int N = 2) const
	{
		for (unsigned int i = 0; i < _nOfInterest.size(); i++) {
			if (_nOfInterest[i] == N) return _EPn[i];
		}
		return -100.;

	}

	void EventPlane::project(const Event& e)
	{
		const Particles& parts = apply<ChargedFinalState>(e, "CFS").particles();

		for (unsigned int i = 0; i < _nOfInterest.size(); i++) {
			double n = static_cast<double>(_nOfInterest[i]);
			double Qx = 0.0;
			double Qy = 0.0;
			foreach (const Particle& p, parts) {
				if (p.pT() < 5. * GeV) {
					Qx += p.pT()*cos(n*p.phi());
					Qy += p.pT()*sin(n*p.phi());
				}
			}
			_EPn[i] = atan2(Qy,Qx)/ n;

		}


	}



}
