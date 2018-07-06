// -*- C++ -*-
#include "Rivet/Projections/ParticleVn.hh"

namespace Rivet
{

	ParticleVn::ParticleVn(const ChargedFinalState& cfs,const std::vector<int>& nOfInterest)
	{

		addProjection(cfs,"CFS");
		addProjection(EventPlane(cfs,nOfInterest),"EventPlane");

		setName("ParticleVn");
		_nOfInterest = nOfInterest;
		_Vn.resize(_nOfInterest.size());
	}


	ParticleVn::ParticleVn(const ChargedFinalState& cfs,int nOfInterest)
	{

		addProjection(cfs,"CFS");
		addProjection(EventPlane(cfs,nOfInterest),"EventPlane");

		setName("ParticleVn");
		_nOfInterest.push_back(nOfInterest);
		_Vn.resize(_nOfInterest.size());
	}
	const double ParticleVn::GetVn(int n = 2) const
	{
		for (unsigned int i = 0; i < _nOfInterest.size(); i++)
			{
				if (_nOfInterest[i] == n) return _Vn[i];
			}
		return -1.;
	}
	void ParticleVn::project(const Event& e)
	{
		const std::vector<double> EP = apply<EventPlane>(e, "EventPlane").EventPlanes();
		const Particles& part = apply<ChargedFinalState>(e,"CFS").particles();
		for (unsigned int i = 0; i < _nOfInterest.size(); i++)
			{
				_Vn[i] = fitToCos(part,EP[i],_nOfInterest[i]);
			}


	}


	double ParticleVn::deltaPhi2(double phi1, double phi2)
	{
		return mapAngle0To2Pi(phi1 - phi2);
	}

	double ParticleVn::fitToCos(const Particles& part, double EPn, int n)
	{
		std::vector<double> phiBins(_n_histbins,0.0);

		double binWidth = 2*pi/_n_histbins;


		foreach (const Particle &p, part)
		{
			if (p.pT() < 5. * GeV)
				{
					phiBins[static_cast<int>(deltaPhi2(p.phi(),EPn)/binWidth)] += p.pT()*GeV;
				}
		}

		double integral = 0.0;
		double midOffset = binWidth * .5;
		double Vn = 0.0;

		foreach (const double &d, phiBins) integral += d;


		for (int i =0; i < _n_histbins; i++)
			Vn += (phiBins[i]) * cos(n * (i * binWidth + midOffset ));


		Vn /= integral;
		return Vn;
	}

}
