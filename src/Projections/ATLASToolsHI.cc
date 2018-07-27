// -*- C++ -*-
#include "Rivet/Projections/ATLASToolsHI.hh"

namespace Rivet {

	ATLASToolsHI::ATLASToolsHI( double trackJetPtQualifier,
	                            double EMTowerQualifier
	                          ) : _tJetPt(trackJetPtQualifier), _towEt(EMTowerQualifier)
	{
		DetectorPixels act(DetectorPixels::ATLAS_MidRapidity);
		setName("ATLASToolsHI");
		addProjection(act, "ACT");

		const Cut &c = act.Cuts();

		FinalState tracks(c && Cuts::pt > 0.15);
		declare(tracks,"TRACKS");

		double minRap = act.RapidityMin();
		double maxRap = act.RapidityMax();
		_nPhi = act.PhiTowers();
		_wPhi = act.PhiWidth();
		_nRap = act.RapTowers();
		_wRap = act.RapWidth();
		_minRap = minRap;
		_maxRap = maxRap;
		_nRap2 = _nRap /2;
		_towerArea = _wPhi * _wRap;


		DetectorPixels fcal(DetectorPixels::ATLAS_FCal);
		declare(fcal,"FCAL");

		ParticleVn fcalVn(fcal,2,_nPhi);
		declare(fcalVn,"FCALVN");


		FastJets jets(act, FastJets::ANTIKT, 0.2, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
		declare(jets, "CJETS");

		FastJets trackjets(tracks, FastJets::ANTIKT, 0.4, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
		declare(trackjets, "TJETS");


	}
	void ATLASToolsHI::project(const Event& e)
	{
		_calorimeterJets.clear();
		_calorimeterTowersInclGhosts = apply<DetectorPixels>(e, "ACT").particles();
		_calorimeterTowersRaw = Particles(_calorimeterTowersInclGhosts);

		// Get Calorimeter Jets
		Jets calorimeterJets(apply<FastJets>(e, "CJETS").jets());
		Jets caloJetsAdj(calorimeterJets);

		// Find seed jets
		_tagSeedJets1(caloJetsAdj);

		// Get v2 from FCAL
		const ParticleVn &fcalVn = apply<ParticleVn>(e, "FCALVN");
		_v2 = fcalVn.GetVn(2);
		_psi2 = fcalVn.EventPlaneN(2);

		// Calculate rhoi(eta) for mid rapidities;
		_calculateRhoi1(caloJetsAdj);

		// Adjust Jet 4 vectors to account for initial background
		caloJetsAdj = calorimeterJets;
		_subtractInterimUE(caloJetsAdj);

		// Find New Seed Jets including Track Seed Jets
		const Jets &trackJets = apply<FastJets>(e, "TJETS").jets(Cuts::pt > _tJetPt * GeV); // per 1208.1967v2
		Jets cSeeds, tSeeds;
		_getSeedJets2(caloJetsAdj, trackJets, cSeeds, tSeeds);

		// Get new rhoi and V2
		_calculateRhoi2( _calorimeterTowersInclGhosts, cSeeds, tSeeds);
		_calcV2Psi2( _calorimeterTowersInclGhosts, cSeeds, tSeeds );

		// Subtrack UE from towers, create final jets
		_subtractUE();
		_finalizeJets(calorimeterJets);
		_calorimeterTowers = filter_select(_calorimeterTowersInclGhosts, Cuts::Et > 1e-5);


	}

	int ATLASToolsHI::compare(const Projection& p) const
	{
		const ATLASToolsHI& other = dynamic_cast<const ATLASToolsHI&>(p);
		if (other._tJetPt != _tJetPt) return UNDEFINED;

		return _towEt == other._towEt ? EQUIVALENT : UNDEFINED;
	}

	void ATLASToolsHI::_tagSeedJets1(const Jets &js)
	{
		_isSeed = std::vector<bool>(js.size(),false);
		for (int i =0; i < (signed)js.size(); ++i)
		{
			const Particles& ps = js[i].particles();
			bool has3gevTower = false;
			int nPart = ps.size();
			double sumEt = 0.;
			double maxEt = 0.;
			foreach (const Particle& p, ps)
			{
				if (p.Et() > 3 * GeV) has3gevTower = true;
				sumEt += p.Et();
				if (p.Et() > maxEt) maxEt = p.Et();
			}
			double avgEt = sumEt / nPart;
			if (has3gevTower && maxEt/avgEt > 4.0) _isSeed[i] = true;
		}

	}

	void ATLASToolsHI::_calculateRhoi1(const Jets& js)
	{
		_rhoi = std::vector<double>(_nRap, 0.0);
		_rhoiEntries = std::vector<int>(_nRap, 0);
		int nJet = js.size();

		for(int i = 0; i < nJet; ++i)
		{
			if (!_isSeed[i])
			{
				const Particles& ps = js[i].constituents();
				foreach (const Particle& p, ps)
				{

					int bin = _rapBin(p.eta());
					_rhoi[bin] += p.Et();
					++_rhoiEntries[bin];

				}

			}
		}

		for (int i = 0; i < _nRap; ++i)
		{
			if (_rhoiEntries[i]) _rhoi[i] /= _rhoiEntries[i] * _towerArea;

		}
	}


	void ATLASToolsHI::_subtractInterimUE(Jets &js)
	{
		foreach (Jet& j, js)
		{
			FourMomentum mom(0,0,0,0);
			foreach (const Particle &p, j.particles())
			{
				mom += _adjMom(p);
			}
			j.setState(mom, j.particles());
		}
	}
	FourMomentum ATLASToolsHI::_adjMom(const Particle& p)
	{
		double theta = p.theta();
		double phi = ghostPhi(p);
		double Et = p.Et();
		double Eadj = 1e-20; // ghost energy
		double EtUE = _towerArea * _rhoi[_rapBin(p.eta())] * (1 + 2 * _v2 * cos(2 * (phi - _psi2)));
		if (Et-EtUE > 0.0) Eadj = (Et - EtUE) / sin(theta);

		return FourMomentum(Eadj,Eadj*cos(phi)*sin(theta),Eadj*sin(phi)*sin(theta),Eadj*cos(theta));

	}

	void ATLASToolsHI::_getSeedJets2(const Jets& caloJets, const Jets& trackJets, Jets& cSeeds, Jets& tSeeds)
	{

		foreach (const Jet& j, caloJets)
		{
			if (j.Et() > 25 * GeV)
			{
				cSeeds.push_back(j);
			}
		}
		foreach (const Jet& j, trackJets)
		{
			if(any(j.particles(),PtGtr(4*GeV)))
			{
				tSeeds.push_back(j);
			}
		}
		_trackJets = tSeeds;
	}


	void ATLASToolsHI::_calculateRhoi2(const Particles& towers, const Jets& cSeeds, const Jets& tSeeds)
	{

		_rhoi = std::vector<double>(_nRap, 0.0);
		_rhoiEntries = std::vector<int>(_nRap, 0);

		foreach (const Particle &p, towers)
		{
			if (!any(cSeeds,DeltaRLess(p,0.4)) && !any(tSeeds,DeltaRLess(p,0.4)))
			{
				int bin = _rapBin(p.eta());
				_rhoi[bin] += p.Et();
				++_rhoiEntries[bin];
			}
		}

		for (int i = 0; i < _nRap; ++i)
		{
			if (_rhoiEntries[i]) _rhoi[i] /= _rhoiEntries[i] * _towerArea;
		}
	}

	void ATLASToolsHI::_calcV2Psi2(const Particles& towers, const Jets& cSeeds, const Jets& tSeeds)
	{
		double Qx = 0.0;
		double Qy = 0.0;
		foreach (const Particle& p, towers)
		{
			if (!any(cSeeds,DeltaRLess(p,0.4)) && !any(tSeeds,DeltaRLess(p,0.4)))
			{
				Qx += p.Et() * cos(2 * p.phi());
				Qy += p.Et() * sin(2 * p.phi());
			}

		}
		_psi2 = 0.5 * atan2(Qy,Qx);

		double EtSum = 0.0;
		std::vector<double> EtPhi(_nPhi, 0.0);

		foreach (const Particle& p, towers)
		{
			if (!any(cSeeds,DeltaRLess(p,0.4)) && !any(tSeeds,DeltaRLess(p,0.4)))
			{
				int bin = static_cast<int>(p.phi() /_wPhi);
				EtPhi[bin] += p.Et();
				EtSum += p.Et();
			}
		}

		_v2 = 0.;
		for (int i = 0; i < _nPhi; ++i)
		{
			double phi = (0.5 + static_cast<double>(i)) * _wPhi;
			_v2 += EtPhi[i] * cos(2*(phi - _psi2));
		}
		_v2 /= EtSum;

	}

	void ATLASToolsHI::_subtractUE()
	{
		foreach (Particle &p, _calorimeterTowersInclGhosts)
		{
			p.setMomentum(_adjMom(p));
		}

	}

	void ATLASToolsHI::_finalizeJets(const Jets& js)
	{


		foreach (const Jet &j, js)
		{


			if (j.particles().size() !=0)
			{
				Particles jPart;
				bool matched = false;
				FourMomentum mom(0,0,0,0);
				foreach (const Particle &p, j.particles())
				{
					const FourMomentum& pMom = _calorimeterTowersInclGhosts[p.pid()].momentum();
					mom += pMom;
					if (pMom.Et() >= _towEt) matched = true;
					jPart.push_back(_calorimeterTowersInclGhosts[p.pid()]);

				}
				Jet calJet(mom,Particles(jPart));
				if (!matched)
				{
					if(any(_trackJets,DeltaRLess(calJet,0.2))) matched = true;
				}
				if (matched) _calorimeterJets.push_back(calJet);
			}
		}
		_theParticles = filter_select(_calorimeterTowersInclGhosts, Cuts::pt > 0.1 * GeV);
	}

	int ATLASToolsHI::_rapBin(double eta)
	{
		if (eta < 0) return static_cast<int>((eta+_maxRap)/_wRap);
		else return _nRap2 + static_cast<int>((eta-_minRap)/_wRap);
	}


}
