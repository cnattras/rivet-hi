// -*- C++ -*-
#ifndef RIVET_HeavyIonInfo_HH
#define RIVET_HeavyIonInfo_HH

#include "Rivet/Projection.hh"
#include "HepMC/HeavyIon.h"

namespace Rivet
{				


	/// @brief Accessor to HeavyIon HepMC object
	///
	/// 
	/// @todo HepMC3 compatibility uses GenHeavyIon object
	class HeavyIonInfo : public Projection
	{
	public:
		HeavyIonInfo() {}
		
		const bool HasHeavyIonData() const { return _isHI; }
		const int HardCollisions() const { return _isHI ? _hepHI.Ncoll_hard() : -1; }
		const int Collisions() const { return _isHI ? _hepHI.Ncoll() : -1; }
		const int ParticipatingNucleonsProjectile() const { return _isHI ? _hepHI.Npart_proj() : -1; }
		const int ParticipatingNucleonsTarget() const { return _isHI ? _hepHI.Npart_targ() : -1; }
		const int ParticipatingNucleonsTotal() const { return _isHI ? _hepHI.Npart_proj() + _hepHI.Npart_targ() : -1; }
		const int SpectatorNeutrons() const { return _isHI ? _hepHI.spectator_neutrons() : -1; }
		const int SpectatorProtons() const { return _isHI ? _hepHI.spectator_protons() : -1; }
		const int WoundedUnwoundedCollisions() const { return _isHI ? _hepHI.N_Nwounded_collisions() + _hepHI.Nwounded_N_collisions(): -1; }
		const int WoundedWoundedCollisions() const { return _isHI ? _hepHI.Nwounded_Nwounded_collisions() : -1; }
		
		const double ImpactParameter() const { return _isHI ? _hepHI.impact_parameter() : -1; }
		const double b() const { return _isHI ? _hepHI.impact_parameter() : -1; }
		const double EventPlaneAngle() const { return _isHI ? _hepHI.event_plane_angle() : -1; }
		const double Eccentricity() const { return _isHI ? _hepHI.eccentricity() : -1; }
		const double NNInelCrossSection() const { return _isHI ? _hepHI.sigma_inel_NN() : -1; }
		
		
		
		DEFAULT_RIVET_PROJ_CLONE(HeavyIonInfo);
		
	protected:
		void project(const Event& e) {
			
			if (e.genEvent()->heavy_ion()) { 
				_hepHI = (*(e.genEvent()->heavy_ion()));
				_isHI = true;
			} else {
				_isHI = false;
			}
		}
		
		int compare(const Projection& p) const {
			return 1;
		}
		
		HepMC::HeavyIon _hepHI;
		bool _isHI;
		
	};
}

#endif
