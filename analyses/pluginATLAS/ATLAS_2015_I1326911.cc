
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/ParticleBase.hh"

#include<cmath>
#include<iostream>
 
namespace Rivet {


  
  class ATLAS_2015_I1326911 : public HeavyIonAnalysis {
  public:

    /// Constructor
    ATLAS_2015_I1326911(): HeavyIonAnalysis("ATLAS_2015_I1326911"){}


    /// Book histograms and initialise projections

    void init() {
		HeavyIonAnalysis::init();
		
	addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 20,"ImpactParameterMethod");

	const FinalState fs;
	addProjection(fs,"FS");
	

	FastJets fj(fs, FastJets::ANTIKT, 0.4);
	addProjection(fj, "FJ");
	
	
	//JetDefinition jet_def(FastJets::JetAlgorithm ANTIKT, 0.2, FastJets::ReconstuctionScheme recomb_scheme = E_scheme);
	
	//ClusterSequence cs(_theParticles, jet_def);
    //vector<PseudoJet> jets1 = sort_by_E(cs.inclusivejets());
	
    /*for (int n=0; n < jets1.size(); n++) {
        double mnt = 0;
        mnt += jets1[n].E()/GeV;
        //if (jets1[n].E()/GeV > 3 && jets1[0]/ /*average of all energies){
            /* extract jet from vector 
        //}
    }
    double avg = (mnt/jets1.size());
    vector<double>newv;
    for(int i = 0; i < jets1.size();i++){
        if (jets1[n].E()/GeV < 3 && (jets1[0].E()/ avg) < 4){
            newv.push_back(jets1[n].E());
        }
        
    }*/
    
        

   // Book histograms
	  
	 //Figure 1

      _h_DDCS_R4_Y1 = bookHisto1D("d02-x01-y01");
      _h_DDCS_R4_Y2 = bookHisto1D("d03-x01-y01");
	  _h_DDCS_R4_Y3 = bookHisto1D("d04-x01-y01");
	  _h_DDCS_R4_Y4 = bookHisto1D("d05-x01-y01");
	  _h_DDCS_R4_Y5 = bookHisto1D("d06-x01-y01");
	  
	  //Figure 2

	  _h_JY_Y21_C1 = bookHisto1D("d07-x01-y01");
	  _h_JY_Y21_C2 = bookHisto1D("d08-x01-y01");
	  _h_JY_Y21_C3 = bookHisto1D("d09-x01-y01");
	  _h_JY_Y21_C4 = bookHisto1D("d10-x01-y01");
	  
	  _h_JY_C1_Y1 = bookHisto1D("d11-x01-y01");
	  _h_JY_C1_Y2 = bookHisto1D("d12-x01-y01");
	  _h_JY_C1_Y3 = bookHisto1D("d13-x01-y01");
	  _h_JY_C1_Y4 = bookHisto1D("d14-x01-y01");
	  
	  //Figure 3

	  RAA_Y1_C1 = bookHisto1D("d15-x01-y01");
	  RAA_Y1_C2 = bookHisto1D("d16-x01-y01");
	  RAA_Y1_C3 = bookHisto1D("d17-x01-y01");
	  
	  RAA_Y2_C1 = bookHisto1D("d18-x01-y01");
	  RAA_Y2_C2 = bookHisto1D("d19-x01-y01");
	  RAA_Y2_C3 = bookHisto1D("d20-x01-y01");
	  
	  RAA_Y3_C1 = bookHisto1D("d21-x01-y01");
	  RAA_Y3_C2 = bookHisto1D("d22-x01-y01");
	  RAA_Y3_C3 = bookHisto1D("d23-x01-y01");
	  
	  //Figure 4

	  RAA_C1 = bookHisto1D("d24-x01-y01");
	  RAA_C2 = bookHisto1D("d25-x01-y01");
	  RAA_C3 = bookHisto1D("d26-x01-y01");
	  
	  RAA_Y1 = bookHisto1D("d27-x01-y01");
	  
    }


    /// Perform the per-event analysis

    void analyze(const Event& event) {
		
	const double c = centrality(event, "ImpactParameterMethod");
      //cout << c << "cent\n";
      if((c < 0.) || (c > 80.)){
        vetoEvent;
      }

	const double weight = event.weight();
	
	//loop for first rapidity range
	
	const Jets& jets = applyProjection<FastJets>(event,"FJ").jetsByPt(Cuts::pT > 32*GeV && Cuts::pT < 500*GeV);
	
	//Table1
	
	foreach (const Jet& j, jets) {
	
		if(inRange(j.absrap(),0.0,2.1)){
	    _h_DDCS_R4_Y1->fill(j.pT()/GeV, weight);
		}
		
		
		if(inRange(j.absrap(),0.0,0.3)){
		_h_DDCS_R4_Y2->fill(j.pT()/GeV, weight);
		}
		   
		if(inRange(j.absrap(),0.3,0.8)){
		_h_DDCS_R4_Y3->fill(j.pT()/GeV, weight);
		}
		
		if(inRange(j.absrap(),0.8,1.2)){
		_h_DDCS_R4_Y4->fill(j.pT()/GeV, weight);
		}
		
		if(inRange(j.absrap(),1.2,2.1)){
		_h_DDCS_R4_Y5->fill(j.pT()/GeV, weight);
		}
		
		
	}
	
	//Table2
	foreach (const Jet& j, jets) {
		if(inRange(j.absrap(),0.0,2.1)){
	
		if (c >= 0 && c <= 10){
			foreach (const Jet& j, jets){
			_h_JY_Y21_C1 ->fill(j.pT()/GeV, weight);
			count1++;
			Npartt1 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_targ(): 0;
			Npartp1 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_proj(): 0;
			}
		}	
	
		if (c >= 20 && c <= 30){
			foreach (const Jet& j, jets) {
			_h_JY_Y21_C2 ->fill(j.pT()/GeV, weight);
			count2++;
			Npartt2 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_targ(): 0;
			Npartp2 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_proj(): 0;
			}
		}
	
		if (c >= 30 && c <= 40){
			foreach (const Jet& j, jets) {
			_h_JY_Y21_C3 ->fill(j.pT()/GeV, weight);
			count3++;
			Npartt3 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_targ(): 0;
			Npartp3 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_proj(): 0;
			}
		}
	
		if (c >= 60 && c <= 80){
			foreach (const Jet& j, jets) {
			_h_JY_Y21_C4 ->fill(j.pT()/GeV, weight);
			count4++;
			Npartt4 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_targ(): 0;
			Npartp4 += event.genEvent() ->heavy_ion() ? event.genEvent() ->heavy_ion() ->Npart_proj(): 0;
			}
		}
		
	}
	}
	int Nevt= count1 + count2 + count3 + count4;
	int NpartTot= Npartt1 + Npartp1 + Npartt2 + Npartp2 + Npartt3 + Npartp3 + Npartp4 + Npartt4;
	
	
	if (c >= 0 && c <= 10){		
		foreach (const Jet& j, jets){
		
		if(inRange(j.absrap(),0.0,0.3)){
		_h_JY_C1_Y1 ->fill(j.pT()/GeV, weight);
		}
		
		if(inRange(j.absrap(),0.3,0.8)){
		_h_JY_C1_Y2 ->fill(j.pT()/GeV, weight);
		}
		
		if(inRange(j.absrap(),0.8,1.2)){
		_h_JY_C1_Y3 ->fill(j.pT()/GeV, weight);
		}
		
		if(inRange(j.absrap(),1.2,2.1)){
		_h_JY_C1_Y4 ->fill(j.pT()/GeV, weight);
		}
	
	}
	
	}
	
	//Table 3
	
	foreach (const Jet& j, jets) {
		//Part1
		if(inRange(j.absrap(),0.0,2.1)){
			if (c >= 0 && c <= 10){	
			RAA_Y1_C1 ->fill(j.pT()/GeV,weight);
			}
			
			if (c >= 30 && c <= 40){
			RAA_Y1_C2 ->fill(j.pT()/GeV, weight);
			}
			
			if (c >= 60 && c <= 80){
			RAA_Y1_C3 ->fill(j.pT()/GeV, weight);
			}
		}
		//Part2
		if(inRange(j.absrap(),0.3,0.8)){
			if (c >= 0 && c <= 10){
			RAA_Y2_C1 ->fill(j.pT()/GeV, weight);
			}
			
			if (c >= 30 && c <= 40){
			RAA_Y2_C2 ->fill(j.pT()/GeV, weight);
			}
			
			if (c>= 60 && c <= 80){
			RAA_Y2_C3 ->fill(j.pT()/GeV, weight);
			}
		}
		//Part3
		if(inRange(j.absrap(),1.2,2.1)){
			if (c >= 0 && c <= 10){
			RAA_Y3_C1 ->fill(j.pT()/GeV, weight);
			}
			
			if (c >= 30 && c <= 40){
			RAA_Y3_C2 ->fill(j.pT()/GeV, weight);
			}
		
			if (c >= 60 && c <= 80){
			RAA_Y3_C3 ->fill(j.pT()/GeV, weight);
			}
		}
	}
	
	//Table 4
	
	foreach (const Jet& j, jets){
		if (j.pT()/GeV > 80 && j.pT()/GeV < 100 ){
			if (c >= 0 && c <= 10){
			RAA_C1 ->fill(j.absrap(), weight);
			}
			
			if (c >= 30 && c <= 40){
			RAA_C2 ->fill(j.absrap(), weight);
			}
			
			if (c >= 60 && c <= 80){
			RAA_C3 ->fill(j.absrap(), weight);
			}
		}
	}
	
	//Table 5
	
	foreach (const Jet& j, jets){
		if (j.pT()/GeV > 80 && j.pT()/GeV < 100 ){
			if(inRange(j.absrap(),0.0,2.1)){
			RAA_Y1 ->fill(NpartTot, weight);
			}
		}
	}
	
	
		
		
	} //end analysis section
	

    
    void finalize() {

      //normalize(_h_YYYY); // normalize to unity
	  
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }




    /// Name Histograms
    // Figure 1
    Histo1DPtr _h_DDCS_R4_Y1;
    Histo1DPtr _h_DDCS_R4_Y2;
    Histo1DPtr _h_DDCS_R4_Y3;
	Histo1DPtr _h_DDCS_R4_Y4;
	Histo1DPtr _h_DDCS_R4_Y5;
    
	//Figure 2
	Histo1DPtr _h_JY_Y21_C1; 
	Histo1DPtr _h_JY_Y21_C2; 
	Histo1DPtr _h_JY_Y21_C3; 
    Histo1DPtr _h_JY_Y21_C4;

	Histo1DPtr _h_JY_C1_Y1; 
	Histo1DPtr _h_JY_C1_Y2;
	Histo1DPtr _h_JY_C1_Y3;
	Histo1DPtr _h_JY_C1_Y4; 

	//Figure 3
	Histo1DPtr RAA_Y1_C1;
	Histo1DPtr RAA_Y1_C2;
	Histo1DPtr RAA_Y1_C3;
	  
	Histo1DPtr RAA_Y2_C1;
	Histo1DPtr RAA_Y2_C2;
	Histo1DPtr RAA_Y2_C3;
	  
	Histo1DPtr RAA_Y3_C1;
	Histo1DPtr RAA_Y3_C2;
	Histo1DPtr RAA_Y3_C3;
	
	//Figure 4
	Histo1DPtr RAA_C1;
	Histo1DPtr RAA_C2;
	Histo1DPtr RAA_C3;
	  
	Histo1DPtr RAA_Y1;
	
	//Counters
	int count1=0;
	int count2=0;
	int count3=0;
	int count4=0;
	
	int Npartt1 =0;
	int Npartp1 =0;
	int Npartt2 =0;
	int Npartp2 =0;
	int Npartt3 =0;
	int Npartp3 =0;
	int Npartt4 =0;
	int Npartp4 =0;
	
	
    
	
    


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1326911);


}
