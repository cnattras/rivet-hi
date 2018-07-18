// -*- C++ -*-
#ifndef RIVET_ReferenceDataLoader_HH
#define RIVET_ReferenceDataLoader_HH

#include <fstream>
#include "YODA/ReaderYODA.h"

namespace Rivet {

  class ReferenceDataLoader {

  public:
    ReferenceDataLoader(std::string fileName){
		loadPPReferenceData(fileName);
	}
	
	void loadPPReferenceData(std::string fileName){
		vAO.clear();
		ifstream ppFile;
		ppFile.open(fileName);
		if (!ppFile.is_open()) { cout << fileName << ": DID NOT OPEN\n"; return; }
		
		YODA::Reader& rdr = YODA::ReaderYODA::create();
			
		rdr.read(ppFile,vAO);

		ppFile.close();
		
	}
	
	Histo1DPtr getPPReferenceHisto1D(std::string histName){
		YODA::Histo1D* hptr;
		for (unsigned int i = 0; i < vAO.size(); i++){
			if (histName==vAO[i]->path()){
				hptr = static_cast<YODA::Histo1D*>(vAO[i]);
				return Histo1DPtr(hptr);
				
			}
		}
		cout << "Not found\n";
		return 0x0;
		
	}

  private:

 
	
	/// vector for pp reference objects
	std::vector<AnalysisObject*> vAO;
  };
}

#endif
