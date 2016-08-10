#include "CaloAnalysis_simple.h"

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/CaloClusterCollection.h"

// ROOT
#include "TObject.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLorentzVector.h"

// STL
#include <vector>
#include <iostream>
#include <bitset>


CaloAnalysis_simple::CaloAnalysis_simple(const double sf, const double x0, const double ENE, const TString particle) 
{

  TH1::AddDirectory(kFALSE);

  SF = sf;
  X0 = x0;
  PARTICLE=particle;
  ENERGY = ENE;

  //Histograms initialization
  histClass = new HistogramClass(SF, X0, ENERGY, PARTICLE);
  histClass->Initialize_histos();
}


CaloAnalysis_simple::~CaloAnalysis_simple() {

  histClass->Delete_histos();
  delete histClass;

}

  

void CaloAnalysis_simple::loop(const std::string filename) {

  //Reset histograms
  histClass->Reset_histos();

  //Open file in the reader
  TString filename_eos;
  auto reader = podio::ROOTReader();
  auto store = podio::EventStore();
  try {
    //filename_eos =  "root://eospublic.cern.ch//eos/fcc/users/n/novaj/June10_ecalShifted/"+filename;
    reader.openFile(filename);
    std::cout << "CaloAnalysis_simple opening file " << filename << std::endl;
  }
  catch(std::runtime_error& err) {
    std::cerr<<err.what()<<". Quitting."<<std::endl;
    exit(1);
  }
  store.setReader(&reader);

  bool verbose = true;

  //Loop over all events
  unsigned nEvents = reader.getEntries();
  std::cout << "Number of events: " << nEvents << std::endl;
  for(unsigned i=0; i<nEvents; ++i) {
    if(i%50==0) std::cout<<"reading event "<<i<<std::endl;
    if(i>11) verbose = false;

    processEvent(store, verbose, reader);

    store.clear();
    reader.endOfEvent();
  }

  std::cout << "Total energy: " << histClass->h_cellEnergy->GetMean() << std::endl;
  std::cout << "End of loop" << std::endl;

  return;
}


void CaloAnalysis_simple::processEvent(podio::EventStore& store, bool verbose,
				podio::ROOTReader& reader) {

  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::CaloHitCollection*     colECalHit(nullptr);
  const fcc::CaloClusterCollection*     colECalCluster(nullptr);

 
  bool colMCParticlesOK = store.get("GenParticles", colMCParticles);
  bool colECalHitOK     = store.get("ECalHits" , colECalHit);
  bool colECalClusterOK     = store.get("ECalClusters" , colECalCluster);


  //Total hit energy per event
  SumE_hit_ecal = 0.;
  SumE_hit_ecal_25X0 = 0.;
  SumE_hit_ecal_27X0 = 0.;
  SumE_hit_ecal_30X0 = 0.;
  SumE_hit_ecal_35X0 = 0.;
  SumE_hit_ecal_40X0 = 0.;
  SumE_hit_ecal_45X0 = 0.;
  SumE_hit_ecal_50X0 = 0.;
  SumE_hit_ecal_55X0 = 0.;
  SumE_hit_ecal_60X0 = 0.;

  
  //Hit collection
  if (colECalHitOK && colECalHitOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #ECalClusters:     " << colECalCluster->size()    << std::endl;
    }
    //Loop through the collection
    for (auto& iecluster=colECalCluster->begin(); iecluster!=colECalCluster->end(); ++iecluster) 
        {
          //if (verbose) std::cout << "ECal hit energy " << iehit->Core().Energy << std::endl;
	  double rho = TMath::Sqrt( TMath::Power(iecluster->Core().position.X,2) + TMath::Power(iecluster->Core().position.Y,2) ) - 2700.0 ;
	  //offset of 2700 is to account for the detector starting at rho=2700mm

          SumE_hit_ecal += iecluster->Core().Energy;

	  
	  if (rho/X0 <= 25.){
	    SumE_hit_ecal_25X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_27X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_30X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_35X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_40X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_45X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;

	  }

	  else if (rho/X0 <= 27.){
	    SumE_hit_ecal_27X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_30X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_35X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_40X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_45X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;

	  }

	  else if (rho/X0 <= 30.){
	    SumE_hit_ecal_30X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_35X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_40X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_45X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;

	  }

	  else if (rho/X0 <= 35.){
	    SumE_hit_ecal_35X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_40X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_45X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;

	  }

	  else if (rho/X0 <= 40.){
	    SumE_hit_ecal_40X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_45X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;
	  }

	  else  if (rho/X0 <= 45.){
	    SumE_hit_ecal_45X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;
	  }

	  else if (rho/X0 <= 50.){
	    SumE_hit_ecal_50X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;
	  }

	  else if (rho/X0 <= 55.){
	    SumE_hit_ecal_55X0 += iecluster->Core().Energy;
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;
	  }

	  else if (rho/X0 <= 60.){
	    SumE_hit_ecal_60X0 += iecluster->Core().Energy;
	  }
	  


	}

    if (verbose) std::cout << "Total hit energy: " << SumE_hit_ecal << " hit collection size: " << colECalCluster->size() << std::endl;

    //Fill histograms
    histClass->h_hitEnergy->Fill(SumE_hit_ecal/GeV);
    histClass->h_cellEnergy->Fill(SumE_hit_ecal*SF/GeV);
    histClass->h_res_sf->Fill(SumE_hit_ecal*SF);
    histClass->h_res_sf_25X0->Fill(SumE_hit_ecal_25X0*SF);
    histClass->h_res_sf_27X0->Fill(SumE_hit_ecal_27X0*SF);
    histClass->h_res_sf_30X0->Fill(SumE_hit_ecal_30X0*SF);
    histClass->h_res_sf_35X0->Fill(SumE_hit_ecal_35X0*SF);
    histClass->h_res_sf_40X0->Fill(SumE_hit_ecal_40X0*SF);
    histClass->h_res_sf_45X0->Fill(SumE_hit_ecal_45X0*SF);
    histClass->h_res_sf_50X0->Fill(SumE_hit_ecal_50X0*SF);
    histClass->h_res_sf_55X0->Fill(SumE_hit_ecal_55X0*SF);
    histClass->h_res_sf_60X0->Fill(SumE_hit_ecal_60X0*SF);



  }
  else {
    if (verbose) {
      std::cout << "No CaloHit or CaloCluster Collection!!!!!" << std::endl;
    }
  }

 
  //MCParticle and Vertices collection 
  if (colMCParticlesOK) {
    if (verbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #MCTruthParticles:     " << colMCParticles->size()    << std::endl;
    }
    //Loop through the collection   
    for (auto& iparticle=colMCParticles->begin(); iparticle!=colMCParticles->end(); ++iparticle) {
      //Fill histogram
      histClass->h_ptGen->Fill( sqrt( pow(iparticle->Core().P4.Px,2)+
				      pow(iparticle->Core().P4.Py,2) ) );
    }
    
  }
  else {
    if (verbose) {
      std::cout << "No MCTruth info available" << std::endl;
    }
  }

}
