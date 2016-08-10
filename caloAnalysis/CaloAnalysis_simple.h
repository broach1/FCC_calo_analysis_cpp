#ifndef __CALOANALYSIS_SIMPLE_H__
#define __CALOANALYSIS_SIMPLE_H__

#include "HistogramClass.h"

#include "TObject.h"
#include "TH1F.h"
#include "TString.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class CaloAnalysis_simple {

 public:
  CaloAnalysis_simple(const double sf, const double x0, const double ENE, const TString particle);
  ~CaloAnalysis_simple();

  void loop(const std::string filename);  //Open the file in the reader and loop through the events
  void processEvent(podio::EventStore& store, bool verbose,
		    podio::ROOTReader& reader);

  HistogramClass* histClass; 

 private:
  const double GeV=1000;
  double SF;               // 1/sampling_fraction
  double X0;
  TString PARTICLE;        // Particle type: e/mu
  double ENERGY;           // Beam energy
  double SumE_hit_ecal;    // Total hit energy per event
  double SumE_hit_ecal_25X0;
  double SumE_hit_ecal_27X0;
  double SumE_hit_ecal_30X0;
  double SumE_hit_ecal_35X0;
  double SumE_hit_ecal_40X0;
  double SumE_hit_ecal_45X0;
  double SumE_hit_ecal_50X0;
  double SumE_hit_ecal_55X0;
  double SumE_hit_ecal_60X0;

};

#endif
