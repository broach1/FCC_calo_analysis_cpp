#ifndef __HISTOGRAMCLASS_H__
#define __HISTOGRAMCLASS_H__

#include "TObject.h"
#include "TH1F.h"
#include "TString.h"

class HistogramClass {

 public:
  HistogramClass(double sf, double x0, double ENE, TString particle);
  ~HistogramClass();

  void Initialize_histos();
  void Delete_histos();
  void Reset_histos();

  TH1F* h_hitEnergy;
  TH1F* h_cellEnergy;
  TH1F* h_ptGen;
  TH1F* h_res_sf;
  TH1F* h_res_sf_25X0;
  TH1F* h_res_sf_27X0;
  TH1F* h_res_sf_30X0;
  TH1F* h_res_sf_35X0;
  TH1F* h_res_sf_40X0;
  TH1F* h_res_sf_45X0;
  TH1F* h_res_sf_50X0;
  TH1F* h_res_sf_55X0;
  TH1F* h_res_sf_60X0;

  std::vector<TH1F*> histVector;

 private:
  double SF;
  TString PARTICLE;
  double ENERGY;

};

#endif
