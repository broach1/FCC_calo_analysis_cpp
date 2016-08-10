#include "HistogramClass.h"

// STL
#include <vector>
#include <iostream>
#include <bitset>

HistogramClass::HistogramClass(double sf, double x0, double ENE, TString particle)
{
  PARTICLE=particle;
  ENERGY = ENE;
}

HistogramClass::~HistogramClass()
{
}


void HistogramClass::Initialize_histos()
{

  h_hitEnergy = new TH1F("h_hitEnergy","h_hitEnergy", 200, 0, ENERGY);
  histVector.push_back(h_hitEnergy);

  if (PARTICLE=="e") {
    h_cellEnergy = new TH1F("h_cellenergy","h_cellEnergy", 100, ENERGY-0.2*ENERGY, ENERGY+0.2*ENERGY);
    h_res_sf = new TH1F("h_res_sf","h_res_SF", 1500e3, 0, 1500e3);
    h_res_sf_25X0 = new TH1F("h_res_sf_25X0","h_res_SF_25X0", 1500e3, 0, 1500e3);
    h_res_sf_27X0 = new TH1F("h_res_sf_27X0","h_res_SF_27X0", 1500e3, 0, 1500e3);

    h_res_sf_30X0 = new TH1F("h_res_sf_30X0","h_res_SF_30X0", 1500e3, 0, 1500e3);
    h_res_sf_35X0 = new TH1F("h_res_sf_35X0","h_res_SF_35X0", 1500e3, 0, 1500e3);
    h_res_sf_40X0 = new TH1F("h_res_sf_40X0","h_res_SF_40X0", 1500e3, 0, 1500e3);
    h_res_sf_45X0 = new TH1F("h_res_sf_45X0","h_res_SF_45X0", 1500e3, 0, 1500e3);
    h_res_sf_50X0 = new TH1F("h_res_sf_50X0","h_res_SF_50X0", 1500e3, 0, 1500e3);
    h_res_sf_55X0 = new TH1F("h_res_sf_55X0","h_res_SF_55X0", 1500e3, 0, 1500e3);
    h_res_sf_60X0 = new TH1F("h_res_sf_60X0","h_res_SF_60X0", 1500e3, 0, 1500e3);
  }
  else {
    if (PARTICLE=="mu") h_cellEnergy = new TH1F("h_cellEnergy","h_cellEnergy", 1000, 0, ENERGY-0.8*ENERGY);
    else {
      std::cout << "WARNING!!! Histogram ranges for " << PARTICLE << " particle not defined!!!" <<std::endl;
      h_cellEnergy = new TH1F("h_cellenergy","h_cellEnergy", 100, 0, ENERGY+0.2*ENERGY);
    }
  }
  histVector.push_back(h_cellEnergy);
  histVector.push_back(h_res_sf);
  histVector.push_back(h_res_sf_25X0);
  histVector.push_back(h_res_sf_27X0);
  histVector.push_back(h_res_sf_30X0);
  histVector.push_back(h_res_sf_35X0);
  histVector.push_back(h_res_sf_40X0);
  histVector.push_back(h_res_sf_45X0);
  histVector.push_back(h_res_sf_50X0);
  histVector.push_back(h_res_sf_55X0);
  histVector.push_back(h_res_sf_60X0);
  

  h_ptGen = new TH1F("h_ptGen","h_ptGen", 100, ENERGY-0.2*ENERGY, ENERGY+0.2*ENERGY);
  histVector.push_back(h_ptGen);

}


void HistogramClass::Reset_histos()
{

  for (auto iterator=histVector.begin(); iterator<histVector.end(); iterator++) {
    (*iterator)->Reset();
    (*iterator)->Sumw2();
  }

} 

void HistogramClass::Delete_histos()
{

  for (auto iterator=histVector.begin(); iterator<histVector.end(); iterator++) {
    delete (*iterator);
  }
  
  histVector.clear();

}
