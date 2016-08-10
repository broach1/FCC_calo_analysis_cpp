// ROOT
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector.h"

// STL
#include <vector>
#include <iostream>
#include <bitset>
#include <string>
#include <sstream>


//Material properties
// Values from http://pdg.lbl.gov/2015/AtomicNuclearProperties/

double rho_LAr = 1.396e3;
double X0_LAr = 14.0/100.;
double lambda_LAr = 85.77/100.;
double RM_LAr = 9.043/100.;
double mip_LAr = 2.105*100.;
double Ec_LAr = 32.84;

double rho_LKr = 2.418e3;
double X0_LKr = 4.703/100.;
double lambda_LKr = 61.8/100.;
double RM_LKr = 5.857/100.;
double mip_LKr = 3.281*100.;
double Ec_LKr = 17.03;

double rho_Sc = 1.032e3;
double X0_Sc = 42.4/100.;
double lambda_Sc = 81.9/1.032/100.;
double RM_Sc = 9.6/100.;
double mip_Sc = 1.936*1.032*100.;
double Ec_Sc = 94.0;

double rho_Si = 2.329e3;
double X0_Si = 9.37/100;
double lambda_Si = 46.52/100;
double RM_Si = 4.944/100;
double mip_Si = 3.876*100;
double Ec_Si = 40.19;

double rho_W = 19.3e3;
double X0_W = 0.3504/100;
double lambda_W = 9.946/100;
double RM_W = 0.9327/100;
double mip_W = 22.10*100;
double Ec_W = 7.97;
double Pr_W = 20.0;

double rho_Pb = 11.35e3;
double X0_Pb = 0.5612/100.0;
double lambda_Pb = 17.59/100.0;
double RM_Pb = 1.602/100.0;
double mip_Pb = 12.74*100.0;
double Ec_Pb = 7.43;
double Pr_Pb = 1.5;

double rho_Cu = 8.960e3;
double X0_Cu = 1.436/100;
double lambda_Cu = 15.32/100;
double RM_Cu = 1.568/100;
double mip_Cu = 12.57*100;
double Ec_Cu = 19.42;
double Pr_Cu = 4.0;

double rho_Fe = 7.874e3;
double X0_Fe = 1.757/100;
double lambda_Fe = 16.77/100;
double RM_Fe = 1.719/100;
double mip_Fe = 11.43*100;
double Ec_Fe = 21.68;
double Pr_Fe = 0.05;

double rho_Al = 2.699e3;
double X0_Al = 8.897/100;
double lambda_Al = 39.7/100;
double RM_Al = 4.419/100;
double mip_Al = 4.358*100;
double Ec_Al = 42.70;
double Pr_Al = 1.4;

double rho_G10 = 1.7e3;
double X0_G10 = 19.4/100;
double lambda_G10 = 90.2/1.7/100;
double mip_G10 = 1.87*1.7*100;
double Ec_G10 = 40.0;




double X0(std::vector<double> vec_W, double W0, std::vector<double> vec_X0){
  double inv_x0 = 0.0;
  for (unsigned i=0; i<vec_W.size(); i++) {
    inv_x0 += (vec_W.at(i)/W0) / vec_X0.at(i);
  }

  return 1/inv_x0;

}


double sig_over_E(double Ec, double E, double tabs, double F, double k, double Eloss_active){
  return TMath::Sqrt( (Ec*tabs/E)*(0.032*0.032/F + 9/(TMath::Power(TMath::Log(k*Eloss_active),2))));
}


double SamplFr(std::vector<double> aa, std::vector<double> mm){
  double dot = 0.0;

  for (unsigned i=0; i<aa.size(); i++){
    dot += aa.at(i)*mm.at(i);
  }

  return aa.at(aa.size()-1)*mm.at(mm.size()-1)/dot;
}

double AvX(std::vector<double> vec_dim, std::vector<double> vec_X0){
  double temp1 = 0.0;
  double temp2 = 0.0;
  for (unsigned i=0; i<vec_dim.size(); i++){
    temp1 += vec_dim.at(i);
    temp2 += vec_dim.at(i)/vec_X0.at(i);
      }
  return temp1/temp2;
}


double AvEc(std::vector<double> vec_dim, std::vector<double> vec_X0, std::vector<double> vec_Ec){
  double temp1 = 0.0;
  double temp2 = 0.0;
  for (unsigned i=0; i<vec_dim.size(); i++){
    temp1 += vec_dim.at(i)*vec_Ec.at(i)/vec_X0.at(i);
    temp2 += vec_dim.at(i)/vec_X0.at(i);
      }

  return temp1/temp2;
}

double ApplSamplTerm(double Ec, double tabs, double Elossact, double f, double kk){
  return TMath::Sqrt( Ec*tabs*(0.032*0.032/f + 0.09/TMath::Power(TMath::Log(kk*Elossact),2)));
}



Double_t fit3(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  std::cout<<mip_LAr<<std::endl;
  Double_t samplfr3 = 3.0*mip_LAr/(xx*mip_Pb+3.0*mip_LAr);
  return par[0]*TMath::Sqrt(3.0/samplfr3);
}

Double_t fit6(Double_t *x2, Double_t *par2){
  Float_t xx = x2[0];
  Double_t samplfr6 = 6.0*mip_LAr/(xx*mip_Pb+6.0*mip_LAr);
  return par2[0]*TMath::Sqrt(6.0/samplfr6);
}


//actually do the plotting
int martin_notebook(){

  //std::vector<double> dim {0.0015, 0.00212*2}; 
  std::vector<double> X0 {X0_Pb, X0_LAr};
  std::vector<double> Ec {Ec_Pb, Ec_LAr};
  double larx = 0.006;
  double larx2 = 0.012;
  double larx3 = 0.003;

  std::vector<double> vec_samp_term_6mm;
  std::vector<double> vec_samp_term_12mm;
  std::vector<double> vec_samp_term_3mm;
  std::vector<double> vec_pbx;

  int nPoints = 30;
  float pbxMin = 0.0;
  float pbxMax = 0.0065;


  for (double pbx=pbxMin; pbx<pbxMax; pbx += (pbxMax-pbxMin)/nPoints){

    std::vector<double> dim6 {pbx, larx};
    std::vector<double> dim12 {pbx, larx2};
    std::vector<double> dim3 {pbx, larx3};
    vec_samp_term_6mm.push_back( 100*ApplSamplTerm( AvEc(dim6, X0, Ec), pbx/X0_Pb, larx*mip_LAr, 1, 13e3));
    vec_samp_term_12mm.push_back( 100*ApplSamplTerm( AvEc(dim12, X0, Ec), pbx/X0_Pb, larx2*mip_LAr, 1, 13e3));
    vec_samp_term_3mm.push_back( 100*ApplSamplTerm( AvEc(dim3, X0, Ec), pbx/X0_Pb, larx3*mip_LAr, 1, 13e3));
    vec_pbx.push_back(pbx*1000);
  }

  double *ar_samp_term_6mm = &vec_samp_term_6mm[0];
  double *ar_samp_term_12mm = &vec_samp_term_12mm[0];
  double *ar_samp_term_3mm = &vec_samp_term_3mm[0];
  double *ar_pbx = &vec_pbx[0];
  

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();

  TLegend *leg1 = new TLegend(0.1,0.8,0.3,0.9);

  
  TGraph *g1 = new TGraph(vec_pbx.size(),ar_pbx,ar_samp_term_6mm);
  g1->Draw("AC");
  g1->GetYaxis()->SetRangeUser(0,20);
  g1->SetLineColor(kWhite);
  g1->SetLineWidth(3);
  g1->SetMarkerColor(kWhite);
  g1->SetLineStyle(7);
  g1->GetXaxis()->SetTitle("Pb layer (mm)");
  g1->GetYaxis()->SetTitle("Stochastic term [%]");
  g1->SetTitle("");
  //leg1->AddEntry(g1,"6mm LAr gap","l");
  
  
  /*
  TGraph *g2 = new TGraph(vec_pbx.size(),ar_pbx,ar_samp_term_12mm);
  g2->Draw("AC SAME");
  g2->SetLineColor(kBlue);
  g2->SetMarkerColor(kBlue);
  g2->SetLineWidth(3);
  g2->SetLineStyle(7);
  leg1->AddEntry(g2,"12mm LAr gap","l");
  */

  /*
  TGraph *g3 = new TGraph(vec_pbx.size(),ar_pbx,ar_samp_term_3mm);
  g3->Draw("SAME");
  g3->SetLineColor(kBlue);
  g3->SetMarkerColor(kBlue);
  g3->SetLineStyle(6);
  g3->SetLineWidth(3);
  leg1->AddEntry(g3,"3mm LAr gap","l");
  c1->Update();
  */

  std::vector<double> pb {2.0, 4.0, 6.0};
  std::vector<double> pb_err {0.0, 0.0, 0.0};
  std::vector<double> stoch_pb3 {8.03, 12.66, 16.29};
  std::vector<double> stoch_pb3err {0.0, 0.0, 0.0};
  std::vector<double> stoch_pb6 {6.53,10.67,13.93};
  std::vector<double> stoch_pb6err {0.0,0.0,0.0};
  
  double *ar_pb = &pb[0];
  double *ar_pb_err = &pb_err[0];
  double *ar_stoch_pb3 = &stoch_pb3[0];
  double *ar_stoch_pb3err = &stoch_pb3err[0];
  double *ar_stoch_pb6 = &stoch_pb6[0];
  double *ar_stoch_pb6err = &stoch_pb6err[0];

  TGraph *g_pb3_points = new TGraphErrors(pb.size(),ar_pb,ar_stoch_pb3,ar_pb_err, ar_stoch_pb3err);
  g_pb3_points->Draw("P");
  g_pb3_points->GetYaxis()->SetRangeUser(0,20);
  g_pb3_points->SetMarkerStyle(4);
  g_pb3_points->SetMarkerSize(1.0);
  g_pb3_points->SetMarkerColor(kBlue);

  
  TGraph *g_pb6_points = new TGraphErrors(pb.size(),ar_pb,ar_stoch_pb6,ar_pb_err,ar_stoch_pb6err);
  g_pb6_points->Draw("P SAME");
  g_pb6_points->SetMarkerStyle(5);
  g_pb6_points->SetMarkerSize(1.0);
  g_pb6_points->SetMarkerColor(kRed);
  c1->Update();
  
  

  TF1 *fit6 = new TF1("fit6","[0]*TMath::Sqrt(6.0+[1]*x/[2])",0,7);
  fit6->SetParameters(2.0,mip_Pb,mip_LAr);
  fit6->FixParameter(1,mip_Pb);
  fit6->FixParameter(2,mip_LAr);
  g_pb6_points->Fit("fit6");
  TF1 *res6 = g_pb6_points->GetFunction("fit6");
  res6->Draw("C SAME");
  res6->SetLineStyle(2);
  res6->SetLineColor(kRed);
  leg1->AddEntry(res6,"6mm LAr");

  TF1 *fit3 = new TF1("fit3","[0]*TMath::Sqrt(3.0+[1]*x/[2])",0,7);
  fit3->SetParameters(2.0,mip_Pb,mip_LAr);
  fit3->FixParameter(1,mip_Pb);
  fit3->FixParameter(2,mip_LAr);
  g_pb3_points->Fit("fit3");
  TF1* res3 = g_pb3_points->GetFunction("fit3");
  res3->Draw("C SAME");
  res3->SetLineStyle(2);
  res3->SetLineColor(kBlue);
  leg1->AddEntry(res3,"3mm LAr");


  leg1->AddEntry(g_pb6_points,"6mm LAr","P");

 
  leg1->AddEntry(g_pb3_points,"3mm LAr","P");
  
  
  leg1->Draw();
  c1->Update();


  
  return 0;

}
