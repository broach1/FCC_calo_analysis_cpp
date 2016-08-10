from ROOT import gSystem
gSystem.Load("libcaloanalysis-myanalysis")
from ROOT import CaloAnalysis_simple, TCanvas, TFile, TF1, gPad
#import os
#import numpy as np

PARTICLE = "e"
X0 = 8.25
ENERGY = 500
PART = [1,2,3,4,5,6,7,8,9,10] 
SUFFIX = ""
SF=18.44






for i in [PART]:

    #filename="root://eospublic//eos/fcc/users/n/novaj/July11_highStat/e"+str(ENERGY)+"_b0_LAr3mm_Lead2mm_rangeCut100mikrons_part"+str(i)+".root" 
    filename="root://eospublic//eos/fcc/users/b/broach/July11/e"+str(ENERGY)+"_n10000_lar"+str(LAR)+"_lead"+str(LEAD)+"_part"+str(PART)+".root"
    #filename="/afs/cern.ch/user/b/broach/FCCSW_updated/FCCSW/output.root"

    print "Processing file ",filename
    ma = CaloAnalysis_simple(SF, X0, ENERGY, PARTICLE)
    ma.loop(filename)
    print "Mean hit energy: ", ma.histClass.h_hitEnergy.GetMean()
    print "1/SF calculated: ", ENERGY/(ma.histClass.h_hitEnergy.GetMean())

    #c1 = TCanvas("c1","c1",1000,1000)
    #c1.Divide(2,2)
    #c1.cd(1)
    #ma.histClass.h_hitEnergy.Draw()
    #c1.cd(2)
    #ma.histClass.h_cellEnergy.Draw()
    #ma.histClass.h_cellEnergy.Fit("gaus")
    #c1.cd(3)
    #ma.histClass.h_ptGen.Draw()
    #c1.cd(4)
    #ma.histClass.h_res_sf.Draw()
    
    f2 = TFile("output-histo-"+PARTICLE+str(ENERGY)+"-b0-"+"lar"+str(LAR)+"-lead"+str(LEAD)+"-part"+str(i)+".root", "recreate")
    ma.histClass.h_res_sf.Write()
    ma.histClass.h_res_sf_25X0.Write()
    ma.histClass.h_res_sf_27X0.Write()
    ma.histClass.h_res_sf_30X0.Write()
    ma.histClass.h_res_sf_35X0.Write()
    ma.histClass.h_res_sf_40X0.Write()
    ma.histClass.h_res_sf_45X0.Write()
    ma.histClass.h_res_sf_50X0.Write()
    ma.histClass.h_res_sf_55X0.Write()
    ma.histClass.h_res_sf_60X0.Write()


#c1.SaveAs("plots_"+PARTICLE+str(ENERGY)+".gif")
#closeInput = raw_input("Press ENTER to exit") 
