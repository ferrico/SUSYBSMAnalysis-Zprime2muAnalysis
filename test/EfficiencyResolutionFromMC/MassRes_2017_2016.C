#include "TFile.h"
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"

//     Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;



void MassRes_2017_2016(){

// 	int mrange[] = {120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
// 	Double_t m_bin[14] = {0};
// 	Double_t m_bin_err[14] = {0};
// 
// 	for(int i = 0; i < 14; i++){
// 		m_bin[i] = (mrange[i] + mrange[i+1]) / 2;
// 		m_bin_err[i] = (mrange[i+1] - mrange[i]) / 2;
// 	}
// 	
// 	Double_t BB_2017[] = {0.0104, 0.0121, 0.0164, 0.0177, 0.0245, 0.024, 0.0267, 0.0299, 0.0314, 0.0344, 0.0348, 0.038, 0.0401, 0.0421};
// 	Double_t BB_2016[] = {0.0117, 0.014, 0.0184, 0.0199, 0.0244, 0.0263, 0.0288, 0.0343, 0.037, 0.0395, 0.0438, 0.0503, 0.0527, 0.0606};
// 	Double_t BB_ratio[14] = {0};
// 	Double_t BB_ratio_err_2017[14] = {0.000113, 0.000139, 0.000227, 0.000181, 0.000325, 0.000246, 0.000405, 0.000306, 0.000383, 0.000341, 0.000367, 0.000345, 0.000418, 0.000328};
// 	Double_t BB_ratio_err_2016[14] = {0.000154, 0.000223, 0.000353, 0.000284, 0.000459, 0.000347, 0.00056, 0.000451, 0.000592, 0.000477, 0.000593, 0.000419, 0.000509, 0.000623};
// 	Double_t BB_ratio_err[14] = {0};
// 	
// 	Double_t BE_2017[] = {0.0153, 0.0163, 0.0199, 0.0214, 0.0259, 0.028, 0.031, 0.0344, 0.0357, 0.0394, 0.0415, 0.0451, 0.0464, 0.0497};
// 	Double_t BE_2016[] = {0.0176, 0.0194, 0.0229, 0.0255, 0.0318, 0.0329, 0.0354, 0.0397, 0.0426, 0.0465, 0.0496, 0.0559, 0.0574, 0.0625};
// 	Double_t BE_ratio[14] = {0};
// 	Double_t BE_ratio_err[14] = {0};
// 	Double_t BE_ratio_err_2017[14] = {0.000173, 0.000171, 0.000232, 0.000118, 0.000339, 0.000226, 0.000443, 0.000356, 0.000482, 0.000464, 0.000567, 0.000549, 0.000674, 0.000582};
// 	Double_t BE_ratio_err_2016[14] = {0.000247, 0.000288, 0.0004, 0.000321, 0.000548, 0.000416, 0.000665, 0.000541, 0.000766, 0.000717, 0.000854, 0.000664, 0.000787, 0.000977};




	int pt_BB[] = {72, 100, 152, 200, 300, 452, 800};
	Double_t pt_bin_BB[6] = {0};
	Double_t pt_bin_err_BB[6] = {0};
	int pt_BE[] = {52, 72, 100, 152, 200, 300, 452};
	Double_t pt_bin_BE[6] = {0};
	Double_t pt_bin_err_BE[6] = {0};


	for(int i = 0; i < 6; i++){
		pt_bin_BB[i] = (pt_BB[i] + pt_BB[i+1]) / 2;
		pt_bin_err_BB[i] = (pt_BB[i+1] - pt_BB[i]) / 2;
		pt_bin_BE[i] = (pt_BE[i] + pt_BE[i+1]) / 2;
		pt_bin_err_BE[i] = (pt_BE[i+1] - pt_BE[i]) / 2;
		std::cout<<pt_bin_BB[i]<<"\t"<<pt_bin_err_BB[i]<<"\t\t\t"<<pt_bin_BE[i]<<"\t"<<pt_bin_err_BE[i]<<std::endl;

	}
	
	Double_t BB_2017[] = {1.75597, 1.9056, 2.07152, 2.15898, 2.37409, 2.03149};
	Double_t BB_2016[] = {1.9, 2.07, 2.25, 2.34, 2.73, 3.26};
	Double_t BB_ratio[6] = {0};
	Double_t BB_ratio_err_2017[] = {0.0108987, 0.0156805, 0.0324627, 0.0470328, 0.0999323, 0.354159};
	Double_t BB_ratio_err_2016[] = {0.00924, 0.015, 0.0304, 0.0434, 0.108, 0.327};
	Double_t BB_ratio_err[6] = {0};
	
	Double_t BE_2017[] = {2.26876, 2.36699, 2.56533, 2.88632, 2.94238, 3.03503};
	Double_t BE_2016[] = {3.02, 3.29, 3.46, 3.84, 4.66, 5.48};
	Double_t BE_ratio[6] = {0};
	Double_t BE_ratio_err[6] = {0};
	Double_t BE_ratio_err_2017[] = {0.0102896, 0.015164, 0.0238536, 0.0556097, 0.0819892, 0.182202};
	Double_t BE_ratio_err_2016[] = {0.029, 0.0486, 0.0812, 0.222, 0.397, 1.5};


// 	for(int i = 0; i < 15; i++){
	for(int i = 0; i < 7; i++){
		BB_ratio[i] = (Double_t)(BB_2017[i] / BB_2016[i]) - 1;
		BE_ratio[i] = (Double_t)(BE_2017[i] / BE_2016[i]) - 1;

		Double_t BB_num =  BB_ratio_err_2017[i] / BE_2016[i];
		Double_t BB_den = BB_2017[i]/pow(BE_2016[i],2) * BB_ratio_err_2016[i];
		BB_ratio_err[i] = sqrt(pow(BB_num,2) + pow(BB_den,2));

		Double_t BE_num =  BE_ratio_err_2017[i] / BE_2016[i];
		Double_t BE_den = BE_2017[i]/pow(BE_2016[i],2) * BE_ratio_err_2016[i];
		BE_ratio_err[i] = sqrt(pow(BE_num,2) + pow(BE_den,2));
	}

		
	TGraphErrors* g_BB = new TGraphErrors(6, pt_bin_BB, BB_ratio, pt_bin_err_BB, BB_ratio_err);
	TGraphErrors* g_BE= new TGraphErrors(6, pt_bin_BE, BE_ratio, pt_bin_err_BE, BE_ratio_err);

// 	TGraphErrors* g_BB = new TGraphErrors(14, m_bin, BB_ratio, m_bin_err, BB_ratio_err);
// 	TGraphErrors* g_BE= new TGraphErrors(14, m_bin, BE_ratio, m_bin_err, BE_ratio_err);


	TCanvas* c_BB = new TCanvas("BB", "BB", 650, 500);
	c_BB->SetGrid();
	g_BB->SetMaximum(1);
	g_BB->SetMinimum(-1);
	g_BB->SetTitle("BB category");
	g_BB->SetMarkerColor(1);
	g_BB->SetMarkerStyle(21);		
	g_BB->Draw("AP");
	g_BB->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	g_BB->GetXaxis()->SetTitle("p_{T} [GeV]");
	g_BB->GetYaxis()->SetTitle("2017 / 2016 - 1");
// 	c_BB->Print("BB_ratio_2017_2016.png");
// 	c_BB->Print("BB_ratio_2017_2016.pdf");
	c_BB->Print("BB_ratio_2017_2016_AtZ.png");
	c_BB->Print("BB_ratio_2017_2016_AtZ.pdf");
	
	TCanvas* c_BE= new TCanvas("BE", "BE", 650, 500);
	c_BE->SetGrid();
	g_BE->SetMaximum(1);
	g_BE->SetMinimum(-1);
	g_BE->SetTitle("BE+EE category");
	g_BE->SetMarkerColor(1);
	g_BE->SetMarkerStyle(21);
	g_BE->Draw("AP");
	g_BE->GetYaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	g_BE->GetXaxis()->SetTitle("p_{T} [GeV]");
	g_BE->GetYaxis()->SetTitle("2017 / 2016 - 1");
// 	c_BE->Print("BE_ratio_2017_2016.png");
// 	c_BE->Print("BE_ratio_2017_2016.pdf");
	c_BE->Print("BE_ratio_2017_2016_AtZ.png");
	c_BE->Print("BE_ratio_2017_2016_AtZ.pdf");

}