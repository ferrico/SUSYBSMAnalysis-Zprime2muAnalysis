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




/////////////////		
/*		int mrange[] = {120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
		Double_t m_bin[14] = {0};
		Double_t m_bin_err[14] = {0};

		for(int i = 0; i < 14; i++){
			m_bin[i] = (mrange[i] + mrange[i+1]) / 2;
			m_bin_err[i] = (mrange[i+1] - mrange[i]) / 2;
		}

		
		Double_t BB_2017[] = 
		Double_t BB_2016[] =

		Double_t BB_err_2017[] =
		Double_t BB_err_2016[] =

		Double_t BE_2017[] =
		Double_t BE_2016[] =

		Double_t BE_err_2017[] =
		Double_t BE_err_2016[] =



			
		Double_t BB_ratio[14] = {0};
		Double_t BB_ratio_err[14] = {0};
		Double_t BE_ratio[14] = {0};
		Double_t BE_ratio_err[14] = {0};*/
/////////////////


/////////////////
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



		Double_t BB_2017[] = {1.73749, 1.87722, 2.05562, 2.19897, 2.74414, 2.20712};
		Double_t BB_2016[] = {1.72425, 1.87468, 2.03128, 2.30353, 2.6332, 3.12305};


		Double_t BB_err_2017[] = {0.0113565, 0.0166641, 0.0360212, 0.0510003, 0.131464, 0.947725};
		Double_t BB_err_2016[] = {0.0103617, 0.0149784, 0.0304594, 0.0473793, 0.121829, 0.775885};

		
		Double_t BE_2017[] = {2.22852, 2.32502, 2.54121, 2.73743, 3.16922, 3.32332};
		Double_t BE_2016[] = {2.18493, 2.28102, 2.46599, 2.79349, 2.93024, 3.05105};


		Double_t BE_err_2017[] = {0.010808, 0.0162842, 0.0250296, 0.0565059, 0.0980819, 0.233123};
		Double_t BE_err_2016[] = {0.00947867, 0.0140517, 0.0216485, 0.0511856, 0.0763764, 0.20141};


		Double_t BB_ratio[6] = {0};
		Double_t BB_ratio_err[6] = {0};
		Double_t BE_ratio[6] = {0};
		Double_t BE_ratio_err[6] = {0};
/////////////////
	
	

// 	for(int i = 0; i < 14; i++){
	for(int i = 0; i < 6; i++){
		BB_ratio[i] = (Double_t)(BB_2017[i] / BB_2016[i]) - 1;
		BE_ratio[i] = (Double_t)(BE_2017[i] / BE_2016[i]) - 1;

		Double_t BB_num =  BB_err_2017[i] / BE_2016[i];
		Double_t BB_den = BB_2017[i]/pow(BE_2016[i],2) * BB_err_2016[i];
		BB_ratio_err[i] = sqrt(pow(BB_num,2) + pow(BB_den,2));

		Double_t BE_num =  BE_err_2017[i] / BE_2016[i];
		Double_t BE_den = BE_2017[i]/pow(BE_2016[i],2) * BE_err_2016[i];
		BE_ratio_err[i] = sqrt(pow(BE_num,2) + pow(BE_den,2));
	}


	bool mass = false;
// 	bool mass = true;
	
	TString first = "2016_Tracker";
	TString second = "2017_94X_Tracker";
	TString axis_name = first + " / " + second + " -1";
	TString save_name_BB_png, save_name_BB_pdf, save_name_BE_png, save_name_BE_pdf;
	if(mass){
		save_name_BB_png = "BB_ratio_" + first + "_" + second + ".png";
		save_name_BB_pdf = "BB_ratio_" + first + "_" + second + ".pdf";
		save_name_BE_png = "BE_ratio_" + first + "_" + second + ".png";
		save_name_BE_pdf = "BE_ratio_" + first + "_" + second + ".pdf";
	}
	else{
		save_name_BB_png = "BB_ratio_" + first + "_" + second + "_AtZ.png";
		save_name_BB_pdf = "BB_ratio_" + first + "_" + second + "_AtZ.pdf";
		save_name_BE_png = "BE_ratio_" + first + "_" + second + "_AtZ.png";
		save_name_BE_pdf = "BE_ratio_" + first + "_" + second + "_AtZ.pdf";
	}

	std::cout<<axis_name<<std::endl;
	std::cout<<save_name_BB_png<<std::endl;
	std::cout<<save_name_BB_pdf<<std::endl;
	
// 	TGraphErrors* g_BB = new TGraphErrors(14, m_bin, BB_ratio, m_bin_err, BB_ratio_err);
// 	TGraphErrors* g_BE = new TGraphErrors(14, m_bin, BE_ratio, m_bin_err, BE_ratio_err);
	TGraphErrors* g_BB = new TGraphErrors(6, pt_bin_BB, BB_ratio, pt_bin_err_BB, BB_ratio_err);
	TGraphErrors* g_BE = new TGraphErrors(6, pt_bin_BE, BE_ratio, pt_bin_err_BE, BE_ratio_err);


	TCanvas* c_BB = new TCanvas("BB", "BB", 650, 500);
	c_BB->SetGrid();
	g_BB->SetMaximum(1);
	g_BB->SetMinimum(-1);
	g_BB->SetTitle("BB category");
	g_BB->SetMarkerColor(1);
	g_BB->SetMarkerStyle(21);		
	g_BB->Draw("AP");
	if(mass) g_BB->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	else g_BB->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	g_BB->GetYaxis()->SetTitle("Zprime_Code_2017 / Santiago_2016 - 1");
	g_BB->GetYaxis()->SetTitle(axis_name);
	c_BB->Print("./Check_Resolution_Code/" + save_name_BB_png);
	c_BB->Print("./Check_Resolution_Code/" + save_name_BB_pdf);

	TCanvas* c_BE= new TCanvas("BE", "BE", 650, 500);
	c_BE->SetGrid();
	g_BE->SetMaximum(1);
	g_BE->SetMinimum(-1);
	g_BE->SetTitle("BE+EE category");
	g_BE->SetMarkerColor(1);
	g_BE->SetMarkerStyle(21);
	g_BE->Draw("AP");
	if(mass) g_BE->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	else g_BE->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	g_BE->GetYaxis()->SetTitle("Zprime_Code_2017 / Santiago_2016 - 1");
	g_BE->GetYaxis()->SetTitle(axis_name);	
	c_BB->Print("./Check_Resolution_Code/" + save_name_BE_png);
	c_BE->Print("./Check_Resolution_Code/" + save_name_BE_pdf);

}