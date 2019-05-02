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

#define Mass_Size 9
#define PtRes_Size 19

//     Int_t  binnum = sizeof(MASS_BINS)/sizeof(float)-1;



void MassRes_2018_2017(){

	TString reco = "";
// 	TString reco = "_Tracker";
// 	TString reco = "_Picky";
// 	TString reco = "_Dyt";
// 	TString reco = "_Std";


/////////////////		

// 		int mrange[] = {120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
		int mrange[] = {50, 120, 200, 400, 800, 1400, 2300, 3500, 4500, 6000};
		float m_bin[Mass_Size] = {0};
		float m_bin_err[Mass_Size] = {0};

		for(int i = 0; i < Mass_Size; i++){
			m_bin[i] = (mrange[i] + mrange[i+1]) / 2;
			m_bin_err[i] = (mrange[i+1] - mrange[i]) / 2;
		}

Double_t BB_NUM[] = {0.00780657, 0.00876999, 0.0115848, 0.0171445, 0.0240651, 0.0307671, 0.0339814, 0.0386358, 0.0425849};
Double_t BB_err_NUM[] = {0.000185573, 0.000170645, 0.000170586, 0.000208264, 0.000263631, 0.000313104, 0.000312262, 0.000342593, 0.0003783};
Double_t BE_NUM[] = {0.0107026, 0.0124983, 0.0154617, 0.0207812, 0.0282268, 0.0359416, 0.0406817, 0.0455708, 0.0466258};
Double_t BE_err_NUM[] = {0.000488789, 0.000219081, 0.000178841, 0.000210167, 0.000263481, 0.000351388, 0.00043464, 0.000558356, 0.000618974};
Double_t EE_NUM[] = {0.0156635, 0.0182229, 0.0223099, 0.027588, 0.0335764, 0.034552, 0.0439049, 0.0513625, 0.0580235};
Double_t EE_err_NUM[] = {0.000429611, 0.00044709, 0.000405259, 0.000453301, 0.000660653, 0.000760658, 0.000985207, 0.00114169, 0.00132394};
Double_t OVER_NUM[] = {0.0110402, 0.0118196, 0.0152128, 0.0204302, 0.0277198, 0.0328343, 0.0375256, 0.0428472, 0.0465368};
Double_t OVER_err_NUM[] = {0.000293469, 0.000243291, 0.000204374, 0.000246446, 0.000307298, 0.000354806, 0.000410498, 0.00049286, 0.000552801};

Double_t BB_DEN[] = {0.00884717, 0.0105042, 0.0135028, 0.0194735, 0.0272876, 0.034704, 0.0385993, 0.0469106, 0.0480378};
Double_t BB_err_DEN[] = {0.000154048, 1.18173, 0.000153767, 0.000199469, 0.000197655, 0.000266289, 0.000357622, 1.05777, 0.000415319};
Double_t BE_DEN[] = {0.0116208, 0.0137291, 0.0169957, 0.022824, 0.0314345, 0.0397074, 0.0458155, 0.0503469, 0.053158};
Double_t BE_err_DEN[] = {0.000421781, 0.000237252, 0.000239558, 0.0003224, 0.000375289, 0.000445728, 0.000555699, 0.000682419, 0.000696949};
Double_t EE_DEN[] = {0.0177427, 0.0203665, 0.024579, 0.0301027, 0.0341128, 0.0373804, 0.0495713, 0.0537728, 0.0628316};
Double_t EE_err_DEN[] = {0.000533759, 0.000619822, 0.000598978, 0.000709128, 0.00117921, 0.0011383, 0.00110793, 0.00177402, 0.00188263};
Double_t OVER_DEN[] = {0.0123626, 0.0130953, 0.0169431, 0.0230984, 0.0303363, 0.0366713, 0.0425724, 0.0472301, 0.0521332};
Double_t OVER_err_DEN[] = {0.000316706, 0.000328518, 0.000236, 0.000394387, 0.000437469, 0.000443543, 0.000471829, 0.000629841, 0.00068492};

		float BB_ratio[Mass_Size] = {0};
		float BB_ratio_err[Mass_Size] = {0};
		float BE_ratio[Mass_Size] = {0};
		float BE_ratio_err[Mass_Size] = {0};
		float OVER_ratio[Mass_Size] = {0};
		float OVER_ratio_err[Mass_Size] = {0};
		float EE_ratio[Mass_Size] = {0};
		float EE_ratio_err[Mass_Size] = {0};
/////////////////


/////////////////
/*		int pt_BB[] = {52, 72, 100, 152, 200, 300, 452, 800};
		float pt_bin_BB[7] = {0};
		float pt_bin_err_BB[7] = {0};
		for(int i = 0; i < 7; i++){
			pt_bin_BB[i] = (pt_BB[i] + pt_BB[i+1]) / 2;
			pt_bin_err_BB[i] = (pt_BB[i+1] - pt_BB[i]) / 2;
		}

		int pt_BE[] = {52, 72, 100, 152, 200, 300, 452};
		float pt_bin_BE[6] = {0};
		float pt_bin_err_BE[6] = {0};
		for(int i = 0; i < 6; i++){
			pt_bin_BE[i] = (pt_BE[i] + pt_BE[i+1]) / 2;
			pt_bin_err_BE[i] = (pt_BE[i+1] - pt_BE[i]) / 2;
		}

		int pt_EE[] = {52, 72, 100, 152, 200, 300, 452};
		float pt_bin_EE[6] = {0};
		float pt_bin_err_EE[6] = {0};
		for(int i = 0; i < 6; i++){
			pt_bin_EE[i] = (pt_EE[i] + pt_EE[i+1]) / 2;
			pt_bin_err_EE[i] = (pt_EE[i+1] - pt_EE[i]) / 2;
		}



		float BB_NUM[] = {1.64921, 1.75854, 1.91713, 1.89008, 2.31855, 2.55133, 1.56992};

		float BB_DEN[] = {1.6577, 1.71697, 1.82817, 1.91262, 2.15577, 2.40231, 3.38848};


		float BB_err_NUM[] = {0.020515, 0.0298464, 0.0449316, 0.0839773, 0.127892, 0.385846, 0.896881};

		float BB_err_DEN[] = {0.00959631, 0.013513, 0.0188022, 0.0368154, 0.0543271, 0.119831, 0.457507};



		float BE_NUM[] = {1.87074, 1.98607, 2.00033, 2.10491, 2.52572, 3.99623};

		float BE_DEN[] = {1.91781, 1.98912, 2.08325, 2.19315, 2.48865, 2.89399};


		float BE_err_NUM[] = {0.0208921, 0.0307667, 0.0440376, 0.093281, 0.160731, 0.660037};

		float BE_err_DEN[] = {0.00954004, 0.013927, 0.0202899, 0.0420597, 0.0697226, 0.176383};


		float EE_NUM[] = {2.42271, 2.58769, 2.81455, 2.86423, 3.08963, 3.10203};

		float EE_DEN[] = {2.4641, 2.63812, 2.79207, 3.04894, 3.32456, 3.76875};


		float EE_err_NUM[] = {0.0372571, 0.0554597, 0.0841685, 0.145758, 0.22425, 0.657327};
		float EE_err_DEN[] = {0.0173261, 0.025645, 0.0383928, 0.0846091, 0.119976, 0.316784};

		float BB_ratio[7] = {0};
		float BB_ratio_err[7] = {0};
		float BE_ratio[6] = {0};
		float BE_ratio_err[6] = {0};
		float EE_ratio[6] = {0};
		float EE_ratio_err[6] = {0};*/
/////////////////


/////////////////		
 /*   Double_t mrange[] = {50, 100, 150, 200, 250, 
	    					300, 350, 400, 450, 500, 
	    					600, 700, 800, 900, 1000, 1250, 1500, 2000, 3000};
		float pt_bin[PtRes_Size-1] = {0};
		float pt_bin_err[PtRes_Size-1] = {0};

		for(int i = 0; i < PtRes_Size-1; i++){
			pt_bin[i] = (mrange[i] + mrange[i+1]) / 2;
			pt_bin_err[i] = (mrange[i+1] - mrange[i]) / 2;
		}




Double_t BB_NUM[] = {0.0783423, 0.0845729, 0.0934588, 0.0942626, 0.0993533, 0.113318, 0.121638, 0.123528, 0.127515, 0.139632, 0.153719, 0.163579, 0.187521, 0.190447, 0.198128, 0.247937, 0.304542, 0.186841};
Double_t BB_err_NUM[] = {0.00125744, 0.00162694, 0.0023455, 0.00200631, 0.00287429, 0.00382489, 0.00330405, 0.00285478, 0.00369721, 0.00365886, 0.00366365, 0.00340155, 0.00488451, 0.00582155, 0.00350454, 0.00666853, 0.00423663, 0.00363002};
Double_t BE_NUM[] = {0.0958728, 0.100751, 0.10828, 0.118991, 0.12739, 0.13254, 0.139119, 0.147368, 0.156832, 0.166742, 0.169049, 0.17502, 0.181605, 0.19708, 0.209821, 0.228817, 0.370666, 0.18585};
Double_t BE_err_NUM[] = {0.00117844, 0.00155211, 0.0017582, 0.00223807, 0.00264215, 0.00267303, 0.0029728, 0.00342201, 0.00411785, 0.00325942, 0.00377025, 0.00465466, 0.00581056, 0.00681065, 0.00565257, 0.00986379, 0.0147792, 0.0455446};
Double_t EE_NUM[] = {0.138122, 0.15562, 0.167419, 0.167024, 0.179602, 0.182408, 0.180119, 0.196708, 0.206247, 0.198728, 0.195613, 0.241186, 0.230743, 0.274382, 0.247871, 0.374909, 0.015906, 3.93989};
Double_t EE_err_NUM[] = {0.00242102, 0.00345162, 0.0050135, 0.00466565, 0.0067063, 0.00869406, 0.00947879, 0.00930194, 0.011471, 0.0104046, 0.0149461, 0.0176595, 0.0241009, 0.0297487, 0.0330282, 0.14853, 0.103073, 567416};
Double_t OVER_NUM[] = {0.0952473, 0.10369, 0.113889, 0.118327, 0.130968, 0.135603, 0.137607, 0.155076, 0.144086, 0.165055, 0.177533, 0.18348, 0.183397, 0.194971, 0.202475, 0.216495, 0.320218, 0.251314};
Double_t OVER_err_NUM[] = {0.00130527, 0.00183657, 0.00230513, 0.00245689, 0.0033498, 0.00377291, 0.00359369, 0.00382344, 0.00434668, 0.00391671, 0.00438713, 0.00468404, 0.00580063, 0.00754376, 0.00508863, 0.00780788, 0.00782604, 0.0121476};

Double_t BB_DEN[] = {0.0791455, 0.0863892, 0.0979498, 0.102924, 0.105726, 0.110074, 0.127396, 0.13657, 0.147381, 0.154752, 0.172928, 0.191975, 0.209692, 0.219499, 0.213951, 0.228017, 0.351691, 0.200262};
Double_t BB_err_DEN[] = {0.00124669, 0.00169444, 0.00239188, 0.00226616, 0.00335203, 0.00401264, 0.00360911, 0.00336434, 0.00432373, 0.00425787, 0.00440065, 0.00425211, 0.00621615, 0.00724004, 0.00387801, 0.00631642, 0.00476495, 0.00413389};
Double_t BE_DEN[] = {0.0968976, 0.10684, 0.118558, 0.129238, 0.134507, 0.145072, 0.156886, 0.162265, 0.16571, 0.175403, 0.184137, 0.195779, 0.199705, 0.212367, 0.225916, 0.233228, 0.407743, -0.0634201};
Double_t BE_err_DEN[] = {0.00122962, 0.00165501, 0.00203607, 0.00254877, 0.0029159, 0.00315525, 0.00359554, 0.004046, 0.00459878, 0.00370259, 0.00431507, 0.00532539, 0.00644067, 0.00780422, 0.00633408, 0.0106735, 0.0196947, 0.0168905};
Double_t EE_DEN[] = {0.142029, 0.159074, 0.168101, 0.181215, 0.186431, 0.196828, 0.220337, 0.202223, 0.206223, 0.200595, 0.222142, 0.21427, 0.239753, 0.264133, 0.273702, 0.136969, 0.00769327, 0.277864};
Double_t EE_err_DEN[] = {0.00255908, 0.00356155, 0.00532286, 0.005289, 0.00738937, 0.00976464, 0.0115167, 0.00950422, 0.0118196, 0.010829, 0.0144528, 0.0140194, 0.0235941, 0.0293551, 0.0281774, 0.0312046, 0.00106075, 3.19157};
Double_t OVER_DEN[] = {0.0985483, 0.113068, 0.118193, 0.132326, 0.141431, 0.154789, 0.150475, 0.16853, 0.16806, 0.182042, 0.189523, 0.185644, 0.202908, 0.215196, 0.225078, 0.249871, 0.397109, 0.288049};
Double_t OVER_err_DEN[] = {0.00136445, 0.00198262, 0.00240835, 0.00277925, 0.00382256, 0.0046085, 0.00396091, 0.00438835, 0.00561301, 0.00473649, 0.00497853, 0.00538327, 0.00755534, 0.0081112, 0.00599332, 0.00967275, 0.0103532, 0.0153752};




		float BB_ratio[PtRes_Size] = {0};
		float BB_ratio_err[PtRes_Size] = {0};
		float BE_ratio[PtRes_Size] = {0};
		float BE_ratio_err[PtRes_Size] = {0};
		float OVER_ratio[PtRes_Size] = {0};
		float OVER_ratio_err[PtRes_Size] = {0};
		float EE_ratio[PtRes_Size] = {0};
		float EE_ratio_err[PtRes_Size] = {0};
*/
/////////////////

	int selection;
	int counter_max;
	
	selection = 0; // vs mass
// 	selection = 1; // vs pt
// 	selection = 2; // pt resolution	
	
	if(selection == 0) counter_max = Mass_Size;
	if(selection == 1) counter_max = 7;
	if(selection == 2) counter_max = PtRes_Size;

	for(int i = 0; i < counter_max; i++){
		BB_ratio[i] = (float)(BB_NUM[i] / BB_DEN[i]) - 1;
		float BB_num =  BB_err_NUM[i] / BB_DEN[i];
		float BB_den = BB_NUM[i]/pow(BB_DEN[i],2) * BB_err_DEN[i];
		BB_ratio_err[i] = sqrt(pow(BB_num,2) + pow(BB_den,2));
	}

	for(int i = 0; i < counter_max; i++){
		BE_ratio[i] = (float)(BE_NUM[i] / BE_DEN[i]) - 1;
		float BE_num =  BE_err_NUM[i] / BE_DEN[i];
		float BE_den = BE_NUM[i]/pow(BE_DEN[i],2) * BE_err_DEN[i];
		BE_ratio_err[i] = sqrt(pow(BE_num,2) + pow(BE_den,2));
	}

	for(int i = 0; i < counter_max; i++){
		OVER_ratio[i] = (float)(OVER_NUM[i] / OVER_DEN[i]) - 1;
		float OVER_num =  OVER_err_NUM[i] / OVER_DEN[i];
		float OVER_den = OVER_NUM[i]/pow(OVER_DEN[i],2) * OVER_err_DEN[i];
		OVER_ratio_err[i] = sqrt(pow(OVER_num,2) + pow(OVER_den,2));
	}
	
	for(int i = 0; i < counter_max; i++){
		EE_ratio[i] = (float)(EE_NUM[i] / EE_DEN[i]) - 1;
		float EE_num =  EE_err_NUM[i] / EE_DEN[i];
		float EE_den = EE_NUM[i]/pow(EE_DEN[i],2) * EE_err_DEN[i];
		EE_ratio_err[i] = sqrt(pow(EE_num,2) + pow(EE_den,2));
	}


	TString numerator = "Crujiff_MC" + reco;;
	TString denominator = "DSCB_MC" + reco;
	TString axis_name = numerator + " / " + denominator + " -1";
	TString save_name_BB_png, save_name_BB_pdf, save_name_BE_png, save_name_BE_pdf, save_name_EE_png, save_name_EE_pdf, save_name_OVER_png, save_name_OVER_pdf;
	if(selection == 0){
		save_name_BB_png = "BB_ratio_" + numerator + "_" + denominator + ".png";
		save_name_BB_pdf = "BB_ratio_" + numerator + "_" + denominator + ".pdf";
		save_name_BE_png = "BE_ratio_" + numerator + "_" + denominator + ".png";
		save_name_BE_pdf = "BE_ratio_" + numerator + "_" + denominator + ".pdf";
		save_name_OVER_png = "OVER_ratio_" + numerator + "_" + denominator + ".png";
		save_name_OVER_pdf = "OVER_ratio_" + numerator + "_" + denominator + ".pdf";
		save_name_EE_png = "EE_ratio_" + numerator + "_" + denominator + ".png";
		save_name_EE_pdf = "EE_ratio_" + numerator + "_" + denominator + ".pdf";
	}
	else if(selection == 1){
		save_name_BB_png = "BB_ratio_" + numerator + "_" + denominator + "._AtZ.png";
		save_name_BB_pdf = "BB_ratio_" + numerator + "_" + denominator + "._AtZ.pdf";
		save_name_BE_png = "BE_ratio_" + numerator + "_" + denominator + "._AtZ.png";
		save_name_BE_pdf = "BE_ratio_" + numerator + "_" + denominator + "._AtZ.pdf";
		save_name_OVER_png = "OVER_ratio_" + numerator + "_" + denominator + "._AtZ.png";
		save_name_OVER_pdf = "OVER_ratio_" + numerator + "_" + denominator + "._AtZ.pdf";
		save_name_EE_png = "EE_ratio_" + numerator + "_" + denominator + "._AtZ.png";
		save_name_EE_pdf = "EE_ratio_" + numerator + "_" + denominator + "._AtZ.pdf";
	}
	else{
		save_name_BB_png = "BB_ratio_" + numerator + "_" + denominator + "_PtRes.png";
		save_name_BB_pdf = "BB_ratio_" + numerator + "_" + denominator + "_PtRes.pdf";
		save_name_BE_png = "BE_ratio_" + numerator + "_" + denominator + "_PtRes.png";
		save_name_BE_pdf = "BE_ratio_" + numerator + "_" + denominator + "_PtRes.pdf";
		save_name_OVER_png = "OVER_ratio_" + numerator + "_" + denominator + "_PtRes.png";
		save_name_OVER_pdf = "OVER_ratio_" + numerator + "_" + denominator + "_PtRes.pdf";
		save_name_EE_png = "EE_ratio_" + numerator + "_" + denominator + "_PtRes.png";
		save_name_EE_pdf = "EE_ratio_" + numerator + "_" + denominator + "_PtRes.pdf";
	}
	
	TGraphErrors* g_BB = new TGraphErrors(Mass_Size, m_bin, BB_ratio, m_bin_err, BB_ratio_err);
	TGraphErrors* g_BE = new TGraphErrors(Mass_Size, m_bin, BE_ratio, m_bin_err, BE_ratio_err);
	TGraphErrors* g_OVER = new TGraphErrors(Mass_Size, m_bin, OVER_ratio, m_bin_err, OVER_ratio_err);
	TGraphErrors* g_EE = new TGraphErrors(Mass_Size, m_bin, EE_ratio, m_bin_err, EE_ratio_err);
// 	TGraphErrors* g_BB = new TGraphErrors(7, pt_bin_BB, BB_ratio, pt_bin_err_BB, BB_ratio_err);
// 	TGraphErrors* g_BE = new TGraphErrors(6, pt_bin_BE, BE_ratio, pt_bin_err_BE, BE_ratio_err);
// 	TGraphErrors* g_OVER = new TGraphErrors(6, pt_bin_OVER, OVER_ratio, pt_bin_err_OVER, OVER_ratio_err);
// 	TGraphErrors* g_EE = new TGraphErrors(6, pt_bin_EE, EE_ratio, pt_bin_err_EE, EE_ratio_err);
// 	TGraphErrors* g_BB = new TGraphErrors(PtRes_Size-1, pt_bin, BB_ratio, pt_bin_err, BB_ratio_err);
// 	TGraphErrors* g_BE = new TGraphErrors(PtRes_Size-1, pt_bin, BE_ratio, pt_bin_err, BE_ratio_err);
// 	TGraphErrors* g_OVER = new TGraphErrors(PtRes_Size-1, pt_bin, OVER_ratio, pt_bin_err, OVER_ratio_err);
// 	TGraphErrors* g_EE = new TGraphErrors(PtRes_Size-3, pt_bin, EE_ratio, pt_bin_err, EE_ratio_err);


	bool save = true;
// 	bool save = false;

	TCanvas* c_BB = new TCanvas("BB", "BB", 650, 500);
	c_BB->SetGrid();
	g_BB->SetMaximum(1);
	g_BB->SetMinimum(-1);
	g_BB->SetTitle("BB category");
	g_BB->SetMarkerColor(1);
	g_BB->SetMarkerStyle(21);		
	g_BB->Draw("AP");
	if(selection == 0) g_BB->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	else g_BB->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	g_BB->GetYaxis()->SetTitle("Zprime_Code_2017 / Santiago_2016 - 1");
	g_BB->GetYaxis()->SetTitle(axis_name);
	if(save) c_BB->Print("./Check_Resolution_Ratio/" + save_name_BB_png);
	if(save) c_BB->Print("./Check_Resolution_Ratio/" + save_name_BB_pdf);

	TCanvas* c_BE= new TCanvas("BE", "BE", 650, 500);
	c_BE->SetGrid();
	g_BE->SetMaximum(1);
	g_BE->SetMinimum(-1);
	g_BE->SetTitle("BE category");
	g_BE->SetMarkerColor(1);
	g_BE->SetMarkerStyle(21);
	g_BE->Draw("AP");
	if(selection == 0) g_BE->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	else g_BE->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	g_BE->GetYaxis()->SetTitle("Zprime_Code_2017 / Santiago_2016 - 1");
	g_BE->GetYaxis()->SetTitle(axis_name);	
	if(save) c_BB->Print("./Check_Resolution_Ratio/" + save_name_BE_png);
	if(save) c_BE->Print("./Check_Resolution_Ratio/" + save_name_BE_pdf);
	
	TCanvas* c_OVER= new TCanvas("OVER", "OVER", 650, 500);
	c_OVER->SetGrid();
	g_OVER->SetMaximum(1);
	g_OVER->SetMinimum(-1);
	g_OVER->SetTitle("OVER category");
	g_OVER->SetMarkerColor(1);
	g_OVER->SetMarkerStyle(21);
	g_OVER->Draw("AP");
	if(selection == 0) g_OVER->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	else g_OVER->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	g_OVER->GetYaxis()->SetTitle("Zprime_Code_2017 / Santiago_2016 - 1");
	g_OVER->GetYaxis()->SetTitle(axis_name);	
	if(save) c_BB->Print("./Check_Resolution_Ratio/" + save_name_OVER_png);
	if(save) c_OVER->Print("./Check_Resolution_Ratio/" + save_name_OVER_pdf);
	
	TCanvas* c_EE= new TCanvas("EE", "EE", 650, 500);
	c_EE->SetGrid();
	g_EE->SetMaximum(1);
	g_EE->SetMinimum(-1);
	g_EE->SetTitle("EE category");
	g_EE->SetMarkerColor(1);
	g_EE->SetMarkerStyle(21);
	g_EE->Draw("AP");
	if(selection == 0) g_EE->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	else g_EE->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	g_EE->GetYaxis()->SetTitle("Zprime_Code_2017 / Santiago_2016 - 1");
	g_EE->GetYaxis()->SetTitle(axis_name);	
	if(save) c_EE->Print("./Check_Resolution_Ratio/" + save_name_EE_png);
	if(save) c_EE->Print("./Check_Resolution_Ratio/" + save_name_EE_pdf);
	

}