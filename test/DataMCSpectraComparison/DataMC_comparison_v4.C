#include "TFile.h"
#include <iostream>
#include <fstream>
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
#include <algorithm> 
#include <TString.h>

TString save_document;
TString name_histo;
void SalvaHisto(TString name, TH1F* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, float Y_min = 0, TString strig_1 = "MC", TString strig_2 = "DATA", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, float Y_min = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, TH2F* h_MC, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic = 1, bool logX = false);
void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA[15], TString save, int a, int b, int c, int d, int e);
void DrawStackPlot(THStack* dimuon_MC_Stack, TH1F* dimuon_MC_2016_clear[45], TH1F* dimuon_MC_2017_clear[16], TH1F* dimuon_MC_2018_clear[9]);

	int c_run = -1;
	int c_lumi = -1;
	int c_event = -1;
	float c_mass = -1;

	int n_run = -1;
	int n_lumi = -1;
	int n_event = -1;
	float n_mass = -1;

  UInt_t event;
  UInt_t run;
  Int_t prev_event = -88;
  UInt_t lumin;

  float dil_mass;
  
  float gen_lep_pt[2];
  float gen_lep_eta[2];
  float gen_lep_phi[2];
  float dil_pt;
  float cos_angle;
//   float vertex_chi2;
//   int dil_chosen;
  float gen_dil_mass;
  float genWeight;
//   int nvertices;
    bool GoodVtx;
    


// 
  float lep_pt[2];
  int lep_id[2];
  float lep_pt_err[2];
  float lep_eta[2];
  float lep_phi[2];

  float lep_tuneP_pt[2];
  float lep_glb_pt[2];
  float lep_picky_pt[2];
  float lep_tpfms_pt[2];
  float lep_dyt_pt[2];
  float lep_std_pt[2];
  float lep_tk_pt[2];
  float lep_cocktail_pt[2];
  
  float lep_tuneP_eta[2];
  float lep_glb_eta[2];
  float lep_picky_eta[2];
  float lep_tpfms_eta[2];
  float lep_dyt_eta[2];
  float lep_std_eta[2];
  float lep_tk_eta[2];
  
  float lep_tuneP_phi[2];
  float lep_glb_phi[2];
  float lep_picky_phi[2];
  float lep_tpfms_phi[2];
  float lep_dyt_phi[2];
  float lep_std_phi[2];
  float lep_tk_phi[2];

  float lep_dB[2];
  float lep_sumPt[2];
  float lep_pfIso[2];
//   float lep_triggerMatchPt[2];
  short lep_glb_numberOfValidTrackerLayers[2]; 
  short lep_glb_numberOfValidPixelHits[2];
  short lep_glb_numberOfValidMuonHits[2];
  short lep_TuneP_numberOfValidMuonHits[2];
  short lep_picky_numberOfValidMuonHits[2];
  short lep_dyt_numberOfValidMuonHits[2];
  short lep_tpfms_numberOfValidMuonHits[2];
  short lep_stanAlone_numberOfValidMuonHits[2];
  short lep_stanAlone_numberOfBadHits[2];
  short lep_stanAlone_numberOfMuonHits[2];
  short lep_numberOfMatchedStations[2];
  unsigned int lep_stationMask[2];
  short lep_numberOfMatchedRPCLayers[2];
  bool lep_isGlobalMuon[2];
  bool lep_isTrackerMuon[2];
//   float gen_lep_pt[2];
  float lep_triggerMatchPt[2];
  float lep_qOverPt[2];
  float vertex_chi2;
  
  float met_pt;
  
 Color_t color;
 float gM;
 float kFactor;
 float kFactor_BB;
 float kFactor_BE;
 
 TLegend * legend = new TLegend(0.85,0.65,0.99,0.95);

	TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

	
    TString samples[45] =  {
     						"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200",
                            "Wjets",
    						"dyInclusive50",
                            "Wantitop", "tW", 
                            "ttbar_lep50to500", "ttbar_lep_500to800", "ttbar_lep_800to1200", "ttbar_lep_1200to1800", "ttbar_lep1800toInf", 							
//                                "ttbar_50to500", "ttbar_500to800", "ttbar_800to1200", "ttbar_1200to1800", "ttbar_1800toInf",                                                      
                         "WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500",
                            "ZZ", "ZZ_ext", 
                            "WZ", "WZ_ext",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
	 						"dyPt0to50", "dyPt50to100", "dyPt100to250", "dyPt250to400", "dyPt400to650", "dyPt650",
                            };
	
	float events[45] = {
						6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
						29705748,
						19385554,  
						6933094, 6952830,
						79092400, 200000, 199800, 200000, 40829, 
						1999000, 200000, 200000, 200000, 38969, 
						990064, 998034, 
						1000000, 2995828,
						2977600, 100000, 100000, 98400, 100000, 95106, 100000, 100000, 100000,
						22782948, (float)39612900, (float)26998200, (float)7190820, 167272, 177101
// 						878212, 1662453, 2833172, 1571199, 48731, 177101					
						};
						
	float sigma[45] = {
						2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
						61526.7,
						6025.2, 
						35.6, 35.6,
						87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
						12.178, 1.385, 0.0566, 0.0035, 0.00005,
					    8.2615, 8.2615, 
					    23.565, 23.565, //16.523, 47.13, 16.523, 47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7,
						5352.58, 363.81, 84.015, 3.2283, 0.43604, 0.04098,
						};

    TString samples_2017[16] =  {
                            "dyInclusive50",
                            "Wantitop", "tW", 
//                             "ttbar",
							"ttbar_CP5",
                            "WW",
                            "ZZ",
                            "WZ",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
//                             "DY1JetsToLL_50_150", "DY1JetsToLL_150_250", "DY1JetsToLL_250_400", "DY1JetsToLL_400_Inf", "DY2JetsToLL_50_150", "DY2JetsToLL_150_250", "DY2JetsToLL_250_400", "DY2JetsToLL_400_Inf",
                            };
	
	float events_2017[16] = {
						18574774,
						7780870, 7581624,
// 						33844772,
						8705576,
						7791498,
						1949768,
						3928630,
						2961000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000,
					};
// 37987087					
// 20220257						
	float sigma_2017[16] = {
						5765.4,
						35.6, 35.6,
// 						831.76,
						72.1,
						118.7,
						16.523,
						47.13,
						2112.905, 20.553, 2.8861, 0.25126, 0.017075, 1.366E-3, 8.178E-5, 3.191E-6, 2.787E-7,

						};         
        
         
    TString samples_2018[9] =  {
//                             "Wantitop", "tW", 
//                             "ttbar",
//                             "WW",
//                             "ZZ",
//                             "WZ",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
                            };
	
	float events_2018[9] = {
// 						7780870, 7581624,
// 						33844772,
// 						7791498,
// 						1949768,
// 						3928630,
						2982000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000,
					};
						
	float sigma_2018[9] = {
// 						35.6, 35.6,
// 						831.76,
// 						118.7,
// 						16.523,
// 						47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7,
						};  

	float LUMINOSITY_2018 = 58542.894;
	float LUMINOSITY_2017 = 42416.950;
	float LUMINOSITY_2016 = 36235.493;
	float LUMINOSITY = LUMINOSITY_2018 + LUMINOSITY_2017 + LUMINOSITY_2016;
	
	float Z_peak = 0.9638;
	float Z_peak_BB = 0.9688;
	float Z_peak_BE = 0.9610;

	float weight[45] = {0};
	float weight_BB[45] = {0};
	float weight_BE[45] = {0};

	float weight_2017[16] = {0};
	float weight_2017_BB[16] = {0};
	float weight_2017_BE[16] = {0};

	float weight_2018[9] = {0};
	float weight_2018_BB[9] = {0};
	float weight_2018_BE[9] = {0};

	const int    NMBINS = 50;
	const double MMIN = 60., MMAX =3500.;
	double logMbins[NMBINS+1];

//     Double_t PT_BINS[] = {200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 600, 750, 1000, 1500};
//     Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;
    Double_t PT_BINS[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 
    					1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800};
    Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;    
    
    Double_t ETA_BINS[] = {-2.4, -1.5, -1.2, -0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9, 1.2, 1.5, 2.4};
//     Double_t ETA_BINS[] = {0, 0.9, 1.2, 1.5, 2.4};//, 5, 5.9, 6.2, 6.5, 7.4};
    Int_t  binnum_eta = sizeof(ETA_BINS)/sizeof(Double_t)-1;
    
//     Double_t PHI_BINS[] = {-3.14, -2.356, -1.57, -0.785, 0, 0.785, 1.57, 2.356, 3.14};
    Double_t PHI_BINS[] = {-3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2};
    Int_t  binnum_phi = sizeof(PHI_BINS)/sizeof(Double_t)-1;
    
    Double_t MET_BINS[] = {0, 25, 50, 75, 100, 200, 350, 500};//, 750, 1000};
    Int_t  binnum_met = sizeof(MET_BINS)/sizeof(Double_t)-1;
    
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 

	TH1F* pt_DATA = new TH1F("DATA_pt", "DATA_pt", binnum_pt, PT_BINS);
	TH1F* pt_DATA_BB = new TH1F("DATA_pt_BB", "DATA_pt_BB",binnum_pt, PT_BINS);
	TH1F* pt_DATA_BE = new TH1F("DATA_pt_BE", "DATA_pt_BE",binnum_pt, PT_BINS);
	TH1F* pt_DATA_plus = new TH1F("pt_DATA_plus", "pt_DATA_plus",binnum_pt, PT_BINS);
	TH1F* pt_DATA_minus = new TH1F("pt_DATA_minus", "pt_DATA_minus",binnum_pt, PT_BINS);
	TH1F* pt_DATA_RunBF = new TH1F("pt_DATA_RunBF", "pt_DATA_RunBF",binnum_pt, PT_BINS);
	TH1F* pt_DATA_RunGH = new TH1F("pt_DATA_RunGH", "pt_DATA_RunGH",binnum_pt, PT_BINS);

	THStack* pt_MC_Stack = new THStack("MC_pt", "MC_pt");
	THStack* pt_BB_MC_Stack = new THStack("MC_pt_BB", "MC_pt_BB");
	THStack* pt_BE_MC_Stack = new THStack("MC_pt_BE", "MC_pt_BE");
	THStack* pt_cumulative_MC_Stack = new THStack("MC_pt_cumulative", "MC_pt_cumulative");
	THStack* pt_cumulative_MC_BB_Stack = new THStack("MC_pt_cumulative_BB", "MC_pt_cumulative_BB");
	THStack* pt_cumulative_MC_BE_Stack = new THStack("MC_pt_cumulative_BE", "MC_pt_cumulative_BE");

	THStack* pt_MC_plus_Stack = new THStack("pt_MC_plus", "pt_MC_plus");
	THStack* pt_MC_minus_Stack = new THStack("pt_MC_minus", "pt_MC_minus");
	THStack* pt_MC_RunBF_Stack = new THStack("pt_MC_RunBF", "pt_MC_RunBF");
	THStack* pt_MC_RunGH_Stack = new THStack("pt_MC_RunGH", "pt_MC_RunGH");

	TH1F* pt_MC_2016_clear[45];
	TH1F* pt_MC_2017_clear[15];
	TH1F* pt_MC_2018_clear[9];
	TH1F* pt_cumulative_MC_2016_clear[45];
	TH1F* pt_cumulative_MC_2017_clear[15];
	TH1F* pt_cumulative_MC_2018_clear[9];
	TH1F* pt_BB_MC_2016_clear[45];
	TH1F* pt_BB_MC_2017_clear[15];
	TH1F* pt_BB_MC_2018_clear[9];
	TH1F* pt_BB_cumulative_MC_2016_clear[45];
	TH1F* pt_BB_cumulative_MC_2017_clear[15];
	TH1F* pt_BB_cumulative_MC_2018_clear[9];
	TH1F* pt_BE_MC_2016_clear[45];
	TH1F* pt_BE_MC_2017_clear[15];
	TH1F* pt_BE_MC_2018_clear[9];
	TH1F* pt_BE_cumulative_MC_2016_clear[45];
	TH1F* pt_BE_cumulative_MC_2017_clear[15];
	TH1F* pt_BE_cumulative_MC_2018_clear[9];
	TH1F* eta_MC_2016_clear[45];
	TH1F* eta_MC_2017_clear[15];
	TH1F* eta_MC_2018_clear[9];
	TH1F* phi_MC_2016_clear[45];
	TH1F* phi_MC_2017_clear[15];
	TH1F* phi_MC_2018_clear[9];
	TH1F* pt_MC_plus_clear;
	TH1F* pt_MC_minus_clear;
	TH1F* pt_MC_RunBF_clear;
	TH1F* pt_MC_RunGH_clear;

	TH1F* pt_cumulative_DATA = new TH1F("DATA_pt_cumulative", "DATA_pt_cumulative", binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_BB = new TH1F("DATA_pt_cumulative_BB", "DATA_pt_cumulative_BB",binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_BE = new TH1F("DATA_pt_cumulative_BE", "DATA_pt_cumulative_BE",binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_plus = new TH1F("DATA_pt_cumulative_plus", "DATA_pt_cumulative_plus",binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_minus = new TH1F("DATA_pt_cumulative_minus", "DATA_pt_cumulative_minus",binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_RunBF = new TH1F("DATA_pt_cumulative_RunBF", "DATA_pt_cumulative_RunBF",binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_RunGH = new TH1F("DATA_pt_cumulative_RunGH", "DATA_pt_cumulative_RunGH",binnum_pt, PT_BINS);

	THStack* pt_cumulative_MC_plus_Stack = new THStack("MC_pt_cumulative_plus", "MC_pt_cumulative_plus");
	THStack* pt_cumulative_MC_minus_Stack = new THStack("MC_pt_cumulative_minus", "MC_pt_cumulative_minus");
	THStack* pt_cumulative_MC_RunBF_Stack = new THStack("MC_pt_cumulative_RunBF", "MC_pt_cumulative_RunBF");
	THStack* pt_cumulative_MC_RunGH_Stack = new THStack("MC_pt_cumulative_RunGH", "MC_pt_cumulative_RunGH");

	TH1F* pt_cumulative_MC_plus_clear;
	TH1F* pt_cumulative_MC_minus_clear;
	TH1F* pt_cumulative_MC_RunBF_clear;
	TH1F* pt_cumulative_MC_RunGH_clear;
	
	TH1F* eta_DATA = new TH1F("DATA_eta", "DATA_eta", binnum_eta, ETA_BINS);
	THStack* eta_MC_Stack = new THStack("MC_eta", "MC_eta");

	TH1F* phi_DATA = new TH1F("DATA_phi", "DATA_phi", binnum_phi, PHI_BINS);
	THStack* phi_MC_Stack = new THStack("MC_phi", "MC_phi");

	TH1F* dB_DATA = new TH1F("DATA_dB", "DATA_dB",  100, 0.005, 0.5);
	THStack* dB_MC_Stack = new THStack("MC_dB", "MC_dB");
	TH1F* dB_MC_clear;

	TH1F* PixelHit_DATA = new TH1F("DATA_PixelHit", "DATA_PixelHit",  15, 0,  15);
	THStack* PixelHit_MC_Stack = new THStack("MC_PixelHit", "MC_PixelHit");
	TH1F* PixelHit_MC_clear;

	TH1F* TkLayer_DATA = new TH1F("DATA_TkLayer", "DATA_TkLayer", 20, 0, 20);
	THStack* TkLayer_MC_Stack = new THStack("MC_TkLayer", "MC_TkLayer");
	TH1F* TkLayer_MC_clear;

	TH1F* Iso_DATA = new TH1F("DATA_Iso", "DATA_Iso", 60, 0, 0.3);
	THStack* Iso_MC_Stack = new THStack("MC_Iso", "MC_Iso");
	TH1F* Iso_MC_clear;

	TH1F* relpTErr_DATA = new TH1F("DATA_relpTErr", "DATA_relpTErr", 50, 10e-3, 0.5);
	THStack* relpTErr_MC_Stack = new THStack("MC_relpTErr", "MC_relpTErr");
	TH1F* relpTErr_MC_clear;

	TH1F* Vtx_DATA = new TH1F("DATA_Vtx", "DATA_Vtx",  60, 0, 30);
	THStack* Vtx_MC_Stack = new THStack("MC_Vtx", "MC_Vtx");
	TH1F* Vtx_MC_clear;

	TH1F* ValidMu_DATA = new TH1F("DATA_ValidMu", "DATA_ValidMu",  55, 0, 55);
	THStack* ValidMu_MC_Stack = new THStack("MC_ValidMu", "MC_ValidMu");
	TH1F* ValidMu_MC_clear;
	
void DataMC_comparison_v4(){

	gStyle->SetOptFit(1);   
	gStyle->SetOptStat(0);        
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "XYZ");
	gStyle->SetLabelSize(0.04, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.06, "XYZ");
	gStyle->SetPadBorderMode(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineColor(0);
	gStyle->SetPaintTextFormat(".3f");
//     gStyle->SetFillColor(0);
//     gStyle->SetPadTopMargin(0.05);
//     gStyle->SetPadBottomMargin(0.13);
//     gStyle->SetPadLeftMargin(0.16);
//     gStyle->SetPadRightMargin(0.05);

	std::vector<TString> variables;
	variables.push_back("mass");
	variables.push_back("pt");
	variables.push_back("pt_bb");
	variables.push_back("pt_be");
	variables.push_back("eta");
	variables.push_back("phi");
	variables.push_back("dB");
	variables.push_back("PixelHi");
	variables.push_back("TkLaye");
	variables.push_back("Iso");
	variables.push_back("relpTErr");
	variables.push_back("Vtx");
	variables.push_back("ValidMu");

	for (int ibin = 0; ibin <= NMBINS; ibin++){
    	logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
    	if(ibin != 0 && logMbins[ibin] <= logMbins[ibin-1]) std::cout<<"PROBLEMA"<<std::endl;
    	std::cout<<logMbins[ibin]<<",";
    }

// 	std::vector<Int_t> nBins;
// 	nBins.push_back(NMBINS);
// 	nBins.push_back(binnum_pt);
// 	nBins.push_back(binnum_pt);
// 	nBins.push_back(binnum_pt);
// 	nBins.push_back(binnum_eta);
// 	nBins.push_back(binnum_phi);
// 	nBins.push_back(100);
// 	nBins.push_back(15);
// 	nBins.push_back(20);
// 	nBins.push_back(60);
// 	nBins.push_back(50);
// 	nBins.push_back(60);
// 	nBins.push_back(55);
// 	
// 	std::vector<Double_t> xMin;
// 	xMin.push_back(0.005);
// 	xMin.push_back(0);
// 	xMin.push_back(0);
// 	xMin.push_back(0);
// 	xMin.push_back(10e-3);
// 	xMin.push_back(0);
// 	xMin.push_back(0);
// 
// 	std::vector<Double_t> xMax;
// 	xMax.push_back(&logMbins[0]);
// 	xMax.push_back(&PT_BINS[0]);
// 	xMax.push_back(&PT_BINS[0]);
// 	xMax.push_back(&PT_BINS[0]);
// 	xMax.push_back(&ETA_BINS[0]);
// 	xMax.push_back(&PHI_BINS[0]);
// 	xMax.push_back(0.5);
// 	xMax.push_back(15);
// 	xMax.push_back(20);
// 	xMax.push_back(0.3);
// 	xMax.push_back(0.5);
// 	xMax.push_back(30);
// 	xMax.push_back(55);
// 
// 	std::vector<TString> type;
// 	type.push_back("DATA");
// 	type.push_back("MC");
// 	type.push_back("GEN");
// 
// 	std::vector<STring> year;
// 	year.push_back("2016");
// 	year.push_back("2017");
// 	year.push_back("2018");
// 
// 
//     
//     for(int v = 0; v < variables.size(); v++){
// 	    for(int t = 0; t < type.size(); t++){
// 			if(t == 0){
// 				TString NOME = variables.at(v) + "_" + type.at(t);
// 		    	if(v < 7){
// 					TH1F* NOME = new TH1F(NOME, NOME, nBins.at(v), xMax.at(v));
// 				}
// 				else{
// 					TH1F* NOME = new TH1F(NOME, NOME, nBins.at(v), xMin.at(v-7), xMax.at(v-7));				
// 				}
// 			}
// 			else{
// 				for(int y = 0; y < year.size(); y++){
// 					TString NOME = variables.at(v) + "_" + type.at(t) + "_" + year.at(y) + "_clear";
// 					if(y == 0) TH1F* NOME[45];// = new TH1F(NOME, NOME, nBins.at(v), nMax.at(v));
// 					if(y == 1) TH1F* NOME[16];
// 					if(y == 2) TH1F* NOME[9];			
// 
// 					TString NOME_Stack = variables.at(v) + "_" + type.at(t) + "_" + year.at(y) + "_Stack";
// 					THStack* NOME_Stack;
// 				}
// 			}
// 	    }    
//     }

	TH1F* dimuon_DATA = new TH1F("DATA_dimuon", "DATA_dimuon", NMBINS, logMbins);
	THStack* dimuon_MC_Stack = new THStack("MC_dimuon", "MC_dimuon");
	TH1F* dimuon_MC_2016_clear[45];
	TH1F* dimuon_MC_2017_clear[16];
	TH1F* dimuon_MC_2018_clear[9];

	TH1F* dimuon_BB_DATA = new TH1F("DATA_dimuon_BB", "DATA_dimuon_BB", NMBINS, logMbins);
	THStack* dimuon_BB_MC_Stack = new THStack("MC_dimuon_BB", "MC_dimuon_BB");
	TH1F* dimuon_BB_MC_2016_clear[45];
	TH1F* dimuon_BB_MC_2017_clear[15];
	TH1F* dimuon_BB_MC_2018_clear[9];

	TH1F* dimuon_BE_DATA = new TH1F("DATA_dimuon_BE", "DATA_dimuon_BE", NMBINS, logMbins);
	THStack* dimuon_BE_MC_Stack = new THStack("MC_dimuon_BE", "MC_dimuon_BE");
	TH1F* dimuon_BE_MC_2016_clear[45];
	TH1F* dimuon_BE_MC_2017_clear[15];
	TH1F* dimuon_BE_MC_2018_clear[9];

	TH1F* dimuon_cumulative_DATA = new TH1F("DATA_dimuon_cumulative", "DATA_dimuon_cumulative", NMBINS, logMbins);
	THStack* dimuon_cumulative_MC_Stack = new THStack("MC_dimuon_cumulative", "MC_dimuon_cumulative");
	TH1F* dimuon_cumulative_MC_2016_clear[45];
	TH1F* dimuon_cumulative_MC_2017_clear[15];
	TH1F* dimuon_cumulative_MC_2018_clear[9];

	TH1F* dimuon_cumulative_BB_DATA = new TH1F("DATA_dimuon_cumulative_BB", "DATA_dimuon_cumulative_BB", NMBINS, logMbins);
	THStack* dimuon_cumulative_BB_MC_Stack = new THStack("MC_dimuon_cumulative_BB", "MC_dimuon_cumulative_BB");
	TH1F* dimuon_BB_cumulative_MC_2016_clear[45];
	TH1F* dimuon_BB_cumulative_MC_2017_clear[15];
	TH1F* dimuon_BB_cumulative_MC_2018_clear[9];

	TH1F* dimuon_cumulative_BE_DATA = new TH1F("DATA_dimuon_cumulative_BE", "DATA_dimuon_cumulative_BE", NMBINS, logMbins);
	THStack* dimuon_cumulative_BE_MC_Stack = new THStack("MC_dimuon_cumulative_BE", "MC_dimuon_cumulative_BE");
	TH1F* dimuon_BE_cumulative_MC_2016_clear[45];
	TH1F* dimuon_BE_cumulative_MC_2017_clear[15];
	TH1F* dimuon_BE_cumulative_MC_2018_clear[9];

	legend->AddEntry(dimuon_DATA, "DATA", "lep");

// 	LUMINOSITY_2017 += LUMINOSITY_2018;
// 	LUMINOSITY_2016 += LUMINOSITY_2018;

	for(int i = 0; i < 45; i++){
	
		weight[i] = LUMINOSITY_2016 * sigma[i] / events[i];
		weight[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
		weight_BB[i] = weight[i] * Z_peak_BB / Z_peak;
		weight_BE[i] = weight[i] * Z_peak_BE/ Z_peak;
	}

	for(int i = 0; i < 16; i++){
	
		weight_2017[i] = LUMINOSITY_2017 * sigma_2017[i] / events_2017[i];
		weight_2017[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
		weight_2017_BB[i] = weight_2017[i] * Z_peak_BB / Z_peak;
		weight_2017_BE[i] = weight_2017[i] * Z_peak_BE/ Z_peak;
	}	

	for(int i = 0; i < 9; i++){
	
		weight_2018[i] = LUMINOSITY_2018 * sigma_2018[i] / events_2018[i];
		weight_2018[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
		weight_2018_BB[i] = weight_2018[i] * Z_peak_BB / Z_peak;
		weight_2018_BE[i] = weight_2018[i] * Z_peak_BE/ Z_peak;

	}	
	
         float ciao_bhu[51] = {60,65.0833,70.5973,76.5784,83.0662,90.1038,
         97.7375,106.018,115,124.743,135.311,146.775,159.21,
         172.699,187.33,203.201,220.417,239.091,259.347,281.319,
         305.153,331.006,359.05,389.469,422.466,458.258,497.082,539.196,
         584.877,634.429,688.179,746.483,809.726,878.327,952.741,1033.46,
         1121.01,1215.99,1319.01,1430.76,1551.97,1683.46,1826.09,1980.8,
         2148.61,2330.65,2528.1,2742.29,2974.62,3226.63,3500};
	
	for(int j = 13; j < 45; j++){
		
      	if(j >= 39) continue; 
//       	if(j < 30) continue;
//       	if(j >= 30) continue;
//       	if(j >= 30 and j <= 38) continue;
//       	if(j == 30) continue; 
//       	if(j == 44) continue; 
//       	if(j >= 21) continue;
//       	if(j < 21) continue;

        std::cout<<"opening.. "<<samples[j]<<" --- "<<events[j]<<std::endl;

        TChain *treeMC = new TChain("SimpleNtupler/t");

        treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/ana_datamc_" + samples[j] + ".root");        

	    Long64_t nentries = treeMC->GetEntries();
	    
    	treeMC->SetBranchAddress("genWeight",&genWeight);
     	    	
      	treeMC->SetBranchAddress("event", &event);
      	treeMC->SetBranchAddress("run", &run);
//       	treeMC->SetBranchAddress("lumi", &lumi);    	
    	treeMC->SetBranchAddress("dil_mass",&dil_mass);
    	
    	treeMC->SetBranchAddress("gen_lep_pt", gen_lep_pt);
    	treeMC->SetBranchAddress("gen_lep_eta", gen_lep_eta);
    	treeMC->SetBranchAddress("gen_lep_phi", gen_lep_phi);
    	treeMC->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
	     treeMC->SetBranchAddress("dil_pt",&dil_pt);
	     treeMC->SetBranchAddress("cos_angle",&cos_angle);
 	    treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeMC->SetBranchAddress("nvertices",&nvertices);
	     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);

	     treeMC->SetBranchAddress("lep_pt",lep_pt);
    	 treeMC->SetBranchAddress("lep_id",lep_id);
	     treeMC->SetBranchAddress("lep_eta",lep_eta);
	     treeMC->SetBranchAddress("lep_phi",lep_phi);
	     treeMC->SetBranchAddress("lep_dB",lep_dB);
	     treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	     treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    	 treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	     treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
	     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
	     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
	     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    	 treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
	     treeMC->SetBranchAddress("lep_pfIso",lep_pfIso);
	     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
    	 treeMC->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    	 treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
    	 treeMC->SetBranchAddress("lep_qOverPt",lep_qOverPt);
    	 

     //treeMC->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeMC->SetBranchAddress("met_pt",&met_pt); 	
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples[j].Data(),j , weight[j], nentries);
         
         dimuon_MC_2016_clear[j] = new TH1F(samples[j] + "MC_dimuon_2016", samples[j] + "MC_dimuon_2016", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BB_MC_2016_clear[j] = new TH1F(samples[j] + "MC_BB_dimuon_2016", samples[j] + "MC_BB_dimuon_2016", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BE_MC_2016_clear[j] = new TH1F(samples[j] + "MC_BE_dimuon_2016", samples[j] + "MC_BE_dimuon_2016", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_cumulative_MC_2016_clear[j] = new TH1F(samples[j] + "MC_dimuon_cumulative_2016", samples[j] + "MC_dimuon_cumulative_2016", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BB_cumulative_MC_2016_clear[j] = new TH1F(samples[j] + "MC_dimuon_BB_cumulative_2016", samples[j] + "MC_dimuon_BB_cumulative_2016", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BE_cumulative_MC_2016_clear[j] = new TH1F(samples[j] + "MC_dimuon_BE_cumulative_2016", samples[j] + "MC_dimuon_BE_cumulative_2016", 50, ciao_bhu);//NMBINS, logMbins);         
		pt_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_2016", samples[j] + "_Lepton_Pt_2016", binnum_pt, PT_BINS);
		pt_BB_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BB_2016", samples[j] + "_Lepton_Pt_BB_2016", binnum_pt, PT_BINS);
		pt_BE_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BE_2016", samples[j] + "_Lepton_Pt_BE_2016", binnum_pt, PT_BINS);
		pt_cumulative_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_cumulative_2016", samples[j] + "_Lepton_Pt_cumulative_2016", binnum_pt, PT_BINS);
		pt_BB_cumulative_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BB_cumulative_2016", samples[j] + "_Lepton_Pt_BB_cumulative_2016", binnum_pt, PT_BINS);
		pt_BE_cumulative_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BE_cumulative_2016", samples[j] + "_Lepton_Pt_BE_cumulative_2016", binnum_pt, PT_BINS);
		eta_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Eta_2016", samples[j] + "_Lepton_Eta_2016", binnum_eta, ETA_BINS);
		phi_MC_2016_clear[j] = new TH1F(samples[j] + "_Lepton_Phi_2016", samples[j] + "_Lepton_Phi_2016", binnum_phi, PHI_BINS);
		pt_MC_plus_clear = new TH1F("MC_pt_plus", "MC_pt_plus", binnum_pt, PT_BINS);
		pt_MC_minus_clear = new TH1F("MC_pt_minus", "MC_pt_minus", binnum_pt, PT_BINS);
		pt_MC_RunBF_clear = new TH1F("MC_pt_RunBF", "MC_pt_RunBF", binnum_pt, PT_BINS);
		pt_MC_RunGH_clear = new TH1F("MC_pt_RunGH", "MC_pt_RunGH", binnum_pt, PT_BINS);

		pt_cumulative_MC_plus_clear = new TH1F("MC_pt_cumulative_plus", "MC_pt_cumulative_plus", binnum_pt, PT_BINS);
		pt_cumulative_MC_minus_clear = new TH1F("MC_pt_cumulative_minus", "MC_pt_cumulative_minus", binnum_pt, PT_BINS);
		pt_cumulative_MC_RunBF_clear = new TH1F("MC_pt_cumulative_RunBF", "MC_pt_cumulative_RunBF", binnum_pt, PT_BINS);
		pt_cumulative_MC_RunGH_clear = new TH1F("MC_pt_cumulative_RunGH", "MC_pt_cumulative_RunGH", binnum_pt, PT_BINS);

			
		dB_MC_clear = new TH1F("MC_dB", "MC_dB",  100, 0.005, 0.5);

		PixelHit_MC_clear = new TH1F("MC_PixelHit", "MC_PixelHit",  15, 0,  15);

		TkLayer_MC_clear = new TH1F("MC_TkLayer", "MC_TkLayer", 20, 0, 20);

		Iso_MC_clear = new TH1F("MC_Iso", "MC_Iso", 60, 0.0, 0.3);

		relpTErr_MC_clear = new TH1F("MC_relpTErr", "MC_relpTErr", 50, 10e-3, 0.5);

		Vtx_MC_clear = new TH1F("MC_Vtx", "MC_Vtx",  60, 0, 30);

		ValidMu_MC_clear = new TH1F("MC_ValidMu", "MC_ValidMu",  55, 0, 55);

    	for(int p=0; p<nentries; p++){
//      	 for(int p=0; p<20000; p++){
//      	 for(int p=0; p<1; p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
		treeMC->GetEntry(p);
				
		if(
// 			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){

				int counting = 0;
				bool trovato = false;
				int successivo = 0;

				if(dil_mass > - 99){			

					c_event = event;
					c_mass = dil_mass;
					int p_1 = p + 1;
// 					std::cout<<"A = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
					treeMC->GetEntry(p_1);
					n_event = event;
					while(c_event == n_event){	
						if(c_mass > 70 && c_mass < 110){
							successivo = 1;
							break;
						}
// 						std::cout<<"B = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
						if(
// 							GoodVtx && 
							fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
							lep_pt[0]>53. && lep_pt[1]>53. && 
							lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
							lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
							fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
							(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
							(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
							lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
							lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
							lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
							lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
							cos_angle>-0.9998 && 
							lep_id[0]*lep_id[1]<0 && 
							(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
							vertex_chi2 < 20
						){					
							n_mass = dil_mass;
// 							counting_OK++;
							if(dil_mass > 70 && dil_mass < 110){
								p = p_1;
// 								std::cout<<"C = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
								trovato = true;
								break;
							}
														
// 							std::cout<<"couting = "<<counting<<std::endl;
						}

						counting++;					
						p_1 = p_1 + 1;
						treeMC->GetEntry(p_1);
						n_event = event;
						
						if(p_1 > nentries)
							break;

					}
				}

				weight[j] *= genWeight;
				weight_BB[j] *= genWeight;
				weight_BE[j] *= genWeight;
				
				gM = gen_dil_mass - 400;
				kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
	 		   	
				if(j >= 30 && j <=38){
					weight[j] *= kFactor;
					weight_BB[j] *= kFactor_BB;
					weight_BE[j] *= kFactor_BE;
				}

					dimuon_MC_2016_clear[j]->Fill(dil_mass,  weight[j]); 
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_MC_2016_clear[j]->Fill(dil_mass,  weight_BB[j]);
					else
						dimuon_BE_MC_2016_clear[j]->Fill(dil_mass,  weight_BE[j]);
						
					for(int i = 1; i <  NMBINS+1; i++){
						float contenuto_mc = dimuon_MC_2016_clear[j]->Integral(i, NMBINS);
						dimuon_cumulative_MC_2016_clear[j]->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BB_MC_2016_clear[j]->Integral(i, NMBINS);
						dimuon_BB_cumulative_MC_2016_clear[j]->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BE_MC_2016_clear[j]->Integral(i, NMBINS);
						dimuon_BE_cumulative_MC_2016_clear[j]->SetBinContent(i, contenuto_mc);
					}

					Vtx_MC_clear->Fill(vertex_chi2,  weight[j]); 
   					
					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_MC_2016_clear[j]->Fill(lep_pt[h],  weight[j]); 
						if(fabs(lep_eta[h]) < 1.2)
							pt_BB_MC_2016_clear[j]->Fill(lep_pt[h],  weight_BB[j]);
						else
							pt_BE_MC_2016_clear[j]->Fill(lep_pt[h],  weight_BE[j]);
							
						if(lep_qOverPt[h] > 0)
							pt_MC_plus_clear->Fill(lep_pt[h],  weight[j]);
						else
							pt_MC_minus_clear->Fill(lep_pt[h],  weight[j]);

						if(fabs(lep_eta[h]) < 1.2){							
							float BF = (weight[j]/LUMINOSITY_2016)*19900;
							float GH = (weight[j]/LUMINOSITY_2016)*16300;
							
							pt_MC_RunBF_clear->Fill(lep_pt[h],  BF);
							pt_MC_RunGH_clear->Fill(lep_pt[h],  GH);
						}
													
						for(int i = 1; i <  binnum_pt+1; i++){
							float contenuto_mc = pt_MC_2016_clear[j]->Integral(i, binnum_pt);
							pt_cumulative_MC_2016_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_BB_MC_2016_clear[j]->Integral(i, binnum_pt);
							pt_BB_cumulative_MC_2016_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_BE_MC_2016_clear[j]->Integral(i, binnum_pt);
							pt_BE_cumulative_MC_2016_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_plus_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_plus_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_minus_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_minus_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_RunBF_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_RunBF_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_RunGH_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_RunGH_clear->SetBinContent(i, contenuto_mc);
						}

						eta_MC_2016_clear[j]->Fill(lep_eta[h],  weight[j]); 

						phi_MC_2016_clear[j]->Fill(lep_phi[h],  weight[j]); 
                                                
						dB_MC_clear->Fill(lep_dB[h],  weight[j]); 

						PixelHit_MC_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]); 

						TkLayer_MC_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]); 

						Iso_MC_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]); 

						relpTErr_MC_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]); 

						ValidMu_MC_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]); 
   	 						    	 						
	     			} // for on muons

					if(j >= 30 && j <=38){
						weight[j] /= kFactor;
						weight_BB[j] /= kFactor_BB;
						weight_BE[j] /= kFactor_BE;
					}
				
					weight[j] /= genWeight;
					weight_BB[j] /= genWeight;
					weight_BE[j] /= genWeight;

					if(!trovato) p += counting;
					p += successivo;		

	       }//end of condition on event

		 }// end loop p
	 
	 	 if(j == 44) color = kMagenta-10;
	 	 if(j == 43) color = kMagenta-8;
	 	 if(j == 42) color = kMagenta-6;
	 	 if(j == 41) color = kMagenta-2;
	 	 if(j == 40) color = kMagenta;
	 	 if(j == 39) color = kMagenta+3;
	 	 if(j >= 30 && j < 39){ 
	 	 	color = 3;
	 	 	if(j == 30) legend->AddEntry(dimuon_MC_2016_clear[j], "DY #rightarrow #mu#mu", "f");
	 	 }
// 	 	 if(j == 30) color = kGreen+3;
	 	 if(j > 27 && j < 30){
	 	 	color = kRed+3;
	 	 	if(j == 28) legend->AddEntry(dimuon_MC_2016_clear[j], "WZ", "f");
	 	 }
	 	 if(j > 25 && j < 28){
	 	 	color = 2;
	 	 	if(j == 26) legend->AddEntry(dimuon_MC_2016_clear[j], "ZZ", "f");
	 	 }
	 	 if(j > 20 && j < 26){
	 	 	color = kRed-10;
	 	 	if(j == 23) legend->AddEntry(dimuon_MC_2016_clear[j], "WW", "f");
	 	 }
	 	 if(j > 15 && j < 21){
	 	 	color = 4;
	 	 	if(j == 18) legend->AddEntry(dimuon_MC_2016_clear[j], "t#bar{t}", "f");
	 	 }
	 	 if(j == 14 || j == 15){
	 	 	color = kBlue+2;
	 	 	if(j == 14) legend->AddEntry(dimuon_MC_2016_clear[j], "SingleTop", "f");
	 	 }
	 	 if(j == 13){
	 	 	color = kYellow;
	 	 	legend->AddEntry(dimuon_MC_2016_clear[j], "DY #rightarrow #tau#tau", "f");
	 	 }
		 dimuon_MC_2016_clear[j]->SetFillColor(color);
		 dimuon_MC_2016_clear[j]->SetLineColor(color);
		 dimuon_BB_MC_2016_clear[j]->SetFillColor(color);
		 dimuon_BB_MC_2016_clear[j]->SetLineColor(color);
		 dimuon_BE_MC_2016_clear[j]->SetFillColor(color);
		 dimuon_BE_MC_2016_clear[j]->SetLineColor(color);
		 dimuon_cumulative_MC_2016_clear[j]->SetFillColor(color);
		 dimuon_cumulative_MC_2016_clear[j]->SetLineColor(color);
		 dimuon_BB_cumulative_MC_2016_clear[j]->SetFillColor(color);
		 dimuon_BB_cumulative_MC_2016_clear[j]->SetLineColor(color);
		 dimuon_BE_cumulative_MC_2016_clear[j]->SetFillColor(color);
		 dimuon_BE_cumulative_MC_2016_clear[j]->SetLineColor(color);
		pt_MC_2016_clear[j]->SetFillColor(color);
		pt_MC_2016_clear[j]->SetLineColor(color);
		pt_BB_MC_2016_clear[j]->SetFillColor(color);
		pt_BB_MC_2016_clear[j]->SetLineColor(color);
		pt_BE_MC_2016_clear[j]->SetFillColor(color);
		pt_BE_MC_2016_clear[j]->SetLineColor(color);
		pt_cumulative_MC_2016_clear[j]->SetFillColor(color);
		pt_cumulative_MC_2016_clear[j]->SetLineColor(color);
		pt_BB_cumulative_MC_2016_clear[j]->SetFillColor(color);
		pt_BB_cumulative_MC_2016_clear[j]->SetLineColor(color);
		pt_BE_cumulative_MC_2016_clear[j]->SetFillColor(color);
		pt_BE_cumulative_MC_2016_clear[j]->SetLineColor(color);
		eta_MC_2016_clear[j]->SetFillColor(color);
		eta_MC_2016_clear[j]->SetLineColor(color);
		phi_MC_2016_clear[j]->SetFillColor(color);
		phi_MC_2016_clear[j]->SetLineColor(color);
		pt_MC_plus_clear->SetFillColor(color);
		pt_MC_plus_clear->SetLineColor(color);
		pt_MC_minus_clear->SetFillColor(color);
		pt_MC_minus_clear->SetLineColor(color);
		pt_MC_RunBF_clear->SetFillColor(color);
		pt_MC_RunBF_clear->SetLineColor(color);
		pt_MC_RunGH_clear->SetFillColor(color);
		pt_MC_RunGH_clear->SetLineColor(color);
		pt_cumulative_MC_plus_clear->SetFillColor(color);
		pt_cumulative_MC_plus_clear->SetLineColor(color);
		pt_cumulative_MC_minus_clear->SetFillColor(color);
		pt_cumulative_MC_minus_clear->SetLineColor(color);
		pt_cumulative_MC_RunBF_clear->SetFillColor(color);
		pt_cumulative_MC_RunBF_clear->SetLineColor(color);
		pt_cumulative_MC_RunGH_clear->SetFillColor(color);
		pt_cumulative_MC_RunGH_clear->SetLineColor(color);
		dB_MC_clear->SetFillColor(color);
		dB_MC_clear->SetLineColor(color);
		PixelHit_MC_clear->SetFillColor(color);
		PixelHit_MC_clear->SetLineColor(color);
		TkLayer_MC_clear->SetFillColor(color);
		TkLayer_MC_clear->SetLineColor(color);
		Iso_MC_clear->SetFillColor(color);
		Iso_MC_clear->SetLineColor(color);
		relpTErr_MC_clear->SetFillColor(color);
		relpTErr_MC_clear->SetLineColor(color);
		Vtx_MC_clear->SetFillColor(color);
		Vtx_MC_clear->SetLineColor(color);
		ValidMu_MC_clear->SetFillColor(color);
		ValidMu_MC_clear->SetLineColor(color);
// 		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[j], "HIST");
// 		dimuon_BB_MC_Stack->Add(dimuon_BB_MC_2016_clear[j], "HIST");
// 		dimuon_BE_MC_Stack->Add(dimuon_BE_MC_2016_clear[j], "HIST");
// 		dimuon_cumulative_MC_Stack->Add(dimuon_cumulative_MC_2016_clear[j], "HIST");
// 		dimuon_cumulative_BB_MC_Stack->Add(dimuon_BB_cumulative_MC_2016_clear[j], "HIST");
// 		dimuon_cumulative_BE_MC_Stack->Add(dimuon_BE_cumulative_MC_2016_clear[j], "HIST");
// 		pt_MC_Stack->Add(pt_MC_2016_clear, "HIST");
// 		pt_MC_BB_Stack->Add(pt_BB_MC_2016_clear, "HIST");
// 		pt_MC_BE_Stack->Add(pt_BE_MC_2016_clear, "HIST");
		pt_MC_plus_Stack->Add(pt_MC_plus_clear, "HIST");
		pt_MC_minus_Stack->Add(pt_MC_minus_clear, "HIST");
		pt_MC_RunBF_Stack->Add(pt_MC_RunBF_clear, "HIST");
		pt_MC_RunGH_Stack->Add(pt_MC_RunGH_clear, "HIST");
// 		pt_cumulative_MC_Stack->Add(pt_cumulative_MC_2016_clear[j], "HIST");
// 		pt_cumulative_MC_BB_Stack->Add(pt_BB_cumulative_MC_2016_clear[j], "HIST");
// 		pt_cumulative_MC_BE_Stack->Add(pt_BE_cumulative_MC_2016_clear[j], "HIST");
		pt_cumulative_MC_plus_Stack->Add(pt_cumulative_MC_plus_clear, "HIST");
		pt_cumulative_MC_minus_Stack->Add(pt_cumulative_MC_minus_clear, "HIST");
		pt_cumulative_MC_RunBF_Stack->Add(pt_cumulative_MC_RunBF_clear, "HIST");
		pt_cumulative_MC_RunGH_Stack->Add(pt_cumulative_MC_RunGH_clear, "HIST");
// 		eta_MC_Stack->Add(eta_MC_2017_clear[45], "HIST");
// 		phi_MC_Stack->Add(phi_MC_2016_clear[j], "HIST");
		dB_MC_Stack->Add(dB_MC_clear, "HIST");
		PixelHit_MC_Stack->Add(PixelHit_MC_clear, "HIST");
		TkLayer_MC_Stack->Add(TkLayer_MC_clear, "HIST");
		Iso_MC_Stack->Add(Iso_MC_clear, "HIST");
		relpTErr_MC_Stack->Add(relpTErr_MC_clear, "HIST");
		Vtx_MC_Stack->Add(Vtx_MC_clear, "HIST");
		ValidMu_MC_Stack->Add(ValidMu_MC_clear, "HIST");

	 	 if(j > 30 && j < 39) color = 3;
	 	 if(j == 30) color = kGreen+3;
	 	 if(j > 27 && j < 30) color = kRed+3;
	 	 if(j > 25 && j < 28) color = 2;
	 	 if(j > 20 && j < 26) color = kRed-10;
	 	 if(j > 15 && j < 21) color = 4;
	 	 if(j == 14 || j == 15) color = kBlue+2;		
	 	 
    }// end loop on MC 2016
       
       
    std::cout<<"PASSO AL 2017"<<std::endl;

	for(int j = 0; j < 16; j++){

        std::cout<<"opening.. "<<samples_2017[j]<<" --- "<<events_2017[j]<<std::endl;

        TChain *treeMC = new TChain("SimpleNtupler/t");
//         treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_Pt_Assignment/ana_datamc_" + samples[j] + ".root");

        treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/ana_datamc_" + samples_2017[j] + ".root");        
//      	TFile *file_MC = new TFile("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/ana_datamc_" + samples[j] + ".root", "READ");
//      	file_MC->cd("Our2016MuonsPlusMuonsMinusHistos");
		
// 		TH1F* dimuon_MC_2016_clear[j] = (TH1F*)gDirectory->Get("DimuonMassVtxConstrainedLog");
// 		TH1F* dimuon_BB_MC_2016_clear[j] = (TH1F*)gDirectory->Get("DimuonMassVtxConstrainedLog_bb");
// 		TH1F* dimuon_BE_MC_2016_clear[j] = (TH1F*)gDirectory->Get("DimuonMassVtxConstrainedLog_be");
// 		TH1F* pt_MC_2016_clear = (TH1F*)gDirectory->Get("LeptonPt");
// 		TH1F* eta_MC_2017_clear[45] = (TH1F*)gDirectory->Get("LeptonEta");
// 		TH1F* phi_MC_2016_clear[j] = (TH1F*)gDirectory->Get("LeptonPhi");
// 		TH1F* PixelHit_MC_clear = (TH1F*)gDirectory->Get("NPxHits");
// 		TH1F* TkLayer_MC_clear = (TH1F*)gDirectory->Get("NTkLayers");
// 		TH1F* Iso_MC_clear = (TH1F*)gDirectory->Get("RelIsoSumPt");
// 		TH1F* relpTErr_MC_clear = (TH1F*)gDirectory->Get("DimuonMuonPtErrOverPt");
// 		TH1F* Vtx_MC_clear = (TH1F*)gDirectory->Get("DimuonMassVtx_chi2");
// 		TH1F* ValidMu_MC_clear = (TH1F*)gDirectory->Get("NMuHits");
// 
// 		dimuon_MC_2016_clear[j]->Scale(weight[j]);
// 		dimuon_BB_MC_2016_clear[j]->Scale(weight[j]*Z_peak_BB/Z_peak);
// 		dimuon_BE_MC_2016_clear[j]->Scale(weight[j]*Z_peak_BE/Z_peak);
// 		pt_MC_2016_clear[j]->Scale(weight[j]);
// 		pt_MC_2016_clear[j]->Rebin(50);
// // 		eta_MC_2017_clear[45]->Scale(weight[j]);
// // 		phi_MC_2016_clear[j]->Scale(weight[j]);
// 		PixelHit_MC_clear->Scale(weight[j]);
// 		TkLayer_MC_clear->Scale(weight[j]);
// 		Iso_MC_clear->Scale(weight[j]);
// 		relpTErr_MC_clear->Scale(weight[j]);
// 		relpTErr_MC_clear->Rebin(5);
// 		Vtx_MC_clear->Scale(weight[j]);
// 		Vtx_MC_clear->Rebin(10);
// 		ValidMu_MC_clear->Scale(weight[j]);

	    Long64_t nentries = treeMC->GetEntries();
	    
    	treeMC->SetBranchAddress("genWeight",&genWeight);
     	    	
      	treeMC->SetBranchAddress("event", &event);
      	treeMC->SetBranchAddress("run", &run);
//       	treeMC->SetBranchAddress("lumi", &lumi);    	
    	treeMC->SetBranchAddress("dil_mass",&dil_mass);
    	
    	treeMC->SetBranchAddress("gen_lep_pt", gen_lep_pt);
    	treeMC->SetBranchAddress("gen_lep_eta", gen_lep_eta);
    	treeMC->SetBranchAddress("gen_lep_phi", gen_lep_phi);
    	treeMC->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
	     treeMC->SetBranchAddress("dil_pt",&dil_pt);
	     treeMC->SetBranchAddress("cos_angle",&cos_angle);
 	    treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeMC->SetBranchAddress("nvertices",&nvertices);
	     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);

	     treeMC->SetBranchAddress("lep_pt",lep_pt);
    	 treeMC->SetBranchAddress("lep_id",lep_id);
	     treeMC->SetBranchAddress("lep_eta",lep_eta);
	     treeMC->SetBranchAddress("lep_phi",lep_phi);
	     treeMC->SetBranchAddress("lep_dB",lep_dB);
	     treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	     treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    	 treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	     treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
	     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
// 	     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
	     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    	 treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
	     treeMC->SetBranchAddress("lep_pfIso",lep_pfIso);
	     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
    	 treeMC->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    	 treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
    	 treeMC->SetBranchAddress("lep_qOverPt",lep_qOverPt);


     //treeMC->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeMC->SetBranchAddress("met_pt",&met_pt); 	
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples_2017[j].Data(),j , weight_2017[j], nentries);

         dimuon_MC_2017_clear[j] = new TH1F(samples[j] + "MC_dimuon_2017", samples[j] + "MC_dimuon_2017", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BB_MC_2017_clear[j] = new TH1F(samples[j] + "MC_BB_dimuon_2017", samples[j] + "MC_BB_dimuon_2017", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BE_MC_2017_clear[j] = new TH1F(samples[j] + "MC_BE_dimuon_2017", samples[j] + "MC_BE_dimuon_2017", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_cumulative_MC_2017_clear[j] = new TH1F(samples[j] + "MC_dimuon_cumulative_2017", samples[j] + "MC_dimuon_cumulative_2017", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BB_cumulative_MC_2017_clear[j] = new TH1F(samples[j] + "MC_dimuon_BB_cumulative_2017", samples[j] + "MC_dimuon_BB_cumulative_2017", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BE_cumulative_MC_2017_clear[j] = new TH1F(samples[j] + "MC_dimuon_BE_cumulative_2017", samples[j] + "MC_dimuon_BE_cumulative_2017", 50, ciao_bhu);//NMBINS, logMbins);         
		pt_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_2017", samples[j] + "_Lepton_Pt_2017", binnum_pt, PT_BINS);
		pt_BB_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BB_2017", samples[j] + "_Lepton_Pt_BB_2017", binnum_pt, PT_BINS);
		pt_BE_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BE_2017", samples[j] + "_Lepton_Pt_BE_2017", binnum_pt, PT_BINS);
		pt_cumulative_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_cumulative_2017", samples[j] + "_Lepton_Pt_cumulative_2017", binnum_pt, PT_BINS);
		pt_BB_cumulative_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BB_cumulative_2017", samples[j] + "_Lepton_Pt_BB_cumulative_2017", binnum_pt, PT_BINS);
		pt_BE_cumulative_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BE_cumulative_2017", samples[j] + "_Lepton_Pt_BE_cumulative_2017", binnum_pt, PT_BINS);
		eta_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Eta_2017", samples[j] + "_Lepton_Eta_2017", binnum_eta, ETA_BINS);
		phi_MC_2017_clear[j] = new TH1F(samples[j] + "_Lepton_Phi_2017", samples[j] + "_Lepton_Phi_2017", binnum_phi, PHI_BINS);
		
		pt_MC_plus_clear = new TH1F("MC_pt_plus", "MC_pt_plus", binnum_pt, PT_BINS);
		pt_MC_minus_clear = new TH1F("MC_pt_minus", "MC_pt_minus", binnum_pt, PT_BINS);
		pt_MC_RunBF_clear = new TH1F("MC_pt_RunBFs", "MC_pt_RunBFs", binnum_pt, PT_BINS);
		pt_MC_RunGH_clear = new TH1F("MC_pt_RunGHs", "MC_pt_RunGHs", binnum_pt, PT_BINS);

		pt_cumulative_MC_plus_clear = new TH1F("MC_pt_cumulative_plus", "MC_pt_cumulative_plus", binnum_pt, PT_BINS);
		pt_cumulative_MC_minus_clear = new TH1F("MC_pt_cumulative_minus", "MC_pt_cumulative_minus", binnum_pt, PT_BINS);
		pt_cumulative_MC_RunBF_clear = new TH1F("MC_pt_cumulative_RunBFs", "MC_pt_cumulative_RunBFs", binnum_pt, PT_BINS);
		pt_cumulative_MC_RunGH_clear = new TH1F("MC_pt_cumulative_RunGHs", "MC_pt_cumulative_RunGHs", binnum_pt, PT_BINS);

		dB_MC_clear = new TH1F("MC_dB", "MC_dB",  100, 0.005, 0.5);

		PixelHit_MC_clear = new TH1F("MC_PixelHit", "MC_PixelHit",  15, 0,  15);

		TkLayer_MC_clear = new TH1F("MC_TkLayer", "MC_TkLayer", 20, 0, 20);

		Iso_MC_clear = new TH1F("MC_Iso", "MC_Iso", 60, 0.0, 0.3);

		relpTErr_MC_clear = new TH1F("MC_relpTErr", "MC_relpTErr", 50, 10e-3, 0.5);

		Vtx_MC_clear = new TH1F("MC_Vtx", "MC_Vtx",  60, 0, 30);

		ValidMu_MC_clear = new TH1F("MC_ValidMu", "MC_ValidMu",  55, 0, 55);

    	for(int p=0; p<nentries; p++){
//      	 for(int p=0; p<10000; p++){
//      	 for(int p=0; p<0; p++){

    	 	if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
    	 	
    	 	treeMC->GetEntry(p);
    	 	
// 			if(j > 30) std::cout<<genWeight;
			
		if(
// 			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
				prev_event = event;

				weight_2017[j] *= genWeight;
				weight_2017_BB[j] *= genWeight;
				weight_2017_BE[j] *= genWeight;
				
// 				gM = gen_dil_mass - 400;
// 				kFactor = 1.047 - 0.000143 * gM + 5.167e-08 * pow(gM,2) - 7.84e-12 * pow(gM,3);
// 			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
// 	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
	 		   	
// 				if(j >= 6 && j <=14){
// 					weight[j] *= kFactor;
// 					weight_BB[j] *= kFactor_BB;
// 					weight_BE[j] *= kFactor_BE;
// 				}

					dimuon_MC_2017_clear[j]->Fill(dil_mass,  weight_2017[j]); 
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_MC_2017_clear[j]->Fill(dil_mass,  weight_2017_BB[j]);
					else
						dimuon_BE_MC_2017_clear[j]->Fill(dil_mass,  weight_2017_BE[j]);
						
					for(int i = 1; i <  NMBINS+1; i++){
						float contenuto_mc = dimuon_MC_2017_clear[j]->Integral(i, NMBINS);
						dimuon_cumulative_MC_2017_clear[j]->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BB_MC_2017_clear[j]->Integral(i, NMBINS);
						dimuon_BB_cumulative_MC_2017_clear[j]->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BE_MC_2017_clear[j]->Integral(i, NMBINS);
						dimuon_BE_cumulative_MC_2017_clear[j]->SetBinContent(i, contenuto_mc);
					}

					Vtx_MC_clear->Fill(vertex_chi2,  weight_2017[j]); 
   					
					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_MC_2017_clear[j]->Fill(lep_pt[h],  weight_2017[j]); 
						if(fabs(lep_eta[h]) < 1.2)
							pt_BB_MC_2017_clear[j]->Fill(lep_pt[h],  weight_2017_BB[j]);
						else
							pt_BE_MC_2017_clear[j]->Fill(lep_pt[h],  weight_2017_BE[j]);
							
						if(lep_qOverPt[h] > 0)
							pt_MC_plus_clear->Fill(lep_pt[h],  weight_2017[j]);
						else
							pt_MC_minus_clear->Fill(lep_pt[h],  weight_2017[j]);
							
						for(int i = 1; i <  binnum_pt+1; i++){
							float contenuto_mc = pt_MC_2017_clear[j]->Integral(i, binnum_pt);
							pt_cumulative_MC_2017_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_BB_MC_2017_clear[j]->Integral(i, binnum_pt);
							pt_BB_cumulative_MC_2017_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_BE_MC_2017_clear[j]->Integral(i, binnum_pt);
							pt_BE_cumulative_MC_2017_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_plus_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_plus_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_minus_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_minus_clear->SetBinContent(i, contenuto_mc);
						}

						eta_MC_2017_clear[j]->Fill(lep_eta[h],  weight_2017[j]); 

						phi_MC_2017_clear[j]->Fill(lep_phi[h],  weight_2017[j]); 
                                                
						dB_MC_clear->Fill(lep_dB[h],  weight_2017[j]); 

						PixelHit_MC_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight_2017[j]); 

						TkLayer_MC_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight_2017[j]); 

						Iso_MC_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight_2017[j]); 

						relpTErr_MC_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight_2017[j]); 

						ValidMu_MC_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight_2017[j]); 
   	 						    	 						
	     			} // for on muons

// 				if(j >= 6 && j <=14){
// 					weight[j] /= kFactor;
// 					weight_BB[j] /= kFactor_BB;
// 					weight_BE[j] /= kFactor_BE;
// 				}
			
			weight_2017[j] /= genWeight;
			weight_2017_BB[j] /= genWeight;
			weight_2017_BE[j] /= genWeight;
	       }//end of condition on event

		 }// end loop p
	 
// 	 	 if(j == 44) color = kMagenta-10;
// 	 	 if(j == 43) color = kMagenta-8;
// 	 	 if(j == 42) color = kMagenta-6;
// 	 	 if(j == 41) color = kMagenta-2;
// 	 	 if(j == 40) color = kMagenta;
// 	 	 if(j == 39) color = kMagenta+3;

	 	 if(j >= 7 && j <= 15) color = 3;	 	 	
	 	 if(j == 6) color = kRed+3;
	 	 if(j == 5)	color = 2;
	 	 if(j == 4)	color = kRed-10;
	 	 if(j == 3) color = 4;
	 	 if(j <= 2) color = kBlue+2;

		 dimuon_MC_2017_clear[j]->SetFillColor(color);
		 dimuon_MC_2017_clear[j]->SetLineColor(color);
		 dimuon_BB_MC_2017_clear[j]->SetFillColor(color);
		 dimuon_BB_MC_2017_clear[j]->SetLineColor(color);
		 dimuon_BE_MC_2017_clear[j]->SetFillColor(color);
		 dimuon_BE_MC_2017_clear[j]->SetLineColor(color);
		 dimuon_cumulative_MC_2017_clear[j]->SetFillColor(color);
		 dimuon_cumulative_MC_2017_clear[j]->SetLineColor(color);
		 dimuon_BB_cumulative_MC_2017_clear[j]->SetFillColor(color);
		 dimuon_BB_cumulative_MC_2017_clear[j]->SetLineColor(color);
		 dimuon_BE_cumulative_MC_2017_clear[j]->SetFillColor(color);
		 dimuon_BE_cumulative_MC_2017_clear[j]->SetLineColor(color);
		pt_MC_2017_clear[j]->SetFillColor(color);
		pt_MC_2017_clear[j]->SetLineColor(color);
		pt_BB_MC_2017_clear[j]->SetFillColor(color);
		pt_BB_MC_2017_clear[j]->SetLineColor(color);
		pt_BE_MC_2017_clear[j]->SetFillColor(color);
		pt_BE_MC_2017_clear[j]->SetLineColor(color);
		pt_cumulative_MC_2017_clear[j]->SetFillColor(color);
		pt_cumulative_MC_2017_clear[j]->SetLineColor(color);
		pt_BB_cumulative_MC_2017_clear[j]->SetFillColor(color);
		pt_BB_cumulative_MC_2017_clear[j]->SetLineColor(color);
		pt_BE_cumulative_MC_2017_clear[j]->SetFillColor(color);
		pt_BE_cumulative_MC_2017_clear[j]->SetLineColor(color);
		eta_MC_2017_clear[j]->SetFillColor(color);
		eta_MC_2017_clear[j]->SetLineColor(color);
		phi_MC_2017_clear[j]->SetFillColor(color);
		phi_MC_2017_clear[j]->SetLineColor(color);
		pt_MC_plus_clear->SetFillColor(color);
		pt_MC_plus_clear->SetLineColor(color);
		pt_MC_minus_clear->SetFillColor(color);
		pt_MC_minus_clear->SetLineColor(color);
		pt_MC_RunBF_clear->SetFillColor(color);
		pt_MC_RunBF_clear->SetLineColor(color);
		pt_MC_RunGH_clear->SetFillColor(color);
		pt_MC_RunGH_clear->SetLineColor(color);
		pt_cumulative_MC_plus_clear->SetFillColor(color);
		pt_cumulative_MC_plus_clear->SetLineColor(color);
		pt_cumulative_MC_minus_clear->SetFillColor(color);
		pt_cumulative_MC_minus_clear->SetLineColor(color);
		pt_cumulative_MC_RunBF_clear->SetFillColor(color);
		pt_cumulative_MC_RunBF_clear->SetLineColor(color);
		pt_cumulative_MC_RunGH_clear->SetFillColor(color);
		pt_cumulative_MC_RunGH_clear->SetLineColor(color);
		dB_MC_clear->SetFillColor(color);
		dB_MC_clear->SetLineColor(color);
		PixelHit_MC_clear->SetFillColor(color);
		PixelHit_MC_clear->SetLineColor(color);
		TkLayer_MC_clear->SetFillColor(color);
		TkLayer_MC_clear->SetLineColor(color);
		Iso_MC_clear->SetFillColor(color);
		Iso_MC_clear->SetLineColor(color);
		relpTErr_MC_clear->SetFillColor(color);
		relpTErr_MC_clear->SetLineColor(color);
		Vtx_MC_clear->SetFillColor(color);
		Vtx_MC_clear->SetLineColor(color);
		ValidMu_MC_clear->SetFillColor(color);
		ValidMu_MC_clear->SetLineColor(color);		
// 		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[j], "HIST");
// 		dimuon_BB_MC_Stack->Add(dimuon_BB_MC_2017_clear[j], "HIST");
// 		dimuon_BE_MC_Stack->Add(dimuon_BE_MC_2017_clear[j], "HIST");
// 		dimuon_cumulative_MC_Stack->Add(dimuon_cumulative_MC_2017_clear[j], "HIST");
// 		dimuon_cumulative_BB_MC_Stack->Add(dimuon_BB_cumulative_MC_2017_clear[j], "HIST");
// 		dimuon_cumulative_BE_MC_Stack->Add(dimuon_BE_cumulative_MC_2017_clear[j], "HIST");
// 		pt_MC_Stack->Add(pt_MC_2017_clear, "HIST");
// 		pt_MC_BB_Stack->Add(pt_BB_MC_2017_clear, "HIST");
// 		pt_MC_BE_Stack->Add(pt_BE_MC_2017_clear, "HIST");
		pt_MC_plus_Stack->Add(pt_MC_plus_clear, "HIST");
		pt_MC_minus_Stack->Add(pt_MC_minus_clear, "HIST");
		pt_MC_RunBF_Stack->Add(pt_MC_RunBF_clear, "HIST");
		pt_MC_RunGH_Stack->Add(pt_MC_RunGH_clear, "HIST");
// 		pt_cumulative_MC_Stack->Add(pt_cumulative_MC_2017_clear[j], "HIST");
// 		pt_cumulative_MC_BB_Stack->Add(pt_BB_cumulative_MC_2017_clear[j], "HIST");
// 		pt_cumulative_MC_BE_Stack->Add(pt_BE_cumulative_MC_2017_clear[j], "HIST");
		pt_cumulative_MC_plus_Stack->Add(pt_cumulative_MC_plus_clear, "HIST");
		pt_cumulative_MC_minus_Stack->Add(pt_cumulative_MC_minus_clear, "HIST");
		pt_cumulative_MC_RunBF_Stack->Add(pt_cumulative_MC_RunBF_clear, "HIST");
		pt_cumulative_MC_RunGH_Stack->Add(pt_cumulative_MC_RunGH_clear, "HIST");
// 		eta_MC_Stack->Add(eta_MC_2017_clear[45], "HIST");
// 		phi_MC_Stack->Add(phi_MC_2017_clear[j], "HIST");
		dB_MC_Stack->Add(dB_MC_clear, "HIST");
		PixelHit_MC_Stack->Add(PixelHit_MC_clear, "HIST");
		TkLayer_MC_Stack->Add(TkLayer_MC_clear, "HIST");
		Iso_MC_Stack->Add(Iso_MC_clear, "HIST");
		relpTErr_MC_Stack->Add(relpTErr_MC_clear, "HIST");
		Vtx_MC_Stack->Add(Vtx_MC_clear, "HIST");
		ValidMu_MC_Stack->Add(ValidMu_MC_clear, "HIST");
			 	 
    }
    // end loop on MC 2017

    std::cout<<"PASSO AL 2018"<<std::endl;

	for(int j = 0; j < 9; j++){

        std::cout<<"opening.. "<<samples_2018[j]<<" --- "<<events_2018[j]<<std::endl;

        TChain *treeMC = new TChain("SimpleNtupler/t");
//         treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_Pt_Assignment/ana_datamc_" + samples[j] + ".root");

        treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2018/MC/ana_datamc_" + samples_2018[j] + ".root");        
//      	TFile *file_MC = new TFile("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/ana_datamc_" + samples[j] + ".root", "READ");
//      	file_MC->cd("Our2016MuonsPlusMuonsMinusHistos");
		
// 		TH1F* dimuon_MC_2017_clear[j] = (TH1F*)gDirectory->Get("DimuonMassVtxConstrainedLog");
// 		TH1F* dimuon_BB_MC_2016_clear[j] = (TH1F*)gDirectory->Get("DimuonMassVtxConstrainedLog_bb");
// 		TH1F* dimuon_BE_MC_2016_clear[j] = (TH1F*)gDirectory->Get("DimuonMassVtxConstrainedLog_be");
// 		TH1F* pt_MC_2016_clear = (TH1F*)gDirectory->Get("LeptonPt");
// 		TH1F* eta_MC_2017_clear[45] = (TH1F*)gDirectory->Get("LeptonEta");
// 		TH1F* phi_MC_2016_clear[j] = (TH1F*)gDirectory->Get("LeptonPhi");
// 		TH1F* PixelHit_MC_clear = (TH1F*)gDirectory->Get("NPxHits");
// 		TH1F* TkLayer_MC_clear = (TH1F*)gDirectory->Get("NTkLayers");
// 		TH1F* Iso_MC_clear = (TH1F*)gDirectory->Get("RelIsoSumPt");
// 		TH1F* relpTErr_MC_clear = (TH1F*)gDirectory->Get("DimuonMuonPtErrOverPt");
// 		TH1F* Vtx_MC_clear = (TH1F*)gDirectory->Get("DimuonMassVtx_chi2");
// 		TH1F* ValidMu_MC_clear = (TH1F*)gDirectory->Get("NMuHits");
// 
// 		dimuon_MC_2017_clear[j]->Scale(weight[j]);
// 		dimuon_BB_MC_2016_clear[j]->Scale(weight[j]*Z_peak_BB/Z_peak);
// 		dimuon_BE_MC_2016_clear[j]->Scale(weight[j]*Z_peak_BE/Z_peak);
// 		pt_MC_2016_clear->Scale(weight[j]);
// 		pt_MC_2016_clear->Rebin(50);
// // 		eta_MC_2017_clear[45]->Scale(weight[j]);
// // 		phi_MC_2016_clear[j]->Scale(weight[j]);
// 		PixelHit_MC_clear->Scale(weight[j]);
// 		TkLayer_MC_clear->Scale(weight[j]);
// 		Iso_MC_clear->Scale(weight[j]);
// 		relpTErr_MC_clear->Scale(weight[j]);
// 		relpTErr_MC_clear->Rebin(5);
// 		Vtx_MC_clear->Scale(weight[j]);
// 		Vtx_MC_clear->Rebin(10);
// 		ValidMu_MC_clear->Scale(weight[j]);

	    Long64_t nentries = treeMC->GetEntries();
	    
    	treeMC->SetBranchAddress("genWeight",&genWeight);
     	    	
      	treeMC->SetBranchAddress("event", &event);
      	treeMC->SetBranchAddress("run", &run);
//       	treeMC->SetBranchAddress("lumi", &lumi);    	
    	treeMC->SetBranchAddress("dil_mass",&dil_mass);
    	
    	treeMC->SetBranchAddress("gen_lep_pt", gen_lep_pt);
    	treeMC->SetBranchAddress("gen_lep_eta", gen_lep_eta);
    	treeMC->SetBranchAddress("gen_lep_phi", gen_lep_phi);
    	treeMC->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
	     treeMC->SetBranchAddress("dil_pt",&dil_pt);
	     treeMC->SetBranchAddress("cos_angle",&cos_angle);
 	    treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeMC->SetBranchAddress("nvertices",&nvertices);
	     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);

	     treeMC->SetBranchAddress("lep_pt",lep_pt);
    	 treeMC->SetBranchAddress("lep_id",lep_id);
	     treeMC->SetBranchAddress("lep_eta",lep_eta);
	     treeMC->SetBranchAddress("lep_phi",lep_phi);
	     treeMC->SetBranchAddress("lep_dB",lep_dB);
	     treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	     treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    	 treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	     treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
	     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
// 	     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
	     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    	 treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
	     treeMC->SetBranchAddress("lep_pfIso",lep_pfIso);
	     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
    	 treeMC->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    	 treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
    	 treeMC->SetBranchAddress("lep_qOverPt",lep_qOverPt);


     //treeMC->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeMC->SetBranchAddress("met_pt",&met_pt); 	
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples_2018[j].Data(),j , weight_2018[j], nentries);

         dimuon_MC_2018_clear[j] = new TH1F(samples[j] + "MC_dimuon_2018", samples[j] + "MC_dimuon_2018", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BB_MC_2018_clear[j] = new TH1F(samples[j] + "MC_BB_dimuon_2018", samples[j] + "MC_BB_dimuon_2018", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BE_MC_2018_clear[j] = new TH1F(samples[j] + "MC_BE_dimuon_2018", samples[j] + "MC_BE_dimuon_2018", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_cumulative_MC_2018_clear[j] = new TH1F(samples[j] + "MC_dimuon_cumulative_2018", samples[j] + "MC_dimuon_cumulative_2018", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BB_cumulative_MC_2018_clear[j] = new TH1F(samples[j] + "MC_dimuon_BB_cumulative_2018", samples[j] + "MC_dimuon_BB_cumulative_2018", 50, ciao_bhu);//NMBINS, logMbins);         
         dimuon_BE_cumulative_MC_2018_clear[j] = new TH1F(samples[j] + "MC_dimuon_BE_cumulative_2018", samples[j] + "MC_dimuon_BE_cumulative_2018", 50, ciao_bhu);//NMBINS, logMbins);         
		pt_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_2018", samples[j] + "_Lepton_Pt_2018", binnum_pt, PT_BINS);
		pt_BB_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BB_2018", samples[j] + "_Lepton_Pt_BB_2018", binnum_pt, PT_BINS);
		pt_BE_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BE_2018", samples[j] + "_Lepton_Pt_BE_2018", binnum_pt, PT_BINS);
		pt_cumulative_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_cumulative_2018", samples[j] + "_Lepton_Pt_cumulative_2018", binnum_pt, PT_BINS);
		pt_BB_cumulative_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BB_cumulative_2018", samples[j] + "_Lepton_Pt_BB_cumulative_2018", binnum_pt, PT_BINS);
		pt_BE_cumulative_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Pt_BE_cumulative_2018", samples[j] + "_Lepton_Pt_BE_cumulative_2018", binnum_pt, PT_BINS);
		eta_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Eta_2018", samples[j] + "_Lepton_Eta_2018", binnum_eta, ETA_BINS);
		phi_MC_2018_clear[j] = new TH1F(samples[j] + "_Lepton_Phi_2018", samples[j] + "_Lepton_Phi_2018", binnum_phi, PHI_BINS);
		
		pt_MC_plus_clear = new TH1F("MC_pt_plus", "MC_pt_plus", binnum_pt, PT_BINS);
		pt_MC_minus_clear = new TH1F("MC_pt_minus", "MC_pt_minus", binnum_pt, PT_BINS);
		pt_MC_RunBF_clear = new TH1F("MC_pt_RunBFs", "MC_pt_RunBFs", binnum_pt, PT_BINS);
		pt_MC_RunGH_clear = new TH1F("MC_pt_RunGHs", "MC_pt_RunGHs", binnum_pt, PT_BINS);

		pt_cumulative_MC_plus_clear = new TH1F("MC_pt_cumulative_plus", "MC_pt_cumulative_plus", binnum_pt, PT_BINS);
		pt_cumulative_MC_minus_clear = new TH1F("MC_pt_cumulative_minus", "MC_pt_cumulative_minus", binnum_pt, PT_BINS);
		pt_cumulative_MC_RunBF_clear = new TH1F("MC_pt_cumulative_RunBFs", "MC_pt_cumulative_RunBFs", binnum_pt, PT_BINS);
		pt_cumulative_MC_RunGH_clear = new TH1F("MC_pt_cumulative_RunGHs", "MC_pt_cumulative_RunGHs", binnum_pt, PT_BINS);

		dB_MC_clear = new TH1F("MC_dB", "MC_dB",  100, 0.005, 0.5);

		PixelHit_MC_clear = new TH1F("MC_PixelHit", "MC_PixelHit",  15, 0,  15);

		TkLayer_MC_clear = new TH1F("MC_TkLayer", "MC_TkLayer", 20, 0, 20);

		Iso_MC_clear = new TH1F("MC_Iso", "MC_Iso", 60, 0.0, 0.3);

		relpTErr_MC_clear = new TH1F("MC_relpTErr", "MC_relpTErr", 50, 10e-3, 0.5);

		Vtx_MC_clear = new TH1F("MC_Vtx", "MC_Vtx",  60, 0, 30);

		ValidMu_MC_clear = new TH1F("MC_ValidMu", "MC_ValidMu",  55, 0, 55);

    	for(int p=0; p<nentries; p++){
//      	 for(int p=0; p<10000; p++){
//      	 for(int p=0; p<20; p++){

    	 	if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
    	 	
    	 	treeMC->GetEntry(p);
    	 	
// 			if(j > 30) std::cout<<genWeight;
			
		if(
// 			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
				prev_event = event;

				weight_2018[j] *= genWeight;
				weight_2018_BB[j] *= genWeight;
				weight_2018_BE[j] *= genWeight;
				
// 				gM = gen_dil_mass - 400;
// 				kFactor = 1.047 - 0.000143 * gM + 5.167e-08 * pow(gM,2) - 7.84e-12 * pow(gM,3);
// 			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
// 	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
	 		   	
// 				if(j >= 6 && j <=14){
// 					weight[j] *= kFactor;
// 					weight_BB[j] *= kFactor_BB;
// 					weight_BE[j] *= kFactor_BE;
// 				}

					dimuon_MC_2018_clear[j]->Fill(dil_mass,  weight_2018[j]); 
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_MC_2018_clear[j]->Fill(dil_mass,  weight_2018_BB[j]);
					else
						dimuon_BE_MC_2018_clear[j]->Fill(dil_mass,  weight_2018_BE[j]);
						
					for(int i = 1; i <  NMBINS+1; i++){
						float contenuto_mc = dimuon_MC_2018_clear[j]->Integral(i, NMBINS);
						dimuon_cumulative_MC_2018_clear[j]->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BB_MC_2018_clear[j]->Integral(i, NMBINS);
						dimuon_BB_cumulative_MC_2018_clear[j]->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BE_MC_2018_clear[j]->Integral(i, NMBINS);
						dimuon_BE_cumulative_MC_2018_clear[j]->SetBinContent(i, contenuto_mc);
					}

					Vtx_MC_clear->Fill(vertex_chi2,  weight_2018[j]); 
   					
					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_MC_2018_clear[j]->Fill(lep_pt[h],  weight_2018[j]); 
						if(fabs(lep_eta[h]) < 1.2)
							pt_BB_MC_2018_clear[j]->Fill(lep_pt[h],  weight_2018_BB[j]);
						else
							pt_BE_MC_2018_clear[j]->Fill(lep_pt[h],  weight_2018_BE[j]);
							
						if(lep_qOverPt[h] > 0)
							pt_MC_plus_clear->Fill(lep_pt[h],  weight_2018[j]);
						else
							pt_MC_minus_clear->Fill(lep_pt[h],  weight_2018[j]);
							
						for(int i = 1; i <  binnum_pt+1; i++){
							float contenuto_mc = pt_MC_2018_clear[j]->Integral(i, binnum_pt);
							pt_cumulative_MC_2018_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_BB_MC_2018_clear[j]->Integral(i, binnum_pt);
							pt_BB_cumulative_MC_2018_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_BE_MC_2018_clear[j]->Integral(i, binnum_pt);
							pt_BE_cumulative_MC_2018_clear[j]->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_plus_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_plus_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_minus_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_minus_clear->SetBinContent(i, contenuto_mc);
						}

						eta_MC_2018_clear[j]->Fill(lep_eta[h],  weight_2018[j]); 

						phi_MC_2018_clear[j]->Fill(lep_phi[h],  weight_2018[j]); 
                                                
						dB_MC_clear->Fill(lep_dB[h],  weight_2018[j]); 

						PixelHit_MC_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight_2018[j]); 

						TkLayer_MC_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight_2018[j]); 

						Iso_MC_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight_2018[j]); 

						relpTErr_MC_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight_2018[j]); 

						ValidMu_MC_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight_2018[j]); 
   	 						    	 						
	     			} // for on muons

// 				if(j >= 6 && j <=14){
// 					weight[j] /= kFactor;
// 					weight_BB[j] /= kFactor_BB;
// 					weight_BE[j] /= kFactor_BE;
// 				}
			
			weight_2018[j] /= genWeight;
			weight_2018_BB[j] /= genWeight;
			weight_2018_BE[j] /= genWeight;
	       }//end of condition on event

		 }// end loop p
	 
// 	 	 if(j == 44) color = kMagenta-10;
// 	 	 if(j == 43) color = kMagenta-8;
// 	 	 if(j == 42) color = kMagenta-6;
// 	 	 if(j == 41) color = kMagenta-2;
// 	 	 if(j == 40) color = kMagenta;
// 	 	 if(j == 39) color = kMagenta+3;

	 	 color = 3;	 	 	
// 	 	 if(j >= 6 && j < 14) 
// 	 	 if(j == 5) color = kRed+3;
// 	 	 if(j == 4)	color = 2;
// 	 	 if(j == 3)	color = kRed-10;
// 	 	 if(j == 2) color = 4;
// 	 	 if(j < 2) color = kBlue+2;

		 dimuon_MC_2018_clear[j]->SetFillColor(color);
		 dimuon_MC_2018_clear[j]->SetLineColor(color);
		 dimuon_BB_MC_2018_clear[j]->SetFillColor(color);
		 dimuon_BB_MC_2018_clear[j]->SetLineColor(color);
		 dimuon_BE_MC_2018_clear[j]->SetFillColor(color);
		 dimuon_BE_MC_2018_clear[j]->SetLineColor(color);
		 dimuon_cumulative_MC_2018_clear[j]->SetFillColor(color);
		 dimuon_cumulative_MC_2018_clear[j]->SetLineColor(color);
		 dimuon_BB_cumulative_MC_2018_clear[j]->SetFillColor(color);
		 dimuon_BB_cumulative_MC_2018_clear[j]->SetLineColor(color);
		 dimuon_BE_cumulative_MC_2018_clear[j]->SetFillColor(color);
		 dimuon_BE_cumulative_MC_2018_clear[j]->SetLineColor(color);
		pt_MC_2018_clear[j]->SetFillColor(color);
		pt_MC_2018_clear[j]->SetLineColor(color);
		pt_BB_MC_2018_clear[j]->SetFillColor(color);
		pt_BB_MC_2018_clear[j]->SetLineColor(color);
		pt_BE_MC_2018_clear[j]->SetFillColor(color);
		pt_BE_MC_2018_clear[j]->SetLineColor(color);
		pt_cumulative_MC_2018_clear[j]->SetFillColor(color);
		pt_cumulative_MC_2018_clear[j]->SetLineColor(color);
		pt_BB_cumulative_MC_2018_clear[j]->SetFillColor(color);
		pt_BB_cumulative_MC_2018_clear[j]->SetLineColor(color);
		pt_BE_cumulative_MC_2018_clear[j]->SetFillColor(color);
		pt_BE_cumulative_MC_2018_clear[j]->SetLineColor(color);
		eta_MC_2018_clear[j]->SetFillColor(color);
		eta_MC_2018_clear[j]->SetLineColor(color);
		phi_MC_2018_clear[j]->SetFillColor(color);
		phi_MC_2018_clear[j]->SetLineColor(color);
		pt_MC_plus_clear->SetFillColor(color);
		pt_MC_plus_clear->SetLineColor(color);
		pt_MC_minus_clear->SetFillColor(color);
		pt_MC_minus_clear->SetLineColor(color);
		pt_MC_RunBF_clear->SetFillColor(color);
		pt_MC_RunBF_clear->SetLineColor(color);
		pt_MC_RunGH_clear->SetFillColor(color);
		pt_MC_RunGH_clear->SetLineColor(color);
		pt_cumulative_MC_plus_clear->SetFillColor(color);
		pt_cumulative_MC_plus_clear->SetLineColor(color);
		pt_cumulative_MC_minus_clear->SetFillColor(color);
		pt_cumulative_MC_minus_clear->SetLineColor(color);
		pt_cumulative_MC_RunBF_clear->SetFillColor(color);
		pt_cumulative_MC_RunBF_clear->SetLineColor(color);
		pt_cumulative_MC_RunGH_clear->SetFillColor(color);
		pt_cumulative_MC_RunGH_clear->SetLineColor(color);
		dB_MC_clear->SetFillColor(color);
		dB_MC_clear->SetLineColor(color);
		PixelHit_MC_clear->SetFillColor(color);
		PixelHit_MC_clear->SetLineColor(color);
		TkLayer_MC_clear->SetFillColor(color);
		TkLayer_MC_clear->SetLineColor(color);
		Iso_MC_clear->SetFillColor(color);
		Iso_MC_clear->SetLineColor(color);
		relpTErr_MC_clear->SetFillColor(color);
		relpTErr_MC_clear->SetLineColor(color);
		Vtx_MC_clear->SetFillColor(color);
		Vtx_MC_clear->SetLineColor(color);
		ValidMu_MC_clear->SetFillColor(color);
		ValidMu_MC_clear->SetLineColor(color);		
// 		dimuon_MC_Stack->Add(dimuon_MC_2018_clear[j], "HIST");
// 		dimuon_BB_MC_Stack->Add(dimuon_BB_MC_2018_clear[j], "HIST");
// 		dimuon_BE_MC_Stack->Add(dimuon_BE_MC_2018_clear[j], "HIST");
// 		dimuon_cumulative_MC_Stack->Add(dimuon_cumulative_MC_2018_clear[j], "HIST");
// 		dimuon_cumulative_BB_MC_Stack->Add(dimuon_BB_cumulative_MC_2018_clear[j], "HIST");
// 		dimuon_cumulative_BE_MC_Stack->Add(dimuon_BE_cumulative_MC_2018_clear[j], "HIST");
// 		pt_MC_Stack->Add(pt_MC_2018_clear, "HIST");
// 		pt_MC_BB_Stack->Add(pt_BB_MC_2018_clear, "HIST");
// 		pt_MC_BE_Stack->Add(pt_BE_MC_2018_clear, "HIST");
		pt_MC_plus_Stack->Add(pt_MC_plus_clear, "HIST");
		pt_MC_minus_Stack->Add(pt_MC_minus_clear, "HIST");
		pt_MC_RunBF_Stack->Add(pt_MC_RunBF_clear, "HIST");
		pt_MC_RunGH_Stack->Add(pt_MC_RunGH_clear, "HIST");
// 		pt_cumulative_MC_Stack->Add(pt_cumulative_MC_2018_clear[j], "HIST");
// 		pt_cumulative_MC_BB_Stack->Add(pt_BB_cumulative_MC_2018_clear[j], "HIST");
// 		pt_cumulative_MC_BE_Stack->Add(pt_BE_cumulative_MC_2018_clear[j], "HIST");
		pt_cumulative_MC_plus_Stack->Add(pt_cumulative_MC_plus_clear, "HIST");
		pt_cumulative_MC_minus_Stack->Add(pt_cumulative_MC_minus_clear, "HIST");
		pt_cumulative_MC_RunBF_Stack->Add(pt_cumulative_MC_RunBF_clear, "HIST");
		pt_cumulative_MC_RunGH_Stack->Add(pt_cumulative_MC_RunGH_clear, "HIST");
// 		eta_MC_Stack->Add(eta_MC_2017_clear[45], "HIST");
// 		phi_MC_Stack->Add(phi_MC_2018_clear[j], "HIST");
		dB_MC_Stack->Add(dB_MC_clear, "HIST");
		PixelHit_MC_Stack->Add(PixelHit_MC_clear, "HIST");
		TkLayer_MC_Stack->Add(TkLayer_MC_clear, "HIST");
		Iso_MC_Stack->Add(Iso_MC_clear, "HIST");
		relpTErr_MC_Stack->Add(relpTErr_MC_clear, "HIST");
		Vtx_MC_Stack->Add(Vtx_MC_clear, "HIST");
		ValidMu_MC_Stack->Add(ValidMu_MC_clear, "HIST");
			 	 
    }
    // end loop on MC 2018
    

DrawStackPlot( dimuon_MC_Stack, dimuon_MC_2016_clear, dimuon_MC_2017_clear, dimuon_MC_2018_clear);
DrawStackPlot( dimuon_BB_MC_Stack, dimuon_BB_MC_2016_clear, dimuon_BB_MC_2017_clear, dimuon_BB_MC_2018_clear);
DrawStackPlot( dimuon_BE_MC_Stack, dimuon_BE_MC_2016_clear, dimuon_BE_MC_2017_clear, dimuon_BE_MC_2018_clear);
DrawStackPlot( dimuon_cumulative_MC_Stack, dimuon_cumulative_MC_2016_clear, dimuon_cumulative_MC_2017_clear, dimuon_cumulative_MC_2018_clear);
DrawStackPlot( dimuon_cumulative_BB_MC_Stack, dimuon_BB_cumulative_MC_2016_clear, dimuon_BB_cumulative_MC_2017_clear, dimuon_BB_cumulative_MC_2018_clear);
DrawStackPlot( dimuon_cumulative_BE_MC_Stack, dimuon_BE_cumulative_MC_2016_clear, dimuon_BE_cumulative_MC_2017_clear, dimuon_BE_cumulative_MC_2018_clear);
DrawStackPlot( pt_MC_Stack, pt_MC_2016_clear, pt_MC_2017_clear, pt_MC_2018_clear);
DrawStackPlot( pt_BB_MC_Stack, pt_BB_MC_2016_clear, pt_BB_MC_2017_clear, pt_BB_MC_2018_clear);
DrawStackPlot( pt_BE_MC_Stack, pt_BE_MC_2016_clear, pt_BE_MC_2017_clear, pt_BE_MC_2018_clear);
DrawStackPlot( pt_cumulative_MC_Stack, pt_cumulative_MC_2016_clear, pt_cumulative_MC_2017_clear, pt_cumulative_MC_2018_clear);
DrawStackPlot( pt_cumulative_MC_BB_Stack, pt_BB_cumulative_MC_2016_clear, pt_BB_cumulative_MC_2017_clear, pt_BB_cumulative_MC_2018_clear);
DrawStackPlot( pt_cumulative_MC_BE_Stack, pt_BE_cumulative_MC_2016_clear, pt_BE_cumulative_MC_2017_clear, pt_BE_cumulative_MC_2018_clear);
DrawStackPlot( eta_MC_Stack, eta_MC_2016_clear, eta_MC_2017_clear, eta_MC_2018_clear);
DrawStackPlot( phi_MC_Stack, phi_MC_2016_clear, phi_MC_2017_clear, phi_MC_2018_clear);


////////      DATA    ///////

	for(int k = 0; k < 3; k++){
		
    TChain *treeDATA = new TChain("SimpleNtupler/t");
    if(k == 0) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/DATA_Pt_Assignment/ana_datamc_data_Pt_Assignment.root");
//     if(k == 1) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA_WITH_OR/ana_datamc_data.root");
    if(k == 1) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA_reReco_31March/ana_datamc_data.root");
    if(k == 2) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2018/DATA_17Sept/ana_datamc_data.root");
     	    	
      	treeDATA->SetBranchAddress("event", &event);
      	treeDATA->SetBranchAddress("run", &run);
      	treeDATA->SetBranchAddress("lumi", &lumin);    	
    	treeDATA->SetBranchAddress("dil_mass",&dil_mass);
    	
	     treeDATA->SetBranchAddress("dil_pt",&dil_pt);
	     treeDATA->SetBranchAddress("cos_angle",&cos_angle);
 	    treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeDATA->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeDATA->SetBranchAddress("nvertices",&nvertices);
	     treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);

	     treeDATA->SetBranchAddress("lep_pt",lep_pt);
    	 treeDATA->SetBranchAddress("lep_id",lep_id);
	     treeDATA->SetBranchAddress("lep_eta",lep_eta);
	     treeDATA->SetBranchAddress("lep_phi",lep_phi);
	     treeDATA->SetBranchAddress("lep_dB",lep_dB);
	     treeDATA->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	     treeDATA->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    	 treeDATA->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	     treeDATA->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	     treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
	     treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	     treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
// 	     treeDATA->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
	     treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    	 treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	     treeDATA->SetBranchAddress("lep_sumPt",lep_sumPt);
	     treeDATA->SetBranchAddress("lep_pfIso",lep_pfIso);
	     treeDATA->SetBranchAddress("lep_pt_err",lep_pt_err);
    	 treeDATA->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    	 treeDATA->SetBranchAddress("lep_tk_pt",lep_tk_pt);
    	 treeDATA->SetBranchAddress("lep_qOverPt",lep_qOverPt);



     //treeDATA->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeDATA->SetBranchAddress("met_pt",&met_pt);     

	Long64_t nentries = treeDATA->GetEntries();
	
	printf("opening... DATA --- %lld\n", nentries);    
		int bhu_data = 0;
	
	for(int p=0; p<nentries; p++){
// 	for(int p=0; p<10000; p++){
// 	for(int p=0; p<1; p++){

		if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
	
		treeDATA->GetEntry(p);
		
		if(
// 			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
				prev_event = event;
				
// 				if(run == 317392 && lumin == 666 && event == 910848904)
// 					std::cout<<lep_pt[0]<<"\t"<<lep_pt[1]<<"\t"<< lep_eta[0]<<"\t"<<lep_eta[1]<<std::endl;
// 				if(run == 317182 && lumin == 238 && event == 301576586)
// 					std::cout<<lep_pt[0]<<"\t"<<lep_pt[1]<<"\t"<< lep_eta[0]<<"\t"<<lep_eta[1]<<std::endl;
// 				if(run == 315690 && lumin == 199 && event == 130272508)
// 					std::cout<<lep_pt[0]<<"\t"<<lep_pt[1]<<"\t"<< lep_eta[0]<<"\t"<<lep_eta[1]<<std::endl;

					dimuon_DATA->Fill(dil_mass);
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_DATA->Fill(dil_mass);
					else
						dimuon_BE_DATA->Fill(dil_mass);
					for(int i = 1; i <  NMBINS+1; i++){
						float contenuto_mc = dimuon_DATA->Integral(i, NMBINS);
						dimuon_cumulative_DATA->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BB_DATA->Integral(i, NMBINS);
						dimuon_cumulative_BB_DATA->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BE_DATA->Integral(i, NMBINS);
						dimuon_cumulative_BE_DATA->SetBinContent(i, contenuto_mc);
					}
					Vtx_DATA->Fill(vertex_chi2);

					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_DATA->Fill(lep_pt[h]);
						if(fabs(lep_eta[h]) < 1.2)
							pt_DATA_BB->Fill(lep_pt[h]);
						else
							pt_DATA_BE->Fill(lep_pt[h]);
							
						if(lep_qOverPt[h] > 0)
							pt_DATA_plus->Fill(lep_pt[h]);
						else
							pt_DATA_minus->Fill(lep_pt[h]);

						if(fabs(lep_eta[h]) < 1.2){
							if(run < 278810)
								pt_DATA_RunBF->Fill(lep_pt[h]);
							else
								pt_DATA_RunGH->Fill(lep_pt[h]);
						}
						for(int i = 1; i <  binnum_pt+1; i++){
							float contenuto_mc = pt_DATA->Integral(i, binnum_pt);
							pt_cumulative_DATA->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_BB->Integral(i, binnum_pt);
							pt_cumulative_DATA_BB->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_BE->Integral(i, binnum_pt);
							pt_cumulative_DATA_BE->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_plus->Integral(i, binnum_pt);
							pt_cumulative_DATA_plus->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_minus->Integral(i, binnum_pt);
							pt_cumulative_DATA_minus->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_RunBF->Integral(i, binnum_pt);
							pt_cumulative_DATA_RunBF->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_RunGH->Integral(i, binnum_pt);
							pt_cumulative_DATA_RunGH->SetBinContent(i, contenuto_mc);
						}

						eta_DATA->Fill(lep_eta[h]);

						phi_DATA->Fill(lep_phi[h]);

						dB_DATA->Fill(lep_dB[h]);

						PixelHit_DATA->Fill(lep_glb_numberOfValidPixelHits[h]);

						TkLayer_DATA->Fill(lep_glb_numberOfValidTrackerLayers[h]);

						Iso_DATA->Fill(lep_sumPt[h]/lep_tk_pt[h]);

						relpTErr_DATA->Fill(lep_pt_err[h]/lep_pt[h]);

						ValidMu_DATA->Fill(lep_glb_numberOfValidMuonHits[h]);
   	 						    	 						
	     			} // for on muons

	       		}//end of condition on event

		 }// end loop p
	} //Run16 and 17 and 18
	
// 	return;
			 
	gROOT->Reset();
	gROOT->SetBatch();
	
// 	TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
// 	dimuon_DATA->Draw();
// 	
// 	TCanvas *c21 = new TCanvas("c21", "21", 500, 500);
// 	dimuon_cumulative_DATA->Draw();
// 	return;
	
	TString dir_save = "./";
// 	TString dir_save = "./20162017/";

	save_document = dir_save + "Variuos_distribution.pdf";

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf[", 0,  0, "p_{T}", 0);

	SalvaHisto("dimuon Mass: stack plot", dimuon_MC_Stack, dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon Mass BB: stack plot", dimuon_BB_MC_Stack, dimuon_BB_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon Mass BE: stack plot", dimuon_BE_MC_Stack, dimuon_BE_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);

	SalvaHisto("dimuon cumulative Mass: stack plot", dimuon_cumulative_MC_Stack, dimuon_cumulative_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon cumulative Mass BB: stack plot", dimuon_cumulative_BB_MC_Stack, dimuon_cumulative_BB_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon cumulative Mass BE: stack plot", dimuon_cumulative_BE_MC_Stack, dimuon_cumulative_BE_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);

	SalvaHisto("Lepton p_{T}", pt_MC_Stack, pt_DATA, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BB", pt_BB_MC_Stack, pt_DATA_BB, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BE", pt_BE_MC_Stack, pt_DATA_BE, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} plus", pt_MC_plus_Stack, pt_DATA_plus, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} minus", pt_MC_minus_Stack, pt_DATA_minus, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} RunBF (BB)", pt_MC_RunBF_Stack, pt_DATA_RunBF, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} RunGH (BB)", pt_MC_RunGH_Stack, pt_DATA_RunGH, save_document, 0,  1, "p_{T}", 10e-5);

	SalvaHisto("Lepton p_{T} cumulative", pt_cumulative_MC_Stack, pt_cumulative_DATA, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BB cumulative", pt_cumulative_MC_BB_Stack, pt_cumulative_DATA_BB, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BE cumulative", pt_cumulative_MC_BE_Stack, pt_cumulative_DATA_BE, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} plus cumulative", pt_cumulative_MC_plus_Stack, pt_cumulative_DATA_plus, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} minus cumulative", pt_cumulative_MC_minus_Stack, pt_cumulative_DATA_minus, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} RunBF cumulative", pt_cumulative_MC_RunBF_Stack, pt_cumulative_DATA_RunBF, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} RunGH cumulative", pt_cumulative_MC_RunGH_Stack, pt_cumulative_DATA_RunGH, save_document, 0,  1, "p_{T}", 10e-5);

	SalvaHisto("Lepton #eta", eta_MC_Stack, eta_DATA, save_document, 0,  1, "#eta", 10);

	SalvaHisto("Lepton #phi", phi_MC_Stack, phi_DATA, save_document, 0,  1, "#phi", 10);

	SalvaHisto("dB", dB_MC_Stack, dB_DATA, save_document, 1,  1, "dB", 10e-5);

	SalvaHisto("Pixel Hits", PixelHit_MC_Stack, PixelHit_DATA, save_document, 0,  1, "PixelHits", 10e-1);

	SalvaHisto("TkLayer", TkLayer_MC_Stack, TkLayer_DATA, save_document, 0,  1, "TkLayer", 10e-1);

	SalvaHisto("Rel Tk Iso", Iso_MC_Stack, Iso_DATA, save_document, 0,  1, "Rel Tk Iso", 10e-1);

	SalvaHisto("relPtErro", relpTErr_MC_Stack, relpTErr_DATA, save_document, 1,  1, "#sigma_{p_T}/p_T", 10e-1);

	SalvaHisto("Vtx", Vtx_MC_Stack, Vtx_DATA, save_document, 0,  1, "Vtx #chi^2", 10e-1);

	SalvaHisto("ValidMuHits", ValidMu_MC_Stack, ValidMu_DATA, save_document, 0,  1, "#mu hits", 10e-1);

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf]", 0,  0, "p_{T}", 0);
	
}// end of function 















void DrawStackPlot(THStack* dimuon_MC_Stack, TH1F* dimuon_MC_2016_clear[45], TH1F* dimuon_MC_2017_clear[16], TH1F* dimuon_MC_2018_clear[9]){

		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[1], "HIST");	
		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[2], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[14], "HIST");	
		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[15], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[3], "HIST");
		for(int i = 16; i <= 20; i++)
			dimuon_MC_Stack->Add(dimuon_MC_2016_clear[i], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[4], "HIST");
		for(int i = 21; i <= 25; i++)
			dimuon_MC_Stack->Add(dimuon_MC_2016_clear[i], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[5], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[26], "HIST");	
		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[27], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2017_clear[6], "HIST");
		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[28], "HIST");	
		dimuon_MC_Stack->Add(dimuon_MC_2016_clear[29], "HIST");
		for(int i = 30; i <= 38; i++)
			dimuon_MC_Stack->Add(dimuon_MC_2016_clear[i], "HIST");
		for(int i = 7; i <= 15; i++)
			dimuon_MC_Stack->Add(dimuon_MC_2017_clear[i], "HIST");
		for(int i = 0; i < 9; i++)
			dimuon_MC_Stack->Add(dimuon_MC_2018_clear[i], "HIST");

}






void SalvaHisto(TString name, TH2F* h_MC, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic, bool logX){

	TCanvas *c1 = new TCanvas(name, name,  500, 500);//210,45,1050,750);
   	c1->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0, 1, 1);
//    	pad11->SetGrid();
   	pad11->SetBottomMargin(0.10);
   	pad11->SetLeftMargin(0.10);
   	pad11->SetRightMargin(0.10);
   	pad11->Draw();
   	pad11->cd();
//    	pad11->SetTicks();
	if(Statistic)
		h_MC->SetStats();
	else
		h_MC->SetStats(0);
		
	if(logX)
		pad11->SetLogx();
		
	h_MC->Draw("COLZ");
// 	h_MC->Draw("SAME TEXT0");
	h_MC->SetTitle(name);
	h_MC->GetYaxis()->SetTitleOffset(0.65);
	h_MC->GetXaxis()->SetTitleOffset(0.6);
	h_MC->GetXaxis()->SetTitle(name_Xaxis);
	h_MC->GetYaxis()->SetTitle(name_Yaxis);
			
	TString save_name_pdf = name + ".pdf";
	TString save_name_png = name + ".png";
	
// 	c1->Print(save_name_pdf);
// 	c1->Print(save_name_png);
	c1->Print(save);
	
}

void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, float Y_min, TString strig_1 = "MC", TString strig_2 = "DATA", TString strig_3, TString strig_4, TString strig_5){
			

	TLatex lat;
	TCanvas *c1 = new TCanvas(name, name,  500, 500);//210,45,1050,750);
   	c1->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();

	if(logx)
		pad11->SetLogx();
	pad11->cd();
	
	if(logy)
		pad11->SetLogy();
	pad11->cd();
	
	Double_t max;
	Double_t min;

	if(logy) min = 10e-5;
	else min = 0;
	
	min = Y_min;

	TH1F* last_hist = (TH1F *)h_MC->GetStack()->Last();	
		
	if(h_DATA->GetMaximum() > last_hist->GetMaximum())
		max = h_DATA->GetMaximum();
	else
		max = last_hist->GetMaximum();
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;

	h_MC->SetTitle(name);
	last_hist->SetTitle(name);
	last_hist->Draw();
	last_hist->SetStats(1);
	c1->Update();

// 	TPaveStats * st_MC = (TPaveStats *)last_hist->GetListOfFunctions()->FindObject("stats");
//     if( st_MC ){ 
// 		st_MC->SetName("Const");
// 		st_MC->SetX1NDC(0.75);
// 		st_MC->SetY1NDC(0.52);
// 		st_MC->SetY2NDC(0.72);
// 		st_MC->SetTextColor(kBlue);
//     }
//     else std::cout << "Null pointer to TPaveStats MC: " << st_MC << std::endl;
  
	h_DATA->SetLineColor(kBlack);
	h_DATA->SetTitle(name);
	h_DATA->Draw();
	h_DATA->SetStats(0);
	TH1F *ratio = (TH1F*) h_DATA->Clone();
	h_DATA->SetStats(1);
	c1->Update();
// 	gPad->Update();

// 	TPaveStats * st_DATA = (TPaveStats *)h_DATA->GetListOfFunctions()->FindObject("stats");
//     if( st_DATA ){ 
//     	st_DATA->SetTextColor(kRed); 
// 		st_DATA->SetName("Const");
// 		st_DATA->SetX1NDC(0.75);
// 		st_DATA->SetY1NDC(0.75);
// 		st_DATA->SetY2NDC(0.95);
// 
//     }
//     else std::cout << "Null pointer to TPaveStats DATA: " << st_DATA << std::endl;

	last_hist->GetYaxis()->SetTitleOffset(1.2);
	h_DATA->GetYaxis()->SetTitleOffset(1.2);

	last_hist->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h_DATA->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
		
    c1->Update();

	h_MC->Draw();
	h_DATA->Draw("samePE");
	h_DATA->SetMarkerStyle(3);
	h_DATA->SetMarkerColor(kBlack);
	h_DATA->SetMarkerSize(0.5);
// 	st_MC->Draw("same");
// 	st_DATA->Draw("same");
//     	
	TPaveLabel *anno = new TPaveLabel(0.60,0.80,0.80,0.88, "2016", "NBNDC");
// 	TPaveLabel *anno = new TPaveLabel(0.60,0.80,0.80,0.88, "2017", "NBNDC");
// 	anno->SetTextAlign(12);
// 	anno->SetTextFont(62);
	anno->SetTextSize(0.6);
// 	anno->SetFillColor(0);
// 	anno->SetFillStyle(0);
// 	anno->SetBorderSize(0);
	anno->Draw();

// 	TLegend *l1 = new TLegend(0.3,0.8,0.5,0.9);
// 	l1->AddEntry(h_MC, strig_1, "l");
// 	l1->AddEntry(h_DATA, strig_2, "l");
// 	l1->Draw();	
	h_MC->GetYaxis()->SetTitleOffset(1.25);
	h_DATA->GetYaxis()->SetTitleOffset(1.25);
	h_MC->GetYaxis()->SetTitle("Entries");
	h_DATA->GetYaxis()->SetTitle("Entries");
	h_MC->GetYaxis()->SetTitleSize(0.035);
	h_DATA->GetYaxis()->SetTitleSize(0.035);
	legend->Draw();
		
	c1->Update();
	c1->cd();
	
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.25);
// 	pad22->SetTopMargin(0);
	pad22->SetTopMargin(0.99);
	pad22->SetBottomMargin(0.35);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();

	if(logx)
		pad22->SetLogx();
	pad22->cd();
	
	ratio->Divide(last_hist);
	ratio->GetYaxis()->SetRangeUser(0, 2);
	ratio->SetTitle("");
	ratio->SetStats(0);	
	
	ratio->GetYaxis()->SetTitleSize(0.1);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("DATA / MC");
	ratio->GetYaxis()->SetTitleOffset(0.50);
	ratio->GetYaxis()->SetNdivisions(506); 
	ratio->GetXaxis()->SetTitle(name_axis);
	ratio->GetXaxis()->SetLabelSize(0.15);
	ratio->GetXaxis()->SetTitleSize(0.2);
	ratio->GetXaxis()->SetTitleOffset(0.75);
	ratio->SetLineColor(kBlack);	
	ratio->Draw();
	c1->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kGreen);
	line->SetLineWidth(1);
	line->Draw();

// 	TF1* f1 = new TF1("f1", "pol1", pad22->GetUxmin(), pad22->GetUxmax());
// 	ratio->Fit("f1","R");
// 
//     TLatex* latexFit = new TLatex();
//     for(int i = 0; i < f1->GetNpar()+1; i++){
//        	latexFit->SetTextSize(0.1);
//     	if(i == 2){
//     		float yPos = 0.8;
// 	    	TString longstring = Form("#chi^{2} = %5.3g", f1->GetChisquare());
//   	       	latexFit->DrawLatex(pad22->GetUxmin()+fabs(pad22->GetUxmin()/(float)10), yPos, longstring);
//   	    }   
//     	float yPos = 1.2 + i*0.5;
//     	TString longstring = Form("%s = %5.3g #pm %5.3g", f1->GetParName(i),f1->GetParameter(i),f1->GetParError(i));
//        	latexFit->DrawLatex(pad22->GetUxmin()+fabs(pad22->GetUxmin()/(float)10), yPos, longstring);
//     }
	TString save_name_pdf = name + ".pdf";
	TString save_name_png = name + ".png";
	
// 	c1->Print(save_name_pdf);
// 	c1->Print(save_name_png);
	c1->Print(save);
	



}

void SalvaHisto(TString name, TH1F* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, float Y_min, TString strig_1, TString strig_2, TString strig_3, TString strig_4, TString strig_5){

	TLatex lat;
	TCanvas *c1 = new TCanvas(name, name,  500, 500);//210,45,1050,750);
	
   	c1->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();

	if(logx)
		pad11->SetLogx();
	pad11->cd();
	
	if(logy)
		pad11->SetLogy();
	pad11->cd();
	
	Double_t max;
	Double_t min;

	if(logy) min = 0.5;
	else min = 0;

	if(h_DATA->GetMaximum() > h_MC->GetMaximum())
		max = h_DATA->GetMaximum();
	else
		max = h_MC->GetMaximum();
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;
	
	h_MC->SetLineColor(kBlue);
	h_MC->SetTitle(name);
	h_MC->Draw();
	h_MC->SetStats(1);
	c1->Update();

// 	TPaveStats * st_MC = (TPaveStats *)h_MC->GetListOfFunctions()->FindObject("stats");
//     if( st_MC ){ 
// 		st_MC->SetName("Const");
// 		st_MC->SetX1NDC(0.75);
// 		st_MC->SetY1NDC(0.52);
// 		st_MC->SetY2NDC(0.72);
// 		st_MC->SetTextColor(kBlue);
//     }
//     else std::cout << "Null pointer to TPaveStats MC: " << st_MC << std::endl;
  
	h_DATA->SetLineColor(kRed);
	h_DATA->SetTitle(name);
	h_DATA->Draw();
// 	h_DATA->SetStats(0);
	TH1F *ratio = (TH1F*) h_DATA->Clone();
	h_DATA->SetStats(1);
	c1->Update();
// 	gPad->Update();

// 	TPaveStats * st_DATA = (TPaveStats *)h_DATA->GetListOfFunctions()->FindObject("stats");
//     if( st_DATA ){ 
//     	st_DATA->SetTextColor(kRed); 
// 		st_DATA->SetName("Const");
// 		st_DATA->SetX1NDC(0.75);
// 		st_DATA->SetY1NDC(0.75);
// 		st_DATA->SetY2NDC(0.95);
// 
//     }
//     else std::cout << "Null pointer to TPaveStats DATA: " << st_DATA << std::endl;

	h_MC->GetYaxis()->SetTitleOffset(1.2);
	h_DATA->GetYaxis()->SetTitleOffset(1.2);

	h_MC->SetMaximum(max);
	h_MC->SetMinimum(min);	
// 	h_MC->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h_DATA->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);

		
    c1->Update();

	h_MC->Draw();
	h_DATA->Draw("sameP");
	h_DATA->SetMarkerStyle(3);
	h_DATA->SetMarkerColor(kRed);
	h_DATA->SetMarkerSize(0.5);
// 	st_MC->Draw("same");
// 	st_DATA->Draw("same");
	
//     	
	TLegend *l1 = new TLegend(0.4,0.8,0.6,0.9);
	l1->AddEntry(h_MC, strig_1, "l");
	l1->AddEntry(h_DATA, strig_2, "l");
	l1->Draw();	
	
	c1->Update();
	c1->cd();
	
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.25);
// 	pad22->SetTopMargin(0);
	pad22->SetTopMargin(0.95);
	pad22->SetBottomMargin(0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();

	if(logx)
		pad22->SetLogx();
	pad22->cd();
	
	ratio->Divide(h_MC);
	ratio->GetYaxis()->SetRangeUser(0, 2);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("DATA / MC");
	ratio->GetYaxis()->SetTitleOffset(0.50);
	ratio->GetYaxis()->SetNdivisions(506); 
	ratio->GetXaxis()->SetTitle(name_axis);
	ratio->GetXaxis()->SetLabelSize(0.15);
	ratio->GetXaxis()->SetTitleSize(0.15);
	ratio->GetXaxis()->SetTitleOffset(0.75);
	ratio->SetLineColor(kBlack);	
	ratio->Draw();
	c1->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kGreen);
	line->SetLineWidth(1);
	line->Draw();
	
		
	TString save_name_pdf = name + ".pdf";
	TString save_name_png = name + ".png";
	
// 	c1->Print(save_name_pdf);
// 	c1->Print(save_name_png);
	c1->Print(save);
}
