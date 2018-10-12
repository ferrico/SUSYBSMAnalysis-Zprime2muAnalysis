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
  float vertex_chi2;
  
  float met_pt;
  
 Color_t color;
 float gM;
 float kFactor;
 float kFactor_BB;
 float kFactor_BE;
 float NNPDF;
 
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

    TString samples_2017[15] =  {
                            "Wantitop", "tW", 
                            "ttbar",
                            "WW",
                            "ZZ",
                            "WZ",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
                            };
	
	float events_2017[15] = {
						7780870, 7581624,
						33844772,
						7791498,
						1949768,
						3928630,
						2961000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000,
					};
						
	float sigma_2017[15] = {
						35.6, 35.6,
						831.76,
						118.7,
						16.523,
						47.13,
						2112.905, 20.553, 2.8861, 0.25126, 0.017075, 1.366E-3, 8.178E-5, 3.191E-6, 2.787E-7,
						};         
         
	float LUMINOSITY_2018 = 27697.269;
	float LUMINOSITY_2017 = 41903.837;
	float LUMINOSITY_2016 = 36235.493;
	float LUMINOSITY = LUMINOSITY_2018 + LUMINOSITY_2017 + LUMINOSITY_2016;
	
	float Z_peak = 0.9638;
	float Z_peak_BB = 0.9688;
	float Z_peak_BE = 0.9610;

	float weight[45] = {0};
	float weight_BB[45] = {0};
	float weight_BE[45] = {0};

	float weight_2017[15] = {0};
	float weight_2017_BB[15] = {0};
	float weight_2017_BE[15] = {0};

	const int    NMBINS = 100;
	const double MMIN = 60., MMAX = 2100.;
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
	THStack* pt_MC_Stack = new THStack("MC_pt", "MC_pt");
	THStack* pt_MC_BB_Stack = new THStack("MC_pt_BB", "MC_pt_BB");
	THStack* pt_MC_BE_Stack = new THStack("MC_pt_BE", "MC_pt_BE");
	TH1F* pt_MC_clear;
	TH1F* pt_MC_BB_clear;
	TH1F* pt_MC_BE_clear;

	TH1F* pt_cumulative_DATA = new TH1F("DATA_pt_cumulative", "DATA_pt_cumulative", binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_BB = new TH1F("DATA_pt_cumulative_BB", "DATA_pt_cumulative_BB",binnum_pt, PT_BINS);
	TH1F* pt_cumulative_DATA_BE = new TH1F("DATA_pt_cumulative_BE", "DATA_pt_cumulative_BE",binnum_pt, PT_BINS);
	THStack* pt_cumulative_MC_Stack = new THStack("MC_pt_cumulative", "MC_pt_cumulative");
	THStack* pt_cumulative_MC_BB_Stack = new THStack("MC_pt_cumulative_BB", "MC_pt_cumulative_BB");
	THStack* pt_cumulative_MC_BE_Stack = new THStack("MC_pt_cumulative_BE", "MC_pt_cumulative_BE");
	TH1F* pt_cumulative_MC_clear;
	TH1F* pt_cumulative_MC_BB_clear;
	TH1F* pt_cumulative_MC_BE_clear;
	
	TH1F* eta_DATA = new TH1F("DATA_eta", "DATA_eta", binnum_eta, ETA_BINS);
	THStack* eta_MC_Stack = new THStack("MC_eta", "MC_eta");
	TH1F* eta_MC_clear;

	TH1F* phi_DATA = new TH1F("DATA_phi", "DATA_phi", binnum_phi, PHI_BINS);
	THStack* phi_MC_Stack = new THStack("MC_phi", "MC_phi");
	TH1F* phi_MC_clear;

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

void DataMC_comparison(){

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


	for (int ibin = 0; ibin <= NMBINS; ibin++)
    	logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);

	TH1F* dimuon_DATA = new TH1F("DATA_dimuon", "DATA_dimuon", NMBINS, logMbins);
	THStack* dimuon_MC_Stack = new THStack("MC_dimuon", "MC_dimuon");
	TH1F* dimuon_MC_clear;

	TH1F* dimuon_BB_DATA = new TH1F("DATA_dimuon_BB", "DATA_dimuon_BB", NMBINS, logMbins);
	THStack* dimuon_BB_MC_Stack = new THStack("MC_dimuon_BB", "MC_dimuon_BB");
	TH1F* dimuon_BB_MC_clear;

	TH1F* dimuon_BE_DATA = new TH1F("DATA_dimuon_BE", "DATA_dimuon_BE", NMBINS, logMbins);
	THStack* dimuon_BE_MC_Stack = new THStack("MC_dimuon_BE", "MC_dimuon_BE");
	TH1F* dimuon_BE_MC_clear;

	TH1F* dimuon_cumulative_DATA = new TH1F("DATA_dimuon_cumulative", "DATA_dimuon_cumulative", NMBINS, logMbins);
	THStack* dimuon_cumulative_MC_Stack = new THStack("MC_dimuon_cumulative", "MC_dimuon_cumulative");
	TH1F* dimuon_cumulative_MC_clear;

	TH1F* dimuon_cumulative_BB_DATA = new TH1F("DATA_dimuon_cumulative_BB", "DATA_dimuon_cumulative_BB", NMBINS, logMbins);
	THStack* dimuon_cumulative_BB_MC_Stack = new THStack("MC_dimuon_cumulative_BB", "MC_dimuon_cumulative_BB");
	TH1F* dimuon_cumulative_BB_MC_clear;

	TH1F* dimuon_cumulative_BE_DATA = new TH1F("DATA_dimuon_cumulative_BE", "DATA_dimuon_cumulative_BE", NMBINS, logMbins);
	THStack* dimuon_cumulative_BE_MC_Stack = new THStack("MC_dimuon_cumulative_BE", "MC_dimuon_cumulative_BE");
	TH1F* dimuon_cumulative_BE_MC_clear;

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

	for(int i = 0; i < 15; i++){
	
		weight_2017[i] = LUMINOSITY_2017 * sigma_2017[i] / events_2017[i];
		weight_2017[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
		weight_2017_BB[i] = weight_2017[i] * Z_peak_BB / Z_peak;
		weight_2017_BE[i] = weight_2017[i] * Z_peak_BE/ Z_peak;
	}	
       
    std::cout<<"PASSO AL 2017"<<std::endl;
       

	for(int j = 0; j < 15; j++){

        std::cout<<"opening.. "<<samples_2017[j]<<" --- "<<events_2017[j]<<std::endl;

        TChain *treeMC = new TChain("SimpleNtupler/t");
//         treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_Pt_Assignment/ana_datamc_" + samples[j] + ".root");

        treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/ana_datamc_" + samples_2017[j] + ".root");        
//      	TFile *file_MC = new TFile("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/ana_datamc_" + samples[j] + ".root", "READ");
//      	file_MC->cd("Our2016MuonsPlusMuonsMinusHistos");
		

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

     //treeMC->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeMC->SetBranchAddress("met_pt",&met_pt); 	
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples_2017[j].Data(),j , weight_2017[j], nentries);

		dimuon_MC_clear = new TH1F("MC_dimuon", "MC_dimuon", NMBINS, logMbins);
		dimuon_BB_MC_clear = new TH1F("MC_dimuon_BB", "MC_dimuon_BB", NMBINS, logMbins);
		dimuon_BE_MC_clear = new TH1F("MC_dimuon_BE", "MC_dimuon_BE", NMBINS, logMbins);

		dimuon_cumulative_MC_clear = new TH1F("MC_dimuon_cumulative", "MC_dimuon_cumulative", NMBINS, logMbins);		
		dimuon_cumulative_BB_MC_clear = new TH1F("MC_dimuon_cumulative_BB", "MC_dimuon_cumulative_BB", NMBINS, logMbins);
		dimuon_cumulative_BE_MC_clear = new TH1F("MC_dimuon_cumulative_BE", "MC_dimuon_cumulative_BE", NMBINS, logMbins);
		
		pt_MC_clear = new TH1F("MC_pt", "MC_pt", binnum_pt, PT_BINS);
		pt_MC_BB_clear = new TH1F("MC_pt_BB", "MC_pt_BB", binnum_pt, PT_BINS);
		pt_MC_BE_clear = new TH1F("MC_pt_BE", "MC_pt_BE", binnum_pt, PT_BINS);

		pt_cumulative_MC_clear = new TH1F("MC_pt_cumulative_cumulative", "MC_pt_cumulative_cumulative", binnum_pt, PT_BINS);
		pt_cumulative_MC_BB_clear = new TH1F("MC_pt_cumulative_BB", "MC_pt_cumulative_BB", binnum_pt, PT_BINS);
		pt_cumulative_MC_BE_clear = new TH1F("MC_pt_cumulative_BE", "MC_pt_cumulative_BE", binnum_pt, PT_BINS);

		eta_MC_clear = new TH1F("MC_eta", "MC_eta", binnum_eta, ETA_BINS);
		
		phi_MC_clear = new TH1F("MC_phi", "MC_phi", binnum_phi, PHI_BINS);
		
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

    	 	if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
    	 	
    	 	treeMC->GetEntry(p);
    	 	
// 			if(j > 30) std::cout<<genWeight;
			
		if(
			GoodVtx && 
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
				
				gM = gen_dil_mass - 400;
				kFactor = 1.047 - 0.000143 * gM + 5.167e-08 * pow(gM,2) - 7.84e-12 * pow(gM,3);
			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
	 		   	NNPDF = 0.9803 -0.0001068 * gM + 1.142e-07 * pow(gM,2) -2.013e-11 * pow(gM,3)-5.904e-15 * pow(gM,4)+1.634e-18 * pow(gM,5);
	 		   	
				if(j >= 6 && j <=14){
					weight[j] *= kFactor*NNPDF;
					weight_BB[j] *= kFactor_BB*NNPDF;
					weight_BE[j] *= kFactor_BE*NNPDF;
				}

					dimuon_MC_clear->Fill(dil_mass,  weight_2017[j]); 
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_MC_clear->Fill(dil_mass,  weight_2017_BB[j]);
					else
						dimuon_BE_MC_clear->Fill(dil_mass,  weight_2017_BE[j]);
						
					for(int i = 1; i <  NMBINS+1; i++){
						float contenuto_mc = dimuon_MC_clear->Integral(i, NMBINS);
						dimuon_cumulative_MC_clear->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BB_MC_clear->Integral(i, NMBINS);
						dimuon_cumulative_BB_MC_clear->SetBinContent(i, contenuto_mc);
						contenuto_mc = dimuon_BE_MC_clear->Integral(i, NMBINS);
						dimuon_cumulative_BE_MC_clear->SetBinContent(i, contenuto_mc);
					}

					Vtx_MC_clear->Fill(vertex_chi2,  weight_2017[j]); 
   					
					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_MC_clear->Fill(lep_pt[h],  weight_2017[j]); 
						if(fabs(lep_eta[h]) < 1.2)
							pt_MC_BB_clear->Fill(lep_pt[h],  weight_2017_BB[j]);
						else
							pt_MC_BE_clear->Fill(lep_pt[h],  weight_2017_BE[j]);
							
						for(int i = 1; i <  binnum_pt+1; i++){
							float contenuto_mc = pt_MC_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_BB_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_BB_clear->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_MC_BE_clear->Integral(i, binnum_pt);
							pt_cumulative_MC_BE_clear->SetBinContent(i, contenuto_mc);
						}

						eta_MC_clear->Fill(lep_eta[h],  weight_2017[j]); 

						phi_MC_clear->Fill(lep_phi[h],  weight_2017[j]); 
                                                
						dB_MC_clear->Fill(lep_dB[h],  weight_2017[j]); 

						PixelHit_MC_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight_2017[j]); 

						TkLayer_MC_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight_2017[j]); 

						Iso_MC_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight_2017[j]); 

						relpTErr_MC_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight_2017[j]); 

						ValidMu_MC_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight_2017[j]); 
   	 						    	 						
	     			} // for on muons

				if(j >= 6 && j <=14){
					weight[j] /= kFactor;
					weight_BB[j] /= kFactor_BB;
					weight_BE[j] /= kFactor_BE;
					weight[j] /= NNPDF;
					weight_BB[j] /= NNPDF;
					weight_BE[j] /= NNPDF;
				}
			
			weight_2017[j] /= genWeight;
			weight_2017_BB[j] /= genWeight;
			weight_2017_BE[j] /= genWeight;
	       }//end of condition on event

		 }// end loop p
	 
	 	 if(j < 2){
	 	 	color = kBlue+2;
	 	 	if(j == 1) legend->AddEntry(dimuon_MC_clear, "SingleTop", "f");
	 	 }
	 	 if(j == 2){
	 	 	color = 4;
	 	 	legend->AddEntry(dimuon_MC_clear, "t#bar{t}", "f");
	 	 }
	 	 if(j == 3){
	 	 	color = kRed-10;
	 	 	legend->AddEntry(dimuon_MC_clear, "WW", "f");
	 	 }
	 	 if(j == 4){
	 	 	color = 2;
	 	 	legend->AddEntry(dimuon_MC_clear, "ZZ", "f");
	 	 }
	 	 if(j == 5){
	 	 	color = kRed+3;
	 	 	legend->AddEntry(dimuon_MC_clear, "WZ", "f");
	 	 }
	 	 if(j >= 6 && j < 14){ 
	 	 	color = 3;
	 	 	if(j == 10) legend->AddEntry(dimuon_MC_clear, "DY #rightarrow #mu#mu", "f");
	 	 }

		 dimuon_MC_clear->SetFillColor(color);
		 dimuon_MC_clear->SetLineColor(color);
		 dimuon_BB_MC_clear->SetFillColor(color);
		 dimuon_BB_MC_clear->SetLineColor(color);
		 dimuon_BE_MC_clear->SetFillColor(color);
		 dimuon_BE_MC_clear->SetLineColor(color);
		 dimuon_cumulative_MC_clear->SetFillColor(color);
		 dimuon_cumulative_MC_clear->SetLineColor(color);
		 dimuon_cumulative_BB_MC_clear->SetFillColor(color);
		 dimuon_cumulative_BB_MC_clear->SetLineColor(color);
		 dimuon_cumulative_BE_MC_clear->SetFillColor(color);
		 dimuon_cumulative_BE_MC_clear->SetLineColor(color);
		pt_MC_clear->SetFillColor(color);
		pt_MC_clear->SetLineColor(color);
		pt_MC_BB_clear->SetFillColor(color);
		pt_MC_BB_clear->SetLineColor(color);
		pt_MC_BE_clear->SetFillColor(color);
		pt_MC_BE_clear->SetLineColor(color);
		pt_cumulative_MC_clear->SetFillColor(color);
		pt_cumulative_MC_clear->SetLineColor(color);
		pt_cumulative_MC_BB_clear->SetFillColor(color);
		pt_cumulative_MC_BB_clear->SetLineColor(color);
		pt_cumulative_MC_BE_clear->SetFillColor(color);
		pt_cumulative_MC_BE_clear->SetLineColor(color);
		eta_MC_clear->SetFillColor(color);
		eta_MC_clear->SetLineColor(color);
		phi_MC_clear->SetFillColor(color);
		phi_MC_clear->SetLineColor(color);
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
		
		dimuon_MC_Stack->Add(dimuon_MC_clear, "HIST");
		dimuon_BB_MC_Stack->Add(dimuon_BB_MC_clear, "HIST");
		dimuon_BE_MC_Stack->Add(dimuon_BE_MC_clear, "HIST");
		dimuon_cumulative_MC_Stack->Add(dimuon_cumulative_MC_clear, "HIST");
		dimuon_cumulative_BB_MC_Stack->Add(dimuon_cumulative_BB_MC_clear, "HIST");
		dimuon_cumulative_BE_MC_Stack->Add(dimuon_cumulative_BE_MC_clear, "HIST");
		pt_MC_Stack->Add(pt_MC_clear, "HIST");
		pt_MC_BB_Stack->Add(pt_MC_BB_clear, "HIST");
		pt_MC_BE_Stack->Add(pt_MC_BE_clear, "HIST");
		pt_cumulative_MC_Stack->Add(pt_cumulative_MC_clear, "HIST");
		pt_cumulative_MC_BB_Stack->Add(pt_cumulative_MC_BB_clear, "HIST");
		pt_cumulative_MC_BE_Stack->Add(pt_cumulative_MC_BE_clear, "HIST");
		eta_MC_Stack->Add(eta_MC_clear, "HIST");
		phi_MC_Stack->Add(phi_MC_clear, "HIST");
		dB_MC_Stack->Add(dB_MC_clear, "HIST");
		PixelHit_MC_Stack->Add(PixelHit_MC_clear, "HIST");
		TkLayer_MC_Stack->Add(TkLayer_MC_clear, "HIST");
		Iso_MC_Stack->Add(Iso_MC_clear, "HIST");
		relpTErr_MC_Stack->Add(relpTErr_MC_clear, "HIST");
		Vtx_MC_Stack->Add(Vtx_MC_clear, "HIST");
		ValidMu_MC_Stack->Add(ValidMu_MC_clear, "HIST");
			 	 
    }// end loop on MC 2017



////////      DATA    ///////

	for(int k = 0; k < 2; k++){
		
    TChain *treeDATA = new TChain("SimpleNtupler/t");
    if(k == 10) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/DATA_Pt_Assignment/ana_datamc_data_Pt_Assignment.root");
    if(k == 1) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA/ana_datamc_data.root");
    if(k == 22) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2018/DATA/ana_datamc_data.root");
     	    	
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


     //treeDATA->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeDATA->SetBranchAddress("met_pt",&met_pt);     

	Long64_t nentries = treeDATA->GetEntries();
	
	printf("opening... DATA --- %lld\n", nentries);    
		int bhu_data = 0;
	
	for(int p=0; p<nentries; p++){
// 	for(int p=0; p<1000000; p++){
// 	for(int p=0; p<10000; p++){
// 	for(int p=0; p<1; p++){

		if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
	
		treeDATA->GetEntry(p);
		
		if(
			GoodVtx && 
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
						for(int i = 1; i <  binnum_pt+1; i++){
							float contenuto_mc = pt_DATA->Integral(i, binnum_pt);
							pt_cumulative_DATA->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_BB->Integral(i, binnum_pt);
							pt_cumulative_DATA_BB->SetBinContent(i, contenuto_mc);
							contenuto_mc = pt_DATA_BE->Integral(i, binnum_pt);
							pt_cumulative_DATA_BE->SetBinContent(i, contenuto_mc);
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

	save_document = dir_save + "Variuos_distribution.pdf";

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf[", 0,  0, "p_{T}", 0);

	SalvaHisto("dimuon Mass: stack plot", dimuon_MC_Stack, dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon Mass BB: stack plot", dimuon_BB_MC_Stack, dimuon_BB_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon Mass BE: stack plot", dimuon_BE_MC_Stack, dimuon_BE_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);

	SalvaHisto("dimuon cumulative Mass: stack plot", dimuon_cumulative_MC_Stack, dimuon_cumulative_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon cumulative Mass BB: stack plot", dimuon_cumulative_BB_MC_Stack, dimuon_cumulative_BB_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);
	SalvaHisto("dimuon cumulative Mass BE: stack plot", dimuon_cumulative_BE_MC_Stack, dimuon_cumulative_BE_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]", 10e-5);

	SalvaHisto("Lepton p_{T}", pt_MC_Stack, pt_DATA, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BB", pt_MC_BB_Stack, pt_DATA_BB, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BE", pt_MC_BE_Stack, pt_DATA_BE, save_document, 0,  1, "p_{T}", 10e-5);

	SalvaHisto("Lepton p_{T} cumulative", pt_cumulative_MC_Stack, pt_cumulative_DATA, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BB cumulative", pt_cumulative_MC_BB_Stack, pt_cumulative_DATA_BB, save_document, 0,  1, "p_{T}", 10e-5);
	SalvaHisto("Lepton p_{T} BE cumulative", pt_cumulative_MC_BE_Stack, pt_cumulative_DATA_BE, save_document, 0,  1, "p_{T}", 10e-5);

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
	TLegend *l1 = new TLegend(0.3,0.8,0.5,0.9);
	l1->AddEntry(h_MC, strig_1, "l");
	l1->AddEntry(h_DATA, strig_2, "l");
// 	l1->Draw();	

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
	
// 	lat.DrawLatex(X_pos_latex, max, strig_1);
// 	lat.DrawLatex(X_pos_latex, 0.9*max, strig_2);
// 	lat.DrawLatex(X_pos_latex, 0.8*max, strig_3);
// 	lat.DrawLatex(X_pos_latex, 0.7*max, strig_4);
// 	lat.DrawLatex(X_pos_latex, 0.6*max, strig_5);
	
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
