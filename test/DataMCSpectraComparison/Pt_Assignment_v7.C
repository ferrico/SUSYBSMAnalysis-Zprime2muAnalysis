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
void SalvaHisto(TString name, TH1F* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, TH2F* h_MC, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic = 1, bool logX = false);
void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA[15], TString save, int a, int b, int c, int d, int e);

  Int_t event;
  Int_t run;
  Int_t prev_event = -88;
//   Int_t lumi;

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

	TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

	
    TString samples[45] =  {"dyInclusive50",
     						"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200",
                            "Wjets",
                            "Wantitop", "tW", 
                            "ttbar_50to500", "ttbar_500to800", "ttbar_800to1200", "ttbar_1200to1800", "ttbar_1800toInf", 							
                            "WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500",
                            "ZZ", "ZZ_ext", 
                            "WZ", "WZ_ext",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
	 						"dyPt0to50", "dyPt50to100", "dyPt100to250", "dyPt250to400", "dyPt400to650", "dyPt650",
                            };
	
	float events[45] = {19385554,  
						6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
						29705748,
						6933094, 6952830,
						79092400, 200000, 199800, 200000, 40829, 
						1999000, 200000, 200000, 200000, 38969, 
						990064, 998034, 
						1000000, 2995828,
						2977600, 100000, 100000, 98400, 100000, 95106, 100000, 100000, 100000,
						22782948, (float)39612900, (float)26998200, (float)7190820, 167272, 177101
// 						878212, 1662453, 2833172, 1571199, 48731, 177101					
						};
						
	float sigma[45] = {6025.2, 
						2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
						61526.7,
						35.6, 35.6,
						87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
						12.178, 1.385, 0.0566, 0.0035, 0.00005,
					    8.2615, 8.2615, 
					    23.565, 23.565, //16.523, 47.13, 16.523, 47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7,
						5352.58, 363.81, 84.015, 3.2283, 0.43604, 0.04098,
						};
         
	float LUMINOSITY = 36235.493;
	
	float Z_peak = 0.9638;

	float weight[45] = {0};
	
	int count_event_MC[6][32] = {0};
	int count_event_DATA[6] = {0};
		
	float lepton_pt[2][8] = {0};
	float lepton_eta[2][8] = {0};
	float lepton_phi[2][8] = {0};
	int numberOfValidMuonHits[2][8];
	int NotGoodQualityMuon_MC = 0;
	int NotGoodQualityMuon_DATA = 0;
	float count_assignment_MC_till400[8][32];
	float count_assignment_MC_400to600[8][32];
	float count_assignment_MC_above600[8][32];
	float count_assignment_DATA_till400[8];
	float count_assignment_DATA_400to600[8];
	float count_assignment_DATA_above600[8];

	float tot_till[8] = {0};
	float tot_middle[8] = {0};
	float tot_above[8] = {0};
	float TOT_MC = 0;
	float TOT_DATA = 0;
	
	TString reconstruction[8] = {"Selected", "global", "picky", "tpfms", "dyt", "tracker", "std", "tuneP"};
	TString pt_name[3] = {"MC_pt_", "DATA_pt_", "GEN_pt_"};
	TString eta_name[3] = {"MC_eta_", "DATA_eta_", "GEN_eta"};
	TString phi_name[3] = {"MC_phi_", "DATA_phi_", "GEN_phi"};
	TString pt_vs_eta_name[2] = {"MC_pt_vs_eta_", "DATA_ptVSeta_"};//DATA_pt_vs_eta_"};
	TString pt_vs_phi_name[2] = {"MC_pt_vs_phi_", "DATA_ptVSphi_"};//DATA_pt_vs_phi_"};
	TString eta_vs_phi_name[2] = {"MC_eta_vs_phi_", "DATA_etaVSphi_"};//DATA_eta_vs_phi_"};
	TString pt_vs_met_name[2] = {"MC_pt_vs_met_", "DATA_ptVSmet_"};//DATA_pt_vs_met_"};
	TH1F* pt_MC[8];
	TH1F* pt_MC_clear[8];
	THStack* pt_MC_Stack[8];
	TH1F* pt_DATA[8];
	TH1F* eta_MC[8];
	TH1F* eta_MC_clear[8];
	THStack* eta_MC_Stack[8];
	TH1F* eta_DATA[8];
	TH1F* phi_MC[8];
	TH1F* phi_MC_clear[8];
	THStack* phi_MC_Stack[8];
	TH1F* phi_DATA[8];
	TH2F* pt_vs_eta_MC[8];
	TH2F* pt_vs_eta_DATA[8];
	TH2F* pt_vs_phi_MC[8];
	TH2F* pt_vs_phi_DATA[8];
	TH2F* eta_vs_phi_MC[8];
	TH2F* eta_vs_phi_DATA[8];
	TH2F* pt_vs_met_MC[8];
	TH2F* pt_vs_met_DATA[8];
	TH1F* Double_Count_MC[6][6];
	TH1F* Double_Count_DATA[6][6];

	TH1F* eta_MC_till400[8];
	TH1F* eta_MC_till400_clear[8];
	THStack* eta_MC_till400_Stack[8];
	TH1F* eta_DATA_till400[8];
	TH1F* eta_MC_400to600[8];
	TH1F* eta_MC_400to600_clear[8];
	THStack* eta_MC_400to600_Stack[8];
	TH1F* eta_DATA_400to600[8];
	TH1F* eta_MC_above600[8];
	TH1F* eta_MC_above600_clear[8];
	THStack* eta_MC_above600_Stack[8];
	TH1F* eta_DATA_above600[8];
	
	TH1F* phi_MC_till400[8];
	TH1F* phi_MC_till400_clear[8];
	THStack* phi_MC_till400_Stack[8];
	TH1F* phi_DATA_till400[8];
	TH1F* phi_MC_400to600[8];
	TH1F* phi_MC_400to600_clear[8];
	THStack* phi_MC_400to600_Stack[8];
	TH1F* phi_DATA_400to600[8];
	TH1F* phi_MC_above600[8];
	TH1F* phi_MC_above600_clear[8];
	THStack* phi_MC_above600_Stack[8];
	TH1F* phi_DATA_above600[8];
	
	TH1D* h_count_assignment_MC_till400[8];
	TH1D* h_count_assignment_DATA_till400[8];
	
	TH1F* pt_CountDouble_MC[5][15];
	TH1F* pt_CountDouble_DATA[5][15];
	
	TH1F* DileptonPt_MC_clear[8];
	THStack* DileptonPt_MC_Stack[8];
	TH1F* DileptonPt_DATA[8];
	
	TH2F* Phi_bin = new TH2F("Check #phi bin", "Check #phi bin", 40, 0, 0.8, 80, 0, 0.8);
	TH2F* Phi_bin_46 = new TH2F("Check #phi bin: 400-600", "Check #phi bin: 400-600", 40, 0, 0.8, 80, 0, 0.8);
	TH2F* Phi_bin_6 = new TH2F("Check #phi bin: above 600", "Check #phi bin: above 600", 40, 0, 0.8, 80, 0, 0.8);
// 	Phi_bin->GetXaxis()->SetTile("reco #phi");
// 	Phi_bin->GetYaxis()->SetTile("gen #phi");
	
	const int    NMBINS = 100;
	const double MMIN = 60., MMAX = 2100.;
	double logMbins[NMBINS+1];

//     Double_t PT_BINS[] = {200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 600, 750, 1000, 1500};
//     Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;
    Double_t PT_BINS[] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 
    					1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800};
    Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;    
    
//     Double_t ETA_BINS[] = {-2.4, -1.5, -1.2, -0.9, 0, 0.9, 1.2, 1.5, 2.4};
    Double_t ETA_BINS[] = {0, 0.9, 1.2, 1.5, 2.4};//, 5, 5.9, 6.2, 6.5, 7.4};
    Int_t  binnum_eta = sizeof(ETA_BINS)/sizeof(Double_t)-1;
    
    Double_t PHI_BINS[] = {-3.14, -2.356, -1.57, -0.785, 0, 0.785, 1.57, 2.356, 3.14};
    Int_t  binnum_phi = sizeof(PHI_BINS)/sizeof(Double_t)-1;
    
    Double_t MET_BINS[] = {0, 25, 50, 75, 100, 200, 350, 500};//, 750, 1000};
    Int_t  binnum_met = sizeof(MET_BINS)/sizeof(Double_t)-1;
   
    int count_double_MC[8][8] = {0};
    int count_double_DATA[8][8] = {0};    
    
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 
	
	TH2F* genW_vs_eta = new TH2F("genW vs eta", "genW vs eta", 3, -1.5, 1.5, binnum_eta, ETA_BINS);
	TH2F* genW_vs_phi = new TH2F("genW vs phi", "genW vs phi", 3, -1.5, 1.5, binnum_phi, PHI_BINS);
	TH2F* genW_vs_pt = new TH2F("genW vs pt", "genW vs pt", 3, -1.5, 1.5, binnum_pt, PT_BINS);

	int countDeltaPhi_MC = 0;
	int countDeltaPhi_DATA = 0;
	int NOcountDeltaPhi_MC = 0;
	int NOcountDeltaPhi_DATA = 0;
	THStack* DeltaPhi_MC_till400_Stack = new THStack(phi_name[0] + reconstruction[7] + "_till400", phi_name[0] + reconstruction[7] + "_till400");
	THStack* DeltaPhi_MC_400to600_Stack;
	THStack* DeltaPhi_MC_above600_Stack;
	TH1F* DeltaPhi_MC_till400_clear = new TH1F(phi_name[0] + reconstruction[7] + "_till400", phi_name[0] + reconstruction[7] + "_till400", 50, 0, 6.28);
	TH1F* DeltaPhi_MC_400to600_clear;
	TH1F* DeltaPhi_MC_above600_clear;
	THStack* DeltaPhi_GEN_till400_Stack = new THStack(phi_name[2] + reconstruction[7] + "_till400", phi_name[2] + reconstruction[7] + "_till400");
	THStack* DeltaPhi_GEN_400to600_Stack;
	THStack* DeltaPhi_GEN_above600_Stack;
	TH1F* DeltaPhi_GEN_till400_clear = new TH1F(phi_name[2] + reconstruction[7] + "_till400", phi_name[2] + reconstruction[7] + "_till400", 50, 0, 6.28);
	TH1F* DeltaPhi_GEN_400to600_clear;
	TH1F* DeltaPhi_GEN_above600_clear;
	TH1F* DeltaPhi_DATA_till400 = new TH1F(phi_name[1] + reconstruction[7] + "_till400", phi_name[1] + reconstruction[7] + "_till400", 50, 0, 6.28);
	TH1F* DeltaPhi_DATA_400to600;
	TH1F* DeltaPhi_DATA_above600;


void Pt_Assignment_v7(){

	gStyle->SetOptFit(1);        
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

	for(int i = 0; i < 45; i++){
		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
	}
	
/*	for(int a = 0; a < 5; a++){
		int b = a + 1;
// 		if(a == 4) b = 6;
		pt_CountDouble_MC[a][0] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][1] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][2] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][3] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][4] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][5] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][6] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][7] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][8] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][9] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);	
		pt_CountDouble_MC[a][10] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][11] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][12] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][13] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_MC[a][14] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][0] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][1] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][2] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][3] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][4] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][5] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][6] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][7] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][8] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][9] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);	
		pt_CountDouble_DATA[a][10] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][11] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][12] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][13] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA[a][14] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
	}
*/	
	THStack *gen_Weight_Stack = new THStack("gen Weight - stack plot", "gen Weight - stack plot");

	TH1F *GenDimuon_MC = new TH1F("GenDimuon mass: MC", "GenDimuon mass: MC", NMBINS, logMbins);
	THStack *GenDimuon_MC_Stack = new THStack("GenDimuon mass: MC - stack plot", "GenDimuon mass: MC - stack plot");
	TH1F *Dimuon_MC = new TH1F("Dimuon mass: MC", "Dimuon mass: MC", NMBINS, logMbins);
	THStack *Dimuon_MC_Stack = new THStack("Dimuon_mass_MC_stack_plot", "Dimuon_mass_MC_stack_plot");

	TH1F *Dimuon_DATA = new TH1F("Dimuon mass: DATA", "Dimuon mass: DATA", NMBINS, logMbins);
// 	TH1F *PickyDimuon_MC = new TH1F("PickyDimuon mass: MC", "PickyDimuon mass: MC", NMBINS, logMbins);
	THStack *PickyDimuon_MC_Stack = new THStack("PickyDimuon mass: MC - stack plot", "PickyDimuon mass: MC - stack plot");
	TH1F *PickyDimuon_DATA = new TH1F("PickyDimuon mass: DATA", "PickyDimuon mass: DATA", NMBINS, logMbins);
// 	TH1F *DYTDimuon_MC = new TH1F("DYTDimuon mass: MC", "DYTDimuon mass: MC", NMBINS, logMbins);
	THStack *DYTDimuon_MC_Stack = new THStack("DYTDimuon mass: MC - stack plot", "DYTDimuon mass: MC - stack plot");
	TH1F *DYTDimuon_DATA = new TH1F("DYTDimuon mass: DATA", "DYTDimuon mass: DATA", NMBINS, logMbins);
// 	TH1F *TrackerDimuon_MC = new TH1F("TrackerDimuon mass: MC", "TrackerDimuon mass: MC", NMBINS, logMbins);
	THStack *TrackerDimuon_MC_Stack = new THStack("TrackerDimuon mass: MC - stack plot", "TrackerDimuon mass: MC - stack plot");
	TH1F *TrackerDimuon_DATA = new TH1F("TrackerDimuon mass: DATA", "TrackerDimuon mass: DATA", NMBINS, logMbins);

	TH2F* DimuonTuneP_vs_Picky = new TH2F("DimuonMass: TuneP vs Picky", "DimuonMass: TuneP vs Picky", 44, 60, 500, 44, 60, 500);
	TH2F* DimuonTuneP_vs_DYT = new TH2F("DimuonMass: TuneP vs DYT", "DimuonMass: TuneP vs DYT", 44, 60, 500, 44, 60, 500);
	TH2F* DimuonTuneP_vs_Tracker = new TH2F("DimuonMass: TuneP vs Tracker", "DimuonMass: TuneP vs Tracker", 44, 60, 500, 44, 60, 500);

		
	THStack* pt_GEN_Stack = new THStack("p_{T} gen", "p_{T} gen");
	THStack* eta_GEN_till400_Stack = new THStack("#eta gen: pT_{reco} < 600 GeV", "#eta gen: pT_{reco} < 600 GeV");
	THStack* phi_GEN_till400_Stack = new THStack("#phi gen: pT_{reco} < 600 GeV", "#phi gen: pT_{reco} < 600 GeV");		
	THStack* eta_GEN_400to600_Stack = new THStack("#eta gen: 400 < pT_{reco} < 600 GeV", "#eta gen: 400 < pT_{reco} < 600 GeV");
	THStack* phi_GEN_400to600_Stack = new THStack("#phi gen: 400 < pT_{reco} < 600 GeV", "#phi gen: 400 < pT_{reco} < 600 GeV");		
	THStack* eta_GEN_above600_Stack = new THStack("#eta gen: pT_{reco} > 600 GeV", "#eta gen: pT_{reco} > 600 GeV");
	THStack* phi_GEN_above600_Stack = new THStack("#phi gen: pT_{reco} > 600 GeV", "#phi gen: pT_{reco} > 600 GeV");		

	for(int i = 0; i < 8; i++){

		pt_MC[i] = new TH1F(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i], binnum_pt, PT_BINS);
		pt_MC_Stack[i] = new THStack(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i]);
		pt_DATA[i] = new TH1F(pt_name[1] + reconstruction[i], pt_name[1] + reconstruction[i], binnum_pt, PT_BINS);
		eta_MC[i] = new TH1F(eta_name[0] + reconstruction[i], eta_name[0] + reconstruction[i], binnum_eta, ETA_BINS);
		eta_MC_Stack[i] = new THStack(eta_name[0] + reconstruction[i], eta_name[0] + reconstruction[i]);
		eta_DATA[i] = new TH1F(eta_name[1] + reconstruction[i], eta_name[1] + reconstruction[i],binnum_eta, ETA_BINS);
		phi_MC[i] = new TH1F(phi_name[0] + reconstruction[i], phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS);
		phi_MC_Stack[i] = new THStack(phi_name[0] + reconstruction[i], phi_name[0] + reconstruction[i]);
		phi_DATA[i] = new TH1F(phi_name[1] + reconstruction[i], phi_name[1] + reconstruction[i],binnum_phi, PHI_BINS);
		pt_vs_eta_MC[i] = new TH2F(pt_vs_eta_name[0] + reconstruction[i], pt_vs_eta_name[0] + reconstruction[i], binnum_eta, ETA_BINS, binnum_pt, PT_BINS);
		pt_vs_eta_DATA[i] = new TH2F(pt_vs_eta_name[1] + reconstruction[i], pt_vs_eta_name[1] + reconstruction[i], binnum_eta, ETA_BINS, binnum_pt, PT_BINS);
		pt_vs_phi_MC[i] = new TH2F(pt_vs_phi_name[0] + reconstruction[i], pt_vs_phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS, binnum_pt, PT_BINS);
		pt_vs_phi_DATA[i] = new TH2F(pt_vs_phi_name[1] + reconstruction[i], pt_vs_phi_name[1] + reconstruction[i], binnum_phi, PHI_BINS, binnum_pt, PT_BINS);
		eta_vs_phi_MC[i] = new TH2F(eta_vs_phi_name[0] + reconstruction[i], eta_vs_phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS, binnum_eta, ETA_BINS);
		eta_vs_phi_DATA[i] = new TH2F(eta_vs_phi_name[1] + reconstruction[i], eta_vs_phi_name[1] + reconstruction[i], binnum_phi, PHI_BINS, binnum_eta, ETA_BINS);

		eta_MC_till400[i] = new TH1F(eta_name[0] + reconstruction[i] + "_till400", eta_name[0] + reconstruction[i] + "_till400", binnum_eta, ETA_BINS);
		eta_MC_till400_Stack[i] = new THStack(eta_name[0] + reconstruction[i] + "_till400", eta_name[0] + reconstruction[i] + "_till400");
		eta_DATA_till400[i] = new TH1F(eta_name[1] + reconstruction[i] + "_till400", eta_name[1] + reconstruction[i] + "_till400",binnum_eta, ETA_BINS);
		eta_MC_400to600[i] = new TH1F(eta_name[0] + reconstruction[i] + "_400to600", eta_name[0] + reconstruction[i] + "_400to600", binnum_eta, ETA_BINS);
		eta_MC_400to600_Stack[i] = new THStack(eta_name[0] + reconstruction[i] + "_400to600", eta_name[0] + reconstruction[i] + "_400to600");
		eta_DATA_400to600[i] = new TH1F(eta_name[1] + reconstruction[i] + "_400to600", eta_name[1] + reconstruction[i] + "_400to600", binnum_eta, ETA_BINS);
		eta_MC_above600[i] = new TH1F(eta_name[0] + reconstruction[i] + "_above600", eta_name[0] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);
		eta_MC_above600_Stack[i] = new THStack(eta_name[0] + reconstruction[i] + "_above600", eta_name[0] + reconstruction[i] + "_above600");
		eta_DATA_above600[i] = new TH1F(eta_name[1] + reconstruction[i] + "_above600", eta_name[1] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);

		phi_MC_till400[i] = new TH1F(phi_name[0] + reconstruction[i] + "_till400", phi_name[0] + reconstruction[i] + "_till400", binnum_phi, PHI_BINS);
		phi_MC_till400_Stack[i] = new THStack(phi_name[0] + reconstruction[i] + "_till400", phi_name[0] + reconstruction[i] + "_till400");
		phi_DATA_till400[i] = new TH1F(phi_name[1] + reconstruction[i] + "_till400", phi_name[1] + reconstruction[i] + "_till400",binnum_phi, PHI_BINS);
		phi_MC_400to600[i] = new TH1F(phi_name[0] + reconstruction[i] + "_400to600", phi_name[0] + reconstruction[i] + "_400to600", binnum_phi, PHI_BINS);
		phi_MC_400to600_Stack[i] = new THStack(phi_name[0] + reconstruction[i] + "_400to600", phi_name[0] + reconstruction[i] + "_400to600");
		phi_DATA_400to600[i] = new TH1F(phi_name[1] + reconstruction[i] + "_400to600", phi_name[1] + reconstruction[i] + "_400to600", binnum_phi, PHI_BINS);
		phi_MC_above600[i] = new TH1F(phi_name[0] + reconstruction[i] + "_above600", phi_name[0] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);
		phi_MC_above600_Stack[i] = new THStack(phi_name[0] + reconstruction[i] + "_above600", phi_name[0] + reconstruction[i] + "_above600");
		phi_DATA_above600[i] = new TH1F(phi_name[1] + reconstruction[i] + "_above600", phi_name[1] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);

		pt_vs_met_MC[i] = new TH2F(pt_vs_met_name[0] + reconstruction[i], pt_vs_met_name[0] + reconstruction[i], binnum_met, MET_BINS, binnum_pt, PT_BINS);
		pt_vs_met_DATA[i] = new TH2F(pt_vs_met_name[1] + reconstruction[i], pt_vs_met_name[1] + reconstruction[i], binnum_met, MET_BINS, binnum_pt, PT_BINS);
		
		DileptonPt_MC_Stack[i] = new THStack(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i]);
		DileptonPt_DATA[i] = new TH1F(pt_name[1] + reconstruction[i], pt_name[1] + reconstruction[i], binnum_pt, PT_BINS);

		pt_MC[i]->GetYaxis()->SetTitle("Entries");
		pt_DATA[i]->GetYaxis()->SetTitle("Entries");
		eta_MC[i]->GetYaxis()->SetTitle("Entries");
		eta_DATA[i]->GetYaxis()->SetTitle("Entries");		
		phi_MC[i]->GetYaxis()->SetTitle("Entries");
		phi_DATA[i]->GetYaxis()->SetTitle("Entries");		

		pt_MC[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		pt_DATA[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		eta_MC[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA[i]->GetXaxis()->SetTitle("#eta");
		phi_MC[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA[i]->GetXaxis()->SetTitle("#phi");

		eta_MC_till400[i]->GetXaxis()->SetTitle("#eta");
		eta_MC_400to600[i]->GetXaxis()->SetTitle("#eta");
		eta_MC_above600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA_till400[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA_400to600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA_above600[i]->GetXaxis()->SetTitle("#eta");

		phi_MC_till400[i]->GetXaxis()->SetTitle("#phi");
		phi_MC_400to600[i]->GetXaxis()->SetTitle("#phi");
		phi_MC_above600[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA_till400[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA_400to600[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA_above600[i]->GetXaxis()->SetTitle("#phi");
		
		pt_vs_met_MC[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		pt_vs_met_MC[i]->GetXaxis()->SetTitle("met [GeV]");
		pt_vs_met_DATA[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		pt_vs_met_DATA[i]->GetXaxis()->SetTitle("met [GeV]");

	}


	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
// 			if(i == j) continue;
			TString vv = "_vs_";
			TString name_MC = pt_name[0] + reconstruction[i+1].Data() + vv + reconstruction[j+1].Data();
			TString name_DATA = pt_name[1] + reconstruction[i+1].Data() + vv + reconstruction[j+1].Data();
// 			std::cout<<i<<j<<"\t"<<name<<std::endl;
			Double_Count_MC[i][j] = new TH1F(name_MC, name_MC, binnum_pt, PT_BINS);
			Double_Count_MC[i][j]->GetXaxis()->SetTitle("p_{T}");
			Double_Count_DATA[i][j] = new TH1F(name_DATA, name_DATA, binnum_pt, PT_BINS);
			Double_Count_DATA[i][j]->GetXaxis()->SetTitle("p_{T}");
		}
	}			
	
	ofstream OutFile;
	OutFile.open ("prova.txt");
	
	for(int j = 14; j < 45; j++){

      	if(j == 44) continue; 
      	if(j > 29 and j < 39) continue;
//       	if(j > 38) continue;

        std::cout<<"opening.. "<<samples[j]<<" --- "<<events[j]<<std::endl;

      	TChain *treeMC = new TChain("SimpleNtupler/t");

     	treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_Pt_Assignment/ana_datamc_" + samples[j] + ".root");     	

	    Long64_t nentries = treeMC->GetEntries();
    	treeMC->SetBranchAddress("genWeight",&genWeight);

//     	float gen_Weight_tot = 0;
//     	for(int p=0; p<nentries; p++){
//     	 	treeMC->GetEntry(p);
//     		gen_Weight_tot += genWeight;
//     	}
    	
//     	std::cout<<"   "<<weight[j]<<" --- "<<gen_Weight_tot<<std::endl;
// 		weight[j] = (float)weight[j] / (float)gen_Weight_tot;
// 		float weight_orig = weight[j] * gen_Weight_tot / events[j];
// 		std::cout<<"      tot = "<<weight[j]<<" ----------- "<<weight_orig<<std::endl;     	
//      	continue;


     	TH1F *gen_Weight_clear = new TH1F("gen Weight", "gen Weight", 3, -1.5, 1.5);
		
		TH1F *GenDimuon_MC_clear = new TH1F("GenDimuon mass: MC", "GenDimuon mass: MC", NMBINS, logMbins);
		TH1F *Dimuon_MC_clear = new TH1F("Dimuon_mass_MC", "Dimuon_mass_MC", NMBINS, logMbins);
		TH1F *PickyDimuon_MC_clear = new TH1F("PickyDimuon mass: MC", "PickyDimuon mass: MC", NMBINS, logMbins);
		TH1F *DYTDimuon_MC_clear = new TH1F("DYTDimuon mass: MC", "DYTDimuon mass: MC", NMBINS, logMbins);
		TH1F *TrackerDimuon_MC_clear = new TH1F("TrackerDimuon mass: MC", "TrackerDimuon mass: MC", NMBINS, logMbins);
		TH1F* pt_GEN_clear = new TH1F("p_{T} gen", "p_{T} gen", binnum_pt, PT_BINS);
		TH1F* eta_GEN_till400_clear = new TH1F("#eta gen: pT_{reco} < 600 GeV", "#eta gen: pT_{reco} < 600 GeV", binnum_eta, ETA_BINS);
		TH1F* eta_GEN_400to600_clear = new TH1F("#eta gen: 400 < pT_{reco} < 600 GeV", "#eta gen: 400 < pT_{reco} < 600 GeV", binnum_eta, ETA_BINS);
		TH1F* eta_GEN_above600_clear = new TH1F("#eta gen: pT_{reco} > 600 GeV", "#eta gen: pT_{reco} > 600 GeV", binnum_eta, ETA_BINS);
		TH1F* phi_GEN_till400_clear = new TH1F("#phi gen: pT_{reco} < 600 GeV", "#phi gen: pT_{reco} < 600 GeV", binnum_phi, PHI_BINS);
		TH1F* phi_GEN_400to600_clear = new TH1F("#phi gen: 400 < pT_{reco} < 600 GeV", "#phi ge: 400 < pT_{reco} < 600 GeV", binnum_phi, PHI_BINS);
		TH1F* phi_GEN_above600_clear = new TH1F("#phi gen: pT_{reco} > 600 GeV", "#phi gen: pT_{reco} > 600 GeV", binnum_phi, PHI_BINS);

		for(int i = 0; i < 8; i++){
			pt_MC_clear[i] = new TH1F(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i], binnum_pt, PT_BINS);
			eta_MC_clear[i] = new TH1F(eta_name[0] + reconstruction[i], eta_name[0] + reconstruction[i], binnum_eta, ETA_BINS);
			phi_MC_clear[i] = new TH1F(phi_name[0] + reconstruction[i], phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS);
// 			pt_MC_clear[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
// 			eta_MC_clear[i]->GetXaxis()->SetTitle("#eta");
// 			phi_MC_clear[i]->GetXaxis()->SetTitle("#phi");
			eta_MC_till400_clear[i] = new TH1F(eta_name[0] + reconstruction[i] + "_till400", eta_name[0] + reconstruction[i] + "_till400", binnum_eta, ETA_BINS);
			eta_MC_400to600_clear[i] = new TH1F(eta_name[0] + reconstruction[i] + "_400to600", eta_name[0] + reconstruction[i] + "_400to600", binnum_eta, ETA_BINS);
			eta_MC_above600_clear[i] = new TH1F(eta_name[0] + reconstruction[i] + "_above600", eta_name[0] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);
			phi_MC_till400_clear[i] = new TH1F(phi_name[0] + reconstruction[i] + "_till400", phi_name[0] + reconstruction[i] + "_till400", binnum_phi, PHI_BINS);
			phi_MC_400to600_clear[i] = new TH1F(phi_name[0] + reconstruction[i] + "_400to600", phi_name[0] + reconstruction[i] + "_400to600", binnum_phi, PHI_BINS);
			phi_MC_above600_clear[i] = new TH1F(phi_name[0] + reconstruction[i] + "_above600", phi_name[0] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);

			DileptonPt_MC_clear[i] = new TH1F(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i], binnum_pt, PT_BINS);
		}
     	    	
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
//      treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeMC->SetBranchAddress("nvertices",&nvertices);
	     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);

	     treeMC->SetBranchAddress("lep_pt",lep_pt);
    	 treeMC->SetBranchAddress("lep_id",lep_id);
	     treeMC->SetBranchAddress("lep_eta",lep_eta);
	     treeMC->SetBranchAddress("lep_phi",lep_phi);
	     treeMC->SetBranchAddress("lep_dB",lep_dB);
//      treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	     treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    	 treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	     treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
	     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
// 	     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
//     	 treeMC->SetBranchAddress("lep_picky_numberOfValidMuonHits",lep_picky_numberOfValidMuonHits);
// 	     treeMC->SetBranchAddress("lep_dyt_numberOfValidMuonHits",lep_dyt_numberOfValidMuonHits);
//     	 treeMC->SetBranchAddress("lep_tpfms_numberOfValidMuonHits",lep_tpfms_numberOfValidMuonHits);
// 	     treeMC->SetBranchAddress("lep_stanAlone_numberOfValidMuonHits",lep_stanAlone_numberOfValidMuonHits);

//      treeMC->SetBranchAddress("lep_stanAlone_numberOfBadHits",lep_stanAlone_numberOfBadHits);
//      treeMC->SetBranchAddress("lep_stanAlone_numberOfMuonHits",lep_stanAlone_numberOfMuonHits);
	     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    	 treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
	     treeMC->SetBranchAddress("lep_pfIso",lep_pfIso);
	     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
    	 treeMC->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    	 treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
	     treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
    	 treeMC->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
	     treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
    	 treeMC->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
    	 treeMC->SetBranchAddress("lep_std_pt",lep_std_pt);
    	 treeMC->SetBranchAddress("lep_cocktail_pt",lep_cocktail_pt);

	     treeMC->SetBranchAddress("lep_glb_eta",lep_glb_eta);
    	 treeMC->SetBranchAddress("lep_dyt_eta",lep_dyt_eta); 
	     treeMC->SetBranchAddress("lep_picky_eta",lep_picky_eta);
    	 treeMC->SetBranchAddress("lep_tpfms_eta",lep_tpfms_eta);
    	 treeMC->SetBranchAddress("lep_std_eta",lep_std_eta);
    	 treeMC->SetBranchAddress("lep_tk_eta",lep_tk_eta);
    	 treeMC->SetBranchAddress("lep_tuneP_eta",lep_tuneP_eta);
    	 
	     treeMC->SetBranchAddress("lep_glb_phi",lep_glb_phi);
    	 treeMC->SetBranchAddress("lep_dyt_phi",lep_dyt_phi); 
	     treeMC->SetBranchAddress("lep_picky_phi",lep_picky_phi);
    	 treeMC->SetBranchAddress("lep_tpfms_phi",lep_tpfms_phi);
    	 treeMC->SetBranchAddress("lep_std_phi",lep_std_phi);
    	 treeMC->SetBranchAddress("lep_tk_phi",lep_tk_phi);
    	 treeMC->SetBranchAddress("lep_tuneP_phi",lep_tuneP_phi);

     //treeMC->SetBranchAddress("gen_lep_pt",gen_lep_pt); 
	     treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);   
	     treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2); 	

	     treeMC->SetBranchAddress("met_pt",&met_pt); 	
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples[j].Data(),j , weight[j], nentries);


    	for(int p=0; p<nentries; p++){
//      	 for(int p=0; p<10000; p++){
//      	 for(int p=0; p<1; p++){

    	 	if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
    	 	
    	 	treeMC->GetEntry(p);
    	 	
// 			if(j > 30) std::cout<<genWeight;
			
	     	if (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
	     		lep_pt[0]>53. && lep_pt[1]>53. && 
// 	     		lep_pt[0]>200. && lep_pt[1]>200. && 
    	 		lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
     			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
     			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
	     		lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
    	 		lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
    	 		(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 

	     		lep_sumPt[0]/lep_tk_pt[0]<0.10 && 
    	 		lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
/*	     		lep_sumPt[0]/lep_tk_pt[0]<0.05 &&  // Rel. trk iso < 0.05
    	 		lep_sumPt[1]/lep_tk_pt[1]<0.05 &&  // Rel. trk iso < 0.05
	     		lep_sumPt[0] < 50 &&  // Abs.trk iso < 5
    	 		lep_sumPt[1] < 50 &&  // Abs.trk iso < 5
	     		lep_pfIso[0]/lep_tk_pt[0]<0.05 &&  // Rel. trk iso < 0.05
    	 		lep_pfIso[1]/lep_tk_pt[1]<0.05 &&  // Rel. trk iso < 0.05
	     		lep_pfIso[0] < 150 &&  // Abs. PF iso < 150
    	 		lep_pfIso[1] < 150 &&  // Abs. PF iso < 150
*/
     			cos_angle>-0.9998 && 
     			lep_id[0]*lep_id[1]<0

				 /////// CON E SENZA VTX CUT ///////				
// 				&& vertex_chi2 < 20
				 /////// CON E SENZA VTX CUT ///////
				) { // this corresponds to the selection of an event: we want two muons, isolated that are at least tracker and global, with good quality tracks, opposite sign, close to beam spot. I let you add the trigger requirement 

					if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
					prev_event = event;
					
					if(dil_mass < 60 || dil_mass > 120) continue;
// 					if(dil_mass < 60 || dil_mass > 400) continue;
// 					if(dil_mass > 800) continue;

    			 	gen_Weight_clear->Fill(genWeight);
    		 	
					weight[j] *= genWeight;
//     				std::cout<<genWeight<<" "<<weight[j]<<std::endl;

     				if (
   	 					(lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 &&(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2))
   	 					 && (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 &&(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))
   						 && lep_pt_err[1]/lep_pt[1]<0.3
   						 && lep_pt_err[0]/lep_pt[0]<0.3
   						 && lep_glb_numberOfValidMuonHits[1] > 0
   						 && lep_glb_numberOfValidMuonHits[0] > 0
   						 && GoodVtx
   						 && vertex_chi2 < 20
     					){  
     						Dimuon_MC_clear->Fill(dil_mass,  weight[j]); 
							GenDimuon_MC_clear->Fill(dil_mass,  weight[j]); 
// 							std::cout<<gen_dil_mass<<" ---- TuneP = "<<dil_mass<<"("<<lep_pt[0]<<", "<<lep_pt[1]<<")";
/*     						lep_1.SetPtEtaPhiM(lep_picky_pt[0], lep_picky_eta[0], lep_picky_phi[0], 0.105);
     						lep_2.SetPtEtaPhiM(lep_picky_pt[1], lep_picky_eta[1], lep_picky_phi[1], 0.105);
     						ZPrime = lep_1 + lep_2;
//      						PickyDimuon_MC->Fill(ZPrime.M(),  weight[j]); 
     						if(lep_picky_pt[0] > -9 || lep_picky_pt[1] > -9){
	     						PickyDimuon_MC_clear->Fill(ZPrime.M(),  weight[j]); 
								DimuonTuneP_vs_Picky->Fill(ZPrime.M(), dil_mass);
							}
							
// 							std::cout<<" --- Picky = "<<ZPrime.M()<<"("<<lep_picky_pt[0]<<", "<<lep_picky_pt[1]<<")";

     						lep_1.SetPtEtaPhiM(lep_dyt_pt[0], lep_dyt_eta[0], lep_dyt_phi[0], 0.105);
     						lep_2.SetPtEtaPhiM(lep_dyt_pt[1], lep_dyt_eta[1], lep_dyt_phi[1], 0.105);
     						ZPrime = lep_1 + lep_2;
//      						DYTDimuon_MC->Fill(ZPrime.M(),  weight[j]); 
     						if(lep_dyt_pt[0] > -9 || lep_dyt_pt[1] > -9){
	     						DYTDimuon_MC_clear->Fill(ZPrime.M(),  weight[j]); 
								DimuonTuneP_vs_DYT->Fill(ZPrime.M(), dil_mass);
							}
							
// 							std::cout<<" --- DYT = "<<ZPrime.M()<<"("<<lep_dyt_pt[0]<<", "<<lep_dyt_pt[1]<<")";

							lep_1.SetPtEtaPhiM(lep_tk_pt[0], lep_tk_eta[0], lep_tk_phi[0], 0.105);
							lep_2.SetPtEtaPhiM(lep_tk_pt[1], lep_tk_eta[1], lep_tk_phi[1], 0.105);			
     						ZPrime = lep_1 + lep_2;
//      						TrackerDimuon_MC->Fill(ZPrime.M(),  weight[j]); 
     						if(lep_tk_pt[0] > -9 || lep_tk_pt[1] > -9){
	     						TrackerDimuon_MC_clear->Fill(ZPrime.M(),  weight[j]); 
    	 						DimuonTuneP_vs_Tracker->Fill(ZPrime.M(), dil_mass);
    	 					}
     						
// 							std::cout<<" --- Tracker = "<<ZPrime.M()<<"("<<lep_tk_pt[0]<<", "<<lep_tk_pt[1]<<")"<<std::endl;
*/
     					}

					int DeltaPhi[2] = {0};

					DileptonPt_MC_clear[7]->Fill(dil_pt, weight[j]);
					
					for(int h = 0; h < 2; h++){  // for on two muons in the event

// 	     				numberOfValidMuonHits[h][0] = lep_dyt_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][1] = lep_glb_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][2] = lep_picky_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][3] = lep_stanAlone_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][4] = lep_tpfms_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][5] = lep_TuneP_numberOfValidMuonHits[h];

	     				if (
    	 					(lep_numberOfMatchedStations[h] > 1 || (lep_numberOfMatchedStations[h] == 1 && !(lep_stationMask[h] == 1 || lep_stationMask[h] == 16)) || (lep_numberOfMatchedStations[h] == 1 &&(lep_stationMask[h] == 1 || lep_stationMask[h] == 16) && lep_numberOfMatchedRPCLayers[h] > 2))
     						 && lep_pt_err[h]/lep_pt[h]<0.3
     						 && lep_glb_numberOfValidMuonHits[h] > 0
     						 && lep_pt[h] > 200.
     					) {
   							
   							DeltaPhi[h] = 1;
   							
     						lepton_pt[h][0] = lep_pt[h];//lep_cocktail_pt[h];
     						lepton_pt[h][1] = lep_glb_pt[h];
     						lepton_pt[h][2] = lep_picky_pt[h];
    	 					lepton_pt[h][3] = lep_tpfms_pt[h];
	     					lepton_pt[h][4] = lep_dyt_pt[h];
     						lepton_pt[h][5] = lep_tk_pt[h];
     						lepton_pt[h][6] = lep_std_pt[h];
     						lepton_pt[h][7] = lep_tuneP_pt[h];
     						
     						lepton_eta[h][0] = lep_eta[h];
     						lepton_eta[h][1] = lep_glb_eta[h];
     						lepton_eta[h][2] = lep_picky_eta[h];
    	 					lepton_eta[h][3] = lep_tpfms_eta[h];
	     					lepton_eta[h][4] = lep_dyt_eta[h];
     						lepton_eta[h][5] = lep_tk_eta[h];
     						lepton_eta[h][6] = lep_std_eta[h];
     						lepton_eta[h][7] = lep_tuneP_eta[h]; 

     						lepton_phi[h][0] = lep_phi[h];
     						lepton_phi[h][1] = lep_glb_phi[h];
     						lepton_phi[h][2] = lep_picky_phi[h];
    	 					lepton_phi[h][3] = lep_tpfms_phi[h];
	     					lepton_phi[h][4] = lep_dyt_phi[h];
     						lepton_phi[h][5] = lep_tk_phi[h];
     						lepton_phi[h][6] = lep_std_phi[h];
     						lepton_phi[h][7] = lep_tuneP_phi[h];     						

//      						if(numberOfValidMuonHits[h][1] == 0) continue;

							int a = 1;
							int b = 2;
							int c = 3;
							int d = 4;
							int e = 5;

							////////////////////////////////////////////
							////////////////////////////////////////////		
							///////////// MULTIPLE COUNTING ////////////
/*							for(int f = 0; f < 5; f++){
								if(lepton_pt[h][a] == lepton_pt[h][7]){
									if(lepton_pt[h][a] == lepton_pt[h][b]){
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][14]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" ABCDE  "<<std::endl; 
												else pt_CountDouble_MC[a-1][10]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" ABCD  "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][11]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<"  ABCE "<<std::endl; 
												else pt_CountDouble_MC[a-1][4]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<"ABC   "<<std::endl; 
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][12]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" ABDE  "<<std::endl; 
												else pt_CountDouble_MC[a-1][5]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<"ABD   "<<std::endl; 
										}
										else if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][6]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" ABE  "<<std::endl; 
										else pt_CountDouble_MC[a-1][0]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" AB  "<<std::endl; 

									}
									else{
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][13]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" ACDE  "<<std::endl; 
												else pt_CountDouble_MC[a-1][7]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<"  ACD "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][8]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" ACE  "<<std::endl; 
												else pt_CountDouble_MC[a-1][1]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" AC  "<<std::endl;
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][9]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<"ADE"<<std::endl; 
											else pt_CountDouble_MC[a-1][2]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" AD  "<<std::endl; 
										}
										else 
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_MC[a-1][3]->Fill(lepton_pt[h][a], weight[j]); //std::cout<<" AE  "<<std::endl; 
									}
								}
// 						else std::cout<<" no tunep"<<std::endl;
								swap(a,b);
								swap(b,c);	
								swap(c,d);	
								swap(d,e);	
// 								std::cout<<"\t\t\t\t\t\t\t\t\t\t\t"<<a<<b<<c<<d<<e<<std::endl;	
							}
*/							///////////// MULTIPLE COUNTING ////////////
							////////////////////////////////////////////
							////////////////////////////////////////////
	     					
// 	     					for(int i = 1; i < 7; i++){ // for on different reconstruction
	     					for(int i = 1; i < 8; i++){ // for on different reconstruction
   	 						    	 						
    	 						if(lep_pt[h] == lepton_pt[h][i]){
    	 						
    	 							bool multiple_reconstruction = false;	
    	 							for(int z = 1; z < 7; z++){
    	 								if(z == i || i == 7) continue;
    	 								if(lepton_pt[h][z] == lepton_pt[h][i]){
    	 									multiple_reconstruction = true;
    	 									break;
    	 								}
    	 							}
    	 							if(multiple_reconstruction) break; //to remove multiple reconstructions
    	 							
//     	 							if(j == 30 && i == 1)
// 			     						std::cout<<"+++++++++++"<<lepton_pt[h][0]<<"\t"<<lepton_pt[h][1]<<"\t"<<lepton_pt[h][2]<<"\t"<<lepton_pt[h][3]
// 			     						<<"\t"<<lepton_pt[h][4]<<"\t"<<lepton_pt[h][5]<<"\t"<<lepton_pt[h][6]<<"\t"<<lepton_pt[h][7]<<std::endl;

									// Select Picky only if Picky = TuneP    	 						
//     	 							if(i == 2 && (lepton_pt[h][i] == lepton_pt[h][3] || lepton_pt[h][i] == lepton_pt[h][5])) continue;
									// Select DYT only if DYT = TuneP    	 						
// 									if(i == 4 && lepton_pt[h][i] == lepton_pt[h][5]) continue;
//      								if(i != 7) count_assignment_MC_till400[i][j-14]++;
	     							pt_MC[i]->Fill(lepton_pt[h][i], weight[j]);
	     							eta_MC[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     							phi_MC[i]->Fill(lepton_phi[h][i], weight[j]);
	     							pt_MC_clear[i]->Fill(lepton_pt[h][i]);//, weight[j]); //for stack plot
	     							eta_MC_clear[i]->Fill(fabs(lepton_eta[h][i]), weight[j]); //for stack plot
	     							phi_MC_clear[i]->Fill(lepton_phi[h][i], weight[j]); //for stack plot
    	 							pt_vs_eta_MC[i]->Fill(lep_eta[h], lep_pt[h], weight[j]);
    	 							pt_vs_phi_MC[i]->Fill(lep_phi[h], lep_pt[h], weight[j]);
    	 							eta_vs_phi_MC[i]->Fill(lep_phi[h], lep_eta[h], weight[j]);
    	 							pt_vs_met_MC[i]->Fill(met_pt, lep_pt[h], weight[j]);

		 	     					if(i != 7) pt_GEN_clear->Fill(gen_lep_pt[h]);//,  weight[j]);
//     	 							pt_MC[7]->Fill(lepton_pt[h][7], weight[j]);
// 	     							eta_MC[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
//      								phi_MC[7]->Fill(lepton_phi[h][7], weight[j]);
//      								pt_MC_clear[7]->Fill(lepton_pt[h][7], weight[j]); //for stack plot
//      								eta_MC_clear[7]->Fill(fabs(lepton_eta[h][7]), weight[j]); //for stack plot
//      								phi_MC_clear[7]->Fill(lepton_phi[h][7], weight[j]); //for stack plot
//    	 								pt_vs_eta_MC[7]->Fill(lep_eta[h], lep_pt[h], weight[j]);
//    	 								pt_vs_phi_MC[7]->Fill(lep_phi[h], lep_pt[h], weight[j]);
//    		 							eta_vs_phi_MC[7]->Fill(lep_phi[h], lep_eta[h], weight[j]);
// 	   	 							pt_vs_met_MC[7]->Fill(met_pt, lep_pt[h], weight[j]);

		     						if(lepton_pt[h][i] < 400){
// 										if(lep_phi[h] > 0 && lep_phi[h] < 0.8) Phi_bin->Fill(lep_phi[h], gen_lep_phi[h]);
	     								count_assignment_MC_till400[i][j-14]++;//=weight[j];
// 	    	 							count_assignment_MC_till400[7][j-14]+=weight[j];
		     							eta_MC_till400[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
		     							phi_MC_till400[i]->Fill(lepton_phi[h][i], weight[j]);
		     							eta_MC_till400_clear[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
		     							phi_MC_till400_clear[i]->Fill(lepton_phi[h][i], weight[j]);
// 		    	 						eta_MC_till400[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
// 	    	 							phi_MC_till400[7]->Fill(lepton_phi[h][7], weight[j]);
// 		     							eta_MC_till400_clear[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
// 		     							phi_MC_till400_clear[7]->Fill(lepton_phi[h][7], weight[j]);
			 	     					if(i != 7) eta_GEN_till400_clear->Fill(fabs(gen_lep_eta[h]),  weight[j]);
			 	     					if(i != 7) phi_GEN_till400_clear->Fill(gen_lep_phi[h],  weight[j]);
		     						}
		     						else if(lepton_pt[h][i] > 400 && lepton_pt[h][i] < 600){
		     						
// 		     							genW_vs_eta->Fill(genWeight, fabs(gen_lep_eta[h]));
// 		     							genW_vs_phi->Fill(genWeight, (gen_lep_phi[h]));
// 		     							genW_vs_pt->Fill(genWeight, (gen_lep_pt[h]));
		     						
// 										if(lep_phi[h] > 0 && lep_phi[h] < 0.8) Phi_bin_46->Fill(lep_phi[h], gen_lep_phi[h]);
	     								count_assignment_MC_400to600[i][j-14]++;//=weight[j];
// 	    	 							count_assignment_MC_400to600[7][j-14]+=weight[j];
	     								eta_MC_400to600[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     								phi_MC_400to600[i]->Fill(lepton_phi[h][i], weight[j]);
	     								eta_MC_400to600_clear[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     								phi_MC_400to600_clear[i]->Fill(lepton_phi[h][i], weight[j]);
// 	     								eta_MC_400to600[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
// 	    	 							phi_MC_400to600[7]->Fill(lepton_phi[h][7], weight[j]);
//     	 								eta_MC_400to600_clear[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
// 	     								phi_MC_400to600_clear[7]->Fill(lepton_phi[h][7], weight[j]);
				     					if(i != 7) eta_GEN_400to600_clear->Fill(fabs(gen_lep_eta[h]),  weight[j]);
				     					if(i != 7) phi_GEN_400to600_clear->Fill(gen_lep_phi[h],  weight[j]);

		     						}
	    	 						else{
// 										if(lep_phi[h] > 0 && lep_phi[h] < 0.8) Phi_bin_6->Fill(lep_phi[h], gen_lep_phi[h]);
	     								count_assignment_MC_above600[i][j-14]++;//=weight[j];
// 	    	 							count_assignment_MC_above600[7][j-14]+=weight[j];
	     								eta_MC_above600[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     								phi_MC_above600[i]->Fill(lepton_phi[h][i], weight[j]);
	     								eta_MC_above600_clear[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     								phi_MC_above600_clear[i]->Fill(lepton_phi[h][i], weight[j]);
// 	     								eta_MC_above600[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
// 	    	 							phi_MC_above600[7]->Fill(lepton_phi[h][7], weight[j]);
//     	 								eta_MC_above600_clear[7]->Fill(fabs(lepton_eta[h][7]), weight[j]);
// 	     								phi_MC_above600_clear[7]->Fill(lepton_phi[h][7], weight[j]);
				     					if(i != 7) eta_GEN_above600_clear->Fill(fabs(gen_lep_eta[h]),  weight[j]);
				     					if(i != 7) phi_GEN_above600_clear->Fill(gen_lep_phi[h],  weight[j]);
	     							}
// 	    	 							std::cout<<lepton_pt[h][i]<<"\t"<<count_assignment_MC_above600[7][j-14]<<"\t"<<count_assignment_MC_till400[7][j-14]<<std::endl;
//     	 							break; // to take only one reconstruction in case of double counting --- DYT chosen because it is the first
    	 						}
     						}
    /*	 					if(lep_pt[h] == lepton_pt[h][7]){
//     	 						count_assignment_MC_till400[7][j-14]++;
     							pt_MC[7]->Fill(lepton_pt[h][7], weight[j]);
     							eta_MC[7]->Fill(lepton_eta[h][7], weight[j]);
     							phi_MC[7]->Fill(lepton_phi[h][7], weight[j]);
     							pt_MC_clear[7]->Fill(lepton_pt[h][7], weight[j]); //for stack plot
     							eta_MC_clear[7]->Fill(lepton_eta[h][7], weight[j]); //for stack plot
     							phi_MC_clear[7]->Fill(lepton_phi[h][7], weight[j]); //for stack plot
   	 							pt_vs_eta_MC[7]->Fill(lep_eta[h], lep_pt[h], weight[j]);
   	 							pt_vs_phi_MC[7]->Fill(lep_phi[h], lep_pt[h], weight[j]);
   	 							eta_vs_phi_MC[7]->Fill(lep_phi[h], lep_eta[h], weight[j]);
   	 							pt_vs_met_MC[7]->Fill(met_pt, lep_pt[h], weight[j]);


		     					if(lepton_pt[h][7] < 600){
	    	 						count_assignment_MC_till400[7][j-14]++;
		     						eta_MC_till400[7]->Fill(lepton_eta[h][7], weight[j]);
	     							phi_MC_till400[7]->Fill(lepton_phi[h][7], weight[j]);
	     							eta_MC_till400_clear[7]->Fill(lepton_eta[h][7], weight[j]);
	     							phi_MC_till400_clear[7]->Fill(lepton_phi[h][7], weight[j]);

		     					}
	    	 					else{
	    	 						count_assignment_MC_above600[7][j-14]++;
	     							eta_MC_above600[7]->Fill(lepton_eta[h][7], weight[j]);
	     							phi_MC_above600[7]->Fill(lepton_phi[h][7], weight[j]);
     								eta_MC_above600_clear[7]->Fill(lepton_eta[h][7], weight[j]);
     								phi_MC_above600_clear[7]->Fill(lepton_phi[h][7], weight[j]);

	     						}
     						} */
//      						std::cout<<count_assignment_MC_till400[0][0]<<"\t"<<count_assignment_MC_till400[1][0]<<"\t"<<count_assignment_MC_till400[2][0]
//      						<<"\t"<<count_assignment_MC_till400[3][0]<<"\t"<<count_assignment_MC_till400[4][0]<<"\t"<<count_assignment_MC_till400[5][0]
//      						<<"\t"<<count_assignment_MC_till400[6][0]<<"\t"<<count_assignment_MC_till400[7][0]
//      						<<"\t\t\t"<<lepton_pt[h][0]<<"\t"<<lepton_pt[h][1]<<"\t"<<lepton_pt[h][2]<<"\t"<<lepton_pt[h][3]
//      						<<"\t"<<lepton_pt[h][4]<<"\t"<<lepton_pt[h][5]<<"\t"<<lepton_pt[h][6]<<"\t"<<lepton_pt[h][7]<<std::endl;

//      						for(int i =0; i < 8; i++)
//      							std::cout<<reconstruction[i].Data()<<"\t"<<count_assignment_MC_till400[i][j-30]<<std::endl;
//      						std::cout<<"-"<<std::endl;
     					}
     					else NotGoodQualityMuon_MC++;     					
     				} // for on muons
			
			if(DeltaPhi[0]*DeltaPhi[1] != 0){
				countDeltaPhi_MC++;
				DeltaPhi_MC_till400_clear->Fill(fabs(lepton_phi[0][7]-lepton_phi[1][7]), weight[j]);
				DeltaPhi_GEN_till400_clear->Fill(fabs(gen_lep_phi[0]-gen_lep_phi[1]), weight[j]);
			}
			else NOcountDeltaPhi_MC++;			

			weight[j] /= genWeight;
	       }//end of condition on event

		 }// end loop p
	 
	 	 Color_t color;
	 	 if(j == 44) color = kMagenta-10;
	 	 if(j == 43) color = kMagenta-8;
	 	 if(j == 42) color = kMagenta-6;
	 	 if(j == 41) color = kMagenta-2;
	 	 if(j == 40) color = kMagenta;
	 	 if(j == 39) color = kMagenta+3;
	 	 if(j > 30 && j < 39) color = 3;
	 	 if(j == 30) color = kGreen+3;
	 	 if(j > 27 && j < 30) color = kRed+3;
	 	 if(j > 25 && j < 28) color = 2;
	 	 if(j > 20 && j < 26) color = kRed-10;
	 	 if(j > 15 && j < 21) color = 4;
	 	 if(j == 14 || j == 15) color = kBlue+2;

		 for(int i = 0; i < 8; i++){

			pt_MC_clear[i]->SetFillColor(color);
			pt_MC_clear[i]->SetLineColor(color);
			eta_MC_clear[i]->SetFillColor(color);
			eta_MC_clear[i]->SetLineColor(color);
			phi_MC_clear[i]->SetFillColor(color);
			phi_MC_clear[i]->SetLineColor(color);
			eta_MC_till400_clear[i]->SetFillColor(color);
			eta_MC_till400_clear[i]->SetLineColor(color);
			eta_MC_400to600_clear[i]->SetFillColor(color);
			eta_MC_400to600_clear[i]->SetLineColor(color);
			eta_MC_above600_clear[i]->SetFillColor(color);
			eta_MC_above600_clear[i]->SetLineColor(color);
			phi_MC_till400_clear[i]->SetFillColor(color);
			phi_MC_till400_clear[i]->SetLineColor(color);
			phi_MC_400to600_clear[i]->SetFillColor(color);
			phi_MC_400to600_clear[i]->SetLineColor(color);
			phi_MC_above600_clear[i]->SetFillColor(color);
			phi_MC_above600_clear[i]->SetLineColor(color);
			pt_MC_Stack[i]->Add(pt_MC_clear[i], "HIST");
			eta_MC_Stack[i]->Add(eta_MC_clear[i], "HIST");
			phi_MC_Stack[i]->Add(phi_MC_clear[i], "HIST");
			eta_MC_till400_Stack[i]->Add(eta_MC_till400_clear[i], "HIST");
			phi_MC_till400_Stack[i]->Add(phi_MC_till400_clear[i], "HIST");
			eta_MC_400to600_Stack[i]->Add(eta_MC_400to600_clear[i], "HIST");
			phi_MC_400to600_Stack[i]->Add(phi_MC_400to600_clear[i], "HIST");
			eta_MC_above600_Stack[i]->Add(eta_MC_above600_clear[i], "HIST");
			phi_MC_above600_Stack[i]->Add(phi_MC_above600_clear[i], "HIST");
			
			DileptonPt_MC_clear[i]->SetFillColor(color);
			DileptonPt_MC_clear[i]->SetLineColor(color);
			DileptonPt_MC_Stack[i]->Add(DileptonPt_MC_clear[i], "HIST");
		}		 

		 GenDimuon_MC_clear->SetFillColor(color);
		 GenDimuon_MC_clear->SetLineColor(color);
		 Dimuon_MC_clear->SetFillColor(color);
		 Dimuon_MC_clear->SetLineColor(color);
		 PickyDimuon_MC_clear->SetFillColor(color);
		 PickyDimuon_MC_clear->SetLineColor(color);
		 DYTDimuon_MC_clear->SetFillColor(color);
		 DYTDimuon_MC_clear->SetLineColor(color);
		 TrackerDimuon_MC_clear->SetFillColor(color);
		 TrackerDimuon_MC_clear->SetLineColor(color);
		 pt_GEN_clear->SetFillColor(color);
		 pt_GEN_clear->SetLineColor(color);
		 eta_GEN_till400_clear->SetFillColor(color);
		 eta_GEN_till400_clear->SetLineColor(color);
		 phi_GEN_till400_clear->SetFillColor(color);
		 phi_GEN_till400_clear->SetLineColor(color);
		 eta_GEN_400to600_clear->SetFillColor(color);
		 eta_GEN_400to600_clear->SetLineColor(color);
		 phi_GEN_400to600_clear->SetFillColor(color);
		 phi_GEN_400to600_clear->SetLineColor(color);
		 eta_GEN_above600_clear->SetFillColor(color);
		 eta_GEN_above600_clear->SetLineColor(color);
		 phi_GEN_above600_clear->SetFillColor(color);
		 phi_GEN_above600_clear->SetLineColor(color);
		 
		 gen_Weight_clear->SetFillColor(color);
		 gen_Weight_clear->SetLineColor(color);
		 gen_Weight_Stack->Add(gen_Weight_clear, "HIST");
		 
		 		 
		 GenDimuon_MC_Stack->Add(GenDimuon_MC_clear, "HIST");
		 Dimuon_MC_Stack->Add(Dimuon_MC_clear, "HIST");
		 PickyDimuon_MC_Stack->Add(PickyDimuon_MC_clear, "HIST");
		 DYTDimuon_MC_Stack->Add(DYTDimuon_MC_clear, "HIST");
		 TrackerDimuon_MC_Stack->Add(TrackerDimuon_MC_clear, "HIST");
		 pt_GEN_Stack->Add(pt_GEN_clear, "HIST");
		 eta_GEN_till400_Stack->Add(eta_GEN_till400_clear, "HIST");
		 phi_GEN_till400_Stack->Add(phi_GEN_till400_clear, "HIST");
		 eta_GEN_400to600_Stack->Add(eta_GEN_400to600_clear, "HIST");
		 phi_GEN_400to600_Stack->Add(phi_GEN_400to600_clear, "HIST");
		 eta_GEN_above600_Stack->Add(eta_GEN_above600_clear, "HIST");
		 phi_GEN_above600_Stack->Add(phi_GEN_above600_clear, "HIST");

		 DeltaPhi_MC_till400_clear->SetFillColor(color);
		 DeltaPhi_MC_till400_clear->SetLineColor(color);		 
		 DeltaPhi_MC_till400_Stack->Add(DeltaPhi_MC_till400_clear, "HIST");
		 DeltaPhi_GEN_till400_clear->SetFillColor(color);
		 DeltaPhi_GEN_till400_clear->SetLineColor(color);		 
		 DeltaPhi_GEN_till400_Stack->Add(DeltaPhi_GEN_till400_clear, "HIST");

	 
     }// end loop on MC
    OutFile.close();
    
////////      DATA    ///////

	for(int i = 0; i < 7; i++)
		for(int j = 0; j < 2; j++)
			lepton_pt[j][i] = {0};
		
    TChain *treeDATA = new TChain("SimpleNtupler/t");
    treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/DATA_Pt_Assignment/ana_datamc_data_Pt_Assignment.root");

     treeDATA->SetBranchAddress("event",&event);
     treeDATA->SetBranchAddress("run",&run);
//      treeDATA->SetBranchAddress("lumi",&lumi);
     	
     treeDATA->SetBranchAddress("dil_mass",&dil_mass);
//      treeDATA->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
     treeDATA->SetBranchAddress("dil_pt",&dil_pt);
	treeDATA->SetBranchAddress("cos_angle",&cos_angle);
//      treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeDATA->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeDATA->SetBranchAddress("nvertices",&nvertices);
    treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);
// 
	treeDATA->SetBranchAddress("lep_pt",lep_pt);
    treeDATA->SetBranchAddress("lep_id",lep_id);
	treeDATA->SetBranchAddress("lep_eta",lep_eta);
	treeDATA->SetBranchAddress("lep_phi",lep_phi);
	treeDATA->SetBranchAddress("lep_dB",lep_dB);
//      treeDATA->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	treeDATA->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    treeDATA->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	treeDATA->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
	treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
    treeDATA->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
    treeDATA->SetBranchAddress("lep_picky_numberOfValidMuonHits",lep_picky_numberOfValidMuonHits);
    treeDATA->SetBranchAddress("lep_dyt_numberOfValidMuonHits",lep_dyt_numberOfValidMuonHits);
    treeDATA->SetBranchAddress("lep_tpfms_numberOfValidMuonHits",lep_tpfms_numberOfValidMuonHits);
    treeDATA->SetBranchAddress("lep_stanAlone_numberOfValidMuonHits",lep_stanAlone_numberOfValidMuonHits);

//      treeDATA->SetBranchAddress("lep_stanAlone_numberOfBadHits",lep_stanAlone_numberOfBadHits);
//      treeDATA->SetBranchAddress("lep_stanAlone_numberOfMuonHits",lep_stanAlone_numberOfMuonHits);
	treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	treeDATA->SetBranchAddress("lep_sumPt",lep_sumPt);
	treeDATA->SetBranchAddress("lep_pfIso",lep_pfIso);
	treeDATA->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    treeDATA->SetBranchAddress("lep_tk_pt",lep_tk_pt);    
	treeDATA->SetBranchAddress("lep_glb_pt",lep_glb_pt);
    treeDATA->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
	treeDATA->SetBranchAddress("lep_picky_pt",lep_picky_pt);
    treeDATA->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
	treeDATA->SetBranchAddress("lep_std_pt",lep_std_pt);
	treeDATA->SetBranchAddress("lep_cocktail_pt",lep_cocktail_pt);
	
    treeDATA->SetBranchAddress("lep_glb_eta",lep_glb_eta);
   	treeDATA->SetBranchAddress("lep_dyt_eta",lep_dyt_eta); 
    treeDATA->SetBranchAddress("lep_picky_eta",lep_picky_eta);
   	treeDATA->SetBranchAddress("lep_tpfms_eta",lep_tpfms_eta);
   	treeDATA->SetBranchAddress("lep_std_eta",lep_std_eta);
   	treeDATA->SetBranchAddress("lep_tk_eta",lep_tk_eta);
  	treeDATA->SetBranchAddress("lep_tuneP_eta",lep_tuneP_eta);
  	
    treeDATA->SetBranchAddress("lep_glb_phi",lep_glb_phi);
   	treeDATA->SetBranchAddress("lep_dyt_phi",lep_dyt_phi); 
    treeDATA->SetBranchAddress("lep_picky_phi",lep_picky_phi);
   	treeDATA->SetBranchAddress("lep_tpfms_phi",lep_tpfms_phi);
   	treeDATA->SetBranchAddress("lep_std_phi",lep_std_phi);
   	treeDATA->SetBranchAddress("lep_tk_phi",lep_tk_phi);
  	treeDATA->SetBranchAddress("lep_tuneP_phi",lep_tuneP_phi);

	treeDATA->SetBranchAddress("lep_pt_err",lep_pt_err);
     //treeDATA->SetBranchAddress("gen_lep_pt",gen_lep_pt); 
	treeDATA->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);   
	treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2); 

	treeDATA->SetBranchAddress("met_pt",&met_pt); 		

	TH1F * hDeltaR = new TH1F("DeltaR", "DeltaR", 200, 0, 10);
	TH1F * pt_DATA_BeforeDeltaRcut = new TH1F("pt_DATA_BeforeDeltaRcut", "pt_DATA_BeforeDeltaRcut", binnum_pt, PT_BINS);
	Long64_t nentries = treeDATA->GetEntries();
	
	printf("opening... DATA --- %lld\n", nentries);    
		int bhu_data = 0;
	
	for(int p=0; p<nentries; p++){
// 	for(int p=0; p<10000; p++){
// 	for(int p=0; p<1; p++){

		if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
	
		treeDATA->GetEntry(p);

	    if (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
	     		lep_pt[0]>53. && lep_pt[1]>53. && 
// 	     		lep_pt[0]>200. && lep_pt[1]>200. && 
    	 		lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
     			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
     			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
	     		lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
    	 		lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
    	 		(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
    	 		
	     		lep_sumPt[0]/lep_tk_pt[0]<0.10 && 
    	 		lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
/*	     		lep_sumPt[0]/lep_tk_pt[0]<0.05 &&  // Rel. trk iso < 0.05
    	 		lep_sumPt[1]/lep_tk_pt[1]<0.05 &&  // Rel. trk iso < 0.05
	     		lep_sumPt[0] < 50 &&  // Abs.trk iso < 5
    	 		lep_sumPt[1] < 50 &&  // Abs.trk iso < 5
	     		lep_pfIso[0]/lep_tk_pt[0]<0.05 &&  // Rel. trk iso < 0.05
    	 		lep_pfIso[1]/lep_tk_pt[1]<0.05 &&  // Rel. trk iso < 0.05
	     		lep_pfIso[0] < 150 &&  // Abs. PF iso < 150
    	 		lep_pfIso[1] < 150 &&  // Abs. PF iso < 150
*/
     			cos_angle>-0.9998 && 
     			lep_id[0]*lep_id[1]<0

				 /////// CON E SENZA VTX CUT ///////				
// 				&& vertex_chi2 < 20
				 /////// CON E SENZA VTX CUT ///////
				) { // this corresponds to the selection of an event: we want two muons, isolated that are at least tracker and global, with good quality tracks, opposite sign, close to beam spot. I let you add the trigger requirement 

					if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt	
					prev_event = event;
					
					if(dil_mass < 60 || dil_mass > 120) continue;	
									
     				if (
   	 					(lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 &&(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2))
   	 					 && (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 &&(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))
   						 && lep_pt_err[1]/lep_pt[1]<0.3
   						 && lep_pt_err[0]/lep_pt[0]<0.3
   						 && lep_glb_numberOfValidMuonHits[1] > 0
   						 && lep_glb_numberOfValidMuonHits[0] > 0
   						 && GoodVtx
   						 && vertex_chi2 < 20
     					) Dimuon_DATA->Fill(dil_mass);

					int DeltaPhi[2] = {0};
					
					DileptonPt_DATA[7]->Fill(dil_pt);
					
					float deltaR = -999;
					float deltaPhi = -999;
					TLorentzVector mu1, mu2;
					mu1.SetPtEtaPhiM(lep_pt[0], lep_eta[0], lep_phi[0], 0.105);
					mu2.SetPtEtaPhiM(lep_pt[1], lep_eta[1], lep_phi[1], 0.105);
					deltaR = mu1.DeltaR(mu2);
					deltaPhi = mu1.DeltaPhi(mu2);
					hDeltaR->Fill(deltaR);
					DeltaPhi_DATA_till400->Fill(deltaPhi);

					for(int h = 0; h < 2; h++){ // for on two muons in the event

// 	     				numberOfValidMuonHits[h][0] = lep_dyt_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][1] = lep_glb_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][2] = lep_picky_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][3] = lep_stanAlone_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][4] = lep_tpfms_numberOfValidMuonHits[h];
// 	     				numberOfValidMuonHits[h][5] = lep_TuneP_numberOfValidMuonHits[h];

	     				if (
    	 					(lep_numberOfMatchedStations[h] > 1 || (lep_numberOfMatchedStations[h] == 1 && !(lep_stationMask[h] == 1 || lep_stationMask[h] == 16)) || (lep_numberOfMatchedStations[h] == 1 &&(lep_stationMask[h] == 1 || lep_stationMask[h] == 16) && lep_numberOfMatchedRPCLayers[h] > 2))
     						 && lep_pt_err[h]/lep_pt[h]<0.3
     						 && lep_glb_numberOfValidMuonHits[h] > 0
     						 && lep_pt[h] > 200
     					) {   	
							
							pt_DATA_BeforeDeltaRcut->Fill(lep_pt[h]);
														
							if(deltaR < 0.1) continue;

							DeltaPhi[h] = 1;     					
// 				     		if(lep_pt[0] < 200. || lep_pt[1] < 200.) continue;
							
							lepton_pt[h][0] = lep_pt[h];//lep_cocktail_pt[h];
     						lepton_pt[h][1] = lep_glb_pt[h];
     						lepton_pt[h][2] = lep_picky_pt[h];
    	 					lepton_pt[h][3] = lep_tpfms_pt[h];
	     					lepton_pt[h][4] = lep_dyt_pt[h];
     						lepton_pt[h][5] = lep_tk_pt[h];
     						lepton_pt[h][6] = lep_std_pt[h];
     						lepton_pt[h][7] = lep_tuneP_pt[h];

     						lepton_eta[h][0] = lep_eta[h];
     						lepton_eta[h][1] = lep_glb_eta[h];
     						lepton_eta[h][2] = lep_picky_eta[h];
    	 					lepton_eta[h][3] = lep_tpfms_eta[h];
	     					lepton_eta[h][4] = lep_dyt_eta[h];
     						lepton_eta[h][5] = lep_tk_eta[h];
     						lepton_eta[h][6] = lep_std_eta[h];
     						lepton_eta[h][7] = lep_tuneP_eta[h];  
     						
     						lepton_phi[h][0] = lep_phi[h];
     						lepton_phi[h][1] = lep_glb_phi[h];
     						lepton_phi[h][2] = lep_picky_phi[h];
    	 					lepton_phi[h][3] = lep_tpfms_phi[h];
	     					lepton_phi[h][4] = lep_dyt_phi[h];
     						lepton_phi[h][5] = lep_tk_phi[h];
     						lepton_phi[h][6] = lep_std_phi[h];
     						lepton_phi[h][7] = lep_tuneP_phi[h];    						

//      						if(numberOfValidMuonHits[h][1] == 0) continue;

							int a = 1;
							int b = 2;
							int c = 3;
							int d = 4;
							int e = 5;

							////////////////////////////////////////////
							////////////////////////////////////////////		
							///////////// MULTIPLE COUNTING ////////////
/*							for(int f = 0; f < 5; f++){
								if(lepton_pt[h][a] == lepton_pt[h][7]){
									if(lepton_pt[h][a] == lepton_pt[h][b]){
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][14]->Fill(lepton_pt[h][a]); //std::cout<<" ABCDE  "<<std::endl; 
												else pt_CountDouble_DATA[a-1][10]->Fill(lepton_pt[h][a]); //std::cout<<" ABCD  "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][11]->Fill(lepton_pt[h][a]); //std::cout<<"  ABCE "<<std::endl; 
												else pt_CountDouble_DATA[a-1][4]->Fill(lepton_pt[h][a]); //std::cout<<"ABC   "<<std::endl; 
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][12]->Fill(lepton_pt[h][a]); //std::cout<<" ABDE  "<<std::endl; 
												else pt_CountDouble_DATA[a-1][5]->Fill(lepton_pt[h][a]); //std::cout<<"ABD   "<<std::endl; 
										}
										else if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][6]->Fill(lepton_pt[h][a]); //std::cout<<" ABE  "<<std::endl; 
										else pt_CountDouble_DATA[a-1][0]->Fill(lepton_pt[h][a]); //std::cout<<" AB  "<<std::endl; 

									}
									else{
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][13]->Fill(lepton_pt[h][a]); //std::cout<<" ACDE  "<<std::endl; 
												else pt_CountDouble_DATA[a-1][7]->Fill(lepton_pt[h][a]); //std::cout<<"  ACD "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][8]->Fill(lepton_pt[h][a]); //std::cout<<" ACE  "<<std::endl; 
												else pt_CountDouble_DATA[a-1][1]->Fill(lepton_pt[h][a]); //std::cout<<" AC  "<<std::endl;
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][9]->Fill(lepton_pt[h][a]); //std::cout<<"ADE"<<std::endl; 
											else pt_CountDouble_DATA[a-1][2]->Fill(lepton_pt[h][a]); //std::cout<<" AD  "<<std::endl; 
										}
										else 
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA[a-1][3]->Fill(lepton_pt[h][a]); //std::cout<<" AE  "<<std::endl; 
									}
								}
// 								else std::cout<<" no tunep"<<std::endl;
								swap(a,b);
								swap(b,c);	
								swap(c,d);	
								swap(d,e);	
							}
*/							///////////// MULTIPLE COUNTING ////////////
							////////////////////////////////////////////
							////////////////////////////////////////////

	     					for(int i = 1; i < 8; i++){ // for on different reconstruction
   	 						    	 						
    	 						if(lep_pt[h] == lepton_pt[h][i]){
    	 						
    	 							bool multiple_reconstruction = false;	
    	 							for(int z = 1; z < 7; z++){
    	 								if(z == i || i == 7) continue;
    	 								if(lepton_pt[h][z] == lepton_pt[h][i]){
    	 									multiple_reconstruction = true;
    	 									break;
    	 								}
    	 							}
    	 							if(multiple_reconstruction) break; //to remove multiple reconstructions
    	 							
//     	 							if(j == 30 && i == 1)
// 			     						std::cout<<"+++++++++++"<<lepton_pt[h][0]<<"\t"<<lepton_pt[h][1]<<"\t"<<lepton_pt[h][2]<<"\t"<<lepton_pt[h][3]
// 			     						<<"\t"<<lepton_pt[h][4]<<"\t"<<lepton_pt[h][5]<<"\t"<<lepton_pt[h][6]<<"\t"<<lepton_pt[h][7]<<std::endl;

									// Select Picky only if Picky = TuneP    	 						
//     	 							if(i == 2 && (lepton_pt[h][i] == lepton_pt[h][3] || lepton_pt[h][i] == lepton_pt[h][5])) continue;
									// Select DYT only if DYT = TuneP    	 						
// 									if(i == 4 && lepton_pt[h][i] == lepton_pt[h][5]) continue;
//      								if(i != 7) count_assignment_DATA_till400[i][j-14]++;
	     							pt_DATA[i]->Fill(lepton_pt[h][i]);
	     							eta_DATA[i]->Fill(fabs(lepton_eta[h][i]));
	     							phi_DATA[i]->Fill(lepton_phi[h][i]);
    	 							pt_vs_eta_DATA[i]->Fill(lep_eta[h], lep_pt[h]);
    	 							pt_vs_phi_DATA[i]->Fill(lep_phi[h], lep_pt[h]);
    	 							eta_vs_phi_DATA[i]->Fill(lep_phi[h], lep_eta[h]);
    	 							pt_vs_met_DATA[i]->Fill(met_pt, lep_pt[h]);

//     	 							pt_DATA[7]->Fill(lepton_pt[h][7]);
// 	     							eta_DATA[7]->Fill(fabs(lepton_eta[h][7]));
//      								phi_DATA[7]->Fill(lepton_phi[h][7]);
//    	 								pt_vs_eta_DATA[7]->Fill(lep_eta[h], lep_pt[h]);
//    	 								pt_vs_phi_DATA[7]->Fill(lep_phi[h], lep_pt[h]);
//    		 							eta_vs_phi_DATA[7]->Fill(lep_phi[h], lep_eta[h]);
// 	   	 							pt_vs_met_DATA[7]->Fill(met_pt, lep_pt[h]);

		     						if(lepton_pt[h][i] < 400){
	     								count_assignment_DATA_till400[i]++;
// 	    	 							count_assignment_DATA_till400[7]++;
		     							eta_DATA_till400[i]->Fill(fabs(lepton_eta[h][i]));
		     							phi_DATA_till400[i]->Fill(lepton_phi[h][i]);
// 		    	 						eta_DATA_till400[7]->Fill(fabs(lepton_eta[h][7]));
// 	    	 							phi_DATA_till400[7]->Fill(lepton_phi[h][7]);

		     						}
	    	 						else if(lepton_pt[h][i] > 400 && lepton_pt[h][i] < 600){
	     								count_assignment_DATA_400to600[i]++;
// 	    	 							count_assignment_DATA_400to600[7]++;
	     								eta_DATA_400to600[i]->Fill(fabs(lepton_eta[h][i]));
	     								phi_DATA_400to600[i]->Fill(lepton_phi[h][i]);
// 	     								eta_DATA_400to600[7]->Fill(fabs(lepton_eta[h][7]));
// 	    	 							phi_DATA_400to600[7]->Fill(lepton_phi[h][7]);
	     							}

	    	 						else{
	     								count_assignment_DATA_above600[i]++;
// 	    	 							count_assignment_DATA_above600[7]++;
	     								eta_DATA_above600[i]->Fill(fabs(lepton_eta[h][i]));
	     								phi_DATA_above600[i]->Fill(lepton_phi[h][i]);
// 	     								eta_DATA_above600[7]->Fill(fabs(lepton_eta[h][7]));
// 	    	 							phi_DATA_above600[7]->Fill(lepton_phi[h][7]);
	     							}
//     	 							break; // to take only one reconstruction in case of double counting --- DYT chosen because it is the first
    	 						}
     						}
    /*	 					if(lep_pt[h] == lepton_pt[h][7]){
//     	 						count_assignment_MC_till400[7]++;
     							pt_DATA[7]->Fill(lepton_pt[h][7]);
     							eta_DATA[7]->Fill(lepton_eta[h][7]);
     							phi_DATA[7]->Fill(lepton_phi[h][7]);
     							pt_DATA_clear[7]->Fill(lepton_pt[h][7]); //for stack plot
     							eta_DATA_clear[7]->Fill(lepton_eta[h][7]); //for stack plot
     							phi_DATA_clear[7]->Fill(lepton_phi[h][7]); //for stack plot
   	 							pt_vs_eta_DATA[7]->Fill(lep_eta[h], lep_pt[h]);
   	 							pt_vs_phi_DATA[7]->Fill(lep_phi[h], lep_pt[h]);
   	 							eta_vs_phi_DATA[7]->Fill(lep_phi[h], lep_eta[h]);
   	 							pt_vs_met_DATA[7]->Fill(met_pt, lep_pt[h]);


		     					if(lepton_pt[h][7] < 600){
	    	 						count_assignment_DATA_till400[7]++;
		     						eta_DATA_till400[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA_till400[7]->Fill(lepton_phi[h][7]);
	     							eta_DATA_till400_clear[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA_till400_clear[7]->Fill(lepton_phi[h][7]);

		     					}
	    	 					else{
	    	 						count_assignment_DATA_above600[7]++;
	     							eta_DATA_above600[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA_above600[7]->Fill(lepton_phi[h][7]);
     								eta_DATA_above600_clear[7]->Fill(lepton_eta[h][7]);
     								phi_DATA_above600_clear[7]->Fill(lepton_phi[h][7]);

	     						}
     						} */


//      						for(int i =0; i < 8; i++)
//      							std::cout<<reconstruction[i].Data()<<"\t"<<count_assignment_DATA_till400[i][j-30]<<std::endl;
//      						std::cout<<"-"<<std::endl;
     					}
     					else NotGoodQualityMuon_DATA++;
     				} // for on muons
// 			if(DeltaPhi[0]*DeltaPhi[1] != 0){ 
// 				countDeltaPhi_DATA++; 
// 				DeltaPhi_DATA_till400->Fill(fabs(lepton_phi[0][7]-lepton_phi[1][7]));
// 			}
// 			else NOcountDeltaPhi_DATA++;

	       }//end of condition on event

		 }// end loop p


std::cout<<"MC: TOT = "<<countDeltaPhi_MC+NOcountDeltaPhi_MC<<"\t SI = "
			<<countDeltaPhi_MC<<"\t"<<(float)countDeltaPhi_MC/(float)(countDeltaPhi_MC+NOcountDeltaPhi_MC)<<"\t NO = "
			<<NOcountDeltaPhi_MC<<"\t"<<(float)NOcountDeltaPhi_MC/(float)(countDeltaPhi_MC+NOcountDeltaPhi_MC)<<"\t"<<std::endl;
std::cout<<"DATA: TOT = "<<countDeltaPhi_DATA+NOcountDeltaPhi_DATA<<"\t SI = "
			<<countDeltaPhi_DATA<<"\t"<<(float)countDeltaPhi_DATA/(float)(countDeltaPhi_DATA+NOcountDeltaPhi_DATA)<<"\t NO = "
			<<NOcountDeltaPhi_DATA<<"\t"<<(float)NOcountDeltaPhi_DATA/(float)(countDeltaPhi_DATA+NOcountDeltaPhi_DATA)<<"\t"<<std::endl;


     float countDATA_err[8] = {0};
     float countDATA_TOT_err[8] = {0};
     float first = 0;
     float second = 0;


     std::cout<<" Till 400 GeV "<<std::endl;
     printf("                --> %15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[7].Data()); 

    for(int i = 0; i < 31; i++){ // samples
 		printf("%15s -->  %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", samples[14+i].Data(), count_assignment_MC_till400[1][i], count_assignment_MC_till400[2][i], count_assignment_MC_till400[3][i], count_assignment_MC_till400[4][i], count_assignment_MC_till400[5][i], count_assignment_MC_till400[7][i]);
	   	for(int j = 1; j < 8; j++) // reconstruction
	    	tot_till[j] += count_assignment_MC_till400[j][i];
     }
     std::cout<<"OK"<<std::endl;
     for(int i = 1; i < 7; i++){
     	TOT_MC+=tot_till[i];
     	TOT_DATA+=count_assignment_DATA_till400[i];
     }     
     for(int i = 1; i< 8; i++){
     	countDATA_err[i] = sqrt(count_assignment_DATA_till400[i]);

		first = countDATA_err[i]/TOT_DATA;
		second = (count_assignment_DATA_till400[i] * sqrt(TOT_DATA)) / pow(TOT_DATA, 2);
		countDATA_TOT_err[i] = sqrt(pow(first, 2) + pow(second, 2));
	}

 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("    MC         --->  %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", tot_till[1], tot_till[2], tot_till[3], tot_till[4],tot_till[5], tot_till[7]);
 	 printf("  DATA         --->  %.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\n", count_assignment_DATA_till400[1], countDATA_err[1], count_assignment_DATA_till400[2], countDATA_err[2], count_assignment_DATA_till400[3], countDATA_err[3], count_assignment_DATA_till400[4], countDATA_err[4], count_assignment_DATA_till400[5], countDATA_err[5], count_assignment_DATA_till400[7], countDATA_err[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
 	 printf("  MC TOTAL = %10.3f\nDATA TOTAL = %f\n", TOT_MC, TOT_DATA); 	 
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("   MC (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)tot_till[1]/TOT_MC, (float)tot_till[2]/TOT_MC,(float)tot_till[3]/TOT_MC, (float)tot_till[4]/TOT_MC, (float)tot_till[5]/TOT_MC, (float)tot_till[7]/TOT_MC);
 	 printf(" DATA (%%)      ---> %.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\n", (float)count_assignment_DATA_till400[1]/TOT_DATA, countDATA_TOT_err[1], (float)count_assignment_DATA_till400[2]/TOT_DATA, countDATA_TOT_err[2], (float)count_assignment_DATA_till400[3]/TOT_DATA, countDATA_TOT_err[3], (float)count_assignment_DATA_till400[4]/TOT_DATA, countDATA_TOT_err[4], (float)count_assignment_DATA_till400[5]/TOT_DATA, countDATA_TOT_err[5], (float)count_assignment_DATA_till400[7]/TOT_DATA, countDATA_TOT_err[7]);

     std::cout<<"No good muon MC   = "<<NotGoodQualityMuon_MC<<std::endl;
     std::cout<<"No good muon DATA = "<<NotGoodQualityMuon_DATA<<std::endl;
     
/*     std::cout<<"\n\n\t\t\t\tMC"<<std::endl;
     printf("                    %15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[0].Data(), reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[6].Data(), reconstruction[7].Data()); 
    for(int i = 0; i < 8; i++){ 	
 		if(i == 7)
 		 	 printf(" --------------------------------------------------------------------------------------------------------------------------------------------- \n");
 		printf("%15s -->  %15d\t%15d\t%15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", reconstruction[i].Data(), count_double_MC[0][i], count_double_MC[1][i], count_double_MC[2][i], count_double_MC[3][i], count_double_MC[4][i], count_double_MC[5][i], count_double_MC[6][i], count_double_MC[7][i]);

 	}
     std::cout<<"\n\t\t\t\tDATA"<<std::endl;
     printf("                    %15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[0].Data(), reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[6].Data(), reconstruction[7].Data()); 
    for(int i = 0; i < 8; i++){ 	
 		if(i == 7)
 		 	 printf(" --------------------------------------------------------------------------------------------------------------------------------------------- \n");
 		printf("%15s -->  %15d\t%15d\t%15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", reconstruction[i].Data(), count_double_DATA[0][i], count_double_DATA[1][i], count_double_DATA[2][i], count_double_DATA[3][i], count_double_DATA[4][i], count_double_DATA[5][i], count_double_DATA[6][i], count_double_DATA[7][i]);

 	}*/
 	

     std::cout<<"Between 400 and 600 GeV "<<std::endl;
     printf("                --> %15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[7].Data()); 

    for(int i = 0; i < 31; i++){ // samples
 		printf("%15s -->  %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", samples[14+i].Data(), count_assignment_MC_400to600[1][i], count_assignment_MC_400to600[2][i], count_assignment_MC_400to600[3][i], count_assignment_MC_400to600[4][i], count_assignment_MC_400to600[5][i], count_assignment_MC_400to600[7][i]);
	   	for(int j = 1; j < 8; j++) // reconstruction
	    	tot_middle[j] += count_assignment_MC_400to600[j][i];
     }
     TOT_MC = 0;
     TOT_DATA = 0;
     for(int i = 1; i < 7; i++){
     	TOT_MC+=tot_middle[i];
     	TOT_DATA+=count_assignment_DATA_400to600[i];
     }
     for(int i = 1; i< 8; i++){
     	countDATA_err[i] = sqrt(count_assignment_DATA_400to600[i]);

		first = countDATA_err[i]/TOT_DATA;
		second = (count_assignment_DATA_400to600[i] * sqrt(TOT_DATA)) / pow(TOT_DATA, 2);
		countDATA_TOT_err[i] = sqrt(pow(first, 2) + pow(second, 2));
	}

 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("    MC         --->  %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", tot_middle[1], tot_middle[2], tot_middle[3], tot_middle[4],tot_middle[5], tot_middle[7]);
 	 printf("  DATA         --->  %.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\n", count_assignment_DATA_400to600[1], countDATA_err[1], count_assignment_DATA_400to600[2], countDATA_err[2], count_assignment_DATA_400to600[3], countDATA_err[3], count_assignment_DATA_400to600[4], countDATA_err[4], count_assignment_DATA_400to600[5], countDATA_err[5], count_assignment_DATA_400to600[7], countDATA_err[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
 	 printf("  MC TOTAL = %10.3f\nDATA TOTAL = %f\n", TOT_MC, TOT_DATA); 	 
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("   MC (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)tot_middle[1]/TOT_MC, (float)tot_middle[2]/TOT_MC,(float)tot_middle[3]/TOT_MC, (float)tot_middle[4]/TOT_MC, (float)tot_middle[5]/TOT_MC, (float)tot_middle[7]/TOT_MC);
 	 printf(" DATA (%%)      ---> %.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\n", (float)count_assignment_DATA_400to600[1]/TOT_DATA, countDATA_TOT_err[1], (float)count_assignment_DATA_400to600[2]/TOT_DATA, countDATA_TOT_err[2], (float)count_assignment_DATA_400to600[3]/TOT_DATA, countDATA_TOT_err[3], (float)count_assignment_DATA_400to600[4]/TOT_DATA, countDATA_TOT_err[4], (float)count_assignment_DATA_400to600[5]/TOT_DATA, countDATA_TOT_err[5], (float)count_assignment_DATA_400to600[7]/TOT_DATA, countDATA_TOT_err[7]);


 	
     std::cout<<"Above 600 GeV "<<std::endl;
     printf("                --> %15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[7].Data()); 

    for(int i = 0; i < 31; i++){ // samples
 		printf("%15s -->  %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", samples[14+i].Data(), count_assignment_MC_above600[1][i], count_assignment_MC_above600[2][i], count_assignment_MC_above600[3][i], count_assignment_MC_above600[4][i], count_assignment_MC_above600[5][i], count_assignment_MC_above600[7][i]);
	   	for(int j = 1; j < 8; j++) // reconstruction
	    	tot_above[j] += count_assignment_MC_above600[j][i];
     }
     TOT_MC = 0;
     TOT_DATA = 0;
     for(int i = 1; i < 7; i++){
     	TOT_MC+=tot_above[i];
     	TOT_DATA+=count_assignment_DATA_above600[i];
     }
     for(int i = 1; i< 8; i++){
     	countDATA_err[i] = sqrt(count_assignment_DATA_above600[i]);

		first = countDATA_err[i]/TOT_DATA;
		second = (count_assignment_DATA_above600[i] * sqrt(TOT_DATA)) / pow(TOT_DATA, 2);
		countDATA_TOT_err[i] = sqrt(pow(first, 2) + pow(second, 2));
	}
	
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("    MC         --->  %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", tot_above[1], tot_above[2], tot_above[3], tot_above[4],tot_above[5], tot_above[7]);
 	 printf("  DATA         --->  %.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\n", count_assignment_DATA_above600[1], countDATA_err[1], count_assignment_DATA_above600[2], countDATA_err[2], count_assignment_DATA_above600[3], countDATA_err[3], count_assignment_DATA_above600[4], countDATA_err[4], count_assignment_DATA_above600[5], countDATA_err[5], count_assignment_DATA_above600[7], countDATA_err[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
 	 printf("  MC TOTAL = %10.3f\nDATA TOTAL = %f\n", TOT_MC, TOT_DATA); 	 
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("   MC (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)tot_above[1]/TOT_MC, (float)tot_above[2]/TOT_MC,(float)tot_above[3]/TOT_MC, (float)tot_above[4]/TOT_MC, (float)tot_above[5]/TOT_MC, (float)tot_above[7]/TOT_MC);
 	 printf(" DATA (%%)      ---> %.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\t%.3f+/-%.3f\n", (float)count_assignment_DATA_above600[1]/TOT_DATA, countDATA_TOT_err[1], (float)count_assignment_DATA_above600[2]/TOT_DATA, countDATA_TOT_err[2], (float)count_assignment_DATA_above600[3]/TOT_DATA, countDATA_TOT_err[3], (float)count_assignment_DATA_above600[4]/TOT_DATA, countDATA_TOT_err[4], (float)count_assignment_DATA_above600[5]/TOT_DATA, countDATA_TOT_err[5], (float)count_assignment_DATA_above600[7]/TOT_DATA, countDATA_TOT_err[7]);

     std::cout<<"No good muon MC   = "<<NotGoodQualityMuon_MC<<std::endl;
     std::cout<<"No good muon DATA = "<<NotGoodQualityMuon_DATA<<std::endl;
     
/*     std::cout<<"\n\n\t\t\t\tMC"<<std::endl;
     printf("                    %15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[0].Data(), reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[6].Data(), reconstruction[7].Data()); 
    for(int i = 0; i < 8; i++){ 	
 		if(i == 7)
 		 	 printf(" --------------------------------------------------------------------------------------------------------------------------------------------- \n");
 		printf("%15s -->  %15d\t%15d\t%15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", reconstruction[i].Data(), count_double_MC[0][i], count_double_MC[1][i], count_double_MC[2][i], count_double_MC[3][i], count_double_MC[4][i], count_double_MC[5][i], count_double_MC[6][i], count_double_MC[7][i]);

 	}
     std::cout<<"\n\t\t\t\tDATA"<<std::endl;
     printf("                    %15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[0].Data(), reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[6].Data(), reconstruction[7].Data()); 
    for(int i = 0; i < 8; i++){ 	
 		if(i == 7)
 		 	 printf(" --------------------------------------------------------------------------------------------------------------------------------------------- \n");
 		printf("%15s -->  %15d\t%15d\t%15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", reconstruction[i].Data(), count_double_DATA[0][i], count_double_DATA[1][i], count_double_DATA[2][i], count_double_DATA[3][i], count_double_DATA[4][i], count_double_DATA[5][i], count_double_DATA[6][i], count_double_DATA[7][i]);

 	}*/
 
    gROOT->Reset();
    gROOT->SetBatch();

 	 	
 	bool save = false;
//  	bool save = true;

	TString dir_save = "./PT_ASSIGNMENT_PLOT_allMC/";

// 	TFile *file = new TFile(dir_save + "Histos.root", "RECREATE");

	int a = 0;
	int b = 1;
	int c = 2;
	int d = 3;
	int e = 4;
// 	for(int k = 0; k < 5; k++){
// 		int j = k + 1;
// 		name_histo = Form("%s_MultipleCounting", reconstruction[j].Data());
// 		SaveMultipleCounting(name_histo, pt_CountDouble_MC[a], pt_CountDouble_DATA[a], dir_save, a, b, c, d, e);
// 		swap(a,b);
// 		swap(b,c);	
// 		swap(c,d);	
// 		swap(d,e);	
// 	}	

//   TCanvas* c4 = new TCanvas("mass", "mass", 200,10,700,500);
//   TH1F* last_hist = (TH1F *)eta_GEN_till400_Stack->GetStack()->Last();	
//   last_hist->Draw();
	SalvaHisto("Dimuon Mass: stack plot", Dimuon_MC_Stack, Dimuon_DATA, dir_save + "DimuonMass_StackPlot.C", 1,  1, "m_{#mu#mu} [GeV]");

	name_histo = Form("#eta %s: p_{T} < 400GeV", reconstruction[7].Data());
	SalvaHisto(name_histo, eta_MC_till400_Stack[7], eta_DATA_till400[7], dir_save + "TuneP_Eta400.C", 0,  0, "#eta");

	name_histo = Form("#eta %s: 400 < p_{T} < 600GeV", reconstruction[7].Data());
	SalvaHisto(name_histo, eta_MC_400to600_Stack[7], eta_DATA_400to600[7], dir_save + "TuneP_Eta400600.C", 0,  0, "#eta");

	name_histo = Form("#eta %s: p_{T} > 600GeV", reconstruction[7].Data());
	SalvaHisto(name_histo, eta_MC_above600_Stack[7], eta_DATA_above600[7], dir_save + "TuneP_Eta600.C", 0,  0, "#eta");

	name_histo = Form("#phi %s: p_{T} < 400GeV", reconstruction[7].Data());
	SalvaHisto(name_histo, phi_MC_till400_Stack[7], phi_DATA_till400[7], dir_save + "TuneP_Phi400.C", 0,  0,  "#phi");

	name_histo = Form("#phi %s: 400 < p_{T} < 600GeV", reconstruction[7].Data());
	SalvaHisto(name_histo, phi_MC_400to600_Stack[7], phi_DATA_400to600[7], dir_save + "TuneP_Phi400600.C", 0,  0,  "#phi");

	name_histo = Form("#phi %s: p_{T} > 600GeV", reconstruction[7].Data());
	SalvaHisto(name_histo, phi_MC_above600_Stack[7], phi_DATA_above600[7], dir_save + "TuneP_Phi600.C", 0,  0,  "#phi");

	name_histo = Form("%s p_{T}", reconstruction[7].Data());
	SalvaHisto(name_histo, pt_MC_Stack[7], pt_DATA[7], dir_save + "TuneP_Pt.C", 0,  1, "p_{T}");

	name_histo = Form("%s p_{T}Cut", reconstruction[7].Data());
	SalvaHisto("pTCut", pt_DATA_BeforeDeltaRcut, pt_DATA[7], dir_save + "PtCut.C", 0,  0, "#DeltaR");



	if(1==1){
	save_document = dir_save + "Variuos_distribution.pdf";

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf[", 0,  0, "p_{T}");
	SalvaHisto("DeltaR", hDeltaR, hDeltaR, dir_save + "Variuos_distribution.pdf[", 0,  0, "#DeltaR");
	SalvaHisto("pTCut", pt_DATA_BeforeDeltaRcut, pt_DATA[7], dir_save + "Variuos_distribution.pdf[", 0,  0, "#DeltaR");

// 	SalvaHisto("DeltaPhi: stack plot", DeltaPhi_MC_till400_Stack, DeltaPhi_DATA_till400, dir_save + "Variuos_distribution.pdf", 0,  0, "#phi");
	SalvaHisto("DeltaPhi: stack plot", DeltaPhi_DATA_till400, DeltaPhi_DATA_till400, dir_save + "Variuos_distribution.pdf", 0,  0, "#phi");
	SalvaHisto("GEN DeltaPhi: stack plot", DeltaPhi_GEN_till400_Stack, DeltaPhi_DATA_till400, dir_save + "Variuos_distribution.pdf", 0,  0, "#phi");

// 	SalvaHisto("gen Weight: stack plot", gen_Weight_Stack, pt_DATA[7], dir_save + "Variuos_distribution.pdf", 0,  1, "weight");
// 	SalvaHisto("genWeight vs eta", genW_vs_eta, dir_save + "Variuos_distribution.pdf", "gen W", "eta", 0, 0);
// 	SalvaHisto("genWeight vs phi", genW_vs_phi, dir_save + "Variuos_distribution.pdf", "gen W", "phi", 0, 0);
// 	SalvaHisto("genWeight vs pt", genW_vs_pt, dir_save + "Variuos_distribution.pdf", "gen W", "pt", 0, 0);
	SalvaHisto("Dilepton", DileptonPt_MC_Stack[7], DileptonPt_DATA[7], dir_save + "Variuos_distribution.pdf[", 0,  1, "p_{T}");
	SalvaHisto("#eta gen distribution (pT_{reco} < 400 GeV): stack plot", eta_GEN_till400_Stack, eta_DATA_till400[7], dir_save + "Variuos_distribution.pdf", 0,  0, "#eta");
	SalvaHisto("#eta gen distribution (400 < pT_{reco} < 600 GeV): stack plot", eta_GEN_400to600_Stack, eta_DATA_400to600[7], dir_save + "Variuos_distribution.pdf", 0,  0, "#eta");
	SalvaHisto("#eta gen distribution (pT_{reco} > 600 GeV): stack plot", eta_GEN_above600_Stack, eta_DATA_above600[7], dir_save + "Variuos_distribution.pdf", 0,  0, "#eta");
	SalvaHisto("#phi gen distribution (pT_{reco} < 400 GeV): stack plot", phi_GEN_till400_Stack, phi_DATA_till400[7], dir_save + "Variuos_distribution.pdf", 0,  0, "#phi");
	SalvaHisto("#phi gen distribution (400 < pT_{reco} < 600 GeV): stack plot", phi_GEN_400to600_Stack, phi_DATA_400to600[7], dir_save + "Variuos_distribution.pdf", 0,  0, "#phi");
	SalvaHisto("#phi gen distribution (pT_{reco} > 600 GeV): stack plot", phi_GEN_above600_Stack, phi_DATA_above600[7], dir_save + "Variuos_distribution.pdf", 0,  0, "#phi");
	SalvaHisto("p_{T} gen distribution: stack plot", pt_GEN_Stack, pt_DATA[7], dir_save + "Variuos_distribution.pdf", 0,  1, "p_{T}");
	SalvaHisto("Gen Dimuon Mass", GenDimuon_MC_Stack, Dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
// 									SalvaHisto("Dimuon Mass", Dimuon_MC, Dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
	SalvaHisto("Dimuon Mass: stack plot", Dimuon_MC_Stack, Dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");

// 	SalvaHisto("PickyDimuon Mass: stack plot", PickyDimuon_MC_Stack, Dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
// 	SalvaHisto("DYTDimuon Mass: stack plot", DYTDimuon_MC_Stack, Dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
// 	SalvaHisto("TrackerDimuon Mass: stack plot", TrackerDimuon_MC_Stack, Dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
// 	SalvaHisto("DimuonMass: TuneP_vs_Picky", DimuonTuneP_vs_Picky, dir_save + "Variuos_distribution.pdf", "Picky", "TuneP", 0, 0);
// 	SalvaHisto("DimuonMass: TuneP_vs_DYT", DimuonTuneP_vs_DYT, dir_save + "Variuos_distribution.pdf", "DYT", "TuneP", 0, 0);
// 	SalvaHisto("DimuonMass: TuneP_vs_Tracker", DimuonTuneP_vs_Tracker, dir_save + "Variuos_distribution.pdf", "Tracker", "TuneP", 0, 0);
// 								SalvaHisto("Phi_check", Phi_bin, dir_save + "Variuos_distribution.pdf", "reco", "gen", 0, 0);
// 								SalvaHisto("Phi_check_400to600", Phi_bin_46, dir_save + "Variuos_distribution.pdf", "reco", "gen", 0, 0);
// 								SalvaHisto("Phi_check_above600", Phi_bin_6, dir_save + "Variuos_distribution.pdf", "reco", "gen", 0, 0);



	for(int i = 1; i < 8; i++){		
		if(i != 2 && i != 4 && i != 5 && i != 7) continue;
	
// 		name_histo = Form("Eta_%s_till400GeV", reconstruction[i].Data());
// 		SalvaHisto(name_histo, eta_MC_till400[i], eta_DATA_till400[i], save_document, 0,  0, "#eta");
		name_histo = Form("#eta %s: p_{T} < 400GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC_till400_Stack[i], eta_DATA_till400[i], save_document, 0,  0, "#eta");

// 		name_histo = Form("Eta_%s_400to600GeV", reconstruction[i].Data());
// 		SalvaHisto(name_histo, eta_MC_400to600[i], eta_DATA_400to600[i], save_document, 0,  0, "#eta");
		name_histo = Form("#eta %s: 400 < p_{T} < 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC_400to600_Stack[i], eta_DATA_400to600[i], save_document, 0,  0, "#eta");

// 		name_histo = Form("Eta_%s_above600GeV", reconstruction[i].Data());
// 		SalvaHisto(name_histo, eta_MC_above600[i], eta_DATA_above600[i], save_document, 0,  0, "#eta");
		name_histo = Form("#eta %s: p_{T} > 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC_above600_Stack[i], eta_DATA_above600[i], save_document, 0,  0, "#eta");

// 		name_histo = Form("Phi_%s_till400GeV", reconstruction[i].Data());
// 		SalvaHisto(name_histo, phi_MC_till400[i], phi_DATA_till400[i], save_document, 0, 0,  "#phi");
		name_histo = Form("#phi %s: p_{T} < 400GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC_till400_Stack[i], phi_DATA_till400[i], save_document, 0,  0,  "#phi");

// 		name_histo = Form("Phi_%s_400to600GeV", reconstruction[i].Data());
// 		SalvaHisto(name_histo, phi_MC_400to600[i], phi_DATA_400to600[i], save_document, 0,  0, "#phi");
		name_histo = Form("#phi %s: 400 < p_{T} < 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC_400to600_Stack[i], phi_DATA_400to600[i], save_document, 0,  0,  "#phi");

// 		name_histo = Form("Phi_%s_above600GeV", reconstruction[i].Data());
// 		SalvaHisto(name_histo, phi_MC_above600[i], phi_DATA_above600[i], save_document, 0,  0, "#phi");
		name_histo = Form("#phi %s: p_{T} > 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC_above600_Stack[i], phi_DATA_above600[i], save_document, 0,  0,  "#phi");

// 		name_histo = Form("%s p_{T}", reconstruction[i].Data());
// 		SalvaHisto(name_histo, pt_MC[i], pt_DATA[i], save_document, 0,  1, "p_{T}");
		name_histo = Form("%s p_{T}", reconstruction[i].Data());
		SalvaHisto(name_histo, pt_MC_Stack[i], pt_DATA[i], save_document, 0,  1, "p_{T}");

// 		name_histo = Form("%s #eta", reconstruction[i].Data());
// 		SalvaHisto(name_histo, eta_MC[i], eta_DATA[i], save_document, 0,  0, "#eta");
		name_histo = Form("%s #eta", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC_Stack[i], eta_DATA[i], save_document, 0,  0, "#eta");

// 		name_histo = Form("%s #phi", reconstruction[i].Data());
// 		SalvaHisto(name_histo, phi_MC[i], phi_DATA[i], save_document, 0,  0, "#phi");
		name_histo = Form("%s #phi", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC_Stack[i], phi_DATA[i], save_document, 0,  0, "#phi");

		}
/*		name_histo = Form("%s p_{T} vs #eta: MC", reconstruction[i].Data());
		if(i != 7){
			pt_vs_eta_MC[i]->Divide(pt_vs_eta_MC[7]);
			pt_vs_eta_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_eta_MC[i], save_document, "#eta", "p_{T}");		
		name_histo = Form("%s p_{T} vs #eta: DATA", reconstruction[i].Data());
		if(i != 7){
			pt_vs_eta_DATA[i]->Divide(pt_vs_eta_DATA[7]);
			pt_vs_eta_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_eta_DATA[i], save_document, "#eta", "p_{T}");
		name_histo = Form("%s p_{T} vs #eta: DATA/MC", reconstruction[i].Data());
		pt_vs_eta_DATA[i]->Divide(pt_vs_eta_MC[i]);
		pt_vs_eta_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, pt_vs_eta_DATA[i], save_document, "#eta", "p_{T}", 0);


		name_histo = Form("%s p_{T} vs #phi: MC", reconstruction[i].Data());
		if(i != 7){
			pt_vs_phi_MC[i]->Divide(pt_vs_phi_MC[7]);
			pt_vs_phi_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_phi_MC[i], save_document, "#phi", "p_{T}");		
		name_histo = Form("%s p_{T} vs #phi: DATA", reconstruction[i].Data());
		if(i != 7){
			pt_vs_phi_DATA[i]->Divide(pt_vs_phi_DATA[7]);
			pt_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_phi_DATA[i], save_document, "#phi", "p_{T}");
		name_histo = Form("%s p_{T} vs #phi: DATA/MC", reconstruction[i].Data());
		pt_vs_phi_DATA[i]->Divide(pt_vs_phi_MC[i]);
		pt_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);
		SalvaHisto(name_histo, pt_vs_phi_DATA[i], save_document, "#phi", "p_{T}", 0);


		name_histo = Form("%s #eta vs #phi: MC", reconstruction[i].Data());
		if(i != 7){
			eta_vs_phi_MC[i]->Divide(eta_vs_phi_MC[7]);
			eta_vs_phi_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, eta_vs_phi_MC[i], save_document, "#phi", "#eta");		
		name_histo = Form("%s #eta vs #phi: DATA", reconstruction[i].Data());
		if(i != 7){
			eta_vs_phi_DATA[i]->Divide(eta_vs_phi_DATA[7]);
			eta_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, eta_vs_phi_DATA[i], save_document, "#phi", "#eta");
		name_histo = Form("%s #eta vs #phi: DATA/MC", reconstruction[i].Data());
		eta_vs_phi_DATA[i]->Divide(eta_vs_phi_MC[i]);
		eta_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, eta_vs_phi_DATA[i], save_document, "#phi", "#eta", 0);


		name_histo = Form("%s p_{T} vs met: MC", reconstruction[i].Data());
		if(i != 7){
			pt_vs_met_MC[i]->Divide(pt_vs_met_MC[7]);
			pt_vs_met_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_met_MC[i], save_document, "met", "p_{T}");
		name_histo = Form("%s p_{T} vs met: DATA", reconstruction[i].Data());
		if(i != 7){
			pt_vs_met_DATA[i]->Divide(pt_vs_met_DATA[7]);
			pt_vs_met_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_met_DATA[i], save_document, "met", "p_{T}");
		name_histo = Form("%s p_{T} vs met: DATA/MC", reconstruction[i].Data());
		pt_vs_met_DATA[i]->Divide(pt_vs_met_MC[i]);
		pt_vs_met_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, pt_vs_met_DATA[i], save_document, "met", "p_{T}", 0);
*/
	}	

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf]", 0,  0, "p_{T}");

/*
	TCanvas *c_Double_Count_MC;
	TLatex lat;
	TString vv = " == ";
	for(int i = 0; i < 6; i++){
		TString cc = "c_";
		TString nome_canvas = cc +  reconstruction[i+1].Data(); 
		c_Double_Count_MC = new TCanvas(nome_canvas, nome_canvas, 1050, 750);
		c_Double_Count_MC->Divide(2,3);
		int z = 0;
		int h = 0;
		for(int j = 0; j < 6; j++){
			z++;
			if(i == j) continue;
			c_Double_Count_MC->cd(z);
			float max = Double_Count_DATA[i][j]->GetMaximum();
			if(max < Double_Count_MC[i][j]->GetMaximum())
				max = Double_Count_MC[i][j]->GetMaximum();
			max = 1.1 * max;
			Double_Count_MC[i][j]->Draw();
			Double_Count_MC[i][j]->SetTitle("");
			Double_Count_MC[i][j]->SetLineColor(kBlue);
			Double_Count_MC[i][j]->GetYaxis()->SetRangeUser(0, max);
			gPad->Update();
			TPaveStats *s_MC = (TPaveStats*)Double_Count_MC[i][j]->GetListOfFunctions()->FindObject("stats");
			s_MC->SetName("Const");
			s_MC->SetX1NDC(0.75);
			s_MC->SetY1NDC(0.52);
			s_MC->SetY2NDC(0.72);
			s_MC->SetTextColor(kBlue);
			Double_Count_DATA[i][j]->Draw();
			Double_Count_DATA[i][j]->SetTitle("");
			Double_Count_DATA[i][j]->SetLineColor(kRed);
			Double_Count_DATA[i][j]->GetYaxis()->SetRangeUser(0, max);
			gPad->Update();
			TPaveStats *st_DATA = (TPaveStats*)Double_Count_DATA[i][j]->GetListOfFunctions()->FindObject("stats");
			st_DATA->SetName("Const");
			st_DATA->SetX1NDC(0.75);
			st_DATA->SetY1NDC(0.75);
			st_DATA->SetY2NDC(0.95);
			st_DATA->SetTextColor(kRed);
// 			pt_MC[i]->GetYaxis()->SetTitleOffset(1.2);
// 			pt_DATA[i]->GetYaxis()->SetTitleOffset(1.2);
	    	Double_Count_MC[i][j]->Draw();
    		Double_Count_DATA[i][j]->Draw("same");	
			s_MC->Draw("same");
			st_DATA->Draw("same");
			TString nome_lat = reconstruction[i+1].Data() + vv + reconstruction[z].Data();
			max = 0.9 * max / 1.1;
			if(max == 0)
				max = 0.9;
			lat.DrawLatex(850, max, nome_lat);

		}
		if(save){
			if(i==0)
		 		c_Double_Count_MC->Print(dir_save + "PT_DoubleCounting_distribution.pdf[");
			c_Double_Count_MC->Print(dir_save + "PT_DoubleCounting_distribution.pdf");
 			if(i==5)
				c_Double_Count_MC->Print(dir_save + "PT_DoubleCounting_distribution.pdf]");
		    c_Double_Count_MC->Write();
		}
	}
*/

// 	DimuonTuneP_vs_Picky->Write();
// 	DimuonTuneP_vs_DYT->Write();
// 	DimuonTuneP_vs_Tracker->Write();
// 	
// 	file->Close();	
	
}// end of function 
























void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA[15], TString save, int a, int b, int c, int d, int e){
	TLatex lat;
	TString nome_lat[5] = {"Global", "Picky", "Tpfms", "DYT", "Tracker"};
	TString inizio = Form("%s%s.pdf[", save.Data(), name.Data());
	TString fine = Form("%s%s.pdf]", save.Data(), name.Data());
	TString mezzo = Form("%s%s.pdf", save.Data(), name.Data());
	
	float max = 0;
	
		SalvaHisto("h_blank", h_blank, h_blank, inizio,  0,  0, "p_{T}");
		SalvaHisto("Multiple Assignment", h[0], h_DATA[0], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b]);
		SalvaHisto("Multiple Assignment", h[1], h_DATA[1], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c]);
		SalvaHisto("Multiple Assignment", h[2], h_DATA[2], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[3], h_DATA[3], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[4], h_DATA[4], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c]);
		SalvaHisto("Multiple Assignment", h[5], h_DATA[5], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[6], h_DATA[6], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[7], h_DATA[7], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[8], h_DATA[8], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[9], h_DATA[9], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[d], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[10], h_DATA[10], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[11], h_DATA[11], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[12], h_DATA[12], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[d], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[13], h_DATA[13], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c], nome_lat[d], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[14], h_DATA[14], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c], nome_lat[d], nome_lat[e]);
		SalvaHisto("h_blank", h_blank, h_blank, fine, 0, 0, "p_{T}");

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
	h_MC->Draw("SAME TEXT0");
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

void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex, TString strig_1, TString strig_2, TString strig_3, TString strig_4, TString strig_5){
			

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

	if(logy) min = 10e-1;
	else min = 0;

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

	TPaveStats * st_MC = (TPaveStats *)last_hist->GetListOfFunctions()->FindObject("stats");
    if( st_MC ){ 
		st_MC->SetName("Const");
		st_MC->SetX1NDC(0.75);
		st_MC->SetY1NDC(0.52);
		st_MC->SetY2NDC(0.72);
		st_MC->SetTextColor(kBlue);
    }
    else std::cout << "Null pointer to TPaveStats MC: " << st_MC << std::endl;
  
	h_DATA->SetLineColor(kBlack);
	h_DATA->SetTitle(name);
	h_DATA->Draw();
// 	h_DATA->SetStats(0);
	TH1F *ratio = (TH1F*) h_DATA->Clone();
	h_DATA->SetStats(1);
	c1->Update();
// 	gPad->Update();

	TPaveStats * st_DATA = (TPaveStats *)h_DATA->GetListOfFunctions()->FindObject("stats");
    if( st_DATA ){ 
    	st_DATA->SetTextColor(kRed); 
		st_DATA->SetName("Const");
		st_DATA->SetX1NDC(0.75);
		st_DATA->SetY1NDC(0.75);
		st_DATA->SetY2NDC(0.95);

    }
    else std::cout << "Null pointer to TPaveStats DATA: " << st_DATA << std::endl;

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
	st_MC->Draw("same");
// 	st_DATA->Draw("same");
//     	
	TLegend *l1 = new TLegend(0.25,0.8,0.35,0.9);
	l1->AddEntry(h_DATA, "DATA", "l");
	l1->AddEntry(h_MC, "MC", "l");
// 	l1->Draw();	
	
	lat.DrawLatex(X_pos_latex, max, strig_1);
	lat.DrawLatex(X_pos_latex, 0.9*max, strig_2);
	lat.DrawLatex(X_pos_latex, 0.8*max, strig_3);
	lat.DrawLatex(X_pos_latex, 0.7*max, strig_4);
	lat.DrawLatex(X_pos_latex, 0.6*max, strig_5);
	
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

void SalvaHisto(TString name, TH1F* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex, TString strig_1, TString strig_2, TString strig_3, TString strig_4, TString strig_5){

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

	if(logy) min = 0.1;
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

	TPaveStats * st_MC = (TPaveStats *)h_MC->GetListOfFunctions()->FindObject("stats");
    if( st_MC ){ 
		st_MC->SetName("Const");
		st_MC->SetX1NDC(0.75);
		st_MC->SetY1NDC(0.52);
		st_MC->SetY2NDC(0.72);
		st_MC->SetTextColor(kBlue);
    }
    else std::cout << "Null pointer to TPaveStats MC: " << st_MC << std::endl;
  
	h_DATA->SetLineColor(kRed);
	h_DATA->SetTitle(name);
	h_DATA->Draw();
// 	h_DATA->SetStats(0);
	TH1F *ratio = (TH1F*) h_DATA->Clone();
	h_DATA->SetStats(1);
	c1->Update();
// 	gPad->Update();

	TPaveStats * st_DATA = (TPaveStats *)h_DATA->GetListOfFunctions()->FindObject("stats");
    if( st_DATA ){ 
    	st_DATA->SetTextColor(kRed); 
		st_DATA->SetName("Const");
		st_DATA->SetX1NDC(0.75);
		st_DATA->SetY1NDC(0.75);
		st_DATA->SetY2NDC(0.95);

    }
    else std::cout << "Null pointer to TPaveStats DATA: " << st_DATA << std::endl;

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
	st_MC->Draw("same");
// 	st_DATA->Draw("same");
	
//     	
	TLegend *l1 = new TLegend(0.25,0.8,0.35,0.9);
	l1->AddEntry(h_DATA, "DATA", "l");
	l1->AddEntry(h_MC, "MC", "l");
	l1->Draw();	
	
	lat.DrawLatex(X_pos_latex, max, strig_1);
	lat.DrawLatex(X_pos_latex, 0.9*max, strig_2);
	lat.DrawLatex(X_pos_latex, 0.8*max, strig_3);
	lat.DrawLatex(X_pos_latex, 0.7*max, strig_4);
	lat.DrawLatex(X_pos_latex, 0.6*max, strig_5);
	
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
