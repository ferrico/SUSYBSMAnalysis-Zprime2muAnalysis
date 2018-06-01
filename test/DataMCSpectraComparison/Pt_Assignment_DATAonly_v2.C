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
#include <algorithm> 
#include <TString.h>

TString save_document;
TString name_histo;
void SalvaHisto(TString name, TH1F* h_DATA16, TH1F* h_DATA17, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, THStack* h_DATA16, TH1F* h_DATA17, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, TH2F* h_DATA16, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic = 1, bool logX = false);
void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA17[15], TString save, int a, int b, int c, int d, int e);

  Int_t event;
  Int_t run;
  Int_t prev_event = -88;
//   Int_t lumi;

  float dil_mass;
  
  TLorentzVector lep_1, lep_2;
  TLorentzVector ZPrime; 
  
//   float dil_pt;
  float cos_angle;
//   float vertex_chi2;
//   int dil_chosen;
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
//   short lep_stanAlone_numberOfBadHits[2];
//   short lep_stanAlone_numberOfMuonHits[2];
  short lep_numberOfMatchedStations[2];
  unsigned int lep_stationMask[2];
  short lep_numberOfMatchedRPCLayers[2];
  bool lep_isGlobalMuon[2];
  bool lep_isTrackerMuon[2];
  float lep_triggerMatchPt[2];
  float vertex_chi2;
  
  float met_pt;

	TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

	float LUMINOSITY_2016 = 36235.493;
	float LUMINOSITY_2017 = 41903.837;

	int count_event_DATA16[6][39] = {0};
	int count_event_DATA17[6] = {0};
		
	float lepton_pt[2][8] = {0};
	float lepton_eta[2][8] = {0};
	float lepton_phi[2][8] = {0};
	int numberOfValidMuonHits[2][8];
	int NotGoodQualityMuon_DATA16 = 0;
	int NotGoodQualityMuon_DATA17 = 0;
	int count_assignment_DATA16_till400[8];
	int count_assignment_DATA16_400to600[8];
	int count_assignment_DATA16_above600[8];
	int count_assignment_DATA17_till400[8];
	int count_assignment_DATA17_400to600[8];
	int count_assignment_DATA17_above600[8];

	float tot_till[8] = {0};
	float tot_middle[8] = {0};
	float tot_above[8] = {0};
	float TOT_DATA16 = 0;
	int TOT_DATA17 = 0;
	
	TString reconstruction[8] = {"Selected", "global", "picky", "tpfms", "dyt", "tracker", "std", "tuneP"};
	TString pt_name[2] = {"DATA16_pt_", "DATA17_pt_"};
	TString eta_name[2] = {"DATA16_eta_", "DATA17_eta_"};
	TString phi_name[2] = {"DATA16_phi_", "DATA17_phi_"};
	TString pt_vs_eta_name[2] = {"DATA16_pt_vs_eta_", "DATA17_ptVSeta_"};//DATA17_pt_vs_eta_"};
	TString pt_vs_phi_name[2] = {"DATA16_pt_vs_phi_", "DATA17_ptVSphi_"};//DATA17_pt_vs_phi_"};
	TString eta_vs_phi_name[2] = {"DATA16_eta_vs_phi_", "DATA17_etaVSphi_"};//DATA17_eta_vs_phi_"};
	TString pt_vs_met_name[2] = {"DATA16_pt_vs_met_", "DATA17_ptVSmet_"};//DATA17_pt_vs_met_"};
	TH1F* pt_DATA16[8];
	TH1F* pt_DATA17[8];
	TH1F* eta_DATA16[8];
	TH1F* eta_DATA17[8];
	TH1F* phi_DATA16[8];
	TH1F* phi_DATA17[8];
	TH2F* pt_vs_eta_DATA16[8];
	TH2F* pt_vs_eta_DATA17[8];
	TH2F* pt_vs_phi_DATA16[8];
	TH2F* pt_vs_phi_DATA17[8];
	TH2F* eta_vs_phi_DATA16[8];
	TH2F* eta_vs_phi_DATA17[8];
	TH2F* pt_vs_met_DATA16[8];
	TH2F* pt_vs_met_DATA17[8];
	TH1F* Double_Count_DATA16[6][6];
	TH1F* Double_Count_DATA17[6][6];

	TH1F* eta_DATA16_till400[8];
	TH1F* eta_DATA17_till400[8];
	TH1F* eta_DATA16_400to600[8];
	TH1F* eta_DATA17_400to600[8];
	TH1F* eta_DATA16_above600[8];
	TH1F* eta_DATA17_above600[8];

	TH1F* phi_DATA16_till400[8];
	TH1F* phi_DATA17_till400[8];
	TH1F* phi_DATA16_400to600[8];
	TH1F* phi_DATA17_400to600[8];
	TH1F* phi_DATA16_above600[8];
	TH1F* phi_DATA17_above600[8];
	
	TH1D* h_count_assignment_DATA16_till400[8];
	TH1D* h_count_assignment_DATA17_till400[8];
	
	TH1F* pt_CountDouble_DATA16[5][15];
	TH1F* pt_CountDouble_DATA17[5][15];
	
	const int    NMBINS = 100;
	const double MMIN = 60., MMAX = 2100.;
	double logMbins[NMBINS+1];

//     Double_t PT_BINS[] = {200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 600, 750, 1000, 1500};
//     Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;
    Double_t PT_BINS[] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 
    					1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800};
    Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;    
    
//     Double_t ETA_BINS[] = {-2.4, -1.5, -1.2, -0.9, 0, 0.9, 1.2, 1.5, 2.4};
    Double_t ETA_BINS[] = {0, 0.9, 1.2, 1.5, 2.4};
    Int_t  binnum_eta = sizeof(ETA_BINS)/sizeof(Double_t)-1;
    
    Double_t PHI_BINS[] = {-3.14, -2.356, -1.57, -0.785, 0, 0.785, 1.57, 2.356, 3.14};
    Int_t  binnum_phi = sizeof(PHI_BINS)/sizeof(Double_t)-1;
    
    Double_t MET_BINS[] = {0, 25, 50, 75, 100, 200, 350, 500};//, 750, 1000};
    Int_t  binnum_met = sizeof(MET_BINS)/sizeof(Double_t)-1;
   
    int count_double_DATA16[8][8] = {0};
    int count_double_DATA17[8][8] = {0};    
    
void Pt_Assignment_DATAonly_v2(){

    gROOT->Reset();
    gROOT->SetBatch();

	gStyle->SetOptFit(1);        
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "XYZ");
	gStyle->SetLabelSize(0.025, "XYZ");
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


	TH1F *Dimuon_DATA16 = new TH1F("Dimuon mass: DATA16", "Dimuon mass: DATA16", NMBINS, logMbins);
	TH1F *Dimuon_Picky_DATA16 = new TH1F("Dimuon mass Picky: DATA16", "Dimuon mass Picky: DATA16", NMBINS, logMbins);
	TH1F *Dimuon_Dyt_DATA16 = new TH1F("Dimuon mass Dyt: DATA16", "Dimuon mass Dyt: DATA16", NMBINS, logMbins);
	TH1F *Dimuon_DATA17 = new TH1F("Dimuon mass: DATA17", "Dimuon mass: DATA17", NMBINS, logMbins);
	TH1F *Dimuon_Picky_DATA17 = new TH1F("Dimuon mass Picky: DATA17", "Dimuon mass Picky: DATA17", NMBINS, logMbins);
	TH1F *Dimuon_Dyt_DATA17 = new TH1F("Dimuon mass Dyt: DATA17", "Dimuon mass Dyt: DATA17", NMBINS, logMbins);
	
	for(int a = 0; a < 5; a++){
		int b = a + 1;
		pt_CountDouble_DATA16[a][0] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][1] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][2] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][3] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][4] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][5] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][6] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][7] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][8] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][9] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);	
		pt_CountDouble_DATA16[a][10] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][11] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][12] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][13] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA16[a][14] = new TH1F(pt_name[0] + reconstruction[b], pt_name[0] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][0] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][1] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][2] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][3] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][4] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][5] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][6] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][7] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][8] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][9] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);	
		pt_CountDouble_DATA17[a][10] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][11] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][12] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][13] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
		pt_CountDouble_DATA17[a][14] = new TH1F(pt_name[1] + reconstruction[b], pt_name[1] + reconstruction[b], binnum_pt, PT_BINS);
	}
	

	for(int i = 0; i < 8; i++){

		pt_DATA16[i] = new TH1F(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i], binnum_pt, PT_BINS);
		pt_DATA17[i] = new TH1F(pt_name[1] + reconstruction[i], pt_name[1] + reconstruction[i], binnum_pt, PT_BINS);
		eta_DATA16[i] = new TH1F(eta_name[0] + reconstruction[i], eta_name[0] + reconstruction[i], binnum_eta, ETA_BINS);
		eta_DATA17[i] = new TH1F(eta_name[1] + reconstruction[i], eta_name[1] + reconstruction[i],binnum_eta, ETA_BINS);
		phi_DATA16[i] = new TH1F(phi_name[0] + reconstruction[i], phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS);
		phi_DATA17[i] = new TH1F(phi_name[1] + reconstruction[i], phi_name[1] + reconstruction[i],binnum_phi, PHI_BINS);
		pt_vs_eta_DATA16[i] = new TH2F(pt_vs_eta_name[0] + reconstruction[i], pt_vs_eta_name[0] + reconstruction[i], binnum_eta, ETA_BINS, binnum_pt, PT_BINS);
		pt_vs_eta_DATA17[i] = new TH2F(pt_vs_eta_name[1] + reconstruction[i], pt_vs_eta_name[1] + reconstruction[i], binnum_eta, ETA_BINS, binnum_pt, PT_BINS);
		pt_vs_phi_DATA16[i] = new TH2F(pt_vs_phi_name[0] + reconstruction[i], pt_vs_phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS, binnum_pt, PT_BINS);
		pt_vs_phi_DATA17[i] = new TH2F(pt_vs_phi_name[1] + reconstruction[i], pt_vs_phi_name[1] + reconstruction[i], binnum_phi, PHI_BINS, binnum_pt, PT_BINS);
		eta_vs_phi_DATA16[i] = new TH2F(eta_vs_phi_name[0] + reconstruction[i], eta_vs_phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS, binnum_eta, ETA_BINS);
		eta_vs_phi_DATA17[i] = new TH2F(eta_vs_phi_name[1] + reconstruction[i], eta_vs_phi_name[1] + reconstruction[i], binnum_phi, PHI_BINS, binnum_eta, ETA_BINS);

		eta_DATA16_till400[i] = new TH1F(eta_name[0] + reconstruction[i] + "_till400", eta_name[0] + reconstruction[i] + "_till400", binnum_eta, ETA_BINS);
		eta_DATA17_till400[i] = new TH1F(eta_name[1] + reconstruction[i] + "_till400", eta_name[1] + reconstruction[i] + "_till400",binnum_eta, ETA_BINS);
		eta_DATA16_400to600[i] = new TH1F(eta_name[0] + reconstruction[i] + "_400to600", eta_name[0] + reconstruction[i] + "_400to600", binnum_eta, ETA_BINS);
		eta_DATA17_400to600[i] = new TH1F(eta_name[1] + reconstruction[i] + "_400to600", eta_name[1] + reconstruction[i] + "_400to600", binnum_eta, ETA_BINS);
		eta_DATA16_above600[i] = new TH1F(eta_name[0] + reconstruction[i] + "_above600", eta_name[0] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);
		eta_DATA17_above600[i] = new TH1F(eta_name[1] + reconstruction[i] + "_above600", eta_name[1] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);

		phi_DATA16_till400[i] = new TH1F(phi_name[0] + reconstruction[i] + "_till400", phi_name[0] + reconstruction[i] + "_till400", binnum_phi, PHI_BINS);
		phi_DATA17_till400[i] = new TH1F(phi_name[1] + reconstruction[i] + "_till400", phi_name[1] + reconstruction[i] + "_till400",binnum_phi, PHI_BINS);
		phi_DATA16_400to600[i] = new TH1F(phi_name[0] + reconstruction[i] + "_400to600", phi_name[0] + reconstruction[i] + "_400to600", binnum_phi, PHI_BINS);
		phi_DATA17_400to600[i] = new TH1F(phi_name[1] + reconstruction[i] + "_400to600", phi_name[1] + reconstruction[i] + "_400to600", binnum_phi, PHI_BINS);
		phi_DATA16_above600[i] = new TH1F(phi_name[0] + reconstruction[i] + "_above600", phi_name[0] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);
		phi_DATA17_above600[i] = new TH1F(phi_name[1] + reconstruction[i] + "_above600", phi_name[1] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);

		pt_vs_met_DATA16[i] = new TH2F(pt_vs_met_name[0] + reconstruction[i], pt_vs_met_name[0] + reconstruction[i], binnum_met, MET_BINS, binnum_pt, PT_BINS);
		pt_vs_met_DATA17[i] = new TH2F(pt_vs_met_name[1] + reconstruction[i], pt_vs_met_name[1] + reconstruction[i], binnum_met, MET_BINS, binnum_pt, PT_BINS);



		pt_DATA16[i]->GetYaxis()->SetTitle("Entries");
		pt_DATA17[i]->GetYaxis()->SetTitle("Entries");
		eta_DATA16[i]->GetYaxis()->SetTitle("Entries");
		eta_DATA17[i]->GetYaxis()->SetTitle("Entries");		
		phi_DATA16[i]->GetYaxis()->SetTitle("Entries");
		phi_DATA17[i]->GetYaxis()->SetTitle("Entries");		

		pt_DATA16[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		pt_DATA17[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		eta_DATA16[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA17[i]->GetXaxis()->SetTitle("#eta");
		phi_DATA16[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA17[i]->GetXaxis()->SetTitle("#phi");


		eta_DATA16_till400[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA16_400to600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA16_above600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA17_till400[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA17_400to600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA17_above600[i]->GetXaxis()->SetTitle("#eta");

		phi_DATA16_till400[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA16_400to600[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA16_above600[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA17_till400[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA17_400to600[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA17_above600[i]->GetXaxis()->SetTitle("#phi");
		
		pt_vs_met_DATA16[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		pt_vs_met_DATA16[i]->GetXaxis()->SetTitle("met [GeV]");
		pt_vs_met_DATA17[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
		pt_vs_met_DATA17[i]->GetXaxis()->SetTitle("met [GeV]");

	}

	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
// 			if(i == j) continue;
			TString vv = "_vs_";
			TString name_DATA16 = pt_name[0] + reconstruction[i+1].Data() + vv + reconstruction[j+1].Data();
			TString name_DATA17 = pt_name[1] + reconstruction[i+1].Data() + vv + reconstruction[j+1].Data();
// 			std::cout<<i<<j<<"\t"<<name<<std::endl;
			Double_Count_DATA16[i][j] = new TH1F(name_DATA16, name_DATA16, binnum_pt, PT_BINS);
			Double_Count_DATA16[i][j]->GetXaxis()->SetTitle("p_{T}");
			Double_Count_DATA17[i][j] = new TH1F(name_DATA17, name_DATA17, binnum_pt, PT_BINS);
			Double_Count_DATA17[i][j]->GetXaxis()->SetTitle("p_{T}");
		}
	}			
     		
   	TChain *treeDATA16 = new TChain("SimpleNtupler/t");

    treeDATA16->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/DATA_Pt_Assignment/ana_datamc_data_Pt_Assignment.root");
 
   	treeDATA16->SetBranchAddress("event", &event);
   	treeDATA16->SetBranchAddress("run", &run);
//       	treeDATA16->SetBranchAddress("lumi", &lumi);    	
  	treeDATA16->SetBranchAddress("dil_mass",&dil_mass);
    	
//      treeDATA16->SetBranchAddress("dil_pt",&dil_pt);
     treeDATA16->SetBranchAddress("cos_angle",&cos_angle);
//      treeDATA16->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeDATA16->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeDATA16->SetBranchAddress("nvertices",&nvertices);
     treeDATA16->SetBranchAddress("GoodVtx",&GoodVtx);

     treeDATA16->SetBranchAddress("lep_pt",lep_pt);
   	 treeDATA16->SetBranchAddress("lep_id",lep_id);
     treeDATA16->SetBranchAddress("lep_eta",lep_eta);
     treeDATA16->SetBranchAddress("lep_phi",lep_phi);
     treeDATA16->SetBranchAddress("lep_dB",lep_dB);
//      treeDATA16->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeDATA16->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
   	 treeDATA16->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeDATA16->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeDATA16->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeDATA16->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

     treeDATA16->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeDATA16->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
   	 treeDATA16->SetBranchAddress("lep_picky_numberOfValidMuonHits",lep_picky_numberOfValidMuonHits);
     treeDATA16->SetBranchAddress("lep_dyt_numberOfValidMuonHits",lep_dyt_numberOfValidMuonHits);
   	 treeDATA16->SetBranchAddress("lep_tpfms_numberOfValidMuonHits",lep_tpfms_numberOfValidMuonHits);
     treeDATA16->SetBranchAddress("lep_stanAlone_numberOfValidMuonHits",lep_stanAlone_numberOfValidMuonHits);

//      treeDATA16->SetBranchAddress("lep_stanAlone_numberOfBadHits",lep_stanAlone_numberOfBadHits);
//      treeDATA16->SetBranchAddress("lep_stanAlone_numberOfMuonHits",lep_stanAlone_numberOfMuonHits);
     treeDATA16->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
   	 treeDATA16->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeDATA16->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeDATA16->SetBranchAddress("lep_pfIso",lep_pfIso);
     treeDATA16->SetBranchAddress("lep_pt_err",lep_pt_err);
   	 treeDATA16->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
   	 treeDATA16->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeDATA16->SetBranchAddress("lep_glb_pt",lep_glb_pt);
   	 treeDATA16->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
     treeDATA16->SetBranchAddress("lep_picky_pt",lep_picky_pt);
   	 treeDATA16->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
   	 treeDATA16->SetBranchAddress("lep_std_pt",lep_std_pt);
   	 treeDATA16->SetBranchAddress("lep_cocktail_pt",lep_cocktail_pt);

     treeDATA16->SetBranchAddress("lep_glb_eta",lep_glb_eta);
   	 treeDATA16->SetBranchAddress("lep_dyt_eta",lep_dyt_eta); 
     treeDATA16->SetBranchAddress("lep_picky_eta",lep_picky_eta);
   	 treeDATA16->SetBranchAddress("lep_tpfms_eta",lep_tpfms_eta);
   	 treeDATA16->SetBranchAddress("lep_std_eta",lep_std_eta);
   	 treeDATA16->SetBranchAddress("lep_tk_eta",lep_tk_eta);
   	 treeDATA16->SetBranchAddress("lep_tuneP_eta",lep_tuneP_eta);
   	 
     treeDATA16->SetBranchAddress("lep_glb_phi",lep_glb_phi);
   	 treeDATA16->SetBranchAddress("lep_dyt_phi",lep_dyt_phi); 
     treeDATA16->SetBranchAddress("lep_picky_phi",lep_picky_phi);
   	 treeDATA16->SetBranchAddress("lep_tpfms_phi",lep_tpfms_phi);
   	 treeDATA16->SetBranchAddress("lep_std_phi",lep_std_phi);
   	 treeDATA16->SetBranchAddress("lep_tk_phi",lep_tk_phi);
   	 treeDATA16->SetBranchAddress("lep_tuneP_phi",lep_tuneP_phi);

     treeDATA16->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);   
     treeDATA16->SetBranchAddress("vertex_chi2",&vertex_chi2); 	

     treeDATA16->SetBranchAddress("met_pt",&met_pt); 	

     Long64_t nentries = treeDATA16->GetEntries();
     
     printf("opening... DATA 2016 --- %lld\n", nentries);
     
     float weight = LUMINOSITY_2017 / LUMINOSITY_2016;

   	 for(int p=0; p<nentries; p++){
//    	 for(int p=0; p<10000; p++){
//    	 for(int p=0; p<1; p++){

   	 	if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
    	 	
   	 	treeDATA16->GetEntry(p);

     	if (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
     		lep_pt[0]>53. && lep_pt[1]>53. && 
//      		lep_pt[0]>200. && lep_pt[1]>200. && 
   	 		lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
   			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
   			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
	   		lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
    		lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
     		(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
    	 		
//    	 		dil_mass > 120 &&
//    	 		(dil_mass < 70 || dil_mass > 110) &&
			dil_mass > 60 && dil_mass < 120 &&

     		lep_sumPt[0]/lep_tk_pt[0]<0.10 && 
   	 		lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
/*     		lep_sumPt[0]/lep_tk_pt[0]<0.05 &&  // Rel. trk iso < 0.05
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
					Dimuon_DATA16->Fill(dil_mass,  weight); 
					lep_1.SetPtEtaPhiM(lep_picky_pt[0], lep_picky_eta[0], lep_picky_phi[0], 0.105);
					lep_2.SetPtEtaPhiM(lep_picky_pt[1], lep_picky_eta[1], lep_picky_phi[1], 0.105);
					ZPrime = lep_1 + lep_2;
					Dimuon_Picky_DATA16->Fill(ZPrime.M(),  weight); 
					lep_1.SetPtEtaPhiM(lep_dyt_pt[0], lep_dyt_eta[0], lep_dyt_phi[0], 0.105);
					lep_2.SetPtEtaPhiM(lep_dyt_pt[1], lep_dyt_eta[1], lep_dyt_phi[1], 0.105);
					ZPrime = lep_1 + lep_2;
					Dimuon_Dyt_DATA16->Fill(ZPrime.M(),  weight); 
				}


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

// 					if(numberOfValidMuonHits[h][1] == 0) continue;

							////////////////////////////////////////////
							////////////////////////////////////////////		
							///////////// MULTIPLE COUNTING ////////////
							int a = 1;
							int b = 2;
							int c = 3;
							int d = 4;
							int e = 5;

							for(int f = 0; f < 5; f++){
								if(lepton_pt[h][a] == lepton_pt[h][7]){
									if(lepton_pt[h][a] == lepton_pt[h][b]){
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][14]->Fill(lepton_pt[h][a], weight); //std::cout<<" ABCDE  "<<std::endl; 
												else pt_CountDouble_DATA16[a-1][10]->Fill(lepton_pt[h][a], weight); //std::cout<<" ABCD  "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][11]->Fill(lepton_pt[h][a], weight); //std::cout<<"  ABCE "<<std::endl; 
												else pt_CountDouble_DATA16[a-1][4]->Fill(lepton_pt[h][a], weight); //std::cout<<"ABC   "<<std::endl; 
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][12]->Fill(lepton_pt[h][a], weight); //std::cout<<" ABDE  "<<std::endl; 
												else pt_CountDouble_DATA16[a-1][5]->Fill(lepton_pt[h][a], weight); //std::cout<<"ABD   "<<std::endl; 
										}
										else if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][6]->Fill(lepton_pt[h][a], weight); //std::cout<<" ABE  "<<std::endl; 
										else pt_CountDouble_DATA16[a-1][0]->Fill(lepton_pt[h][a], weight); //std::cout<<" AB  "<<std::endl; 

									}
									else{
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][13]->Fill(lepton_pt[h][a], weight); //std::cout<<" ACDE  "<<std::endl; 
												else pt_CountDouble_DATA16[a-1][7]->Fill(lepton_pt[h][a], weight); //std::cout<<"  ACD "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][8]->Fill(lepton_pt[h][a], weight); //std::cout<<" ACE  "<<std::endl; 
												else pt_CountDouble_DATA16[a-1][1]->Fill(lepton_pt[h][a], weight); //std::cout<<" AC  "<<std::endl;
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][9]->Fill(lepton_pt[h][a], weight); //std::cout<<"ADE"<<std::endl; 
											else pt_CountDouble_DATA16[a-1][2]->Fill(lepton_pt[h][a], weight); //std::cout<<" AD  "<<std::endl; 
										}
										else 
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA16[a-1][3]->Fill(lepton_pt[h][a], weight); //std::cout<<" AE  "<<std::endl; 
									}
								}
// 						else std::cout<<" no tunep"<<std::endl;
								swap(a,b);
								swap(b,c);	
								swap(c,d);	
								swap(d,e);	
// 								std::cout<<"\t\t\t\t\t\t\t\t\t\t\t"<<a<<b<<c<<d<<e<<std::endl;	
							}
							///////////// MULTIPLE COUNTING ////////////
							////////////////////////////////////////////
							////////////////////////////////////////////
	     					
	     			for(int i = 1; i < 7; i++){ // for on different reconstruction
   	 						    	 						
    	 				if(lep_pt[h] == lepton_pt[h][i]){
    	 					
    	 					bool multiple_reconstruction = false;	
    	 					for(int z = 1; z < 7; z++){
    	 						if(z == i) continue;
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
//      								if(i != 7) count_assignment_DATA16_till400[i][j-14]++;
   							pt_DATA16[i]->Fill(lepton_pt[h][i], weight);
   							eta_DATA16[i]->Fill(fabs(lepton_eta[h][i]), weight);
   							phi_DATA16[i]->Fill(lepton_phi[h][i], weight);
 							pt_vs_eta_DATA16[i]->Fill(lep_eta[h], lep_pt[h], weight);
 							pt_vs_phi_DATA16[i]->Fill(lep_phi[h], lep_pt[h], weight);
 							eta_vs_phi_DATA16[i]->Fill(lep_phi[h], lep_eta[h], weight);
 							pt_vs_met_DATA16[i]->Fill(met_pt, lep_pt[h], weight);

 							pt_DATA16[7]->Fill(lepton_pt[h][7], weight);
   							eta_DATA16[7]->Fill(fabs(lepton_eta[h][7]), weight);
							phi_DATA16[7]->Fill(lepton_phi[h][7], weight);
							pt_vs_eta_DATA16[7]->Fill(lep_eta[h], lep_pt[h], weight);
   	 						pt_vs_phi_DATA16[7]->Fill(lep_phi[h], lep_pt[h], weight);
   		 					eta_vs_phi_DATA16[7]->Fill(lep_phi[h], lep_eta[h], weight);
	   	 					pt_vs_met_DATA16[7]->Fill(met_pt, lep_pt[h], weight);
	   	 							
	   	 							
		     				if(lepton_pt[h][i] < 400){
   								count_assignment_DATA16_till400[i]++;
	    	 					count_assignment_DATA16_till400[7]++;
	   							eta_DATA16_till400[i]->Fill(fabs(lepton_eta[h][i]), weight);
	   							phi_DATA16_till400[i]->Fill(lepton_phi[h][i], weight);
	   							eta_DATA16_till400[7]->Fill(fabs(lepton_eta[h][7]), weight);
	    	 					phi_DATA16_till400[7]->Fill(lepton_phi[h][7], weight);
     						}
     						else if(lepton_pt[h][i] > 400 && lepton_pt[h][i] < 600){
   								count_assignment_DATA16_400to600[i]++;
	    	 					count_assignment_DATA16_400to600[7]++;
	     						eta_DATA16_400to600[i]->Fill(fabs(lepton_eta[h][i]), weight);
	     						phi_DATA16_400to600[i]->Fill(lepton_phi[h][i], weight);
	     						eta_DATA16_400to600[7]->Fill(fabs(lepton_eta[h][7]), weight);
	    	 					phi_DATA16_400to600[7]->Fill(lepton_phi[h][7], weight);
		     				}
	    	 				else{
	     						count_assignment_DATA16_above600[i]++;
	    	 					count_assignment_DATA16_above600[7]++;
	     						eta_DATA16_above600[i]->Fill(fabs(lepton_eta[h][i]), weight);
	     						phi_DATA16_above600[i]->Fill(lepton_phi[h][i], weight);
	     						eta_DATA16_above600[7]->Fill(fabs(lepton_eta[h][7]), weight);
	    	 					phi_DATA16_above600[7]->Fill(lepton_phi[h][7], weight);
    	 					}
// 	    	 							std::cout<<lepton_pt[h][i]<<"\t"<<count_assignment_DATA16_above600[7][j-14]<<"\t"<<count_assignment_DATA16_till400[7][j-14]<<std::endl;
//     	 							break; // to take only one reconstruction in case of double counting --- DYT chosen because it is the first
    	 				}
     				}
    /*	 					if(lep_pt[h] == lepton_pt[h][7]){
//     	 						count_assignment_DATA16_till400[7][j-14]++;
     							pt_DATA16[7]->Fill(lepton_pt[h][7], weight);
     							eta_DATA16[7]->Fill(lepton_eta[h][7], weight);
     							phi_DATA16[7]->Fill(lepton_phi[h][7], weight);
     							pt_DATA16_clear[7]->Fill(lepton_pt[h][7], weight); //for stack plot
     							eta_DATA16_clear[7]->Fill(lepton_eta[h][7], weight); //for stack plot
     							phi_DATA16_clear[7]->Fill(lepton_phi[h][7], weight); //for stack plot
   	 							pt_vs_eta_DATA16[7]->Fill(lep_eta[h], lep_pt[h], weight);
   	 							pt_vs_phi_DATA16[7]->Fill(lep_phi[h], lep_pt[h], weight);
   	 							eta_vs_phi_DATA16[7]->Fill(lep_phi[h], lep_eta[h], weight);
   	 							pt_vs_met_DATA16[7]->Fill(met_pt, lep_pt[h], weight);


		     					if(lepton_pt[h][7] < 600){
	    	 						count_assignment_DATA16_till400[7][j-14]++;
		     						eta_DATA16_till400[7]->Fill(lepton_eta[h][7], weight);
	     							phi_DATA16_till400[7]->Fill(lepton_phi[h][7], weight);
	     							eta_DATA16_till400_clear[7]->Fill(lepton_eta[h][7], weight);
	     							phi_DATA16_till400_clear[7]->Fill(lepton_phi[h][7], weight);

		     					}
	    	 					else{
	    	 						count_assignment_DATA16_above600[7][j-14]++;
	     							eta_DATA16_above600[7]->Fill(lepton_eta[h][7], weight);
	     							phi_DATA16_above600[7]->Fill(lepton_phi[h][7], weight);
     								eta_DATA16_above600_clear[7]->Fill(lepton_eta[h][7], weight);
     								phi_DATA16_above600_clear[7]->Fill(lepton_phi[h][7], weight);

	     						}
     						} */
//      						std::cout<<count_assignment_DATA16_till400[0][0]<<"\t"<<count_assignment_DATA16_till400[1][0]<<"\t"<<count_assignment_DATA16_till400[2][0]
//      						<<"\t"<<count_assignment_DATA16_till400[3][0]<<"\t"<<count_assignment_DATA16_till400[4][0]<<"\t"<<count_assignment_DATA16_till400[5][0]
//      						<<"\t"<<count_assignment_DATA16_till400[6][0]<<"\t"<<count_assignment_DATA16_till400[7][0]
//      						<<"\t\t\t"<<lepton_pt[h][0]<<"\t"<<lepton_pt[h][1]<<"\t"<<lepton_pt[h][2]<<"\t"<<lepton_pt[h][3]
//      						<<"\t"<<lepton_pt[h][4]<<"\t"<<lepton_pt[h][5]<<"\t"<<lepton_pt[h][6]<<"\t"<<lepton_pt[h][7]<<std::endl;

//      						for(int i =0; i < 8; i++)
//      							std::cout<<reconstruction[i].Data()<<"\t"<<count_assignment_DATA16_till400[i][j-30]<<std::endl;
//      						std::cout<<"-"<<std::endl;
     			}
     			else NotGoodQualityMuon_DATA16++;
			} // for on muons
      }// end loop p
	}// end loop on DATA16
    
    
////////      DATA17    ///////

	for(int i = 0; i < 7; i++)
		for(int j = 0; j < 2; j++)
			lepton_pt[j][i] = {0};
		
    TChain *treeDATA17 = new TChain("SimpleNtupler/t");
     treeDATA17->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA/ana_datamc_data.root");     	


     treeDATA17->SetBranchAddress("event",&event);
     treeDATA17->SetBranchAddress("run",&run);
//      treeDATA17->SetBranchAddress("lumi",&lumi);
     	
     treeDATA17->SetBranchAddress("dil_mass",&dil_mass);
//      treeDATA17->SetBranchAddress("dil_pt",&dil_pt);
	treeDATA17->SetBranchAddress("cos_angle",&cos_angle);
//      treeDATA17->SetBranchAddress("vertex_chi2",&vertex_chi2);
//      treeDATA17->SetBranchAddress("dil_chosen",&dil_chosen);
//      treeDATA17->SetBranchAddress("nvertices",&nvertices);
    treeDATA17->SetBranchAddress("GoodVtx",&GoodVtx);
// 
	treeDATA17->SetBranchAddress("lep_pt",lep_pt);
    treeDATA17->SetBranchAddress("lep_id",lep_id);
	treeDATA17->SetBranchAddress("lep_eta",lep_eta);
	treeDATA17->SetBranchAddress("lep_phi",lep_phi);
	treeDATA17->SetBranchAddress("lep_dB",lep_dB);
//      treeDATA17->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
	treeDATA17->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
    treeDATA17->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
	treeDATA17->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
	treeDATA17->SetBranchAddress("lep_stationMask",&lep_stationMask);
	treeDATA17->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	treeDATA17->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
    treeDATA17->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
    treeDATA17->SetBranchAddress("lep_picky_numberOfValidMuonHits",lep_picky_numberOfValidMuonHits);
    treeDATA17->SetBranchAddress("lep_dyt_numberOfValidMuonHits",lep_dyt_numberOfValidMuonHits);
    treeDATA17->SetBranchAddress("lep_tpfms_numberOfValidMuonHits",lep_tpfms_numberOfValidMuonHits);
    treeDATA17->SetBranchAddress("lep_stanAlone_numberOfValidMuonHits",lep_stanAlone_numberOfValidMuonHits);

//      treeDATA17->SetBranchAddress("lep_stanAlone_numberOfBadHits",lep_stanAlone_numberOfBadHits);
//      treeDATA17->SetBranchAddress("lep_stanAlone_numberOfMuonHits",lep_stanAlone_numberOfMuonHits);
	treeDATA17->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
    treeDATA17->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
	treeDATA17->SetBranchAddress("lep_sumPt",lep_sumPt);
	treeDATA17->SetBranchAddress("lep_pfIso",lep_pfIso);
	treeDATA17->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
    treeDATA17->SetBranchAddress("lep_tk_pt",lep_tk_pt);    
	treeDATA17->SetBranchAddress("lep_glb_pt",lep_glb_pt);
    treeDATA17->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
	treeDATA17->SetBranchAddress("lep_picky_pt",lep_picky_pt);
    treeDATA17->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
	treeDATA17->SetBranchAddress("lep_std_pt",lep_std_pt);
	treeDATA17->SetBranchAddress("lep_cocktail_pt",lep_cocktail_pt);
	
    treeDATA17->SetBranchAddress("lep_glb_eta",lep_glb_eta);
   	treeDATA17->SetBranchAddress("lep_dyt_eta",lep_dyt_eta); 
    treeDATA17->SetBranchAddress("lep_picky_eta",lep_picky_eta);
   	treeDATA17->SetBranchAddress("lep_tpfms_eta",lep_tpfms_eta);
   	treeDATA17->SetBranchAddress("lep_std_eta",lep_std_eta);
   	treeDATA17->SetBranchAddress("lep_tk_eta",lep_tk_eta);
  	treeDATA17->SetBranchAddress("lep_tuneP_eta",lep_tuneP_eta);
  	
    treeDATA17->SetBranchAddress("lep_glb_phi",lep_glb_phi);
   	treeDATA17->SetBranchAddress("lep_dyt_phi",lep_dyt_phi); 
    treeDATA17->SetBranchAddress("lep_picky_phi",lep_picky_phi);
   	treeDATA17->SetBranchAddress("lep_tpfms_phi",lep_tpfms_phi);
   	treeDATA17->SetBranchAddress("lep_std_phi",lep_std_phi);
   	treeDATA17->SetBranchAddress("lep_tk_phi",lep_tk_phi);
  	treeDATA17->SetBranchAddress("lep_tuneP_phi",lep_tuneP_phi);

	treeDATA17->SetBranchAddress("lep_pt_err",lep_pt_err);
	treeDATA17->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);   
	treeDATA17->SetBranchAddress("vertex_chi2",&vertex_chi2); 

	treeDATA17->SetBranchAddress("met_pt",&met_pt); 		

	nentries = treeDATA17->GetEntries();
	
	printf("opening... DATA17 --- %lld\n", nentries);    
	
	for(int p=0; p<nentries; p++){
// 	for(int p=0; p<10000; p++){
// 	for(int p=0; p<1; p++){

		if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
	
		treeDATA17->GetEntry(p);

	    if (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
	     		lep_pt[0]>53. && lep_pt[1]>53. && 
// 	     		lep_pt[0]>200. && lep_pt[1]>200. && 
    	 		lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
     			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
     			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
	     		lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
    	 		lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
    	 		(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
    	 		
//     	 		dil_mass > 120 &&
//     	 		(dil_mass < 70 || dil_mass > 110) &&
				dil_mass > 60 && dil_mass < 120 &&
    	 		
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
							Dimuon_DATA17->Fill(dil_mass,  weight); 
							lep_1.SetPtEtaPhiM(lep_picky_pt[0], lep_picky_eta[0], lep_picky_phi[0], 0.105);
							lep_2.SetPtEtaPhiM(lep_picky_pt[1], lep_picky_eta[1], lep_picky_phi[1], 0.105);
							ZPrime = lep_1 + lep_2;
							Dimuon_Picky_DATA17->Fill(ZPrime.M(),  weight); 
							lep_1.SetPtEtaPhiM(lep_dyt_pt[0], lep_dyt_eta[0], lep_dyt_phi[0], 0.105);
							lep_2.SetPtEtaPhiM(lep_dyt_pt[1], lep_dyt_eta[1], lep_dyt_phi[1], 0.105);
							ZPrime = lep_1 + lep_2;
							Dimuon_Dyt_DATA17->Fill(ZPrime.M(),  weight); 
						}

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
							for(int f = 0; f < 5; f++){
								if(lepton_pt[h][a] == lepton_pt[h][7]){
									if(lepton_pt[h][a] == lepton_pt[h][b]){
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][14]->Fill(lepton_pt[h][a]); //std::cout<<" ABCDE  "<<std::endl; 
												else pt_CountDouble_DATA17[a-1][10]->Fill(lepton_pt[h][a]); //std::cout<<" ABCD  "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][11]->Fill(lepton_pt[h][a]); //std::cout<<"  ABCE "<<std::endl; 
												else pt_CountDouble_DATA17[a-1][4]->Fill(lepton_pt[h][a]); //std::cout<<"ABC   "<<std::endl; 
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][12]->Fill(lepton_pt[h][a]); //std::cout<<" ABDE  "<<std::endl; 
												else pt_CountDouble_DATA17[a-1][5]->Fill(lepton_pt[h][a]); //std::cout<<"ABD   "<<std::endl; 
										}
										else if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][6]->Fill(lepton_pt[h][a]); //std::cout<<" ABE  "<<std::endl; 
										else pt_CountDouble_DATA17[a-1][0]->Fill(lepton_pt[h][a]); //std::cout<<" AB  "<<std::endl; 

									}
									else{
										if(lepton_pt[h][a] == lepton_pt[h][c]){
											if(lepton_pt[h][a] == lepton_pt[h][d]){
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][13]->Fill(lepton_pt[h][a]); //std::cout<<" ACDE  "<<std::endl; 
												else pt_CountDouble_DATA17[a-1][7]->Fill(lepton_pt[h][a]); //std::cout<<"  ACD "<<std::endl; 
											}
											else{
												if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][8]->Fill(lepton_pt[h][a]); //std::cout<<" ACE  "<<std::endl; 
												else pt_CountDouble_DATA17[a-1][1]->Fill(lepton_pt[h][a]); //std::cout<<" AC  "<<std::endl;
											}
										}
										else if(lepton_pt[h][a] == lepton_pt[h][d]){
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][9]->Fill(lepton_pt[h][a]); //std::cout<<"ADE"<<std::endl; 
											else pt_CountDouble_DATA17[a-1][2]->Fill(lepton_pt[h][a]); //std::cout<<" AD  "<<std::endl; 
										}
										else 
											if(lepton_pt[h][a] == lepton_pt[h][e]) pt_CountDouble_DATA17[a-1][3]->Fill(lepton_pt[h][a]); //std::cout<<" AE  "<<std::endl; 
									}
								}
// 								else std::cout<<" no tunep"<<std::endl;
								swap(a,b);
								swap(b,c);	
								swap(c,d);	
								swap(d,e);	
							}
							///////////// MULTIPLE COUNTING ////////////
							////////////////////////////////////////////
							////////////////////////////////////////////

	     					for(int i = 1; i < 7; i++){ // for on different reconstruction
   	 						    	 						
    	 						if(lep_pt[h] == lepton_pt[h][i]){
    	 						
    	 							bool multiple_reconstruction = false;	
    	 							for(int z = 1; z < 7; z++){
    	 								if(z == i) continue;
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
//      								if(i != 7) count_assignment_DATA17_till400[i][j-14]++;
	     							pt_DATA17[i]->Fill(lepton_pt[h][i]);
	     							eta_DATA17[i]->Fill(fabs(lepton_eta[h][i]));
	     							phi_DATA17[i]->Fill(lepton_phi[h][i]);
    	 							pt_vs_eta_DATA17[i]->Fill(lep_eta[h], lep_pt[h]);
    	 							pt_vs_phi_DATA17[i]->Fill(lep_phi[h], lep_pt[h]);
    	 							eta_vs_phi_DATA17[i]->Fill(lep_phi[h], lep_eta[h]);
    	 							pt_vs_met_DATA17[i]->Fill(met_pt, lep_pt[h]);

    	 							pt_DATA17[7]->Fill(lepton_pt[h][7]);
	     							eta_DATA17[7]->Fill(fabs(lepton_eta[h][7]));
     								phi_DATA17[7]->Fill(lepton_phi[h][7]);
   	 								pt_vs_eta_DATA17[7]->Fill(lep_eta[h], lep_pt[h]);
   	 								pt_vs_phi_DATA17[7]->Fill(lep_phi[h], lep_pt[h]);
   		 							eta_vs_phi_DATA17[7]->Fill(lep_phi[h], lep_eta[h]);
	   	 							pt_vs_met_DATA17[7]->Fill(met_pt, lep_pt[h]);

		     						if(lepton_pt[h][i] < 400){
	     								count_assignment_DATA17_till400[i]++;
	    	 							count_assignment_DATA17_till400[7]++;
		     							eta_DATA17_till400[i]->Fill(fabs(lepton_eta[h][i]));
		     							phi_DATA17_till400[i]->Fill(lepton_phi[h][i]);
		    	 						eta_DATA17_till400[7]->Fill(fabs(lepton_eta[h][7]));
	    	 							phi_DATA17_till400[7]->Fill(lepton_phi[h][7]);

		     						}
	    	 						else if(lepton_pt[h][i] > 400 && lepton_pt[h][i] < 600){
	     								count_assignment_DATA17_400to600[i]++;
	    	 							count_assignment_DATA17_400to600[7]++;
	     								eta_DATA17_400to600[i]->Fill(fabs(lepton_eta[h][i]));
	     								phi_DATA17_400to600[i]->Fill(lepton_phi[h][i]);
	     								eta_DATA17_400to600[7]->Fill(fabs(lepton_eta[h][7]));
	    	 							phi_DATA17_400to600[7]->Fill(lepton_phi[h][7]);
	     							}

	    	 						else{
	     								count_assignment_DATA17_above600[i]++;
	    	 							count_assignment_DATA17_above600[7]++;
	     								eta_DATA17_above600[i]->Fill(fabs(lepton_eta[h][i]));
	     								phi_DATA17_above600[i]->Fill(lepton_phi[h][i]);
	     								eta_DATA17_above600[7]->Fill(fabs(lepton_eta[h][7]));
	    	 							phi_DATA17_above600[7]->Fill(lepton_phi[h][7]);
	     							}
//     	 							break; // to take only one reconstruction in case of double counting --- DYT chosen because it is the first
    	 						}
     						}
    /*	 					if(lep_pt[h] == lepton_pt[h][7]){
//     	 						count_assignment_DATA16_till400[7]++;
     							pt_DATA17[7]->Fill(lepton_pt[h][7]);
     							eta_DATA17[7]->Fill(lepton_eta[h][7]);
     							phi_DATA17[7]->Fill(lepton_phi[h][7]);
     							pt_DATA17_clear[7]->Fill(lepton_pt[h][7]); //for stack plot
     							eta_DATA17_clear[7]->Fill(lepton_eta[h][7]); //for stack plot
     							phi_DATA17_clear[7]->Fill(lepton_phi[h][7]); //for stack plot
   	 							pt_vs_eta_DATA17[7]->Fill(lep_eta[h], lep_pt[h]);
   	 							pt_vs_phi_DATA17[7]->Fill(lep_phi[h], lep_pt[h]);
   	 							eta_vs_phi_DATA17[7]->Fill(lep_phi[h], lep_eta[h]);
   	 							pt_vs_met_DATA17[7]->Fill(met_pt, lep_pt[h]);


		     					if(lepton_pt[h][7] < 600){
	    	 						count_assignment_DATA17_till400[7]++;
		     						eta_DATA17_till400[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA17_till400[7]->Fill(lepton_phi[h][7]);
	     							eta_DATA17_till400_clear[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA17_till400_clear[7]->Fill(lepton_phi[h][7]);

		     					}
	    	 					else{
	    	 						count_assignment_DATA17_above600[7]++;
	     							eta_DATA17_above600[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA17_above600[7]->Fill(lepton_phi[h][7]);
     								eta_DATA17_above600_clear[7]->Fill(lepton_eta[h][7]);
     								phi_DATA17_above600_clear[7]->Fill(lepton_phi[h][7]);

	     						}
     						} */


//      						for(int i =0; i < 8; i++)
//      							std::cout<<reconstruction[i].Data()<<"\t"<<count_assignment_DATA17_till400[i][j-30]<<std::endl;
//      						std::cout<<"-"<<std::endl;
     					}
     					else NotGoodQualityMuon_DATA17++;
     				} // for on muons

	       }//end of condition on event

		 }// end loop p



     std::cout<<" Till 400 GeV "<<std::endl;
     int TOT_DATA16 = 0;
     int TOT_DATA17 = 0;
     for(int i = 1; i < 7; i++){
     	TOT_DATA16+=count_assignment_DATA16_till400[i];
     	TOT_DATA17+=count_assignment_DATA17_till400[i];
     }
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("  DATA16         --->  %15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA16_till400[1], count_assignment_DATA16_till400[2], count_assignment_DATA16_till400[3], count_assignment_DATA16_till400[4], count_assignment_DATA16_till400[5], count_assignment_DATA16_till400[7]);
 	 printf("  DATA17         --->  %15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA17_till400[1], count_assignment_DATA17_till400[2], count_assignment_DATA17_till400[3], count_assignment_DATA17_till400[4], count_assignment_DATA17_till400[5], count_assignment_DATA17_till400[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf(" DATA16 (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA16_till400[1]/TOT_DATA16, (float)count_assignment_DATA16_till400[2]/TOT_DATA16, (float)count_assignment_DATA16_till400[3]/TOT_DATA16, (float)count_assignment_DATA16_till400[4]/TOT_DATA16, (float)count_assignment_DATA16_till400[5]/TOT_DATA16, (float)count_assignment_DATA16_till400[7]/TOT_DATA16);
 	 printf(" DATA17 (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA17_till400[1]/TOT_DATA17, (float)count_assignment_DATA17_till400[2]/TOT_DATA17, (float)count_assignment_DATA17_till400[3]/TOT_DATA17, (float)count_assignment_DATA17_till400[4]/TOT_DATA17, (float)count_assignment_DATA17_till400[5]/TOT_DATA17, (float)count_assignment_DATA17_till400[7]/TOT_DATA17);


     std::cout<<"Between 400 and 600 GeV "<<std::endl;

     TOT_DATA16 = 0;
     TOT_DATA17 = 0;
     for(int i = 1; i < 7; i++){
     	TOT_DATA16+=count_assignment_DATA16_400to600[i];
     	TOT_DATA17+=count_assignment_DATA17_400to600[i];
     }
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("  DATA16         --->  %15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA16_400to600[1], count_assignment_DATA16_400to600[2], count_assignment_DATA16_400to600[3], count_assignment_DATA16_400to600[4], count_assignment_DATA16_400to600[5], count_assignment_DATA16_400to600[7]);
 	 printf("  DATA17         --->  %15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA17_400to600[1], count_assignment_DATA17_400to600[2], count_assignment_DATA17_400to600[3], count_assignment_DATA17_400to600[4], count_assignment_DATA17_400to600[5], count_assignment_DATA17_400to600[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf(" DATA16 (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA16_400to600[1]/TOT_DATA16, (float)count_assignment_DATA16_400to600[2]/TOT_DATA16, (float)count_assignment_DATA16_400to600[3]/TOT_DATA16, (float)count_assignment_DATA16_400to600[4]/TOT_DATA16, (float)count_assignment_DATA16_400to600[5]/TOT_DATA16, (float)count_assignment_DATA16_400to600[7]/TOT_DATA16);
 	 printf(" DATA17 (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA17_400to600[1]/TOT_DATA17, (float)count_assignment_DATA17_400to600[2]/TOT_DATA17, (float)count_assignment_DATA17_400to600[3]/TOT_DATA17, (float)count_assignment_DATA17_400to600[4]/TOT_DATA17, (float)count_assignment_DATA17_400to600[5]/TOT_DATA17, (float)count_assignment_DATA17_400to600[7]/TOT_DATA17); 	
 	


     std::cout<<"Above 600 GeV "<<std::endl;
     TOT_DATA16 = 0;
     TOT_DATA17 = 0;
     for(int i = 1; i < 7; i++){
     	TOT_DATA16+=count_assignment_DATA16_above600[i];
     	TOT_DATA17+=count_assignment_DATA17_above600[i];
     }
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("  DATA16         --->  %15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA16_above600[1], count_assignment_DATA16_above600[2], count_assignment_DATA16_above600[3], count_assignment_DATA16_above600[4], count_assignment_DATA16_above600[5], count_assignment_DATA16_above600[7]);
 	 printf("  DATA17         --->  %15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA17_above600[1], count_assignment_DATA17_above600[2], count_assignment_DATA17_above600[3], count_assignment_DATA17_above600[4], count_assignment_DATA17_above600[5], count_assignment_DATA17_above600[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf(" DATA16 (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA16_above600[1]/TOT_DATA16, (float)count_assignment_DATA16_above600[2]/TOT_DATA16, (float)count_assignment_DATA16_above600[3]/TOT_DATA16, (float)count_assignment_DATA16_above600[4]/TOT_DATA16, (float)count_assignment_DATA16_above600[5]/TOT_DATA16, (float)count_assignment_DATA16_above600[7]/TOT_DATA16);
 	 printf(" DATA17 (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA17_above600[1]/TOT_DATA17, (float)count_assignment_DATA17_above600[2]/TOT_DATA17, (float)count_assignment_DATA17_above600[3]/TOT_DATA17, (float)count_assignment_DATA17_above600[4]/TOT_DATA17, (float)count_assignment_DATA17_above600[5]/TOT_DATA17, (float)count_assignment_DATA17_above600[7]/TOT_DATA17);


 	 	
//  	bool save = false;
 	bool save = true;

	TString dir_save = "./PT_ASSIGNMENT_PLOT_DATAonly/";

	int a = 0;
	int b = 1;
	int c = 2;
	int d = 3;
	int e = 4;
	for(int k = 0; k < 5; k++){
		int j = k + 1;
		name_histo = Form("%s_MultipleCounting", reconstruction[j].Data());
		SaveMultipleCounting(name_histo, pt_CountDouble_DATA16[a], pt_CountDouble_DATA17[a], dir_save, a, b, c, d, e);
		swap(a,b);
		swap(b,c);	
		swap(c,d);	
		swap(d,e);	
	}	
	save_document = dir_save + "Variuos_distribution.pdf";
	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf[", 0,  0, "p_{T}");
	SalvaHisto("Dimuon Mass: TuneP", Dimuon_DATA16, Dimuon_DATA17, dir_save + "Variuos_distribution.pdf[", 1,  1, "m_{#mu#mu} [GeV]");
	SalvaHisto("Dimuon Mass: Picky", Dimuon_Picky_DATA16, Dimuon_Picky_DATA17, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
	SalvaHisto("Dimuon Mass: Dyt", Dimuon_Dyt_DATA16, Dimuon_Dyt_DATA17, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
	
	for(int i = 1; i < 8; i++){
	
		if(i != 2 && i != 4 && i != 7) continue;
	
		name_histo = Form("#eta %s: p_{T} < 400GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_DATA16_till400[i], eta_DATA17_till400[i], save_document, 0,  0, "#eta");

		name_histo = Form("#eta %s: 400 < p_{T} < 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_DATA16_400to600[i], eta_DATA17_400to600[i], save_document, 0,  0, "#eta");

		name_histo = Form("#eta %s: p_{T} > 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_DATA16_above600[i], eta_DATA17_above600[i], save_document, 0,  0, "#eta");


		name_histo = Form("#phi %s: p_{T} < 400GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_DATA16_till400[i], phi_DATA17_till400[i], save_document, 0, 0,  "#phi");

		name_histo = Form("#phi %s: 400 < p_{T} < 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_DATA16_400to600[i], phi_DATA17_400to600[i], save_document, 0,  0, "#phi");

		name_histo = Form("#phi %s: p_{T} > 600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_DATA16_above600[i], phi_DATA17_above600[i], save_document, 0,  0, "#phi");


		name_histo = Form("%s p_{T}", reconstruction[i].Data());
		SalvaHisto(name_histo, pt_DATA16[i], pt_DATA17[i], save_document, 0,  1, "p_{T}");

		name_histo = Form("%s #eta", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_DATA16[i], eta_DATA17[i], save_document, 0,  0, "#eta");

		name_histo = Form("%s #phi", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_DATA16[i], phi_DATA17[i], save_document, 0,  0, "#phi");


/*		name_histo = Form("%s p_{T} vs #eta: DATA16", reconstruction[i].Data());
		if(i != 7){
			pt_vs_eta_DATA16[i]->Divide(pt_vs_eta_DATA16[7]);
			pt_vs_eta_DATA16[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_eta_DATA16[i], save_document, "#eta", "p_{T}");		
		name_histo = Form("%s p_{T} vs #eta: DATA17", reconstruction[i].Data());
		if(i != 7){
			pt_vs_eta_DATA17[i]->Divide(pt_vs_eta_DATA17[7]);
			pt_vs_eta_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_eta_DATA17[i], save_document, "#eta", "p_{T}");
		name_histo = Form("%s p_{T} vs #eta: DATA17/DATA16", reconstruction[i].Data());
		pt_vs_eta_DATA17[i]->Divide(pt_vs_eta_DATA16[i]);
		pt_vs_eta_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, pt_vs_eta_DATA17[i], save_document, "#eta", "p_{T}", 0);


		name_histo = Form("%s p_{T} vs #phi: DATA16", reconstruction[i].Data());
		if(i != 7){
			pt_vs_phi_DATA16[i]->Divide(pt_vs_phi_DATA16[7]);
			pt_vs_phi_DATA16[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_phi_DATA16[i], save_document, "#phi", "p_{T}");		
		name_histo = Form("%s p_{T} vs #phi: DATA17", reconstruction[i].Data());
		if(i != 7){
			pt_vs_phi_DATA17[i]->Divide(pt_vs_phi_DATA17[7]);
			pt_vs_phi_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_phi_DATA17[i], save_document, "#phi", "p_{T}");
		name_histo = Form("%s p_{T} vs #phi: DATA17/DATA16", reconstruction[i].Data());
		pt_vs_phi_DATA17[i]->Divide(pt_vs_phi_DATA16[i]);
		pt_vs_phi_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);
		SalvaHisto(name_histo, pt_vs_phi_DATA17[i], save_document, "#phi", "p_{T}", 0);


		name_histo = Form("%s #eta vs #phi: DATA16", reconstruction[i].Data());
		if(i != 7){
			eta_vs_phi_DATA16[i]->Divide(eta_vs_phi_DATA16[7]);
			eta_vs_phi_DATA16[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, eta_vs_phi_DATA16[i], save_document, "#phi", "#eta");		
		name_histo = Form("%s #eta vs #phi: DATA17", reconstruction[i].Data());
		if(i != 7){
			eta_vs_phi_DATA17[i]->Divide(eta_vs_phi_DATA17[7]);
			eta_vs_phi_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, eta_vs_phi_DATA17[i], save_document, "#phi", "#eta");
		name_histo = Form("%s #eta vs #phi: DATA17/DATA16", reconstruction[i].Data());
		eta_vs_phi_DATA17[i]->Divide(eta_vs_phi_DATA16[i]);
		eta_vs_phi_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, eta_vs_phi_DATA17[i], save_document, "#phi", "#eta", 0);


		name_histo = Form("%s p_{T} vs met: DATA16", reconstruction[i].Data());
		if(i != 7){
			pt_vs_met_DATA16[i]->Divide(pt_vs_met_DATA16[7]);
			pt_vs_met_DATA16[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_met_DATA16[i], save_document, "met", "p_{T}");
		name_histo = Form("%s p_{T} vs met: DATA17", reconstruction[i].Data());
		if(i != 7){
			pt_vs_met_DATA17[i]->Divide(pt_vs_met_DATA17[7]);
			pt_vs_met_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_met_DATA17[i], save_document, "met", "p_{T}");
		name_histo = Form("%s p_{T} vs met: DATA17/DATA16", reconstruction[i].Data());
		pt_vs_met_DATA17[i]->Divide(pt_vs_met_DATA16[i]);
		pt_vs_met_DATA17[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, pt_vs_met_DATA17[i], save_document, "met", "p_{T}", 0);
*/
	}	

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf]", 0,  0, "p_{T}");

	TCanvas *c_Double_Count_DATA16;
	TLatex lat;
	TString vv = " == ";
	for(int i = 0; i < 6; i++){
		TString cc = "c_";
		TString nome_canvas = cc +  reconstruction[i+1].Data(); 
		c_Double_Count_DATA16 = new TCanvas(nome_canvas, nome_canvas, 1050, 750);
		c_Double_Count_DATA16->Divide(2,3);
		int z = 0;
		int h = 0;
		for(int j = 0; j < 6; j++){
			z++;
			if(i == j) continue;
			c_Double_Count_DATA16->cd(z);
			float max = Double_Count_DATA17[i][j]->GetMaximum();
			if(max < Double_Count_DATA16[i][j]->GetMaximum())
				max = Double_Count_DATA16[i][j]->GetMaximum();
			max = 1.1 * max;
			Double_Count_DATA16[i][j]->Draw();
			Double_Count_DATA16[i][j]->SetTitle("");
			Double_Count_DATA16[i][j]->SetLineColor(kBlue);
			Double_Count_DATA16[i][j]->GetYaxis()->SetRangeUser(0, max);
			gPad->Update();
			TPaveStats *s_DATA16 = (TPaveStats*)Double_Count_DATA16[i][j]->GetListOfFunctions()->FindObject("stats");
			s_DATA16->SetName("Const");
			s_DATA16->SetX1NDC(0.75);
			s_DATA16->SetY1NDC(0.52);
			s_DATA16->SetY2NDC(0.72);
			s_DATA16->SetTextColor(kBlue);
			Double_Count_DATA17[i][j]->Draw();
			Double_Count_DATA17[i][j]->SetTitle("");
			Double_Count_DATA17[i][j]->SetLineColor(kRed);
			Double_Count_DATA17[i][j]->GetYaxis()->SetRangeUser(0, max);
			gPad->Update();
			TPaveStats *st_DATA17 = (TPaveStats*)Double_Count_DATA17[i][j]->GetListOfFunctions()->FindObject("stats");
			st_DATA17->SetName("Const");
			st_DATA17->SetX1NDC(0.75);
			st_DATA17->SetY1NDC(0.75);
			st_DATA17->SetY2NDC(0.95);
			st_DATA17->SetTextColor(kRed);
// 			pt_DATA16[i]->GetYaxis()->SetTitleOffset(1.2);
// 			pt_DATA17[i]->GetYaxis()->SetTitleOffset(1.2);
	    	Double_Count_DATA16[i][j]->Draw();
    		Double_Count_DATA17[i][j]->Draw("same");	
			s_DATA16->Draw("same");
			st_DATA17->Draw("same");
			TString nome_lat = reconstruction[i+1].Data() + vv + reconstruction[z].Data();
			max = 0.9 * max / 1.1;
			if(max == 0)
				max = 0.9;
			lat.DrawLatex(850, max, nome_lat);

		}
		if(save){
			if(i==0)
		 		c_Double_Count_DATA16->Print(dir_save + "PT_DoubleCounting_distribution.pdf[");
			c_Double_Count_DATA16->Print(dir_save + "PT_DoubleCounting_distribution.pdf");
 			if(i==5)
				c_Double_Count_DATA16->Print(dir_save + "PT_DoubleCounting_distribution.pdf]");
		    c_Double_Count_DATA16->Write();
		}
	}
	
}// end of function 


void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA17[15], TString save, int a, int b, int c, int d, int e){
	TLatex lat;
	TString nome_lat[5] = {"Global", "Picky", "Tpfms", "DYT", "Tracker"};
	TString inizio = Form("%s%s.pdf[", save.Data(), name.Data());
	TString fine = Form("%s%s.pdf]", save.Data(), name.Data());
	TString mezzo = Form("%s%s.pdf", save.Data(), name.Data());
	
	float max = 0;
	
		SalvaHisto("h_blank", h_blank, h_blank, inizio,  0,  0, "p_{T}");
		SalvaHisto("Multiple Assignment", h[0], h_DATA17[0], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b]);
		SalvaHisto("Multiple Assignment", h[1], h_DATA17[1], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c]);
		SalvaHisto("Multiple Assignment", h[2], h_DATA17[2], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[3], h_DATA17[3], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[4], h_DATA17[4], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c]);
		SalvaHisto("Multiple Assignment", h[5], h_DATA17[5], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[6], h_DATA17[6], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[7], h_DATA17[7], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[8], h_DATA17[8], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[9], h_DATA17[9], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[d], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[10], h_DATA17[10], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c], nome_lat[d]);
		SalvaHisto("Multiple Assignment", h[11], h_DATA17[11], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[12], h_DATA17[12], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[d], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[13], h_DATA17[13], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[c], nome_lat[d], nome_lat[e]);
		SalvaHisto("Multiple Assignment", h[14], h_DATA17[14], mezzo, 0, 0, "p_{T}", 250, nome_lat[a], nome_lat[b], nome_lat[c], nome_lat[d], nome_lat[e]);
		SalvaHisto("h_blank", h_blank, h_blank, fine, 0, 0, "p_{T}");

}

void SalvaHisto(TString name, TH2F* h_DATA16, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic, bool logX){

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
		h_DATA16->SetStats();
	else
		h_DATA16->SetStats(0);
		
	if(logX)
		pad11->SetLogx();
		
	h_DATA16->Draw("COLZ");
// 	h_DATA16->Draw("SAME TEXT0");
	h_DATA16->SetTitle(name);
	h_DATA16->GetYaxis()->SetTitleOffset(0.65);
	h_DATA16->GetXaxis()->SetTitleOffset(0.6);
	h_DATA16->GetXaxis()->SetTitle(name_Xaxis);
	h_DATA16->GetYaxis()->SetTitle(name_Yaxis);
			
	TString save_name_pdf = name + ".pdf";
	TString save_name_png = name + ".png";
	
// 	c1->Print(save_name_pdf);
// 	c1->Print(save_name_png);
	c1->Print(save);
	
}

void SalvaHisto(TString name, TH1F* h_DATA16, TH1F* h_DATA17, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex, TString strig_1, TString strig_2, TString strig_3, TString strig_4, TString strig_5){

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

	if(h_DATA17->GetMaximum() > h_DATA16->GetMaximum())
		max = h_DATA17->GetMaximum();
	else
		max = h_DATA16->GetMaximum();
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;
	
	h_DATA16->SetLineColor(kBlue);
	h_DATA16->SetTitle(name);
	h_DATA16->Draw();
	h_DATA16->SetStats(1);
	c1->Update();

	TPaveStats * st_DATA16 = (TPaveStats *)h_DATA16->GetListOfFunctions()->FindObject("stats");
    if( st_DATA16 ){ 
		st_DATA16->SetName("Const");
		st_DATA16->SetX1NDC(0.75);
		st_DATA16->SetY1NDC(0.52);
		st_DATA16->SetY2NDC(0.72);
		st_DATA16->SetTextColor(kBlue);
    }
    else std::cout << "Null pointer to TPaveStats DATA16: " << st_DATA16 << std::endl;
  
	h_DATA17->SetLineColor(kRed);
	h_DATA17->SetTitle(name);
	h_DATA17->Draw();
// 	h_DATA17->SetStats(0);
	TH1F *ratio = (TH1F*) h_DATA17->Clone();
	h_DATA17->SetStats(1);
	c1->Update();
// 	gPad->Update();

	TPaveStats * st_DATA17 = (TPaveStats *)h_DATA17->GetListOfFunctions()->FindObject("stats");
    if( st_DATA17 ){ 
    	st_DATA17->SetTextColor(kRed); 
		st_DATA17->SetName("Const");
		st_DATA17->SetX1NDC(0.75);
		st_DATA17->SetY1NDC(0.75);
		st_DATA17->SetY2NDC(0.95);

    }
    else std::cout << "Null pointer to TPaveStats DATA17: " << st_DATA17 << std::endl;

	h_DATA16->GetYaxis()->SetTitleOffset(1.2);
	h_DATA17->GetYaxis()->SetTitleOffset(1.2);

	h_DATA16->SetMaximum(max);
	h_DATA16->SetMinimum(min);	
// 	h_DATA16->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h_DATA17->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);

		
    c1->Update();

	h_DATA16->Draw();
	h_DATA17->Draw("sameP");
	h_DATA17->SetMarkerStyle(3);
	h_DATA17->SetMarkerColor(kRed);
	h_DATA17->SetMarkerSize(0.5);
	st_DATA16->Draw("same");
// 	st_DATA17->Draw("same");
	h_DATA16->GetXaxis()->SetTitle(name_axis);
	h_DATA16->GetYaxis()->SetTitleSize(0.04);
	h_DATA16->GetYaxis()->SetLabelSize(0.04);
	h_DATA16->GetXaxis()->SetTitleSize(0.04);
	h_DATA16->GetXaxis()->SetLabelSize(0.04);
	h_DATA17->GetXaxis()->SetTitle(name_axis);
	h_DATA17->GetYaxis()->SetTitleSize(0.04);
	h_DATA17->GetYaxis()->SetLabelSize(0.04);
	h_DATA17->GetXaxis()->SetTitleSize(0.04);
	h_DATA17->GetXaxis()->SetLabelSize(0.04);
	
//     	
	TLegend *l1 = new TLegend(0.4,0.8,0.6,0.9);
	l1->AddEntry(h_DATA17, "DATA 2017", "l");
	l1->AddEntry(h_DATA16, "DATA 2016", "l");
	l1->Draw();	
	
	lat.DrawLatex(X_pos_latex, max, strig_1);
	lat.DrawLatex(X_pos_latex, 0.9*max, strig_2);
	lat.DrawLatex(X_pos_latex, 0.8*max, strig_3);
	lat.DrawLatex(X_pos_latex, 0.7*max, strig_4);
	lat.DrawLatex(X_pos_latex, 0.6*max, strig_5);
	
	c1->Update();
	c1->cd();
	
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.25);
	pad22->SetTopMargin(0.95);
	pad22->SetBottomMargin(0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();

	if(logx)
		pad22->SetLogx();
	pad22->cd();
	
	ratio->Divide(h_DATA16);
	ratio->GetYaxis()->SetRangeUser(0, 2);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("DATA17 / DATA16");
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
