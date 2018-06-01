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
void SalvaHisto(TString name, TH1F* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, TH2F* h_MC, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic = 1, bool logX = false);
void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA[15], TString save, int a, int b, int c, int d, int e);

  Int_t event;
  Int_t run;
  Int_t prev_event = -88;
//   Int_t lumi;

  float dil_mass;
//   float dil_pt;
  float cos_angle;
//   float vertex_chi2;
//   int dil_chosen;
//   float gen_dil_mass;
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
//   float gen_lep_pt[2];
  float lep_triggerMatchPt[2];
  float vertex_chi2;
  
  float met_pt;

	TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);
         
	float LUMINOSITY_2016 = 36235.493;
	float LUMINOSITY_2017_reReco = 41903.837;
	
	float weight[39] = {0};
	
	int count_event_MC[6][39] = {0};
	int count_event_DATA[6] = {0};
		
	float lepton_pt[2][8] = {0};
	float lepton_eta[2][8] = {0};
	float lepton_phi[2][8] = {0};
	int numberOfValidMuonHits[2][8];
	int NotGoodQualityMuon_MC = 0;
	int NotGoodQualityMuon_DATA = 0;
	int count_assignment_MC[8];
	int count_assignment_DATA[8];

	int tot[8] = {0};
	int TOT_MC = 0, TOT_DATA = 0;
	
	TString reconstruction[8] = {"Selected", "global", "picky", "tpfms", "dyt", "tracker", "std", "tuneP"};
	TString pt_name[2] = {"2017_RR_pt_", "2016_pt_"};
	TString eta_name[2] = {"2017_RR_eta_", "2016_eta_"};
	TString phi_name[2] = {"2017_RR_phi_", "2016_phi_"};
	TString pt_vs_eta_name[2] = {"2017_RR_pt_vs_eta_", "2016_ptVSeta_"};//2016_pt_vs_eta_"};
	TString pt_vs_phi_name[2] = {"2017_RR_pt_vs_phi_", "2016_ptVSphi_"};//2016_pt_vs_phi_"};
	TString eta_vs_phi_name[2] = {"2017_RR_eta_vs_phi_", "2016_etaVSphi_"};//2016_eta_vs_phi_"};
	TString pt_vs_met_name[2] = {"2017_RR_pt_vs_met_", "2016_ptVSmet_"};//2016_pt_vs_met_"};
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
	TH1F* eta_MC_till600[8];
	TH1F* eta_DATA_till600[8];
	TH1F* eta_MC_above600[8];
	TH1F* eta_DATA_above600[8];
	TH1F* phi_MC_till600[8];
	TH1F* phi_DATA_till600[8];
	TH1F* phi_MC_above600[8];
	TH1F* phi_DATA_above600[8];

	TH1D* h_count_assignment_MC[8];
	TH1D* h_count_assignment_DATA[8];
	
	TH1F* pt_CountDouble_MC[5][15];
	TH1F* pt_CountDouble_DATA[5][15];
	
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
   
    int count_double_MC[8][8] = {0};
    int count_double_DATA[8][8] = {0};
    
void Pt_Assignment_DATAonly(){

// 	gROOT->LoadMacro("/afs/cern.ch/work/f/ferrico/private/Ratioplot.C+");

    gROOT->Reset();
    gROOT->SetBatch();

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
	
	TH1F *Dimuon_MC = new TH1F("Dimuon mass: DATA 2017 RR", "Dimuon mass: DATA 2017 RR", NMBINS, logMbins);
	THStack *Dimuon_MC_Stack = new THStack("Dimuon mass: MC - stack plot", "Dimuon mass: MC - stack plot");
	TH1F *Dimuon_DATA = new TH1F("Dimuon mass: DATA", "Dimuon mass: DATA", NMBINS, logMbins);
	
     
	for(int i = 0; i < 39; i++){
// 		weight[i] = LUMINOSITY * sigma[i] / events[i];
// 		weight[i] *= Z_peak;
		weight[i] = LUMINOSITY_2016 / LUMINOSITY_2017_reReco;
	}
	
	for(int a = 0; a < 5; a++){
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
		eta_MC_till600[i] = new TH1F(eta_name[0] + reconstruction[i] + "_till600", eta_name[0] + reconstruction[i] + "_till600", binnum_eta, ETA_BINS);
		eta_DATA_till600[i] = new TH1F(eta_name[1] + reconstruction[i] + "_till600", eta_name[1] + reconstruction[i] + "_till600",binnum_eta, ETA_BINS);
		eta_MC_above600[i] = new TH1F(eta_name[0] + reconstruction[i] + "_above600", eta_name[0] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);
		eta_DATA_above600[i] = new TH1F(eta_name[1] + reconstruction[i] + "_above600", eta_name[1] + reconstruction[i] + "_above600", binnum_eta, ETA_BINS);
		phi_MC_till600[i] = new TH1F(phi_name[0] + reconstruction[i] + "_till600", phi_name[0] + reconstruction[i] + "_till600", binnum_phi, PHI_BINS);
		phi_DATA_till600[i] = new TH1F(phi_name[1] + reconstruction[i] + "_till600", phi_name[1] + reconstruction[i] + "_till600",binnum_phi, PHI_BINS);
		phi_MC_above600[i] = new TH1F(phi_name[0] + reconstruction[i] + "_above600", phi_name[0] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);
		phi_DATA_above600[i] = new TH1F(phi_name[1] + reconstruction[i] + "_above600", phi_name[1] + reconstruction[i] + "_above600", binnum_phi, PHI_BINS);

		pt_vs_met_MC[i] = new TH2F(pt_vs_met_name[0] + reconstruction[i], pt_vs_met_name[0] + reconstruction[i], binnum_met, MET_BINS, binnum_pt, PT_BINS);
		pt_vs_met_DATA[i] = new TH2F(pt_vs_met_name[1] + reconstruction[i], pt_vs_met_name[1] + reconstruction[i], binnum_met, MET_BINS, binnum_pt, PT_BINS);



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


		eta_MC_till600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA_till600[i]->GetXaxis()->SetTitle("#eta");
		eta_MC_above600[i]->GetXaxis()->SetTitle("#eta");
		eta_DATA_above600[i]->GetXaxis()->SetTitle("#eta");
		phi_MC_till600[i]->GetXaxis()->SetTitle("#phi");
		phi_DATA_till600[i]->GetXaxis()->SetTitle("#phi");
		phi_MC_above600[i]->GetXaxis()->SetTitle("#phi");
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
	
// 	TH2F* std_vs_tuneP_MC = new TH2F("Std vs TuneP: MC", "Std vs TuneP: MC", 50, 200, 205, 50, 200, 205);//binnum_pt, PT_BINS, binnum_pt, PT_BINS);
// 	TH2F* std_vs_tuneP_DATA = new TH2F("Std vs TuneP: DATA", "Std vs TuneP: DATA", 1800, 200, 1800, 1800, 200, 2000);


     for(int j = 38; j <39; j++){
     
		TH1F *Dimuon_MC_clear = new TH1F("Dimuon mass: MC", "Dimuon mass: MC", NMBINS, logMbins);
		for(int i = 0; i < 8; i++){
			pt_MC_clear[i] = new TH1F(pt_name[0] + reconstruction[i], pt_name[0] + reconstruction[i], binnum_pt, PT_BINS);
			eta_MC_clear[i] = new TH1F(eta_name[0] + reconstruction[i], eta_name[0] + reconstruction[i], binnum_eta, ETA_BINS);
			phi_MC_clear[i] = new TH1F(phi_name[0] + reconstruction[i], phi_name[0] + reconstruction[i], binnum_phi, PHI_BINS);
			pt_MC_clear[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
			eta_MC_clear[i]->GetXaxis()->SetTitle("#eta");
			phi_MC_clear[i]->GetXaxis()->SetTitle("#phi");
		}

     	TChain *treeMC = new TChain("SimpleNtupler/t");
//      	treeMC->Add("./mc_Pt_assignment/ana_datamc_"+samples[j]+".root");

     	treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2018/DATA/ana_datamc_data.root");     	
 
      	treeMC->SetBranchAddress("event", &event);
      	treeMC->SetBranchAddress("run", &run);
//       	treeMC->SetBranchAddress("lumi", &lumi);    	
    	treeMC->SetBranchAddress("dil_mass",&dil_mass);
//      treeMC->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
//      treeMC->SetBranchAddress("dil_pt",&dil_pt);
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
	     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);	
    	 treeMC->SetBranchAddress("lep_picky_numberOfValidMuonHits",lep_picky_numberOfValidMuonHits);
	     treeMC->SetBranchAddress("lep_dyt_numberOfValidMuonHits",lep_dyt_numberOfValidMuonHits);
    	 treeMC->SetBranchAddress("lep_tpfms_numberOfValidMuonHits",lep_tpfms_numberOfValidMuonHits);
	     treeMC->SetBranchAddress("lep_stanAlone_numberOfValidMuonHits",lep_stanAlone_numberOfValidMuonHits);

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

	     Long64_t nentries = treeMC->GetEntries();
	     
		printf("opening... 2017 DATA --- %lld --- %.5f\n", nentries, weight[j]);    
	
		int bhu_mc = 0;

    	 for(int p=0; p<nentries; p++){
//     	 for(int p=0; p<10000; p++){

			if(p % 100000 == 0) std::cout<<p<<" su "<<nentries<<std::endl;		
  
    	 	treeMC->GetEntry(p);
			
//     	 	std::cout<<lep_std_pt[0]<<"\t"<<lep_std_pt[1]<<std::endl;

	     	if (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
	     		lep_pt[0]>53. && lep_pt[1]>53. && 
// 	     		lep_pt[0]>200. && lep_pt[1]>200. && 
    	 		lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
     			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
     			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
	     		lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
    	 		lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
    	 		(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
    	 		
    	 		dil_mass > 120 &&

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
     					){  Dimuon_MC->Fill(dil_mass,  weight[j]); Dimuon_MC_clear->Fill(dil_mass,  weight[j]); }

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
    	 									continue;
    	 								}
    	 							}
    	 							if(multiple_reconstruction) continue; //to remove multiple reconstructions
    	 						
									// Select Picky only if Picky = TuneP    	 						
    	 							if(i == 2 && (lepton_pt[h][i] == lepton_pt[h][3] || lepton_pt[h][i] == lepton_pt[h][5])) continue;
									// Select DYT only if DYT = TuneP    	 						
									if(i == 4 && lepton_pt[h][i] == lepton_pt[h][5]) continue;
//      								if(i != 7) count_assignment_MC[i]++;
	     							pt_MC[i]->Fill(lepton_pt[h][i], weight[j]);
	     							eta_MC[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     							phi_MC[i]->Fill(lepton_phi[h][i], weight[j]);
    	 							pt_vs_eta_MC[i]->Fill(lep_eta[h], lep_pt[h], weight[j]);
    	 							pt_vs_phi_MC[i]->Fill(lep_phi[h], lep_pt[h], weight[j]);
    	 							eta_vs_phi_MC[i]->Fill(lep_phi[h], lep_eta[h], weight[j]);
    	 							pt_vs_met_MC[i]->Fill(met_pt, lep_pt[h], weight[j]);

		     						if(lepton_pt[h][i] < 600){
		     							eta_MC_till600[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
		     							phi_MC_till600[i]->Fill(lepton_phi[h][i], weight[j]);
		     						}
	    	 						else{
	     								if(i != 7) count_assignment_MC[i]++;
	     								eta_MC_above600[i]->Fill(fabs(lepton_eta[h][i]), weight[j]);
	     								phi_MC_above600[i]->Fill(lepton_phi[h][i], weight[j]);
	     							}

    	 							break; // to take only one reconstruction in case of double counting --- DYT chosen because it is the first
    	 						}
     						}
    	 					if(lep_pt[h] == lepton_pt[h][7]){
//     	 						count_assignment_MC[7]++;
     							pt_MC[7]->Fill(lepton_pt[h][7], weight[j]);
     							eta_MC[7]->Fill(lepton_eta[h][7], weight[j]);
     							phi_MC[7]->Fill(lepton_phi[h][7], weight[j]);
   	 							pt_vs_eta_MC[7]->Fill(lep_eta[h], lep_pt[h], weight[j]);
   	 							pt_vs_phi_MC[7]->Fill(lep_phi[h], lep_pt[h], weight[j]);
   	 							eta_vs_phi_MC[7]->Fill(lep_phi[h], lep_eta[h], weight[j]);
   	 							pt_vs_met_MC[7]->Fill(met_pt, lep_pt[h], weight[j]);


		     					if(lepton_pt[h][7] < 600){
		     						eta_MC_till600[7]->Fill(lepton_eta[h][7], weight[j]);
	     							phi_MC_till600[7]->Fill(lepton_phi[h][7], weight[j]);
		     					}
	    	 					else{
	    	 						count_assignment_MC[7]++;
	     							eta_MC_above600[7]->Fill(lepton_eta[h][7], weight[j]);
	     							phi_MC_above600[7]->Fill(lepton_phi[h][7], weight[j]);
	     						}
     						}

//      						for(int i =0; i < 8; i++)
//      							std::cout<<reconstruction[i].Data()<<"\t"<<count_assignment_MC[i][j-30]<<std::endl;
//      						std::cout<<"-"<<std::endl;
     					}
     					else NotGoodQualityMuon_MC++;
     				} // for on muons
	       }//end of condition on event

		 }// end loop p
		 
     }// end loop on MC
    
    
////////      DATA    ///////

	for(int i = 0; i < 7; i++)
		for(int j = 0; j < 2; j++)
			lepton_pt[j][i] = {0};
		
    TChain *treeDATA = new TChain("SimpleNtupler/t");
//      	treeDATA->Add("./mc_Pt_assignment/ana_datamc_"+samples[j]+".root");
    treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/DATA_Pt_Assignment/ana_datamc_data_Pt_Assignment.root");

     treeDATA->SetBranchAddress("event",&event);
     treeDATA->SetBranchAddress("run",&run);
//      treeDATA->SetBranchAddress("lumi",&lumi);
     	
     treeDATA->SetBranchAddress("dil_mass",&dil_mass);
//      treeDATA->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
//      treeDATA->SetBranchAddress("dil_pt",&dil_pt);
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

	Long64_t nentries = treeDATA->GetEntries();
	
	printf("opening... 2016 DATA --- %lld\n", nentries);    
		int bhu_data = 0;
	
	for(int p=0; p<nentries; p++){
// 	for(int p=0; p<10000; p++){

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
    	 		
    	 		dil_mass > 120 &&

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
     					) Dimuon_DATA->Fill(dil_mass);

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
    	 									continue;
    	 								}
    	 							}
    	 							if(multiple_reconstruction) continue; //to remove multiple reconstructions
    	 						
									// Select Picky only if Picky = TuneP    	 						
    	 							if(i == 2 && (lepton_pt[h][i] == lepton_pt[h][3] || lepton_pt[h][i] == lepton_pt[h][5])) continue;
									// Select DYT only if DYT = TuneP    	 						
									if(i == 4 && lepton_pt[h][i] == lepton_pt[h][5]) continue;
//      								if(i != 7) count_assignment_DATA[i]++;
	     							pt_DATA[i]->Fill(lepton_pt[h][i]);
	     							eta_DATA[i]->Fill(fabs(lepton_eta[h][i]));
	     							phi_DATA[i]->Fill(lepton_phi[h][i]);
    	 							pt_vs_eta_DATA[i]->Fill(lep_eta[h], lep_pt[h]);
    	 							pt_vs_phi_DATA[i]->Fill(lep_phi[h], lep_pt[h]);
    	 							eta_vs_phi_DATA[i]->Fill(lep_phi[h], lep_eta[h]);
    	 							pt_vs_met_DATA[i]->Fill(met_pt, lep_pt[h]);

		    	 					if(lepton_pt[h][i] < 600){
		     							eta_DATA_till600[i]->Fill(fabs(lepton_eta[h][i]));
		     							phi_DATA_till600[i]->Fill(lepton_phi[h][i]);
		     						}
	    	 						else{
	     								if(i != 7) count_assignment_DATA[i]++;
	     								eta_DATA_above600[i]->Fill(fabs(lepton_eta[h][i]));
	     								phi_DATA_above600[i]->Fill(lepton_phi[h][i]);
	     							}

    	 							break; // to take only one reconstruction in case of double counting --- DYT chosen because it is the first
    	 						}
     						}
    	 					if(lep_pt[h] == lepton_pt[h][7]){ 
//     	 						count_assignment_DATA[7]++;
    	 						pt_DATA[7]->Fill(lepton_pt[h][7]);
    	 						eta_DATA[7]->Fill(lepton_eta[h][7]);
    	 						phi_DATA[7]->Fill(lepton_phi[h][7]);
   	 							pt_vs_eta_DATA[7]->Fill(lep_eta[h], lep_pt[h]);
   	 							pt_vs_phi_DATA[7]->Fill(lep_phi[h], lep_pt[h]);
   	 							eta_vs_phi_DATA[7]->Fill(lep_phi[h], lep_eta[h]);
   	 							pt_vs_met_DATA[7]->Fill(met_pt, lep_pt[h]);

		     					if(lepton_pt[h][7] < 600){
		     						eta_DATA_till600[7]->Fill(lepton_eta[h][7]);
		     						phi_DATA_till600[7]->Fill(lepton_phi[h][7]);
		     					}
	    	 					else{
	    	 						count_assignment_DATA[7]++;
	     							eta_DATA_above600[7]->Fill(lepton_eta[h][7]);
	     							phi_DATA_above600[7]->Fill(lepton_phi[h][7]);
	     						}

    	 					}


//      						for(int i =0; i < 8; i++)
//      							std::cout<<reconstruction[i].Data()<<"\t"<<count_assignment_DATA[i][j-30]<<std::endl;
//      						std::cout<<"-"<<std::endl;
     					}
     					else NotGoodQualityMuon_DATA++;
     				} // for on muons

	       }//end of condition on event

		 }// end loop p



     printf("                --> %15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t\t|%15s\n", reconstruction[0].Data(), reconstruction[1].Data(), reconstruction[2].Data(), reconstruction[3].Data(), reconstruction[4].Data(), reconstruction[5].Data(), reconstruction[6].Data(), reconstruction[7].Data()); 


     for(int i = 1; i < 7; i++){
     	TOT_MC+=count_assignment_MC[i];
     	TOT_DATA+=count_assignment_DATA[i];
     }

 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf("  DATA 2017 (94X)         --->  %15d\t%15d\t%15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_MC[0], count_assignment_MC[1], count_assignment_MC[2], count_assignment_MC[3], count_assignment_MC[4], count_assignment_MC[5], count_assignment_MC[6], count_assignment_MC[7]);
 	 printf("  DATA 2016               --->  %15d\t%15d\t%15d\t%15d\t%15d\t%15d\t%15d\t\t|%15d\n", count_assignment_DATA[0], count_assignment_DATA[1], count_assignment_DATA[2], count_assignment_DATA[3], count_assignment_DATA[4], count_assignment_DATA[5], count_assignment_DATA[6], count_assignment_DATA[7]);
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
 	 printf("  2017 TOTAL = %10d\2016 TOTAL = %10d\n", TOT_MC, TOT_DATA); 	 
 	 printf(" ----------------------------------------------------------------------------------------------------------------------------------- \n");
     printf(" DATA 2017 (94X) (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_MC[0]/TOT_MC, (float)count_assignment_MC[1]/TOT_MC, (float)count_assignment_MC[2]/TOT_MC, (float)count_assignment_MC[3]/TOT_MC, (float)count_assignment_MC[4]/TOT_MC, (float)count_assignment_MC[5]/TOT_MC, (float)count_assignment_MC[6]/TOT_MC, (float)count_assignment_MC[7]/TOT_MC);
 	 printf(" DATA 2016       (%%)      ---> %15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t%15.3f\t\t|%15.3f\n", (float)count_assignment_DATA[0]/TOT_DATA, (float)count_assignment_DATA[1]/TOT_DATA, (float)count_assignment_DATA[2]/TOT_DATA, (float)count_assignment_DATA[3]/TOT_DATA, (float)count_assignment_DATA[4]/TOT_DATA, (float)count_assignment_DATA[5]/TOT_DATA, (float)count_assignment_DATA[6]/TOT_DATA, (float)count_assignment_DATA[7]/TOT_DATA);

//      std::cout<<"No good muon MC   = "<<NotGoodQualityMuon_MC<<std::endl;
//      std::cout<<"No good muon DATA = "<<NotGoodQualityMuon_DATA<<std::endl;
     
 	bool save = false;
//  	bool save = true;

	TString dir_save = "./PT_ASSIGNMENT_PLOT_DATAonly/";

	int a = 0;
	int b = 1;
	int c = 2;
	int d = 3;
	int e = 4;
	for(int k = 0; k < 5; k++){
		int j = k + 1;
		name_histo = Form("%s_MultipleCounting", reconstruction[j].Data());
		SaveMultipleCounting(name_histo, pt_CountDouble_MC[a], pt_CountDouble_DATA[a], dir_save, a, b, c, d, e);
		swap(a,b);
		swap(b,c);	
		swap(c,d);	
		swap(d,e);	
	}
	
	save_document = dir_save + "Variuos_distribution.pdf";
	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf[", 0,  0, "p_{T}");
	SalvaHisto("Dimuon Mass", Dimuon_MC, Dimuon_DATA, dir_save + "Variuos_distribution.pdf[", 1,  1, "m_{#mu#mu} [GeV]");
	for(int i = 1; i < 8; i++){
	
		name_histo = Form("Eta_%s_till600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC_till600[i], eta_DATA_till600[i], save_document, 0,  0, "#eta");
		name_histo = Form("Eta_%s_above600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC_above600[i], eta_DATA_above600[i], save_document, 0,  0, "#eta");
		name_histo = Form("Phi_%s_till600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC_till600[i], phi_DATA_till600[i], save_document, 0, 0,  "#phi");
		name_histo = Form("Phi_%s_above600GeV", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC_above600[i], phi_DATA_above600[i], save_document, 0,  0, "#phi");

		name_histo = Form("%s p_{T}", reconstruction[i].Data());
		SalvaHisto(name_histo, pt_MC[i], pt_DATA[i], save_document, 0,  1, "p_{T}");
		name_histo = Form("%s #eta", reconstruction[i].Data());
		SalvaHisto(name_histo, eta_MC[i], eta_DATA[i], save_document, 0,  0, "#eta");
		name_histo = Form("%s #phi", reconstruction[i].Data());
		SalvaHisto(name_histo, phi_MC[i], phi_DATA[i], save_document, 0,  0, "#phi");

		name_histo = Form("%s p_{T} vs #eta: DATA 2017 RR", reconstruction[i].Data());
		if(i != 7){
			pt_vs_eta_MC[i]->Divide(pt_vs_eta_MC[7]);
			pt_vs_eta_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_eta_MC[i], save_document, "#eta", "p_{T}");		
		name_histo = Form("%s p_{T} vs #eta: DATA 2016", reconstruction[i].Data());
		if(i != 7){
			pt_vs_eta_DATA[i]->Divide(pt_vs_eta_DATA[7]);
			pt_vs_eta_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_eta_DATA[i], save_document, "#eta", "p_{T}");
		name_histo = Form("%s p_{T} vs #eta: DATA 2016 / DATA 2017 RR", reconstruction[i].Data());
		pt_vs_eta_DATA[i]->Divide(pt_vs_eta_MC[i]);
		pt_vs_eta_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, pt_vs_eta_DATA[i], save_document, "#eta", "p_{T}", 0);


		name_histo = Form("%s p_{T} vs #phi: DATA 2017 RR", reconstruction[i].Data());
		if(i != 7){
			pt_vs_phi_MC[i]->Divide(pt_vs_phi_MC[7]);
			pt_vs_phi_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_phi_MC[i], save_document, "#phi", "p_{T}");		
		name_histo = Form("%s p_{T} vs #phi: DATA 2016", reconstruction[i].Data());
		if(i != 7){
			pt_vs_phi_DATA[i]->Divide(pt_vs_phi_DATA[7]);
			pt_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_phi_DATA[i], save_document, "#phi", "p_{T}");
		name_histo = Form("%s p_{T} vs #phi: DATA 2016 / DATA 2017 RR", reconstruction[i].Data());
		pt_vs_phi_DATA[i]->Divide(pt_vs_phi_MC[i]);
		pt_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);
		SalvaHisto(name_histo, pt_vs_phi_DATA[i], save_document, "#phi", "p_{T}", 0);


		name_histo = Form("%s #eta vs #phi: DATA 2017 RR", reconstruction[i].Data());
		if(i != 7){
			eta_vs_phi_MC[i]->Divide(eta_vs_phi_MC[7]);
			eta_vs_phi_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, eta_vs_phi_MC[i], save_document, "#phi", "#eta");		
		name_histo = Form("%s #eta vs #phi: DATA 2016", reconstruction[i].Data());
		if(i != 7){
			eta_vs_phi_DATA[i]->Divide(eta_vs_phi_DATA[7]);
			eta_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, eta_vs_phi_DATA[i], save_document, "#phi", "#eta");
		name_histo = Form("%s #eta vs #phi: DATA 2016 / DATA 2017 RR", reconstruction[i].Data());
		eta_vs_phi_DATA[i]->Divide(eta_vs_phi_MC[i]);
		eta_vs_phi_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, eta_vs_phi_DATA[i], save_document, "#phi", "#eta", 0);


		name_histo = Form("%s p_{T} vs met: DATA 2017 RR", reconstruction[i].Data());
		if(i != 7){
			pt_vs_met_MC[i]->Divide(pt_vs_met_MC[7]);
			pt_vs_met_MC[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_met_MC[i], save_document, "met", "p_{T}");
		name_histo = Form("%s p_{T} vs met: DATA 2016", reconstruction[i].Data());
		if(i != 7){
			pt_vs_met_DATA[i]->Divide(pt_vs_met_DATA[7]);
			pt_vs_met_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		}
		SalvaHisto(name_histo, pt_vs_met_DATA[i], save_document, "met", "p_{T}");
		name_histo = Form("%s p_{T} vs met: DATA 2016 / DATA 2017 RR", reconstruction[i].Data());
		pt_vs_met_DATA[i]->Divide(pt_vs_met_MC[i]);
		pt_vs_met_DATA[i]->GetZaxis()->SetRangeUser(0.0, 1.0);			
		SalvaHisto(name_histo, pt_vs_met_DATA[i], save_document, "met", "p_{T}", 0);

	}	
	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf]", 0,  0, "p_{T}");

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
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0);
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

	TH1F* last_hist = (TH1F *)h_MC->GetStack()->Last();	
		
	if(h_DATA->GetMaximum() > last_hist->GetMaximum())
		max = h_DATA->GetMaximum();
	else
		max = last_hist->GetMaximum();
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;

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

	last_hist->GetYaxis()->SetTitleOffset(1.2);
	h_DATA->GetYaxis()->SetTitleOffset(1.2);

	last_hist->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h_DATA->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
		
    c1->Update();

	h_MC->Draw();
	h_DATA->Draw("sameP");
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
	
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.3);
	pad22->SetTopMargin(0);
	pad22->SetBottomMargin(0.25);
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
	ratio->GetYaxis()->SetTitle("DATA 2016 / DATA 2017 RR");
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
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0);
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
// 	h_MC->SetMarkerStyle(20);
// 	h_MC->SetMarkerColor(kBlue);
// 	h_MC->SetMarkerSize(0.5);
	h_DATA->Draw("sameP");
	h_DATA->SetMarkerStyle(3);
	h_DATA->SetMarkerColor(kRed);
	h_DATA->SetMarkerSize(0.5);
	st_MC->Draw("same");
// 	st_DATA->Draw("same");
	
//     	
	TLegend *l1 = new TLegend(0.2,0.8,0.4,0.9);
	l1->AddEntry(h_DATA, "DATA 2016", "l");
	l1->AddEntry(h_MC, "DATA 2017 RR", "l");
	l1->Draw();	
	
	lat.DrawLatex(X_pos_latex, max, strig_1);
	lat.DrawLatex(X_pos_latex, 0.9*max, strig_2);
	lat.DrawLatex(X_pos_latex, 0.8*max, strig_3);
	lat.DrawLatex(X_pos_latex, 0.7*max, strig_4);
	lat.DrawLatex(X_pos_latex, 0.6*max, strig_5);
	
	c1->Update();
	c1->cd();
	
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.1, 1, 0.3);
	pad22->SetTopMargin(0);
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
	ratio->GetYaxis()->SetTitle("DATA 2016 / DATA 2017 RR");
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
