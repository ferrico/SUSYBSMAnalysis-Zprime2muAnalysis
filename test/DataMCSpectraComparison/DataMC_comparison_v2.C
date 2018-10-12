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
void SalvaHisto(TString name, TH1F* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "MC", TString strig_2 = "DATA", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex = 0, TString strig_1 = "", TString strig_2 = "", TString strig_3 = "", TString strig_4 = "", TString strig_5 = "");
void SalvaHisto(TString name, TH2F* h_MC, TString save, TString name_Xaxis, TString name_Yaxis, int Statistic = 1, bool logX = false);
void SaveMultipleCounting(TString name, TH1F* h[15], TH1F* h_DATA[15], TString save, int a, int b, int c, int d, int e);

  UInt_t event;
  UInt_t run;
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
         
// 	float LUMINOSITY = 78139.33;
// 	float LUMINOSITY = 41903.837;
	float LUMINOSITY = 36235.493;
	
	float Z_peak = 0.9638;
	float Z_peak_BB = 0.9688;
	float Z_peak_BE = 0.9610;

	float weight[45] = {0};

	const int    NMBINS = 100;
	const double MMIN = 60., MMAX = 2100.;
	double logMbins[NMBINS+1];

//     Double_t PT_BINS[] = {200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 600, 750, 1000, 1500};
//     Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;
    Double_t PT_BINS[] = {0, 50, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 
    					1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800};
    Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;    
    
//     Double_t ETA_BINS[] = {-2.4, -1.5, -1.2, -0.9, 0, 0.9, 1.2, 1.5, 2.4};
    Double_t ETA_BINS[] = {0, 0.9, 1.2, 1.5, 2.4};//, 5, 5.9, 6.2, 6.5, 7.4};
    Int_t  binnum_eta = sizeof(ETA_BINS)/sizeof(Double_t)-1;
    
//     Double_t PHI_BINS[] = {-3.14, -2.356, -1.57, -0.785, 0, 0.785, 1.57, 2.356, 3.14};
    Double_t PHI_BINS[] = {-3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2};
    Int_t  binnum_phi = sizeof(PHI_BINS)/sizeof(Double_t)-1;
    
    Double_t MET_BINS[] = {0, 25, 50, 75, 100, 200, 350, 500};//, 750, 1000};
    Int_t  binnum_met = sizeof(MET_BINS)/sizeof(Double_t)-1;
    
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 

	TH1F* pt_DATA = new TH1F("DATA_pt", "DATA_pt",binnum_pt, PT_BINS);
	TH1F* pt_DATA_BB = new TH1F("DATA_pt_BB", "DATA_pt_BB",binnum_pt, PT_BINS);
	TH1F* pt_DATA_BE = new TH1F("DATA_pt_BE", "DATA_pt_BE",binnum_pt, PT_BINS);
	THStack* pt_MC_Stack = new THStack("MC_pt", "MC_pt");
	THStack* pt_MC_BB_Stack = new THStack("MC_pt_BB", "MC_pt_BB");
	THStack* pt_MC_BE_Stack = new THStack("MC_pt_BE", "MC_pt_BE");
	TH1F* pt_MC_clear = new TH1F("MC_pt", "MC_pt", binnum_pt, PT_BINS);
	TH1F* pt_MC_BB_clear = new TH1F("MC_pt_BB", "MC_pt_BB", binnum_pt, PT_BINS);
	TH1F*pt_MC_BE_clear = new TH1F("MC_pt_BE", "MC_pt_BE", binnum_pt, PT_BINS);

	TH1F* eta_DATA = new TH1F("DATA_eta", "DATA_eta",binnum_eta, ETA_BINS);
	THStack* eta_MC_Stack = new THStack("MC_eta", "MC_eta");
	TH1F* eta_MC_clear = new TH1F("MC_eta", "MC_eta", binnum_eta, ETA_BINS);

	TH1F* phi_DATA = new TH1F("DATA_phi", "DATA_phi",binnum_phi, ETA_BINS);
	THStack* phi_MC_Stack = new THStack("MC_phi", "MC_phi");
	TH1F* phi_MC_clear = new TH1F("MC_phi", "MC_phi", binnum_phi, ETA_BINS);

	TH1F* dB_DATA = new TH1F("DATA_dB_TuneP", "DATA_dB_TuneP",  100, 0.0, 1);
	THStack* dB_MC_Stack = new THStack("MC_dB_TuneP", "MC_dB_TuneP");
	TH1F* dB_MC_clear = new TH1F("MC_dB_TuneP", "MC_dB_TuneP",  100, 0.0, 1);

	TH1F* PixelHit_DATA = new TH1F("DATA_PixelHit_TuneP", "DATA_PixelHit_TuneP",  100, 0, 100);
	THStack* PixelHit_MC_Stack = new THStack("MC_PixelHit_TuneP", "MC_PixelHit_TuneP");
	TH1F* PixelHit_MC_clear = new TH1F("MC_PixelHit_TuneP", "MC_PixelHit_TuneP",  100, 0, 100);

	TH1F* TkLayer_DATA = new TH1F("DATA_TkLayer_TuneP", "DATA_TkLayer_TuneP",  100, 0, 100);
	THStack* TkLayer_MC_Stack = new THStack("MC_TkLayer_TuneP", "MC_TkLayer_TuneP");
	TH1F* TkLayer_MC_clear = new TH1F("MC_TkLayer_TuneP", "MC_TkLayer_TuneP",  100, 0, 100);

	TH1F* Iso_DATA = new TH1F("DATA_Iso_TuneP", "DATA_Iso_TuneP",  100, 0.0, 1);
	THStack* Iso_MC_Stack = new THStack("MC_Iso_TuneP", "MC_Iso_TuneP");
	TH1F* Iso_MC_clear = new TH1F("MC_Iso_TuneP", "MC_Iso_TuneP",  100, 0.0, 1);

	TH1F* relpTErr_DATA = new TH1F("DATA_relpTErr_TuneP", "DATA_relpTErr_TuneP",  100, 0.0, 1);
	THStack* relpTErr_MC_Stack = new THStack("MC_relpTErr_TuneP", "MC_relpTErr_TuneP");
	TH1F* relpTErr_MC_clear = new TH1F("MC_relpTErr_TuneP", "MC_relpTErr_TuneP",  100, 0.0, 1);

	TH1F* Vtx_DATA = new TH1F("DATA_Vtx_TuneP", "DATA_Vtx_TuneP",  100, 0, 100);
	THStack* Vtx_MC_Stack = new THStack("MC_Vtx_TuneP", "MC_Vtx_TuneP");
	TH1F* Vtx_MC_clear = new TH1F("MC_Vtx_TuneP", "MC_Vtx_TuneP",  100, 0, 100);

	TH1F* ValidMu_DATA = new TH1F("DATA_ValidMu_TuneP", "DATA_ValidMu_TuneP",  100, 0, 100);
	THStack* ValidMu_MC_Stack = new THStack("MC_ValidMu_TuneP", "MC_ValidMu_TuneP");
	TH1F* ValidMu_MC_clear = new TH1F("MC_ValidMu_TuneP", "MC_ValidMu_TuneP",  100, 0, 100);

void DataMC_comparison(){

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
    	
	TH1F* dimuon_DATA = new TH1F("DATA_dimuon", "DATA_dimuon", MBINS, logMbins);
	THStack* dimuon_MC_Stack = new THStack("MC_dimuon", "MC_dimuon");
	TH1F* dimuon_MC_clear = new TH1F("MC_dimuon", "MC_dimuon", MBINS, logMbins);

	TH1F* dimuon_BB_DATA = new TH1F("DATA_dimuon_BB", "DATA_dimuon_BB", MBINS, logMbins);
	THStack* dimuon_BB_MC_Stack = new THStack("MC_dimuon_BB", "MC_dimuon_BB");
	TH1F* dimuon_BB_MC_clear = new TH1F("MC_dimuon_BB", "MC_dimuon_BB", MBINS, logMbins);

	TH1F* dimuon_BE_DATA = new TH1F("DATA_dimuon_BE", "DATA_dimuon_BE", MBINS, logMbins);
	THStack* dimuon_BE_MC_Stack = new THStack("MC_dimuon_BE", "MC_dimuon_BE");
	TH1F* dimuon_BE_MC_clear = new TH1F("MC_dimuon_BE", "MC_dimuon_BE", MBINS, logMbins);

	for(int i = 0; i < 45; i++){
		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
	}	
	
	for(int j = 14; j < 45; j++){
		
      	if(j > 38) continue; 
//       	if(j < 30) continue;
//       	if(j >= 30) continue;
//       	if(j >= 30 and j <= 38) continue;
//       	if(j == 30) continue; 
//       	if(j == 44) continue; 
//       	if(j >= 21) continue;
//       	if(j < 21) continue;

        std::cout<<"opening.. "<<samples[j]<<" --- "<<events[j]<<std::endl;

      	TChain *treeMC = new TChain("SimpleNtupler/t");

     	treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_Pt_Assignment/ana_datamc_" + samples[j] + ".root");     	

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
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples[j].Data(),j , weight[j], nentries);


    	for(int p=0; p<nentries; p++){
//      	 for(int p=0; p<10000; p++){
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
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
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
					
					dimuon_MC_clear->Fill(dil_mass,  weight[j]); 
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_MC_clear->Fill(dil_mass,  weight[j]*Z_peak_BB/Z_peak); 
					else
						dimuon_BE_MC_clear->Fill(dil_mass,  weight[j])*Z_peak_BE/Z_peak); 
					Vtx_MC_clear->Fill(vertex_chi2,  weight[j]); 
   					
					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_MC_clear->Fill(lep_pt[h],  weight[j]); 
						pt_MC_BB_clear->Fill(lep_pt[h],  weight[j]*Z_peak_BB/Z_peak); 
						pt_MC_BE_clear->Fill(lep_pt[h],  weight[j]*Z_peak_BE/Z_peak); 

						eta_MC_clear->Fill(lep_eta[h],  weight[j]); 

						phi_MC_clear->Fill(lep_phi[h],  weight[j]); 

						dB_MC_clear->Fill(lep_dB[h],  weight[j]); 

						PixelHit_MC_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]); 

						TkLayer_MC_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]); 

						Iso_MC_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]); 

						relpTErr_MC_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]); 

						ValidMu_MC_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]); 
   	 						    	 						
	     			} // for on muons

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

		 dimuon_MC_clear->SetFillColor(color);
		 dimuon_MC_clear->SetLineColor(color);
		 dimuon_BB_MC_clear->SetFillColor(color);
		 dimuon_BB_MC_clear->SetLineColor(color);
		 dimuon_BE_MC_clear->SetFillColor(color);
		 dimuon_BE_MC_clear->SetLineColor(color);
		pt_MC_clear->SetFillColor(color);
		pt_MC_clear->SetLineColor(color);
		pt_MC_BB_clear->SetFillColor(color);
		pt_MC_BB_clear->SetLineColor(color);
		pt_MC_BE_clear->SetFillColor(color);
		pt_MC_BE_clear->SetLineColor(color);
		Leading_pt_MC_clear->SetFillColor(color);
		Leading_pt_MC_clear->SetLineColor(color);
		Subleading_pt_MC_clear->SetFillColor(color);
		Subleading_pt_MC_clear->SetLineColor(color);
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
		pt_MC_Stack[i]->Add(pt_MC_clear, "HIST");
		pt_MC_BB_Stack[i]->Add(pt_MC_BB_clear, "HIST");
		pt_MC_BE_Stack[i]->Add(pt_MC_BE_clear, "HIST");
		Leading_pt_MC_Stack[i]->Add(Leading_pt_MC_clear, "HIST");
		Subleading_pt_MC_Stack[i]->Add(Subleading_pt_MC_clear, "HIST");
		eta_MC_Stack[i]->Add(eta_MC_clear, "HIST");
		phi_MC_Stack[i]->Add(phi_MC_clear, "HIST");
		dB_MC_Stack[i]->Add(dB_MC_clear, "HIST");
		PixelHit_MC_Stack[i]->Add(PixelHit_MC_clear, "HIST");
		TkLayer_MC_Stack[i]->Add(TkLayer_MC_clear, "HIST");
		Iso_MC_Stack[i]->Add(Iso_MC_clear, "HIST");
		relpTErr_MC_Stack[i]->Add(relpTErr_MC_clear, "HIST");
		Vtx_MC_Stack[i]->Add(Vtx_MC_clear, "HIST");
		ValidMu_MC_Stack[i]->Add(ValidMu_MC_clear, "HIST");
	 	 
    }// end loop on MC
    OutFile.close();
       
////////      DATA    ///////

	for(int k = 0; k < 1; k++){
	for(int i = 0; i < 7; i++)
		for(int j = 0; j < 2; j++)
			lepton_pt[j][i] = {0};
		
    TChain *treeDATA = new TChain("SimpleNtupler/t");
    if(k == 0) treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/DATA_Pt_Assignment/ana_datamc_data_Pt_Assignment.root");
    else  treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA/ana_datamc_data.root");

    	treeDATA->SetBranchAddress("genWeight",&genWeight);
     	    	
      	treeDATA->SetBranchAddress("event", &event);
      	treeDATA->SetBranchAddress("run", &run);
//       	treeDATA->SetBranchAddress("lumi", &lumi);    	
    	treeDATA->SetBranchAddress("dil_mass",&dil_mass);
    	
    	treeDATA->SetBranchAddress("gen_lep_pt", gen_lep_pt);
    	treeDATA->SetBranchAddress("gen_lep_eta", gen_lep_eta);
    	treeDATA->SetBranchAddress("gen_lep_phi", gen_lep_phi);
    	treeDATA->SetBranchAddress("gen_dil_mass",&gen_dil_mass);
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
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
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
					
					dimuon_DATA->Fill(dil_mass);
					if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2)
						dimuon_BB_DATA->Fill(dil_mass);
					else
						dimuon_BE_DATA->Fill(dil_mass);
					Vtx_DATA->Fill(vertex_chi2,  weight[j]); 
					Vtx_DATA->Fill(vertex_chi2);

					for(int h = 0; h < 2; h++){ // for on two muons in the event									    						

						pt_DATA->Fill(lep_pt[h]);
						pt_MC_BB_clear->Fill(lep_pt[h]);
						pt_MC_BE_clear->Fill(lep_pt[h]);

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
	} //Run16 and 17
	
		 
	gROOT->Reset();
	gROOT->SetBatch();

	TString dir_save = "./20162017/";

	save_document = dir_save + "Variuos_distribution.pdf";

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf[", 0,  0, "p_{T}");

	SalvaHisto("dimuon Mass: stack plot", dimuon_MC_Stack, dimuon_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
	SalvaHisto("dimuon Mass BB: stack plot", dimuon_BB_MC_Stack, dimuon_BB_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");
	SalvaHisto("dimuon Mass BE: stack plot", dimuon_BE_MC_Stack, dimuon_BE_DATA, dir_save + "Variuos_distribution.pdf", 1,  1, "m_{#mu#mu} [GeV]");

	name_histo = Form("%s p_{T}", reconstruction[i].Data());
	SalvaHisto(name_histo, pt_MC_Stack[i], pt_DATA[i], save_document, 0,  1, "p_{T}");

	name_histo = Form("%s p_{T} BB", reconstruction[i].Data());
	SalvaHisto(name_histo, pt_MC_BB_Stack[i], pt_DATA_BB[i], save_document, 0,  1, "p_{T}");

	name_histo = Form("%s p_{T} BE", reconstruction[i].Data());
	SalvaHisto(name_histo, pt_MC_BE_Stack[i], pt_DATA_BE[i], save_document, 0,  1, "p_{T}");

	name_histo = Form("%s #eta", reconstruction[i].Data());
	SalvaHisto(name_histo, eta_MC_Stack[i], eta_DATA[i], save_document, 0,  0, "#eta");

	name_histo = Form("%s #phi", reconstruction[i].Data());
	SalvaHisto(name_histo, phi_MC_Stack[i], phi_DATA[i], save_document, 0,  0, "#phi");

	name_histo = Form("dB", reconstruction[i].Data());
	SalvaHisto(name_histo, dB_MC_Stack[i], dB_DATA[i], save_document, 0,  1, "dB}");

	name_histo = Form("Pixel Hits", reconstruction[i].Data());
	SalvaHisto(name_histo, PixelHit_MC_Stack[i], PixelHit_DATA[i], save_document, 0,  1, "PixelHits");

	name_histo = Form("TkLayer", reconstruction[i].Data());
	SalvaHisto(name_histo, TkLayer_MC_Stack[i], TkLayer_DATA[i], save_document, 0,  1, "TkLayer");

	name_histo = Form("Rel Tk Iso", reconstruction[i].Data());
	SalvaHisto(name_histo, Iso_MC_Stack[i], Iso_DATA[i], save_document, 0,  0, "Rel Tk Iso");

	name_histo = Form("relPtErro", reconstruction[i].Data());
	SalvaHisto(name_histo, relpTErr_MC_Stack[i], relpTErr_DATA[i], save_document, 0,  0, "#sigma_{p_T}/p_T");

	name_histo = Form("Vtx", reconstruction[i].Data());
	SalvaHisto(name_histo, Vtx_MC_Stack[i], Vtx_DATA[i], save_document, 0,  0, "Vtx #chi^2");

	name_histo = Form("ValidMuHits", reconstruction[i].Data());
	SalvaHisto(name_histo, ValidMu_MC_Stack[i], ValidMu_DATA[i], save_document, 0,  0, "#mu hits");

	SalvaHisto("h_blank", h_blank, h_blank, dir_save + "Variuos_distribution.pdf]", 0,  0, "p_{T}");
	
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

void SalvaHisto(TString name, THStack* h_MC, TH1F* h_DATA, TString save, bool logx, bool logy, TString name_axis, int X_pos_latex, TString strig_1 = "MC", TString strig_2 = "DATA", TString strig_3, TString strig_4, TString strig_5){
			

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
	TLegend *l1 = new TLegend(0.3,0.8,0.5,0.9);
	l1->AddEntry(h_MC, strig_1, "l");
	l1->AddEntry(h_DATA, strig_2, "l");
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
