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
#include "TPaveStats.h"
#include "TPad.h"

#define n_bins 19
#define mass_bin 9

void Draw_MAss(TString nminus1, TH1F* MC_No_Cut_mass, TH1F* MC_Cut_mass, TH1F* DATA_No_Cut_mass, TH1F* DATA_Cut_mass, TLegend *legend_16);

void Nminus1_2017(){

gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetLegendTextSize(0.0275);


// 	    Double_t MASS_BINS[] = {120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
		Double_t MASS_BINS[] = {50, 120, 200, 400, 800, 1400, 2300, 3500};
	    Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;
	    
	    float Mass_BB_sigma[mass_bin] = {0};
	    float Mass_BB_sigma_err[mass_bin] = {0};
	    float Mass_BE_sigma[mass_bin] = {0};
	    float Mass_BE_sigma_err[mass_bin] = {0};
	    float Mass_EE_sigma[mass_bin] = {0};
	    float Mass_EE_sigma_err[mass_bin] = {0};
	    float Mass_OVER_sigma[mass_bin] = {0};
	    float Mass_OVER_sigma_err[mass_bin] = {0};

	    Double_t PT_BINS[] = {50, 100, 150, 200, 250, 
	    					300, 350, 400, 450, 500, 
	    					600, 700, 800, 900, 1000, 1250, 1500, 2000, 3000};
// 	    Double_t PT_BINS[] = {(float)1/(float)3000, (float)1/(float)2000, (float)1/(float)1500, (float)1/(float)1250, (float)1/(float)1000,
// 	    					(float)1/(float)900, (float)1/(float)800, (float)1/(float)700, (float)1/(float)600, (float)1/(float)500,
// 	    					(float)1/(float)450, (float)1/(float)400, (float)1/(float)350, (float)1/(float)300, (float)1/(float)250, 
// 	    					(float)1/(float)200, (float)1/(float)150, (float)1/(float)100, (float)1/(float)50};
	    Int_t  pt_binnum = sizeof(PT_BINS)/sizeof(Double_t)-1;
	    
	    float Pt_BB_sigma[n_bins] = {0};
	    float Pt_BB_sigma_err[n_bins] = {0};
	    float Pt_BE_sigma[n_bins] = {0};
	    float Pt_BE_sigma_err[n_bins] = {0};
	    float Pt_EE_sigma[n_bins] = {0};
	    float Pt_EE_sigma_err[n_bins] = {0};
	    float Pt_OVER_sigma[n_bins] = {0};
	    float Pt_OVER_sigma_err[n_bins] = {0};

	    
          std::vector<float> SIGMA;
          std::vector<float> SIGMA_ERR;
		
    TString samples[15] =  {
                            "Wantitop", "tW", 
                            "ttbar",
                            "WW",
                            "ZZ",
                            "WZ",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
                            };
	
	float events[15] = {
						7780870, 7581624,
						33844772,
						7791498,
						1949768,
						3928630,
						2961000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000,
					};
						
	float sigma[15] = {
						35.6, 35.6,
						831.76,
						118.7,
						16.523,
						47.13,
						2112.905, 20.553, 2.8861, 0.25126, 0.017075, 1.366E-3, 8.178E-5, 3.191E-6, 2.787E-7,
						};         
         
	float LUMINOSITY = 41903.837;


	float weight[15] = {0};
	
	TString NOME;
	float  fit_min, fit_max;
	float yPos;
	TString longstring;
	

	
	
	Long64_t ne;
	Long64_t Nne;
  float dil_mass;
  float cos_angle;
  float vertex_chi2;
  int dil_chosen;
  
  Int_t event;
  Int_t run;
  unsigned lumi;

  float lep_pt[2];
  float lep_phi[2];
  int lep_id[2];
  float lep_pt_err[2];
  float lep_eta[2];
  float lep_tk_pt[2];
  float lep_tk_eta[2];
  float lep_tk_phi[2];
  float lep_std_pt[2];
  float lep_std_eta[2];
  float lep_std_phi[2];
  float lep_picky_pt[2];
  float lep_picky_eta[2];
  float lep_picky_phi[2];
  float lep_dyt_pt[2];
  float lep_dyt_eta[2];
  float lep_dyt_phi[2];
  float lep_glb_pt[2];
  float lep_tpfms_pt[2];
  float lep_dB[2];
  float lep_sumPt[2];
  float lep_triggerMatchPt[2];
  short lep_glb_numberOfValidTrackerLayers[2]; 
  short lep_glb_numberOfValidPixelHits[2];
  short lep_glb_numberOfValidMuonHits[2];
  short lep_numberOfMatchedStations[2];
  bool lep_isGlobalMuon[2];
  bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool GoodVtx;
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    float vertex_m;
    float gen_dil_mass;
    float gen_lep_qOverPt[2];
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
	
	float MASS = -1;
	float MASS_SCALE = -1;
	float MASS_GEN = -1;
	int i = 0;
	int n_count;
	int p_count;
	
	float gen_eta_first = -999;
	float gen_eta_second = -999;
	
	float gen_pt_first = -999;
	float gen_pt_second = -999;
	
	float pt_first = -999;
	float pt_second = -999;
	float max_pt = -999;
	
	float gen_pt_plus = -999;
	float gen_pt_minus = -999;
		
	int pp;
	bool next;
	
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 
	
  Int_t prev_event = -88;
  
	float count_event[9] = {0};
	float count_event_BB[9] = {0};
	float count_event_BE[9] = {0};

			
	for(int i = 0; i < 15; i++){
		weight[i] = LUMINOSITY * sigma[i] / events[i];
// 		weight[i] = 1;
	}


	TH1F * MC_Cut_BELOW = new TH1F("MC_Cut_BELOW", "MC_Cut_BELOW", 11, -0.5, 10.5);
	TH1F * MC_No_Cut_BELOW = new TH1F("MC_No_Cut_BELOW", "MC_No_Cut_BELOW", 11, -0.5, 10.5);
	TH1F * DATA_Cut_BELOW = new TH1F("DATA_Cut_BELOW", "DATA_Cut_BELOW", 11, -0.5, 10.5);
	TH1F * DATA_No_Cut_BELOW = new TH1F("DATA_No_Cut_BELOW", "DATA_No_Cut_BELOW", 11, -0.5, 10.5);
	TH1F * MC_Cut_ABOVE = new TH1F("MC_Cut_ABOVE", "MC_Cut_ABOVE", 11, -0.5, 10.5);
	TH1F * MC_No_Cut_ABOVE = new TH1F("MC_No_Cut_ABOVE", "MC_No_Cut_ABOVE", 11, -0.5, 10.5);
	TH1F * DATA_Cut_ABOVE = new TH1F("DATA_Cut_ABOVE", "DATA_Cut_ABOVE", 11, -0.5, 10.5);
	TH1F * DATA_No_Cut_ABOVE = new TH1F("DATA_No_Cut_ABOVE", "DATA_No_Cut_ABOVE", 11, -0.5, 10.5);

	TH1F * MC_No_Cut_mass = new TH1F("MC_No_Cut_mass", "MC_No_Cut_mass", binnum, MASS_BINS);
	TH1F * DATA_No_Cut_mass = new TH1F("DATA_No_Cut_mass", "DATA_No_Cut_mass", binnum, MASS_BINS);

	TH1F * MC_Cut_mass_PT = new TH1F("MC_Cut_mass_PT", "MC_Cut_mass_PT", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_DB = new TH1F("MC_Cut_mass_DB", "MC_Cut_mass_DB", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_ISO = new TH1F("MC_Cut_mass_ISO", "MC_Cut_mass_ISO", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_TK = new TH1F("MC_Cut_mass_TK", "MC_Cut_mass_TK", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_PX = new TH1F("MC_Cut_mass_PX", "MC_Cut_mass_PX", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_HIT = new TH1F("MC_Cut_mass_HIT", "MC_Cut_mass_HIT", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_STATION = new TH1F("MC_Cut_mass_STATION", "MC_Cut_mass_STATION", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_VTX = new TH1F("MC_Cut_mass_VTX", "MC_Cut_mass_VTX", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_BACK = new TH1F("MC_Cut_mass_BACK", "MC_Cut_mass_BACK", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_SIGMA = new TH1F("MC_Cut_mass_SIGMA", "MC_Cut_mass_SIGMA", binnum, MASS_BINS);
	TH1F * MC_Cut_mass_TRIG = new TH1F("MC_Cut_mass_TRIG", "MC_Cut_mass_TRIG", binnum, MASS_BINS);

	TH1F * DATA_Cut_mass_PT = new TH1F("DATA_Cut_mass_PT", "DATA_Cut_mass_PT", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_DB = new TH1F("DATA_Cut_mass_DB", "DATA_Cut_mass_DB", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_ISO = new TH1F("DATA_Cut_mass_ISO", "DATA_Cut_mass_ISO", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_TK = new TH1F("DATA_Cut_mass_TK", "DATA_Cut_mass_TK", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_PX = new TH1F("DATA_Cut_mass_PX", "DATA_Cut_mass_PX", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_HIT = new TH1F("DATA_Cut_mass_HIT", "DATA_Cut_mass_HIT", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_STATION = new TH1F("DATA_Cut_mass_STATION", "DATA_Cut_mass_STATION", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_VTX = new TH1F("DATA_Cut_mass_VTX", "DATA_Cut_mass_VTX", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_BACK = new TH1F("DATA_Cut_mass_BACK", "DATA_Cut_mass_BACK", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_SIGMA = new TH1F("DATA_Cut_mass_SIGMA", "DATA_Cut_mass_SIGMA", binnum, MASS_BINS);
	TH1F * DATA_Cut_mass_TRIG = new TH1F("DATA_Cut_mass_TRIG", "DATA_Cut_mass_TRIG", binnum, MASS_BINS);
	
			
  for(int j=0; j < 15; j++){  
  
	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    

     TChain *treeMC = new TChain("SimpleNtupler/t");

//      treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_Pt_Assignment/ana_datamc_"+samples[j]+".root");
     treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/ana_datamc_" + samples[j] + ".root");        

     treeMC->SetBranchAddress("event",&event);
     treeMC->SetBranchAddress("run",&run);
     treeMC->SetBranchAddress("lumi",&lumi);
     treeMC->SetBranchAddress("dil_mass",&dil_mass);
     treeMC->SetBranchAddress("cos_angle",&cos_angle);
     treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
     treeMC->SetBranchAddress("lep_pt",lep_pt);
     treeMC->SetBranchAddress("lep_id",lep_id);
     treeMC->SetBranchAddress("lep_eta",lep_eta);
     treeMC->SetBranchAddress("lep_phi",lep_phi);
     treeMC->SetBranchAddress("lep_dB",lep_dB);
     treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
     treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeMC->SetBranchAddress("lep_tk_eta",lep_tk_eta);
     treeMC->SetBranchAddress("lep_tk_phi",lep_tk_phi);
//      treeMC->SetBranchAddress("lep_std_pt",lep_std_pt);
//      treeMC->SetBranchAddress("lep_std_eta",lep_std_eta);
//      treeMC->SetBranchAddress("lep_std_phi",lep_std_phi);
     treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeMC->SetBranchAddress("lep_picky_eta",lep_picky_eta);
     treeMC->SetBranchAddress("lep_picky_phi",lep_picky_phi);
//      treeMC->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
//      treeMC->SetBranchAddress("lep_dyt_eta",lep_dyt_eta);
//      treeMC->SetBranchAddress("lep_dyt_phi",lep_dyt_phi);
     treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeMC->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeMC->SetBranchAddress("vertex_m",&vertex_m);
     treeMC->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);
     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
     treeMC->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
     treeMC->SetBranchAddress("gen_lep_qOverPt", gen_lep_qOverPt);
     treeMC->SetBranchAddress("gen_lep_eta", gen_lep_eta);
     treeMC->SetBranchAddress("gen_lep_pt", gen_lep_pt);	

	
	ne = treeMC->GetEntries();
	std::cout<<"START"<<std::endl;
	for ( int p=0; p < ne ;p++){
// 	for ( int p=0; p<100000 ;p++){
// 	for ( int p=0; p<0 ;p++){
	
		if(p % 100000 == 0) std::cout<<p<<" su "<<ne<<std::endl;		
			
		treeMC->GetEntry(p);

// 		if(dil_mass < 60 || dil_mass > 120) continue;
				
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120){
					MC_No_Cut_BELOW->Fill(0.0, weight[j]);
					MC_No_Cut_BELOW->Fill(1, weight[j]);
					MC_No_Cut_BELOW->Fill(2, weight[j]);
					MC_No_Cut_BELOW->Fill(3, weight[j]);
					MC_No_Cut_BELOW->Fill(4, weight[j]);
					MC_No_Cut_BELOW->Fill(5, weight[j]);
					MC_No_Cut_BELOW->Fill(6, weight[j]);
					MC_No_Cut_BELOW->Fill(7, weight[j]);
					MC_No_Cut_BELOW->Fill(8, weight[j]);
					MC_No_Cut_BELOW->Fill(9, weight[j]);
					MC_No_Cut_BELOW->Fill(10, weight[j]);
				}
				if(dil_mass > 120){
					MC_No_Cut_ABOVE->Fill(0.0, weight[j]);
					MC_No_Cut_ABOVE->Fill(1, weight[j]);
					MC_No_Cut_ABOVE->Fill(2, weight[j]);
					MC_No_Cut_ABOVE->Fill(3, weight[j]);
					MC_No_Cut_ABOVE->Fill(4, weight[j]);
					MC_No_Cut_ABOVE->Fill(5, weight[j]);
					MC_No_Cut_ABOVE->Fill(6, weight[j]);
					MC_No_Cut_ABOVE->Fill(7, weight[j]);
					MC_No_Cut_ABOVE->Fill(8, weight[j]);
					MC_No_Cut_ABOVE->Fill(9, weight[j]);
					MC_No_Cut_ABOVE->Fill(10, weight[j]);
				}
				MC_No_Cut_mass->Fill(dil_mass, weight[j]);
				

			}
			
		// PT	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
// 			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(0.0, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(0.0, weight[j]);
				MC_Cut_mass_PT->Fill(dil_mass, weight[j]);
				

			}


		// DB
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
// 			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(1, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(1, weight[j]);
				MC_Cut_mass_DB->Fill(dil_mass, weight[j]);

			}
			
		// ISOL	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
// 			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(2, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(2, weight[j]);
				MC_Cut_mass_ISO->Fill(dil_mass, weight[j]);
				
			}
			
		// TRK LAYERS	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
// 			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(3, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(3, weight[j]);
				MC_Cut_mass_TK->Fill(dil_mass, weight[j]);

			}
			
		// PIXEL HITS	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
// 			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(4, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(4, weight[j]);
				MC_Cut_mass_PX->Fill(dil_mass, weight[j]);
			}
			
		// MUON HITS	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
// 			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(5, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(5, weight[j]);
				MC_Cut_mass_HIT->Fill(dil_mass, weight[j]);
			}
			
		// MATCHED STATIONS
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
// 			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(6, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(6, weight[j]);
				MC_Cut_mass_STATION->Fill(dil_mass, weight[j]);
			}
			
		// VTX
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) //&& 
// 			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(7, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(7, weight[j]);
				MC_Cut_mass_VTX->Fill(dil_mass, weight[j]);
			}
			
		// BACK - TO - BACK	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
// 			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(8, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(8, weight[j]);
				MC_Cut_mass_BACK->Fill(dil_mass, weight[j]);
			}
			
		// sigmaPT	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
// 			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(9, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(9, weight[j]);
				MC_Cut_mass_SIGMA->Fill(dil_mass, weight[j]);

			}
			
		// TRIGGER	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
// 			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) MC_Cut_BELOW->Fill(10, weight[j]);
				if(dil_mass > 120) 					MC_Cut_ABOVE->Fill(10, weight[j]);
				MC_Cut_mass_TRIG->Fill(dil_mass, weight[j]);
			}			
	} // for entries
	
	
	std::cout<<"ESCO DAL SAMPLES"<<std::endl;
	
 } // for samples
 
     TChain *treeDATA = new TChain("SimpleNtupler/t");
     treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA/ana_datamc_data.root");

     treeDATA->SetBranchAddress("event",&event);
     treeDATA->SetBranchAddress("run",&run);
     treeDATA->SetBranchAddress("lumi",&lumi);
     treeDATA->SetBranchAddress("dil_mass",&dil_mass);
     treeDATA->SetBranchAddress("cos_angle",&cos_angle);
     treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeDATA->SetBranchAddress("dil_chosen",&dil_chosen);
     treeDATA->SetBranchAddress("lep_pt",lep_pt);
     treeDATA->SetBranchAddress("lep_id",lep_id);
     treeDATA->SetBranchAddress("lep_eta",lep_eta);
     treeDATA->SetBranchAddress("lep_phi",lep_phi);
     treeDATA->SetBranchAddress("lep_dB",lep_dB);
     treeDATA->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeDATA->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
     treeDATA->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeDATA->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeDATA->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeDATA->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeDATA->SetBranchAddress("lep_tk_eta",lep_tk_eta);
     treeDATA->SetBranchAddress("lep_tk_phi",lep_tk_phi);
     treeDATA->SetBranchAddress("lep_std_pt",lep_std_pt);
     treeDATA->SetBranchAddress("lep_std_eta",lep_std_eta);
     treeDATA->SetBranchAddress("lep_std_phi",lep_std_phi);
     treeDATA->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeDATA->SetBranchAddress("lep_picky_eta",lep_picky_eta);
     treeDATA->SetBranchAddress("lep_picky_phi",lep_picky_phi);
     treeDATA->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
     treeDATA->SetBranchAddress("lep_dyt_eta",lep_dyt_eta);
     treeDATA->SetBranchAddress("lep_dyt_phi",lep_dyt_phi);
     treeDATA->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeDATA->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeDATA->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeDATA->SetBranchAddress("vertex_m",&vertex_m);
     treeDATA->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);
     treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	Long64_t nentries = treeDATA->GetEntries();
	
	printf("opening... DATA --- %lld\n", nentries);    
		
	for(int p=0; p<nentries; p++){
// 	for(int p=0; p<100000; p++){
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

				if(dil_mass > 60 && dil_mass < 120){
					DATA_No_Cut_BELOW->Fill(0.0);
					DATA_No_Cut_BELOW->Fill(1);
					DATA_No_Cut_BELOW->Fill(2);
					DATA_No_Cut_BELOW->Fill(3);
					DATA_No_Cut_BELOW->Fill(4);
					DATA_No_Cut_BELOW->Fill(5);
					DATA_No_Cut_BELOW->Fill(6);
					DATA_No_Cut_BELOW->Fill(7);
					DATA_No_Cut_BELOW->Fill(8);
					DATA_No_Cut_BELOW->Fill(9);
					DATA_No_Cut_BELOW->Fill(10);
				}
				if(dil_mass > 120){
					DATA_No_Cut_ABOVE->Fill(0.0);
					DATA_No_Cut_ABOVE->Fill(1);
					DATA_No_Cut_ABOVE->Fill(2);
					DATA_No_Cut_ABOVE->Fill(3);
					DATA_No_Cut_ABOVE->Fill(4);
					DATA_No_Cut_ABOVE->Fill(5);
					DATA_No_Cut_ABOVE->Fill(6);
					DATA_No_Cut_ABOVE->Fill(7);
					DATA_No_Cut_ABOVE->Fill(8);
					DATA_No_Cut_ABOVE->Fill(9);
					DATA_No_Cut_ABOVE->Fill(10);
				}
				
				DATA_No_Cut_mass->Fill(dil_mass);
			}
		
		
		// PT	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
// 			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(0.0);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(0.0);
				DATA_Cut_mass_PT->Fill(dil_mass);		
			}


		// DB
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
// 			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(1);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(1);
				DATA_Cut_mass_DB->Fill(dil_mass);
			}
			
		// ISOL	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
// 			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(2);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(2);
				DATA_Cut_mass_ISO->Fill(dil_mass);

			}
			
		// TRK LAYERS	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
// 			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(3);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(3);
				DATA_Cut_mass_TK->Fill(dil_mass);

			}
			
		// PIXEL HITS	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
// 			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(4);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(4);
				DATA_Cut_mass_PX->Fill(dil_mass);
			}
			
		// MUON HITS	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
// 			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(5);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(5);
				DATA_Cut_mass_HIT->Fill(dil_mass);
			}
			
		// MATCHED STATIONS
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
// 			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
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
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(6);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(6);
				DATA_Cut_mass_STATION->Fill(dil_mass);
			}
			
		// VTX
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) //&& 
// 			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(7);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(7);
				DATA_Cut_mass_VTX->Fill(dil_mass);
			}
			
		// BACK - TO - BACK	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
// 			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(8);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(8);
				DATA_Cut_mass_BACK->Fill(dil_mass);
			}
			
		// sigmaPT	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
// 			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(9);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(9);
				DATA_Cut_mass_SIGMA->Fill(dil_mass);
			}
			
		// TRIGGER	
		if( 
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
// 			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;
				
				if(dil_mass > 60 && dil_mass < 120) DATA_Cut_BELOW->Fill(10);
				if(dil_mass > 120) 					DATA_Cut_ABOVE->Fill(10);
				DATA_Cut_mass_TRIG->Fill(dil_mass);
			}	

	}
 
 

// gROOT->Reset();
// gROOT->SetBatch();

TLegend *legend_16 = new TLegend(0.6, 0.2, 0.9, 0.4);

TEfficiency* MC_eff_BELOW = 0;
TEfficiency* DATA_eff_BELOW = 0;
MC_Cut_BELOW->Sumw2();
MC_No_Cut_BELOW->Sumw2();
MC_eff_BELOW = new TEfficiency(*MC_No_Cut_BELOW, *MC_Cut_BELOW);
DATA_Cut_BELOW->Sumw2();
DATA_No_Cut_BELOW->Sumw2();
DATA_eff_BELOW = new TEfficiency(*DATA_No_Cut_BELOW, *DATA_Cut_BELOW);
TCanvas *c3 = new TCanvas("BELOW Z peak", "60 < m_{ll} < 120", 500, 500);
c3->cd();
TPad* pad11 = new TPad("pad1", "pad1", 0, 0.05, 1, 1);
pad11->SetBottomMargin(0.15);
pad11->SetTicks();
pad11->Draw();
pad11->cd();
MC_eff_BELOW->SetMarkerStyle(20);
MC_eff_BELOW->SetMarkerSize(1);
MC_eff_BELOW->SetMarkerColor(kGreen+3);
MC_eff_BELOW->SetFillColor(kGreen+3);
MC_eff_BELOW->SetFillStyle(1001);
MC_eff_BELOW->Draw("AP2");
DATA_eff_BELOW->SetMarkerStyle(20);
DATA_eff_BELOW->SetMarkerSize(1);
DATA_eff_BELOW->SetMarkerColor(kBlack);
DATA_eff_BELOW->Draw("same P");
MC_eff_BELOW->SetTitle("60 < m_{ll} < 120; ciao ; #epsilon"); 
DATA_eff_BELOW->SetTitle("60 < m_{ll} < 120; ciao ; #epsilon"); 
MC_No_Cut_BELOW->SetMinimum(0.85);
MC_No_Cut_BELOW->SetMaximum(1.1);
MC_Cut_BELOW->SetMinimum(0.85);
MC_Cut_BELOW->SetMaximum(1.1);
DATA_No_Cut_BELOW->SetMinimum(0.85);
DATA_No_Cut_BELOW->SetMaximum(1.1);
DATA_Cut_BELOW->SetMinimum(0.85);
DATA_Cut_BELOW->SetMaximum(1.1);
legend_16->AddEntry(MC_eff_BELOW, "Simulation", "f");
legend_16->AddEntry(DATA_eff_BELOW, "2016 Data 36.2 /fb", "lep");
legend_16->Draw();

TEfficiency* MC_eff_ABOVE = 0;
TEfficiency* DATA_eff_ABOVE = 0;
MC_Cut_ABOVE->Sumw2();
MC_No_Cut_ABOVE->Sumw2();
MC_eff_ABOVE = new TEfficiency(*MC_No_Cut_ABOVE, *MC_Cut_ABOVE);
DATA_Cut_ABOVE->Sumw2();
DATA_No_Cut_ABOVE->Sumw2();
DATA_eff_ABOVE = new TEfficiency(*DATA_No_Cut_ABOVE, *DATA_Cut_ABOVE);
TCanvas *c4 = new TCanvas("ABOVE Z peak", "m_{ll} > 120", 500, 500);
c4->cd();
TPad* pad12 = new TPad("pad1", "pad1", 0, 0.05, 1, 1);
pad12->SetBottomMargin(0.15);
pad12->SetTicks();
pad12->Draw();
pad12->cd();
MC_eff_ABOVE->SetMarkerStyle(20);
MC_eff_ABOVE->SetMarkerSize(1);
MC_eff_ABOVE->SetMarkerColor(kGreen+3);
MC_eff_ABOVE->SetFillColor(kGreen+3);
MC_eff_ABOVE->SetFillStyle(1001);
MC_eff_ABOVE->Draw("AP2");
DATA_eff_ABOVE->SetMarkerStyle(20);
DATA_eff_ABOVE->SetMarkerSize(1);
DATA_eff_ABOVE->SetMarkerColor(kBlack);
DATA_eff_ABOVE->Draw("same P");
MC_eff_ABOVE->SetTitle("m_{ll} > 120; ciao ; #epsilon"); 
DATA_eff_ABOVE->SetTitle("m_{ll} > 120; ciao ; #epsilon"); 
legend_16->Draw();

Draw_MAss("PT", MC_No_Cut_mass, MC_Cut_mass_PT, DATA_No_Cut_mass, DATA_Cut_mass_PT, legend_16);
Draw_MAss("DB", MC_No_Cut_mass, MC_Cut_mass_DB, DATA_No_Cut_mass, DATA_Cut_mass_DB, legend_16);
Draw_MAss("ISO", MC_No_Cut_mass, MC_Cut_mass_ISO, DATA_No_Cut_mass, DATA_Cut_mass_ISO, legend_16);
Draw_MAss("TK", MC_No_Cut_mass, MC_Cut_mass_TK, DATA_No_Cut_mass, DATA_Cut_mass_TK, legend_16);
Draw_MAss("PX", MC_No_Cut_mass, MC_Cut_mass_PX, DATA_No_Cut_mass, DATA_Cut_mass_PX, legend_16);
Draw_MAss("HIT", MC_No_Cut_mass, MC_Cut_mass_HIT, DATA_No_Cut_mass, DATA_Cut_mass_HIT, legend_16);
Draw_MAss("STATION", MC_No_Cut_mass, MC_Cut_mass_STATION, DATA_No_Cut_mass, DATA_Cut_mass_STATION, legend_16);
Draw_MAss("VTX", MC_No_Cut_mass, MC_Cut_mass_VTX, DATA_No_Cut_mass, DATA_Cut_mass_VTX, legend_16);
Draw_MAss("BACK", MC_No_Cut_mass, MC_Cut_mass_BACK, DATA_No_Cut_mass, DATA_Cut_mass_BACK, legend_16);
Draw_MAss("SIGMA", MC_No_Cut_mass, MC_Cut_mass_SIGMA, DATA_No_Cut_mass, DATA_Cut_mass_SIGMA, legend_16);
Draw_MAss("TRIG", MC_No_Cut_mass, MC_Cut_mass_TRIG, DATA_No_Cut_mass, DATA_Cut_mass_TRIG, legend_16);


return;

} // main function


void Draw_MAss(TString nminus1, TH1F* MC_No_Cut_mass, TH1F* MC_Cut_mass, TH1F* DATA_No_Cut_mass, TH1F* DATA_Cut_mass, TLegend *legend_16){
	TEfficiency* MC_mass = 0;
	MC_No_Cut_mass->Sumw2();
	MC_Cut_mass->Sumw2();
	MC_mass = new TEfficiency(*MC_No_Cut_mass, *MC_Cut_mass);
	TEfficiency* DATA_mass = 0;
	DATA_No_Cut_mass->Sumw2();
	DATA_Cut_mass->Sumw2();
	DATA_mass = new TEfficiency(*DATA_No_Cut_mass, *DATA_Cut_mass);
	TCanvas *c2 = new TCanvas(nminus1, nminus1, 500, 500);
	MC_mass->SetMarkerStyle(20);
	MC_mass->SetMarkerSize(1);
	MC_mass->SetMarkerColor(kGreen+3);
	MC_mass->SetFillColor(kGreen+3);
	MC_mass->SetFillStyle(1001);
	MC_mass->Draw("AP2");
	DATA_mass->SetMarkerStyle(20);
	DATA_mass->SetMarkerSize(1);
	DATA_mass->SetMarkerColor(kBlack);
	DATA_mass->Draw("same P");
	MC_mass->SetTitle("N-1 vs Mass; m_{ll}; #epsilon"); 
	DATA_mass->SetTitle("N-1 vs Mass; m_{ll}; #epsilon"); 
	legend_16->Draw();
}
