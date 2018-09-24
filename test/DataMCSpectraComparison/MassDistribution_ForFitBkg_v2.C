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


void MassDistribution_ForFitBkg_v2(){

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetLegendTextSize(0.03);
gROOT->LoadMacro("cruijff.C+");


	TH1F *mass = new TH1F("mass", "mass", 6000, 0, 6000);
	TH1F *BB_mass = new TH1F("BB_mass", "BB_mass", 6000, 0, 6000);
	TH1F *BE_mass = new TH1F("BE_mass", "BE_mass", 6000, 0, 6000);
	
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
	float weight_BB[15] = {0};
	float weight_BE[15] = {0};
	
 float gM;
 float kFactor;
 float kFactor_BB;
 float kFactor_BE;

	
	Long64_t ne;
	Long64_t Nne;
  float dil_mass;
  float cos_angle;
  float vertex_chi2;
  int dil_chosen;
  
  Int_t event;
  Int_t run;
  unsigned lumi;
  Int_t prev_event = -88;

  float lep_pt[2];
  float lep_phi[2];
  int lep_id[2];
  float lep_pt_err[2];
  float lep_eta[2];
  float lep_tk_pt[2];
  float lep_glb_pt[2];
  float lep_picky_pt[2];
  float lep_tpfms_pt[2];
  float lep_dB[2];
  float lep_sumPt[2];
  float lep_triggerMatchPt[2];
  short lep_glb_numberOfValidTrackerLayers[2]; 
  short lep_glb_numberOfValidPixelHits[2];
  short lep_glb_numberOfValidMuonHits[2];
  short lep_TuneP_numberOfValidMuonHits[2];
  short lep_numberOfMatchedStations[2];
  bool lep_isGlobalMuon[2];
  bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool GoodVtx;
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    float vertex_m;
    float gen_dil_mass;
  float genWeight;    
    float gen_lep_qOverPt[2];
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
//     float M[562242];
//     float MS[562242];
    
    
	TLorentzVector c_base;
	TLorentzVector c_daughter_0, c_daughter_1;

	TLorentzVector n_base;
	TLorentzVector n_daughter_0, n_daughter_1;

	int p_doppioni = 0;
	int n_doppioni = 0;	
	bool ok;
	int cont = -1;

	int c_run = -1;
	int c_lumi = -1;
	int c_event = -1;
	float c_mass = -1;
	float c_mass_scale = -1;
	
	int p_run = -1;
	int p_lumi = -1;
	int p_event = -1;
	float p_mass = -1;
	float p_mass_scale = -1;

	int n_run = -1;
	int n_lumi = -1;
	int n_event = -1;
	float n_mass = -1;
	float n_mass_scale = -1;

	int n_n_event = -1;
	float n_n_mass = -1;

	
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
	
	bool reco;
	bool noIB;
	
	reco = false;
	noIB = false;
			
	for(int i = 0; i < 15; i++){
		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight_BB[i] = weight[i];
		weight_BE[i] = weight[i];
	}
	
                        	
  for(int j=5; j < 6; j++){   


	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    

     TChain *treeMC = new TChain("SimpleNtupler/t");

     treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/ana_datamc_" + samples[j] + ".root");        
     treeMC->SetBranchAddress("genWeight",&genWeight);
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
     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
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

	int p_chosen = -99;

	std::cout<<"START"<<std::endl;

// 	for ( int p=0; p<ne ;p++){
	for ( int p=0; p<1000 ;p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
		treeMC->GetEntry(p);
		
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
				
				int counting = 0;
				int counting_OK = 0;
			
				do{
					treeMC->GetEntry(p);	
					c_event = event;
					c_mass = dil_mass;
// 					std::cout<<"A = "<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
					int p_1 = p + 1;
					treeMC->GetEntry(p_1);
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
						n_mass = dil_mass;
						counting_OK++;
// 						std::cout<<"B = "<<p_1<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
// 						std::cout<<"couting = "<<counting<<std::endl;
					}

					n_event = event;
					counting++;					
					p = p + 1;	
				}
				while(c_event == n_event);
				
				p -= counting;
				treeMC->GetEntry(p);	
				
				if(counting_OK > 0){
					std::vector<int> number_zed;
					for(int i = 0; i <= counting_OK; i++){
						treeMC->GetEntry(p+i);
						if(dil_mass > 70 && dil_mass < 110)
							number_zed.push_back(p+i);
						std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<std::endl;
					}
					if(number_zed.size() != 0){
						float max_pt = -99;
						for(int i = 0; i < number_zed.size(); i++){
							treeMC->GetEntry(number_zed.at(i));
							if(lep_pt[0] + lep_pt[1] > max_pt){
								max_pt = lep_pt[0] + lep_pt[1];
								p = number_zed.at(i);
							}
						}
					}
				}
				
				treeMC->GetEntry(p);	
// 				std::cout<<"Selected = "<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				
				weight[j] *= genWeight;
				weight_BB[j] *= genWeight;
				weight_BE[j] *= genWeight;
				
				gM = gen_dil_mass - 400;
				kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);

				if(j >= 6 && j <=14){
					weight[j] *= kFactor;
					weight_BB[j] *= kFactor_BB;
					weight_BE[j] *= kFactor_BE;
				}
				
				mass->Fill(vertex_m, weight[j]);
				
				if(fabs(lep_eta[0])<=1.2 && fabs(lep_eta[1]) <=1.2)				
					BB_mass->Fill(vertex_m, weight_BB[j]);
				else
					BE_mass->Fill(vertex_m, weight_BE[j]);
			
				if(j >= 6 && j <=14){
					weight[j] /= kFactor;
					weight_BB[j] /= kFactor_BB;
					weight_BE[j] /= kFactor_BE;
				}

				weight[j] /= genWeight;
				weight_BB[j] /= genWeight;
				weight_BE[j] /= genWeight;
				
				if(counting > 1)
					p += counting;


	       }//end of condition on event
	
		} // for event
	}
	
// 	TCanvas *c1 = new TCanvas("c1", "c1", 500 ,500);
// 	BB_mass->Draw();
// 	TCanvas *c2 = new TCanvas("c2", "c2", 500 ,500);
// 	BE_mass->Draw();
	TCanvas *c3 = new TCanvas("c3", "c3", 500 ,500);
	mass->Draw();
	
// 	TFile *f1 = new TFile("./MassDistributionForFit.root", "RECREATE");
// 	f1->cd();
// 	mass->Write();
// 	BB_mass->Write();
// 	BE_mass->Write();
// 	f1->Close();
	
}
