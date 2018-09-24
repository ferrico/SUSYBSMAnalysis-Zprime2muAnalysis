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

Double_t MASS_BINS[] = {50, 120, 200, 400, 800, 1400, 2300, 3500, 4500, 6000};
// Double_t MASS_BINS[] = {0, 100, 200, 300, 400, 600, 800, 1000, 1400, 1800, 2200, 2800, 3400, 4000, 5000, 6000};
// Double_t MASS_BINS[] = {0, 120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;

TString NOME;

float  fit_min, fit_max;

void ExtractResolution(TH2F *h_2, TString NOME, TString save, TH1F* res, TH1F* mean);
void SalvaRatio(TString name, TH1F* h_NUM, TH1F* h_DEN, TString save);

void MassScale_MC_v1(){

gROOT->LoadMacro("cruijff.C+");

gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetLegendTextSize(0.03);
	
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
    
	int c_run = -1;
	int c_lumi = -1;
	int c_event = -1;
	float c_mass = -1;

	int n_run = -1;
	int n_lumi = -1;
	int n_event = -1;
	float n_mass = -1;

	TLorentzVector c_base;
	TLorentzVector c_daughter_0, c_daughter_1;
	float MASS = -1;
	float MASS_SCALE = -1;
	float MASS_GEN = -1;
	
	for(int i = 0; i < 15; i++){
		weight[i] = LUMINOSITY * sigma[i] / events[i];
// 		weight[i] = 1;
		weight_BB[i] = weight[i];
		weight_BE[i] = weight[i];
	}
	
	TH2F *res_vs_mass_BB_NS = new TH2F("Res vs mass BB", "Res vs mass BB", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_BB_S = new TH2F("Res vs mass BB: scale", "Res vs mass BB: scale", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_BE_NS = new TH2F("Res vs mass BE+EE", "Res vs mass BE+EE", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_BE_S = new TH2F("Res vs mass BE+EE: scale", "Res vs mass BE+EE: scale", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_09 = new TH2F("Res vs mass 09", "Res vs mass 09", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_EE = new TH2F("Res vs mass EE", "Res vs mass EE", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_OVER = new TH2F("Res vs mass Overlap", "Res vs mass Overlap", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F *res_vs_mass_OTHER = new TH2F("Res vs mass BE", "Res vs mass BE", binnum, MASS_BINS, 400, -1., 1.0);

	TH1F *resolution_BB_NS = new TH1F("Resolution BB", "Resolution BB", binnum, MASS_BINS);
	TH1F *resolution_BB_S = new TH1F("Resolution BB: scale", "Resolution BB: scale", binnum, MASS_BINS);
	TH1F *resolution_BE_NS = new TH1F("Resolution BE+EE", "Resolution BE+EE", binnum, MASS_BINS);
	TH1F *resolution_BE_S = new TH1F("Resolution BE+EE: scale", "Resolution BE+EE: scale", binnum, MASS_BINS);
	TH1F *resolution_09 = new TH1F("Resolution 09", "Resolution 09", binnum, MASS_BINS);
	TH1F *resolution_EE = new TH1F("Resolution EE", "Resolution EE", binnum, MASS_BINS);
	TH1F *resolution_OVER = new TH1F("Resolution Overlap", "Resolution Overlap", binnum, MASS_BINS);
	TH1F *resolution_OTHER = new TH1F("Resolution BE", "Resolution BE", binnum, MASS_BINS);

	TH1F *mean_BB_NS = new TH1F("Mean BB", "Mean BB", binnum, MASS_BINS);
	TH1F *mean_BB_S = new TH1F("Mean BB: scale", "Mean BB: scale", binnum, MASS_BINS);
	TH1F *mean_BE_NS = new TH1F("Mean BE+EE", "Mean BE+EE", binnum, MASS_BINS);
	TH1F *mean_BE_S = new TH1F("Mean BE+EE: scale", "Mean BE+EE: scale", binnum, MASS_BINS);
	TH1F *mean_09 = new TH1F("Mean 09", "Mean 09", binnum, MASS_BINS);
	TH1F *mean_EE = new TH1F("Mean EE", "Mean EE", binnum, MASS_BINS);
	TH1F *mean_OVER = new TH1F("Mean Overlap", "Mean Overlap", binnum, MASS_BINS);
	TH1F *mean_OTHER = new TH1F("Mean BE", "Mean BE", binnum, MASS_BINS);
	
                        	
  for(int j=6; j < 15; j++){   

	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    

     TChain *treeSCALE = new TChain("SimpleNtupler/t");

     treeSCALE->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/DY_kFactor_Scale/ana_datamc_" + samples[j] + ".root");        

     treeSCALE->SetBranchAddress("genWeight",&genWeight);
     treeSCALE->SetBranchAddress("event",&event);
     treeSCALE->SetBranchAddress("run",&run);
     treeSCALE->SetBranchAddress("lumi",&lumi);
     treeSCALE->SetBranchAddress("dil_mass",&dil_mass);
     treeSCALE->SetBranchAddress("cos_angle",&cos_angle);
     treeSCALE->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeSCALE->SetBranchAddress("dil_chosen",&dil_chosen);
     treeSCALE->SetBranchAddress("lep_pt",lep_pt);
     treeSCALE->SetBranchAddress("lep_id",lep_id);
     treeSCALE->SetBranchAddress("lep_eta",lep_eta);
     treeSCALE->SetBranchAddress("lep_phi",lep_phi);
     treeSCALE->SetBranchAddress("lep_dB",lep_dB);
     treeSCALE->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeSCALE->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
     treeSCALE->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeSCALE->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeSCALE->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeSCALE->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
     treeSCALE->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeSCALE->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeSCALE->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeSCALE->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeSCALE->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeSCALE->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeSCALE->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeSCALE->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeSCALE->SetBranchAddress("vertex_m",&vertex_m);
     treeSCALE->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeSCALE->SetBranchAddress("GoodVtx",&GoodVtx);
     treeSCALE->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeSCALE->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
     treeSCALE->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
     treeSCALE->SetBranchAddress("gen_lep_qOverPt", gen_lep_qOverPt);
     treeSCALE->SetBranchAddress("gen_lep_eta", gen_lep_eta);
     treeSCALE->SetBranchAddress("gen_lep_pt", gen_lep_pt);
		
	ne = treeSCALE->GetEntries();

	int p_chosen = -99;

	std::cout<<"START: "<<ne<<std::endl;

	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<10000; p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
		treeSCALE->GetEntry(p);
				
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
					treeSCALE->GetEntry(p_1);
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
						treeSCALE->GetEntry(p_1);
						n_event = event;
						
						if(p_1 > ne)
							break;

					}
				}

				treeSCALE->GetEntry(p);	

// 				std::cout<<"Selected = "<<p<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

// 				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				
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
				
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				
				MASS_SCALE = dil_mass;
				MASS_GEN = gen_dil_mass;

// 				res_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
// 				res_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
			
				if(fabs(gen_lep_eta[0])<1.2 && fabs(gen_lep_eta[1])<1.2){
					res_vs_mass_BB_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
				}
				else{
					res_vs_mass_BE_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
				}
	  

				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
	
				if(j >= 6 && j <=14){
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
	
		} // for event
	}

  for(int j=6; j < 15; j++){   

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

	std::cout<<"START: "<<ne<<std::endl;

	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<10000; p++){

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
						
						if(p_1 > ne)
							break;

					}
				}

				treeMC->GetEntry(p);	

// 				std::cout<<"Selected = "<<p<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

// 				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				
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
				
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 

				MASS = dil_mass;
				MASS_GEN = gen_dil_mass;
				
				weight[j] = 1;

// 				res_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
// 				res_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
			
				if(fabs(gen_lep_eta[0])<1.2 && fabs(gen_lep_eta[1])<1.2){
					res_vs_mass_BB_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				}
				else{
					res_vs_mass_BE_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				}
				
				if (fabs(gen_lep_eta[0]) < 0.9 && fabs(gen_lep_eta[1]) < 0.9){
					res_vs_mass_09->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				}
				else if (fabs(gen_lep_eta[0]) > 1.2 && fabs(gen_lep_eta[1]) > 1.2){
					res_vs_mass_EE->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				}
				else if((fabs(gen_lep_eta[0]) > 0.9 && fabs(gen_lep_eta[0]) < 1.2) || (fabs(gen_lep_eta[1]) > 0.9 && fabs(gen_lep_eta[1]) < 1.2)){
					res_vs_mass_OVER->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				}
				else{
					res_vs_mass_OTHER->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				}
		  		  		  

				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
	
				if(j >= 6 && j <=14){
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
	
		} // for event
	}



	gROOT->Reset();
	gROOT->SetBatch();

	ExtractResolution(res_vs_mass_BB_NS, "Res_BB_NS", "Res_BB_NS", resolution_BB_NS, mean_BB_NS);
	ExtractResolution(res_vs_mass_BB_S, "Res_BB_S", "Res_BB_S", resolution_BB_S, mean_BB_S);
	ExtractResolution(res_vs_mass_BE_NS, "Res_BE_NS", "Res_BE_NS", resolution_BE_NS, mean_BE_NS);
	ExtractResolution(res_vs_mass_BE_S, "Res_BE_S", "Res_BE_S", resolution_BE_S, mean_BE_S);
	ExtractResolution(res_vs_mass_09, "Below_09", "Below_09", resolution_09, mean_09);
	ExtractResolution(res_vs_mass_EE, "EE", "EE", resolution_EE, mean_EE);
	ExtractResolution(res_vs_mass_OVER, "OVER", "OVER", resolution_OVER, mean_OVER);
	ExtractResolution(res_vs_mass_OTHER, "OTHER", "OTHER", resolution_OTHER, mean_OTHER);

	SalvaRatio("BB", resolution_BB_S, resolution_BB_NS, "./MassScale/MC/Resolution_W_WO_scale_BB");
	SalvaRatio("BE", resolution_BE_S, resolution_BE_NS, "./MassScale/MC/Resolution_W_WO_scale_BE");

/*
	TCanvas	*c1	= new TCanvas ("c1", "c1", 500, 500);
	resolution_BB_NS->SetTitle("BB category");
	resolution_BB_NS->SetLineColor(kBlue+1);
	resolution_BB_NS->SetMarkerStyle(20);
	resolution_BB_NS->SetMarkerSize(0.5);
	resolution_BB_NS->SetMarkerColor(kBlue+1);
	resolution_BB_NS->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_BB_NS->Draw();
	resolution_BB_S->SetTitle("");
	resolution_BB_S->SetLineColor(kRed+1);
	resolution_BB_S->SetMarkerStyle(20);
	resolution_BB_S->SetMarkerSize(0.5);
	resolution_BB_S->SetMarkerColor(kRed);
	resolution_BB_S->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_BB_S->Draw("same");
	TLegend	*legend_1	=	new	TLegend(0.15,0.7,0.35,0.85);
	legend_1->AddEntry(resolution_BB_NS, "Without Scale", "lep");
	legend_1->AddEntry(resolution_BB_S,"With scale", "lep");
	legend_1->Draw(); 
	c1->Print("./MassScale/MC/Resolution_W_WO_scale_BB.png");
	c1->Print("./MassScale/MC/Resolution_W_WO_scale_BB.pdf");

	TCanvas	*c2	= new TCanvas ("c2", "c2", 500, 500);
	resolution_BE_NS->SetTitle("BE+EE category");
	resolution_BE_NS->SetLineColor(kBlue+1);
	resolution_BE_NS->SetMarkerStyle(20);
	resolution_BE_NS->SetMarkerSize(0.5);
	resolution_BE_NS->SetMarkerColor(kBlue+1);
	resolution_BE_NS->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_BE_NS->Draw();
	resolution_BE_S->SetTitle("");
	resolution_BE_S->SetLineColor(kRed+1);
	resolution_BE_S->SetMarkerStyle(20);
	resolution_BE_S->SetMarkerSize(0.5);
	resolution_BE_S->SetMarkerColor(kRed+1);
	resolution_BE_S->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_BE_S->Draw("same");
	TLegend	*legend_2	=	new	TLegend(0.15,0.7,0.35,0.85);
	legend_2->AddEntry(resolution_BE_NS, "Without Scale", "lep");
	legend_2->AddEntry(resolution_BE_S,"With scale", "lep");
	legend_2->Draw(); 
	c2->Print("./MassScale/MC/Resolution_W_WO_scale_BE.png");
	c2->Print("./MassScale/MC/Resolution_W_WO_scale_BE.pdf");
	*/
	TCanvas	*c3	= new TCanvas ("c3", "c3", 500, 500);
	resolution_09->SetTitle("");
	resolution_09->SetLineColor(kRed+1);
	resolution_09->SetMarkerStyle(20);
	resolution_09->SetMarkerSize(0.5);
	resolution_09->SetMarkerColor(kRed+1);
	resolution_09->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_09->SetStats(0);
	resolution_09->Draw();
	resolution_OTHER->SetTitle("");
	resolution_OTHER->SetLineColor(kGreen+1);
	resolution_OTHER->SetMarkerStyle(20);
	resolution_OTHER->SetMarkerSize(0.5);
	resolution_OTHER->SetMarkerColor(kGreen+1);
	resolution_OTHER->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_OTHER->SetStats(0);
	resolution_OTHER->Draw("same");
	resolution_EE->SetTitle("");
	resolution_EE->SetLineColor(kBlue+1);
	resolution_EE->SetMarkerStyle(20);
	resolution_EE->SetMarkerSize(0.5);
	resolution_EE->SetMarkerColor(kBlue+1);
	resolution_EE->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_EE->SetStats(0);	
	resolution_EE->Draw("same");
	resolution_OVER->SetTitle("");
	resolution_OVER->SetLineColor(kBlack);
	resolution_OVER->SetMarkerStyle(20);
	resolution_OVER->SetMarkerSize(0.5);
	resolution_OVER->SetMarkerColor(kBlack);
	resolution_OVER->GetYaxis()->SetRangeUser(0,	0.15);
	resolution_OVER->SetStats(0);
	resolution_OVER->Draw("same");
	TLegend	*legend_3	=	new	TLegend(0.15,0.7,0.88,0.85);
	legend_3->AddEntry(resolution_09, "BB category: both #mu |#eta| < 0.9", "lep");
	legend_3->AddEntry(resolution_OVER,"OVERLAP category: at least one #mu 0.9 < |#eta| < 1.2", "lep");
	legend_3->AddEntry(resolution_EE,"EE category: both #mu |#eta| > 1.2", "lep");
	legend_3->AddEntry(resolution_OTHER,"BE category", "lep");
	legend_3->Draw(); 
	c3->Print("./MassScale/MC/Resolution.png");
	c3->Print("./MassScale/MC/Resolution.pdf");

}














void ExtractResolution(TH2F *h_2, TString NOME, TString save, TH1F* res, TH1F* mean){

	save = "./MassScale/MC/" + save + ".pdf";		
	
	for(int i=1; i <= h_2->GetNbinsX(); i++){
		
		NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
		TCanvas *canvas = new TCanvas(NOME, NOME, 500, 500);
		TH1D* proiezione = h_2->ProjectionY(NOME,i,i);
		fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
	 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();

// 		proiezione->Draw();
	
	 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
	 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());

	    proiezione->Fit("gaus","M0R+");
    
	    TF1* f1 = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
	    f1->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 1.4, 1.4);// #15, 0.001);             
	    f1->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
	    f1->SetLineColor(kBlue);
	    f1->SetLineWidth(2);
	    proiezione->Fit("cruijff","MR+"); 
	    TString title = Form("BB: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    
		res->SetBinContent(i, f1->GetParameter(2));
		res->SetBinError(i, f1->GetParError(2));

		mean->SetBinContent(i, f1->GetParameter(1));
		mean->SetBinError(i, f1->GetParError(1));
			
		std::cout<<" ---------------------------------------------------------------------------- "<<MASS_BINS[i-1]<<" "<<MASS_BINS[i]<<"      "<<f1->GetParameter(2)<<std::endl;


        if(i==1){
        	TString save_init = save + "[";
              canvas->Print(save_init);
        }

        canvas->Print(save);

        if(i==h_2->GetNbinsX()){
        	TString save_final = save + "]";
        	canvas->Print(save_final);
        }
              
        canvas->Write();
    }
}

void SalvaRatio(TString name, TH1F* h_SCALE, TH1F* h_NOSCALE, TString save){

	TCanvas *c1 = new TCanvas(name, name,  500, 500);
	
   	c1->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
	
	Double_t max = 0.1;
	Double_t min = 0.0;
	
	h_NOSCALE->SetLineColor(kBlue);
	h_NOSCALE->SetMarkerColor(kBlue);
	h_NOSCALE->SetTitle(name);
	h_NOSCALE->Draw();
	h_NOSCALE->SetStats(0);
	c1->Update();
  
	h_SCALE->SetLineColor(kRed);
	h_SCALE->SetMarkerColor(kRed);
	h_SCALE->SetTitle(name);
	h_SCALE->Draw("same");
	h_SCALE->SetStats(0);

	h_NOSCALE->GetYaxis()->SetTitle("Mass resolution");
	h_SCALE->GetYaxis()->SetTitle("Mass resolution");
	
	TH1F *ratio = (TH1F*) h_SCALE->Clone();
	c1->Update();
// 	gPad->Update();

	h_NOSCALE->GetYaxis()->SetTitleOffset(1.2);
	h_SCALE->GetYaxis()->SetTitleOffset(1.2);

	h_NOSCALE->SetMaximum(max);
	h_NOSCALE->SetMinimum(min);	
// 	h_NOSCALE->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
// 	h_SCALE->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);

	h_SCALE->SetMarkerStyle(3);
	h_SCALE->SetMarkerColor(kRed);
	h_SCALE->SetMarkerSize(0.5);

	TLegend *l1 = new TLegend(0.4,0.8,0.6,0.9);
	l1->AddEntry(h_NOSCALE, "Without Scale", "l");
	l1->AddEntry(h_SCALE, "With Scale", "l");
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
	
	ratio->Add(h_NOSCALE, -1);
	ratio->Divide(h_NOSCALE);
	ratio->GetYaxis()->SetRangeUser(-0.5, 0.5);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("(With / Without) - 1");
	ratio->GetYaxis()->SetTitleOffset(0.50);
	ratio->GetYaxis()->SetNdivisions(506); 
	ratio->GetXaxis()->SetTitle("m_{ll} (GeV)");
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
		
	TString save_name_pdf = save + ".pdf";
	TString save_name_png = save + ".png";
	
	c1->Print(save_name_pdf);
	c1->Print(save_name_png);
}
