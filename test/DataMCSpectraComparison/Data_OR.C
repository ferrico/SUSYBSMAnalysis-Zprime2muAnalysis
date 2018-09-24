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


void Data_OR(){

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
	
	const int    NMBINS = 110;//100;
	const double MMIN = 60., MMAX = 3000;//2100.;
	double logMbins[NMBINS+1];

	for (int ibin = 0; ibin <= NMBINS; ibin++){
	   	logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
// 		std::cout<<ibin<<"\t"<<logMbins[ibin]<<std::endl; 	
	}

	TH1F* OR = new TH1F("WITH OR", "WITH OR", NMBINS, logMbins);
	TH1F* OLD = new TH1F("WITHOUT OR", "WITHOUT OR", NMBINS, logMbins);
	
	
     TChain *treeOR = new TChain("SimpleNtupler/t");

     treeOR->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA_WITH_OR/ana_datamc_data.root");        

     treeOR->SetBranchAddress("event",&event);
     treeOR->SetBranchAddress("run",&run);
     treeOR->SetBranchAddress("lumi",&lumi);
     treeOR->SetBranchAddress("dil_mass",&dil_mass);
     treeOR->SetBranchAddress("cos_angle",&cos_angle);
     treeOR->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeOR->SetBranchAddress("dil_chosen",&dil_chosen);
     treeOR->SetBranchAddress("lep_pt",lep_pt);
     treeOR->SetBranchAddress("lep_id",lep_id);
     treeOR->SetBranchAddress("lep_eta",lep_eta);
     treeOR->SetBranchAddress("lep_phi",lep_phi);
     treeOR->SetBranchAddress("lep_dB",lep_dB);
     treeOR->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeOR->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
     treeOR->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeOR->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeOR->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeOR->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
     treeOR->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeOR->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeOR->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeOR->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeOR->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeOR->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeOR->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeOR->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeOR->SetBranchAddress("vertex_m",&vertex_m);
     treeOR->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeOR->SetBranchAddress("GoodVtx",&GoodVtx);
     treeOR->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeOR->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
		
	ne = treeOR->GetEntries();

	int p_chosen = -99;

	std::cout<<"START: "<<ne<<std::endl;

	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<1000 ;p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
		treeOR->GetEntry(p);
				
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
				bool trovato = false;
				int successivo = 0;

				if(dil_mass > - 99){			

					c_event = event;
					c_mass = dil_mass;
					int p_1 = p + 1;
// 					std::cout<<"A = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
					treeOR->GetEntry(p_1);
					n_event = event;
					while(c_event == n_event){	
						if(c_mass > 70 && c_mass < 110){
							successivo = 1;
							break;
						}
// 						std::cout<<"B = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
						if(
							GoodVtx && 
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
						treeOR->GetEntry(p_1);
						n_event = event;
						
						if(p_1 > ne)
							break;

					}
				}

				treeOR->GetEntry(p);	

// 				std::cout<<"Selected = "<<p<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

// 				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				

				
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////

				OR->Fill(dil_mass);

				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////

				
				if(!trovato) p += counting;
				p += successivo;

	       }//end of condition on event
	
		} // for event
		
		
		
     TChain *treeOLD = new TChain("SimpleNtupler/t");

     treeOLD->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA/ana_datamc_data.root");        

     treeOLD->SetBranchAddress("event",&event);
     treeOLD->SetBranchAddress("run",&run);
     treeOLD->SetBranchAddress("lumi",&lumi);
     treeOLD->SetBranchAddress("dil_mass",&dil_mass);
     treeOLD->SetBranchAddress("cos_angle",&cos_angle);
     treeOLD->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeOLD->SetBranchAddress("dil_chosen",&dil_chosen);
     treeOLD->SetBranchAddress("lep_pt",lep_pt);
     treeOLD->SetBranchAddress("lep_id",lep_id);
     treeOLD->SetBranchAddress("lep_eta",lep_eta);
     treeOLD->SetBranchAddress("lep_phi",lep_phi);
     treeOLD->SetBranchAddress("lep_dB",lep_dB);
     treeOLD->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeOLD->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
     treeOLD->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeOLD->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeOLD->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeOLD->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
     treeOLD->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeOLD->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeOLD->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeOLD->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeOLD->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeOLD->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeOLD->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeOLD->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeOLD->SetBranchAddress("vertex_m",&vertex_m);
     treeOLD->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeOLD->SetBranchAddress("GoodVtx",&GoodVtx);
     treeOLD->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeOLD->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
		
	ne = treeOLD->GetEntries();

	p_chosen = -99;

	std::cout<<"START: "<<ne<<std::endl;

	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<1000 ;p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
		treeOLD->GetEntry(p);
				
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
				bool trovato = false;
				int successivo = 0;

				if(dil_mass > - 99){			

					c_event = event;
					c_mass = dil_mass;
					int p_1 = p + 1;
// 					std::cout<<"A = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
					treeOLD->GetEntry(p_1);
					n_event = event;
					while(c_event == n_event){	
						if(c_mass > 70 && c_mass < 110){
							successivo = 1;
							break;
						}
// 						std::cout<<"B = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
						if(
							GoodVtx && 
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
						treeOLD->GetEntry(p_1);
						n_event = event;
						
						if(p_1 > ne)
							break;

					}
				}

				treeOLD->GetEntry(p);	

// 				std::cout<<"Selected = "<<p<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

// 				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				

				
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////

				OLD->Fill(dil_mass);

				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////

				
				if(!trovato) p += counting;
				p += successivo;

	       }//end of condition on event
	
		} // for event
		
		
		
		
	TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
   	c1->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();

	pad11->SetLogy();
	pad11->cd();

	float max;
	if(OR->GetMaximum() > OLD->GetMaximum())
		max = OR->GetMaximum();
	else
		max = OLD->GetMaximum();
	max = 1.1 * max;
	
	OR->SetLineColor(kBlue);
// 	OR->SetTitle(name);
	OR->Draw();
	OR->SetStats(1);
	c1->Update();

	TPaveStats * st_OR = (TPaveStats *)OR->GetListOfFunctions()->FindObject("stats");
    if( st_OR ){ 
		st_OR->SetName("Const");
		st_OR->SetX1NDC(0.75);
		st_OR->SetY1NDC(0.52);
		st_OR->SetY2NDC(0.72);
		st_OR->SetTextColor(kBlue);
    }
    else std::cout << "Null pointer to TPaveStats OR: " << st_OR << std::endl;
  
	OLD->SetLineColor(kRed);
// 	OLD->SetTitle(name);
	OLD->Draw();
// 	OLD->SetStats(0);
	TH1F *ratio = (TH1F*) OLD->Clone();
	OLD->SetStats(1);
	c1->Update();
// 	gPad->Update();

	TPaveStats * st_OLD = (TPaveStats *)OLD->GetListOfFunctions()->FindObject("stats");
    if( st_OLD ){ 
    	st_OLD->SetTextColor(kRed); 
		st_OLD->SetName("Const");
		st_OLD->SetX1NDC(0.75);
		st_OLD->SetY1NDC(0.75);
		st_OLD->SetY2NDC(0.95);

    }
    else std::cout << "Null pointer to TPaveStats OLD: " << st_OLD << std::endl;

	OR->GetYaxis()->SetTitleOffset(1.2);
	OLD->GetYaxis()->SetTitleOffset(1.2);

	OR->SetMaximum(max);
	OLD->SetMaximum(max);
// 	OR->SetMinimum(min);	
// 	OR->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
// 	OLD->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);

		
    c1->Update();

	OR->Draw();
	OLD->Draw("sameP");
	OLD->SetMarkerStyle(3);
	OLD->SetMarkerColor(kRed);
	OLD->SetMarkerSize(0.5);
// 	st_OR->Draw("same");
// 	st_OLD->Draw("same");
	
//     	
	TLegend *l1 = new TLegend(0.4,0.8,0.6,0.9);
	l1->AddEntry(OR, "Data W OR", "l");
	l1->AddEntry(OLD, "Data W/O OR", "l");
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

	ratio->Divide(OR);
	ratio->GetYaxis()->SetRangeUser(0, 2);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("(W/O OR) / (W OR)");
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
		
// 	c1->Print(save_name_pdf);
// 	c1->Print(save_name_png);
	c1->Print("Data_WITH_WITHOUT_OR.png");	
	c1->Print("Data_WITH_WITHOUT_OR.pdf");	
		
		
		
		
}
