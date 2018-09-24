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


void MassScale_DATA_v1(){

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetLegendTextSize(0.03);
		
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
 
    Double_t MASS_BINS[] = {50, 120, 200, 400, 800, 1400, 3500};
    Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;
 
	TH2F *res = new TH2F("Resolution", "Resolution", binnum, MASS_BINS, 240, -0.3, 0.3);
	res->GetXaxis()->SetTitle("(Mass - Mass_{scale}) / Mass"); 
	res->GetYaxis()->SetTitle("Entries"); 
	res->SetTitle("Mass residuals");
	TH2F *res_BB = new TH2F("Resolution BB", "Resolution BB", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_BB->GetXaxis()->SetTitle("(Mass - Mass_{scale}) / Mass");
	res_BB->GetYaxis()->SetTitle("Entries"); 
	res_BB->SetTitle("BB mass residuals");
	TH2F *res_BE = new TH2F("Resolution BE + EE", "Resolution BE + EE", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_BE->GetXaxis()->SetTitle("(Mass - Mass_{scale}) / Mass"); 
	res_BE->GetYaxis()->SetTitle("Entries"); 
	res_BE->SetTitle("BE + EE mass residuals"); 
                        	
     TChain *treeDATA = new TChain("SimpleNtupler/t");

     treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA_SCALE/ana_datamc_data.root");        

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
     treeDATA->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
     treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeDATA->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeDATA->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeDATA->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeDATA->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeDATA->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeDATA->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeDATA->SetBranchAddress("vertex_m",&vertex_m);
     treeDATA->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);
     treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
		
	ne = treeDATA->GetEntries();
	
	std::cout<<"START ---- "<<ne<<std::endl;

	int p_chosen = -99;

	std::cout<<"START: "<<ne<<std::endl;

	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<1000 ;p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
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

				int counting = 0;
				bool trovato = false;
				int successivo = 0;

				if(dil_mass > - 99){			

					c_event = event;
					c_mass = dil_mass;
					int p_1 = p + 1;
// 					std::cout<<"A = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
					treeDATA->GetEntry(p_1);
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
						treeDATA->GetEntry(p_1);
						n_event = event;
						
						if(p_1 > ne)
							break;

					}
				}

				treeDATA->GetEntry(p);	

// 				std::cout<<"Selected = "<<p<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

// 				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				
				
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////

				c_daughter_0.SetPtEtaPhiM(lep_pt[0], lep_eta[0], lep_phi[0], 0.10566);
				c_daughter_1.SetPtEtaPhiM(lep_pt[1], lep_eta[1], lep_phi[1], 0.10566);
				c_base = c_daughter_0 + c_daughter_1;

				MASS = c_base.M();
				MASS_SCALE = dil_mass;

				res->Fill(MASS, (MASS - MASS_SCALE)/MASS);
			
				if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1])<1.2){
					res_BB->Fill(MASS, (MASS - MASS_SCALE)/MASS);
				}
				else{
					res_BE->Fill(MASS, (MASS - MASS_SCALE)/MASS);
				}

				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				
				if(!trovato) p += counting;
				p += successivo;

	       }//end of condition on event
	
		} // for event


gROOT->Reset();
gROOT->SetBatch();

	TString NOME;
	
	TH1F* Mean_BB = new TH1F("Mean BB", "Mean BB", binnum, MASS_BINS);
	TH1F* Mean_BE = new TH1F("Mean BE+EE", "Mean BE+EE", binnum, MASS_BINS);
	
	TH1F* RMS_BB = new TH1F("RMS BB", "RMS BB", binnum, MASS_BINS);
	TH1F* RMS_BE = new TH1F("RMS BE+EE", "RMS BE+EE", binnum, MASS_BINS);

	
    for(int i=1; i <= res_BB->GetNbinsX(); i++){

          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BB = res_BB->ProjectionY(NOME,i,i);
          
          proiezione_BB->Draw();
          
          Mean_BB->SetBinContent(i, proiezione_BB->GetMean());
          Mean_BB->SetBinError(i, proiezione_BB->GetMeanError());

          RMS_BB->SetBinContent(i, proiezione_BB->GetRMS());
          RMS_BB->SetBinError(i, proiezione_BB->GetRMSError());
	

        if(i==1)
              canvas->Print("./MassScale/DATA/BB.pdf[");

        canvas->Print("./MassScale/DATA/BB.pdf");

        if(i==res_BB->GetNbinsX())
              canvas->Print("./MassScale/DATA/BB.pdf]");

        canvas->Write();

      }

    for(int i=1; i <= res_BE->GetNbinsX(); i++){

          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BE = res_BE->ProjectionY(NOME,i,i);

          proiezione_BE->Draw();
          
          Mean_BE->SetBinContent(i, proiezione_BE->GetMean());
          Mean_BE->SetBinError(i, proiezione_BE->GetMeanError());
          
          RMS_BE->SetBinContent(i, proiezione_BE->GetRMS());
          RMS_BE->SetBinError(i, proiezione_BE->GetRMSError());


        if(i==1)
              canvas->Print("./MassScale/DATA/BE.pdf[");

        canvas->Print("./MassScale/DATA/BE.pdf");

        if(i==res_BE->GetNbinsX())
              canvas->Print("./MassScale/DATA/BE.pdf]");

        canvas->Write();

      }

	Mean_BB->SetStats(0);
	Mean_BE->SetStats(0);
	RMS_BB->SetStats(0);
	RMS_BE->SetStats(0);

	TCanvas *BB = new TCanvas("Mean: (Scale - NoScale) / NoScale; BB", "Mean: (Scale - NoScale) / NoScale; BB", 210,45,750,500);
	Mean_BB->Draw();
	BB->Print("./MassScale/DATA/Mean_vs_mass_BB.pdf");
	BB->Print("./MassScale/DATA/Mean_vs_mass_BB.png");


	TCanvas *BE = new TCanvas("Mean: (Scale - NoScale) / NoScale; BE", "Mean: (Scale - NoScale) / NoScale; BE", 210,45,750,500);
	Mean_BE->Draw();
	BE->Print("./MassScale/DATA/Mean_vs_mass_BE.pdf");
	BE->Print("./MassScale/DATA/Mean_vs_mass_BE.png");

	TCanvas *BB_RMS = new TCanvas("RMS: (Scale - NoScale) / NoScale; BB", "RMS: (Scale - NoScale) / NoScale; BB", 210,45,750,500);
	RMS_BB->Draw();
	BB_RMS->Print("./MassScale/DATA/RMS_vs_mass_BB.pdf");
	BB_RMS->Print("./MassScale/DATA/RMS_vs_mass_BB.png");


	TCanvas *BE_RMS = new TCanvas("RMS: (Scale - NoScale) / NoScale; BE", "RMS: (Scale - NoScale) / NoScale; BE", 210,45,750,500);
	RMS_BE->Draw();
	BE_RMS->Print("./MassScale/DATA/RMS_vs_mass_BE.pdf");
	BE_RMS->Print("./MassScale/DATA/RMS_vs_mass_BE.png");







	
}
