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

void MassScale_DATA(){

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetLegendTextSize(0.03);
gROOT->Reset();
gROOT->SetBatch();
gROOT->LoadMacro("./EfficiencyResolutionFromMC/cruijff.C");

	    Double_t MASS_BINS[] = {50, 120, 200, 400, 800, 1400, 2300, 3500, 4500, 6000};
// 	    Double_t MASS_BINS[] = {0, 100, 200, 300, 400, 600, 800, 1000, 1400, 1800, 2200, 2800, 3400, 4000, 5000, 6000};
// 	    Double_t MASS_BINS[] = {0, 120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
	    Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;
	    
	      TString NOME;
          std::vector<float> SIGMA;
          std::vector<float> SIGMA_ERR;

	TH1F *h_scale = new TH1F("Dilepton mass W scale correction", "Dilepton mass W scale correction", 2500, 0, 2500);
	h_scale->GetYaxis()->SetTitle("Entries/20 GeV"); 
	h_scale->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV]");
	TH1F *h_Nscale = new TH1F("Dilepton mass W/O scale correction", "Dilepton mass W/O scale correction", 2500, 0, 2500);
	h_Nscale->GetYaxis()->SetTitle("Entries/20 GeV"); 
	 
	TH2F *res_be = new TH2F("Resolution BE + EE", "Resolution BE + EE", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_be->GetXaxis()->SetTitle("(Mass - Mass_{scale}) / Mass"); 
	res_be->GetYaxis()->SetTitle("Entries"); 
	res_be->SetTitle("BE + EE mass residuals");
	TH2F *res_bb = new TH2F("Resolution BB", "Resolution BB", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_bb->GetXaxis()->SetTitle("(Mass - Mass_{scale}) / Mass");
	res_bb->GetYaxis()->SetTitle("Entries"); 
	res_bb->SetTitle("BB mass residuals");
	TH2F *res = new TH2F("Resolution", "Resolution", binnum, MASS_BINS, 240, -0.3, 0.3);
	res->GetXaxis()->SetTitle("(Mass - Mass_{scale}) / Mass"); 
	res->GetYaxis()->SetTitle("Entries"); 
	res->SetTitle("Mass residuals");
	
	TH2F *res_be_NS = new TH2F("Resolution(NoScale) BE + EE", "Resolution(NoScale) BE + EE", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_be_NS->GetYaxis()->SetTitle("(Mass_{GEN} - Mass) / Mass_{GEN}"); 
	res_be_NS->GetXaxis()->SetTitle("Mass"); 
	res_be_NS->SetTitle("BE + EE mass residuals: No Mass scale wrt Gen");
	TH2F *res_bb_NS = new TH2F("Resolution(NoScale) BB", "Resolution(NoScale) BB", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_bb_NS->GetYaxis()->SetTitle("(Mass_{GEN} - Mass) / Mass_{GEN}"); 
	res_bb_NS->GetXaxis()->SetTitle("Mass"); 
	res_bb_NS->SetTitle("BB mass residuals: No Mass scale wrt Gen");
	TH2F *res_NS = new TH2F("Resolution(NoScale)", "Resolution(NoScale)", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_NS->GetYaxis()->SetTitle("(Mass_{GEN} - Mass) / Mass_{GEN}"); 
	res_NS->GetXaxis()->SetTitle("Mass"); 
	res_NS->SetTitle("Mass residuals: No Mass scale wrt Gen");
	
	TH2F *res_be_S = new TH2F("Resolution(Scale) BE + EE", "Resolution(Scale) BE + EE", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_be_S->GetYaxis()->SetTitle("(Mass_{GEN} - Mass_{scale}) / Mass_{GEN}"); 
	res_be_S->GetXaxis()->SetTitle("Mass"); 
	res_be_S->SetTitle("BE + EE mass residuals: Mass scale wrt Gen");
	TH2F *res_bb_S = new TH2F("Resolution(Scale) BB", "Resolution(Scale) BB", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_bb_S->GetYaxis()->SetTitle("(Mass_{GEN} - Mass_{scale}) / Mass_{GEN}"); 
	res_bb_S->GetXaxis()->SetTitle("Mass"); 
	res_bb_S->SetTitle("BB mass residuals: Mass scale wrt Gen");
	TH2F *res_S = new TH2F("Resolution(Scale)", "Resolution(Scale)", binnum, MASS_BINS, 240, -0.3, 0.3);
	res_S->GetYaxis()->SetTitle("(Mass_{GEN} - Mass_{scale}) / Mass_{GEN}"); 
	res_S->GetXaxis()->SetTitle("Mass"); 
	res_S->SetTitle("Mass residuals: Mass scale wrt Gen");
	
	
	TH1F  *resolution_bb_NS = new TH1F("Resolution vs Mass (BB, Noscale)", "Resolution vs Mass (BB, Noscale)", binnum, MASS_BINS);
	resolution_bb_NS->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	resolution_bb_NS->GetYaxis()->SetTitle("Mass resolution");
	TH1F  *resolution_bb_S = new TH1F("Resolution vs Mass (BB, Scale)", "Resolution vs Mass (BB, Scale)", binnum, MASS_BINS);
	resolution_bb_S->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	resolution_bb_S->GetYaxis()->SetTitle("Mass resolution");
	TH1F  *resolution_be_NS = new TH1F("Resolution vs Mass (BE + EE, Noscale)", "Resolution vs Mass (BE + EE, Noscale)", binnum, MASS_BINS);
	resolution_be_NS->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	resolution_be_NS->GetYaxis()->SetTitle("Mass resolution");
	TH1F  *resolution_be_S = new TH1F("Resolution vs Mass (BE + EE, Scale)", "Resolution vs Mass (BE + EE Scale)", binnum, MASS_BINS);
	resolution_be_S->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	resolution_be_S->GetYaxis()->SetTitle("Mass resolution");
	
	
	
	TH1F  *mean_bb_NS = new TH1F("Mean vs Mass (BB, Noscale)", "Mean vs Mass (BB, Noscale)", binnum, MASS_BINS);
	mean_bb_NS->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	mean_bb_NS->GetYaxis()->SetTitle("Mean");
	TH1F  *mean_bb_S = new TH1F("Mean vs Mass (BB, Scale)", "Mean vs Mass (BB, Scale)", binnum, MASS_BINS);
	mean_bb_S->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	mean_bb_S->GetYaxis()->SetTitle("Mean");
	TH1F  *mean_be_NS = new TH1F("Mean vs Mass (BE + EE, Noscale)", "Mean vs Mass (BE + EE, Noscale)", binnum, MASS_BINS);
	mean_be_NS->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	mean_be_NS->GetYaxis()->SetTitle("Mean");
	TH1F  *mean_be_S = new TH1F("Mean vs Mass (BE + EE, Scale)", "Mean vs Mass (BE + EE Scale)", binnum, MASS_BINS);
	mean_be_S->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	mean_be_S->GetYaxis()->SetTitle("Mean");

	
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
// 		weight[i] = 1;
	}
	
                        	
  for(int j=6; j < 15; j++){   


	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j );    

     TChain *treeDATA = new TChain("SimpleNtupler/t");

     treeDATA->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA/");        
            
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
	std::cout<<"START"<<std::endl;
	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<1000 ;p++){
		if(p % 100000 == 0) std::cout<<p<<std::endl;		

		pp = p+1;
// 		p = p + count;
		
		n_count = 0;
		p_count = 0;

		// next event
		
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

				c_daughter_0.SetPtEtaPhiM(lep_pt[0], lep_eta[0], lep_phi[0], 0.10566);
				c_daughter_1.SetPtEtaPhiM(lep_pt[1], lep_eta[1], lep_phi[1], 0.10566);
				c_base = c_daughter_0 + c_daughter_1;

				MASS = c_base.M();
				MASS_SCALE = dil_mass;

				res->Fill(MASS, (MASS - MASS_SCALE)/MASS););
			
				if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1])<1.2){
					res_bb->Fill(MASS, (MASS - MASS_SCALE)/MASS);
				}
				else
					res_be->Fill(MASS, (MASS - MASS_SCALE)/MASS);
				}

		} // if selection
	
	} // for event
	

  } // for samples
	
	std::cout<<"STOP"<<std::endl;
	
    
    TCanvas *canvas = new TCanvas("canvas", "canvas", 210,45,1000,700);
    float  min, max;
    float fattore = 1.5;
  
    for(int i=1; i <= res_bb->GetNbinsX(); i++){

          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BB_NS = res_bb->ProjectionY(NOME,i,i);
          
          mean_bb->SetBinContent(i, proiezione_BB_NS->GetMean());
          mean_bb->SetBinError(i, proiezione_BB_NS->GetMeanError());
	

        if(i==1)
              canvas->Print("./MassScale/DATA/BB.pdf[");

        canvas->Print("./MassScale/DATA/BB.pdf");

        if(i==res_bb->GetNbinsX())
              canvas->Print("./MassScale/DATA/BB.pdf]");

        canvas->Write();

      }

    for(int i=1; i <= res_be->GetNbinsX(); i++){

          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BE_NS = res_be->ProjectionY(NOME,i,i);
          
          mean_be->SetBinContent(i, proiezione_BE_NS->GetMean());
          mean_be->SetBinError(i, proiezione_BE_NS->GetMeanError());
	

        if(i==1)
              canvas->Print("./MassScale/DATA/BE.pdf[");

        canvas->Print("./MassScale/DATA/BE.pdf");

        if(i==res_be->GetNbinsX())
              canvas->Print("./MassScale/DATA/BE.pdf]");

        canvas->Write();

      }


	TCanvas *BB = new TCanvas("Mean: (Scale - NoScale) / NoScale; BB", "Mean: (Scale - NoScale) / NoScale; BB", 210,45,750,500);
	mean_bb->Draw();
	BB->Print("Mean_vs_mass_BB.pdf");
	BB->Print("Mean_vs_mass_BB.png");


	TCanvas *BE = new TCanvas("Mean: (Scale - NoScale) / NoScale; BE", "Mean: (Scale - NoScale) / NoScale; BE", 210,45,750,500);
	mean_bb->Draw();
	BE->Print("Mean_vs_mass_BE.pdf");
	BE->Print("Mean_vs_mass_BE.png");

		
}
