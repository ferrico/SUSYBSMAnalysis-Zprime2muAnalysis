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

void DimuonMassRes_checkAtZ(){

gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetLegendTextSize(0.03);

gROOT->LoadMacro("cruijff.C+");

gROOT->Reset();
gROOT->SetBatch();

	    Double_t MASS_BINS[] = {72, 100, 152, 200, 300, 452, 800};

	    Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;
	    
          std::vector<float> SIGMA;
          std::vector<float> SIGMA_ERR;
		
     TString samples[39] =  {"dyInclusive50", 
	 						"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200", 
	 						"Wantitop", "tW", 
	 						"ttbar_lep50to500", "ttbar_lep_500to800", "ttbar_lep_800to1200", "ttbar_lep_1200to1800", "ttbar_lep1800toInf", 
	 						"Wjets", 
	 						"WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500", 
	 						"ZZ", "WZ", "ZZ_ext", "WZ_ext", 
	 						"dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000"};


	float events[39] = {19385554,  6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
						6933094, 6952830,
						79092400, 200000, 199800, 200000, 40829, 
						29705748,
						1999000, 200000, 200000, 200000, 38969, 
						990064, 1000000, 998034, 2995828,
						2977600, 100000, 100000, 98400, 100000, 100000, 100000, 100000, 100000
						};
						
	float sigma[39] = {6025.2, 2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
						35.6, 35.6,
						87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
						61526.7,
						12.178, 1.385, 0.0566, 0.0035, 0.00005,
						16.523, 47.13, 16.523, 47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7
						};

	float LUMINOSITY = 39484;

	float weight[39] = {0};
	
	Double_t BB_NS[15] = {0};
	Double_t BB_S[15] = {0};
	Double_t BB_NS_err[15] = {0};
	Double_t BB_S_err[15] = {0};

	Double_t BE_NS[15] = {0};
	Double_t BE_S[15] = {0};
	Double_t BE_NS_err[15] = {0};
	Double_t BE_S_err[15] = {0};
	
	Double_t BB[15] = {0};
	Double_t BE[15] = {0};
	Double_t BB_err[15] = {0};
	Double_t BE_err[15] = {0};
	
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
  float lep_glb_pt[2];
  float lep_picky_pt[2];
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
    
    float M[562242];
    float MS[562242];
    
    
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
	
			
	for(int i = 0; i < 39; i++){
// 		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight[i] = 1;
	}
	
	double mass;
	double rdil;
	TH2F* DileptonMass_2d_vsPt = new TH2F("DileptonMass_2d_vsPt", "DileptonMass_2d_vsPt", 200, 50., 150., 500., 0., 2000.);
	TH2F* DileptonMass_2d_vsPt_BB = new TH2F("DileptonMass_2d_vsPt: BB", "DileptonMass_2d_vsPt: BB", 200, 50., 150., 500., 0., 2000.);
	TH2F* DileptonMass_2d_vsPt_BE = new TH2F("DileptonMass_2d_vsPt: BE+EE", "DileptonMass_2d_vsPt: BE+EE", 200, 50., 150., 500., 0., 2000.);
	TH2F* DileptonMass_2d_vsPt_EE = new TH2F("DileptonMass_2d_vsPt: EE", "DileptonMass_2d_vsPt: EE", 200, 50., 150., 500., 0., 2000.);

	TH1F *res = new TH1F("Resolution", "Resolution", binnum, MASS_BINS);
	res->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res->GetYaxis()->SetTitle("Entries"); 
	res->SetTitle("Mass residuals");
	TH1F *res_bb = new TH1F("Resolution BB", "Resolution BB", binnum, MASS_BINS);
	res_bb->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res_bb->GetYaxis()->SetTitle("Entries"); 
	res_bb->SetTitle("BB mass residuals");
	TH1F *res_be = new TH1F("Resolution BE + EE", "Resolution BE + EE", binnum, MASS_BINS);
	res_be->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res_be->GetYaxis()->SetTitle("Entries"); 
	res_be->SetTitle("BE + EE mass residuals");
	TH1F *res_ee = new TH1F("Resolution EE", "Resolution EE", binnum, MASS_BINS);
	res_ee->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res_ee->GetYaxis()->SetTitle("Entries"); 
	res_ee->SetTitle("EE mass residuals");
	
	
                       	
// 	TFile *scale = new TFile("/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/MCMCSpectraComparison/MC/YesScale_YesEtaCut/Run2016MuonsOnly/", "READ");

     TChain *treeMC = new TChain("SimpleNtupler/t");
     treeMC->Add("../DataMCSpectraComparison/data/Run2017MuonsOnly/ana_datamc_data.root");
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
	std::cout<<"START"<<std::endl;
	for ( int p=0; p < ne ;p++){
// 	for ( int p=0; p<1000 ;p++){
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
		
			mass         = dil_mass; //vertex_m			
			DileptonMass_2d_vsPt->Fill(mass, lep_pt[0]);
		    DileptonMass_2d_vsPt->Fill(mass, lep_pt[1]);
  
		  if (abs(lep_eta[0]) < 1.2 && abs(lep_eta[1]) < 1.2) { 
	    	DileptonMass_2d_vsPt_BB->Fill(mass, lep_pt[0]);
	    	DileptonMass_2d_vsPt_BB->Fill(mass, lep_pt[1]);
		  }
		  else { 
    		DileptonMass_2d_vsPt_BE->Fill(mass, lep_pt[0]);
		    DileptonMass_2d_vsPt_BE->Fill(mass, lep_pt[1]);
		  }	
		  if((lep_eta[0] > 1.2 || lep_eta[0] < -1.2) && (lep_eta[1] > 1.2 || lep_eta[1] < -1.2)) {
		    DileptonMass_2d_vsPt_BE->Fill(mass, lep_pt[0]);
		    DileptonMass_2d_vsPt_BE->Fill(mass, lep_pt[1]);
		  }                                                                        

			
			
			
			
		} //if selection
	} // for entries

 

 for(int i = 1; i <= DileptonMass_2d_vsPt_BB->GetNbinsX(); i++){
 	
 	NOME = Form("DimuonMass Resolution BB: %.0f < m < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_BB->ProjectionY(NOME,i,i);
//  	proiezione->Draw();
 	
 	fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BB: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    proiezione->SetTitle(title);
        
    res_bb->SetBinContent(i, funct->GetParameter(2));
    res_bb->SetBinError(i, funct->GetParError(2));
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(-0.85, yPos, longstring);
    }
 	 	
//  	if(i==1)
//  		canvas->Print("./BB_resolution.pdf[");
// 	canvas->Print("./BB_resolution.pdf");
//  	if(i==DileptonMass_2d_vsPt_BB->GetNbinsX())
//  		canvas->Print("./BB_resolution.pdf]");
//     canvas->Write();

}

 for(int i = 1; i <= DileptonMass_2d_vsPt_BE->GetNbinsX(); i++){
 	
 	NOME = Form("DimuonMass Resolution BE+EE: %.0f < m < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_BE->ProjectionY(NOME,i,i);
//  	proiezione->Draw();
 	
 	fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BE+EE: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    proiezione->SetTitle(title); 
       
    res_be->SetBinContent(i, funct->GetParameter(2));
    res_be->SetBinError(i, funct->GetParError(2));
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(-0.85, yPos, longstring);
    }
 	 	
//  	if(i==1)
//  		canvas->Print("./BE_resolution.pdf[");
// 	canvas->Print("./BE_resolution.pdf");
//  	if(i==DileptonMass_2d_vsPt_BE->GetNbinsX())
//  		canvas->Print("./BE_resolution.pdf]");
//     canvas->Write();

}

 for(int i = 1; i <= DileptonMass_2d_vsPt_EE->GetNbinsX(); i++){
 	
 	NOME = Form("DimuonMass Resolution EE: %.0f < m < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_EE->ProjectionY(NOME,i,i);
//  	proiezione->Draw();
 	
 	fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("EE: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    proiezione->SetTitle(title);
    
    res_ee->SetBinContent(i, funct->GetParameter(2));
    res_ee->SetBinError(i, funct->GetParError(2));
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(-0.85, yPos, longstring);
    }
    
 	 	
//  	if(i==1)
//  		canvas->Print("./EE_resolution.pdf[");
// 	canvas->Print("./EE_resolution.pdf");
//  	if(i==DileptonMass_2d_vsPt_EE->GetNbinsX())
//  		canvas->Print("./EE_resolution.pdf]");
//     canvas->Write();

}

 
 TCanvas* c1 = new TCanvas("res", "res", 500, 500);
 c1->SetGrid();
 res_bb->SetTitle("");
 res_bb->SetLineColor(kRed+1);
 res_bb->SetMarkerStyle(20);
 res_bb->SetMarkerSize(0.5);
 res_bb->SetMarkerColor(kRed+1);
 res_bb->GetYaxis()->SetRangeUser(0, 0.15);
 res_bb->Draw();
 res_be->SetTitle("");
 res_be->SetLineColor(kGreen+1);
 res_be->SetMarkerStyle(20);
 res_be->SetMarkerSize(0.5);
 res_be->SetMarkerColor(kGreen+1);
 res_be->GetYaxis()->SetRangeUser(0, 0.15);
 res_be->Draw("same");
 res_ee->SetTitle("");
 res_ee->SetLineColor(kBlue+1);
 res_ee->SetMarkerStyle(20);
 res_ee->SetMarkerSize(0.5);
 res_ee->SetMarkerColor(kBlue+1);
 res_ee->GetYaxis()->SetRangeUser(0, 0.15);
 res_ee->Draw("same");
 TLegend *legend_3 = new TLegend(0.6,0.75,0.85,0.85);
 legend_3->AddEntry(res_bb, "BB category", "lep");
 legend_3->AddEntry(res_be,"BE+EE category", "lep");
 legend_3->AddEntry(res_ee,"EE category", "lep");
 legend_3->Draw(); 
//  c1->Print("Resolution.png");
//  c1->Print("Resolution.pdf");

//  DileptonMass_2d_vsPt->Draw();
//  TCanvas* c2 = new TCanvas("res_BB", "res_BB", 210,45,1050,750);
//  DileptonMass_2d_vsPt_BB->Draw();
//  TCanvas* c3 = new TCanvas("res_BE", "res_BE", 210,45,1050,750);
//  DileptonMass_2d_vsPt_BE->Draw();
//  TCanvas* c4 = new TCanvas("res_EE", "res_EE", 210,45,1050,750);
//  DileptonMass_2d_vsPt_EE->Draw();
 
 
} // main function



