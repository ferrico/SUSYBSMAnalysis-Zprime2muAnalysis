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

	double cruijff(double *x,double *par) {
// 		double arg = 0;
		double sigma;
		double alpha;
		
		if (x[0] >= par[0]){
			sigma = par[2];
			alpha = par[4];
		}
		else{
			sigma = par[1];
			alpha = par[3];
		}
		double delta_x = x[0] - par[0];
		double delta_x2 = delta_x * delta_x;
	
		return exp(- delta_x2 / (2 * sigma * sigma + alpha * delta_x2));
	
	}



void MassScale_MC(){

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetLegendTextSize(0.03);
gROOT->Reset();
gROOT->SetBatch();
// gROOT->LoadMacro("cruijff.C");


	    Double_t MASS_BIN[15] = {0};//{100, 200, 300, 400, 600, 800, 1000, 1400, 1800, 2200, 2800, 3400, 4000, 5000, 6000};

// 	    Double_t MASS_BINS[] = {50, 120, 200, 400, 800, 1400, 2300, 3500, 4500, 6000};
	    Double_t MASS_BINS[] = {0, 100, 200, 300, 400, 600, 800, 1000, 1400, 1800, 2200, 2800, 3400, 4000, 5000, 6000};
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



	
// 	TH1F *pt_plus_res_be_S = new TH1F("p_{T}+ resolution(Scale) BE + EE", "p_{T}+ resolution(Scale) BE + EE", 100, -0.2, 0.2);
// 	pt_plus_res_be_S->GetXaxis()->SetTitle("(pT_{GEN} - pT_{scale}) / pT_{GEN}"); 
// 	pt_plus_res_be_S->GetYaxis()->SetTitle("Entries"); 
// 	pt_plus_res_be_S->SetTitle("BE + EE p_{T}+ residuals: p_{T}+ scale wrt Gen");
// 	TH1F *pt_plus_res_bb_S = new TH1F("p_{T}+ resolution(Scale) BB", "p_{T}+ resolution(Scale) BB", 100, -0.2, 0.2);
// 	pt_plus_res_bb_S->GetXaxis()->SetTitle("(pT_{GEN} - pT_{scale}) / pT_{GEN}"); 
// 	pt_plus_res_bb_S->GetYaxis()->SetTitle("Entries"); 
// 	pt_plus_res_bb_S->SetTitle("BB mass residuals: p_{T}+ scale wrt Gen");
// 	TH1F *pt_plus_res_S = new TH1F("p_{T}+ resolution(Scale)", "p_{T}+ resolution(Scale)", 100, -0.2, 0.2);
// 	pt_plus_res_S->GetXaxis()->SetTitle("(pT_{GEN} - pT_{scale}) / pT_{GEN}"); 
// 	pt_plus_res_S->GetYaxis()->SetTitle("Entries"); 
// 	pt_plus_res_S->SetTitle("p_{T}+ residuals: p_{T}+ scale wrt Gen");
// 
// 	TH1F *pt_plus_res_be_NS = new TH1F("p_{T}+ resolution(NoScale) BE + EE", "p_{T}+ resolution(NoScale) BE + EE", 100, -0.2, 0.2);
// 	pt_plus_res_be_NS->GetXaxis()->SetTitle("(pT_{GEN} - pT) / pT_{GEN}"); 
// 	pt_plus_res_be_NS->GetYaxis()->SetTitle("Entries"); 
// 	pt_plus_res_be_NS->SetTitle("BE + EE p_{T}+ residuals: p_{T}+ NoScale wrt Gen");
// 	TH1F *pt_plus_res_bb_NS = new TH1F("p_{T}+ resolution(NoScale) BB", "p_{T}+ resolution(NoScale) BB", 100, -0.2, 0.2);
// 	pt_plus_res_bb_NS->GetXaxis()->SetTitle("(pT_{GEN} - pT) / pT_{GEN}"); 
// 	pt_plus_res_bb_NS->GetYaxis()->SetTitle("Entries"); 
// 	pt_plus_res_bb_NS->SetTitle("BB p_{T}+ residuals: p_{T}+ NoScale wrt Gen");
// 	TH1F *pt_plus_res_NS = new TH1F("p_{T}+ resolution(NoScale)", "p_{T}+ resolution(NoScale)", 100, -0.2, 0.2);
// 	pt_plus_res_NS->GetXaxis()->SetTitle("(pT_{GEN} - pT) / pT_{GEN}"); 
// 	pt_plus_res_NS->GetYaxis()->SetTitle("Entries"); 
// 	pt_plus_res_NS->SetTitle("p_{T}+ residuals: p_{T}+ NoScale wrt Gen");
// 	
// 	TH1F *pt_minus_res_be_S = new TH1F("p_{T}- resolution(Scale) BE + EE", "p_{T}- resolution(Scale) BE + EE", 100, -0.2, 0.2);
// 	pt_minus_res_be_S->GetXaxis()->SetTitle("(pT_{GEN} - pT_{scale}) / pT_{GEN}"); 
// 	pt_minus_res_be_S->GetYaxis()->SetTitle("Entries"); 
// 	pt_minus_res_be_S->SetTitle("BE + EE p_{T}- residuals: p_{T}- scale wrt Gen");
// 	TH1F *pt_minus_res_bb_S = new TH1F("p_{T}- resolution(Scale) BB", "p_{T}- resolution(Scale) BB", 100, -0.2, 0.2);
// 	pt_minus_res_bb_S->GetXaxis()->SetTitle("(pT_{GEN} - pT_{scale}) / pT_{GEN}"); 
// 	pt_minus_res_bb_S->GetYaxis()->SetTitle("Entries"); 
// 	pt_minus_res_bb_S->SetTitle("BB mass residuals: p_{T}- scale wrt Gen");
// 	TH1F *pt_minus_res_S = new TH1F("p_{T}- resolution(Scale)", "p_{T}- resolution(Scale)", 100, -0.2, 0.2);
// 	pt_minus_res_S->GetXaxis()->SetTitle("(pT_{GEN} - pT_{scale}) / pT_{GEN}"); 
// 	pt_minus_res_S->GetYaxis()->SetTitle("Entries"); 
// 	pt_minus_res_S->SetTitle("p_{T}- residuals: p_{T}- scale wrt Gen");
// 
// 	TH1F *pt_minus_res_be_NS = new TH1F("p_{T}- resolution(NoScale) BE + EE", "p_{T}- resolution(NoScale) BE + EE", 100, -0.2, 0.2);
// 	pt_minus_res_be_NS->GetXaxis()->SetTitle("(pT_{GEN} - pT) / pT_{GEN}"); 
// 	pt_minus_res_be_NS->GetYaxis()->SetTitle("Entries"); 
// 	pt_minus_res_be_NS->SetTitle("BE + EE p_{T}- residuals: p_{T}- NoScale wrt Gen");
// 	TH1F *pt_minus_res_bb_NS = new TH1F("p_{T}- resolution(NoScale) BB", "p_{T}- resolution(NoScale) BB", 100, -0.2, 0.2);
// 	pt_minus_res_bb_NS->GetXaxis()->SetTitle("(pT_{GEN} - pT) / pT_{GEN}"); 
// 	pt_minus_res_bb_NS->GetYaxis()->SetTitle("Entries"); 
// 	pt_minus_res_bb_NS->SetTitle("BB p_{T}- residuals: p_{T}- NoScale wrt Gen");
// 	TH1F *pt_minus_res_NS = new TH1F("p_{T}- resolution(NoScale)", "p_{T}- resolution(NoScale)", 100, -0.2, 0.2);
// 	pt_minus_res_NS->GetXaxis()->SetTitle("(pT_{GEN} - pT) / pT_{GEN}"); 
// 	pt_minus_res_NS->GetYaxis()->SetTitle("Entries"); 
// 	pt_minus_res_NS->SetTitle("p_{T}- residuals: p_{T}- NoScale wrt Gen");

	
	
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
	float LUMINOSITY = 36295.39;

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
	
	bool reco;
	bool noIB;
	
	reco = false;
	noIB = false;
			
	for(int i = 0; i < 39; i++){
// 		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight[i] = 1;
	}
	
                        	
// 	TFile *scale = new TFile("/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/MCMCSpectraComparison/MC/YesScale_YesEtaCut/Run2016MuonsOnly/", "READ");
  for(int j=30; j < 39; j++){   

// 	if(j != 35) continue;  // 30 per dy50to120 ---- 33 per dy400to800 ---- 35 per dy1400to2300 //

	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    
     TChain *treeMC = new TChain("SimpleNtupler/t");
     if(j == 35)
     	TChain *treeMC = new TChain("SimpleNtupler/t;3");
//      treeMC->Add("mc/mc_YesEtaCut_NoScale/ana_datamc_"+samples[j]+".root");
//       treeMC->Add("mc/YesScale/ana_datamc_"+samples[j]+".root");
    if(noIB)
      treeMC->Add("mc/DY_scale_genInfo_NOIB/ana_datamc_"+samples[j]+".root");
    else
      treeMC->Add("mc/DY_scale_genInfo_status1/ana_datamc_"+samples[j]+".root");
      
      

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
	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<1000 ;p++){
		if(p % 100000 == 0) std::cout<<p<<std::endl;		

		pp = p+1;
// 		p = p + count;
		
		n_count = 0;
		p_count = 0;

		// next event
		
		treeMC->GetEntry(p);
// 		std::cout<<p<<"   BASE       "<<run<<"   "<<lumi<<"   "<<event<<"   "<<dil_mass<<std::endl;
		ok = false;
		
		next = false;
		
		if(fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20){
// 		if(dil_mass > 900 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20){

			treeMC->GetEntry(p);
			c_daughter_0.SetPtEtaPhiM(lep_pt[0], lep_eta[0], lep_phi[0], 0.10566);
			c_daughter_1.SetPtEtaPhiM(lep_pt[1], lep_eta[1], lep_phi[1], 0.10566);
			c_base = c_daughter_0 + c_daughter_1;
			c_run = run;
			c_lumi = lumi;
			c_event = event;
			if(reco){
				gen_eta_first = lep_eta[0];
				gen_eta_second = lep_eta[1];
			}
			else{
			gen_eta_first = gen_lep_eta[0];
			gen_eta_second = gen_lep_eta[1];
			}			
			pt_first = lep_pt[0];
			pt_second = lep_pt[1];
			gen_pt_first = gen_lep_pt[0];
			gen_pt_second = gen_lep_pt[1];
			max_pt = pt_first+pt_second;
			MASS = c_base.M();
			MASS_SCALE = dil_mass;
			MASS_GEN = gen_dil_mass;
// 			if(lep_id[0] > 0){
// 				pt_minus = lep_pt[0];
// 				pt_plus = lep_pt[1];
// 			}
// 			else{
// 				pt_minus = lep_pt[1];
// 				pt_plus = lep_pt[0];
// 			}
// 			if(gen_lep_qOverPt[0] > 0){
// 				gen_pt_plus	= 1/gen_lep_qOverPt[0];
// 				gen_pt_minus	= 1/gen_lep_qOverPt[1];
// 			}
// 			else{
// 				gen_pt_plus	= 1/gen_lep_qOverPt[1];
// 				gen_pt_minus	= 1/gen_lep_qOverPt[0];
// 			}			

			
			if(p < ne-1){
				treeMC->GetEntry(pp);				
				int ddd = 0;
				while(c_run == run && c_lumi == lumi && c_event == event){
// 					if(ddd > 0) std::cout<<ddd<<std::endl;
					if(fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20){
						n_daughter_0.SetPtEtaPhiM(lep_pt[0], lep_eta[0], lep_phi[0], 0.10566);
						n_daughter_1.SetPtEtaPhiM(lep_pt[1], lep_eta[1], lep_phi[1], 0.10566);
						n_base = n_daughter_0 + n_daughter_1;
// 						if(MASS_GEN < gen_dil_mass){
						if(max_pt < (lep_pt[0]+lep_pt[1])){
							MASS = n_base.M();
							MASS_SCALE = dil_mass;
							MASS_GEN = gen_dil_mass;
							if(reco){
								gen_eta_first = lep_eta[0];
								gen_eta_second = lep_eta[1];
							}
							else{
								gen_eta_first = gen_lep_eta[0];
								gen_eta_second = gen_lep_eta[1];
							}	
							pt_first = lep_pt[0];
							pt_second = lep_pt[1];
							max_pt = pt_first + pt_second;
							gen_pt_first = gen_lep_pt[0];
							gen_pt_second = gen_lep_pt[1];
// 							if(lep_id[0] > 0){
// 								pt_minus = lep_pt[0];
// 								pt_plus = lep_pt[1];
// 							}
// 							else{
// 								pt_minus = lep_pt[1];
// 								pt_plus = lep_pt[0];
// 							}
// 							if(gen_lep_qOverPt[0] > 0){
// 								gen_pt_plus	= 1/gen_lep_qOverPt[0];
// 								gen_pt_minus	= 1/gen_lep_qOverPt[1];
// 							}
// 							else{
// 								gen_pt_plus	= 1/gen_lep_qOverPt[1];
// 								gen_pt_minus	= 1/gen_lep_qOverPt[0];
// 							}	
						}
						n_doppioni++;
					}
						
					pp++;
					p++;
					ddd++;
					treeMC->GetEntry(pp);					
				}
			}				

// 			std::cout<<p<<") "<<c_run<<"   "<<c_lumi<<"   "<<c_event<<"   "<<MASS<<"   "<<MASS_SCALE<<"   "<<std::endl;
// 			std::cout<<p<<") "<<run<<"   "<<lumi<<"   "<<event<<"   "<<MASS<<"   "<<MASS_SCALE<<"   "<<p_count<<"   "<<n_count<<std::endl;

			res->Fill(MASS, (MASS - MASS_SCALE)/MASS, weight[j]);
			res_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
			res_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
			
			if(fabs(gen_eta_first)>1.2 || fabs(gen_eta_second)>1.2){
				res_be->Fill(MASS, (MASS - MASS_SCALE)/MASS, weight[j]);
				res_be_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				res_be_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
			}
			if(fabs(gen_eta_first)<1.2 && fabs(gen_eta_second)<1.2){
				res_bb->Fill(MASS, (MASS - MASS_SCALE)/MASS, weight[j]);
				res_bb_NS->Fill(MASS_GEN, (MASS_GEN - MASS)/MASS_GEN, weight[j]);
				res_bb_S->Fill(MASS_GEN, (MASS_GEN - MASS_SCALE)/MASS_GEN, weight[j]);
			}

// 			res->Fill(MASS, (MASS - MASS_SCALE)/MASS);
// 			res_NS->Fill(MASS, (MASS_GEN - MASS)/MASS_GEN);
// 			res_S->Fill(MASS_SCALE, (MASS_GEN - MASS_SCALE)/MASS_GEN);
// 			
// 			if(fabs(gen_eta_first)>1.2 || fabs(gen_eta_second)>1.2){
// 				res_be->Fill(MASS, (MASS - MASS_SCALE)/MASS);
// 				res_be_NS->Fill(MASS, (MASS_GEN - MASS)/MASS_GEN);
// 				res_be_S->Fill(MASS_SCALE, (MASS_GEN - MASS_SCALE)/MASS_GEN);
// 			}
// 			if(fabs(gen_eta_first)<1.2 && fabs(gen_eta_second)<1.2){
// 				res_bb->Fill(MASS, (MASS - MASS_SCALE)/MASS);
// 				res_bb_NS->Fill(MASS, (MASS_GEN - MASS)/MASS_GEN);
// 				res_bb_S->Fill(MASS_SCALE, (MASS_GEN - MASS_SCALE)/MASS_GEN);
// 			}


			

		} // if selection
	
	} // for event
	
	std::cout<<"number of samples: "<<j-29<<std::endl;
// 	std::cout<<res_bb_NS->GetStdDev()<<std::endl;
// 	std::cout<<res_bb_S->GetStdDev()<<std::endl;
// 	std::cout<<res_be_NS->GetStdDev()<<std::endl;
// 	std::cout<<res_be_S->GetStdDev()<<std::endl;
	

//     float  sigma, alpha;
       
// 	res_bb_NS->Reset();
// 	res_bb_S->Reset();
// 	res_be_NS->Reset();
// 	res_be_S->Reset();
// 	
// 	 

  } // for samples
	
	std::cout<<"STOP"<<std::endl;
	
	std::cout<<"Doppioni = "<<p_doppioni<<" --- "<<n_doppioni<<std::endl;
		
    
    TCanvas *canvas = new TCanvas("canvas", "canvas", 210,45,1000,700);
//    	TString nome_file;   	
    float  min, max;
//     TF1 *func;
    float fattore = 1.5;
  
    for(int i=1; i <= res_bb_NS->GetNbinsX(); i++){

          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BB_NS = res_bb_NS->ProjectionY(NOME,i,i);
//           resolution_bb_NS->GetXaxis()->SetTitle("Mass resolution");
//           resolution_bb_NS->GetYaxis()->SetTitle("entries");
//           resolution_bb_NS->SetTitle("Z' Mass Res (Z' mass function) - GEM");
//           MEAN_PT.push_back(resolution_bb_NS->GetMean());
//           RMS_PT.push_back(resolution_bb_NS->GetRMS());
//           MIN = MEAN_PT.at(i-1) - 1.5 * RMS_PT.at(i-1);
//           MAX = MEAN_PT.at(i-1) + 1.5 * RMS_PT.at(i-1);

           min = -fattore * proiezione_BB_NS->GetRMS();
           max = fattore * proiezione_BB_NS->GetRMS();
          TF1 *f1 = new TF1("f1","gaus",min,max);
          proiezione_BB_NS->Fit("f1","R");
//           RMS_RES_PT_GEM.push_back(resolution_bb_NS->GetRMS());
//           RMS_RES_PT_ERR_GEM.push_back(resolution_bb_NS->GetRMSError());
		resolution_bb_NS->SetBinContent(i, f1->GetParameter(2));
		resolution_bb_NS->SetBinError(i, f1->GetParError(2));

		mean_bb_NS->SetBinContent(i, f1->GetParameter(1));
		mean_bb_NS->SetBinError(i, f1->GetParError(1));

		
		BB_NS[i-1] = f1->GetParameter(2);
		BB_NS_err[i-1] = f1->GetParError(2);
		
		std::cout<<" ---------------------------------------------------------------------------- BB NO SCALE"<<MASS_BINS[i-1]<<" "<<MASS_BINS[i]<<"      "<<f1->GetParameter(2)<<std::endl;
		
//           PT_GEM.push_back(res_bb_NS->GetXaxis()->GetBinCenter(i));
//           PT_ERR_GEM.push_back((MASS_BINS[i]-MASS_BINS[i-1])/2);

          if(i==1)
              canvas->Print("./MassScale/MC/BB_NS.pdf[");

              canvas->Print("./MassScale/MC/BB_NS.pdf");

          if(i==res_bb_NS->GetNbinsX())
              canvas->Print("./MassScale/MC/BB_NS.pdf]");

          canvas->Write();

      }

    for(int i=1; i <= res_bb_S->GetNbinsX(); i++){
//           cout<<" ---------------------------------------------------------------------------- "<<i<<endl;
          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BB_S = res_bb_S->ProjectionY(NOME,i,i);

           min = -fattore * proiezione_BB_S->GetRMS();
           max = fattore * proiezione_BB_S->GetRMS();
          TF1 *f1 = new TF1("f1","gaus",min,max);
          proiezione_BB_S->Fit("f1","R");
		resolution_bb_S->SetBinContent(i, f1->GetParameter(2));
		resolution_bb_S->SetBinError(i, f1->GetParError(2));

		mean_bb_S->SetBinContent(i, f1->GetParameter(1));
		mean_bb_S->SetBinError(i, f1->GetParError(1));
		
		BB_S[i-1] = f1->GetParameter(2);
		BB_S_err[i-1] = f1->GetParError(2);
         
         		std::cout<<" ---------------------------------------------------------------------------- BE SCALE"<<MASS_BINS[i-1]<<" "<<MASS_BINS[i]<<"      "<<f1->GetParameter(2)<<std::endl;
         
          if(i==1)
              canvas->Print("./MassScale/MC/BB_S.pdf[");

              canvas->Print("./MassScale/MC/BB_S.pdf");

          if(i==res_bb_S->GetNbinsX())
              canvas->Print("./MassScale/MC/BB_S.pdf]");

          canvas->Write();

      }
      
      
    for(int i=1; i <= res_be_NS->GetNbinsX(); i++){
//           cout<<" ---------------------------------------------------------------------------- "<<i<<endl;
          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BE_NS = res_be_NS->ProjectionY(NOME,i,i);

           min = -fattore * proiezione_BE_NS->GetRMS();
           max = fattore * proiezione_BE_NS->GetRMS();
          TF1 *f1 = new TF1("f1","gaus",min,max);
          proiezione_BE_NS->Fit("f1","R");

		resolution_be_NS->SetBinContent(i, f1->GetParameter(2));
		resolution_be_NS->SetBinError(i, f1->GetParError(2));

		mean_be_NS->SetBinContent(i, f1->GetParameter(1));
		mean_be_NS->SetBinError(i, f1->GetParError(1));

		BE_NS[i-1] = f1->GetParameter(2);
		BE_NS_err[i-1] = f1->GetParError(2);
		
          		std::cout<<" ---------------------------------------------------------------------------- BE NO SCALE"<<MASS_BINS[i-1]<<" "<<MASS_BINS[i]<<"      "<<f1->GetParameter(2)<<std::endl;
          
          if(i==1)
              canvas->Print("./MassScale/MC/BE_NS.pdf[");

              canvas->Print("./MassScale/MC/BE_NS.pdf");

          if(i==res_be_NS->GetNbinsX())
              canvas->Print("./MassScale/MC/BE_NS.pdf]");

          canvas->Write();

      }

    for(int i=1; i <= res_be_S->GetNbinsX(); i++){
//           cout<<" ---------------------------------------------------------------------------- "<<i<<endl;
          NOME = Form("Mass Resolution: %.0f < mass < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
          TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
          TH1D* proiezione_BE_S = res_be_S->ProjectionY(NOME,i,i);
//           cout<<" ------------------------------------------------------------------------------------------ Integral =  "<<proiezione_BE_S->Integral()<<" --------------- "<<res_be_S->Integral(i-1, i)<<endl;
           min = -fattore * proiezione_BE_S->GetRMS();
           max = fattore * proiezione_BE_S->GetRMS();
          TF1 *f1 = new TF1("f1","gaus",min,max);
          proiezione_BE_S->Fit("f1","R");
//           if(noIB){
//           		std::cout<<" -----------------------------------------------------------------------------------------cruijff"<<std::endl;
// 	          TF1 *f2 = new TF1("f2",cruijff, -0.3, 0.3, 5);
// 	          f2->SetParNames("mean", "#sigma_{R}", "#alpha_{R}", "#sigma_{L}", "#alpha_{L}");
// 	          f2->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), 0., 0.);
// 	          proiezione_BE_S->Fit("f2","R");
// 	          resolution_be_S->SetBinContent(i, f2->GetParameter(2));
// 	          resolution_be_S->SetBinError(i, f2->GetParError(2));
// 	          std::cout<<f2->GetParameter(0)<<std::endl;
// 	          std::cout<<f2->GetParameter(1)<<std::endl;
// 	          std::cout<<f2->GetParameter(2)<<std::endl;
// 	          std::cout<<f2->GetParameter(3)<<std::endl;
// 	          std::cout<<f2->GetParameter(4)<<std::endl;
//           }
		resolution_be_S->SetBinContent(i, f1->GetParameter(2));
		resolution_be_S->SetBinError(i, f1->GetParError(2));
		
		mean_be_S->SetBinContent(i, f1->GetParameter(1));
		mean_be_S->SetBinError(i, f1->GetParError(1));
		
		MASS_BIN[i] = res_be_S->GetXaxis()->GetBinCenter(i+1);
		
				BE_S[i-1] = f1->GetParameter(2);
				BE_S_err[i-1] = f1->GetParError(2);
          
          		std::cout<<" ---------------------------------------------------------------------------- BE SCALE"<<MASS_BINS[i-1]<<" "<<MASS_BINS[i]<<"      "<<f1->GetParameter(2)<<std::endl;
          
          if(i==1)
              canvas->Print("./MassScale/MC/BE_S.pdf[");

              canvas->Print("./MassScale/MC/BE_S.pdf");

          if(i==res_be_S->GetNbinsX())
              canvas->Print("./MassScale/MC/BE_S.pdf]");

          canvas->Write();

      }     


	Double_t MASS_BINS_err[15] = {0};
	MASS_BIN[0] = 50;
	for(int i = 1; i < 16; i++){
			MASS_BINS_err[i-1] = (MASS_BINS[i] - MASS_BINS[i-1])/2;
// 			std::cout<<MASS_BINS_err[i]<<std::endl;
	}

	
	for(i = 0; i < 15; i++){
		BB[i] = BB_S[i] / BB_NS[i] - 1;
		BE[i] = BE_S[i] / BE_NS[i] - 1;
// 		if(i < 14)
			std::cout<<i<<"   centro = "<<MASS_BIN[i]<<"   estremi ="<<MASS_BINS[i]<<" "<<MASS_BINS[i+1]<<" errore = "<<MASS_BINS_err[i]<<std::endl;
		std::cout<<i<<"  "<<BB_S[i]<<"   "<<BB_NS[i]<<"  "<<BB[i]<<std::endl;
		std::cout<<i<<"  "<<BE_S[i]<<"   "<<BE_NS[i]<<"  "<<BE[i]<<std::endl;
		
		BB_err[i] = sqrt( (BB_S_err[i]/(BB_NS[i])) * (BB_S_err[i]/(BB_NS[i])) + (BB_S[i] * BB_NS_err[i] / (BB_NS[i] * BB_NS[i])) * (BB_S[i] * BB_NS_err[i] / (BB_NS[i] * BB_NS[i])));
		BE_err[i] = sqrt( (BE_S_err[i]/(BE_NS[i])) * (BE_S_err[i]/(BE_NS[i])) + (BE_S[i] * BE_NS_err[i] / (BE_NS[i] * BE_NS[i])) * (BE_S[i] * BE_NS_err[i] / (BE_NS[i] * BE_NS[i])));
	}



	TCanvas *bb = new TCanvas("Ratio: NoScale / Scale; BB", "Ratio: NoScale / Scale; BB", 210,45,750,500);
	TPad *pad1 = new TPad("pad1","This is pad1", 0.05, 0.33, 1, 1.0);
	TPad *pad2 = new TPad("pad2","This is pad2", 0.05, 0.05, 1, 0.33);

	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->Draw();

	pad1->cd();
	resolution_bb_S->Draw("E");
	resolution_bb_NS->Draw("E same");
	resolution_bb_S->SetTitle("Mass resolution vs mass: BB");
	resolution_bb_NS->SetTitle("Mass resolution vs mass: BB");
	resolution_bb_NS->GetYaxis()->SetTitleSize(10);
	resolution_bb_NS->SetMarkerStyle(20);
	resolution_bb_S->SetMarkerStyle(20);
	resolution_bb_NS->SetMarkerSize(0.75);
	resolution_bb_S->SetMarkerSize(0.75);
	resolution_bb_S->SetMarkerColor(kRed);
	resolution_bb_S->SetLineColor(kRed);
	resolution_bb_NS->SetMarkerColor(kBlue);
	resolution_bb_NS->SetLineColor(kBlue);
	resolution_bb_NS->SetStats(0);
	resolution_bb_S->SetStats(0);
// 	resolution_bb_S->GetYaxis()->SetLabelSize(0.1);
// 	resolution_bb_S->GetYaxis()->SetTitleSize(15);
 	TLegend *legend_3 = new TLegend(0.25,0.8,0.5,0.9);
 	legend_3->AddEntry(resolution_bb_NS, "No Scale", "lep");
 	legend_3->AddEntry(resolution_bb_S,"Scale", "lep");
 	legend_3->Draw();

	bb->cd();
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw(); 	
 	pad2->cd();
	gr_bb = new TGraphErrors(15, MASS_BIN, BB, MASS_BINS_err, BB_err);
// 	gr_bb->SetTitle("Ratio: NoScale / Scale - BB");
	gr_bb->SetTitle(" ");
	gr_bb->SetMarkerStyle(21);
	gr_bb->SetMarkerSize(0.5);
	gr_bb->GetYaxis()->SetTitleSize(12);
	gr_bb->GetYaxis()->SetTitleFont(43);
	gr_bb->GetYaxis()->SetLabelSize(0.1);
	gr_bb->GetXaxis()->SetLabelSize(0.07);	
	gr_bb->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	gr_bb->GetXaxis()->SetTitleSize(12);
	gr_bb->GetXaxis()->SetTitleFont(43);
	gr_bb->GetXaxis()->SetRangeUser(0, 6000);

	gr_bb->GetYaxis()->SetTitle("(Scale - NoScale) / NoScale");
	gr_bb->Draw("AE*");
	
	bb->Print("./MassScale/MC/BB_vs_Mass.pdf");
	bb->Print("./MassScale/MC/BB_vs_Mass.png");
	
	
	TCanvas *be = new TCanvas("Ratio: NoScale / Scale; BE", "Ratio: NoScale / Scale; BE", 210,45,750,500);
	TPad *pad11 = new TPad("pad11","This is pad11", 0.05, 0.33, 1, 1.0);
	TPad *pad22 = new TPad("pad22","This is pad22", 0.05, 0.05, 1, 0.33);

	pad11->SetBottomMargin(0);
	pad11->SetGridx();
	pad11->Draw();

	pad11->cd();
	resolution_be_S->Draw("E");
	resolution_be_NS->Draw("E same");
	resolution_be_S->SetTitle("Mass resolution vs mass: BE");
	resolution_be_NS->SetTitle("Mass resolution vs mass: BE");
	resolution_be_NS->GetYaxis()->SetTitleSize(10);
	resolution_be_NS->SetMarkerStyle(20);
	resolution_be_S->SetMarkerStyle(20);
	resolution_be_NS->SetMarkerSize(0.75);
	resolution_be_S->SetMarkerSize(0.75);
	resolution_be_S->SetMarkerColor(kRed);
	resolution_be_S->SetLineColor(kRed);
	resolution_be_NS->SetMarkerColor(kBlue);
	resolution_be_NS->SetLineColor(kBlue);
	resolution_be_NS->SetStats(0);
	resolution_be_S->SetStats(0);
// 	resolution_be_S->GetYaxis()->SetLabelSize(0.1);
// 	resolution_be_S->GetYaxis()->SetTitleSize(15);
 	TLegend *legend_be_3 = new TLegend(0.25,0.8,0.5,0.9);
 	legend_be_3->AddEntry(resolution_be_NS, "No Scale", "lep");
 	legend_be_3->AddEntry(resolution_be_S,"Scale", "lep");
 	legend_be_3->Draw();

	be->cd();
	pad22->SetTopMargin(0);
	pad22->SetBottomMargin(0.2);
	pad22->SetGridx();
	pad22->Draw(); 	
 	pad22->cd();
	gr_be = new TGraphErrors(15, MASS_BIN, BE, MASS_BINS_err, BE_err);
// 	gr_be->SetTitle("Ratio: NoScale / Scale - BE");
	gr_be->SetTitle(" ");
	gr_be->SetMarkerStyle(21);
	gr_be->SetMarkerSize(0.5);
	gr_be->GetYaxis()->SetTitleSize(12);
	gr_be->GetYaxis()->SetTitleFont(43);
	gr_be->GetYaxis()->SetLabelSize(0.1);
	gr_be->GetXaxis()->SetLabelSize(0.07);	
	gr_be->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
	gr_be->GetXaxis()->SetTitleSize(12);
	gr_be->GetXaxis()->SetTitleFont(43);
	gr_be->GetXaxis()->SetRangeUser(0, 6000);

	gr_be->GetYaxis()->SetTitle("(Scale - NoScale) / NoScale");
	gr_be->Draw("AE*");
	
	be->Print("./MassScale/MC/BE_vs_Mass.pdf");
	be->Print("./MassScale/MC/BE_vs_Mass.png");


	TCanvas *prova = new TCanvas("Mean vs mass: BB", "Mean vs mass: BB", 210,45,750,500);
	mean_bb_S->Draw("E");
	mean_bb_NS->Draw("E same");
	mean_bb_S->SetTitle("Mean vs mass: BB");
	mean_bb_NS->SetTitle("Mean vs mass: BB");
	mean_bb_NS->SetMarkerStyle(20);
	mean_bb_S->SetMarkerStyle(20);
	mean_bb_NS->SetMarkerSize(0.75);
	mean_bb_S->SetMarkerSize(0.75);
	mean_bb_S->SetMarkerColor(kRed);
	mean_bb_S->SetLineColor(kRed);
	mean_bb_NS->SetMarkerColor(kBlue);
	mean_bb_NS->SetLineColor(kBlue);
	mean_bb_NS->SetStats(0);
	mean_bb_S->SetStats(0);
 	TLegend *legend = new TLegend(0.65,0.8,0.9,0.9);
 	legend->AddEntry(mean_bb_NS, "No Scale", "lep");
 	legend->AddEntry(mean_bb_S,"Scale", "lep");
 	legend->Draw();
 	
 	prova->Print("./MassScale/MC/Mean_MC_vs_mass_BB.png");
 	prova->Print("./MassScale/MC/Mean_MC_vs_mass_BB.pdf");
 	

	TCanvas *prova_be = new TCanvas("Mean vs mass: BE+EE", "Mean vs mass: BE+EE", 210,45,750,500);
	mean_be_S->Draw("E");
	mean_be_NS->Draw("E same");
	mean_be_S->SetTitle("Mean vs mass: BE+EE");
	mean_be_NS->SetTitle("Mean vs mass: BE+EE");
	mean_be_NS->SetMarkerStyle(20);
	mean_be_S->SetMarkerStyle(20);
	mean_be_NS->SetMarkerSize(0.75);
	mean_be_S->SetMarkerSize(0.75);
	mean_be_S->SetMarkerColor(kRed);
	mean_be_S->SetLineColor(kRed);
	mean_be_NS->SetMarkerColor(kBlue);
	mean_be_NS->SetLineColor(kBlue);
	mean_be_NS->SetStats(0);
	mean_be_S->SetStats(0);
 	TLegend *legend_be = new TLegend(0.65,0.8,0.9,0.9);
 	legend_be->AddEntry(mean_be_NS, "No Scale", "lep");
 	legend_be->AddEntry(mean_be_S,"Scale", "lep");	
 	legend_be->Draw();
 	
 	prova_be->Print("./MassScale/MC/Mean_MC_vs_mass_BE.png");
 	prova_be->Print("./MassScale/MC/Mean_MC_vs_mass_BE.pdf");

      
	TCanvas *prova_2 = new TCanvas("Resolution vs mass: BB", "Resolution vs mass: BB", 210,45,750,500);
	resolution_bb_S->Draw("E");
	resolution_bb_NS->Draw("E same");
	resolution_bb_S->SetTitle("Mass resolution vs mass: BB");
	resolution_bb_NS->SetTitle("Mass resolution vs mass: BB");
	resolution_bb_NS->SetMarkerStyle(20);
	resolution_bb_S->SetMarkerStyle(20);
	resolution_bb_NS->SetMarkerSize(0.75);
	resolution_bb_S->SetMarkerSize(0.75);
	resolution_bb_S->SetMarkerColor(kRed);
	resolution_bb_S->SetLineColor(kRed);
	resolution_bb_NS->SetMarkerColor(kBlue);
	resolution_bb_NS->SetLineColor(kBlue);
	resolution_bb_NS->SetStats(0);
	resolution_bb_S->SetStats(0);
 	TLegend *legend_2 = new TLegend(0.25,0.8,0.5,0.9);
 	legend_2->AddEntry(resolution_bb_NS, "No Scale", "lep");
 	legend_2->AddEntry(resolution_bb_S,"Scale", "lep");
 	legend_2->Draw();
 	
 	prova_2->Print("./MassScale/MC/Resolution_MC_vs_mass_BB.png");
 	prova_2->Print("./MassScale/MC/Resolution_MC_vs_mass_BB.pdf");
 	

	TCanvas *prova_be_2 = new TCanvas("Resolution vs mass: BE+EE", "Resolution vs mass: BE+EE", 210,45,750,500);
	resolution_be_S->Draw("E");
	resolution_be_NS->Draw("E same");
	resolution_be_S->SetTitle("Mass resolution vs mass: BE+EE");
	resolution_be_NS->SetTitle("Mass resolution vs mass: BE+EE");
	resolution_be_NS->SetMarkerStyle(20);
	resolution_be_S->SetMarkerStyle(20);
	resolution_be_NS->SetMarkerSize(0.75);
	resolution_be_S->SetMarkerSize(0.75);
	resolution_be_S->SetMarkerColor(kRed);
	resolution_be_S->SetLineColor(kRed);
	resolution_be_NS->SetMarkerColor(kBlue);
	resolution_be_NS->SetLineColor(kBlue);
	resolution_be_NS->SetStats(0);
	resolution_be_S->SetStats(0);
 	TLegend *legend_be_2 = new TLegend(0.25,0.8,0.5,0.9);
 	legend_be_2->AddEntry(resolution_be_NS, "No Scale", "lep");
 	legend_be_2->AddEntry(resolution_be_S,"Scale", "lep");	
 	legend_be_2->Draw();
 	
 	prova_be_2->Print("./MassScale/MC/Resolution_MC_vs_mass_BE.png");
 	prova_be_2->Print("./MassScale/MC/Resolution_MC_vs_mass_BE.pdf");


	///////////////   W/O wrt W
// 	TCanvas *c_3 = new TCanvas("MC W & W/O scale correction", "MC W & W/O scale correction", 210,45,750,500);
// 	res->Draw();
// 	c_3->Print("./MassScale/MC/MassScale_Resolution_MC.png");
// 
// 	TCanvas *c = new TCanvas("MC W & W/O scale correction BE+EE", "MC W & W/O scale correction BE+EE", 210,45,750,500);
// 	res_be->Draw();
// 	c->Print("./MassScale/MC/MassScale_Resolution_BE_MC.png");
// 
// 	TCanvas *c_2 = new TCanvas("MC W & W/O scale correction BB", "MC W & W/O scale correction BB", 210,45,750,500);
// 	res_bb->Draw();
// 	c_2->Print("./MassScale/MC/MassScale_Resolution_BB_MC.png");
	
	
	////////////////    W/O wrt gen
	TCanvas *c_33 = new TCanvas("MC W scale correction", "MC W scale correction", 210,45,750,500);
	res_S->Draw("COLZ");
	c_33->Print("./MassScale/MC/MassScale_Resolution_W_MC.png");

	TCanvas *cc = new TCanvas("MC W correction BE+EE", "MC W correction BE+EE", 210,45,750,500);
	res_be_S->Draw("COLZ");
	cc->Print("./MassScale/MC/MassScale_Resolution_BE_W_MC.png");

	TCanvas *c_22 = new TCanvas("MC W correction BB", "MC W correction BB", 210,45,750,500);
	res_bb_S->Draw("COLZ");
	res_bb_S->SetMaximum(500);
	c_22->Print("./MassScale/MC/MassScale_Resolution_BB_W_MC.png");
	
	
	///////////////    W wrt gen
	TCanvas *c_333 = new TCanvas("MC W/O correction", "MC W/O correction", 210,45,750,500);
	res_NS->Draw("COLZ");
	c_333->Print("./MassScale/MC/MassScale_Resolution_WO_MC.png");

	TCanvas *ccc = new TCanvas("MC W/O correction BE+EE", "MC W/O correction BE+EE", 210,45,750,500);
	res_be_NS->Draw("COLZ");
	ccc->Print("./MassScale/MC/MassScale_Resolution_BE_WO_MC.png");

	TCanvas *c_222 = new TCanvas("MC W/O correction BB", "MC W/O correction BB", 210,45,750,500);
	res_bb_NS->Draw("COLZ");
	c_222->Print("./MassScale/MC/MassScale_Resolution_BB_WO_MC.png");


















// 	//////////////////////  plus W/O wrt gen
// 	TCanvas *c_33_plus = new TCanvas("MC W pt plus scale correction", "MC W pt plus scale correction", 210,45,750,500);
// 	pt_plus_res_S->Draw();
// 	c_33_plus->Print("./MassScale/MC/PtPlus_Resolution_W_MC.png");
// 
// 	TCanvas *cc_plus = new TCanvas("MC W pt plus correction BE+EE", "MC W pt plus correction BE+EE", 210,45,750,500);
// 	pt_plus_res_be_S->Draw();
// 	cc_plus->Print("./MassScale/MC/PtPlus_Resolution_BE_W_MC.png");
// 
// 	TCanvas *c_22_plus = new TCanvas("MC W pt plus correction BB", "MC W pt plus correction BB", 210,45,750,500);
// 	pt_plus_res_bb_S->Draw();
// 	c_22_plus->Print("./MassScale/MC/PtPlus_Resolution_BB_W_MC.png");
// 	
// 	
// 	////////////////////// 	  plus W wrt gen
// 	TCanvas *c_333_plus = new TCanvas("MC W/O pt plus correction", "MC W/O pt plus correction", 210,45,750,500);
// 	pt_plus_res_NS->Draw();
// 	c_333_plus->Print("./MassScale/MC/PtPlus_Resolution_WO_MC.png");
// 
// 	TCanvas *ccc_plus = new TCanvas("MC W/O pt plus correction BE+EE", "MC W/O pt plus correction BE+EE", 210,45,750,500);
// 	pt_plus_res_be_NS->Draw();
// 	ccc->Print("./MassScale/MC/PtPlus_Resolution_BE_WO_MC.png");
// 
// 	TCanvas *c_222_plus = new TCanvas("MC W/O pt plus correction BB", "MC W/O pt plus correction BB", 210,45,750,500);
// 	pt_plus_res_bb_NS->Draw();
// 	c_222_plus->Print("./MassScale/MC/PtPlus_Resolution_BB_WO_MC.png");
// 	
// 	
// 	
// 	////////////////////// 	  minus W/O wrt gen
// 	TCanvas *c_33_minus = new TCanvas("MC W pt minus scale correction", "MC W pt minus scale correction", 210,45,750,500);
// 	pt_minus_res_S->Draw();
// 	c_33_minus->Print("./MassScale/MC/PtPlus_Resolution_W_MC.png");
// 
// 	TCanvas *cc_minus = new TCanvas("MC W pt minus correction BE+EE", "MC W pt minus correction BE+EE", 210,45,750,500);
// 	pt_minus_res_be_S->Draw();
// 	cc_minus->Print("./MassScale/MC/PtPlus_Resolution_BE_W_MC.png");
// 
// 	TCanvas *c_22_minus = new TCanvas("MC W pt minus correction BB", "MC W pt minus correction BB", 210,45,750,500);
// 	pt_minus_res_bb_S->Draw();
// 	c_22_minus->Print("./MassScale/MC/PtPlus_Resolution_BB_W_MC.png");
// 	
// 	
// 	////////////////////// 	  minus W wrt gen
// 	TCanvas *c_333_minus = new TCanvas("MC W/O pt minus correction", "MC W/O pt minus correction", 210,45,750,500);
// 	pt_minus_res_NS->Draw();
// 	c_333_minus->Print("./MassScale/MC/PtPlus_Resolution_WO_MC.png");
// 
// 	TCanvas *ccc_minus = new TCanvas("MC W/O pt minus correction BE+EE", "MC W/O pt minus correction BE+EE", 210,45,750,500);
// 	pt_minus_res_be_NS->Draw();
// 	ccc->Print("./MassScale/MC/PtPlus_Resolution_BE_WO_MC.png");
// 
// 	TCanvas *c_222_minus = new TCanvas("MC W/O pt minus correction BB", "MC W/O pt minus correction BB", 210,45,750,500);
// 	pt_minus_res_bb_NS->Draw();
// 	c_222_minus->Print("./MassScale/MC/PtPlus_Resolution_BB_WO_MC.png");






// 	c->SetLogy();
// 
// 	c->Modified();
// 	c->Update();
// 		
// 	h_scale->Draw();
// 	gPad->Update();
// 	TPaveStats *tps_scale = (TPaveStats*)c->GetPrimitive("stats");
// 	tps_scale->SetName("Scale correction");
// 	tps_scale->SetTextColor(kBlue);
// 	tps_scale->SetX1NDC(0.7);
// 	tps_scale->SetY1NDC(0.7);
// 	tps_scale->SetY2NDC(0.8);
// 
// 	
// 	h_Nscale->Draw();
// 	gPad->Update();
// 	TPaveStats *tps_Nscale = (TPaveStats*)c->GetPrimitive("stats");
// 	tps_Nscale->SetName("NO Scale correction");
// 	tps_Nscale->SetTextColor(kRed);
// 	tps_Nscale->SetX1NDC(0.7);
// 	tps_Nscale->SetY1NDC(0.8);
// 	tps_Nscale->SetY2NDC(0.9);
// 
// 	
// 	h_Nscale->SetLineColor(kBlue);
// 	h_scale->SetLineColor(kRed);
// 
// 	h_scale->Rebin(rebin);
// 	h_Nscale->Rebin(rebin);
// 	
// 	h_scale->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	h_Nscale->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	h_scale->GetYaxis()->SetTitle("Event / 20 GeV");
// 	h_Nscale->GetYaxis()->SetTitle("Event / 20 GeV");
// 	h_scale->SetTitle("Comparison MC W - W/O scale correction");
// 	h_Nscale->SetTitle("Comparison MC W - W/O scale correction");
// 	
// 	h_scale->SetMinimum(0.01);
// 	h_scale->SetMaximum(5*10e5);
// 	h_Nscale->SetMinimum(0.01);
// 	h_Nscale->SetMaximum(5*10e5);
// 	
// 	h_Nscale->Draw();
// 	h_scale->Draw("same");
// // 	tps_scale->Draw("same");
// // 	tps_Nscale->Draw("same");
// 
// 	TLegend *l = new TLegend(0.6,0.75,0.9,0.9);
// 	l->AddEntry(h_scale,"MC W scale correction", "l");	
// 	l->AddEntry(h_Nscale,"MC W/O scale correction", "l");
// 	l->Draw();
// 		
// 	c->Print("Mass_scale_Resolution_BB.png");		
}
