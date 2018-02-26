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

// gROOT->Reset();
// gROOT->SetBatch();

	    Double_t PT_BINS_BB[] = {72, 100, 152, 200, 300, 452, 800};
	    Int_t  binnum_BB = sizeof(PT_BINS_BB)/sizeof(Double_t)-1;

	    Double_t PT_BINS_BE[] = {52, 72, 100, 152, 200, 300, 452};
	    Int_t  binnum_BE = sizeof(PT_BINS_BB)/sizeof(Double_t)-1;
	    
	    float BB_sigma[6] = {0};
	    float BB_sigma_err[6] = {0};
	    float BE_sigma[6] = {0};
	    float BE_sigma_err[6] = {0};
	    
          std::vector<float> SIGMA;
          std::vector<float> SIGMA_ERR;
		
     TString samples[43] =  {"dyInclusive50", 
	 						"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200", 
	 						"Wantitop", "tW", 
	 						"ttbar_lep50to500", "ttbar_lep_500to800", "ttbar_lep_800to1200", "ttbar_lep_1200to1800", "ttbar_lep1800toInf", 
	 						"Wjets", 
	 						"WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500", 
	 						"ZZ", "WZ", "ZZ_ext", "WZ_ext", 
	 						"dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000", //};
	 						"dyPt100to250", "dyPt250to400", "dyPt400to650", "dyPt650"};


	float events[43] = {19385554,  6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
						6933094, 6952830,
						79092400, 200000, 199800, 200000, 40829, 
						29705748,
						1999000, 200000, 200000, 200000, 38969, 
						990064, 1000000, 998034, 2995828,
						2990000, 100000, 100000, 98000, 96613, 100000, 100000, 100000, 100000,
						2699820, 4854240, 167272, 177101
						};
						
	float sigma[43] = {6025.2, 2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
						35.6, 35.6,
						87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
						61526.7,
						12.178, 1.385, 0.0566, 0.0035, 0.00005,
						16.523, 47.13, 16.523, 47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7,
						83.12, 3.047, 0.3921, 0.03636
						};

	float LUMINOSITY = 39484;

	float weight[43] = {0};
	
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
  float dil_lep_pt[2];
  float lep_phi[2];
  int lep_id[2];
  float lep_pt_err[2];
  float lep_eta[2];
  float dil_lep_eta[2];
  float lep_tk_pt[2];
  float lep_tk_eta[2];
  float lep_tk_phi[2];
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
	
	
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 
	
  Int_t prev_event = -88;
	
			
	for(int i = 0; i < 43; i++){
// 		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight[i] = 1;
	}
	
	double mass;
	double rdil;
	TH2F* DileptonMass_2d_vsPt_MC = new TH2F("DileptonMass_2d_vsPt - MC", "DileptonMass_2d_vsPt - MC", 80, 75, 105., binnum_BB, PT_BINS_BB);
	TH2F* DileptonMass_2d_vsPt_BB_MC = new TH2F("DileptonMass_2d_vsPt: BB - MC", "DileptonMass_2d_vsPt: BB - MC", 80, 75, 105., binnum_BB, PT_BINS_BB);
	TH2F* DileptonMass_2d_vsPt_BE_MC = new TH2F("DileptonMass_2d_vsPt: BE+EE - MC", "DileptonMass_2d_vsPt: BE+EE - MC", 80, 75, 105., binnum_BE, PT_BINS_BE);
	TH2F* DileptonMass_2d_vsPt_DATA = new TH2F("DileptonMass_2d_vsPt - DATA", "DileptonMass_2d_vsPt - DATA", 80, 75, 105., binnum_BB, PT_BINS_BB);
	TH2F* DileptonMass_2d_vsPt_BB_DATA = new TH2F("DileptonMass_2d_vsPt: BB - DATA", "DileptonMass_2d_vsPt: BB - DATA", 80, 75, 105., binnum_BB, PT_BINS_BB);
	TH2F* DileptonMass_2d_vsPt_BE_DATA = new TH2F("DileptonMass_2d_vsPt: BE+EE - DATA", "DileptonMass_2d_vsPt: BE+EE - DATA", 80, 75, 105., binnum_BE, PT_BINS_BE);
	
	TH1F* Mass_dist = new TH1F("Mass_dist", "Mass_dist", 300, 70, 6070);
	TH1F* Mass_dist_BB = new TH1F("Mass_dist_BB", "Mass_dist_BB", 300, 70, 6070);
	TH1F* Mass_dist_BE = new TH1F("Mass_dist_BE", "Mass_dist_BE", 300, 70, 6070);

// 	TH2F* pt_vs_eta_MC[9];
// 	for(int i = 0; i < 9; i++)
// 		pt_vs_eta_MC[i] = new TH2F(samples[30+i], samples[30+i], 80, 75, 105., binnum_BB, PT_BINS_BB);


// 	TH2F* DileptonMass_2d_vsPt_EE = new TH2F("DileptonMass_2d_vsPt: EE", "DileptonMass_2d_vsPt: EE", 80, 75, 105., binnum_BE, PT_BINS_BE);

	TH1F *MC_res = new TH1F("Resolution - MC", "Resolution - MC", binnum_BB, PT_BINS_BB);
	MC_res->GetXaxis()->SetTitle("p_{T} [GeV]");
	MC_res->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	MC_res->SetTitle("Mass residuals");
	TH1F *MC_res_bb = new TH1F("Resolution BB - MC", "Resolution BB - MC", binnum_BB, PT_BINS_BB);
	MC_res_bb->GetXaxis()->SetTitle("p_{T} [GeV]");
	MC_res_bb->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	MC_res_bb->SetTitle("BB mass residuals");
	TH1F *MC_res_be = new TH1F("Resolution BE + EE - MC", "Resolution BE + EE - MC", binnum_BE, PT_BINS_BE);
	MC_res_be->GetXaxis()->SetTitle("p_{T} [GeV]");
	MC_res_be->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	MC_res_be->SetTitle("BE + EE mass residuals");

	TH1F *DATA_res = new TH1F("Resolution - DATA", "Resolution - DATA", binnum_BB, PT_BINS_BB);
	DATA_res->GetXaxis()->SetTitle("p_{T} [GeV]");
	DATA_res->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	DATA_res->SetTitle("Mass residuals");
	TH1F *DATA_res_bb = new TH1F("Resolution BB - DATA", "Resolution BB - DATA", binnum_BB, PT_BINS_BB);
	DATA_res_bb->GetXaxis()->SetTitle("p_{T} [GeV]");
	DATA_res_bb->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	DATA_res_bb->SetTitle("BB mass residuals");
	TH1F *DATA_res_be = new TH1F("Resolution BE + EE - DATA", "Resolution BE + EE - DATA", binnum_BE, PT_BINS_BE);
	DATA_res_be->GetXaxis()->SetTitle("p_{T} [GeV]");
	DATA_res_be->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	DATA_res_be->SetTitle("BE + EE mass residuals");

// 	TH1F *res_ee = new TH1F("Resolution EE", "Resolution EE", binnum_BE, PT_BINS_BE);
// 	res_ee->GetXaxis()->SetTitle("p_{T} [GeV]");
// 	res_ee->GetYaxis()->SetTitle("#sigma (Z peak)"); 
// 	res_ee->SetTitle("EE mass residuals");
	
	
                       	
  for(int j=43; j < 43; j++){   

	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    

     TChain *treeMC = new TChain("SimpleNtupler/t");
//      treeMC->Add("../DataMCSpectraComparison/mc/ana_datamc_"+samples[j]+".root");
     treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/ana_datamc_"+samples[j]+".root");
     treeMC->SetBranchAddress("event",&event);
     treeMC->SetBranchAddress("run",&run);
     treeMC->SetBranchAddress("lumi",&lumi);
     treeMC->SetBranchAddress("dil_mass",&dil_mass);
     treeMC->SetBranchAddress("cos_angle",&cos_angle);
     treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
     treeMC->SetBranchAddress("lep_pt",lep_pt);
     treeMC->SetBranchAddress("dil_lep_pt",dil_lep_pt);
     treeMC->SetBranchAddress("lep_id",lep_id);
     treeMC->SetBranchAddress("lep_eta",lep_eta);
//      treeMC->SetBranchAddress("dil_lep_eta",dil_lep_eta);
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
     treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeMC->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeMC->SetBranchAddress("vertex_m",&vertex_m);
     treeMC->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);
     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);

	ne = treeMC->GetEntries();
	std::cout<<"START"<<std::endl;
	for ( int p=0; p < ne ;p++){
// 	for ( int p=0; p<100000 ;p++){
		if(p % 100000 == 0) std::cout<<p<<" su "<<ne<<std::endl;		
		
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

			if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
			prev_event = event;

		
			mass         = dil_mass; //vertex_m			
			DileptonMass_2d_vsPt_MC->Fill(mass, dil_lep_pt[0], weight[j]);
		    DileptonMass_2d_vsPt_MC->Fill(mass, dil_lep_pt[1], weight[j]);
		    Mass_dist->Fill(mass, weight[j]);
  
		  if (abs(lep_eta[0]) < 1.2 && abs(lep_eta[1]) < 1.2) { 
	    	DileptonMass_2d_vsPt_BB_MC->Fill(mass, dil_lep_pt[0], weight[j]);
	    	DileptonMass_2d_vsPt_BB_MC->Fill(mass, dil_lep_pt[1], weight[j]);
		    Mass_dist_BB->Fill(mass, weight[j]);
		  }
		  else { 
    		DileptonMass_2d_vsPt_BE_MC->Fill(mass, dil_lep_pt[0], weight[j]);
		    DileptonMass_2d_vsPt_BE_MC->Fill(mass, dil_lep_pt[1], weight[j]);
		    Mass_dist_BE->Fill(mass, weight[j]);
		  }	
		  
// 		  if((dil_lep_eta[0] > 1.2 || dil_lep_eta[0] < -1.2) && (dil_lep_eta[1] > 1.2 || dil_lep_eta[1] < -1.2)) {
// 		    DileptonMass_2d_vsPt_EE->Fill(mass, dil_lep_pt[0]);
// 		    DileptonMass_2d_vsPt_EE->Fill(mass, dil_lep_pt[1]);
// 		  }                                                                        

// 		  pt_vs_eta_MC[j-30]->Fill(mass, dil_lep_pt[0]);
// 		  pt_vs_eta_MC[j-30]->Fill(mass, dil_lep_pt[1]);
					
		} //if selection
	} // for entries
 } // for samples
 
 int counting_data[4] = {0};

//////////////////////////
////////// DATA //////////
//////////////////////////

     TChain *treeDATA = new TChain("SimpleNtupler/t");
     treeDATA->Add("../DataMCSpectraComparison/data_OK/Run2017MuonsOnly/ana_datamc_data.root");
     treeDATA->SetBranchAddress("event",&event);
     treeDATA->SetBranchAddress("run",&run);
     treeDATA->SetBranchAddress("lumi",&lumi);
     treeDATA->SetBranchAddress("dil_mass",&dil_mass);
     treeDATA->SetBranchAddress("cos_angle",&cos_angle);
     treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeDATA->SetBranchAddress("dil_chosen",&dil_chosen);
     treeDATA->SetBranchAddress("lep_pt",lep_pt);
     treeDATA->SetBranchAddress("dil_lep_pt",dil_lep_pt);
     treeDATA->SetBranchAddress("lep_id",lep_id);
     treeDATA->SetBranchAddress("lep_eta",lep_eta);
//      treeDATA->SetBranchAddress("dil_lep_eta",dil_lep_eta);
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
	for ( int p=0; p < ne ;p++){
// 	for ( int p=0; p<100000 ;p++){
		if(p % 100000 == 0) std::cout<<p<<" su "<<ne<<std::endl;		
		
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

			
// 			if(mass > 900)
// 				counting_data[3]++;
// 			else{
// 				if(mass > 600)
// 					counting_data[2]++;
// 				else{
// 					if(mass > 400) 
// 						counting_data[1]++;
// 					else
// 						counting_data[0]++;
// 				}
// 			}
		
// 			lep_1.SetPtEtaPhiM(lep_tk_pt[0], lep_tk_eta[0], lep_tk_phi[0], 0.105);
// 			lep_2.SetPtEtaPhiM(lep_tk_pt[1], lep_tk_eta[1], lep_tk_phi[1], 0.105);
// 			ZPrime = lep_1 + lep_2;
// 
// 			mass 		 = ZPrime.M(); // resolution using Tracker track		
// 			DileptonMass_2d_vsPt_DATA->Fill(mass, lep_tk_pt[0]);
// 		    DileptonMass_2d_vsPt_DATA->Fill(mass, lep_tk_pt[1]);
//   
// 		  if (abs(lep_eta[0]) < 1.2 && abs(lep_eta[1]) < 1.2) { 
// 	    	DileptonMass_2d_vsPt_BB_DATA->Fill(mass, lep_tk_pt[0]);
// 	    	DileptonMass_2d_vsPt_BB_DATA->Fill(mass, lep_tk_pt[1]);
// 		  }
// 		  else { 
//     		DileptonMass_2d_vsPt_BE_DATA->Fill(mass, lep_tk_pt[0]);
// 		    DileptonMass_2d_vsPt_BE_DATA->Fill(mass, lep_tk_pt[1]);
// 		  }	
    
			mass         = dil_mass; //vertex_m	
		
			DileptonMass_2d_vsPt_DATA->Fill(mass, dil_lep_pt[0]);
		    DileptonMass_2d_vsPt_DATA->Fill(mass, dil_lep_pt[1]);
  
		  if (abs(lep_eta[0]) < 1.2 && abs(lep_eta[1]) < 1.2) { 
	    	DileptonMass_2d_vsPt_BB_DATA->Fill(mass, dil_lep_pt[0]);
	    	DileptonMass_2d_vsPt_BB_DATA->Fill(mass, dil_lep_pt[1]);
		  }
		  else { 
    		DileptonMass_2d_vsPt_BE_DATA->Fill(mass, dil_lep_pt[0]);
		    DileptonMass_2d_vsPt_BE_DATA->Fill(mass, dil_lep_pt[1]);
		  }	
// 		  if((dil_lep_eta[0] > 1.2 || dil_lep_eta[0] < -1.2) && (dil_lep_eta[1] > 1.2 || dil_lep_eta[1] < -1.2)) {
// 		    DileptonMass_2d_vsPt_EE->Fill(mass, dil_lep_pt[0]);
// 		    DileptonMass_2d_vsPt_EE->Fill(mass, dil_lep_pt[1]);
// 		  }                                                                         

			
			
			
			
		} //if selection
	} // for entries

//////////////////////////
////////// DATA //////////
//////////////////////////


bool save = true;
	
// TH2F* tmp;
// TH2F* prova = new TH2F("prova", "prova", 80, 75, 105., binnum_BB, PT_BINS_BB);
// for(int i = 0; i < 9; i++){	
// 		 TH2F* tmp = (TH2F*)pt_vs_eta_MC[i]->Clone();
// 	   	 if(i==0){
// 	   	 	tmp->Scale(weight[30+i]);
// 	   	 	prova = tmp;
// 	   	 }
// 	   	 else{
// 	   	 	tmp->Scale(weight[30+i]);
// 	   	 	prova->Add(tmp);
// 	   	 }
// 
//  	NOME = Form("%s - MC", samples[30+i].Data());
//  	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
//  	if(i==0)
//  		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./ciao.pdf[");
// 	canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./ciao.pdf");
//  	if(i==8)
//  		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./ciao.pdf]");
//     canvas->Write();	   	 	   	
// }
//    TCanvas* c44 = new TCanvas("prova", "prova", 200,10,700,500);
//    prova->Draw("COLZ");



/*
 for(int i = 1; i <= DileptonMass_2d_vsPt_BB_MC->GetNbinsY(); i++){
//  for(int i = 1; i <= DileptonMass_2d_vsPt_BB_MC->GetNbinsY(); i++){
 	
 	NOME = Form("DimuonMass Resolution at Z BB: %.0f < m < %.0f - MC", PT_BINS_BB[i-1], PT_BINS_BB[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_BB_MC->ProjectionX(NOME,i,i);
//  	TH1D* proiezione = DileptonMass_2d_vsPt_BB_MC->ProjectionX(NOME,i,i);
 	proiezione->Rebin(2);
 	proiezione->Draw();
 	
 	fit_min = 75;//proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = 105;//proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 1.4,1.3);
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BB: Mass for %.0f < p_{T} <%.0f", PT_BINS_BB[i-1], PT_BINS_BB[i]);
    proiezione->SetTitle(title);
    proiezione->GetXaxis()->SetRangeUser(fit_min, fit_max);
    proiezione->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
        
    MC_res_bb->SetBinContent(i, funct->GetParameter(2));
    MC_res_bb->SetBinError(i, funct->GetParError(2));
    
    BB_sigma[i-1] = funct->GetParameter(2);
    BB_sigma_err[i-1] = funct->GetParError(2);
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(80, yPos, longstring);
    }
 	 
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BB_resolutionAtZ_MC.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BB_resolutionAtZ_MC.pdf");
 	if(i==DileptonMass_2d_vsPt_BB_MC->GetNbinsY())
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BB_resolutionAtZ_MC.pdf]");
    canvas->Write();

}

 for(int i = 1; i <= DileptonMass_2d_vsPt_BE_MC->GetNbinsY(); i++){
 	
 	NOME = Form("DimuonMass Resolution at Z BE+EE: %.0f < m < %.0f - MC", PT_BINS_BE[i-1], PT_BINS_BE[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_BE_MC->ProjectionX(NOME,i,i);
 	proiezione->Draw();

 	fit_min = 75;//proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = 105;//proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BE+EE: Mass for %.0f < p_{T} <%.0f", PT_BINS_BE[i-1], PT_BINS_BE[i]);
    proiezione->SetTitle(title); 
    proiezione->GetXaxis()->SetRangeUser(fit_min, fit_max);
    proiezione->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
       
    MC_res_be->SetBinContent(i, funct->GetParameter(2));
    MC_res_be->SetBinError(i, funct->GetParError(2));
    
    BE_sigma[i-1] = funct->GetParameter(2);
    BE_sigma_err[i-1] = funct->GetParError(2);

    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(80, yPos, longstring);
    }

 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BE_resolutionAtZ_MC.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BE_resolutionAtZ_MC.pdf");
 	if(i==DileptonMass_2d_vsPt_BE_MC->GetNbinsY())
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BE_resolutionAtZ_MC.pdf]");
    canvas->Write();

}
*/

 for(int i = 1; i <= DileptonMass_2d_vsPt_BB_DATA->GetNbinsY(); i++){
 	
 	NOME = Form("DimuonMass Resolution at Z BB: %.0f < m < %.0f - DATA", PT_BINS_BB[i-1], PT_BINS_BB[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_BB_DATA->ProjectionX(NOME,i,i);
//  	proiezione->Draw();
 	
 	fit_min = 75;//proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = 105;//proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 1.4,1.3);
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BB: Mass for %.0f < p_{T} <%.0f", PT_BINS_BB[i-1], PT_BINS_BB[i]);
    proiezione->SetTitle(title);
    proiezione->GetXaxis()->SetRangeUser(fit_min, fit_max);
    proiezione->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
        
    DATA_res_bb->SetBinContent(i, funct->GetParameter(2));
    DATA_res_bb->SetBinError(i, funct->GetParError(2));
    
    BB_sigma[i-1] = funct->GetParameter(2);
    BB_sigma_err[i-1] = funct->GetParError(2);
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(80, yPos, longstring);
    }

 if(save){ 	 	
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BB_resolutionAtZ_DATA.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BB_resolutionAtZ_DATA.pdf");
 	if(i==DileptonMass_2d_vsPt_BB_DATA->GetNbinsY())
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BB_resolutionAtZ_DATA.pdf]");
    canvas->Write();
}

}

 for(int i = 1; i <= DileptonMass_2d_vsPt_BE_DATA->GetNbinsY(); i++){
 	
 	NOME = Form("DimuonMass Resolution at Z BE+EE: %.0f < m < %.0f - DATA", PT_BINS_BE[i-1], PT_BINS_BE[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_BE_DATA->ProjectionX(NOME,i,i);
//  	proiezione->Draw();
 	
 	fit_min = 75;//proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = 105;//proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BE+EE: Mass for %.0f < p_{T} <%.0f", PT_BINS_BE[i-1], PT_BINS_BE[i]);
    proiezione->SetTitle(title); 
    proiezione->GetXaxis()->SetRangeUser(fit_min, fit_max);
    proiezione->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
       
    DATA_res_be->SetBinContent(i, funct->GetParameter(2));
    DATA_res_be->SetBinError(i, funct->GetParError(2));
    
    BE_sigma[i-1] = funct->GetParameter(2);
    BE_sigma_err[i-1] = funct->GetParError(2);

    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(80, yPos, longstring);
    }

 if(save){ 	 	
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BE_resolutionAtZ_DATA.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BE_resolutionAtZ_DATA.pdf");
 	if(i==DileptonMass_2d_vsPt_BE_DATA->GetNbinsY())
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./BE_resolutionAtZ_DATA.pdf]");
    canvas->Write();
}

}


//////////////////////////
///// BE INUTILE /////////
//////////////////////////
/*
 for(int i = 1; i <= DileptonMass_2d_vsPt_EE->GetNbinsY(); i++){
 	
 	NOME = Form("DimuonMass Resolution at Z EE: %.0f < m < %.0f", PT_BINS_BB[i-1], PT_BINS_BB[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMass_2d_vsPt_EE->ProjectionX(NOME,i,i);
//  	proiezione->Draw();
 	
 	fit_min = 75;//proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = 105;//proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("EE: Mass for %.0f < p_{T} <%.0f", PT_BINS_BB[i-1], PT_BINS_BB[i]);
    proiezione->SetTitle(title);
    proiezione->GetXaxis()->SetRangeUser(fit_min, fit_max);
    
    res_ee->SetBinContent(i, funct->GetParameter(2));
    res_ee->SetBinError(i, funct->GetParError(2));
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
    for(int i = 0; i < funct->GetNpar(); i++){
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(75, yPos, longstring);
    }
    
 	 	
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./EE_resolutionAtZ.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./EE_resolutionAtZ.pdf");
 	if(i==DileptonMass_2d_vsPt_EE->GetNbinsY())
 		canvas->Print("./Check_Resolution_Code_2017_94X_AtZ/./EE_resolutionAtZ.pdf]");
    canvas->Write();

}
*/
//////////////////////////
///// BE INUTILE /////////
//////////////////////////


// TLegend *legend_MC = new TLegend(0.6,0.75,0.85,0.85);
// legend_MC->AddEntry(DATA_res_bb, "DATA", "lep");
// legend_MC->AddEntry(MC_res_bb,"DY", "lep");
 
 TCanvas* c_BB = new TCanvas("res_BB", "res_BB", 500, 500);
//  c_BB->SetGrid();
 TPad* pad_1_bb = new TPad("pad1_bb", "pad1_bb", 0.05, 0.25, 1, 1);
 pad_1_bb->SetBottomMargin(0);
 pad_1_bb->SetGridx();
 pad_1_bb->Draw();
 pad_1_bb->cd();
 DATA_res_bb->SetTitle("");
 DATA_res_bb->SetLineColor(kBlack);
 DATA_res_bb->SetMarkerStyle(20);
 DATA_res_bb->SetMarkerSize(0.5);
 DATA_res_bb->SetMarkerColor(kBlack);
 DATA_res_bb->GetYaxis()->SetRangeUser(0, 8);
 DATA_res_bb->Draw();
//  MC_res_bb->SetTitle("");
//  MC_res_bb->SetLineColor(kRed+1);
//  MC_res_bb->SetMarkerStyle(20);
//  MC_res_bb->SetMarkerSize(0.5);
//  MC_res_bb->SetMarkerColor(kRed+1);
//  MC_res_bb->GetYaxis()->SetRangeUser(0, 8);
//  MC_res_bb->Draw("same");
//  legend_MC->Draw();

 c_BB->Update();
 c_BB->cd();
/*
 TPad* pad_2_bb = new TPad("pad2_bb", "pad2_bb", 0.05, 0.05, 1, 0.25);
 pad_2_bb->SetTopMargin(0);
 pad_2_bb->SetBottomMargin(0.2);
 pad_2_bb->SetGridx();
 pad_2_bb->Draw();
 pad_2_bb->cd();
 TH1F* ratio_bb = (TH1F*)DATA_res_bb->Clone();
 ratio_bb->Divide(MC_res_bb);
 ratio_bb->SetLineColor(kBlack);
 ratio_bb->SetTitle("");
 ratio_bb->GetYaxis()->SetRangeUser(0, 2);
 ratio_bb->GetXaxis()->SetLabelFont(42);
//  ratio_bb->GetXaxis()->SetLabelOffset(0.007);
 ratio_bb->GetXaxis()->SetLabelSize(0.1);
 ratio_bb->GetXaxis()->SetTitle("p_{T} [GeV]");
 ratio_bb->GetXaxis()->SetTitleSize(1);
 ratio_bb->GetYaxis()->SetLabelSize(0.1);
 ratio_bb->GetYaxis()->SetTitle("DATA / MC");
 ratio_bb->GetYaxis()->SetTitleSize(1);
 TLine* l_bb = new TLine(PT_BINS_BB[0], 1, PT_BINS_BB[sizeof(PT_BINS_BB)/sizeof(Double_t)-1], 1);
 l_bb->SetLineColor(kRed);
 l_bb->SetLineWidth(2);
 ratio_bb->Draw();
 l_bb->Draw();
 c_BB->Update();*/
 if(save){
 c_BB->Print("./Check_Resolution_Code_2017_94X_AtZ/ResolutionAtZ_BB.png");
 c_BB->Print("./Check_Resolution_Code_2017_94X_AtZ/ResolutionAtZ_BB.pdf");
	}
 TCanvas* c_BE = new TCanvas("res_BE", "res_BE", 500, 500);
//  c_BE->SetGrid();
 TPad* pad_1_be = new TPad("pad1_be", "pad1_be", 0.05, 0.25, 1, 1);
 pad_1_be->SetBottomMargin(0);
 pad_1_be->SetGridx();
 pad_1_be->Draw();
 pad_1_be->cd();
 DATA_res_be->SetTitle("");
 DATA_res_be->SetLineColor(kBlack);
 DATA_res_be->SetMarkerStyle(20);
 DATA_res_be->SetMarkerSize(0.5);
 DATA_res_be->SetMarkerColor(kBlack);
 DATA_res_be->GetYaxis()->SetRangeUser(0, 8);
 DATA_res_be->Draw();
//  MC_res_be->SetTitle("");
//  MC_res_be->SetLineColor(kRed+1);
//  MC_res_be->SetMarkerStyle(20);
//  MC_res_be->SetMarkerSize(0.5);
//  MC_res_be->SetMarkerColor(kRed+1);
//  MC_res_be->GetYaxis()->SetRangeUser(0, 8);
//  MC_res_be->Draw("same");
//  legend_MC->Draw();
/* 
 c_BE->Update();
 c_BE->cd();
 
 TPad* pad_2_be = new TPad("pad2_be", "pad2_be", 0.05, 0.05, 1, 0.25);
 pad_2_be->SetTopMargin(0);
 pad_2_be->SetBottomMargin(0.2);
 pad_2_be->SetGridx();
 pad_2_be->Draw();
 pad_2_be->cd();
 TH1F* ratio_be = (TH1F*)DATA_res_be->Clone();
 ratio_be->Divide(MC_res_be);
 ratio_be->SetLineColor(kBlack);
 ratio_be->SetTitle("");
 ratio_be->GetYaxis()->SetRangeUser(0, 2);
 ratio_be->GetXaxis()->SetLabelFont(42);
//  ratio_be->GetXaxis()->SetLabelOffset(0.007);
 ratio_be->GetXaxis()->SetLabelSize(0.1);
 ratio_be->GetXaxis()->SetTitle("p_{T} [GeV]");
 ratio_be->GetXaxis()->SetTitleSize(3);
 ratio_be->GetYaxis()->SetLabelSize(0.1);
 ratio_be->GetYaxis()->SetTitle("DATA / MC");
 ratio_be->GetYaxis()->SetTitleSize(3);
 TLine* l_be = new TLine(PT_BINS_BE[0], 1, PT_BINS_BE[sizeof(PT_BINS_BE)/sizeof(Double_t)-1], 1);
 l_be->SetLineColor(kRed);
 l_be->SetLineWidth(2);
 ratio_be->Draw();
 l_be->Draw();*/
 if(save){
 c_BE->Print("./Check_Resolution_Code_2017_94X_AtZ/ResolutionAtZ_BE.png");
 c_BE->Print("./Check_Resolution_Code_2017_94X_AtZ/ResolutionAtZ_BE.pdf");
}



/*
 TCanvas* c11 = new TCanvas("res", "res", 200,10,700,500);
 DileptonMass_2d_vsPt_MC->Draw("COLZ");
 TCanvas* c2 = new TCanvas("res_BB", "res_BB", 200,10,700,500);
 DileptonMass_2d_vsPt_BB_MC->Draw("COLZ");
 TCanvas* c3 = new TCanvas("res_BE", "res_BE", 200,10,700,500);
 DileptonMass_2d_vsPt_BE_MC->Draw("COLZ");
//  TCanvas* c4 = new TCanvas("res_EE", "res_EE", 200,10,700,500);
//  DileptonMass_2d_vsPt_EE->Draw();
 */
 
//   TCanvas* c4 = new TCanvas("mass", "mass", 200,10,700,500);
//  Mass_dist->Draw();
//   TCanvas* c5 = new TCanvas("mass_BB", "mass_BB", 200,10,700,500);
//  Mass_dist_BB->Draw();
//   TCanvas* c6 = new TCanvas("mass_BE", "mass_BE", 200,10,700,500);
//  Mass_dist_BE->Draw();
 


 for(int i = 0; i < 6; i++)
	std::cout<<BB_sigma[i]<<", ";
 std::cout<<""<<std::endl;
 for(int i = 0; i < 6; i++)
	std::cout<<BB_sigma_err[i]<<", ";	
 std::cout<<""<<std::endl;
 for(int i = 0; i < 6; i++)
	std::cout<<BE_sigma[i]<<", ";
 std::cout<<""<<std::endl;
 for(int i = 0; i < 6; i++)
	std::cout<<BE_sigma_err[i]<<", ";	    
 std::cout<<""<<std::endl;

} // main function



