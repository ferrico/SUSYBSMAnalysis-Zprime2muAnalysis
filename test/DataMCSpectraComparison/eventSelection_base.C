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

  UInt_t event;
  UInt_t run;
  Int_t prev_event = -88;
  UInt_t lumin;

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
  float lep_qOverPt[2];
  float vertex_chi2;
  
  float met_pt;
  
 Color_t color;
 float gM;
 float kFactor;
 float kFactor_BB;
 float kFactor_BE;
 
 TLegend * legend = new TLegend(0.85,0.65,0.99,0.95);

	TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

	
    TString samples[45] =  {
     						"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200",
                            "Wjets",
    						"dyInclusive50",
                            "Wantitop", "tW", 
                            "ttbar_lep50to500", "ttbar_lep_500to800", "ttbar_lep_800to1200", "ttbar_lep_1200to1800", "ttbar_lep1800toInf", 							
//                                "ttbar_50to500", "ttbar_500to800", "ttbar_800to1200", "ttbar_1200to1800", "ttbar_1800toInf",                                                      
                         "WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500",
                            "ZZ", "ZZ_ext", 
                            "WZ", "WZ_ext",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
	 						"dyPt0to50", "dyPt50to100", "dyPt100to250", "dyPt250to400", "dyPt400to650", "dyPt650",
                            };
	
	float events[45] = {
						6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
						29705748,
						19385554,  
						6933094, 6952830,
						79092400, 200000, 199800, 200000, 40829, 
						1999000, 200000, 200000, 200000, 38969, 
						990064, 998034, 
						1000000, 2995828,
						2977600, 100000, 100000, 98400, 100000, 95106, 100000, 100000, 100000,
						22782948, (float)39612900, (float)26998200, (float)7190820, 167272, 177101
// 						878212, 1662453, 2833172, 1571199, 48731, 177101					
						};
						
	float sigma[45] = {
						2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
						61526.7,
						6025.2, 
						35.6, 35.6,
						87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
						12.178, 1.385, 0.0566, 0.0035, 0.00005,
					    8.2615, 8.2615, 
					    23.565, 23.565, //16.523, 47.13, 16.523, 47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7,
						5352.58, 363.81, 84.015, 3.2283, 0.43604, 0.04098,
						};
           
	float LUMINOSITY_2016 = 36235.493;
	
	float Z_peak = 0.9638;
	float Z_peak_BB = 0.9688;
	float Z_peak_BE = 0.9610;

	float weight[45] = {0};
	float weight_BB[45] = {0};
	float weight_BE[45] = {0};

	const int    NMBINS = 50;
	const double MMIN = 60., MMAX =3500.;
	double logMbins[NMBINS+1];

    Double_t PT_BINS[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 
    					1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800};
    Int_t  binnum_pt = sizeof(PT_BINS)/sizeof(Double_t)-1;    
    
    Double_t ETA_BINS[] = {-2.4, -1.5, -1.2, -0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9, 1.2, 1.5, 2.4};
    Int_t  binnum_eta = sizeof(ETA_BINS)/sizeof(Double_t)-1;
    
    Double_t PHI_BINS[] = {-3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2};
    Int_t  binnum_phi = sizeof(PHI_BINS)/sizeof(Double_t)-1;
    
    Double_t MET_BINS[] = {0, 25, 50, 75, 100, 200, 350, 500};//, 750, 1000};
    Int_t  binnum_met = sizeof(MET_BINS)/sizeof(Double_t)-1;
    
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 

void DataMC_comparison_v3(){

	gStyle->SetOptFit(1);   
	gStyle->SetOptStat(0);        
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

	for(int i = 0; i < 45; i++){
	
		weight[i] = LUMINOSITY_2016 * sigma[i] / events[i];
		weight[i] *= Z_peak;
// 		weight[i] = 1;
// 		if(i>13) std::cout<<weight[i]<<std::endl;
		weight_BB[i] = weight[i] * Z_peak_BB / Z_peak;
		weight_BE[i] = weight[i] * Z_peak_BE/ Z_peak;
	}
	
	for(int j = 13; j < 45; j++){
		
      	if(j >= 39) continue; 
//       	if(j < 30) continue;
//       	if(j >= 30) continue;
//       	if(j >= 30 and j <= 38) continue;
//       	if(j == 30) continue; 
//       	if(j == 44) continue; 
//       	if(j >= 21) continue;
//       	if(j < 21) continue;

        std::cout<<"opening.. "<<samples[j]<<" --- "<<events[j]<<std::endl;

        TChain *treeMC = new TChain("SimpleNtupler/t");

        treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2016/MC_OK/ana_datamc_" + samples[j] + ".root");        

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
    	 treeMC->SetBranchAddress("lep_qOverPt",lep_qOverPt);
    	 

     //treeMC->SetBranchAddress("gen_lep_pt",gen_lep_pt); 

	     treeMC->SetBranchAddress("met_pt",&met_pt); 	
	     
         printf("opening.. %s %i  --- %.5f ---- %lld\n",samples[j].Data(),j , weight[j], nentries);

//     	for(int p=0; p<nentries; p++){
//      	 for(int p=0; p<20000; p++){
     	 for(int p=0; p<1; p++){

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

				weight[j] *= genWeight;
				weight_BB[j] *= genWeight;
				weight_BE[j] *= genWeight;
				
				gM = gen_dil_mass - 400;
				kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
	 		   	
				if(j >= 30 && j <=38){
					weight[j] *= kFactor;
					weight_BB[j] *= kFactor_BB;
					weight_BE[j] *= kFactor_BE;
				}
   					








				if(j >= 30 && j <=38){
					weight[j] /= kFactor;
					weight_BB[j] /= kFactor_BB;
					weight_BE[j] /= kFactor_BE;
				}

				weight[j] /= genWeight;
				weight_BB[j] /= genWeight;
				weight_BE[j] /= genWeight;
				
	       }//end of condition on event

		 }// end loop p
	 
	 	 if(j == 44) color = kMagenta-10;
	 	 if(j == 43) color = kMagenta-8;
	 	 if(j == 42) color = kMagenta-6;
	 	 if(j == 41) color = kMagenta-2;
	 	 if(j == 40) color = kMagenta;
	 	 if(j == 39) color = kMagenta+3;
	 	 if(j >= 30 && j < 39){ 
	 	 	color = 3;
	 	 	if(j == 30) legend->AddEntry(dimuon_MC_clear, "DY #rightarrow #mu#mu", "f");
	 	 }
// 	 	 if(j == 30) color = kGreen+3;
	 	 if(j > 27 && j < 30){
	 	 	color = kRed+3;
	 	 	if(j == 28) legend->AddEntry(dimuon_MC_clear, "WZ", "f");
	 	 }
	 	 if(j > 25 && j < 28){
	 	 	color = 2;
	 	 	if(j == 26) legend->AddEntry(dimuon_MC_clear, "ZZ", "f");
	 	 }
	 	 if(j > 20 && j < 26){
	 	 	color = kRed-10;
	 	 	if(j == 23) legend->AddEntry(dimuon_MC_clear, "WW", "f");
	 	 }
	 	 if(j > 15 && j < 21){
	 	 	color = 4;
	 	 	if(j == 18) legend->AddEntry(dimuon_MC_clear, "t#bar{t}", "f");
	 	 }
	 	 if(j == 14 || j == 15){
	 	 	color = kBlue+2;
	 	 	if(j == 14) legend->AddEntry(dimuon_MC_clear, "SingleTop", "f");
	 	 }
	 	 if(j == 13){
	 	 	color = kYellow;
	 	 	legend->AddEntry(dimuon_MC_clear, "DY #rightarrow #tau#tau", "f");
	 	 }

	 	 if(j > 30 && j < 39) color = 3;
	 	 if(j == 30) color = kGreen+3;
	 	 if(j > 27 && j < 30) color = kRed+3;
	 	 if(j > 25 && j < 28) color = 2;
	 	 if(j > 20 && j < 26) color = kRed-10;
	 	 if(j > 15 && j < 21) color = 4;
	 	 if(j == 14 || j == 15) color = kBlue+2;		
	 	 
    }// end loop on MC 2016
       
// 	return;
			 
// 	gROOT->Reset();
// 	gROOT->SetBatch();
	
	
}// end of function 
