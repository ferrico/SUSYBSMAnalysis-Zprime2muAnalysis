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



void skim(){


	TH1F *h_zz = new TH1F("DimuonMassVertexConstrained", "DimuonMassVertexConstrained", 20000, 0, 20000);
// 	TH1F *h_zz_ext = new TH1F("DimuonMassVertexConstrained", "DimuonMassVertexConstrained", 20000, 0, 20000); 
// 	TH1F *h_wz = new TH1F("DimuonMassVertexConstrained", "DimuonMassVertexConstrained", 20000, 0, 20000);
	
	TH1F *h_zz_bb = new TH1F("DimuonMassVertexConstrained_bb", "DimuonMassVertexConstrained_bb", 20000, 0, 20000);
// 	TH1F *h_zz_ext_bb = new TH1F("DimuonMassVertexConstrained_bb", "DimuonMassVertexConstrained_bb", 20000, 0, 20000); 
// 	TH1F *h_wz_bb = new TH1F("DimuonMassVertexConstrained_bb", "DimuonMassVertexConstrained_bb", 20000, 0, 20000);
// 	
	TH1F *h_zz_be = new TH1F("DimuonMassVertexConstrained_be", "DimuonMassVertexConstrained_be", 20000, 0, 20000);
// 	TH1F *h_zz_ext_be = new TH1F("DimuonMassVertexConstrained_be", "DimuonMassVertexConstrained_be", 20000, 0, 20000); 
// 	TH1F *h_wz_be = new TH1F("DimuonMassVertexConstrained_be", "DimuonMassVertexConstrained_be", 20000, 0, 20000); 
	
	Long64_t ne;
	Long64_t Nne;
    float   DimuonMassVertexConstrained, DimuonMassVertexConstrained_bb, DimuonMassVertexConstrained_be;
    
    	
	TFile *f_zz = new TFile("/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/mc/mc_YesEtaCut_NoScale/MC_OK/ana_datamc_ZZ_ext.root", "READ");
	
	f_zz->cd("Our2016MuonsPlusMuonsMinusHistos");
	TH1F *hist = (TH1F*)gDirectory->Get("DimuonMassVertexConstrained");
	TH1F *hist_bb = (TH1F*)gDirectory->Get("DimuonMassVertexConstrained_bb");
	TH1F *hist_be = (TH1F*)gDirectory->Get("DimuonMassVertexConstrained_be");

	for (int p=0; p<hist->GetNbinsX() ;p++){
// 		std::cout<<hist->GetBinLowEdge(p)<<" --- "<<hist->GetBinContent(p)<<std::endl;
		if(hist->GetBinLowEdge(p) > 1300 && hist->GetBinLowEdge(p) < 1330 && hist->GetBinContent(p) > 0){
			std::cout<<" ===================== "<<std::endl;
			h_zz->SetBinContent(p, 0);
			h_zz_bb->SetBinContent(p, 0);
			h_zz_be->SetBinContent(p, 0);		
			continue;
		}
			
		h_zz->SetBinContent(p,  hist->GetBinContent(p));
		h_zz_bb->SetBinContent(p,  hist_bb->GetBinContent(p));
		h_zz_be->SetBinContent(p,  hist_be->GetBinContent(p));
	}

	TFile *f_zz_new = new TFile("/afs/cern.ch/work/f/ferrico/private/ZPrime_code/CMSSW_8_0_21/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/mc/ana_datamc_ZZ_ext_skim.root", "RECREATE");	
	TDirectory *dir = f_zz_new->mkdir("Our2016MuonsPlusMuonsMinusHistos");
	dir->cd();
// 	f_zz_new->cd("Our2016MuonsPlusMuonsMinusHistos");
	h_zz->Write();
	h_zz_bb->Write();
	h_zz_be->Write();
	f_zz_new->Close();
// 	
	std::cout<<h_zz->Integral()<<std::endl;
	std::cout<<h_zz_bb->Integral()<<std::endl;
	std::cout<<h_zz_be->Integral()<<std::endl;

}
