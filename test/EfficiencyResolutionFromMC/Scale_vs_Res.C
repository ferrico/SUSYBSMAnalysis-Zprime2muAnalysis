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



void Scale_vs_Res(){

	float bias = gRandom->Gaus(0,0.025);
	float MuonPt = 1500;
	float res_pt;
	float scale_pt;


	for(int i = 0; i < 3500; i++){

		scale_pt = MuonPt/1000.;
		scale_pt = 1/scale_pt;
		scale_pt = scale_pt + bias;
		scale_pt = 1/scale_pt;
		scale_pt = scale_pt*1000.;
		
		res_pt = gRandom->Gaus(MuonPt, 0.01 * MuonPt);
		
		std::cout<<"Nominal = "<<MuonPt<<"\t scale = "<<scale_pt<<"\t resolution = "<<res_pt<<
		"\t\t Nominal / Scale = "<<MuonPt/scale_pt - 1<<"\t\t Nominal / resolution = "<<MuonPt/res_pt -1 <<std::endl;
		
		MuonPt++;
	
	}

}