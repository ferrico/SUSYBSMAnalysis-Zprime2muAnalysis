#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

#path = 'zp2mu_histos.root'

# path = './data/Run2016MuonsOnly/ana_datamc_data.root'
path = './data_OK/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root'

# path = 'zp2mu_histos_TunePNew_Scale.root'
# path = 'zp2mu_histos_TunePNew_NoScale.root'
# path = './mc/ana_datamc_ZZ_ext.root'

#path = 'zp2mu_ana_datamc_data_RunG_AOD/Run2016G/Run2016MuonsOnly/ana_datamc_data.root'
#path = 'zp2mu_ana_datamc_data_RunG/Run2016MuonsOnly/ana_datamc_data.root'
#path = 'zp2mu_ana_datamc_data_RunG_ReReco/Run2016MuonsOnly/ana_datamc_data.root'
#tmp_fn = 'micro_ntuple.temp_runG_ReReco.txt'
#tmp_fn = 'micro_ntuple.temp_runG.txt'
#tmp_fn = 'micro_ntuple.temp_runH.txt'
#tmp_fn = 'micro_ntuple.temp_MINIAOD.txt'
#tmp_fn = 'micro_ntuple.temp_AOD.txt'

# tmp_fn = 'micro_ntuple.temp_NoScale_TUTTO.txt'
tmp_fn = 'cancel.txt'
#branch_spec = 'run:lumi:event:dil_mass:vertex_chi2:lep_triggerMatchPt[0]:lep_triggerMatchPt[1]'

# branch_spec = 'vertex_m:dil_mass:run:lumi:event:lep_pt[0]:lep_pt[1]:lep_eta[0]:lep_eta[1]:lep_phi[0]:lep_phi[1]'
# branch_spec = 'run:lumi:event:dil_mass:lep_pt[1]:lep_eta[1]:lep_phi[1]:lep_qOverPt[1]:lep_pt[0]:lep_eta[0]:lep_phi[0]:lep_qOverPt[0]'
branch_spec = 'run:lumi:event:dil_mass:lep_eta[0]:lep_eta[1]'
#branch_spec = 'run:lumi:event:vertex_m' 

#branch_spec = 'run:lumi:event:vertex_m:lep_pt[1]:lep_pt[0]:lep_sumPt[1]:lep_tk_pt[1]:lep_sumPt[0]:lep_tk_pt[0]'#dil_mass
#branch_spec = 'run:lumi:event:vertex_m'
#branch_spec = 'run:lumi:event:dil_mass:lep_pt[0]:lep_pt[1]:lep_eta[0]:lep_eta[1]:triggerMatched:trigger_match_0:trigger_match_1:lep_tk_numberOfValidPixelHits[0]:lep_tk_numberOfValidPixelHits[1]:lep_glb_numberOfValidPixelHits[0]:lep_glb_numberOfValidPixelHits[1]:lep_triggerMatchEta[0]:lep_triggerMatchEta[1]'
#branch_spec = 'run:lumi:event:dil_mass:lep_pt[0]:lep_pt[1]:lep_eta[0]:lep_eta[1]:triggerMatched:trigger_match_0:trigger_match_1'
#cut = 'vertex_m>400'
#cut = 'event ==  56306745 && lep_isGlobalMuon[0] && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_isTrackerMuon[0] && lep_pt[0] > 45 && lep_pt[1] > 45 && abs(lep_dB[0]) < 0.2 &&  abs(lep_dB[1]) < 0.2  && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 &&  GoodData && OppSign && vertex_m>50'
########## 2 TeV
#cut = 'event ==  56306745 && lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 45 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 45 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && ((lep_triggerMatchPt[0] > 45) || (lep_triggerMatchPt[1] > 45)) && cos_angle > -0.9998   && vertex_chi2 < 10 && GoodData && OppSign && vertex_m>50 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 '

#cut = 'event==429190704'
#cut = 'lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 53 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 53 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[1] / lep_tk_pt[1] < 0.1  && (lep_id[0]*lep_id[1]<0)  && cos_angle > -0.9998 && (GoodDataRan  && GoodVtx) && vertex_m > 900 && ((lep_triggerMatchPt[0] > 50) || (lep_triggerMatchPt[1] > 50)) && vertex_chi2 < 20'

#cut = 'lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 53 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 53 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[1] / lep_tk_pt[1] < 0.1  && (lep_id[0]*lep_id[1]<0)  && cos_angle > -0.9998 && (GoodDataRan  && GoodVtx) && dil_mass > 900 && ((lep_triggerMatchPt[0] > 50) || (lep_triggerMatchPt[1] > 50)) && vertex_chi2 >= 20'

#2016
# cut = 'lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && dil_chosen==0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20'

#2012
#cut = 'lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && lep_numberOfMatchedStations[0] > 1 && lep_numberOfMatchedStations[1] > 1 && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && dil_chosen==0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 900'
# cut='fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20'
cut='event <= 278808'
# (fabs(lep_eta[0])>1.2 || fabs(lep_eta[1])>1.2) && 


### OR ov trigger filters ### offline cut only pt>45
#cut = 'run==251643 && lep_isGlobalMuon[0] && lep_isTrackerMuon[0] && lep_pt[0] > 50 && abs(lep_dB[0]) < 0.2 && lep_glb_numberOfValidTrackerLayers[0] > 5 && lep_glb_numberOfValidPixelHits[0] >= 1 && lep_glb_numberOfValidMuonHits[0] > 0 && lep_numberOfMatchedStations[0] > 1 && lep_pt_err[0] / lep_pt[0] < 0.3 && lep_sumPt[0] / lep_tk_pt[0] < 0.1 && lep_isGlobalMuon[1] && lep_isTrackerMuon[1] && lep_pt[1] > 50 && abs(lep_dB[1]) < 0.2 && lep_glb_numberOfValidTrackerLayers[1] > 5 && lep_glb_numberOfValidPixelHits[1] >= 1 && lep_glb_numberOfValidMuonHits[1] > 0 && lep_numberOfMatchedStations[1] > 1 && lep_pt_err[1] / lep_pt[1] < 0.3 && lep_sumPt[1] / lep_tk_pt[1] < 0.1 && ((lep_triggerMatchPt[0] > 50 ) || (lep_triggerMatchPt[1] > 50 )) && cos_angle > -0.9998   && vertex_chi2 < 10 && GoodData && OppSign && vertex_m>50'
#&& vertex_m>60 && vertex_m<120
#cut = 'loose_2012_0 && loose_2012_1'
#cut = 'OurSel2012 && vertex_m>60'

f = ROOT.TFile(path)
#t = f.SimpleNtuplertunepnew.Get('t')
#t = f.SimpleNtuplerstartup.Get('t')
t = f.SimpleNtupler.Get('t')
t.GetPlayer().SetScanRedirect(True)
t.GetPlayer().SetScanFileName(tmp_fn)
t.Scan(branch_spec, cut,  "colsize=15 precision=15")
t.GetPlayer().SetScanRedirect(False)
f.Close()

lines = [line.split(' *')[1:] for line in open(tmp_fn).readlines() if ' * ' in line and 'Row' not in line]
lines = ['\t'.join(y.strip() for y in x) for x in lines]
open(tmp_fn, 'wt').write('\n'.join(lines))

# f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'RECREATE')
# t = ROOT.TTree('t','')
# t.ReadFile(tmp_fn, branch_spec)
# f.Write()
# f.Close()

#os.remove(tmp_fn)
