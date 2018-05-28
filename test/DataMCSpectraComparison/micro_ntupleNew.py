#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

ROOT.gStyle.SetPaintTextFormat("%.3f");

# path = '/afs/cern.ch/work/f/ferrico/private/PickEvent_2017/CMSSW_9_2_0/src/pickevents_RunB3_MINIAOD.root'

path = './data/Run2017MuonsOnly/ana_datamc_data.root'
# path = './mc/mc_YesEtaCut_NoScale/MC_OK/ana_datamc_ttbar_lep.root'

# path = 'zp2mu_histos.root'

tmp_fn = 'cancella.txt'
# tmp_fn = 'List_400_600.txt'
# tmp_fn = 'OurListNoDop_BE_Scaled.txt'
# tmp_fn = 'OurListNoDop_BE_NoScaled.txt'

branch_spec = 'run:lumi:event:vertex_m:max(abs(lep_eta[0]), abs(lep_eta[1])):lep_pt[0]+lep_pt[1]'


# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 120 && vertex_m < 400'
# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 400 && vertex_m < 600'
# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 600 && vertex_m < 900'
# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 900 && vertex_m < 1300'
# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 1300 && vertex_m < 1800'
# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 1800 && vertex_m < 8000'
# cut='dil_chosen==0 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 780 && vertex_m < 785'
cut='GoodVtx && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20'# && vertex_m > 400'# && vertex_m < 900'


# cut='fabs(lep_eta[0])<1.2 && fabs(lep_eta[1])<1.2 && fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20'
#(fabs(lep_eta[0])>1.2 || fabs(lep_eta[1])>1.2)

#cut='lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && lep_numberOfMatchedStations[0] > 1 && lep_numberOfMatchedStations[1] > 1 && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m>900'

f = ROOT.TFile(path)
t = f.SimpleNtupler.Get('t')
t.GetPlayer().SetScanRedirect(True)
t.GetPlayer().SetScanFileName(tmp_fn)
#- colsize=10 is needed in order not to currupt the event number
t.Scan(branch_spec, cut, "colsize=15 precision=10")
#t.Scan('*', cut)
t.GetPlayer().SetScanRedirect(False)
f.Close()

lines = [line.split(' *')[1:] for line in open(tmp_fn).readlines() if ' * ' in line and 'Row' not in line]
#- This loop is to keep only the first (highest-rank) dimuon in an event
prev_run = prev_lumi = prev_event = sum_pt = 0
cleaned_lines = []
for x in lines:
    cleaned_line = []
    curr_run   = x[0]
    curr_lumi  = x[1]
    curr_event = x[2]
    mass = x[3]
    eta = x[4]
    curr_sum_pt = x[5]


    if curr_run != prev_run or curr_lumi != prev_lumi or curr_event != prev_event:
        cleaned_line.append(mass)
        cleaned_lines.append(cleaned_line)
        prev_run   = curr_run
        prev_lumi  = curr_lumi
        prev_event = curr_event
        prev_sum_pt = curr_sum_pt
        if(float(mass) > 200 and float(mass) < 400):
        	print mass, curr_run, curr_lumi, curr_event, eta, curr_sum_pt
#        	print mass, curr_run, curr_lumi, curr_event, eta, curr_sum_pt
    else:
    	if(curr_sum_pt > prev_sum_pt):
	        print 'delete curr: ',curr_run,' ',curr_lumi,' ',curr_event,' mass ',mass, curr_sum_pt
    	else:
	    	print 'delete pre: ',curr_run,' ',curr_lumi,' ',curr_event,' mass ',mass, prev_sum_pt
lines = ['\t'.join(y.strip() for y in x) for x in cleaned_lines]
open(tmp_fn, 'wt').write('\n'.join(lines))

# f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'CREATE')
# t = ROOT.TTree('t','')
# t.ReadFile(tmp_fn, 'dil_mass')
# f.Write()
# f.Close()

#os.remove(tmp_fn)
