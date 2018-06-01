#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.roottools import ROOT

path = './data_OK/NoScale_YesEtaCut_Run2016MuonsOnly/ana_datamc_data.root'

tmp_fn = 'ca.txt'
# tmp_fn = 'OurListNoDop_BE_Scaled.txt'
# tmp_fn = 'OurListNoDop_BE_NoScaled.txt'

branch_spec = 'vertex_m:lep_pt[0]:lep_pt[1]:lep_eta[0]:lep_phi[0]:lep_eta[1]:lep_phi[1]:run:lumi:event'#:vertex_chi2:lep_isTrackerMuon[0]:lep_isTrackerMuon[1]:lep_isGlobalMuon[0]:lep_isGlobalMuon[1]:GoodDataRan:GoodVtx:lep_dB[0]:lep_dB[1]:lep_triggerMatchPt[0]:lep_triggerMatchPt[1]:lep_glb_numberOfValidMuonHits[0]:lep_glb_numberOfValidMuonHits[1]:lep_glb_numberOfValidPixelHits[0]:lep_glb_numberOfValidPixelHits[1]:lep_glb_numberOfValidTrackerLayers[0]:lep_glb_numberOfValidTrackerLayers[1]:lep_id[0]:lep_id[1]:lep_sumPt[0]/lep_tk_pt[0]:lep_sumPt[1]/lep_tk_pt[1]:lep_pt_err[0]/lep_pt[0]:lep_pt_err[1]/lep_pt[1]:cos_angle'
 

cut='fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && vertex_m > 780 && vertex_m < 785'
# cut = 'vertex_m < 1061.1 && vertex_m > 1061.'


# cut='fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && GoodDataRan && GoodVtx && lep_pt[0]>53. && lep_pt[1]>53. && lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && ( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && cos_angle>-0.9998 && lep_id[0]*lep_id[1]<0 && (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && vertex_chi2 < 20 && 1299 < vertex_m && vertex_m < 1300.5'
#fabs(lep_eta[0])<1.2 && fabs(lep_eta[1])<1.2 && 
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
prev_run = prev_lumi = prev_event = 0
cleaned_lines = []
for x in lines:
    cleaned_line = []
    curr_run   = x[0]
    curr_lumi  = x[1]
    curr_event = x[2]
    mass = x[3]
    pt_1 = x[4]
    pt_2 = x[5]
    eta_1 = x[6]
    phi_1 = x[7]
    eta_2 = x[8]
    phi_2 = x[9]
    #    print curr_run, curr_lumi, curr_event, mass
    if curr_run != prev_run or curr_lumi != prev_lumi or curr_event != prev_event:
#         cleaned_line.append( curr_run, curr_lumi, curr_event, mass)
        cleaned_line.append(mass, pt_1, pt_2, eta_1, phi_1, eta_2, phi_2, curr_run, curr_lumi, curr_event)
        cleaned_lines.append(cleaned_line)
        prev_run   = curr_run
        prev_lumi  = curr_lumi
        prev_event = curr_event
    else:
        print 'delete line: ',curr_run,' ',curr_lumi,' ',curr_event,' mass ',mass
lines = ['\t'.join(y.strip() for y in x) for x in cleaned_lines]
open(tmp_fn, 'wt').write('\n'.join(lines))

# f = ROOT.TFile(path.replace('.root', '.microtuple.root'), 'CREATE')
# t = ROOT.TTree('t','')
# t.ReadFile(tmp_fn, 'dil_mass')
# f.Write()
# f.Close()

#os.remove(tmp_fn)
