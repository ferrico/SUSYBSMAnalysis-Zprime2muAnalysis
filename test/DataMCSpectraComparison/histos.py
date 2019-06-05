#!/usr/bin/env python


miniAOD = True
Electrons = False

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =[#'file:./pat.root'
#         '/store/mc/RunIIAutumn18MiniAOD/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/120000/078DB2B1-40DD-634D-A3CF-D2E377CAFA48.root'

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/F6872BE7-B6FF-E811-B523-3417EBE528B5.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/E661E021-B7FF-E811-9221-A4BF0112BBDE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/B8F9E019-A702-E911-A9B2-0025905A613C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/A860E574-B7FF-E811-8F4A-0425C5DE7BEE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/A6CF84E0-7E01-E911-BECD-28924A33AF26.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/5483381F-B9FF-E811-9939-90B11CBCFF75.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/46A97DA6-B6FF-E811-94A8-A0369FE2C07E.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/34CEC3CA-B6FF-E811-8239-0CC47A7EEE96.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/22033CCC-A802-E911-9C36-44A84225C893.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/FCDAFAF4-E9FF-E811-9067-1866DAEECF18.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/D21267A3-0DFF-E811-B019-0CC47A4C8E0E.root',
'/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/844FF8D1-4FFF-E811-A7FD-00000931FE80.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/22A4ECE2-52FF-E811-98A0-0CC47A745294.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/ECB1A1F1-AF01-E911-8152-0CC47A4C8F30.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/CA5E0C0A-EAFF-E811-B933-001EC94B4EF3.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/96BCFB1B-C201-E911-B21C-0025904C66A2.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/5454F8C1-B001-E911-A13D-842B2B7682C7.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/50E40CAC-B201-E911-BEAB-008CFA197A70.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/3E2FE142-AF01-E911-816C-EC0D9A8221FE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/389722EA-AF01-E911-87FE-0CC47A13D2A4.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/361A4781-AF01-E911-8DC4-A0369FE2C0DE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/0AB8932F-B101-E911-AE3D-1CB72C1B6568.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/20000/0A0BBEF9-AF01-E911-A569-A4BF0112BD7C.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/D2652C30-F9FF-E811-9A82-A0369FC5B56C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/CC4C5B02-F9FF-E811-A257-0CC47A13D3A8.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/B4156155-F9FF-E811-AB7F-D4AE526DF7FF.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/905A2D74-FCFF-E811-A291-0025904C6620.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/904527DE-F8FF-E811-9ED8-E0071B74AC00.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/72892446-FAFF-E811-89D2-0242AC1C0504.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/38FE9F6F-FCFF-E811-94F8-0CC47A4C8F08.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/10F317C6-FCFF-E811-B8F5-20CF307C98DC.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/06E7EDC2-F9FF-E811-9364-0CC47A5FC67D.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/024AE64F-FCFF-E811-B06B-A4BF0108B90A.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/FCC3FBEB-FAFE-E811-8520-0CC47AF9B32A.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/BA8792B6-F9FE-E811-9871-48FD8EE73AE9.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/A63DB95A-FBFE-E811-840E-6C3BE5B52368.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/A24F67EA-FAFE-E811-A1E9-1CB72C1B2D88.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/7E1B86D8-FAFE-E811-BD29-90B11CBCFFEA.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/7CCFCCB8-F9FE-E811-97C4-24BE05C4D821.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6EDAE4F8-F9FE-E811-95EF-0CC47A4C8F18.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/4A403B97-FAFE-E811-A10D-C4346BC80410.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/445B98A3-F9FE-E811-A8A5-0CC47A0091C6.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/06BE50F2-FAFE-E811-B663-A4BF0112BD8E.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/AEE1011E-8AFE-E811-A460-801844DEEC30.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/9AA7BE87-89FE-E811-B1F8-001D0968F760.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/960B6F6B-8BFE-E811-ABB2-A0B3CCDFB624.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/8673FA6E-8BFE-E811-8D7E-0CC47A0093EC.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/7EDC2367-8BFE-E811-9B80-A0369FC5EEB0.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/24BCC472-8BFE-E811-8B26-0025907277CE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/1CA1D868-8BFE-E811-BE70-AC1F6B0DE0C4.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/0265FE66-8BFE-E811-997E-0CC47A5FBE21.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/D4C21D4D-8808-E911-B2E4-0025905B85FE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/B83883DD-8808-E911-BFA3-90E2BAD5729C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/98505B04-8808-E911-8E5B-A4BF0108B83A.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/76CF84F7-8608-E911-9976-0CC47A57D164.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/700F2FD2-8608-E911-ACF7-001E67A3EF48.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/6878F5D2-8608-E911-9AA5-0CC47AD98F72.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/5AD73B3C-8608-E911-8928-003048C446B8.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/4AC6D584-8808-E911-9E7B-484D7E8DF09F.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/38284F43-8808-E911-80A0-0CC47AFCC65E.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/26B8090A-8908-E911-8905-D4AE527EDC59.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/D6CAF73C-14FE-E811-ACA4-001E67792672.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/D6C5F95A-2BFE-E811-8542-A4BF0112BDE8.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/CCFEB521-24FE-E811-8D48-0025905A609E.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/A2B50D10-24FE-E811-9BE2-0CC47A4DEE50.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/00000/700FCFFD-24FE-E811-B0BE-0242AC130002.root',

# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/FA79AB8D-87FE-E811-88F0-246E96D14B5C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/F422CC95-39FE-E811-8E41-20040FE8E824.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/F04611E5-5BFF-E811-8E6C-20CF307C98DC.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E6911657-A3FF-E811-ACF8-008CFA0646D4.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E60C0CB5-11FF-E811-BDE1-008CFAF754BE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E28C7349-700A-E911-BB17-003048FFD75A.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/E01A6566-6C0A-E911-9F98-AC1F6B0DE30C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/DAB88C7A-700A-E911-9D73-3417EBE51CDF.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/D86BDF7A-0800-E911-AA11-A0369FC52374.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/D0D883C3-C1FE-E811-9A15-00238BAA1C4C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/D0B118B5-C0FF-E811-836D-0CC47AD98F78.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/D0718C63-B2FE-E811-A4EE-A4BF0112BD24.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/CEE8F0F9-F900-E911-B4C6-002590E7DE70.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/CA39B76A-41FE-E811-963E-A4BF0112BC7A.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/C4D2EE05-6F01-E911-9C9B-1866DAEA8808.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/BCEB6F2E-B802-E911-9F42-44A842CFD633.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/B8102909-E8FF-E811-8F50-0CC47A4C8F2C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/AC67079B-BAFE-E811-99E9-A4BF01125958.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/AC5A2808-B1FF-E811-847B-AC1F6B1AF186.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/AC29FD8D-4CFE-E811-B9D6-ECB1D7B67E10.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/9ED94F11-BBFE-E811-AAB9-801844DEEFD8.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/9E1C7516-790A-E911-98C6-A0369FC5DCB8.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/9AD5B91E-B5FE-E811-A8A5-A4BF01125B38.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/98B6134B-37FE-E811-A47F-0025905C43EC.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/98A67742-6C0A-E911-9A1F-0025905C53F2.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/947512A1-12FF-E811-A009-003048FFD7AA.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/92D45F5F-6C0A-E911-A0CF-28924A35059A.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/906A5E9D-C0FE-E811-AE46-A4BF01125A30.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/900B627D-AFFE-E811-9CF0-AC1F6B0F6758.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/8C282682-AFFE-E811-897A-A4BF0112BD24.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/8C0C50E4-790A-E911-9970-0242AC130002.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/8A39B5D2-4001-E911-A734-1CB72C1B2D80.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/80DE1CD8-71FF-E811-BFBB-0CC47A4DEDA2.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/765C8FDF-8200-E911-93EA-44A84225C911.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/760929D9-740A-E911-A93F-0242AC1C0500.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/72F26E53-8BFE-E811-87F8-A4BF0112BD4E.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6E13E852-37FE-E811-BFFB-0CC47AFCC3A2.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6AC2CC2B-8DFE-E811-9B2B-0025904CF766.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/6AAFEBA7-5EFF-E811-8444-0CC47AF9B302.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/68AF7BC1-790A-E911-AAE8-AC1F6B1AF186.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/684506A2-C9FF-E811-81C9-B8CA3A709648.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/680A798D-730A-E911-BE53-B8CA3A708F98.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/54172ACD-A4FE-E811-9F47-14187741121F.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/50A3D4D2-C6FE-E811-AE32-24BE05C6D731.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/4C9309E5-A000-E911-9A6E-001E67A3E8F9.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/40DB9224-B1FF-E811-8DD4-549F3525DEFC.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/4019248E-87FF-E811-9994-0CC47A7EEE0E.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/3AF92D6D-43FE-E811-A161-A4BF0112DC34.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/3A05E184-790A-E911-8DA5-6C3BE5B50170.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/2C556A19-3EFF-E811-A79F-002590FD5E88.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/2C47F6A6-00FF-E811-A9DB-001E67E6F88C.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/2AA736F4-3DFE-E811-AC11-001E677926FA.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/28F43855-CAFF-E811-8217-002481CFE708.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/280907CF-E1FE-E811-AB77-B083FED424C3.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/1A956BF3-790A-E911-9FAA-0025904B02AE.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/1A8CD5BA-7400-E911-A8A4-008CFA1C6458.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/0E0AE564-6C0A-E911-A649-0CC47A6C1874.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/0C098951-54FE-E811-8CE1-A0369FD0B130.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/06F28D30-47FE-E811-BDC7-001E67792810.root',
# '/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/06CCD709-53FF-E811-B5BF-28924A33BBAA.root',
			   ]
process.maxEvents.input = -1
process.GlobalTag.globaltag ='94X_mc2017_realistic_v14' #94X_mc2017_realistic_MuonTrackFix_01' # 2017 DY TRK FIXED
# process.GlobalTag.globaltag ='102X_upgrade2018_realistic_v12'
# process.GlobalTag.globaltag ='102X_upgrade2018_realistic_v15' # DY 102
# process.GlobalTag.globaltag ='102X_dataRun2_Prompt_v11' # DATA D
# process.GlobalTag.globaltag ='102X_dataRun2_Sep2018Rereco_v1' # reReco A,B,C
#process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold

# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale

process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

if miniAOD:
	from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
	electrons_miniAOD(process)

	from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
	HistosFromPAT.leptonsFromDileptons = True
	####################################
	####################################
	####################################
	
	HistosFromPAT.usekFactor = False #### Set TRUE to use K Factor #####

	####################################
	####################################
	####################################
	
	
	####################################
	####################################
	####################################
	
	ZSkim = False #### Set TRUE to skim dy50to120 with a Z pt < 100 GeV #####

	####################################
	####################################
	####################################
	
else:
	from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
	HistosFromPAT.leptonsFromDileptons = True
	
# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
#import SUSYBSMAnalysis.Zprime2muAnalysis.VBTFSelection_cff as VBTFSelection
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionOld_cff as OurSelectionOld
#import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2011EPS_cff as OurSelection2011EPS
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionNew_cff as OurSelectionNew
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff as OurSelectionDec2012
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2018_cff as OurSelection2018





# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
	('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
	('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
	('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
	]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
# 	'Our2012'  : OurSelectionDec2012,
	'Our2016'  : OurSelection2018,
	#'OurNoIso' : OurSelectionDec2012,
	'Simple'   : OurSelection2018, # The selection cuts in the module listed here are ignored below.
	#'OurMuPrescaledNew'  : OurSelectionNew,
	#'OurMuPrescaled2012' : OurSelectionDec2012
	}

if miniAOD and Electrons:
	dils = [('MuonsPlusMuonsMinus',          '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
		('MuonsPlusMuonsPlus',           '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',         'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
		('MuonsMinusMuonsMinus',         '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
		('MuonsSameSign',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
		('MuonsAllSigns',                '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',         ''),
		('MuonsPlusElectronsMinus',      '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == -2'),
		('MuonsMinusElectronsPlus',      '%(leptons_name)s:muons@- %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == 2'),
		('MuonsPlusElectronsPlus',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     'daughter(0).pdgId() + daughter(1).pdgId() == -24'),
		('MuonsMinusElectronsMinus',     '%(leptons_name)s:muons@- %(leptons_name)s:electrons@-',     'daughter(0).pdgId() + daughter(1).pdgId() == 24'),
		('MuonsElectronsOppSign',        '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',     ''),
		('MuonsElectronsSameSign',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
		('MuonsElectronsAllSigns',       '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',     ''),
		]
	
	cuts = {'Our2012'  : OurSelectionDec2012,
		'Our2016'  : OurSelection2018,
		'EmuVeto'  : OurSelectionDec2012, # this switches on the dRMuEl veto
		'Simple'   : OurSelectionDec2012, # The selection cuts in the module listed here are ignored below.
		}
	

# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
	# Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present.
    if miniAOD:
	    path_list.append(process.egmGsfElectronIDSequence)
	    
   
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'MuPrescaled' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut
    if miniAOD:
    	leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
    else:
    	leptons = process.leptons.clone(muon_cuts = muon_cuts)

    if  miniAOD and Electrons:
	    if cut_name == 'EmuVeto':
		    leptons.electron_muon_veto_dR = 0.1

    # Keep using old TuneP for past selections
#     if 'Dec2012' not in Selection.__file__:
#         leptons.muon_track_for_momentum = cms.string('TunePNew')
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the MuPrescaled paths, we don't care about e-mu events.
        if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
            continue

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)

        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        # Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        if cut_name == 'Simple':
            alldil.electron_cut_mask = cms.uint32(0)
            #alldil.loose_cut = 'isGlobalMuon && pt > 20.'#to be changed for first runs
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.sort_by_pt = True
            dil.do_remove_overlap = False
            if hasattr(dil, 'back_to_back_cos_angle_min'):
                delattr(dil, 'back_to_back_cos_angle_min')
            if hasattr(dil, 'vertex_chi2_max'):
                delattr(dil, 'vertex_chi2_max')
            if hasattr(dil, 'dpt_over_pt_max'):
                delattr(dil, 'dpt_over_pt_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
        elif 'MuPrescaled' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
            assert alldil.tight_cut == trigger_match
            alldil.tight_cut = prescaled_trigger_match

    # Histos now just needs to know which leptons and dileptons to use.
      
	histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)


       #define the list of MC samples to be read here. be careful that if WWinclusive or tautau sample are not commented it will apply the filters when running locally.

    samples = [
 		########## DY 2017 TRK NOT FIXED ###########
#           ('dy50to120', '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'),
#       ('dy120to200',    '/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#           ('dy200to400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#       ('dy400to800', '/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#           ('dy800to1400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#       ('dy1400to2300', '/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#       ('dy2300to3500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#       ('dy3500to4500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
#       ('dy4500to6000', '/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
		########## DY 2017 TRK NOT FIXED ###########
#		
		########## DY 2017 TRK FIXED ###########
		('dy50to120', '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
		('dy120to200', '/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
# 		('dy200to400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
# 		('dy400to800', '/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
		('dy800to1400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
# 		('dy1400to2300', '/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
# 		('dy2300to3500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
# 		('dy3500to4500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
# 		('dy4500to6000', '/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/RunIIFall17MiniAODv2-MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'),
		########## DY 2017 TRK FIXED ###########


# 	    ('dy50to120', '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#     	('dy120to200', '/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
# 	    ('dy200to400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#     	('dy400to800', '/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
# 	    ('dy800to1400', '/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'),
#     	('dy1400to2300', '/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#     	('dy2300to3500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#     	('dy3500to4500', '/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#     	('dy4500to6000', '/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),

#	    ('ZZ', '/ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#    	('WZ', '/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM'),
#		('WW', '/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'),
#	    ('Wantitop', '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'),
#    	('tW', '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'),
#    	('ttbar_lep', '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'),
#    	('dyInclusive50', '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'),


# 		('dyPt100to250', '/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext5-v1/MINIAODSIM'),
# 		('dyPt250to400', '/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext5-v1/MINIAODSIM'),
# 		('dyPt400to650', '/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 		('dyPt650', '/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),

#     	('WZ_ext', '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'),
#     	
#     	('ZZ_ext', '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'),
# 
#     	('WWinclusive', '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('WW200to600', '/WWTo2L2Nu_Mll_200To600_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('WW600to1200', '/WWTo2L2Nu_Mll_600To1200_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('WW1200to2500', '/WWTo2L2Nu_Mll_1200To2500_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('WW2500', '/WWTo2L2Nu_Mll_2500ToInf_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
# 
# 	    ('Wjets', '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
#     	('ttbar_lep50to500', '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 		('ttbar_lep_500to800', '/TTToLL_MLL_500To800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 		('ttbar_lep_800to1200', '/TTToLL_MLL_800To1200_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 		('ttbar_lep_1200to1800', '/TTToLL_MLL_1200To1800_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 		('ttbar_lep1800toInf', '/TTToLL_MLL_1800ToInf_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 
# 
# ###	    ('qcd50to80', '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('qcd80to120', '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('qcd120to170', '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('qcd170to300', '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('qcd300to470', '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('qcd470to600', '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('qcd600to800', '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('qcd800to1000', '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('qcd1000to1400', '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('qcd1400to1800', '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('qcd1800to2400', '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
#     	('qcd2400to3200', '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'),
# 	    ('qcd3200', '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM'),
		
        ]

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    if miniAOD:
        process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')
        process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
	pobj = process.EventCounter * process.dileptonPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)
    else:
    	pobj = process.muonPhotonMatch * reduce(lambda x,y: x*y, path_list)


    for name, ana_dataset in samples:
    
	if 'dy50to120' in name:
		if ZSkim:
			if miniAOD:
				process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
				process.DYGenMassFilter = cms.EDFilter('DyPt_ZSkim',
							       src = cms.InputTag('prunedGenParticles'),
							       min_mass = cms.double(0),
							       max_mass = cms.double(100), 
							       )
			else:
				process.DYGenMassFilter = cms.EDFilter('DyPt_ZSkim',
							       src = cms.InputTag('prunedMCLeptons'),
							       min_mass = cms.double(0),
							       max_mass = cms.double(100), 
							       )   
			pobj = process.DYGenMassFilter * pobj

	if 'ttbar_lep50to500' in name:
			if miniAOD:
				process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
				process.DYGenMassFilter = cms.EDFilter('TTbarSelection',
							       src = cms.InputTag('prunedGenParticles'),
							       min_mass = cms.double(50),
							       max_mass = cms.double(500), 
							       )
			else:
				process.DYGenMassFilter = cms.EDFilter('TTbarSelection',
							       src = cms.InputTag('prunedMCLeptons'),
							       min_mass = cms.double(50),
							       max_mass = cms.double(500), 
							       )   
			pobj = process.DYGenMassFilter * pobj
		
	if 'WWinclusive' in name:
		if miniAOD:
			process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
			process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
							       src = cms.InputTag('prunedGenParticles'),
							       min_mass = cms.double(50),
							       max_mass = cms.double(200), 
							       )
		else:
			process.DYGenMassFilter = cms.EDFilter('DibosonGenMass',
							       src = cms.InputTag('prunedMCLeptons'),
							       min_mass = cms.double(50),
							       max_mass = cms.double(200),
							       )   
		pobj = process.DYGenMassFilter * pobj
		     
	if 'dyInclusive50' in name:
		if miniAOD:
			process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
			process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
							       src = cms.InputTag('prunedGenParticles'),                                      
							       )
		else:
			process.DYGenMassFilter = cms.EDFilter('TauTauSelection',
							       src = cms.InputTag('prunedMCLeptons'),                                      
							       ) 
			
		pobj = process.DYGenMassFilter * pobj


	
    if 'VBTF' not in cut_name and cut_name != 'Simple':
        if not miniAOD:
		pobj = process.goodDataFilter * pobj
	else:
		process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
		for dataFilter in goodDataFiltersMiniAOD:
			#setattr(process,dataFilter 
			pobj = dataFilter * pobj


    if 'MuPrescaled' in cut_name: ####### Now it seams that there are no prescaled path ########
        if miniAOD:
            pobj = process.PrescaleToCommonMiniAOD * pobj ####### Now it seams that there are no prescaled path ########
        else:
            pobj = process.PrescaleToCommon * pobj ####### Now it seams that there are no prescaled path ########
    path = cms.Path(pobj)
    setattr(process, pathname, path)


def ntuplify(process, fill_gen_info=True):


    if miniAOD:
	process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
    	obj = process.prunedMCLeptons
	obj.src = cms.InputTag('prunedGenParticles')

	process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD',
                                           dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
						met_src = cms.InputTag("slimmedMETs"),
						jet_src = cms.InputTag("slimmedJets"),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
								TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
# 								TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
                                           genEventInfo = cms.untracked.InputTag('generator'),
                                           metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter"))
                                           )
 
    else:
	process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler',
                                           dimu_src = cms.InputTag('SimpleMuonsAllSigns'),
					   met_src = cms.InputTag("patMETsPF"),
					   jet_src = cms.InputTag("cleanPatJets"),
                                           beamspot_src = cms.InputTag('offlineBeamSpot'),
                                           vertices_src = cms.InputTag('offlinePrimaryVertices'),
					   TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),
                                           genEventInfo = cms.untracked.InputTag('generator')
                                           )
    if miniAOD and Electrons:
	    process.SimpleNtuplerEmu = process.SimpleNtupler.clone(dimu_src = cms.InputTag('SimpleMuonsElectronsAllSigns'))

    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        process.SimpleNtupler.hardInteraction = hardInteraction
        
    if hasattr(process, 'pathSimple'):
	if miniAOD and fill_gen_info:
        	process.pathSimple *=obj * process.SimpleNtupler 
		if Electrons:
			process.pathSimple *=obj * process.SimpleNtupler * process.SimpleNtuplerEmu
	else:
		process.pathSimple *= process.SimpleNtupler 
		if Electrons:
			process.pathSimple *= process.SimpleNtupler * process.SimpleNtuplerEmu

ntuplify(process) #to have ntuples also running in interactive way

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.triggerSummaryAnalyzerAOD

    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.pathSimple *= process.PrintOriginalMuons

    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.PrintEventSimple

    #- 2011-2012 selection (Nlayers > 8)
    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu

    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
    if hasattr(process, 'pathOur2012'):
        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
    process.CheckPrescale.trigger_paths = cms.vstring(*trigger_paths)
    process.pCheckPrescale = cms.Path(process.CheckPrescale)

def for_data(process):
	process.GlobalTag.globaltag ='94X_mc2017_realistic_v14' #94X_mc2017_realistic_MuonTrackFix_01' # 2017 DY TRK FIXED
# 	process.GlobalTag.globaltag ='102X_upgrade2018_realistic_v12'
# 	process.GlobalTag.globaltag ='102X_upgrade2018_realistic_v15' # DY 102
# 	process.GlobalTag.globaltag ='102X_dataRun2_Prompt_v11' # DATA D
# 	process.GlobalTag.globaltag ='102X_dataRun2_Sep2018Rereco_v1' # reReco A,B,C
	ntuplify(process)
    #check_prescale(process, trigger_paths) ####### Now it seams that there are no prescaled path ########

def for_mc(process, hlt_process_name, fill_gen_info):
    ntuplify(process, fill_gen_info)
    switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)

def get_dataset(run):
    #JMTBAD common with dataset_details in submit below, make a DataSamples.py?
    run = int(run)
    if 190450 <= run <= 191284:
        return '/SingleMu/tucker-datamc_SingleMuRun2012A_Prompt_190450_191284_20120418134612-57b19813ab8f2ab142c4566dc6738156/USER'
    else:
        raise ValueError('dunno how to do run %i' % run)

if 'int_data' in sys.argv:
    for_data(process)
    printify(process)
    
if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', False)
    printify(process)
    
if 'gogo' in sys.argv:
    for_data(process)
    printify(process)
    
    n = sys.argv.index('gogo')
    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
    print run, lumi, event
    run = int(run)
    lumi = int(lumi)
    event = int(event)
    filename = [x for x in sys.argv if x.endswith('.root')]
    if filename:
        filename = filename[0]
    else:
        dataset = get_dataset(run)
        print dataset
        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
        print repr(output)
        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
    print filename
    process.source.fileNames = [filename]
    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
    set_events_to_process(process, [(run, event)])

f = file('outfile_histos1', 'w')
f.write(process.dumpPython())
f.close()

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'ana_datamc_%(name)s'
config.General.workArea = 'crab_2017_DYTrkFixed_smeared_v3'
#config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
#config.JobType.priority = 1
config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'global'
job_control
# config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/ferrico'
config.Data.ignoreLocality = True 
config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_IT_Bari'
#config.Debug.scheddName = 'crab3-5@vocms059.cern.ch'
'''
    
    just_testing = 'testing' in sys.argv
        
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [
# 						('SingleMuonRun2018A_v1-PromptReco', '/SingleMuon/Run2018A-PromptReco-v1/MINIAOD'),
# 						('SingleMuonRun2018A_v2-PromptReco', '/SingleMuon/Run2018A-PromptReco-v2/MINIAOD'),
# 						('SingleMuonRun2018A_v3-PromptReco', '/SingleMuon/Run2018A-PromptReco-v3/MINIAOD'),
# 						('SingleMuonRun2018B_v1-PromptReco', '/SingleMuon/Run2018B-PromptReco-v1/MINIAOD'),
# 						('SingleMuonRun2018B_v2-PromptReco', '/SingleMuon/Run2018B-PromptReco-v2/MINIAOD'),
# 						('SingleMuonRun2018C_v1-PromptReco', '/SingleMuon/Run2018C-PromptReco-v1/MINIAOD'),
# 						('SingleMuonRun2018C_v2-PromptReco', '/SingleMuon/Run2018C-PromptReco-v2/MINIAOD'),
# 						('SingleMuonRun2018C_v3-PromptReco', '/SingleMuon/Run2018C-PromptReco-v3/MINIAOD'),

# 						('SingleMuonRun2018A_v2-17Sep2018', '/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD'),
# 						('SingleMuonRun2018B_v1-17Sep2018', '/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD'),
# 						('SingleMuonRun2018C_v1-17Sep2018', '/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD'),
						('SingleMuonRun2018D_v2-PromptReco', '/SingleMuon/Run2018D-PromptReco-v2/MINIAOD'),

            ]

        lumi_lists = [
		# 'NoLumiMask'
		#           'DCSOnly',
		#            'Run2012PlusDCSOnlyMuonsOnly',
		'Run2018MuonsOnly',
		# 'Run2015',
		]

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:
            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name

            new_py = open('histos.py').read()
            new_py += "\nfor_data(process)\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'LumiBased'
#config.Data.runRange = '256843-257490'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 150
config.Data.lumiMask = '%(lumi_mask)s' #######
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 100000
    ''')

       
        for name, ana_dataset in samples:
            print name

            new_py = open('histos.py').read()
            
            open('histos_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % locals())
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

#         if not just_testing:
# 		os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')

