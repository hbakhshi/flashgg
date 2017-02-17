#import flashgg.Taggers.globalVariables as globalVars

#import FWCore.ParameterSet.Config as cms

#from globalVariables_cff import globalVariables

vtx_variables=[
    "vtxprob                := diPhotonMVA.vtxprob",
    "ptbal                  := diPhoton.ptBal",
    "ptasym                 := diPhoton.ptAsym",
    "logspt2                := diPhoton.logSumPt2",
    "p2conv                 := diPhoton.pullConv",
    "nconv                  := diPhoton.nConv",
    "vtxmva                 := diPhoton.vtxProbMVA",
    "vtxdz                  := diPhoton.dZ1",
    "vtx_x                  := diPhoton.vtx.x", 
    "vtx_y                  := diPhoton.vtx.y", 
    "vtx_z                  := diPhoton.vtx.z", 
    "gv_x                   := diPhoton.genPV.x", 
    "gv_y                   := diPhoton.genPV.y", 
    "gv_z                   := diPhoton.genPV.z"

]

dipho_variables=[
    "dipho_sumpt            := diPhoton.sumPt",
    "dipho_cosphi           := abs(cos(diPhoton.leadingPhoton.phi - diPhoton.subLeadingPhoton.phi))",
    "dipho_mass             := diPhoton.mass",
    "dipho_pt               := diPhoton.pt",
    "dipho_phi              := diPhoton.phi",
    "dipho_eta              := diPhoton.eta",
    "dipho_PtoM             := diPhoton.pt/diPhoton.mass",
    "cosphi                 := diPhotonMVA.CosPhi",
    "sigmaMrvoM             := diPhotonMVA.sigmarv",
    "sigmaMwvoM             := diPhotonMVA.sigmawv",
]

photon_variables=[
    "dipho_leadPt           := diPhoton.leadingPhoton.pt",
    "dipho_leadEt           := diPhoton.leadingPhoton.et",
    "dipho_leadEta          := diPhoton.leadingPhoton.eta",
    "dipho_leadPhi          := diPhoton.leadingPhoton.phi",
    "dipho_lead_sieie       := diPhoton.leadingPhoton.sigmaIetaIeta",
    "dipho_lead_hoe         := diPhoton.leadingPhoton.hadronicOverEm",
    "dipho_lead_sigmaEoE    := diPhoton.leadingPhoton.sigEOverE",
    "dipho_lead_ptoM        := diPhoton.leadingPhoton.pt/diPhoton.mass",
    "dipho_leadR9           := diPhoton.leadingPhoton.r9",
    "dipho_leadIDMVA        := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
    "dipho_lead_elveto      := diPhoton.leadingPhoton.passElectronVeto",
    "dipho_lead_prompt      := diPhoton.leadingPhoton.genMatchType",
    "dipho_lead_chiso       := diPhoton.leadingView.pfChIso03WrtChosenVtx",
    "dipho_lead_chisow      := diPhoton.leadingPhoton.pfChgIsoWrtWorstVtx04",
    "dipho_lead_phoiso      := diPhoton.leadingPhoton.pfPhoIso03",
    "dipho_lead_phoiso04    := diPhoton.leadingPhoton.pfPhoIso04",
    "dipho_lead_neutiso      := diPhoton.leadingPhoton.pfNeutIso03",
    "dipho_lead_ecaliso03   := diPhoton.leadingPhoton.ecalRecHitSumEtConeDR03",
    "dipho_lead_hcaliso03   := diPhoton.leadingPhoton.hcalTowerSumEtConeDR03",
    "dipho_lead_pfcluecal03 := diPhoton.leadingPhoton.ecalPFClusterIso",
    "dipho_lead_pfcluhcal03 := diPhoton.leadingPhoton.hcalPFClusterIso",
    "dipho_lead_trkiso03    := diPhoton.leadingPhoton.trkSumPtHollowConeDR03",
    "dipho_lead_pfchiso2    := diPhoton.leadingView.pfChIso02WrtChosenVtx",
    "dipho_lead_haspixelseed:= diPhoton.leadingPhoton.hasPixelSeed",
    "dipho_lead_sieip       := diPhoton.leadingPhoton.sieip",
    "dipho_lead_etawidth    := diPhoton.leadingPhoton.superCluster.etaWidth",
    "dipho_lead_phiwidth    := diPhoton.leadingPhoton.superCluster.phiWidth",
    "dipho_lead_regrerr     := diPhoton.leadingPhoton.sigEOverE * diPhoton.leadingPhoton.energy",
    "dipho_lead_s4ratio     :=  diPhoton.leadingPhoton.s4",
    "dipho_lead_effSigma    :=  diPhoton.leadingPhoton.esEffSigmaRR",
    "dipho_lead_scraw       :=  diPhoton.leadingPhoton.superCluster.rawEnergy",
    "dipho_lead_ese         :=  diPhoton.leadingPhoton.superCluster.preshowerEnergy",

    "dipho_subleadPt        := diPhoton.subLeadingPhoton.pt",
    "dipho_subleadEt        := diPhoton.subLeadingPhoton.et",
    "dipho_subleadEta       := diPhoton.subLeadingPhoton.eta",
    "dipho_subleadPhi       := diPhoton.subLeadingPhoton.phi",
    "dipho_sublead_sieie    := diPhoton.subLeadingPhoton.sigmaIetaIeta",
    "dipho_sublead_hoe      := diPhoton.subLeadingPhoton.hadronicOverEm",
    "dipho_sublead_sigmaEoE := diPhoton.subLeadingPhoton.sigEOverE",
    "dipho_sublead_ptoM     := diPhoton.subLeadingPhoton.pt/diPhoton.mass",
    "dipho_subleadR9        := diPhoton.subLeadingPhoton.r9",
    "dipho_subleadIDMVA     := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
    "dipho_sublead_elveto   := diPhoton.subLeadingPhoton.passElectronVeto",
    "dipho_sulead_prompt    := diPhoton.subLeadingPhoton.genMatchType",
    "dipho_sublead_chiso       := diPhoton.leadingView.pfChIso03WrtChosenVtx",
    "dipho_sublead_chisow   := diPhoton.subLeadingPhoton.pfChgIsoWrtWorstVtx04",
    "dipho_sublead_phoiso   := diPhoton.subLeadingPhoton.pfPhoIso03",
    "dipho_sublead_phoiso04 := diPhoton.subLeadingPhoton.pfPhoIso04",
    "dipho_sublead_neutiso   := diPhoton.subLeadingPhoton.pfNeutIso03",
    "dipho_sublead_ecaliso03:= diPhoton.subLeadingPhoton.ecalRecHitSumEtConeDR03",
    "dipho_sublead_hcaliso03:= diPhoton.subLeadingPhoton.hcalTowerSumEtConeDR03",
    "dipho_sublead_pfcluecal03 := diPhoton.subLeadingPhoton.ecalPFClusterIso",
    "dipho_sublead_pfcluhcal03 := diPhoton.subLeadingPhoton.hcalPFClusterIso",
    "dipho_sublead_trkiso03    := diPhoton.subLeadingPhoton.trkSumPtHollowConeDR03",
    "dipho_sublead_pfchiso2    := diPhoton.subLeadingView.pfChIso02WrtChosenVtx",
    "dipho_sublead_haspixelseed:= diPhoton.subLeadingPhoton.hasPixelSeed",
    "dipho_sublead_sieip       := diPhoton.subLeadingPhoton.sieip",
    "dipho_sublead_etawidth    := diPhoton.subLeadingPhoton.superCluster.etaWidth",
    "dipho_sublead_phiwidth    := diPhoton.subLeadingPhoton.superCluster.phiWidth",
    "dipho_sublead_regrerr     := diPhoton.subLeadingPhoton.sigEOverE * diPhoton.subLeadingPhoton.energy",
    "dipho_sublead_s4ratio     :=  diPhoton.subLeadingPhoton.s4",
    "dipho_sublead_effSigma    :=  diPhoton.subLeadingPhoton.esEffSigmaRR",
    "dipho_sublead_scraw       :=  diPhoton.subLeadingPhoton.superCluster.rawEnergy",
    "dipho_sublead_ese         :=  diPhoton.subLeadingPhoton.superCluster.preshowerEnergy",
        
]

lepton_variables=[
    "LeptonType             :=getLeptonType()",
    "n_ele                  := electrons.size",
    "ele1_pt                := ?(electrons.size>0)? electrons.at(0).pt : -999",
    "ele2_pt                := ?(electrons.size>1)? electrons.at(1).pt : -999",
    "ele1_eta               := ?(electrons.size>0)? electrons.at(0).superCluster().eta : -999",
    "ele2_eta               := ?(electrons.size>1)? electrons.at(1).superCluster().eta : -999",
    "ele1_phi               := ?(electrons.size>0)? electrons.at(0).superCluster().phi : -999",
    "ele2_phi               := ?(electrons.size>1)? electrons.at(1).superCluster().phi : -999",
    "ele1_ch                := ?(electrons.size>0)? electrons.at(0).charge : -999",
    "ele2_ch                := ?(electrons.size>1)? electrons.at(1).charge : -999",
    "ele1_sigmaIetaIeta     := ?(electrons.size>0)? electrons.at(0).full5x5_sigmaIetaIeta : -999",
    "ele2_sigmaIetaIeta     := ?(electrons.size>1)? electrons.at(1).full5x5_sigmaIetaIeta : -999",
    "ele1_dEtaInSeed        := ?(electrons.size>0)? electrons.at(0).deltaEtaSuperClusterTrackAtVtx() - electrons.at(0).superCluster().eta() + electrons.at(0).superCluster().seed().eta  : -999",
    "ele2_dEtaInSeed        := ?(electrons.size>1)? electrons.at(1).deltaEtaSuperClusterTrackAtVtx() - electrons.at(1).superCluster().eta() + electrons.at(1).superCluster().seed().eta  : -999",
    "ele1_dPhiIn            := ?(electrons.size>0)? electrons.at(0).deltaPhiSuperClusterTrackAtVtx : -999",
    "ele2_dPhiIn            := ?(electrons.size>1)? electrons.at(1).deltaPhiSuperClusterTrackAtVtx : -999",
    "ele1_hOverE            := ?(electrons.size>0)? electrons.at(0).hadronicOverEm : -999",
    "ele2_hOverE            := ?(electrons.size>1)? electrons.at(1).hadronicOverEm : -999",
    "ele1_RelIsoEA          := ?(electrons.size>0)? ( electrons.at(0).pfIsolationVariables().sumChargedHadronPt + max(0.0, electrons.at(0).pfIsolationVariables().sumNeutralHadronEt + electrons.at(0).pfIsolationVariables().sumPhotonEt - getElecAlpha(0)*getrho) ) / electrons.at(0).pt : -999",
    "ele2_RelIsoEA          := ?(electrons.size>1)? ( electrons.at(1).pfIsolationVariables().sumChargedHadronPt + max(0.0, electrons.at(1).pfIsolationVariables().sumNeutralHadronEt + electrons.at(1).pfIsolationVariables().sumPhotonEt - getElecAlpha(1)*getrho) ) / electrons.at(1).pt : -999",
    "ele1_ooEmooP           := ?(electrons.size>0)? abs(1.0 - electrons.at(0).eSuperClusterOverP)*(1./electrons.at(0).ecalEnergy) : -999",
    "ele2_ooEmooP           := ?(electrons.size>1)? abs(1.0 - electrons.at(1).eSuperClusterOverP)*(1./electrons.at(1).ecalEnergy) : -999",
    "ele1_dxy               := ?(electrons.size>0)? getLeadingLeptonVertexDxy(\"electron\"): -999",
    "ele2_dxy               := ?(electrons.size>0)? getSubLeadingLeptonVertexDxy(\"electron\"):  -999",
    "ele1_diphodxy          := ?(electrons.size>0)? getLeadingLeptonDiphoVertexDxy(\"electron\"):  -999",
    "ele2_diphodxy          := ?(electrons.size>1)? getSubLeadingLeptonDiphoVertexDxy(\"electron\"):  -999",
    "ele1_dz                := ?(electrons.size>0)? getLeadingLeptonVertexDz(\"electron\") : -999",
    "ele2_dz                := ?(electrons.size>1)? getSubLeadingLeptonVertexDz(\"electron\"):  -999",
    "ele1_diphodz           := ?(electrons.size>0)? getLeadingLeptonDiphoVertexDz(\"electron\"):  -999",
    "ele2_diphodz           := ?(electrons.size>1)? getSubLeadingLeptonDiphoVertexDz(\"electron\"):  -999",
    "ele1_misHits           := ?(electrons.size>0)? getLeadingElectronMisHits: -999",
    "ele2_misHits           := ?(electrons.size>0)? getSubLeadingElectronMisHits: -999",
    "ele1_ConversionVeto    := ?(electrons.size>0)? electrons.at(0).passConversionVeto : -999",
    "ele2_ConversionVeto    := ?(electrons.size>1)? electrons.at(1).passConversionVeto : -999",
    "ele1_ChargedHadronPt   := ?(electrons.size>0)? electrons.at(0).pfIsolationVariables().sumChargedHadronPt: -999",
    "ele2_ChargedHadronPt   := ?(electrons.size>1)? electrons.at(1).pfIsolationVariables().sumChargedHadronPt: -999",
    "ele2_NeutralHadronEt   := ?(electrons.size>1)? electrons.at(1).pfIsolationVariables().sumNeutralHadronEt: -999",
    "ele1_NeutralHadronEt   := ?(electrons.size>0)? electrons.at(0).pfIsolationVariables().sumNeutralHadronEt: -999",
    "ele1_PhotonEt          := ?(electrons.size>0)? electrons.at(0).pfIsolationVariables().sumPhotonEt: -999",
    "ele2_PhotonEt          := ?(electrons.size>1)? electrons.at(1).pfIsolationVariables().sumPhotonEt: -999",
    "n_muons                := muons.size",
    "muon1_pt               := ?(muons.size>0)? muons.at(0).pt : -999",
    "muon2_pt               := ?(muons.size>1)? muons.at(1).pt : -999",
    "muon1_eta              := ?(muons.size>0)? muons.at(0).eta : -999",
    "muon2_eta              := ?(muons.size>1)? muons.at(1).eta : -999",
    "muon1_phi              := ?(muons.size>0)? muons.at(0).phi : -999",
    "muon2_phi              := ?(muons.size>1)? muons.at(1).phi : -999",
    "muon1_ch               := ?(muons.size>0)? muons.at(0).charge : -999",
    "muon2_ch               := ?(muons.size>1)? muons.at(1).charge : -999",
    "muon1_iso              := ?(muons.size>0)? (muons.at(0).pfIsolationR04().sumChargedHadronPt+ max( 0.,muons.at(0).pfIsolationR04().sumNeutralHadronEt + muons.at(0).pfIsolationR04().sumPhotonEt - 0.5 * muons.at(0).pfIsolationR04().sumPUPt)) / ( muons.at(0).pt ) : -999.",
    "muon2_iso              := ?(muons.size>1)? (muons.at(1).pfIsolationR04().sumChargedHadronPt+ max( 0.,muons.at(1).pfIsolationR04().sumNeutralHadronEt + muons.at(1).pfIsolationR04().sumPhotonEt - 0.5 * muons.at(1).pfIsolationR04().sumPUPt)) / ( muons.at(1).pt ) : -999.",
    "muon1_chi2             := ?(muons.size>0)? muons.at(0).innerTrack().normalizedChi2 : -999",
    "muon2_chi2             := ?(muons.size>1)? muons.at(1).innerTrack().normalizedChi2 : -999",
    "muon1_mHits            := ?(muons.size>0)? muons.at(0).innerTrack().hitPattern().numberOfValidMuonHits : -999",
    "muon2_mHits            := ?(muons.size>1)? muons.at(1).innerTrack().hitPattern().numberOfValidMuonHits : -999",
    "muon1_mStations        := ?(muons.size>0)? muons.at(0).numberOfMatchedStations : -999",
    "muon2_mStations        := ?(muons.size>1)? muons.at(1).numberOfMatchedStations : -999",
    "muon1_dxy              := ?(muons.size>0)? getLeadingLeptonVertexDxy(\"muon\"): -999",
    "muon2_dxy              := ?(muons.size>0)? getSubLeadingLeptonVertexDxy(\"muon\"):  -999",
    "muon1_diphodxy         := ?(muons.size>0)? getLeadingLeptonDiphoVertexDxy(\"muon\"):  -999",
    "muon2_diphodxy         := ?(muons.size>1)? getSubLeadingLeptonDiphoVertexDxy(\"muon\"):  -999",
    "muon1_dz               := ?(muons.size>0)? getLeadingLeptonVertexDz(\"muon\") : -999",
    "muon2_dz               := ?(muons.size>1)? getSubLeadingLeptonVertexDz(\"muon\"):  -999",
    "muon1_diphodz          := ?(muons.size>0)? getLeadingLeptonDiphoVertexDz(\"muon\"):  -999",
    "muon2_diphodz          := ?(muons.size>1)? getSubLeadingLeptonDiphoVertexDz(\"muon\"):  -999",
    "muon1_pxHits           := ?(muons.size>0)? muons.at(0).innerTrack().hitPattern().numberOfValidPixelHits : -999",
    "muon2_pxHits           := ?(muons.size>1)? muons.at(1).innerTrack().hitPattern().numberOfValidPixelHits : -999",
    "muon1_tkLayers         := ?(muons.size>0)? muons.at(0).innerTrack().hitPattern().trackerLayersWithMeasurement : -999",
    "muon2_tkLayers         := ?(muons.size>1)? muons.at(1).innerTrack().hitPattern().trackerLayersWithMeasurement : -999"

    

]
jet_variables=[
    
    "n_fwdjets              := Jets_EtaSorted.size",
    "fwdjet1_pt             := ?Jets_EtaSorted.size>0? Jets_EtaSorted.at(0).pt : -999",
    "fwdjet2_pt             := ?Jets_EtaSorted.size>1? Jets_EtaSorted.at(1).pt : -999",
    "fwdjet1_eta            := ?Jets_EtaSorted.size>0? Jets_EtaSorted.at(0).eta: -999.",
    "fwdjet2_eta            := ?Jets_EtaSorted.size>1? Jets_EtaSorted.at(1).eta:-999.",
    "fwdjet1_phi            := ?Jets_EtaSorted.size>0? Jets_EtaSorted.at(0).phi: -999.",
    "fwdjet2_phi            := ?Jets_EtaSorted.size>1? Jets_EtaSorted.at(1).phi: -999.",

    "n_M_bjets                := nMedium_bJets",
    "n_L_bjets                := nLoose_bJets",
    "n_T_bjets                := nTight_bJets",

    "n_bjets                := bJets.size",
    "bjet1_pt               := ?bJets.size>0? bJets.at(0).pt : -999",
    "bjet2_pt               := ?bJets.size>1? bJets.at(1).pt : -999",
    "bjet1_eta              := ?bJets.size>0? bJets.at(0).eta: -999.",
    "bjet2_eta              := ?bJets.size>1? bJets.at(1).eta:-999.",
    "bjet1_phi              := ?bJets.size>0? bJets.at(0).phi: -999.",
    "bjet2_phi              := ?bJets.size>1? bJets.at(1).phi: -999.",

    # new variables
    "n_jets                 := jets.size",
    "jet1_pt                := ?jets.size>0? jets.at(0).pt : -999",
    "jet2_pt                := ?jets.size>1? jets.at(1).pt : -999",
    "jet3_pt                := ?jets.size>2? jets.at(2).pt : -999",
    "jet1_eta               := ?jets.size>0? jets.at(0).eta : -999",
    "jet2_eta               := ?jets.size>1? jets.at(1).eta : -999",
    "jet3_eta               := ?jets.size>2? jets.at(2).eta : -999",
    "jet1_phi               := ?jets.size>0? jets.at(0).phi : -999",
    "jet2_phi               := ?jets.size>1? jets.at(1).phi : -999",
    "jet3_phi               := ?jets.size>2? jets.at(2).phi : -999",
    
]

thqmva_variables=[
    "bTagWeight             := bTagWeight",
    "photonWeights          := photonWeights",
    "FoxWolf                :=getFoxWolframMoment_ONE",
    "Aplanarity             :=getAplanarity()",
    "MET                    :=getMET()",
    "METPhi                 :=getMET_Phi()"
]

for label in ["HighestBTagVal", "Medium" , "Loose" , "Tight"]:
    thqmva_variables.append('fwdJetEta_{0}             := ?thqleptonicMvaRes("{0}")>-10.? getFwdJet("{0}").eta : -999'.format(label) )
    thqmva_variables.append('MVA_{0}                   := thqleptonicMvaRes("{0}")'.format(label) )
    thqmva_variables.append('bJetPt_{0}                := ?thqleptonicMvaRes("{0}")>-10.? getbJet("{0}").pt : -999'.format(label) )


truth_variables=[
    ]
