minimalVariables = ["CMS_hgg_mass[160,100,180]:=diPhoton().mass",
                    "dZ[40,-20.,20.]:=(tagTruth().genPV().z-diPhoton().vtx().z)", # store actual value
                                                                               #when doing systematics, variables need to have a binning
                                                                               #specified, otherwise the rooDataHist end up empty.
            								       #an assert in the code prevents you from doing this.
                    "centralObjectWeight[1,-999999.,999999.] := centralWeight"]

minimalHistograms = []

minimalNonSignalVariables = ["CMS_hgg_mass[160,100,180]:=diPhoton().mass"]#,"centralObjectWeight[1,-999999.,999999.] := centralWeight"]

defaultVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass", 
                                    "leadPt                   :=diPhoton().leadingPhoton.pt",
                                    "subleadPt                :=diPhoton().subLeadingPhoton.pt",
                                    "leadEta                  :=diPhoton().leadingPhoton.eta",
                                    "subleadEta               :=diPhoton().subLeadingPhoton.eta",
                                    "diphoMVA                 :=diPhotonMVA().result",
                                    "leadIDMVA                :=leadingView.phoIdMvaWrtChosenVtx",
                                    "subleadIDMVA             :=subLeadingView.phoIdMvaWrtChosenVtx",
                                    "maxEta                   :=max(abs(diPhoton().leadingPhoton.superCluster.eta),abs(diPhoton().leadingPhoton.superCluster.eta))",
                                    "vtxZ           :=diPhoton().vtx().z",
                                     "bjet1_pt := ?bJets.size>0? bJets.at(0).pt : -1",
                        "bjet2_pt := ?bJets.size>1? bJets.at(1).pt : -1",
                        "bjet1_eta:= ?bJets.size>0? bJets.at(0).eta: -1",
                        "bjet2_eta:= ?bJets.size>1? bJets.at(1).eta:-1",
                        "bjet1_abseta := ?bJets.size>0? abs(bJets.at(0).eta): -1",
                        "bjet2_abseta := ?bJets.size>1? abs(bJets.at(1).eta): -1",
                        "jet1_pt  := ?(jets.size>0)? jets.at(0).pt : -1",
                        "jet2_pt  := ?(jets.size>1)? jets.at(1).pt : -1",
                        "jet1_eta := ?(jets.size>0)? jets.at(0).eta : -1",
                        "jet2_eta := ?(jets.size>1)? jets.at(1).eta: -1",
                        "jet1_abseta := ?(jets.size>0)? abs(jets.at(0).eta): -1",
                        "jet2_abseta := ?(jets.size>1)? abs(jets.at(1).eta): -1",
                        "electronjet1deta := ?(electrons.size>0)? abs(electrons.at(0).eta - jets.at(0).eta) : -1",
                        "electronjet2deta := ?(electrons.size>0)? abs(electrons.at(0).eta - jets.at(1).eta) : -1",
                        "muonjet1deta := ?(muons.size>0)? abs(muons.at(0).eta - jets.at(0).eta) : -1",
"muonjet2deta := ?(muons.size>0)? abs(muons.at(0).eta - jets.at(1).eta) : -1",
                                    ]


defaultHistograms=["CMS_hgg_mass>>mass(160,100,180)",
                                     "subleadPt:leadPt>>ptLeadvsSub(180,20,200:180,20,200)",
                                     "diphoMVA>>diphoMVA(50,0,1)",
                                     "maxEta>>maxEta[0.,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4442,1.566,1.7,1.8,2.,2.2,2.3,2.5]"
                                     ]

systematicVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass"]#,"centralObjectWeight[1,-999999.,999999.] := centralWeight"]
systematicHistograms=["CMS_hgg_mass>>mass(160,100,180)"]
