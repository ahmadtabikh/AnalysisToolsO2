// To avoid your local _inputs.h files to be replaced by the template in the git repository, you should ask git to ignore them with git update-index --assume-unchanged (for more details see https://stackoverflow.com/questions/3319479/can-i-git-commit-a-file-and-ignore-its-content-changes)

// TFile* file_AliAnalysis = new TFile("../AnalysisResults_Run2_merged_Jaime.root");
TFile* file_AliAnalysis; //dummy

// //////// -------- Track Sel check rough -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_GlobalWithOptionalTpc", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalWithOptionalTpc", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";




// //////// -------- Track Sel check Detailed-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_intermTrackSel2", "LHC23zzh_intermTrackSel1bis", "LHC23zzh_intermTrackSel1", "LHC23zzh_GlobalTpcOptional", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "intermTrackSel2", "intermTrackSel1bis", "intermTrackSel1", "GlobalTpcOptional", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Track Sel check small LHCzzh -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMTracks", "LHC23zzh_GlobalTracks"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";



// //////// -------- Apass vs Cpasses comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Cpass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass1", "LHC23zzh_cpass11", "LHC23zzh_cpass10", "LHC23zzh_cpass8", "LHC23zzh_cpass1", "LHC23zzh_cpass0"};
// const TString DatasetsNames[nDatasets] = {"apass1", "cpass11", "cpass10", "cpass8", "cpass1", "cpass0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- apass3 vs apass 2 with sel8Full -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("apass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8Full", "LHC23zzh_apass3_sel8Full"};
// const TString DatasetsNames[nDatasets] = {"apass2", "apass3"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- sel8 vs sel8Full -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("EvtSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8", "LHC23zzh_apass2_sel8Full"};
// const TString DatasetsNames[nDatasets] = {"sel8", "sel8Full"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Track Sel check small LHC23k6d -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23k6d_MMTracks", "LHC23k6d_GlobalTracksTpcOptional", "LHC23k6d_GlobalTracks"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTpcOptional", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- LHC24b1 vs ppJetAnchored -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh_apass3");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass3", "LHC23zzh_apass3"};
// const TString DatasetsNames[nDatasets] = {"globalTracks", "uniformTracks"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency_uniformTracks"
//                                           };


// //////// -------- Run 2 data - O2Physics vs Aliphysics -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Software");
// const TString* texDatasetsComparisonCommonDenominator = new TString("Run2Data");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass3", "LHC23zzh_apass3"};
// const TString DatasetsNames[nDatasets] = {"Run3GlobalTracks", "Run3UniformTracks"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency_uniformTracks"
//                                           };

// //////// -------- Run 3 data - sigmapt mc vs data -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("sim/data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_global", "LHC23zzh_apass2_uniform", "LHC23k6d_global", "LHC23k6d_uniform"};
// const TString DatasetsNames[nDatasets] = {"data global", "data uniform", "mc global", "mc uniform"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency"
//                                           };


// //////// -------- Run 3 data apass3 1)-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("global/uniform");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform"};
// const TString DatasetsNames[nDatasets] = {"data global", "data uniform"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency_uniformTracks"
//                                           };

// //////// -------- Run 3 data apass3 2) -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("global/uniform");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 4;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// const TString DatasetsNames[nDatasets] = {"global", "uniform", "itsOnly", "globalWithBadTPCCrossedRows"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency"
//                                           };

// const TString trainId = "";
            

// //////// -------- Run 3 data apass4 - merged LHC23zzh -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("globalTracks");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 1;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_areaLeadingCuts_centralityWindowsAvailable"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};



// //////// -------- MC pp anchored PbPB vs pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("anchor");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet sim");
// const int nDatasets = 3;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_local"};
// const TString DatasetsNames[nDatasets] = {"PbPbAnchored", "ppAnchored", "LHC24f3 pp"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };
// // const TString analysisWorkflow[nDatasets] = {"track-efficiency",
// //                                             "track-efficiency"
// //                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true, false};

// jaime testJetAnchoredGapGen-20240801-155530


// //////// -------- Run 3 data apass4 - run by run -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_areaLeadingCuts/run544116", "LHC23zzh_apass4_areaLeadingCuts/run544123"};
// const TString DatasetsNames[nDatasets] = {"run 544116", "run 544123"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                             "track-efficiency"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false};



// //////// -------- MC pp anchored PbPB vs pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("includePartonEvent");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet sim");
// const int nDatasets = 2;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetPbPb_includePartonEvent"};
// const TString DatasetsNames[nDatasets] = {"includePartonEvent off", "includePartonEvent on"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true};

// // jaime testJetAnchoredGapGen-20240801-155530


// //////// -------- LHC23zzh - apass 4 withOUT area and leadingtrackpt cuts-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_areaLeadingCuts_AND_NoAreaNorLeadingCuts"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh_apass4_NoAreaNorLeadingCuts"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_noAreaLeadingJetCut"
// const TString analysisWorkflow[nDatasets] = {"track-efficiency"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};


// //////// -------- Run 3 data apass4 - run by run -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train254040/run544116", "LHC23zzh_apass4_train254040/run544121", "LHC23zzh_apass4_train254040/run544123"};
// const TString DatasetsNames[nDatasets] = {"run544116", "run544121", "run544123"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false};
// const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true};



// //////// -------- LHC23zzh - apass 4 occupancy analysis - sel8FullPbPb-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("occupancy");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8FullPbPb");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"sel8FullPbPb/LHC23zzh_apass4_occupancy01000", "sel8FullPbPb/LHC23zzh_apass4_occupancy05000", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-5k occupancy", "5k-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false, false, false, false};

// const std::string histDrawColorsOption = "";


// //////// -------- LHC23zzh - apass 4 occupancy analysis - sel8Full -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("occupancy");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8Full");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"sel8Full/LHC23zzh_apass4_occupancy01000", "sel8Full/LHC23zzh_apass4_occupancy05000", "sel8Full/LHC23zzh_apass4_occupancy10000", "sel8Full/LHC23zzh_apass4_occupancy20000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-5k occupancy", "5k-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false, false, false, false};

// const std::string histDrawColorsOption = "";




// //////// -------- LHC23zzh - apass 4 occupancy analysis - sel8Full vs sel8FullPbPb -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("sel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"sel8Full/LHC23zzh_apass4_occupancy01000", "sel8FullPbPb/LHC23zzh_apass4_occupancy01000", "sel8Full/LHC23zzh_apass4_occupancy05000", "sel8FullPbPb/LHC23zzh_apass4_occupancy05000", "sel8Full/LHC23zzh_apass4_occupancy10000", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000", "sel8Full/LHC23zzh_apass4_occupancy20000", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k sel8Full", "0k-1k sel8FullPbPb", "1k-5k sel8Full", "1k-5k sel8FullPbPb", "5k-10k sel8Full", "5k-10k sel8FullPbPb", "10-20k sel8Full", "10-20k sel8FullPbPb"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false, false, false, false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false, false, false, false, false, false, false, false};

// const std::string histDrawColorsOption = "colorPairs";


// //////// -------- LHC23zzh - apass 4 occupancy analysis - sel8Full vs sel8FullPbPb -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("sel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train254040/merged"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh apass4"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false};

// const std::string histDrawColorsOption = "colorPairs";



// //////// -------- LHC23zzh - apass 4 occupancy analysis - sel8FullPbPb-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("occupancy");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8FullPbPb");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"sel8FullPbPb/LHC23zzh_apass4_occupancy01000", "sel8FullPbPb/LHC23zzh_apass4_occupancy02000", "sel8FullPbPb/LHC23zzh_apass4_occupancy06000", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-2k occupancy", "2k-6k occupancy", "6-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false, false, false, false, false};

// const std::string histDrawColorsOption = "";

// //////// -------- local pp test sim anchored to PbPb comparison to pp anchor -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("ppAnchorPbPb");
// const TString* texDatasetsComparisonCommonDenominator = new TString("sim");
// const int nDatasets = 2;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"pp_sim_anchored_to_PbPb_5360GeV", "pp_sim_anchored_to_pp_1360GeV_jaime"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"PbPb anchor (5.36TeV)", "pp anchor (13.6TeV)"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa",
//                                             };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true};
// const bool trackHistsObsoleteVersion[nDatasets] = {false, false};
// const std::string histDrawColorsOption = "";


// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 6000occupancy - run(IR) comparison-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23 golden runs");
// const int nDatasets = 27;
// const TString Datasets[nDatasets] = {"run544124","run544123","run544477","run544392","run544492","run544391","run544122","run544476","run544390","run544098","run544454","run544475","run544121","run544032","run544491","run544095","run544389","run544451","run544510","run544474","run544185","run544091","run544028","run544184","run544116","run544508","run544490"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"6.6 kHz - run544124","11 kHz - run544123","13 kHz - run544477","13 kHz - run544392","15 kHz - run544492","15 kHz - run544391","16 kHz - run544122","16 kHz - run544476","18 kHz - run544390","18 kHz - run544098","19 kHz - run544454","19 kHz - run544475","22 kHz - run544121","23 kHz - run544032","24 kHz - run544491","25 kHz - run544095","27 kHz - run544389","28 kHz - run544451","29 kHz - run544510","29 kHz - run544474","29 kHz - run544185","29 kHz - run544091","30 kHz - run544028","32 kHz - run544184","38 kHz - run544116","39 kHz - run544508","43 kHz - run544490"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[10]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[11]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[12]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[13]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[14]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[15]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[16]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[17]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[18]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[19]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[20]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[21]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[22]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[23]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[24]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[25]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[26]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
// const std::string histDrawColorsOption = "";




// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_sel8FullPbPb_train345079","LHC23zzh_sel8FullPbPb_train345079"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"00-10%","50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_id24141",
//                                           "jet-finder-charged-qa_id24142"
//                                           };

// const TString wagonId[nDatasets] = {"_id24141",
//                                     "_id24142"
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false,false};
// const std::string histDrawColorsOption = "";
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false};



// //////// -------- LHC23zzh - apass 4 occupancy analysis sel8FullPbPb-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8FullPbPb");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"sel8FullPbPb/LHC23zzh_apass4_occupancy01000_train297793", "sel8FullPbPb/LHC23zzh_apass4_occupancy02000_train297794", "sel8FullPbPb/LHC23zzh_apass4_occupancy06000_train302311", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000_train297796", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000_train297797"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-2k occupancy", "2k-6k occupancy", "6k-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString wagonId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false};
// const std::string histDrawColorsOption = "";
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false};


// //////// -------- LHC23zzh - apass 4 occupancy analysis sel8FullPbPb-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8FullPbPb");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"sel8FullPbPb/LHC23zzh_apass4_occupancy01000", "sel8FullPbPb/LHC23zzh_apass4_occupancy05000", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-5k occupancy", "5k-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString wagonId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const std::string histDrawColorsOption = "";
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false};
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- LHC23zzh - apass 4 - 1000 occupancy cut with new event sel, and tighter eta (0.8 for tracks, 0.6 for jets - central 00-10% -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHCzzh_apass4_jetspectra_train365500","LHCzzh_apass4_jetspectra_train365500"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"00-10%","50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"jet-spectra-charged_central",
//                                              "jet-spectra-charged_peripheral"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false,false};
// const std::string histDrawColorsOption = "";
// const bool datasetsAreSubsetsofId0 = false;
// const bool trackHistsObsoleteVersion[nDatasets] = {false,false};



// //////// -------- LHC25b6 - pp sim anchored to PbPb 10% ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_10percent_train370115", "LHC25b4a_pp_ref_Tracks_train371362"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet MC Pb-Pb anchor", "ppRef gen.purp. MC"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false};
// const std::string histDrawColorsOption = "";
// const bool datasetsAreSubsetsofId0 = false;
// const bool trackHistsObsoleteVersion[nDatasets] = {true,true};




// /////////// -------- Run 3 LHC24_ppref_pass1 - run by run - Data/MC comparaison LHC24_ppref_pass1 vs LHC25b4ab6_train549962 ------ ////////
TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV"); 
const TString* texDatasetsComparisonType = new TString("run");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24_ppref_pass1");
const int nDatasets = 8;
const TString Datasets[nDatasets] = {"Run559456", "Run559437", "Run559444", "Run559443", "Run559423", "Run559408", "Run559410", "Run559409"};

const TString DatasetsNames[nDatasets] = {
                                        // "All runs wo 559348 and 559385",
                                        "1095 kHz - Run559456",
                                        "1101 kHz - Run559437",
                                        "1112 kHz - Run559444",
                                        "1112 kHz - Run559443",
                                        "1115 kHz - Run559423",
                                        "1118 kHz - Run559408",
                                        "1120 kHz - Run559410",
                                        "1146 kHz - Run559409"
                                        };

// excute the macro in the "yyyymmdd_discription" directroies (/home/tabikh/alice/MyWork/yyyymmdd_discription)
// TFile* file_O2Analysis_list[nDatasets] = {
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[0]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[1]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[2]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[3]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[4]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[5]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[6]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[7]+"/AnalysisResults.root")
//                                 };
// const TString analysisWorkflow[nDatasets] = {
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency"
//                                             };

// ////////////////////////////////// Use these **instead of the above** for taking the eta and phi track spectrum from the 2D et-phi map in the jet-sprectra-charged workflow ////////
TFile* file_O2Analysis_list[nDatasets] = { 
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[0]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[1]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[2]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[3]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[4]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[5]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[6]+"/AnalysisResults.root"),
                                new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[7]+"/AnalysisResults.root")
                                };
const TString analysisWorkflow[nDatasets] = {
                                            "jet-spectra-charged",
                                            "jet-spectra-charged",
                                            "jet-spectra-charged",
                                            "jet-spectra-charged",
                                            "jet-spectra-charged",
                                            "jet-spectra-charged",
                                            "jet-spectra-charged",
                                            "jet-spectra-charged"
                                            };                                           

// ############################################################################ MONTE CARLO DATASETS ############################################################################
//////////////////////////////for MC LHC25b4ab6_train549962 ////////////////////////////////
const TString MC_Datasets[nDatasets] = {"Run559456", "Run559437", "Run559444", "Run559443", "Run559423", "Run559408", "Run559410", "Run559409"};

TFile* file_O2AnalysisMC_list[nDatasets] = {new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[0]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[1]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[2]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[3]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[4]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[5]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[6]+"/AnalysisResults.root"), 
                                            new TFile("../Datasets/LHC25b4b5_train533385/"+MC_Datasets[7]+"/AnalysisResults.root")
                                          };

// const TString analysisWorkflowMC[nDatasets] = {"track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency"
//                                             };

const TString analysisWorkflowMC[nDatasets] = {"jet-spectra-charged","jet-spectra-charged","jet-spectra-charged","jet-spectra-charged","jet-spectra-charged","jet-spectra-charged","jet-spectra-charged","jet-spectra-charged"
                                            };
const TString wagonId[nDatasets] = {"", "", "", "", "", "", "", ""};
// const TString wagonId[nDatasets] = {"", "", ""};
const TString trainId = "";
const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false, false, false, false};
// const bool isDatasetWeighted[nDatasets] = {false, false, false};

const std::string histDrawColorsOption = "";
const bool datasetsAreSubsetsofId0 = false;
const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true, true ,true, true, true, true};
// const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true};

const std::string histDatasetComparisonStructure = "";








// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                            

// // //////// -------- Run 3 LHC24_ppref_medium test outout Test -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24_ppref_medium_test output");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {""};

// const TString DatasetsNames[nDatasets] = {"All runs"};

// // excute the macro in the "yyyymmdd_discription" directroies (/home/tabikh/alice/MyWork/yyyymmdd_discription)
// TFile* file_O2Analysis_list[nDatasets] = {
//                                 new TFile("../Datasets/LHC24_ppref_outputTestTracksQA/"+Datasets[0]+"/AnalysisResults.root")};

// const TString analysisWorkflow[nDatasets] = {
//                                             "track-efficiency"};
// const TString wagonId[nDatasets] = {""};
// // const TString wagonId[nDatasets] = {"", "", ""};

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// // const bool isDatasetWeighted[nDatasets] = {false, false, false};

// const std::string histDrawColorsOption = "";
// const bool datasetsAreSubsetsofId0 = false;
// const bool trackHistsObsoleteVersion[nDatasets] = {true};
// // const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true};

// const std::string histDatasetComparisonStructure = "";




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// // // /////////// -------- Run 3 LHC24_ppref_pass1 - run by run - Data/MC comparaison LHC24_ppref_pass1 vs LHC25b4a5_train558017_tracksQA ------ ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24_ppref_pass1");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"Run559348", "Run559361", "Run559362", "Run559385", "Run559387"
//                                      }; 

// const TString DatasetsNames[nDatasets] = {
//                                         // "All runs wo 559348 and 559385",
//                                         "490 kHz - Run559348",
//                                         "1109 kHz - Run559361",
//                                         "1107 kHz - Run559362",
//                                         "1101 kHz - Run559385",
//                                         "1107 kHz - Run559387"
//                                         };

// // excute the macro in the "yyyymmdd_discription" directroies (/home/tabikh/alice/MyWork/yyyymmdd_discription)
// TFile* file_O2Analysis_list[nDatasets] = {
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[0]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[1]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[2]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[3]+"/AnalysisResults.root"),
//                                 new TFile("../Datasets/LHC24_ppref_pass1_train387507_tracksQA/"+Datasets[4]+"/AnalysisResults.root")
//                                 };
// // TFile* file_O2Analysis_list[nDatasets] = {
// //                                 new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[0]+"/AnalysisResults.root"),
// //                                 new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[1]+"/AnalysisResults.root"),
// //                                 new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[2]+"/AnalysisResults.root"),
// //                                 new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[3]+"/AnalysisResults.root"),
// //                                 new TFile("../Datasets/LHC24_ppref_pass1_train380686/"+Datasets[4]+"/AnalysisResults.root")
// //                                 };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency"
//                                             };
// // const TString analysisWorkflow[nDatasets] = {"jet-spectra-charged", "jet-spectra-charged", "jet-spectra-charged", "jet-spectra-charged", "jet-spectra-charged"
// //                                             };  

// //////////////////////////////for MC LHC25b4a5_train558017_tracksQA ////////////////////////////////
// const TString MC_Datasets[nDatasets] = {"Run559348", "Run559361", "Run559362", "Run559385", "Run559387"
//                                         };                     
// TFile* file_O2AnalysisMC_list[nDatasets] = {new TFile("../Datasets/LHC25b4a5_train558017_tracksQA/"+MC_Datasets[0]+"/AnalysisResults.root"), 
//                                             new TFile("../Datasets/LHC25b4a5_train558017_tracksQA/"+MC_Datasets[1]+"/AnalysisResults.root"), 
//                                             new TFile("../Datasets/LHC25b4a5_train558017_tracksQA/"+MC_Datasets[2]+"/AnalysisResults.root"),
//                                             new TFile("../Datasets/LHC25b4a5_train558017_tracksQA/"+MC_Datasets[3]+"/AnalysisResults.root"), 
//                                             new TFile("../Datasets/LHC25b4a5_train558017_tracksQA/"+MC_Datasets[4]+"/AnalysisResults.root")
//                                           };

// // TFile* file_O2AnalysisMC_list[nDatasets] = {new TFile("../Datasets/LHC25b4a5_train535170/"+MC_Datasets[0]+"/AnalysisResults.root"), 
// //                                             new TFile("../Datasets/LHC25b4a5_train535170/"+MC_Datasets[1]+"/AnalysisResults.root"), 
// //                                             new TFile("../Datasets/LHC25b4a5_train535170/"+MC_Datasets[2]+"/AnalysisResults.root"),
// //                                             new TFile("../Datasets/LHC25b4a5_train535170/"+MC_Datasets[3]+"/AnalysisResults.root"), 
// //                                             new TFile("../Datasets/LHC25b4a5_train535170/"+MC_Datasets[4]+"/AnalysisResults.root")
// //                                           };

// const TString analysisWorkflowMC[nDatasets] = {"track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency"
//                                             };
// // const TString analysisWorkflowMC[nDatasets] = {"jet-spectra-charged", "jet-spectra-charged", "jet-spectra-charged", "jet-spectra-charged", "jet-spectra-charged"
// //                                             };                                            
// const TString wagonId[nDatasets] = {"", "", "", "", ""};
// // const TString wagonId[nDatasets] = {"", "", ""};
// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false};
// // const bool isDatasetWeighted[nDatasets] = {false, false, false};

// const std::string histDrawColorsOption = "";
// const bool datasetsAreSubsetsofId0 = false;
// const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true, true ,true};
// // const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true};

// const std::string histDatasetComparisonStructure = "";