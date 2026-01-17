
// To avoid your local _inputs.h files to be replaced by the template in the git repository, you should ask git to ignore them with git update-index --assume-unchanged (for more details see https://stackoverflow.com/questions/3319479/can-i-git-commit-a-file-and-ignore-its-content-changes)

#ifndef JETSPECTRUM_INPUTS_H
#define JETSPECTRUM_INPUTS_H

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////       file access choice       ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
TFile* file_O2Analysis_run2ComparisonFileHannaBossiLaura = new TFile("Datasets/Run2_Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange/Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange.root");
TFile* file_O2Analysis_run2ComparisonFileMLPaper = new TFile("Datasets/Run2_Unfolding_MachineLearningMethod_R020/Ch-jetSuppression_PbPb502TeV.root");


//////// -------- pp LHC24ap_pass1_small_train379582 as data with LHC25b4a as MC  ///////
TString* texEnergy = new TString("#sqrt{#it{s}} = 5.36 TeV"); 

TString* texCollisionDataType = new TString("Data pp"); 
TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+" "+(TString)*texEnergy); 

TString* texCollisionMCType = new TString("MC pp PYTHIA"); 
TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+" "+(TString)*texEnergy); 

const TString* texDatasetsComparisonType = new TString("0");
const TString* texDatasetsComparisonCommonDenominator = new TString("");
const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC24_ppref_pass1_train380686"};
// const TString DatasetsMC[1] = {"LHC25b4b5_train533385"};
const TString DatasetsMC[1] = {"LHC25b4b5_train533385"};
// const TString Datasets[nDatasets] = {"LHC24_ppref_train587361"};
const TString DatasetsNames[nDatasets] = {"ppref"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("../Datasets/"+Datasets[0]+"/AnalysisResults.root")};

// TFile* file_O2Analysis_MCfileForMatrix = nullptr;
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4ab6_train564527_trackLoss3perCent/AnalysisResults.root");
TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4b5_train533385/AnalysisResults.root"); // previously i used this LHC25b6_train380971 is MCweighted --> in settings set : const bool mcIsWeighted = true;
TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("../Datasets/LHC25b4a5_train535170/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)


// TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b6_train405048/AnalysisResults.root"); // previously i used this LHC25b6_train380971 is MCweighted --> in settings set : const bool mcIsWeighted = true;
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("../Datasets/LHC25b6_train405048/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// ///// !!!!!!! if controlMC is true use this ; if not comment it and uncomment the two line above 
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b6_McControl_train405048DividedByTwo/merged_AnalysisResults1.root"); 
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("../Datasets/LHC25b6_McControl_train405048DividedByTwo/merged_AnalysisResults2.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// const TString trainId = "_id12832";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString trainId = "_id12436";
// const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
const TString trainIdData = "";
const TString analysisWorkflowData = "jet-spectra-charged"+trainIdData;

const TString trainIdBkg = "";
const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;

const TString trainIdUnfoldingControl = "";
const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

const TString trainIdMC = "";
const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;

const bool oldformatNEventsGenInO2Analysis = true;
// bool oldformatNEventsGenInO2Analysis;
const bool etaCutOnMatchedJetsIsObsoleteVersion = false;


// add two matching matrix from different MC : LHC25b4a & LHC25b4b
const bool addTwoMatchingMatrixFromDifferentMC = false; 
//the first one will be the : file_O2Analysis_MCfileForMatrix
TFile* file_O2Analysis_MCfileForSecondMatrix = new TFile("../Datasets/LHC25b4b_train403772_MC/AnalysisResults.root");


// //////// -------- pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o pass7");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC22o_pass7_train261733"};
// const TString DatasetsNames[nDatasets] = {"LHC22o_pass7"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/LHC24f3b_train261768/half1/AnalysisResults.root");

// // TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24f3_sel8MC_train240962/AnalysisResults.root")};
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC24f3b_train261768/half2/AnalysisResults.root")};

// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/LHC24b1b_sel8Full_train239409/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted/AnalysisResults.root");

// const TString trainId = "";
// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";
// // const TString analysisWorkflowMC = "jet-finder-charged-qa_global_CollMatch";
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;





// //////// -------- pp spectrum with Joonsuk files tests -------- ////////
// TString* texEnergy = new TString("#sqrt{#it{s}} = 13.6 TeV"); 
// TString* texCollisionDataType = new TString("Data pp"); 
// TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+" "+(TString)*texEnergy); 
// TString* texCollisionMCType = new TString("MC pp PYTHIA"); 
// TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+" "+(TString)*texEnergy); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"Data_LHC22o_minBias_small_Joonsuk"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC24f3b_withTrackTuner_Joonsuk/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/ppSim_LHC24f3b_withTrackTuner_Joonsuk/AnalysisResults.root")};

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;

// const TString trainId = "";
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;
// const TString analysisWorkflow_unfoldingControl = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;







// //////// -------- Pb-Pb spectrum with Wenhui Angantyr files tests -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"Data_halfMCAngantyr_train361349"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr_train361349/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/MC_halfMCAngantyr_train361349/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdData = "";
// const TString analysisWorkflowData = "jet-spectra-charged"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;








// //////// -------- LHC25b6 - pp sim anchored to PbPb 10%  - just looking at eff and purity ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV - Pb-Pb anchor"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_10percent_train366181_jetSpectrumWorkflow"};
// const TString DatasetsNames[nDatasets] = {"LHC25b6"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_10percent_train366181_jetSpectrumWorkflow/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_10percent_train366181_jetSpectrumWorkflow/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdData = "";
// const TString analysisWorkflowData = "jet-spectra-charged"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;




// //////// -------- test LHC25b6 efficiency finderQA vs spectraCharged ///////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC25b6_localWorkflow"};
// const TString DatasetsNames[nDatasets] = {"spectraCharged"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix = new TFile("Datasets/LHC25b6_localWorkflow/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl = {new TFile("Datasets/LHC25b6_localWorkflow/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdData = "_id26156";
// const TString analysisWorkflowData = "jet-spectra-charged_central"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged_noOccupancyCut"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = false;


#endif
