#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
#include "THnSparse.h"
// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./TrackQC_settings.h"
#include "./TrackQC_inputs.h"
#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramPlotting.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

// Plot Utilities
TString contextTrackDatasetComp(std::string options);

//////////// QC plot functions
// Dataset comparison
void Draw_Pt_DatasetComparison(std::string options);
void Draw_Eta_DatasetComparison(float* ptRange, std::string options);
void Draw_Phi_DatasetComparison(float* ptRange, std::string options);

// void Draw_Eta_DatasetComparison_trackSelComp();
// void Draw_Phi_DatasetComparison_trackSelComp();


void Draw_Pt_CentralityComparison(int iDataset);
void Draw_Eta_CentralityComparison(int iDataset);
void Draw_Phi_CentralityComparison(int iDataset);

void Draw_Pt_Run2Run3Comparison_0010Cent(int iDataset);
void Draw_Eta_Run2Run3Comparison_0010Cent(int iDataset);
void Draw_Phi_Run2Run3Comparison_0010Cent(int iDataset);

void Draw_Sigmapt_vs_pt_DatasetComp();

void Draw_Sigmapt_vs_pt_DatasetComp_legacyTH3();
void Draw_Sigmapt_nonGlobal_uniformTracks_legacyTH3();
void Draw_Sigmapt_nonGlobal_uniformTracks_fromSubtraction_legacyTH3();
void Draw_Sigmapt_nonGlobal_uniformTracks_centralEta_legacyTH3();
void Draw_Sigmapt_vs_pt_DatasetComp_centralEta_legacyTH3();

void Draw_SelectedMultiplicity_DatasetComp();
void Draw_SelectedMultiplicity_CentralityComp(int iDataset);

void Draw_Mean_Ntrack_vs_Dataset();
void Draw_Mean_Pt_vs_Dataset();
void Draw_Mean_Eta_vs_Dataset();
void Draw_DcaXY_DatasetComp();

void Draw_tracks_2D_profile(std::string options, int refIndex);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void TrackQC() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");



  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    // Draw_Pt_CentralityComparison(iDataset);
    // Draw_Eta_CentralityComparison(iDataset);
    // Draw_Phi_CentralityComparison(iDataset);
  }
  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_Pt_Run2Run3Comparison_0010Cent(iDataset);
  //   Draw_Eta_Run2Run3Comparison_0010Cent(iDataset);
  //   Draw_Phi_Run2Run3Comparison_0010Cent(iDataset);
  // }

  // const int nPtBins = 10;
  const int nPtBins = 1;
  float jetPtMinCut, jetPtMaxCut;
  // float jetPtMinCutArray[nPtBins+1] = {0, 1, 2, 4, 6, 8, 10, 15, 20, 30, 200};
  float jetPtMinCutArray[nPtBins+1] = {0, 200};

  // Draw_tracks_2D_profile("etaphi", 0) ;
  Draw_Pt_DatasetComparison("evtNorm");
  // Draw_Pt_DatasetComparison("entriesNorm");
  for(int iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    jetPtMinCut = jetPtMinCutArray[iPtBin];
    jetPtMaxCut = jetPtMinCutArray[iPtBin+1];

    float ptRange[2] = {jetPtMinCut, jetPtMaxCut};
    // Draw_Eta_DatasetComparison(ptRange, "evtNorm");
    // Draw_Eta_DatasetComparison(ptRange, "entriesNorm");
    // Draw_Phi_DatasetComparison(ptRange, "evtNorm");
    // Draw_Phi_DatasetComparison(ptRange, "entriesNorm");

  // Draw_Eta_DatasetComparison_trackSelComp();
  // Draw_Phi_DatasetComparison_trackSelComp();
  }

  // Draw_Sigmapt_vs_pt_DatasetComp();
  // Draw_Sigmapt_nonGlobal_uniformTracks();
  // Draw_Sigmapt_nonGlobal_uniformTracks_fromSubtraction();
  // Draw_Sigmapt_nonGlobal_uniformTracks_centralEta();

  // Draw_SelectedMultiplicity_DatasetComp();
  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_SelectedMultiplicity_CentralityComp(iDataset);
  // }

  // Draw_Mean_Pt_vs_Dataset();
  // Draw_Mean_Eta_vs_Dataset();
  // Draw_Mean_Ntrack_vs_Dataset();
}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs() {
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  // gSystem->Load("libMinuit");
  // gSystem->Load("libSTEERBase");
  // gSystem->Load("libESD");
  // gSystem->Load("libAOD");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libCORRFW");
  // gSystem->Load("libPWGTools");
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLineScalePS(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextTrackDatasetComp(std::string options){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  //  if (options.find("track") == std::string::npos) {
  texcontextDatasetCompAndRadiusAndVarRange = *texDatasetsComparisonCommonDenominator;
  return texcontextDatasetCompAndRadiusAndVarRange;
}

// TString contextPtRange(float* PtRange){
//   std::stringstream ss;
//   ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
//   TString textContext((TString)ss.str());
//   // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
//   return textContext;
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw_Pt_DatasetComparison(std::string options) {
  TH2D* H2D_centrality_track[nDatasets];
  TH3D* H3D_track[nDatasets];

  TH1D* H1D_trackPt[nDatasets];
  TH1D* H1D_trackPt_rebinned[nDatasets];
  
  TH1D* H1D_trackPt_rebinned_ratios[nDatasets];
  TH1D* H1D_spectrum_sum = nullptr;
  TH1D* H1D_spectrum_sum_exclude = nullptr;
  TH1D* H1D_trackPt_all_avg = nullptr;
  TH1D* H1D_trackPt_all_avg_3 = nullptr;



  bool divideSuccess = false;
  std::vector<int> excludedDatasets = {1, 4, 13}; 
  double Nevents; 

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    if (trackHistsObsoleteVersion[iDataset]) {
      H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset], iCollSystem == 0 ? 1 : 0, iCollSystem == 0 ? H2D_centrality_track[iDataset]->GetNbinsX() : -1, "e");
    } else {
      H3D_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt[iDataset] = (TH1D*)H3D_track[iDataset]->ProjectionX("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H3D_track[iDataset]->GetNbinsY(), 1, H3D_track[iDataset]->GetNbinsZ(), "e");
    }

    // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(2.,"trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    cout << "########################################################################################################"<< endl;
    cout << "number of bins x "<< H1D_trackPt[iDataset]->GetNbinsX()<< ",,,,,"<<H1D_trackPt[iDataset]->GetBinCenter(1)<<",,,"<< H1D_trackPt[iDataset]->GetXaxis()->GetXmax()<<endl;
    cout << "########################################################################################################"<< endl;

    int nBinsLog = 20;
    std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(H1D_trackPt[iDataset]->GetBinCenter(1), H1D_trackPt[iDataset]->GetXaxis()->GetXmax(), nBinsLog);
    // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(0.5, 100, nBinsLog);
    double* O2ptLogBins = &O2H1DPtLogBinsVector[0];
    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(nBinsLog, "trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset], O2ptLogBins);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////// add spectrum then divide by nDatasets /////////////////////////////////////////////////////////////////////////////////////////////
    TH1D* rebinnedClone = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_clone_" + Datasets[iDataset] + DatasetsNames[iDataset]);
    double Neventstemp = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
    NormaliseYieldToNEvents(rebinnedClone, Neventstemp);
    if (H1D_spectrum_sum == nullptr) {
        H1D_spectrum_sum = (TH1D*)rebinnedClone->Clone("H1D_spectrum_sum");
    } else {
        H1D_spectrum_sum->Add(rebinnedClone);
    }

    //////////add spectrum except desired ones/////////////////////////////
    if (std::find(excludedDatasets.begin(), excludedDatasets.end(), iDataset) == excludedDatasets.end()) {
      TH1D* rebinnedCloneExclude = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_clone_wExclusion" + Datasets[iDataset] + DatasetsNames[iDataset]);
      double Neventstemp2 = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
      NormaliseYieldToNEvents(rebinnedCloneExclude, Neventstemp2);
      if (H1D_spectrum_sum_exclude == nullptr) {
          H1D_spectrum_sum_exclude = (TH1D*)rebinnedCloneExclude->Clone("H1D_spectrum_sum_exculde");
      } else {
          H1D_spectrum_sum_exclude->Add(rebinnedCloneExclude);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);

    if (options.find("evtNorm") != std::string::npos) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]), analysisWorkflow[iDataset];
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_trackPt_rebinned[iDataset], Nevents);
      cout << "Dataset " << iDataset << ": Nevents = " << Nevents << endl;
    }
    if (options.find("entriesNorm") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackPt_rebinned[iDataset]);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////suite///////////////////////////////////////////////////////
  if (H1D_spectrum_sum && nDatasets > 0) {
    H1D_spectrum_sum->Scale(1.0 / nDatasets);
  }   // H1D_spectrum_sum is the sum of all run's normalized tracks spctra, then divided bu the number of runs

  int nUsedDatasets = nDatasets - excludedDatasets.size();
  if (H1D_spectrum_sum_exclude && nUsedDatasets > 0) {
    H1D_spectrum_sum_exclude->Scale(1.0 / nUsedDatasets);
  }   // H1D_spectrum_sum-exclude is the sum of all runs normalized by the tracks' spectra (excluding problematic ones), and then divided by the number of used runs.
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  TString DatasetsNamesPairRatio[nDatasets];
  int nHistPairRatio = (int)nDatasets / 2;;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      if (iDataset < nHistPairRatio) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
        H1D_trackPt_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPt_rebinned[2*iDataset]->Clone("trackPt_rebinned_ratios"+Datasets[2*iDataset]+DatasetsNames[2*iDataset]);
        H1D_trackPt_rebinned_ratios[iDataset]->Reset("M");
        divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[2*iDataset], H1D_trackPt_rebinned[2*iDataset+1], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      }
    } else {
      H1D_trackPt_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt_rebinned_ratios[iDataset]->Reset("M");
      // divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_trackPt_rebinned[0], 1., 1., datasetsAreSubsetsofId0 ? "b" : ""); // ratio in pt = run/AllRuns
      // divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_spectrum_sum, 1., 1., datasetsAreSubsetsofId0 ? "b" : "");        // ratio in pt = run/AvgOfRuns
      divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_spectrum_sum_exclude, 1., 1., datasetsAreSubsetsofId0 ? "b" : "");   // ratio in pt = run/AvgOfUsedRuns

      
    }
  }
  TString textContext(contextTrackDatasetComp(""));

  TString* textYaxis;
  TString pdfNameNorm;
  if (options.find("evtNorm") != std::string::npos) {
    textYaxis = texTrackPtYield_EntriesNorm;
    pdfNameNorm = (TString)"_EventNorm";
  }
  if (options.find("entriesNorm") != std::string::npos) {
    textYaxis = texTrackPtYield_EntriesNorm;
    pdfNameNorm = (TString)"_EntriesNorm";
  }

  TString* pdfName = new TString("track_Pt_DataComp"+pdfNameNorm);
  TString* pdfName_ratio = new TString("track_Pt_DataComp"+pdfNameNorm+"_ratio");


  Draw_TH1_Histograms(H1D_trackPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, textYaxis, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy"+histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackPt_rebinned_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texPtX, texRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,zoomToOneMedium1");
    } else {
      Draw_TH1_Histograms(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst,logx"+histDatasetComparisonStructure);
    }
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }

  //////////////////////////////////////////////////////////////////////////
  //////////////////suite//////////////////////////////////////////////////
  // TString* pdfNametemp = new TString("track_Pt_after_Sum_all_normspectrum_thenDivide_nDataset");
  // Draw_TH1_Histogram(H1D_spectrum_sum, textContext, pdfNametemp, texPtX, textYaxis, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx, logy");

  // TString* pdfNametemp2 = new TString("ratio_track_Pt_AllRuns_div_averagedSumDataset");
  // TString* y = new TString("ratio Allruns/Avrg");
  bool divideSuccess1 = false;
  H1D_trackPt_all_avg = (TH1D*)H1D_trackPt_rebinned[0]->Clone("trackPt_rebinned_ratios_all_avg");
  H1D_trackPt_all_avg->Reset("M");
  divideSuccess1 = H1D_trackPt_all_avg->Divide(H1D_trackPt_rebinned[0], H1D_spectrum_sum, 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
  // if (divideSuccess == true) {
  //   Draw_TH1_Histogram(H1D_trackPt_all_avg, textContext, pdfNametemp2, texPtX, y, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  // }

  // TString* pdfNametemp3 = new TString("ratio_track_Pt_AllRuns_div_averagedUsedSumDataset");
  // TString* Y = new TString("ratio Allruns/AvrgOfUsed");
  bool divideSuccess3 = false;
  H1D_trackPt_all_avg_3 = (TH1D*)H1D_trackPt_rebinned[0]->Clone("trackPt_rebinned_ratios_all_Usedavg");
  H1D_trackPt_all_avg_3->Reset("M");
  divideSuccess3 = H1D_trackPt_all_avg_3->Divide(H1D_trackPt_rebinned[0], H1D_spectrum_sum_exclude, 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
  // if (divideSuccess3 == true) {
  //   Draw_TH1_Histogram(H1D_trackPt_all_avg_3, textContext, pdfNametemp3, texPtX, Y, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  // }

  TString* pdfNametemp4 = new TString("Ratio");
  TString* Yaxislabel = new TString("ratio Allruns/Avrg");
  const TString Names[2] = {"All Runs/averaged Runs",
                            "All Runs/averaged Used Runs"};
  TH1D** denominatorRatios = new TH1D*[2];
  denominatorRatios[0] = H1D_trackPt_all_avg;
  denominatorRatios[1] = H1D_trackPt_all_avg_3;
  Draw_TH1_Histograms(denominatorRatios, Names, 2, textContext, pdfNametemp4, texPtX, Yaxislabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
}


void Draw_Eta_DatasetComparison(float* ptRange, std::string options) {
  TH2D* H2D_centrality_track[nDatasets];
  TH3D* H3D_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];
  TH1D* H1D_spectrum_sum = nullptr;
  TH1D* H1D_spectrum_sum_exclude = nullptr;
  TH1D* H1D_trackPt_all_avg = nullptr;
  TH1D* H1D_trackPt_all_avg_3 = nullptr;

  bool divideSuccess = false;
  std::vector<int> excludedDatasets = {1, 4, 13}; 
  double Nevents; 

  float ptCutLow = ptRange[0];
  float ptCutHigh = ptRange[1];


  for(int iDataset = 0; iDataset < nDatasets; iDataset++){


    if (trackHistsObsoleteVersion[iDataset]) {
      cout << "WARNING: OBSOLETE track histogram version selected, cannot cut on pT" << endl;
      H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackEta[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset]+DatasetsNames[iDataset],  iCollSystem == 0 ? 1 : 0,  iCollSystem == 0 ? H2D_centrality_track[iDataset]->GetNbinsX() : -1, "e");
    } else {
      H3D_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);

      int ibinPt_low = H3D_track[iDataset]->GetXaxis()->FindBin(ptCutLow);
      int ibinPt_high = H3D_track[iDataset]->GetXaxis()->FindBin(ptCutHigh);
      if (ibinPt_low == 0) 
        cout << "WARNING: Eta_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
      if (ibinPt_high == H3D_track[iDataset]->GetXaxis()->GetNbins()+1) 
        cout << "WARNING: Eta_DatasetComparison is counting the overflow with the chosen PtRange" << endl;

      H1D_trackEta[iDataset] = (TH1D*)H3D_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset]+DatasetsNames[iDataset], ibinPt_low, ibinPt_high, 1, H3D_track[iDataset]->GetNbinsZ(), "e");
    }

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////// add spectrum then divide by nDatasets /////////////////////////////////////////////////////////////////////////////////////////////
    TH1D* rebinnedClone = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_clone_" + Datasets[iDataset] + DatasetsNames[iDataset]);
    double Neventstemp = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
    NormaliseYieldToNEvents(rebinnedClone, Neventstemp);
    if (H1D_spectrum_sum == nullptr) {
        H1D_spectrum_sum = (TH1D*)rebinnedClone->Clone("H1D_spectrum_sum");
    } else {
        H1D_spectrum_sum->Add(rebinnedClone);
    }

    //////////add spectrum except desired ones/////////////////////////////
    if (std::find(excludedDatasets.begin(), excludedDatasets.end(), iDataset) == excludedDatasets.end()) {
      TH1D* rebinnedCloneExclude = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_clone_wExclusion" + Datasets[iDataset] + DatasetsNames[iDataset]);
      double Neventstemp2 = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
      NormaliseYieldToNEvents(rebinnedCloneExclude, Neventstemp2);
      if (H1D_spectrum_sum_exclude == nullptr) {
          H1D_spectrum_sum_exclude = (TH1D*)rebinnedCloneExclude->Clone("H1D_spectrum_sum_exculde");
      } else {
          H1D_spectrum_sum_exclude->Add(rebinnedCloneExclude);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    if (options.find("evtNorm") != std::string::npos) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]), analysisWorkflow[iDataset];
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_trackEta_rebinned[iDataset], Nevents);
    }
    if (options.find("entriesNorm") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackEta_rebinned[iDataset]);
    }
  }


  ////////////////////////////////////////////////////////////////////////
  ////////////suite///////////////////////////////////////////////////////
  if (H1D_spectrum_sum && nDatasets > 0) {
    H1D_spectrum_sum->Scale(1.0 / nDatasets);
  }   // H1D_spectrum_sum is the sum of all run's normalized tracks spctra, then divided bu the number of runs

  int nUsedDatasets = nDatasets - excludedDatasets.size();
  if (H1D_spectrum_sum_exclude && nUsedDatasets > 0) {
    H1D_spectrum_sum_exclude->Scale(1.0 / nUsedDatasets);
  }   // H1D_spectrum_sum-exclude is the sum of all runs normalized by the tracks' spectra (excluding problematic ones), and then divided by the number of used runs.
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  TString DatasetsNamesPairRatio[nDatasets];
  int nHistPairRatio = (int)nDatasets / 2;;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      if (iDataset < nHistPairRatio) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
        H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[2*iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[2*iDataset]+DatasetsNames[2*iDataset]);
        H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
        divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[2*iDataset], H1D_trackEta_rebinned[2*iDataset+1], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      }
    } else {
      H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
      divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      // divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_spectrum_sum, 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      // divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_spectrum_sum_exclude, 1., 1., datasetsAreSubsetsofId0 ? "b" : "");

    }
  }

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextPtRange(ptRange), ""));

  TString* textYaxis;
  TString pdfNameNorm;
  if (options.find("evtNorm") != std::string::npos) {
    textYaxis = texTrackEtaYield_EventNorm;
    pdfNameNorm = (TString)"_EventNorm";
  }
  if (options.find("entriesNorm") != std::string::npos) {
    textYaxis = texTrackEtaYield_EntriesNorm;
    pdfNameNorm = (TString)"_EntriesNorm";
  }

  TString* pdfName = new TString((TString)"track_Eta_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]"+pdfNameNorm);
  TString* pdfName_ratio = new TString((TString)"track_Eta_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]"+pdfNameNorm+"_ratio");

  std::array<std::array<float, 2>, 2> drawnWindowEta = {{{-1, 1}, {260, 390}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnWindowEtaZoom = {{{-1, 1}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
  Draw_TH1_Histograms(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texEtaX, textYaxis, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, ""+histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackEta_rebinned_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texEtaX, texRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, ",zoomToOneExtra");
    } else {
      Draw_TH1_Histograms(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst"+histDatasetComparisonStructure);
    }
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}

void Draw_Phi_DatasetComparison(float* ptRange, std::string options) { 
  TH3D* H3D_track[nDatasets];
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents;

  float ptCutLow = ptRange[0];
  float ptCutHigh = ptRange[1];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    if (trackHistsObsoleteVersion[iDataset]) {
      cout << "WARNING: OBSOLETE track histogram version selected, cannot cut on pT" << endl;
      H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPhi[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPhi_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");
    } else {
      H3D_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);

      int ibinPt_low = H3D_track[iDataset]->GetXaxis()->FindBin(ptCutLow);
      int ibinPt_high = H3D_track[iDataset]->GetXaxis()->FindBin(ptCutHigh);
      if (ibinPt_low == 0) 
        cout << "WARNING: Eta_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
      if (ibinPt_high == H3D_track[iDataset]->GetXaxis()->GetNbins()+1) 
        cout << "WARNING: Eta_DatasetComparison is counting the overflow with the chosen PtRange" << endl;

      H1D_trackPhi[iDataset] = (TH1D*)H3D_track[iDataset]->ProjectionZ("trackPhi_"+Datasets[iDataset]+DatasetsNames[iDataset], ibinPt_low, ibinPt_high, 1, H3D_track[iDataset]->GetNbinsY(), "e");
    }

    H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    cout << "H1D_trackPhi_rebinned[iDataset]->GetEntries() preNorm = " << H1D_trackPhi_rebinned[iDataset]->Integral() << endl;

    if (options.find("evtNorm") != std::string::npos) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]), analysisWorkflow[iDataset];
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
      }
      NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iDataset], Nevents);
      cout << "H1D_trackPhi_rebinned[iDataset]->GetEntries() postNorm = " << H1D_trackPhi_rebinned[iDataset]->Integral() << ", Nevents = " << Nevents << endl;
    }
    if (options.find("entriesNorm") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackPhi_rebinned[iDataset]);
    }
  }


  TString DatasetsNamesPairRatio[nDatasets];
  int nHistPairRatio = (int)nDatasets / 2;;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      if (iDataset < nHistPairRatio) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
        H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[2*iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[2*iDataset]+DatasetsNames[2*iDataset]);
        H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
        divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[2*iDataset], H1D_trackPhi_rebinned[2*iDataset+1], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      }
    } else {
      H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
      divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[iDataset], H1D_trackPhi_rebinned[0], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
    }
  }

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextPtRange(ptRange), ""));

  TString* textYaxis;
  TString pdfNameNorm;
  if (options.find("evtNorm") != std::string::npos) {
    textYaxis = texTrackPhiYield_EventNorm;
    pdfNameNorm = (TString)"_EventNorm";
  }
  if (options.find("entriesNorm") != std::string::npos) {
    textYaxis = texTrackPhiYield_EntriesNorm;
    pdfNameNorm = (TString)"_EntriesNorm";
  }

  TString* pdfName = new TString((TString)"track_Phi_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]"+pdfNameNorm);
  TString* pdfName_ratio = new TString((TString)"track_Phi_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]"+pdfNameNorm+"_ratio");
  
  std::array<std::array<float, 2>, 2> legendPlacementCustom = {{{0.2, 0.2}, {0.4, 0.45}}}; // {{{x1, y1}, {x2, y2}}}

  Draw_TH1_Histograms(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, textYaxis, texCollisionDataInfo, drawnWindowAuto, legendPlacementCustom, contextPlacementAuto, "histWithLine"+histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackPhi_rebinned_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texPhiX, texRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "zoomToOneExtraExtra");
    } else {
      Draw_TH1_Histograms(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst"+histDatasetComparisonStructure);
    }
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}

void Draw_Pt_CentralityComparison(int iDataset) {

  TH3D* H3D_trackPttrackCent;
  TH2D* H2D_trackPttrackCent;
  TH1D* H1D_trackPt[nCentralityBins];
  TH1D* H1D_trackPt_rebinned[nCentralityBins];
  
  if (trackHistsObsoleteVersion[iDataset]) {
    H2D_trackPttrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_CentralityComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_trackPttrackCent->Sumw2();
  } else {
    H3D_trackPttrackCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_eta"))->Clone("Draw_Pt_CentralityComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H3D_trackPttrackCent->Sumw2();
  }

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    if (trackHistsObsoleteVersion[iDataset]) {
      ibinCent_low = H2D_trackPttrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
      ibinCent_high = H2D_trackPttrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
      H1D_trackPt[iCentralityBin] = (TH1D*)H2D_trackPttrackCent->ProjectionY("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    } else {
      ibinCent_low = H3D_trackPttrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
      ibinCent_high = H3D_trackPttrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
      H1D_trackPt[iCentralityBin] = (TH1D*)H3D_trackPttrackCent->ProjectionY("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, 1, H3D_trackPttrackCent->GetNbinsZ(), "e");
    }
      
    H1D_trackPt_rebinned[iCentralityBin] = (TH1D*)H1D_trackPt[iCentralityBin]->Rebin(1.,"trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    NormaliseYieldToNEvents(H1D_trackPt_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_Pt");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_trackPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
  TString* pdfName2 = new TString("track_CentralityComp_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_Pt_logx");
  Draw_TH1_Histograms(H1D_trackPt_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
}



void Draw_Eta_CentralityComparison(int iDataset) {

  TH3D* H3D_trackEtatrackCent;
  TH2D* H2D_trackEtatrackCent;
  TH1D* H1D_trackEta[nCentralityBins];
  TH1D* H1D_trackEta_rebinned[nCentralityBins];
  
  if (trackHistsObsoleteVersion[iDataset]) {
    H2D_trackEtatrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_CentralityComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_trackEtatrackCent->Sumw2();
  } else {
    H3D_trackEtatrackCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_eta"))->Clone("Draw_Eta_CentralityComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H3D_trackEtatrackCent->Sumw2();
  }

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    if (trackHistsObsoleteVersion[iDataset]) {
      cout << "WARNING: OBSOLETE track histogram version selected, cannot cut on pT" << endl;
      ibinCent_low = H2D_trackEtatrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
      ibinCent_high = H2D_trackEtatrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
      H1D_trackEta[iCentralityBin] = (TH1D*)H2D_trackEtatrackCent->ProjectionY("trackEta_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    } else {
      ibinCent_low = H3D_trackEtatrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
      ibinCent_high = H3D_trackEtatrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
      H1D_trackEta[iCentralityBin] = (TH1D*)H3D_trackEtatrackCent->ProjectionZ("trackEta_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, 1, H3D_trackEtatrackCent->GetNbinsY(), "e");
    }
    H1D_trackEta_rebinned[iCentralityBin] = (TH1D*)H1D_trackEta[iCentralityBin]->Rebin(1.,"trackEta_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    // NormaliseYieldToIntegral(H1D_trackEta_rebinned[iCentralityBin]);
    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_Eta");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_trackEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texEtaX, texTrackEtaYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  // Draw_TH1_Histograms(H1D_trackEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texTrackEtaYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");

}

void Draw_Phi_CentralityComparison(int iDataset) {

  TH2D* H2D_trackPhitrackCent;
  TH3D* H3D_trackPhitrackCent;
  TH1D* H1D_trackPhi[nCentralityBins];
  TH1D* H1D_trackPhi_rebinned[nCentralityBins];
  
  if (trackHistsObsoleteVersion[iDataset]) {
    H2D_trackPhitrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_CentralityComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_trackPhitrackCent->Sumw2();
  } else {
    H3D_trackPhitrackCent = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_phi"))->Clone("Draw_Phi_CentralityComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H3D_trackPhitrackCent->Sumw2();
  }

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){

    if (trackHistsObsoleteVersion[iDataset]) {
      cout << "WARNING: OBSOLETE track histogram version selected, cannot cut on pT" << endl;
      ibinCent_low = H2D_trackPhitrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
      ibinCent_high = H2D_trackPhitrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
      H1D_trackPhi[iCentralityBin] = (TH1D*)H2D_trackPhitrackCent->ProjectionY("trackPhi_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
    } else {
      ibinCent_low = H3D_trackPhitrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
      ibinCent_high = H3D_trackPhitrackCent->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
      H1D_trackPhi[iCentralityBin] = (TH1D*)H3D_trackPhitrackCent->ProjectionZ("trackPhi_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, 1, H3D_trackPhitrackCent->GetNbinsY(), "e");
    }
    H1D_trackPhi_rebinned[iCentralityBin] = (TH1D*)H1D_trackPhi[iCentralityBin]->Rebin(1.,"trackPhi_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_Phi");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_trackPhi_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPhiX, texTrackPhiYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}



void Draw_Pt_Run2Run3Comparison_0010Cent(int iDataset) {

  int ibinCent_low, ibinCent_high;

  cout << "o2physics" << endl;
  // O2Physics Run 3

  TH2D* H3D_trackPttrackCent;
  TH1D* H1D_trackPt;
  TH1D* H1D_trackPt_rebinned;

  H3D_trackPttrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_CentralityComparison"+DatasetsNames[iDataset]);
  H3D_trackPttrackCent->Sumw2();

  ibinCent_low = H3D_trackPttrackCent->GetXaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_trackPttrackCent->GetXaxis()->FindBin((double)10 - GLOBAL_epsilon);
  H1D_trackPt = (TH1D*)H3D_trackPttrackCent->ProjectionY("trackPt_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
  H1D_trackPt_rebinned = (TH1D*)H1D_trackPt->Rebin(1.,"trackPt_rebinned_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  NormaliseYieldToNEvents(H1D_trackPt_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], 0, 10));

  cout << "aliPhysics" << endl;
  // AliPhysics Run 2

  TH1D* H1D_run2_trackPt;
  TH1D* H1D_run2_trackPt_rebinned;
  H1D_run2_trackPt = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("tracks"))->FindObject("histTrackPt_0");
  std::vector<double> O2H1DPtbinsVector = GetTH1Bins(H1D_trackPt_rebinned);
  double* O2ptBins = &O2H1DPtbinsVector[0];
  H1D_run2_trackPt_rebinned = (TH1D*)H1D_run2_trackPt->Rebin(H1D_trackPt_rebinned->GetNbinsX(), "H1D_run2_trackPt_rebinned", O2ptBins);

  double nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_trackPt_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  cout << "ratio" << endl;
  // ratio

  TH1D* H1D_trackVar_rebinned_ratio;
  bool divideSuccess = false;

  H1D_trackVar_rebinned_ratio = (TH1D*)H1D_trackPt_rebinned->Clone("trackVar_rebinned_ratios"+DatasetsNames[iDataset]);
  H1D_trackVar_rebinned_ratio->Reset("M");
  divideSuccess = H1D_trackVar_rebinned_ratio->Divide(H1D_trackPt_rebinned, H1D_run2_trackPt_rebinned);

  cout << "plotting" << endl;
  // plotting
  TH1D* hist_list[2] = {H1D_trackPt_rebinned, H1D_run2_trackPt_rebinned};
  TString RunCompLegend[2];
  std::stringstream ss;

  ss << "run 3";
  RunCompLegend[0] = (TString)ss.str();
  ss.str("");
  ss.clear();
  ss << "run 2";
  RunCompLegend[1] = (TString)ss.str();
  ss.str("");
  ss.clear();

  TString* pdfName = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Pt");
  TString* pdfName_ratio = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Pt_ratio");

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms(hist_list, RunCompLegend, 2, textContext, pdfName, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
  TString* pdfName2 = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Pt_logx");
  Draw_TH1_Histograms(hist_list, RunCompLegend, 2, textContext, pdfName2, texPtX, texTrackPtYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
  if (divideSuccess == true) {
    Draw_TH1_Histogram(H1D_trackVar_rebinned_ratio, textContext, pdfName_ratio, texPtX, texRatioRun3Run2, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, logx, ratioLine");
  }
  else {
    cout << "Divide failed in Draw_Pt_Run2Run3Comparison_0010Cent" << endl;
  }
}

void Draw_Eta_Run2Run3Comparison_0010Cent(int iDataset) {

  int ibinCent_low, ibinCent_high;


  // O2Physics Run 3

  TH2D* H3D_trackEtatrackCent;
  TH1D* H1D_trackEta;
  TH1D* H1D_trackEta_rebinned;

  H3D_trackEtatrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_CentralityComparison"+DatasetsNames[iDataset]);
  H3D_trackEtatrackCent->Sumw2();

  ibinCent_low = H3D_trackEtatrackCent->GetXaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_trackEtatrackCent->GetXaxis()->FindBin((double)10 - GLOBAL_epsilon);
  H1D_trackEta = (TH1D*)H3D_trackEtatrackCent->ProjectionY("trackEta_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
  H1D_trackEta_rebinned = (TH1D*)H1D_trackEta->Rebin(1.,"trackEta_rebinned_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  NormaliseYieldToNEvents(H1D_trackEta_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], 0, 10));

  // AliPhysics Run 2

  TH1D* H1D_run2_trackEta;
  TH1D* H1D_run2_trackEta_rebinned;
  H1D_run2_trackEta = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("tracks"))->FindObject("histTrackEta_0");
  std::vector<double> O2H1DEtabinsVector = GetTH1Bins(H1D_trackEta_rebinned);
  double* O2etaBins = &O2H1DEtabinsVector[0];
  H1D_run2_trackEta_rebinned = (TH1D*)H1D_run2_trackEta->Rebin(H1D_trackEta_rebinned->GetNbinsX(), "H1D_run2_trackEta_rebinned", O2etaBins);

  double nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_trackEta_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  // plotting
  TH1D* hist_list[2] = {H1D_trackEta_rebinned, H1D_run2_trackEta_rebinned};
  TString RunCompLegend[2];
  std::stringstream ss;

  ss << "run 3";
  RunCompLegend[0] = (TString)ss.str();
  ss.str("");
  ss.clear();
  ss << "run 2";
  RunCompLegend[1] = (TString)ss.str();
  ss.str("");
  ss.clear();



  TString* pdfName = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Eta");

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms(hist_list, RunCompLegend, 2, textContext, pdfName, texEtaX, texTrackEtaYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}

void Draw_Phi_Run2Run3Comparison_0010Cent(int iDataset) {

  int ibinCent_low, ibinCent_high;


  // O2Physics Run 3

  TH2D* H3D_trackPhitrackCent;
  TH1D* H1D_trackPhi;
  TH1D* H1D_trackPhi_rebinned;

  H3D_trackPhitrackCent = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_phi"))->Clone("Draw_Phi_CentralityComparison"+DatasetsNames[iDataset]);
  H3D_trackPhitrackCent->Sumw2();

  ibinCent_low = H3D_trackPhitrackCent->GetXaxis()->FindBin((double)0 + GLOBAL_epsilon);
  ibinCent_high = H3D_trackPhitrackCent->GetXaxis()->FindBin((double)10 - GLOBAL_epsilon);
  H1D_trackPhi = (TH1D*)H3D_trackPhitrackCent->ProjectionY("trackPhi_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");
  H1D_trackPhi_rebinned = (TH1D*)H1D_trackPhi->Rebin(1.,"trackPhi_rebinned_"+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

  NormaliseYieldToNEvents(H1D_trackPhi_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], 0, 10));

  // AliPhysics Run 2

  TH1D* H1D_run2_trackPhi;
  TH1D* H1D_run2_trackPhi_rebinned;
  H1D_run2_trackPhi = (TH1D*)((TObject*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("tracks"))->FindObject("histTrackPhi_0");
  std::vector<double> O2H1DPhibinsVector = GetTH1Bins(H1D_trackPhi_rebinned);
  double* O2phiBins = &O2H1DPhibinsVector[0];
  H1D_run2_trackPhi_rebinned = (TH1D*)H1D_run2_trackPhi->Rebin(H1D_trackPhi_rebinned->GetNbinsX(), "H1D_run2_trackPhi_rebinned", O2phiBins);

  double nEvents_CentWindow;
  TH1D* H1D_Centrality_Run2= (TH1D*)((TObject*)file_AliAnalysis->Get("JetCore_JetCore_57_2050_04_histos"))->FindObject("fHistCentrality");
  int iBinCent_low = H1D_Centrality_Run2->GetXaxis()->FindBin(0 + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_Run2->GetXaxis()->FindBin(10 - GLOBAL_epsilon);
  nEvents_CentWindow = H1D_Centrality_Run2->Integral(iBinCent_low, iBinCent_high);

  H1D_run2_trackPhi_rebinned->Scale(1./nEvents_CentWindow,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

  // plotting
  TH1D* hist_list[2] = {H1D_trackPhi_rebinned, H1D_run2_trackPhi_rebinned};
  TString RunCompLegend[2];
  std::stringstream ss;

  ss << "run 3";
  RunCompLegend[0] = (TString)ss.str();
  ss.str("");
  ss.clear();
  ss << "run 2";
  RunCompLegend[1] = (TString)ss.str();
  ss.str("");
  ss.clear();



  TString* pdfName = new TString("track_Run2_vs_Run3_"+DatasetsNames[iDataset]+"_@cent[00,10]_Phi");

  TString textContext(contextCustomOneField(RunCompLegend[0]+" vs "+RunCompLegend[1], ""));

  Draw_TH1_Histograms(hist_list, RunCompLegend, 2, textContext, pdfName, texPhiX, texTrackPhiYield_EventNorm_CentWindow, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_Sigmapt_vs_pt_DatasetComp() {
  TH3D* HDsparse_cent_sigmapt_pt[nDatasets];
  TH3D* H3D_cent_sigmapt_pt[nDatasets];

  TH2D* H2D_sigmapt_pt_concatenated[nDatasets];
  TH2D* H2D_sigmapt_pt[nDatasets];
  TH2D* H2D_sigmapt_pt_high[nDatasets];
  TH2D* H2D_sigmapt_pt_high_rebinnedX[nDatasets];

  TH1D* H1D_sigmapt_pt_X_forMedian[nDatasets];

  TH1D* H1D_sigmapt_pt_mean_withProfile[nDatasets];
  TH1D* H1D_sigmapt_pt_high_mean_withProfile[nDatasets];
  TH1D* H1D_sigmapt_pt_median[nDatasets];

  TH1D* H1D_sigmapt_pt_mean_withProfile_concatenated[nDatasets];

  double rebin_xMin = 0;
  double rebin_xMiddle = 10;
  double rebin_xMax = 100;
  int rebin_nLeft = 50;
  int rebin_nRight = 10;
  int nbinsX;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_sigmapt_pt[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_track_pt_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_DatasetComp_left"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_sigmapt_pt_high[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_track_pt_high_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_high_DatasetComp_right"+Datasets[iDataset]+DatasetsNames[iDataset]);


    H2D_sigmapt_pt_high_rebinnedX[iDataset] = (TH2D*)H2D_sigmapt_pt_high[iDataset]->RebinX(1.,"H2D_sigmapt_pt_high_rebinnedX"+Datasets[iDataset]+DatasetsNames[iDataset]);

    // H2D_sigmapt_pt[iDataset] = (TH2D*)HDsparse_cent_sigmapt_pt[iDataset]->Project3D(dummyName[0]+Form("%d", iDataset)+"_sigmapt_e_zy"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    // H2D_sigmapt_pt[iDataset]->Sumw2();

    // H2D_sigmapt_pt_rebinned_simple[iDataset] = (TH2D*)H2D_sigmapt_pt[iDataset]->RebinX(1.,"H2D_sigmapt_pt_rebinnedSimple"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_sigmapt_pt_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_sigmapt_pt_rebinned[iDataset]->GetNbinsX()-10);
    // H2D_rhoCentrality_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_sigmapt_pt_rebinned[iDataset]->FindLastBinAbove(1, 2)); //(asks for the last bin on the y axis (axis number 2) to have strictly more than 1 entry)

    // H1D_sigmapt_pt_median[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionX("H1D_sigmapt_pt_median"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    // H1D_sigmapt_pt_median[iDataset]->Reset("M");
    // cout << "testC" << endl;
 
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins((TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionX("H1D_sigmapt_pt_left"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e"));
    std::vector<double> xbinsVectorRight = GetTH1Bins((TH1D*)H2D_sigmapt_pt_high_rebinnedX[iDataset]->ProjectionX("H1D_sigmapt_pt_right"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e"));
    xbinsVectorRight.erase(xbinsVectorRight.begin());
    cout << "xbinsVectorLeft.front() = " << xbinsVectorLeft.front() << ", xbinsVectorLeft.back() = " << xbinsVectorLeft.back() << ", xbinsVectorRight.front() = " << xbinsVectorRight.front() << ", xbinsVectorRight.back() = " << xbinsVectorRight.back() << endl;
    std::vector<double> xbinsVectorCombination = xbinsVectorLeft;
    xbinsVectorCombination.insert( xbinsVectorCombination.end(), xbinsVectorRight.begin(), xbinsVectorRight.end() );
    double* xbins_new = &xbinsVectorCombination[0];
    cout << "xbinsVectorCombination.size() = " << xbinsVectorCombination.size() << endl;



    H1D_sigmapt_pt_mean_withProfile[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProfileX("H1D_sigmapt_pt_rebinned_mean_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    H1D_sigmapt_pt_high_mean_withProfile[iDataset] = (TH1D*)H2D_sigmapt_pt_high_rebinnedX[iDataset]->ProfileX("H1D_sigmapt_pt_high_rebinned_mean_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");

    TH1D H1D_sigmapt_pt_mean_withProfile_concatenated_temp("H1D_sigmapt_pt_mean_withProfile_concatenated_temp"+Datasets[iDataset]+DatasetsNames[iDataset], "H1D_sigmapt_pt_mean_withProfile_concatenated_temp"+Datasets[iDataset]+DatasetsNames[iDataset], xbinsVectorCombination.size()-1, xbins_new);
    H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset] = (TH1D*)H1D_sigmapt_pt_mean_withProfile_concatenated_temp.Clone("H1D_sigmapt_pt_mean_withProfile_concatenated"+Datasets[iDataset]+DatasetsNames[iDataset]);

    // y-axis
    std::vector<double> ybinsVector = GetTH1Bins((TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionY("H1D_sigmapt_pt_y"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e"));
    double* ybins_new = &ybinsVector[0];

    TH2D H2D_sigmapt_pt_concatenated_temp("H2D_sigmapt_pt_concatenated_temp"+Datasets[iDataset]+DatasetsNames[iDataset], "H2D_sigmapt_pt_concatenated_temp"+Datasets[iDataset]+DatasetsNames[iDataset], xbinsVectorCombination.size()-1, xbins_new, ybinsVector.size()-1, ybins_new);
    H2D_sigmapt_pt_concatenated[iDataset] = (TH2D*)H2D_sigmapt_pt_concatenated_temp.Clone("H2D_sigmapt_pt_concatenated"+Datasets[iDataset]+DatasetsNames[iDataset]);


    double binContent, binError;
    int iBinLowY_OldHist, iBinHighY_OldHist;
    for(int iBinX = 1; iBinX <= H2D_sigmapt_pt[iDataset]->GetNbinsX(); iBinX++){
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinContent(iBinX, H1D_sigmapt_pt_mean_withProfile[iDataset]->GetBinContent(iBinX));
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinError(iBinX, H1D_sigmapt_pt_mean_withProfile[iDataset]->GetBinError(iBinX));
      for(int iBinY = 1; iBinY <= H2D_sigmapt_pt[iDataset]->GetNbinsY(); iBinY++){
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinContent(iBinX, iBinY));
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinError(iBinX, iBinY));
        if (iDataset ==15){
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinContent(iBinX, iBinY) - H2D_sigmapt_pt[0]->GetBinContent(iBinX, iBinY));
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(iBinX, iBinY, sqrt(H2D_sigmapt_pt_concatenated[iDataset]->GetBinContent(iBinX, iBinY)));
        } 
      }
    }

    for(int iBinX = 1; iBinX <= H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetNbinsX(); iBinX++){
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, H1D_sigmapt_pt_high_mean_withProfile[iDataset]->GetBinContent(iBinX));
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, H1D_sigmapt_pt_high_mean_withProfile[iDataset]->GetBinError(iBinX));
      for(int iBinY = 1; iBinY <= H2D_sigmapt_pt[iDataset]->GetNbinsY(); iBinY++){
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetBinContent(iBinX, iBinY));
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetBinError(iBinX, iBinY));
        if (iDataset ==15){
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetBinContent(iBinX, iBinY) - H2D_sigmapt_pt_high_rebinnedX[0]->GetBinContent(iBinX, iBinY));
          H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, sqrt(H2D_sigmapt_pt_concatenated[iDataset]->GetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY)));
        } 
      }
    }

    H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset] = (TH1D*)H2D_sigmapt_pt_concatenated[iDataset]->ProfileX("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");



    // // median instead of mean:
    // double x, q;
    // q = 0.5;
    // for(int iBin = 1; iBin <= H1D_sigmapt_pt_median[iDataset]->GetNbinsX(); iBin++){
    //   // H2D_sigmapt_pt[iDataset]->GetXaxis()->SetRange(iBin,iBin);
    //   H1D_sigmapt_pt_X_forMedian[iDataset] = (TH1D*)H2D_sigmapt_pt[iDataset]->ProjectionY("H2D_sigmapt_pt_X"+Datasets[iDataset]+DatasetsNames[iDataset], iBin, iBin, "e");

    //   H1D_sigmapt_pt_X_forMedian[iDataset]->GetQuantiles(1, &x, &q);
    //   H1D_sigmapt_pt_median[iDataset]->SetBinContent(iBin, x);
    //   H1D_sigmapt_pt_median[iDataset]->SetBinError(iBin, 0.0001); // no idea how to get the error on the median calculation
    // }

    H2D_sigmapt_pt_concatenated[iDataset]->Scale(1./H2D_sigmapt_pt_concatenated[iDataset]->Integral(1, H2D_sigmapt_pt_concatenated[iDataset]->GetNbinsX(), 1, H2D_sigmapt_pt_concatenated[iDataset]->GetNbinsY(), "width"));
  }

  TString* pdfName = new TString("track_sigmapt_vs_pt_DataComp");
  TString* pdfName_logy = new TString("track_sigmapt_vs_pt_DataComp_logy");

  TString* pdfName_mean_withProfile = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile");
  TString* pdfName_mean_withProfile_logy = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile_logy");

  // TString* pdfName_median = new TString("track_sigmapt_median_vs_pt_DataComp");
  // TString* pdfName_median_logy = new TString("track_sigmapt_median_vs_pt_DataComp_logy");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, *texDatasetsComparisonType, ""));

  std::array<std::array<float, 3>, 3> drawnWindow2DSigma = {{{0.1, 100}, {0.001, 100}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
  // Draw_TH2_Histograms(H2D_sigmapt_pt_concatenated, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowSigma, th2ContoursNone, contourNumberNone, "logx,logz,autoRangeSame"); // ?
  // Draw_TH2_Histograms(H2D_sigmapt_pt_concatenated, DatasetsNames, nDatasets, textContext, pdfName_logy, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindow2DSigma, th2ContoursNone, contourNumberNone, "logx,logy,logz,autoRangeSame"); // ?

  // Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx");
  std::array<std::array<float, 2>, 2> drawnWindowSigmaAverage = {{{0.1, 100}, {0.001, 10}}}; // {{xmin, xmax}, {ymin, ymax}}
  Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile_logy, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowSigmaAverage, legendPlacementAuto, contextPlacementAuto, "logx,logy");

  // Draw_TH1_Histograms(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx");
  // Draw_TH1_Histograms(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median_logy, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
}


void Draw_SelectedMultiplicity_DatasetComp() {

  TH2D* H2D_rhoMult[nDatasets];
  TH1D* H1D_mult[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoMult[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_ntracks"))->Clone("Draw_NtracksPerEvent_DatasetComp"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_rhoMult[iDataset]->Sumw2();

    // H2D_rhoMult_rebinned[iDataset] = (TH2D*)H2D_rhoMult[iDataset]->RebinX(1.,"H2D_rhoMult_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_rhoMult_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
    // H2D_rhoMult_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 2));
    // H2D_rhoMult_rebinned[iDataset]->Scale(1./H2D_rhoMult_rebinned[iDataset]->Integral(1, H2D_rhoMult_rebinned[iDataset]->GetNbinsX(), 1, H2D_rhoMult_rebinned[iDataset]->GetNbinsY(), "width"));

    H1D_mult[iDataset] = H2D_rhoMult[iDataset]->ProjectionY("H1D_mult"+(TString)Datasets[iDataset], 1, H2D_rhoMult[iDataset]->GetNbinsX());
  }

  TString* pdfName = new TString("track_DataComp_multiplicitySelected");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_mult, DatasetsNames, nDatasets, textContext, pdfName, texSelectedMultiplicity, texCount, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
}


void Draw_SelectedMultiplicity_CentralityComp(int iDataset, std::string options) {

  TH2D* H2D_multCentrality;

  TH1D* H1D_mult[nCentralityBins];
  TH1D* H1D_mult_rebinned[nCentralityBins];

  H2D_multCentrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_ntracks"))->Clone("Draw_SelectedMultiplicity_CentralityProjection"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H2D_multCentrality->Sumw2();

  // float FluctuLow = -60;
  // float FluctuHigh = 60;
  // int ibinFluctu_low  = H2D_multCentrality->GetYaxis()->FindBin(FluctuLow);
  // int ibinFluctu_high = H2D_multCentrality->GetYaxis()->FindBin(FluctuHigh);
  // H2D_multCentrality->GetYaxis()->SetRange(ibinFluctu_low, ibinFluctu_high);

  int ibinCent_low, ibinCent_high;
  TString CentralityLegend[nCentralityBins];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ibinCent_low = H2D_multCentrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin]);
    ibinCent_high = H2D_multCentrality->GetXaxis()->FindBin(arrayCentralityBinning[iCentralityBin+1])-1;
    // cout << "ibinCent_low = " << ibinCent_low << ", ibinCent_high = " << ibinCent_high << endl;
    H1D_mult[iCentralityBin] = (TH1D*)H2D_multCentrality->ProjectionY("multCentrality_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ibinCent_low, ibinCent_high, "e");

    H1D_mult_rebinned[iCentralityBin] = (TH1D*)H1D_mult[iCentralityBin]->Rebin(1.,"multCentrality_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]");

    yAxisLabel = texNoNorm_SelMultYield;
    if (options.find("normEntries") != std::string::npos) {
      NormaliseYieldToNEntries(H1D_mult_rebinned[iCentralityBin]);
      yAxisLabel = texEntriesNorm_selMultYield;
    }
    if (options.find("normEvents") != std::string::npos) {
      NormaliseYieldToNEvents(H1D_mult_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]));
      yAxisLabel = texCollNorm_selMultYield;
    }
    if (options.find("normEventsCentrality") != std::string::npos) {
      NormaliseYieldToNEvents(H1D_mult_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));
      yAxisLabel = texCollNorm_selMultYield_CentWindow;
    }
    
    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+DatasetsNames[iDataset]+"_selMult_logy");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_mult_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texSelectedMultiplicity, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoXrange,logy");
  TString* pdfName2 = new TString("track_CentralityComp_"+DatasetsNames[iDataset]+"_selMult");
  Draw_TH1_Histograms(H1D_mult_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texSelectedMultiplicity, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoXrange");
}


void Draw_Mean_Pt_vs_Dataset() {
  TH2D* H2D_centrality_track[nDatasets];
  TH3D* H3D_centrality_track[nDatasets];

  TH1D* H1D_trackPt[nDatasets];
  
  TH1D* H1D_mean_vs_dataset = new TH1D("H1D_mean_vs_dataset", "H1D_mean_vs_dataset", nDatasets, 0, (float)nDatasets);

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    if (trackHistsObsoleteVersion[iDataset]) {
      H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");
    } else {
      H3D_centrality_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_eta"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt[iDataset] = (TH1D*)H3D_centrality_track[iDataset]->ProjectionY("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H3D_centrality_track[iDataset]->GetNbinsX(), 1, H3D_centrality_track[iDataset]->GetNbinsZ(), "e");
    }
    H1D_trackPt[iDataset]->ResetStats();

    H1D_mean_vs_dataset->SetBinContent(iDataset+1, H1D_trackPt[iDataset]->GetMean());
    H1D_mean_vs_dataset->SetBinError(iDataset+1, H1D_trackPt[iDataset]->GetMean(11));
    H1D_mean_vs_dataset->GetXaxis()->SetBinLabel(iDataset+1, DatasetsNames[iDataset]);
  }

  TString* pdfName = new TString("track_mean_Pt_vs_dataset");

  TString textContext(contextTrackDatasetComp(""));
  TString* emptyAxisName = new TString("");

  Draw_TH1_Histogram(H1D_mean_vs_dataset, textContext, pdfName, emptyAxisName, texMeanPt, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "datasetXaxisBinLabels");
  cout << " has memory leaks somewhere" << endl;
}

void Draw_Mean_Eta_vs_Dataset() {
  TH3D* H3D_centrality_track[nDatasets];
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  
  TH1D* H1D_mean_vs_dataset = new TH1D("H1D_mean_vs_dataset", "H1D_mean_vs_dataset", nDatasets, 0, (float)nDatasets);

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    if (trackHistsObsoleteVersion[iDataset]) {
      H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackEta[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");
    } else {
      H3D_centrality_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_centrality_track_pt_track_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackEta[iDataset] = (TH1D*)H3D_centrality_track[iDataset]->ProjectionZ("trackEta_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H3D_centrality_track[iDataset]->GetNbinsX(), 1, H3D_centrality_track[iDataset]->GetNbinsY(), "e");
    }
    H1D_trackEta[iDataset]->ResetStats();

    H1D_mean_vs_dataset->SetBinContent(iDataset+1, H1D_trackEta[iDataset]->GetMean());
    H1D_mean_vs_dataset->SetBinError(iDataset+1, H1D_trackEta[iDataset]->GetMean(11));
    H1D_mean_vs_dataset->GetXaxis()->SetBinLabel(iDataset+1, DatasetsNames[iDataset]);
  }

  TString* pdfName = new TString("track_mean_Eta_vs_dataset");

  TString textContext(contextTrackDatasetComp(""));
  TString* emptyAxisName = new TString("");

  Draw_TH1_Histogram(H1D_mean_vs_dataset, textContext, pdfName, emptyAxisName, texMeanEta, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "datasetXaxisBinLabels");
}

void Draw_Mean_Ntrack_vs_Dataset() {
  TH2D* H2D_centrality_ntrack[nDatasets];

  TH1D* H1D_ntrack[nDatasets];
  
  TH1D* H1D_mean_vs_dataset = new TH1D("H1D_mean_vs_dataset", "H1D_mean_vs_dataset", nDatasets, 0, (float)nDatasets);

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centrality_ntrack[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_ntracks"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H1D_ntrack[iDataset] = (TH1D*)H2D_centrality_ntrack[iDataset]->ProjectionY("ntracks_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_centrality_ntrack[iDataset]->GetNbinsX(), "e");
    H1D_ntrack[iDataset]->ResetStats();

    H1D_mean_vs_dataset->SetBinContent(iDataset+1, H1D_ntrack[iDataset]->GetMean());
    H1D_mean_vs_dataset->SetBinError(iDataset+1, H1D_ntrack[iDataset]->GetMean(11));
    H1D_mean_vs_dataset->GetXaxis()->SetBinLabel(iDataset+1, DatasetsNames[iDataset]);
  }

  TString* pdfName = new TString("track_mean_Ntracks_vs_dataset");

  TString textContext(contextTrackDatasetComp(""));
  TString* emptyAxisName = new TString("");

  Draw_TH1_Histogram(H1D_mean_vs_dataset, textContext, pdfName, emptyAxisName, texMeanNtracks, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "datasetXaxisBinLabels");
}

void Draw_DcaXY_DatasetComp();



void Draw_tracks_2D_profile(std::string options, int refIndex = 0){
  //only used when "trackHistsObsoleteVersion" is false in the settings, and the analysisWorkflow name is "jet-finder-charged-qa" because the 3D histogram "h3_track_pt_track_eta_track_phi" exists there
  TH2D* H2D_etaPhi[nDatasets];
  TH3D* H3D_track[nDatasets];

  // Create a TCanvas to hold the plots
  TString* yAxisLabel;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    if (analysisWorkflow[iDataset].Contains("jet-finder-charged-qa") == true){
      H3D_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi"))->Clone("Draw_2D_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    }

    TString projectionAxis = "yz";  // Default is eta-phi
    if (options.find("pteta") != std::string::npos) {
      projectionAxis = "xy";
    } else if (options.find("ptphi") != std::string::npos){
      projectionAxis = "xz";
    } 

    H2D_etaPhi[iDataset] = dynamic_cast<TH2D*>(H3D_track[iDataset]->Project3D(projectionAxis));
    H2D_etaPhi[iDataset]->SetName(Form("TrackProjection_%s_%d", projectionAxis.Data(), iDataset));
    H2D_etaPhi[iDataset]->SetTitle(Form("Tracks %s Projection - Run %s", projectionAxis.Data(), DatasetsNames[iDataset].Data()));
    
    double Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset]);
    H2D_etaPhi[iDataset]->Scale(1.0 / Nevents);
  }

  //Draw in one canvas 

  TCanvas* c = new TCanvas("c", "H2D_TracketaPhi collection", 1200, 800);

  // Compute rows and columns for dividing canvas
  int nCols = ceil(sqrt(nDatasets));
  int nRows = ceil((double)nDatasets / nCols);
  c->Divide(nCols, nRows);

  // Loop over the histograms and draw them
  for (int i = 0; i < nDatasets; ++i) {
      c->cd(i+1); // Canvas pads are 1-indexed
      H2D_etaPhi[i]->Draw("COLZ");
        // Add title to the top-left corner of each pad
      TLatex latex;
      latex.SetNDC(); // Use Normalized Device Coordinates
      latex.SetTextSize(0.04); // Adjust size as needed
      latex.DrawLatex(0.1, 0.92, DatasetsNames[i]); // x, y, and text
  }

  c->Update();  
  //Plot_2D_Ratio(H2D_etaPhi, "etaphi", nDatasets, refIndex);


}

