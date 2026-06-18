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
#include "TProfile.h"
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
#include <iostream>
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

void Draw_Eta_DatasetComparison_trackSelComp();
void Draw_Phi_DatasetComparison_trackSelComp();


void Draw_Pt_CentralityComparison(int iDataset);
void Draw_Eta_CentralityComparison(int iDataset, std::string options);
void Draw_Phi_CentralityComparison(int iDataset, std::string options);

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
void Draw_SelectedMultiplicity_CentralityComp(int iDataset, std::string options);

void Draw_Mean_Ntrack_vs_Dataset();
void Draw_Mean_Pt_vs_Dataset();
void Draw_Mean_Eta_vs_Dataset();
void Draw_DcaXY_DatasetComp();

void Draw_TrackSelDistribs_DatasetComparison(float* ptRange);
void Draw_TrackSelDistribs_PtComparison(int iDataset, float* ptBins, int nBinsPt);

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
    // Draw_Eta_CentralityComparison(iDataset, "entriesNorm");
    // Draw_Phi_CentralityComparison(iDataset, "entriesNorm");
  }
  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_Pt_Run2Run3Comparison_0010Cent(iDataset);
  //   Draw_Eta_Run2Run3Comparison_0010Cent(iDataset);
  //   Draw_Phi_Run2Run3Comparison_0010Cent(iDataset);
  // }

  // const int nPtBins = 6;
  const int nPtBins = 1;
  float jetPtMinCut, jetPtMaxCut;
  // float jetPtMinCutArray[nPtBins+1] = {0.15, 1, 2, 4, 6, 20, 100};
  float jetPtMinCutArray[nPtBins+1] = {0, 200};


  // Draw_Pt_DatasetComparison("evtNorm");
  // Draw_Pt_DatasetComparison("entriesNorm");
  Draw_Pt_DatasetComparison("20GevNorm");
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

    // Draw_TrackSelDistribs_DatasetComparison(ptRange);

  }
  int iDataset = 0;

  // Draw_Sigmapt_vs_pt_DatasetComp();
  // Draw_Sigmapt_nonGlobal_uniformTracks();
  // Draw_Sigmapt_nonGlobal_uniformTracks_fromSubtraction();
  // Draw_Sigmapt_nonGlobal_uniformTracks_centralEta();

  // Draw_SelectedMultiplicity_DatasetComp();
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    // Draw_SelectedMultiplicity_CentralityComp(iDataset, "normEvents");
  }

  // Draw_Mean_Pt_vs_Dataset();
  // Draw_Mean_Eta_vs_Dataset();
  // Draw_Mean_Ntrack_vs_Dataset();
  // Draw_DcaXY_DatasetComp();

  // Draw_TrackSelDistribs_PtComparison(iDataset, jetPtMinCutArray, nPtBins);

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

  bool divideSuccess = false;

  double Nevents; 

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    if (trackHistsObsoleteVersion[iDataset]) {
      H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_track_pt"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset], iCollSystem == 0 ? 1 : 0, iCollSystem == 0 ? H2D_centrality_track[iDataset]->GetNbinsX() : -1, "e");
    } else {
      H3D_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackPt[iDataset] = (TH1D*)H3D_track[iDataset]->ProjectionX("trackPt_"+Datasets[iDataset]+DatasetsNames[iDataset], 1, H3D_track[iDataset]->GetNbinsY(), 1, H3D_track[iDataset]->GetNbinsZ(), "e");
    }

    // ===================== Basic rebinning =====================
    // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(2.,"trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    
    // int nBinsLogRough = 100;
    // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(H1D_trackPt[iDataset], nBinsLogRough);
    // int nBinsLogResult = O2H1DPtLogBinsVector.size()-1;
    // // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(0.5, 100, nBinsLogRough);
    // double* O2ptLogBins = &O2H1DPtLogBinsVector[0];
    // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(nBinsLogResult, "trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset], O2ptLogBins);

    // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);

    // ===================== Advanced rebinning =====================
    std::vector<double> xbinsVector = GetTH1Bins(H1D_trackPt[iDataset]);
    double* xbinsInitial = &xbinsVector[0];
    double ptBinsNew[500]; // 500 to be safe
    int iBinNew = 0;
    int iBinInitial = 0;
    double binWidth = H1D_trackPt[iDataset]->GetXaxis()->GetBinWidth(1); //to be use in case we want to have larger bins in mid-pt range

    int increment;

    while (iBinInitial < H1D_trackPt[iDataset]->GetNbinsX()) {
        double pt = xbinsInitial[iBinInitial];

        if (pt < 0.3) {
            increment = 1;
        }
        else if (pt < 25.0) {
            increment = 1; 
            // increment = 1 + pt / binWidth / 8; //in case we want to have larger bins in mid-pt range
        }
        else {
            //no need of increment because we want to stop at 100
            ptBinsNew[iBinNew++] = pt;   // left edge of the large bin
            ptBinsNew[iBinNew++] = 100.0; // right edge of the large bin
            break;
      }
        ptBinsNew[iBinNew++] = pt;
        iBinInitial += increment;
    }

    // security to ensure the last bin edge is 100 GeV
    if (ptBinsNew[iBinNew-1] < 100.0) ptBinsNew[iBinNew++] = 100.0;

    int nBinsNew = iBinNew - 1;

    // ======= Creation of the rebinned histogram =======
    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(nBinsNew,"H1D_trackPt_rebinned" + Datasets[iDataset] + DatasetsNames[iDataset], ptBinsNew);

    // ===== Display of the initial histogram bins (for debbug) =====
    std::cout << "=== Bins de l'histogramme initial ===" << std::endl;
    for (int i = 1; i <= H1D_trackPt[iDataset]->GetNbinsX(); i++) {
        double xlow  = H1D_trackPt[iDataset]->GetBinLowEdge(i);
        double xhigh = H1D_trackPt[iDataset]->GetBinLowEdge(i) + H1D_trackPt[iDataset]->GetBinWidth(i);
        double content = H1D_trackPt[iDataset]->GetBinContent(i);
        std::cout << "Bin " << i << ": [" << xlow << ", " << xhigh << "] -> " << content << std::endl;
    }

    // ===== Affichage des bins du nouvel histogramme rebinned (for debbug) =====
    std::cout << "=== Bins du nouvel histogramme rebinned ===" << std::endl;
    for (int i = 1; i <= H1D_trackPt_rebinned[iDataset]->GetNbinsX(); i++) {
        double xlow  = H1D_trackPt_rebinned[iDataset]->GetBinLowEdge(i);
        double xhigh = H1D_trackPt_rebinned[iDataset]->GetBinLowEdge(i) + H1D_trackPt_rebinned[iDataset]->GetBinWidth(i);
        double content = H1D_trackPt_rebinned[iDataset]->GetBinContent(i);
        std::cout << "Bin " << i << ": [" << xlow << ", " << xhigh << "] -> " << content << std::endl;
    }
    // ===================== End rebinning =====================


    if (options.find("evtNorm") != std::string::npos) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]);
      }
      NormaliseAsYieldToInputN(H1D_trackPt_rebinned[iDataset], Nevents);
      cout << "Dataset " << iDataset << ": Nevents = " << Nevents << endl;
    }
    if (options.find("entriesNorm") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackPt_rebinned[iDataset]);
    }
    if (options.find("20GevNorm") != std::string::npos) {
      double pt20GeV = 20;
      int iBinAt20GeV = H1D_trackPt_rebinned[iDataset]->GetXaxis()->FindBin(pt20GeV);
      int nTracksAt20GeV = H1D_trackPt_rebinned[iDataset]->GetBinContent(iBinAt20GeV);
      NormaliseAsYieldToInputN(H1D_trackPt_rebinned[iDataset], nTracksAt20GeV);
    }
  }
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
      divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_trackPt_rebinned[0], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
    }
  }
  TString textContext(contextTrackDatasetComp(""));

  TString* textYaxis;
  TString pdfNameNorm;
  if (options.find("evtNorm") != std::string::npos) {
    textYaxis = texTrackPtYield_EventNorm;
    pdfNameNorm = (TString)"_EventNorm";
  }
  if (options.find("entriesNorm") != std::string::npos) {
    textYaxis = texTrackPtYield_EntriesNorm;
    pdfNameNorm = (TString)"_EntriesNorm";
  }

  if (options.find("20GevNorm") != std::string::npos) {
    textYaxis = texTrackPtYield_20GevNorm;
    pdfNameNorm = (TString)"_20GevNorm";
  }

  TString* pdfName = new TString("track_Pt_DataComp"+pdfNameNorm);
  TString* pdfName_ratio = new TString("track_Pt_DataComp"+pdfNameNorm+"_ratio");
  TString* pdfName_ratio_zoom = new TString("track_Pt_DataComp"+pdfNameNorm+"_ratio_zoom");

  std::array<std::array<float, 2>, 2> drawnWindowCustomRatio_zoom = {{{0.1, 100}, {0.9, 1.1}}}; // {{xmin, xmax}, {ymin, ymax}}
  const std::array<std::array<float, 2>, 2> drawnWindowPt = {{{-999, -999}, {-999, -999}}}; // {{{xmin, xmax}, {ymin, ymax}}}
  const std::array<std::array<float, 2>, 2> legendPlacementPt = {{{0.7, 0.65}, {0.85, 0.85}}}; // {{{x1, y1}, {x2, y2}}}
  const std::array<std::array<float, 2>, 2> drawnWindowPtRatio = {{{-999, -999}, {0.3, 1.8}}}; // {{{xmin, xmax}, {ymin, ymax}}}
  const std::array<std::array<float, 2>, 2> legendPlacementPtRatio = {{{0.17, 0.72}, {0.45, 0.81}}}; // {{{x1, y1}, {x2, y2}}}



  Draw_TH1_Histograms(H1D_trackPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtX, textYaxis, texCollisionDataInfo, drawnWindowPt, legendPlacementPt, contextPlacementAuto, "logx,logy"+histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackPt_rebinned_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texPtX, texRatio, texCollisionDataInfo, drawnWindowPtRatio, legendPlacementPtRatio, contextPlacementAuto, "logx,zoomToOneMedium1");
    } else {
      Draw_TH1_Histograms(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst,logx"+histDatasetComparisonStructure);
      Draw_TH1_Histograms(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio_zoom, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowCustomRatio_zoom, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst,logx,zoomToOneMedium2"+histDatasetComparisonStructure);
    }
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}


void Draw_Eta_DatasetComparison(float* ptRange, std::string options) {
  TH2D* H2D_centrality_track[nDatasets];
  TH3D* H3D_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

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

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(4.,"trackEta_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);


    if (options.find("evtNorm") != std::string::npos) {
      if (isDatasetWeighted[iDataset]) {
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]);
      }
      NormaliseAsYieldToInputN(H1D_trackEta_rebinned[iDataset], Nevents);
    }
    if (options.find("entriesNorm") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackEta_rebinned[iDataset]);
    }
  }

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
  TString* pdfName_ratio_zoom = new TString((TString)"track_Eta_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]"+pdfNameNorm+"_ratio_zoom");
  TString* pdfName_ratio_zoom2 = new TString((TString)"track_Eta_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]"+pdfNameNorm+"_ratio_zoom2");

  std::array<std::array<float, 2>, 2> drawnWindowEta = {{{-1, 1}, {260, 390}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnWindowEtaZoom = {{{-1, 1}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnWindowEtaTwoByTwoRatio = {{{-1, 1}, {0.5, 1.3}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> legendPlacementCustom = {{{0.65, 0.65}, {0.85, 0.85}}}; // {{{x1, y1}, {x2, y2}}}
  std::array<std::array<float, 2>, 2> legendPlacementCustom2 = {{{0.2, 0.2}, {0.7, 0.38}}}; // {{{x1, y1}, {x2, y2}}}



  Draw_TH1_Histograms(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texEtaX, textYaxis, texCollisionDataInfo, drawnWindowEtaZoom, legendPlacementCustom, contextPlacementAuto, ""+histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackEta_rebinned_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texEtaX, texRatio, texCollisionDataInfo, drawnWindowEtaTwoByTwoRatio, legendPlacementCustom2, contextPlacementAuto, ",zoomToOneExtra");
    } else {
      Draw_TH1_Histograms(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst"+histDatasetComparisonStructure);
      Draw_TH1_Histograms(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio_zoom, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst,zoomToOneMedium2"+histDatasetComparisonStructure);
      Draw_TH1_Histograms(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio_zoom2, texEtaX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "noMarkerFirst,zoomToOneExtraExtra"+histDatasetComparisonStructure);
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
      H1D_trackPhi[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPhi_"+Datasets[iDataset]+DatasetsNames[iDataset], iCollSystem == 0 ? 1 : 0,  iCollSystem == 0 ? H2D_centrality_track[iDataset]->GetNbinsX() : -1, "e");
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
        Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]);
      } else {
        Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]);
      }
      NormaliseAsYieldToInputN(H1D_trackPhi_rebinned[iDataset], Nevents);
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
  
  std::array<std::array<float, 2>, 2> drawnWindowEtaTwoByTwoRatio = {{{-1, 7}, {0.5, 1.7}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> legendPlacementCustomRatio = {{{0.2, 0.2}, {0.5, 0.3}}}; // {{{x1, y1}, {x2, y2}}}
  std::array<std::array<float, 2>, 2> legendPlacementCustom = {{{0.65, 0.68}, {0.85, 0.85}}}; // {{{x1, y1}, {x2, y2}}}
  std::array<std::array<float, 2>, 2> drawnWindowEta = {{{-999, -999}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiX, textYaxis, texCollisionDataInfo, drawnWindowEta, legendPlacementCustom, contextPlacementAuto, "histWithLine"+histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackPhi_rebinned_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texPhiX, texRatio, texCollisionDataInfo, drawnWindowEtaTwoByTwoRatio, legendPlacementCustomRatio, contextPlacementAuto, "zoomToOneExtraExtra");
    } else {
      Draw_TH1_Histograms(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiX, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementCustom, contextPlacementAuto, "noMarkerFirst"+histDatasetComparisonStructure);
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
      
        // ===================== Advanced rebinning =====================
    std::vector<double> xbinsVector = GetTH1Bins(H1D_trackPt[iCentralityBin]);
    double* xbinsInitial = &xbinsVector[0];
    double ptBinsNew[500]; // 500 to be safe
    int iBinNew = 0;
    int iBinInitial = 0;
    double binWidth = H1D_trackPt[iCentralityBin]->GetXaxis()->GetBinWidth(1); //to be use in case we want to have larger bins in mid-pt range

    int increment;

    while (iBinInitial < H1D_trackPt[iCentralityBin]->GetNbinsX()) {
        double pt = xbinsInitial[iBinInitial];

        if (pt < 025) {
            increment = 1;
        }
        else if (pt < 100.0) {
            // increment = 1; 
            increment = 1 + pt / binWidth / 8; //in case we want to have larger bins in mid-pt range
        }
        else {
            //no need of increment because we want to stop at 100
            ptBinsNew[iBinNew++] = pt;   // left edge of the large bin
            ptBinsNew[iBinNew++] = 100.0; // right edge of the large bin
            break;
      }
        ptBinsNew[iBinNew++] = pt;
        iBinInitial += increment;
    }

    // security to ensure the last bin edge is 100 GeV
    if (ptBinsNew[iBinNew-1] < 100.0) ptBinsNew[iBinNew++] = 100.0;

    int nBinsNew = iBinNew - 1;

    // ======= Creation of the rebinned histogram =======
    H1D_trackPt_rebinned[iCentralityBin] = (TH1D*)H1D_trackPt[iCentralityBin]->Rebin(nBinsNew,"trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_@cent["+Form("%.1d", ibinCent_low)+","+Form("%.1d", ibinCent_high)+"]", ptBinsNew);

    NormaliseAsYieldToInputN(H1D_trackPt_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));

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



void Draw_Eta_CentralityComparison(int iDataset, std::string options) {

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

  TString* yAxisLabel;
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


    yAxisLabel = texNoNorm_SelMultYield;
    if (options.find("entriesNorm") != std::string::npos) {
      NormaliseYieldToIntegral(H1D_trackEta_rebinned[iCentralityBin]);
      yAxisLabel = texTrackEtaYield_EntriesNorm;
    }
    if (options.find("evtNorm") != std::string::npos) {
      NormaliseAsYieldToInputN(H1D_trackEta_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));
      yAxisLabel = texTrackEtaYield_EventNorm_CentWindow;
    }

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_Eta");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  std::array<std::array<float, 2>, 2> drawnWindow = {{{-999, -999}, {22, 33}}}; // {{xmin, xmax}, {ymin, ymax}}

  Draw_TH1_Histograms(H1D_trackEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  // Draw_TH1_Histograms(H1D_trackEta_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, texTrackEtaYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");

}

void Draw_Phi_CentralityComparison(int iDataset, std::string options) {

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
  
  TString* yAxisLabel;

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

    yAxisLabel = texNoNorm_SelMultYield;
    if (options.find("entriesNorm") != std::string::npos) {
      NormaliseYieldToIntegral(H1D_trackPhi_rebinned[iCentralityBin]);
      yAxisLabel = texTrackPhiYield_EntriesNorm;
    }
    if (options.find("evtNorm") != std::string::npos) {
      NormaliseAsYieldToInputN(H1D_trackPhi_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));
      yAxisLabel = texTrackPhiYield_EventNorm_CentWindow;
    }

    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_Phi");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_trackPhi_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "zoom");
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

  NormaliseAsYieldToInputN(H1D_trackPt_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], 0, 10));

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

  NormaliseAsYieldToInputN(H1D_trackEta_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], 0, 10));

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

  NormaliseAsYieldToInputN(H1D_trackPhi_rebinned, GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], 0, 10));

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
  TH1D* H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[nDatasets];

  double rebin_xMin = 0;
  double rebin_xMiddle = 10;
  double rebin_xMax = 100;
  int rebin_nLeft = 50;
  int rebin_nRight = 10;
  int nbinsX, nbinsY;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_sigmapt_pt[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_track_pt_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_DatasetComp_left"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_sigmapt_pt_high[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_track_pt_high_track_sigmapt"))->Clone("Draw_Sigmapt_vs_pt_high_DatasetComp_right"+Datasets[iDataset]+DatasetsNames[iDataset]);


    H2D_sigmapt_pt_high_rebinnedX[iDataset] = (TH2D*)H2D_sigmapt_pt_high[iDataset]->RebinX(1.,"H2D_sigmapt_pt_high_rebinnedX"+Datasets[iDataset]+DatasetsNames[iDataset]);

    // gustavo
    // nbinsY = H2D_sigmapt_pt[iDataset]->GetNbinsY();
    // for(int j = 1; j <= H2D_sigmapt_pt[iDataset]->GetNbinsX(); j++) {
    //   TH1D* temp1 = h2D->ProjectionY(Form("Bin%d",j),j,j);
    //   float scale1 = temp1 -> Integral(-1,-1);
      
    //   for(int i = 1; i <= nbinsy; i++) {
    //     if ( scale1 > 0 ) {
    //       h2D->SetBinContent(j,i, h2D->GetBinContent(j,i)/scale1);
    //       h2D->SetBinError  (j,i, h2D->GetBinError  (j,i)/scale1);
    //     } else {
    //       h2D->SetBinContent(j,i, 0);
    //       h2D->SetBinError  (j,i, 0);       
    //     }
    //   } // y bin loop 
    // } // x bin loop
  
    NormaliseXSlicesToOne(H2D_sigmapt_pt[iDataset]);
    NormaliseXSlicesToOne(H2D_sigmapt_pt_high_rebinnedX[iDataset]);


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
        // if (iDataset ==3){ //hardcoded for some old check
        //   H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(iBinX, iBinY, H2D_sigmapt_pt[iDataset]->GetBinContent(iBinX, iBinY) - H2D_sigmapt_pt[0]->GetBinContent(iBinX, iBinY));
        //   H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(iBinX, iBinY, sqrt(H2D_sigmapt_pt_concatenated[iDataset]->GetBinContent(iBinX, iBinY)));
        // } 
      }
    }

    for(int iBinX = 1; iBinX <= H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetNbinsX(); iBinX++){
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, H1D_sigmapt_pt_high_mean_withProfile[iDataset]->GetBinContent(iBinX));
      // H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, H1D_sigmapt_pt_high_mean_withProfile[iDataset]->GetBinError(iBinX));
      for(int iBinY = 1; iBinY <= H2D_sigmapt_pt[iDataset]->GetNbinsY(); iBinY++){
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetBinContent(iBinX, iBinY));
        H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetBinError(iBinX, iBinY));
        // if (iDataset ==3){//hardcoded for some old check
        //   H2D_sigmapt_pt_concatenated[iDataset]->SetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, H2D_sigmapt_pt_high_rebinnedX[iDataset]->GetBinContent(iBinX, iBinY) - H2D_sigmapt_pt_high_rebinnedX[0]->GetBinContent(iBinX, iBinY));
        //   H2D_sigmapt_pt_concatenated[iDataset]->SetBinError(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY, sqrt(H2D_sigmapt_pt_concatenated[iDataset]->GetBinContent(H2D_sigmapt_pt[iDataset]->GetNbinsX()+iBinX, iBinY)));
        // } 
      }
    }

    H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset] = (TH1D*)((TProfile*)H2D_sigmapt_pt_concatenated[iDataset]->ProfileX("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e"))->ProjectionX("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile_th1"+Datasets[iDataset]+DatasetsNames[iDataset]);



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
  TString* pdfName_mean_withProfile_ratio = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile_ratio");
  TString* pdfName_mean_withProfile_ratio_zoom = new TString("track_sigmapt_mean_vs_pt_DataComp_withProfile_ratio_zoom");

  // TString* pdfName_median = new TString("track_sigmapt_median_vs_pt_DataComp");
  // TString* pdfName_median_logy = new TString("track_sigmapt_median_vs_pt_DataComp_logy");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, *texDatasetsComparisonType, ""));
  std::array<std::array<float, 2>, 2> legendPlacement = {{{0.55, 0.65}, {0.65, 0.75}}}; // {{{x1, y1}, {x2, y2}}}

  std::array<std::array<float, 2>, 3> drawnWindow2DSigma = {{{0.1, 100}, {0.001, 100}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
  // Draw_TH2_Histograms(H2D_sigmapt_pt_concatenated, DatasetsNames, nDatasets, textContext, pdfName, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindowSigma, th2ContoursNone, contourNumberNone, "logx,logz,autoRangeSame"); // ?
  Draw_TH2_Histograms(H2D_sigmapt_pt_concatenated, DatasetsNames, nDatasets, textContext, pdfName_logy, texPtX, texSigmaPt, texCollisionDataInfo, drawnWindow2DSigma, th2ContoursNone, contourNumberNone, "logx,logy,logz,autoRangeSame"); // ?

  // Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx");
  std::array<std::array<float, 2>, 2> drawnWindowSigmaAverage = {{{0.1, 100}, {0.005, 2}}}; // {{xmin, xmax}, {ymin, ymax}}
  Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile_logy, texPtX, texSigmaPtMean, texCollisionDataInfo, drawnWindowSigmaAverage, legendPlacementAuto, contextPlacementAuto, "logx,logy"+histDatasetComparisonStructure);

  std::array<std::array<float, 2>, 2> legendPlacementCustom = {{{0.2, 0.15}, {0.4, 0.4}}}; // {{{x1, y1}, {x2, y2}}}

  bool divideSuccess;
  TString DatasetsNamesPairRatio[nDatasets];
  int nHistPairRatio = (int)nDatasets / 2;;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      if (iDataset < nHistPairRatio) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
        H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset] = (TH1D*)H1D_sigmapt_pt_mean_withProfile_concatenated[2*iDataset]->Clone("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile_ratios"+Datasets[2*iDataset]+DatasetsNames[2*iDataset]);
        H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset]->Reset("M");
        divideSuccess = H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset]->Divide(H1D_sigmapt_pt_mean_withProfile_concatenated[2*iDataset], H1D_sigmapt_pt_mean_withProfile_concatenated[2*iDataset+1], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      }
    } else if (histDatasetComparisonStructure.find("twoRatioSetsInOne") != std::string::npos) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[iDataset]+(TString)"/"+DatasetsNames[iDataset%2];
        H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset] = (TH1D*)H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->Clone("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]);
        H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset]->Reset("M");
        divideSuccess = H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset]->Divide(H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset], H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset%2], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
    } else {
      H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset] = (TH1D*)H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset]->Clone("H1D_sigmapt_pt_concatenated_rebinned_mean_withProfile_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset]->Reset("M");
      divideSuccess = H1D_sigmapt_pt_mean_withProfile_concatenated_ratios[iDataset]->Divide(H1D_sigmapt_pt_mean_withProfile_concatenated[iDataset], H1D_sigmapt_pt_mean_withProfile_concatenated[0], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
    }
  }
  if (divideSuccess == true) {
    std::array<std::array<float, 2>, 2> drawnWindowCustomRatio = {{{0.1, 100}, {0, 2}}}; // {{xmin, xmax}, {ymin, ymax}}
    std::array<std::array<float, 2>, 2> drawnWindowCustomRatio_zoom = {{{0.1, 100}, {0.7, 1.3}}}; // {{xmin, xmax}, {ymin, ymax}}
    std::array<std::array<float, 2>, 2> drawnWindowCustomRatio_twoRatioSetsInOne_zoom = {{{0.1, 100}, {0.92, 1.12}}}; // {{xmin, xmax}, {ymin, ymax}}
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_mean_withProfile_ratio, texPtX, texRatio, texCollisionDataInfo, drawnWindowCustomRatio, legendPlacementCustom, contextPlacementAuto, "logx"); //no histDatasetComparisonStructure because we actually don't want to plot the ratio with two by two colors etc, that's just for the non ratio plots
      Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_mean_withProfile_ratio_zoom, texPtX, texRatio, texCollisionDataInfo, drawnWindowCustomRatio_zoom, legendPlacementCustom, contextPlacementAuto, "logx"); //no histDatasetComparisonStructure because we actually don't want to plot the ratio with two by two colors etc, that's just for the non ratio plots
    } else if (histDatasetComparisonStructure.find("twoRatioSetsInOne") != std::string::npos){
      Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated_ratios, DatasetsNamesPairRatio, nDatasets, textContext, pdfName_mean_withProfile_ratio, texPtX, texRatio, texCollisionDataInfo, drawnWindowCustomRatio, legendPlacementCustom, contextPlacementAuto, "logx, noMarkerFirst, noMarkerSecond, ratioLine"+histDatasetComparisonStructure+"twoByTwoDatasetPairs");
      Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated_ratios, DatasetsNamesPairRatio, nDatasets, textContext, pdfName_mean_withProfile_ratio_zoom, texPtX, texRatio, texCollisionDataInfo, drawnWindowCustomRatio_twoRatioSetsInOne_zoom, legendPlacementCustom, contextPlacementAuto, "logx, noMarkerFirst, noMarkerSecond, ratioLine"+histDatasetComparisonStructure+"twoByTwoDatasetPairs");
    } else {
      Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated_ratios, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowCustomRatio, legendPlacementCustom, contextPlacementAuto, "logx, noMarkerFirst"+histDatasetComparisonStructure);
      Draw_TH1_Histograms(H1D_sigmapt_pt_mean_withProfile_concatenated_ratios, DatasetsNames, nDatasets, textContext, pdfName_mean_withProfile_ratio_zoom, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowCustomRatio_zoom, legendPlacementCustom, contextPlacementAuto, "logx, noMarkerFirst"+histDatasetComparisonStructure);
    }
  }
  else {
    cout << "Divide failed in Draw_Sigmapt_vs_pt_DatasetComp" << endl;
  }

  // Draw_TH1_Histograms(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx");
  // Draw_TH1_Histograms(H1D_sigmapt_pt_median, DatasetsNames, nDatasets, textContext, pdfName_median_logy, texPtX, texSigmaPtMedian, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
}


void Draw_SelectedMultiplicity_DatasetComp() {

  TH2D* H2D_rhoMult[nDatasets];
  TH1D* H1D_mult[nDatasets];
  TH1D* H1D_mult_rebinned[nDatasets];

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H2D_rhoMult[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_ntracks"))->Clone("Draw_NtracksPerEvent_DatasetComp"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_rhoMult[iDataset]->Sumw2();

    // H2D_rhoMult_rebinned[iDataset] = (TH2D*)H2D_rhoMult[iDataset]->RebinX(1.,"H2D_rhoMult_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_rhoMult_rebinned[iDataset]->GetXaxis()->SetRange(0,H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 1)); //(asks for the last bin on the x axis (axis number 1) to have strictly more than 1 entry)
    // H2D_rhoMult_rebinned[iDataset]->GetYaxis()->SetRange(1, H2D_rhoMult_rebinned[iDataset]->FindLastBinAbove(1, 2));
    // H2D_rhoMult_rebinned[iDataset]->Scale(1./H2D_rhoMult_rebinned[iDataset]->Integral(1, H2D_rhoMult_rebinned[iDataset]->GetNbinsX(), 1, H2D_rhoMult_rebinned[iDataset]->GetNbinsY(), "width"));

    H1D_mult[iDataset] = H2D_rhoMult[iDataset]->ProjectionY("H1D_mult"+(TString)Datasets[iDataset], iCollSystem == 0 ? 1 : 0, iCollSystem == 0 ? H2D_rhoMult[iDataset]->GetNbinsX() : -1);
    H1D_mult_rebinned[iDataset] = (TH1D*)H1D_mult[iDataset]->Rebin(1,"H1D_Mult_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }

  TString* pdfName = new TString("track_DataComp_multiplicitySelected");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_mult_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texSelectedMultiplicityTracks, texCount, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
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
      NormaliseYieldToIntegral(H1D_mult_rebinned[iCentralityBin]);
      yAxisLabel = texEntriesNorm_selMultYield;
    }
    if (options.find("normEvents") != std::string::npos) {
      int nEvents = H2D_multCentrality->Integral(1, H2D_multCentrality->GetNbinsX(), 1, H2D_multCentrality->GetNbinsY());
      // NormaliseAsYieldToInputN(H1D_mult_rebinned[iCentralityBin], GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], collHistIsObsolete[iDataset]));
      NormaliseAsYieldToInputN(H1D_mult_rebinned[iCentralityBin], nEvents);
      yAxisLabel = texCollNorm_selMultYield;
    }
    if (options.find("normEventsCentrality") != std::string::npos) {
      NormaliseAsYieldToInputN(H1D_mult_rebinned[iCentralityBin], GetNEventsSelectedCentrality_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflow[iDataset], arrayCentralityBinning[iCentralityBin], arrayCentralityBinning[iCentralityBin+1]));
      yAxisLabel = texCollNorm_selMultYield_CentWindow;
    }
  
    ss << "Cent " << arrayCentralityBinning[iCentralityBin] << " - " << arrayCentralityBinning[iCentralityBin+1] << " ";
    CentralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }

  TString* pdfName = new TString("track_CentralityComp_"+DatasetsNames[iDataset]+"_selMult_logy");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_mult_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName, texSelectedMultiplicityTracks, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoXrange,logy");
  TString* pdfName2 = new TString("track_CentralityComp_"+DatasetsNames[iDataset]+"_selMult");
  Draw_TH1_Histograms(H1D_mult_rebinned, CentralityLegend, nCentralityBins, textContext, pdfName2, texSelectedMultiplicityTracks, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoXrange");
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

void Draw_DcaXY_DatasetComp() { 

  TH1D* H1D_trackDcaXY[nDatasets];  
  TH1D* H1D_trackDcaXY_ratios[nDatasets];

  bool divideSuccess = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H1D_trackDcaXY[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/TrackPar/dcaXY"))->Clone("Draw_DcaXY_DatasetComp"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H1D_trackDcaXY[iDataset]->Sumw2();
  }

  cout << "tesst0" << endl;

  TString DatasetsNamesPairRatio[nDatasets];
  int nHistPairRatio = (int)nDatasets / 2;;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    cout << "tesst1, " << iDataset << endl;
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      if (iDataset < nHistPairRatio) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
        H1D_trackDcaXY_ratios[iDataset] = (TH1D*)H1D_trackDcaXY[2*iDataset]->Clone("trackDcaXY_rebinned_ratios"+Datasets[2*iDataset]+DatasetsNames[2*iDataset]);
        H1D_trackDcaXY_ratios[iDataset]->Reset("M");
        divideSuccess = H1D_trackDcaXY_ratios[iDataset]->Divide(H1D_trackDcaXY[2*iDataset], H1D_trackDcaXY[2*iDataset+1], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
      }
    } else {
      H1D_trackDcaXY_ratios[iDataset] = (TH1D*)H1D_trackDcaXY[iDataset]->Clone("trackDcaXY_rebinned_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H1D_trackDcaXY_ratios[iDataset]->Reset("M");
      divideSuccess = H1D_trackDcaXY_ratios[iDataset]->Divide(H1D_trackDcaXY[iDataset], H1D_trackDcaXY[0], 1., 1., datasetsAreSubsetsofId0 ? "b" : "");
    }
  }

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* pdfName = new TString((TString)"track_DcaXY_DataComp");
  TString* pdfName_ratio = new TString((TString)"track_DcaXY_DataComp_ratio");
  
  std::array<std::array<float, 2>, 2> drawnWindowCustom = {{{-999, -999}, {0.99, 1.01}}}; // {{{x1, y1}, {x2, y2}}}
  cout << "tesst2" << endl;

  Draw_TH1_Histograms(H1D_trackDcaXY, DatasetsNames, nDatasets, textContext, pdfName, texDCAXY, texCount, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, histDatasetComparisonStructure);
  if (divideSuccess == true) {
    if (histDatasetComparisonStructure.find("twoByTwoDatasetPairs") != std::string::npos) {
      Draw_TH1_Histograms(H1D_trackDcaXY_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texDCAXY, texRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "zoomToOneExtraExtra");
    } else {
      Draw_TH1_Histograms(H1D_trackDcaXY_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texDCAXY, texRatioDatasets, texCollisionDataInfo, drawnWindowCustom, legendPlacementAuto, contextPlacementAuto, "zoomToOneExtraExtra, noMarkerFirst"+histDatasetComparisonStructure);
    }
  }
  else {
    cout << "Divide failed in Draw_DcaXY_DatasetComp" << endl;
  }
  cout << "tesst3" << endl;

}



void Draw_TrackSelDistribs_DatasetComparison(float* ptRange) {

  TString textContext(contextPtRange(ptRange));

  TString* textYaxis;
  TString pdfNameNorm;

  float ptCutLow = ptRange[0];
  float ptCutHigh = ptRange[1];

  double Nevents; 

  const int nHistVariables = 6;
  std::array<TString, nHistVariables> histVariableName{"tpccrossedrows", "tpccrossedrowsoverfindable", "chi2ncls_tpc", "chi2ncls_its", "dcaxy", "dcaz"};

  TH2D* H2D_trackselplot[nHistVariables][nDatasets];
  TH1D* H1D_trackselplot[nHistVariables][nDatasets];
  TH1D* H1D_trackselplot_rebinnedX[nHistVariables][nDatasets];
  
  TString* texXaxis[nHistVariables] = {texCrossedRows, texCrossedRowsOverFindable, texChi2NclsTPC, texChi2NclsITS, texDcaXY, texDcaZ};
  TString* texYaxis[nHistVariables] = {texTrackCrossedRowsYield_EntriesNorm, texTrackCrossedRowsOverFindableYield_EntriesNorm, texTrackChi2NclsTPCYield_EntriesNorm, texTrackChi2NclsITSYield_EntriesNorm, texTrackDcaXYYield_EntriesNorm, texTrackDcaZYield_EntriesNorm};

  for (int iHistVariable = 0; iHistVariable < nHistVariables; iHistVariable++){ 
    for(int iDataset = 0; iDataset < nDatasets; iDataset++){
      H2D_trackselplot[iHistVariable][iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_trackselplot_pt_"+histVariableName[iHistVariable]))->Clone("Draw_htrackselplot_"+histVariableName[iHistVariable]+"_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
      H2D_trackselplot[iHistVariable][iDataset]->Sumw2();


      int ibinPt_low = H2D_trackselplot[iHistVariable][iDataset]->GetXaxis()->FindBin(ptCutLow);
      int ibinPt_high = H2D_trackselplot[iHistVariable][iDataset]->GetXaxis()->FindBin(ptCutHigh);
      if (ibinPt_low == 0) 
        cout << "WARNING: Eta_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
      if (ibinPt_high == H2D_trackselplot[iHistVariable][iDataset]->GetXaxis()->GetNbins()+1) 
        cout << "WARNING: Eta_DatasetComparison is counting the overflow with the chosen PtRange" << endl;

      H1D_trackselplot[iHistVariable][iDataset] = (TH1D*)H2D_trackselplot[iHistVariable][iDataset]->ProjectionY("Draw_htrackselplot_"+histVariableName[iHistVariable]+Datasets[iDataset]+DatasetsNames[iDataset], ibinPt_low, ibinPt_high, "e");
      H1D_trackselplot_rebinnedX[iHistVariable][iDataset] = (TH1D*)H1D_trackselplot[iHistVariable][iDataset]->RebinX(2.,"Draw_htrackselplot_rebinnedX"+histVariableName[iHistVariable]+Datasets[iDataset]+DatasetsNames[iDataset]);

      // ===================== Basic rebinning =====================
      // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(2.,"trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]);
      
      // int nBinsLogRough = 100;
      // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(H1D_trackPt[iDataset], nBinsLogRough);
      // int nBinsLogResult = O2H1DPtLogBinsVector.size()-1;
      // // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(0.5, 100, nBinsLogRough);
      // double* O2ptLogBins = &O2H1DPtLogBinsVector[0];
      // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(nBinsLogResult, "trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset], O2ptLogBins);

      // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);


      NormaliseYieldToIntegral(H1D_trackselplot_rebinnedX[iHistVariable][iDataset]);
    }

    TString* pdfName = new TString("track_selectionDistrib_"+histVariableName[iHistVariable]+"_EntriesNorm_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]");
    TString* pdfName_logy = new TString("track_selectionDistrib_"+histVariableName[iHistVariable]+"_EntriesNorm_DataComp_@pT["+Form("%03.0f", ptCutLow)+","+Form("%03.0f", ptCutHigh)+"]_logy");

    std::array<std::array<float, 2>, 2> drawnWindowCustomRatio_zoom = {{{0.1, 100}, {0.9, 1.1}}}; // {{xmin, xmax}, {ymin, ymax}}
    const std::array<std::array<float, 2>, 2> drawnWindow = {{{-999, -999}, {-999, -999}}}; // {{{xmin, xmax}, {ymin, ymax}}}
    const std::array<std::array<float, 2>, 2> legendPlacement = {{{0.7, 0.65}, {0.85, 0.85}}}; // {{{x1, y1}, {x2, y2}}}
    const std::array<std::array<float, 2>, 2> legendPlacementLog = {{{0.7, 0.40}, {0.85, 0.55}}}; // {{{x1, y1}, {x2, y2}}}



    Draw_TH1_Histograms(H1D_trackselplot_rebinnedX[iHistVariable], DatasetsNames, nDatasets, textContext, pdfName, texXaxis[iHistVariable], texYaxis[iHistVariable], texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacementAuto, "histWithLine, smallMarkers"+histDatasetComparisonStructure);
    Draw_TH1_Histograms(H1D_trackselplot_rebinnedX[iHistVariable], DatasetsNames, nDatasets, textContext, pdfName_logy, texXaxis[iHistVariable], texYaxis[iHistVariable], texCollisionDataInfo, drawnWindow, legendPlacementLog, contextPlacementAuto, "histWithLine, smallMarkers, logy"+histDatasetComparisonStructure);
  }
}

void Draw_TrackSelDistribs_PtComparison(int iDataset, float* ptBins, int nBinsPt) {

  TString textContext(contextTrackDatasetComp("")+" "+DatasetsNames[iDataset]);

  TString* textYaxis;
  TString pdfNameNorm;

  double Nevents; 

  const int nHistVariables = 6;
  std::array<TString, nHistVariables> histVariableName{"tpccrossedrows", "tpccrossedrowsoverfindable", "chi2ncls_tpc", "chi2ncls_its", "dcaxy", "dcaz"};

  TH2D* H2D_trackselplot[nHistVariables];
  TH1D* H1D_trackselplot[nHistVariables][nBinsPt];
  TString* texXaxis[nHistVariables] = {texCrossedRows, texCrossedRowsOverFindable, texChi2NclsTPC, texChi2NclsITS, texDcaXY, texDcaZ};
  TString* texYaxis[nHistVariables] = {texTrackCrossedRowsYield_EntriesNorm, texTrackCrossedRowsOverFindableYield_EntriesNorm, texTrackChi2NclsTPCYield_EntriesNorm, texTrackChi2NclsITSYield_EntriesNorm, texTrackDcaXYYield_EntriesNorm, texTrackDcaZYield_EntriesNorm};

  TString ptLegend[nBinsPt];
  std::stringstream ss;

  for (int iHistVariable = 0; iHistVariable < nHistVariables; iHistVariable++){ 
    H2D_trackselplot[iHistVariable] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_trackselplot_pt_"+histVariableName[iHistVariable]))->Clone("Draw_htrackselplot_"+histVariableName[iHistVariable]+"_DatasetComparison"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_trackselplot[iHistVariable]->Sumw2();
    for(int iBinPt = 0; iBinPt < nBinsPt; iBinPt++){

      int ibinPt_low = H2D_trackselplot[iHistVariable]->GetXaxis()->FindBin(ptBins[iBinPt]);
      int ibinPt_high = H2D_trackselplot[iHistVariable]->GetXaxis()->FindBin(ptBins[iBinPt+1]);
      if (ibinPt_low == 0) 
        cout << "WARNING: Eta_DatasetComparison is counting the underflow with the chosen PtRange" << endl;
      if (ibinPt_high == H2D_trackselplot[iHistVariable]->GetXaxis()->GetNbins()+1) 
        cout << "WARNING: Eta_DatasetComparison is counting the overflow with the chosen PtRange" << endl;

      H1D_trackselplot[iHistVariable][iBinPt] = (TH1D*)H2D_trackselplot[iHistVariable]->ProjectionY("Draw_htrackselplot_"+histVariableName[iHistVariable]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%03.0f", ptBins[iBinPt])+","+Form("%03.0f", ptBins[iBinPt+1]), ibinPt_low, ibinPt_high, "e");

      if (iHistVariable == 5) {
        H1D_trackselplot[iHistVariable][iBinPt]->Rebin(2.);
      }
      // ===================== Basic rebinning =====================
      // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(2.,"trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset]);
      
      // int nBinsLogRough = 100;
      // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(H1D_trackPt[iDataset], nBinsLogRough);
      // int nBinsLogResult = O2H1DPtLogBinsVector.size()-1;
      // // std::vector<double> O2H1DPtLogBinsVector = MakeVariableBinning_logarithmic(0.5, 100, nBinsLogRough);
      // double* O2ptLogBins = &O2H1DPtLogBinsVector[0];
      // H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(nBinsLogResult, "trackPt_rebinned_"+Datasets[iDataset]+DatasetsNames[iDataset], O2ptLogBins);

      // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);


      NormaliseYieldToIntegral(H1D_trackselplot[iHistVariable][iBinPt]);

      ss << "pt " << ptBins[iBinPt] << " - " << ptBins[iBinPt+1] << " GeV";
      ptLegend[iBinPt] = (TString)ss.str();
      ss.str("");
      ss.clear();
    }

    TString* pdfName = new TString("track_selectionDistrib_"+histVariableName[iHistVariable]+"_EntriesNorm_PtComp_"+DatasetsNames[iDataset]);
    TString* pdfName_logy = new TString("track_selectionDistrib_"+histVariableName[iHistVariable]+"_EntriesNorm_PtComp_"+DatasetsNames[iDataset]+"_logy");

    std::array<std::array<float, 2>, 2> drawnWindowCustomRatio_zoom = {{{0.1, 100}, {0.9, 1.1}}}; // {{xmin, xmax}, {ymin, ymax}}
    const std::array<std::array<float, 2>, 2> drawnWindow = {{{-999, -999}, {-999, -999}}}; // {{{xmin, xmax}, {ymin, ymax}}}
    const std::array<std::array<float, 2>, 2> legendPlacement = {{{0.7, 0.65}, {0.85, 0.85}}}; // {{{x1, y1}, {x2, y2}}}
    const std::array<std::array<float, 2>, 2> legendPlacementLog = {{{0.7, 0.65}, {0.85, 0.85}}}; // {{{x1, y1}, {x2, y2}}}
    const std::array<std::array<float, 2>, 2> legendPlacementLog_crossedRows = {{{0.7, 0.30}, {0.85, 0.45}}}; // {{{x1, y1}, {x2, y2}}}
    

    Draw_TH1_Histograms(H1D_trackselplot[iHistVariable], ptLegend, nBinsPt, textContext, pdfName, texXaxis[iHistVariable], texYaxis[iHistVariable], texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacementAuto, "histWithLine, smallMarkers"+histDatasetComparisonStructure);
    Draw_TH1_Histograms(H1D_trackselplot[iHistVariable], ptLegend, nBinsPt, textContext, pdfName_logy, texXaxis[iHistVariable], texYaxis[iHistVariable], texCollisionDataInfo, drawnWindow, (iHistVariable == 0 || iHistVariable == 1) ? legendPlacementLog_crossedRows : legendPlacementLog, contextPlacementAuto, "histWithLine, smallMarkers, logy"+histDatasetComparisonStructure);
  }
}