#ifndef JETSPECTRUM_DRAWINGFUNCTIONS_C
#define JETSPECTRUM_DRAWINGFUNCTIONS_C

#include "TStyle.h"
#include "TGraph.h"
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
#include "TMultiGraph.h"
#include <RooUnfold.h> // one should likely do `aliBuild build RooUnfold` then `alienv enter RooUnfold/latest` as alidist roounfold version is usually quite old
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"

//My Libraries
#include "./JetSpectrum_settings.h"
#include "./JetSpectrum_inputs.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_Unfolding.h"
#include "./JetSpectrum_Unfolding.C"
#include "./JetSpectrum_EfficiencyPurityGetters.h"
#include "./JetSpectrum_EfficiencyPurityGetters.C"
#include "./JetSpectrum_systematics.C"

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 
#include "../Utilities/Fits.C"
#include "../Utilities/Fits.h" 

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <stdlib.h>     /* abort, NULL */
#include <fstream>
using namespace std;


/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

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

void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step){
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin)/step) + 1;
  std::stringstream ss;
  ss.precision(2);
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    ss << "k_{unfold} = " << unfoldIterationMax - iUnfoldIteration * step; 
    iterationLegend[iUnfoldIteration] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// Spectrum plotting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_raw;

  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_raw, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_raw, iDataset, iRadius, options);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_raw");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots) {
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    yAxisLabel = texCount;
    *pdfName = *pdfName+(TString)"_noEventNormNorBinWidthScaling";
  }

  Draw_TH1_Histogram(H1D_jetPt_raw, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
}

void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcp_genBinning;
  TH1D* H1D_jetPt_mcp_recBinning;
  TH1D* H1D_jetPt_mcp_collection[2];

  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_genBinning, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_recBinning, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_genBinning, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_recBinning, iDataset, iRadius, options);
  }
  
  H1D_jetPt_mcp_collection[0] = H1D_jetPt_mcp_genBinning;
  H1D_jetPt_mcp_collection[1] = H1D_jetPt_mcp_recBinning;


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcp");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots) {
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    yAxisLabel = texCount;
    *pdfName = *pdfName+(TString)"_noEventNormNorBinWidthScaling";
  }
  TString genVsRecBinningLegend[2] = {"gen binning", "rec binning"};

  Draw_TH1_Histograms(H1D_jetPt_mcp_collection, genVsRecBinningLegend, 2, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy"); 
}

void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcdMatched;
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_mcdMatched_genBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);
  }


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcdMatched");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots) {
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    yAxisLabel = texCount;
    *pdfName = *pdfName+(TString)"_noEventNormNorBinWidthScaling";
  }

  Draw_TH1_Histogram(H1D_jetPt_mcdMatched, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
}


void Draw_Pt_efficiency_jets(int iRadius, std::string options) {
  TH1D* H1D_jetEfficiency[nDatasets];
  bool divideSuccess[nDatasets];
  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    if (useFineBinningTest) {
      divideSuccess[iDataset] = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency[iDataset], iDataset, iRadius, options);
    } else {
      divideSuccess[iDataset] = Get_Pt_JetEfficiency(H1D_jetEfficiency[iDataset], iDataset, iRadius, options);
    }
  }

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_efficiency");
  if (std::all_of(std::begin(divideSuccess), std::end(divideSuccess), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccess are true
    Draw_TH1_Histograms(H1D_jetEfficiency, DatasetsNames, nDatasets, textContext, pdfName, texPtJetGen, texJetEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");

    if (writeOutputRootFile_efficiency) {
      cout << "######################### FILE CREATED IN PRINCIPLE##################" << endl; 
      TFile* outFile_eff = new TFile("output.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
      TString histoName = Form("H1D_jetEfficiency_%s", MC_Datasets[0].Data());
      H1D_jetEfficiency[0]->Write(histoName);
      outFile_eff->Close();
      cout << "######################### HISTO SAVED IN PRINCIPLE ##################" << endl; 
    }
  }
}

void Draw_kinematicEfficiency(int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations[nDatasets];
  TH2D* H2D_jetPtResponseMatrix_detectorResponse[nDatasets];
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[nDatasets];
  TH1D* H1D_kinematicEfficiency[nDatasets];

  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

    Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse[iDataset], iDataset, iRadius, "");

    Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations[iDataset], iDataset, iRadius, "");
    Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[iDataset], H2D_jetPtResponseMatrix_detectorResponse[iDataset], H2D_jetPtResponseMatrix_fluctuations[iDataset], iDataset, iRadius, options);

    Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency[iDataset], H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[iDataset], name_H1D_kinematicEfficiency, iRadius);
  }
  TString priorInfo = (TString)unfoldingPrior;

  TString partialUniqueSpecifier = (TString)"R="+Form("%.1f",arrayRadius[iRadius]);
  TString* pdfName = new TString("kinematicEfficiency_"+partialUniqueSpecifier+priorInfo);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_kinematicEfficiency, DatasetsNames, nDatasets, textContext, pdfName, texPtJetGen, texJetKinematicEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_FakeRatio(int iRadius, std::string options) {
  TH1D* H1D_fakeRatio[nDatasets];

  TString partialUniqueSpecifier = (TString)" R="+Form("%.1f",arrayRadius[iRadius]);

  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    Get_Pt_JetFakes(H1D_fakeRatio[iDataset], iDataset, iRadius, options);
  }
  TString* pdfName = new TString("fakeRatio_"+partialUniqueSpecifier);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_fakeRatio, DatasetsNames, nDatasets, textContext, pdfName, texPtJetRec, texFakeRatio, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "");

  if (writeOutputRootFile_efficiency) {
      cout << "######################### FILE CREATED IN PRINCIPLE##################" << endl; 
      TFile* outFile_eff = new TFile("output.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
      TString histoName = Form("H1D_fakeRatio%s", MC_Datasets[0].Data());
      H1D_fakeRatio[0]->Write(histoName);
      outFile_eff->Close();
      cout << "######################### HISTO SAVED IN PRINCIPLE ##################" << endl; 
    }
}

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");

  TString priorInfo = (TString)unfoldingPrior;

  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/ResponseMatrices", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/ResponseMatrices", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/ResponseMatrices", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
  //     mkdir("pdfFolder/ResponseMatrices", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
  //     mkdir("pngFolder/ResponseMatrices", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");
  // TString* pdfNameFullRes_logz = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"FullRes_logz");
  TString* pdfName = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  // TString* pdfNameFullRes = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_FullRes");


  TString texCombinedMatrix = contextCustomOneField((TString)"ALICE Performance", ""); // Response matrix - "+(TString)*texEnergy
  TString textContextMatrixDetails = contextCustomFiveFields((TString)"Bkg. fluctuation response ", "", (TString)*texCollisionDataType, (TString)*texEnergyPbPb, contextJetRadius(arrayRadius[iRadius]), "");

  // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
  TH2D* MatrixResponse;
  TString* xLabel;
  TString* yLabel;
  if (transposeResponseHistogramsInDrawing) {
    MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_fluctuations).Clone("Draw_ResponseMatrices_Fluctuations"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetFluctCorrectedX;
    yLabel = texPtJetRec;
  } else {
    MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_fluctuations->Clone("Draw_ResponseMatrices_Fluctuations"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetRec;
    yLabel = texPtJetFluctCorrectedX;
  }

  double th2ContourCustom[1] = {0.000001}; // hardcoded at 10-6 for now
  int contourNumberCustom = 1;

  // std::array<std::array<float, 2>, 3> drawnWindowYaxianRequest = {{{-999, -999}, {-999, -999}, {1e-6, 1e-1}}}; // {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}} /// put AUTO again after perf figure is done

  // Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName_logz, texPtJetRec, texPtJetFluctCorrectedX, texCollisionDataInfo, drawnWindow2DAuto, th2ContourCustom, contourNumberCustom, "logz");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
}

void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  cout << "Draw_ResponseMatrices_detectorResponse 1" << endl;
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "");
  cout << "Draw_ResponseMatrices_detectorResponse 2" << endl;

  TString priorInfo = (TString)unfoldingPrior;


  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/ResponseMatrices", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/ResponseMatrices", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/ResponseMatrices", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
  //     mkdir("pdfFolder/ResponseMatrices", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
  //     mkdir("pngFolder/ResponseMatrices", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString* pdfName = new TString("ResponseMatrices/responseMatrix_detectorEffects_"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_detectorEffects_"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

  TString texCombinedMatrix = contextCustomOneField((TString)"ALICE Simulation", ""); // Response matrix - "+(TString)*texEnergy
  TString textContextMatrixDetails = contextCustomFiveFields((TString)"Detector response ", "", (TString)*texCollisionMCType, (TString)*texEnergy, (TString)contextJetRadius(arrayRadius[iRadius]), "");


  // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
  TH2D* MatrixResponse;
  TString* xLabel;
  TString* yLabel;
  if (transposeResponseHistogramsInDrawing) {
    MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorResponse).Clone("Draw_ResponseMatrices_detectorResponse"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetGen;
    yLabel = texPtJetRec;
  } else {    MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorResponse->Clone("Draw_ResponseMatrices_detectorResponse"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetRec;
    yLabel = texPtJetGen;
  }
  // std::array<std::array<float, 2>, 3> drawnWindowRaymondRequest = {{{-999, -999}, {-999, -999}, {1e-5, 6e-1}}}; // {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}} /// put AUTO again after perf figure is done

  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
}

void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;


  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "");
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix_priorAndNormYslicesAndMergeBins(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  TString priorInfo = (TString)unfoldingPrior;


  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/ResponseMatrices", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/ResponseMatrices", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/ResponseMatrices", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
  //     mkdir("pdfFolder/ResponseMatrices", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
  //     mkdir("pngFolder/ResponseMatrices", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

  TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
  TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

  // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
  TH2D* MatrixResponse;
  TString* xLabel;
  TString* yLabel;
  if (transposeResponseHistogramsInDrawing) {
    MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined).Clone("Draw_ResponseMatrices_DetectorAndFluctuationsCombined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetGen;
    yLabel = texPtJetRec;
  } else {
    MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("Draw_ResponseMatrices_DetectorAndFluctuationsCombined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetRec;
    yLabel = texPtJetGen;
  }

  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");

  // save matrix in .root file
  TString fileName = "jetPtResponseMatrixCombined";
  Save_PtResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, fileName);
}

void Draw_Pt_spectrum_unfolded_singleDataset(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  bool splitTestControlMC = true;

  TH1D* H1D_jetPt_measured;
  TH1D* H1D_jetPt_measured_genBinning;
  TH1D* H1D_jetPt_unfolded;
  TH1D* H1D_jetPt_unfoldedThenRefolded;
  TH1D* H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod;
  TH1D* H1D_jetPt_mcpFolded;
  TH1D* H1D_jetPt_mcpFolded2;
  TH1D* H1D_jetPt_mcpFoldedThenUnfolded;
  TH1D* H1D_jetPt_unfolded_mcpComp[2];
  TH1D* H1D_jetPt_unfolded_run2Comp_fitRebin[3];
  TH1D* H1D_jetPt_unfolded_run2Comp_shapeComp[2];
  TH1D* H1D_jetPt_unfolded_run2Comp[3];
  TH1D* H1D_jetPt_unfolded_measuredComp[2];
  TH1D* H1D_jetPt_unfolded_refoldedComp[3];
  TH1D* H1D_jetPt_unfolded_mcpFoldedComp[2];
  TH1D* H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[2];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcp_recBinControl;
  TH1D* H1D_jetPt_run2_HannaBossiLauraFile;
  TGraph* Graph_jetPt_run2_MLPaperFile;
  TH1D* H1D_jetPt_run2_MLPaperFile = new TH1D("H1D_jetPt_run2_MLPaperFile", "H1D_jetPt_run2_MLPaperFile", nBinPtJetsGen_run2[iRadius], ptBinsJetsGen_run2[iRadius]);
  TH1D* H1D_jetPt_run2_MLPaperFile_rebinned;
  std::vector<TGraphErrors*> TGraph_jetPt_run2_MLPaperFile_fit = {};
  std::vector<TGraphErrors*> TGraph_jetPt_unfolded_run2Comp_fits = {};
  TH1D* H1D_jetPt_ratio_mcp;
  TH1D* H1D_jetPt_ratio_run2_fitRebin[2];
  TH1D* H1D_jetPt_ratio_run2_shapeComp[2];
  TH1D* H1D_jetPt_ratio_run2[2];
  TH1D* H1D_jetPt_unfolded_run2Comp_xT[2];
  TH1D* H1D_jetPt_ratio_run2Comp_xT;
  TH1D* H1D_jetPt_unfolded_run2Comp_fits[2];
  TH1D* H1D_jetPt_ratio_run2Comp_fits;
  TH1D* H1D_jetPt_ratio_measured;
  TH1D* H1D_jetPt_ratio_measuredRefolded[2];
  TH1D* H1D_jetPt_ratio_mcpFoldedMcp;
  TH1D* H1D_jetPt_ratio_mcpFoldedUnfoldedMcp;

  TH1D* measuredInput_mcSplitInput;
  TH1D* H1D_jetPt_mcp_mcSplitInput;
  TH1D* H1D_jetPt_mcp_mcSplitInput_recBinning;
  TH1D* H1D_jetPt_unfolded_inputSplitClosure;
  TH1D* H1D_jetPt_unfolded_mcdSplitClosure[2];
  TH1D* H1D_jetPt_ratio_mcdSplitClosure;
  // RUN 2 settings
  if (comparePbPbWithRun2) {
    H1D_jetPt_run2_HannaBossiLauraFile = (TH1D*)((TH1D*)file_O2Analysis_run2ComparisonFileHannaBossiLaura->Get("Bayesian_Unfoldediter15"))->Clone("H1D_jetPt_run2_HannaBossiLauraFile");
    int NcollRun2 = 4619963; // central (see Laura discussion mattermost) 
    H1D_jetPt_run2_HannaBossiLauraFile->Scale(1./NcollRun2);

    double Ncoll;
    if (centralityRange[0] == 00 && centralityRange[1] == 10) {
      // Ncoll = (1780.9+1387.0)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (1956+1722+1521+1346)/4; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else if (centralityRange[0] == 50 && centralityRange[1] == 70) {
      // Ncoll = (103.7+46.1)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (89.8+39.8)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else {
      cout << "comparison with run2: Ncoll hasn't been calculated for this centrality interval" << endl;
    }
    double sigmaNN = 67.6; // value for sqrt(s) = 5.02 TeV https://arxiv.org/abs/1710.07098
    double T_AA = Ncoll / sigmaNN;
    Graph_jetPt_run2_MLPaperFile = ((TGraph*)((TDirectoryFile*)file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObjectAny("Graph1D_y1")); // https://doi.org/10.1016/j.physletb.2023.138412
    // H1D_jetPt_run2_MLPaperFile = (TH1D*)((TH1D*)(file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObject("Graph1D_y1"))->Clone("H1D_jetPt_run2_MLPaperFile");
    int Ngraph = Graph_jetPt_run2_MLPaperFile->GetN();
    for (int i=0; i < Ngraph; ++i) // setting bin contents to the TGraph values
    {
      double x,y;
      Graph_jetPt_run2_MLPaperFile->GetPoint(i, x, y);
      H1D_jetPt_run2_MLPaperFile->Fill(x, y); // uncertainties are of course screwed up
      int iHist = H1D_jetPt_run2_MLPaperFile->GetXaxis()->FindBin(x);
      H1D_jetPt_run2_MLPaperFile->SetBinError(iHist, Graph_jetPt_run2_MLPaperFile->GetErrorY(i));
    }
    H1D_jetPt_run2_MLPaperFile->Scale(T_AA);
    // now the spectre from the file is 1/N d2N/dpTdeta, instead of 1/T_AA 1/N d2N/dpTdeta
  }

  bool divideSuccessMcp;
  bool divideSuccessRun2_fitRebin[2];
  bool divideSuccessRun2_shapeComp;
  bool divideSuccessRun2[2];
  bool divideSuccessRun2_xt;
  bool divideSuccessRun2_fits;  
  bool divideSuccessMeasured;
  bool divideSuccessMeasuredRefolded[2];
  bool divideSuccessMcpFoldedMcp;
  bool divideSuccessMcpFoldedUnfoldedMcp;
  bool divideSuccessMcdSplitClosure;
  TString partialUniqueSpecifier;

  if (!useFineBinningTest) {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
    if (doClosure_splitMC_mcdUnfoldedVsGen) {
      Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
    }
  } else {
    Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
    if (doClosure_splitMC_mcdUnfoldedVsGen) {
      Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
    }
  }
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  int unfoldParameter, unfoldParameter_mcSplitClosure;

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    } else {
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    }
  } else{
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    } else {
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    }
  }

  unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).first;
  // TH1D* H1D_jetPt_unfolded2 = (TH1D*)H1D_jetPt_unfolded->Clone(H1D_jetPt_unfolded->GetName()+(TString)"H1D_jetPt_unfolded2");
  ////////////////////////////////////////////////////////////
  
  if (writeOutputRootFile) {
    cout << "######################### FILE CREATED IN PRINCIPLE##################" << endl; 
    TFile* outFile = new TFile("output.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
    TString histoName = Form("H1D_Pt_Unfolded_w_%s", mcIsWeighted ? "JJ" : "MB");
    H1D_jetPt_unfolded->Write(histoName);
    outFile->Close();
    cout << "######################### HISTO SAVED IN PRINCIPLE ##################" << endl; 
  }
  ////////////////////////////////////////////////////////////

  cout << "comparison with measured" << endl; 
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);
  }
  H1D_jetPt_unfolded_measuredComp[0] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_measuredComp[1] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measured = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
  // divideSuccessMeasured = H1D_jetPt_ratio_measured->Divide(H1D_jetPt_measured_genBinning);
  divideSuccessMeasured = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measured, H1D_jetPt_measured_genBinning);
  

  cout << "comparison with mcp truth" << endl; 
  H1D_jetPt_unfolded_mcpComp[0] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_mcp_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_mcpComp[1] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_mcp = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
  divideSuccessMcp = H1D_jetPt_ratio_mcp->Divide(H1D_jetPt_mcp);

  cout << "comparison with run2" << endl; 
  std::vector<double> xtBinningVectorRun2 = {};
  std::vector<double> xtBinningVectorRun3 = {};
  if (comparePbPbWithRun2) {
    // comparison with run2 results rebinned using a fit (errors look way underestimated; tsallis function not great aboe 100+ GeV; where should one eval the function inside a bin? probably not just the center)
    H1D_jetPt_unfolded_run2Comp_fitRebin[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp_fitRebin"+partialUniqueSpecifier);
    double fitPtRange[2] = {ptBinsJetsGen_run2[iRadius][0], ptMaxFit}; //-1 because Tsallis shape only accurate until 120GeV or so
    std::tuple<TF1*, TMatrixDSym, TFitResultPtr> exponentialFitWithLogTransfoResult_run2 = ExponentialFitWithLogTransfo(H1D_jetPt_run2_MLPaperFile, fitPtRange);
    TString histName_run2 = "rebinWithFit_run2"+partialUniqueSpecifier;
    std::pair<TH1D*, TGraphErrors*> pairResult_run2 = RebinWithFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange, histName_run2, exponentialFitWithLogTransfoResult_run2);
    H1D_jetPt_run2_MLPaperFile_rebinned = pairResult_run2.first;
    TGraph_jetPt_run2_MLPaperFile_fit.push_back(pairResult_run2.second);
    H1D_jetPt_unfolded_run2Comp_fitRebin[1] = (TH1D*)H1D_jetPt_run2_MLPaperFile_rebinned->Clone("H1D_jetPt_unfolded_run2_rebinned_fitRebin"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_fitRebin[2] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2_fitRebin"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_unfolded_run2Comp_HannaBossiLauraFile"+partialUniqueSpecifier);
    H1D_jetPt_ratio_run2_fitRebin[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_ratio_run2_fitRebin"+partialUniqueSpecifier);
    // H1D_jetPt_ratio_run2[1] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_ratio_run2_HannaBossiLauraFile"+partialUniqueSpecifier);
    divideSuccessRun2_fitRebin[0] = H1D_jetPt_ratio_run2_fitRebin[0]->Divide(H1D_jetPt_run2_MLPaperFile_rebinned);
    // divideSuccessRun2[1] = H1D_jetPt_ratio_run2[1]->Divide(H1D_jetPt_unfolded);


    // comparison with run2 results by rebinning Run3 into Run2 bins
    H1D_jetPt_unfolded_run2Comp[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp[1] = (TH1D*)H1D_jetPt_unfolded->Rebin(nBinPtJetsGen_run2[iRadius],"H1D_jetPt_unfolded_run2Comp_run2Rebin"+partialUniqueSpecifier, ptBinsJetsGen_run2[iRadius]);
    
    int scalingFactorRebin[9] = {2, 2, 2, 2, 2, 3, 3, 1, 1}; //width of run 2 histogram
    for (auto i = 1; i <= H1D_jetPt_unfolded_run2Comp[1]->GetNbinsX(); i++) {
      H1D_jetPt_unfolded_run2Comp[1]->SetBinContent(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[1]->GetBinContent(i));
      H1D_jetPt_unfolded_run2Comp[1]->SetBinError(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[1]->GetBinError(i));
    }
    H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_unfolded_run2Comp_HannaBossiLauraFile"+partialUniqueSpecifier);
    H1D_jetPt_ratio_run2[0] = (TH1D*)H1D_jetPt_unfolded_run2Comp[1]->Clone("H1D_jetPt_ratio_run2"+partialUniqueSpecifier);
    // H1D_jetPt_ratio_run2[1] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_ratio_run2_HannaBossiLauraFile"+partialUniqueSpecifier);
    divideSuccessRun2[0] = H1D_jetPt_ratio_run2[0]->Divide(H1D_jetPt_run2_MLPaperFile);
    // divideSuccessRun2[1] = H1D_jetPt_ratio_run2[1]->Divide(H1D_jetPt_unfolded);

    H1D_jetPt_unfolded_run2Comp_shapeComp[0] = (TH1D*)H1D_jetPt_unfolded_run2Comp[1]->Clone("H1D_jetPt_unfolded_run2Comp_run3rebinned_shapeComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_shapeComp[1] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_run2rescaled_shapeComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_shapeComp[1]->Scale(H1D_jetPt_unfolded_run2Comp_shapeComp[0]->GetBinContent(1)/H1D_jetPt_unfolded_run2Comp_shapeComp[1]->GetBinContent(1));

    H1D_jetPt_ratio_run2_shapeComp[0] = (TH1D*)H1D_jetPt_unfolded_run2Comp_shapeComp[0]->Clone("H1D_jetPt_unfolded_run2Comp_run2rescaled_shapeComp_ratio"+partialUniqueSpecifier);
    divideSuccessRun2_shapeComp = H1D_jetPt_ratio_run2_shapeComp[0]->Divide(H1D_jetPt_unfolded_run2Comp_shapeComp[1]);

    ///////////////////////////////////////
    // xT comparison: xT=2pT/sqrt(s) //////
    ///////////////////////////////////////

    double sqrtS_run2 = 5020; // same unit as pT
    double sqrtS_run3 = 5360; // same unit as pT

    // H1D_jetPt_unfolded_run2Comp_xT[0] = H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_xT_run2"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_run2Comp_xT[1] = H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp_xT_run3"+partialUniqueSpecifier);

    for (int iBin = 1; iBin <= H1D_jetPt_run2_MLPaperFile->GetNbinsX()+1; iBin++) {
      xtBinningVectorRun2.push_back(2./sqrtS_run2*H1D_jetPt_run2_MLPaperFile->GetXaxis()->GetBinLowEdge(iBin));
    }
    double* xtBinningRun2 = &xtBinningVectorRun2[0];
    for (int iBin = 1; iBin <= H1D_jetPt_unfolded->GetNbinsX()+1; iBin++) {
      xtBinningVectorRun3.push_back(2./sqrtS_run3*H1D_jetPt_unfolded->GetXaxis()->GetBinLowEdge(iBin));
    }
    double* xtBinningRun3 = &xtBinningVectorRun3[0];

    H1D_jetPt_unfolded_run2Comp_xT[0] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run2", "H1D_jetPt_unfolded_run2Comp_xT_run2", H1D_jetPt_run2_MLPaperFile->GetNbinsX(), xtBinningRun2);
    H1D_jetPt_unfolded_run2Comp_xT[1] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run3", "H1D_jetPt_unfolded_run2Comp_xT_run3", H1D_jetPt_unfolded->GetNbinsX(), xtBinningRun3);

    double dpt_dxt_run2=sqrtS_run2/2.;
    double dpt_dxt_run3=sqrtS_run3/2.;
    for (int iBin = 1; iBin <= H1D_jetPt_run2_MLPaperFile->GetNbinsX(); iBin++) {
      H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinContent(iBin, dpt_dxt_run2*H1D_jetPt_run2_MLPaperFile->GetBinContent(iBin));
      H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinError(iBin, dpt_dxt_run2*H1D_jetPt_run2_MLPaperFile->GetBinError(iBin));
    }
    for (int iBin = 1; iBin <= H1D_jetPt_unfolded->GetNbinsX(); iBin++) {
      H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinContent(iBin, dpt_dxt_run3*H1D_jetPt_unfolded->GetBinContent(iBin));
      H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinError(iBin, dpt_dxt_run3*H1D_jetPt_unfolded->GetBinError(iBin));
    }



    // ////////// using fit: commented for now as the fits aren't great //////////
    // // make xT binning
    // double binWidth = 0.001;
    // float maxPt = ptWindowDisplay[1];
    // float maxXt = 2*maxPt/sqrtS_run2; // sqrtS_run3 is larger than sqrtS_run2 so maxXtRun3 is smaller than maxXtRun2
    // for (int iBin = 1; iBin <= maxXt/binWidth; iBin++) { //pt bins of 1GeV are travelled
    //   xtBinningVector.push_back(binWidth * iBin);
    // }
    // int nBinsXt = xtBinningVector.size() - 1;
    // double* xtBinning = &xtBinningVector[0];

    // // get fit functions
    // std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tupleFitResult_run2 = TsallisFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange);
    // TF1* TF1_jetPt_run2_fit = std::get<0>(tupleFitResult_run2);
    // TMatrixDSym covMatrix_run2_fit = std::get<1>(tupleFitResult_run2);
    // double fitPtRange_run3[2] = {ptBinsJetsGen[iRadius][0], ptBinsJetsGen[iRadius][nBinPtJetsGen[iRadius]]};
    // std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tupleFitResult_run3 = TsallisFit(H1D_jetPt_unfolded, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange);
    // TF1* TF1_jetPt_run3_fit = std::get<0>(tupleFitResult_run3);
    // TMatrixDSym covMatrix_run3_fit = std::get<1>(tupleFitResult_run3);


    // // TF1_jetPt_run2_fit->SetParameters()

    // double parfitFunctionRun2[2];
    // double parfitFunctionRun3[2];
    // TF1_jetPt_run2_fit->GetParameters(&parfitFunctionRun2[0]);
    // TF1_jetPt_run3_fit->GetParameters(&parfitFunctionRun3[0]);

    // // TF1* initial function = new TF1("dNdptFunction_run2", "x*(1+1/([0]*[1])*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]);
    // // (FoG)'(x)=F'oG(x)*G'(x)
    // // [2] is dpt/dxt
    // // [3] is replacing x with xt=2pt/sqrtS
    // // TF1* dNdxtFunction_run2 = new TF1("dNdxtFunction_run2", "[2]*[3]*x*(1+1/([0]*[1])*[3]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]); 
    // TF1* dNdxtFunction_run2 = new TF1("dNdxtFunction_run2", "x*(1+1/([0]*[1])*[2]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]); //first [2]*[3] term, with [2] being dpt/dxt=sqrtS/2 and [3] being 2./sqrtS, cancel each other 
    // dNdxtFunction_run2->SetParameters(parfitFunctionRun2[0], parfitFunctionRun2[1], 2./sqrtS_run2); //transfor of pt -> xt in variable: xt=2*pt/sqrtS ; pt=sqrtS/2*xt ; dN/dxt(xt) = dN/dpt*dpt/dxt = dN/dpt(2*pt/sqrtS)*sqrtS/2
    // // TF1* dNdxtFunction_run3 = new TF1("dNdxtFunction_run3", "[2]*[3]*x*(1+1/([0]*[1])*[3]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]);
    // TF1* dNdxtFunction_run3 = new TF1("dNdxtFunction_run3", "x*(1+1/([0]*[1])*[2]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]);
    // dNdxtFunction_run3->SetParameters(parfitFunctionRun3[0], parfitFunctionRun3[1], 2./sqrtS_run2); //transfor of pt -> xt in variable: xt=2*pt/sqrtS ; pt=sqrtS/2*xt ; dN/dxt(xt) = dN/dpt*dpt/dxt = dN/dpt(2*pt/sqrtS)*sqrtS/2
    // // 2./sqrtS_run2


    // int nRowsNew = 3; // for parameters [0], [1], [2]
    // TMatrixDSym newCovMatrix_run2 = TMatrixDSym(nRowsNew);
    // TMatrixDSym newCovMatrix_run3 = TMatrixDSym(nRowsNew);
    // std::vector<double> initialisationVector;
    // for (int i = 0; i < nRowsNew*nRowsNew; i++) {
    //   initialisationVector.push_back(0);
    // }
    // newCovMatrix_run2.SetMatrixArray(&initialisationVector[0]);
    // newCovMatrix_run2.SetSub(0, 0, covMatrix_run2_fit); // inserts covMatrix_run2_fit as submatrix at row0 column0; other errors are 0 given [2] and[3] are constants
    // newCovMatrix_run3.SetMatrixArray(&initialisationVector[0]);
    // newCovMatrix_run3.SetSub(0, 0, covMatrix_run3_fit); // inserts covMatrix_run2_fit as submatrix at row0 column0; other errors are 0 given [2] and[3] are constants

    // Double_t *pData = newCovMatrix_run2.GetMatrixArray();
    // for (int i = 0; i < nRowsNew*nRowsNew; i++) {
    //   cout << "test i = " << i << ", newCovMatrix_run2[i]" << pData[i] << endl;
    // }

    // double xtRange[2] = {xtBinningVector.front(), xtBinningVector.back()};
    // TGraphErrors* fitFunctionTGraphErrors_run2 = getFunctionTGraphErrorsFromCovMatrix(xtRange, dNdxtFunction_run2, &newCovMatrix_run2, nBinsXt);    
    // TGraphErrors* fitFunctionTGraphErrors_run3 = getFunctionTGraphErrorsFromCovMatrix(xtRange, dNdxtFunction_run3, &newCovMatrix_run3, nBinsXt);    

    // // define xT histograms
    // H1D_jetPt_unfolded_run2Comp_xT[0] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run2", "H1D_jetPt_unfolded_run2Comp_xT_run2", nBinsXt, xtBinning);
    // H1D_jetPt_unfolded_run2Comp_xT[1] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run3", "H1D_jetPt_unfolded_run2Comp_xT_run3", nBinsXt, xtBinning);
    // H1D_jetPt_unfolded_run2Comp_xT[0]->Sumw2();
    // H1D_jetPt_unfolded_run2Comp_xT[1]->Sumw2();


    // // fill the histograms
    // double H1D_jetPt_unfolded_run2Comp_xT_errorRun2, H1D_jetPt_unfolded_run2Comp_xT_errorRun3;
    // double xtAtCenterOfBin;
    // for (int iBin = 1; iBin <= nBinsXt; iBin++) {
    //   xtAtCenterOfBin = (xtBinning[iBin-1]+xtBinning[iBin])/2;
    //   H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinContent(iBin, dNdxtFunction_run2->Eval(xtAtCenterOfBin));
    //   H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinError(iBin, fitFunctionTGraphErrors_run2->GetErrorY(iBin-1));
    //   H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinContent(iBin, dNdxtFunction_run3->Eval(xtAtCenterOfBin));
    //   H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinError(iBin, fitFunctionTGraphErrors_run3->GetErrorY(iBin-1));
    //   cout << "test fitFunctionTGraphErrors_run2->GetErrorY(iBin) = " << fitFunctionTGraphErrors_run2->GetErrorY(iBin-1) << endl;
    //   cout << "errors of xt graph look way too small (1E-8 or 1E-10)" << endl;
    // }

    // //ratio
    // H1D_jetPt_ratio_run2Comp_xT = (TH1D*)H1D_jetPt_unfolded_run2Comp_xT[1]->Clone("H1D_jetPt_ratio_run2Comp_xT"+partialUniqueSpecifier);
    // divideSuccessRun2_xt = H1D_jetPt_ratio_run2Comp_xT->Divide(H1D_jetPt_unfolded_run2Comp_xT[0]); // run3/run2




    // comparison with run2 results with fits only
    H1D_jetPt_unfolded_run2Comp_fits[0] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_fits_run2"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_fits[1] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp_fits_run3"+partialUniqueSpecifier);

    std::tuple<TF1*, TMatrixDSym, TFitResultPtr> exponentialFitWithLogTransfoResult_run3 = ExponentialFitWithLogTransfo(H1D_jetPt_unfolded, fitPtRange);
    TString histName_run3 = "rebinWithFit_run3"+partialUniqueSpecifier;
    std::pair<TH1D*, TGraphErrors*> pairResult_run3 = RebinWithFit(H1D_jetPt_unfolded, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange, histName_run3, exponentialFitWithLogTransfoResult_run3);

    TGraph_jetPt_unfolded_run2Comp_fits.push_back(pairResult_run2.second);
    TGraph_jetPt_unfolded_run2Comp_fits.push_back(pairResult_run3.second);
    

    double ptAtCenterOfBin;
    int nPoints = TGraph_jetPt_unfolded_run2Comp_fits.at(0)->GetN();
    TH1D* H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[2];
    H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0] = new TH1D("H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted_run2", "H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted_run2", nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]);
    H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1] = new TH1D("H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted_run3", "H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted_run3", nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]);
    H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0]->Sumw2();
    H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1]->Sumw2();
    for (int iBin = 1; iBin <= nBinPtJetsFine[iRadius]; iBin++) {
      ptAtCenterOfBin = (ptBinsJetsFine[iRadius][iBin-1]+ptBinsJetsFine[iRadius][iBin])/2;
      H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0]->SetBinContent(iBin, TGraph_jetPt_unfolded_run2Comp_fits.at(0)->GetPointY(ptAtCenterOfBin));
      H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0]->SetBinError(iBin, TGraph_jetPt_unfolded_run2Comp_fits.at(0)->GetErrorY(iBin-1));
      H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1]->SetBinContent(iBin, TGraph_jetPt_unfolded_run2Comp_fits.at(1)->GetPointY(ptAtCenterOfBin));
      H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1]->SetBinError(iBin, TGraph_jetPt_unfolded_run2Comp_fits.at(1)->GetErrorY(iBin-1));
      // cout << "run2 at bin "<< iBin << ":" << H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0]->GetBinContent(iBin) << endl;
      // cout << "run3 at bin "<< iBin << ":" << H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1]->GetBinContent(iBin) << endl;
    }

    //ratio
    H1D_jetPt_ratio_run2Comp_fits = (TH1D*)H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1]->Clone("H1D_jetPt_ratio_run2Comp_fits"+partialUniqueSpecifier);
    divideSuccessRun2_fits = H1D_jetPt_ratio_run2Comp_fits->Divide(H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0]); // run3/run2
  }

  cout << "comparison with refolded" << endl; 
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  }
  Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);
  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod(H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod, measuredInput, iDataset, iRadius, unfoldParameterInput, options);
  H1D_jetPt_unfolded_refoldedComp[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[1] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod->Clone("H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[2] = (TH1D*)H1D_jetPt_measured->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measuredRefolded[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measuredRefolded[1] = (TH1D*)H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod->Clone("H1D_jetPt_ratio_refoldedComp_RooUnfoldMethod"+partialUniqueSpecifier);
  // divideSuccessMeasuredRefolded[0] = H1D_jetPt_ratio_measuredRefolded[0]->Divide(H1D_jetPt_measured);
  // divideSuccessMeasuredRefolded[1] = H1D_jetPt_ratio_measuredRefolded[1]->Divide(H1D_jetPt_measured);
  divideSuccessMeasuredRefolded[0] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measuredRefolded[0], H1D_jetPt_measured);
  divideSuccessMeasuredRefolded[1] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measuredRefolded[1], H1D_jetPt_measured);

  // Cross section 
  if (doClosure_splitMC_mcdUnfoldedVsGen) {
    std::string detRespOption = useFactorisedMatrixInMcdUnfoldedClosure ? ", useIdentityForFluctResp" : "";
    unfoldParameter_mcSplitClosure = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded_inputSplitClosure, measuredInput_mcSplitInput, iDataset, iRadius, unfoldParameterInput, options+(std::string)detRespOption, splitTestControlMC).first;

    H1D_jetPt_unfolded_mcdSplitClosure[0] = (TH1D*)H1D_jetPt_mcp_mcSplitInput->Clone("H1D_jetPt_mcp_mcSplitClosure"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcdSplitClosure[1] = (TH1D*)H1D_jetPt_unfolded_inputSplitClosure->Clone("H1D_jetPt_unfolded_inputSplitClosure"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcdSplitClosure = (TH1D*)H1D_jetPt_unfolded_inputSplitClosure->Clone("H1D_jetPt_ratio_mcSplitClosure"+partialUniqueSpecifier);
    // divideSuccessMcdSplitClosure = H1D_jetPt_ratio_mcdSplitClosure->Divide(H1D_jetPt_mcp_mcSplitInput);
    divideSuccessMcdSplitClosure = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_mcdSplitClosure, H1D_jetPt_mcp_mcSplitInput);

  }

  if (doClosure_splitMC_mcpFoldedWithFluct) {
    cout << "comparison mcp folded with fluctuations vs mcp" << endl; 
    Get_Pt_spectrum_mcpFoldedWithFluctuations(H1D_jetPt_mcpFolded, iDataset, iRadius, options); // 
    Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_mcSplitInput_recBinning, iDataset, iRadius, options);

    H1D_jetPt_unfolded_mcpFoldedComp[0] = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_mcpFoldedComp_mcpFolded"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpFoldedComp[1] = (TH1D*)H1D_jetPt_mcp_mcSplitInput_recBinning->Clone("H1D_jetPt_mcpFoldedComp_mcp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcpFoldedMcp = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_ratio_mcpFoldedMcp"+partialUniqueSpecifier);
    // divideSuccessMcpFoldedMcp = H1D_jetPt_ratio_mcpFoldedMcp->Divide(H1D_jetPt_mcp_mcSplitInput_recBinning);
    divideSuccessMcpFoldedMcp = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_mcpFoldedMcp,H1D_jetPt_mcp_mcSplitInput_recBinning);

    
    // cout << "Integral mcp folded: " << H1D_jetPt_mcpFolded->Integral(1, H1D_jetPt_mcpFolded->GetNbinsX()) << endl;
    // cout << "Integral mcp       : " << H1D_jetPt_mcp_mcSplitInput->Integral(1, H1D_jetPt_mcp_mcSplitInput->GetNbinsX()) << endl;
  
    cout << "comparison mcp folded with fluctuations then unfolded vs mcp" << endl; 
    if (!normGenAndMeasByNEvtsForUnfoldingInput) {
      Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    } else{
      Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    }  
    Get_Pt_spectrum_unfolded(H1D_jetPt_mcpFoldedThenUnfolded, H1D_jetPt_mcpFolded2, iDataset, iRadius, unfoldParameterInput, options+", noKineEff, noPurity, noEff, inputIsMC, inputIsMCPFoldedTest, useIdentityForDetResp"); // input is mcp with fluctuations smearing: there are no fake jets, and the ptbinrange is the gen one so no kine efficiency
    H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[0] = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcpFoldedUnfolded"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcpFoldedUnfoldedMcp = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_ratio_mcpFoldedUnfoldedMcp"+partialUniqueSpecifier);
    // divideSuccessMcpFoldedUnfoldedMcp = H1D_jetPt_ratio_mcpFoldedUnfoldedMcp->Divide(H1D_jetPt_mcp); // divided by mcp Get_Pt_spectrum_mcp_genBinning
    divideSuccessMcpFoldedUnfoldedMcp = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_mcpFoldedUnfoldedMcp, H1D_jetPt_mcp); // divided by mcp Get_Pt_spectrum_mcp_genBinning
  }

  TString unfoldingCode;
  if (useManualRespMatrixSettingMethod){
    unfoldingCode = "myUnfold";
  } else {
    unfoldingCode = "joUnfold";
  }
  TString unfoldingInfo = (TString)unfoldingMethod+"-k="+Form("%i", unfoldParameter)+"-"+(TString)unfoldingPrior+"-"+unfoldingCode+"-matrixTransfo"+matrixTransformationOrder;

  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/IterationsDump", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/IterationsDump", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/IterationsDump", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/IterationsDump", &st1) == -1) {
  //     mkdir("pdfFolder/IterationsDump", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/IterationsDump", &st2) == -1) {
  //     mkdir("pngFolder/IterationsDump", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString dummyLegend[1] = {""};

  TString* yAxisLabel = texJet_dNdeta;
  TString* yAxisLabelXt = texJet_dNdeta;
  if (normaliseDistribsInComparisonPlots || normaliseUnfoldingResultsAtEnd) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
    yAxisLabelXt = texJet_d2Ndxtdeta_EventNorm;
  }

  TString pdfTitleBase = (TString)"IterationsDump/jet_"+unfoldingInfo;//+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_";
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  // comparison with measured
  TString unfoldedMeasuredCompLegend[2] = {"measured (gen binning)", "unfolded data"};
  TString* pdfName_measuredComp = new TString(pdfTitleBase+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldedMeasuredCompLegend, 2, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasured){
    TString* pdfName_ratio_measured = new TString(pdfTitleBase+"_measuredComp_ratio");
    Draw_TH1_Histogram(H1D_jetPt_ratio_measured, textContext, pdfName_ratio_measured, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    // TString* pdfName_ratio_measured_zoom = new TString(pdfTitleBase+"_measuredComp_ratio_zoom");
    // Draw_TH1_Histogram(H1D_jetPt_ratio_measured, textContext, pdfName_ratio_measured_zoom, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine,zoomToOneMedium2");
  }

  // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"mcp truth", "unfolded data"};
  TString* pdfName_mcpComp = new TString(pdfTitleBase+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldedTruthCompLegend, 2, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMcp){
    TString* pdfName_ratio_mcp = new TString(pdfTitleBase+"_mcpComp_ratio");
    Draw_TH1_Histogram(H1D_jetPt_ratio_mcp, textContext, pdfName_ratio_mcp, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_mcp_zoom = new TString(pdfTitleBase+"_mcpComp_ratio_zoom");
    Draw_TH1_Histogram(H1D_jetPt_ratio_mcp, textContext, pdfName_ratio_mcp_zoom, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with refolded
  TString unfoldedRefoldedCompLegend[3] = {"refolded manually", "refolded roounfold (noErrors)", "measured"};
  TString* pdfName_refoldedComp = new TString(pdfTitleBase+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldedRefoldedCompLegend, 3, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasuredRefolded[0] && divideSuccessMeasuredRefolded[1]) {
    TString* pdfName_ratio_refoldedComp = new TString(pdfTitleBase+"_RefoldedComp_ratio");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldedRefoldedCompLegend, 2, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_refoldedComp_zoom = new TString(pdfTitleBase+"_RefoldedComp_ratio_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldedRefoldedCompLegend, 2, textContext, pdfName_ratio_refoldedComp_zoom, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with Run 2
  if (comparePbPbWithRun2) {
    TString unfoldedRun2CompLegend_fitRebin[3] = {"unfolded Run3", "unfolded Run2 ML rebinned", "unfolded Run2 ML initial"};
    TString* pdfName_run2Comp_fitRebin = new TString(pdfTitleBase+"_run2Comp_fitRebin");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_fitRebin, unfoldedRun2CompLegend_fitRebin, 3, textContext, pdfName_run2Comp_fitRebin, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy,fitSingle", TGraph_jetPt_run2_MLPaperFile_fit);
    if (divideSuccessRun2_fitRebin[0] || divideSuccessRun2_fitRebin[1]) {
      TString* pdfName_ratio_run2_fitRebin = new TString(pdfTitleBase+"_run2Comp_fitRebin_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2_fitRebin[0], textContext, pdfName_ratio_run2_fitRebin, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    }

    TString unfoldedRun2CompLegend[3] = {"unfolded Run3", "unfolded Run3 rebinned", "unfolded Run2"};
    TString* pdfName_run2Comp = new TString(pdfTitleBase+"_run2Comp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp, unfoldedRun2CompLegend, 3, textContext, pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessRun2[0] || divideSuccessRun2[1]) {
      TString* pdfName_ratio_run2 = new TString(pdfTitleBase+"_run2Comp_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2[0], textContext, pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    }

    TString unfoldedRun2CompLegend_shapeComp[2] = {"unfolded Run3 rebinned", "unfolded Run2 scaled up to Run 3"};
    TString* pdfName_run2Comp_shapeComp = new TString(pdfTitleBase+"_run2Comp_shapeComp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_shapeComp, unfoldedRun2CompLegend_shapeComp, 2, textContext, pdfName_run2Comp_shapeComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessRun2[0] || divideSuccessRun2[1]) {
      TString* pdfName_ratio_shapeComp_run2 = new TString(pdfTitleBase+"_run2Comp_shapeComp_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2_shapeComp[0], textContext, pdfName_ratio_shapeComp_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    }

    TString unfoldedRun2CompLegend_xtComp[2] = {"Run2", "Run3"};
    TString* pdfName_run2Comp_xtComp = new TString(pdfTitleBase+"_run2Comp_xtComp");
    std::array<std::array<float, 2>, 2> drawnWindowXt = {{{(float)xtBinningVectorRun3.front(), (float)xtBinningVectorRun3.back()}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_xT, unfoldedRun2CompLegend_xtComp, 2, textContext, pdfName_run2Comp_xtComp, texXtX, yAxisLabelXt, texCollisionDataInfo, drawnWindowXt, legendPlacementAuto, contextPlacementAuto, "logy");
    // if (divideSuccessRun2_xt) {
    //   TString* pdfName_ratio_xtComp_run2 = new TString(pdfTitleBase+"_run2Comp_xtComp_ratio");
    //   Draw_TH1_Histogram(H1D_jetPt_ratio_run2Comp_xT, textContext, pdfName_ratio_xtComp_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowXt, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    // }

    TString unfoldedRun2CompLegend_fits[3] = {"unfolded Run2 ML", "unfolded Run3"};
    TString* pdfName_run2Comp_fits = new TString(pdfTitleBase+"_run2Comp_fits");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_fits, unfoldedRun2CompLegend_fits, 2, textContext, pdfName_run2Comp_fits, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy,fitCollection", TGraph_jetPt_unfolded_run2Comp_fits);
    if (divideSuccessRun2_fits) {
      TString* pdfName_ratio_run2_fits = new TString(pdfTitleBase+"_run2Comp_fits_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2Comp_fits, textContext, pdfName_ratio_run2_fits, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneExtraExtra, ratioLine");
    }
  }

  if (doClosure_splitMC_mcpFoldedWithFluct) {
    // comparison mcp folded with fluctuations vs mcp
    TString unfoldedMcpFoldedCheckLegend[2] = {"mcp-folded", "mcp"};
    TString* pdfName_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedComp, unfoldedMcpFoldedCheckLegend, 2, textContext, pdfName_McpFoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (H1D_jetPt_ratio_mcpFoldedMcp) {
      TString* pdfName_ratio_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedMcp, textContext, pdfName_ratio_McpFoldedCheck, texPtX, texRatioMcpFoldedVsMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }

    // comparison mcp folded with fluctuations then unfolded vs mcp
    TString unfoldedMcpFoldedUnfoldedCheckLegend[2] = {"mcp-folded unfolded", "mcp"};
    TString* pdfName_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedUnfoldedComp, unfoldedMcpFoldedUnfoldedCheckLegend, 2, textContext, pdfName_McpFoldedUnfoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMcpFoldedUnfoldedMcp) {
      TString* pdfName_ratio_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedUnfoldedMcp, textContext, pdfName_ratio_McpFoldedUnfoldedCheck, texPtX, texRatioMcpFoldedUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }
  }

  if (doClosure_splitMC_mcdUnfoldedVsGen) {
    // comparison mcd from controlMC unfolded vs mcp from controlMC
    TString unfoldedMcdClosureCheckLegend[2] = {"mcp", "mcd unfolded"};
    TString* pdfName_McdUnfoldedClosureCheck = new TString(pdfTitleBase+"_McdUnfoldedClosureCheck");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_mcdSplitClosure, unfoldedMcdClosureCheckLegend, 2, textContext, pdfName_McdUnfoldedClosureCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMcdSplitClosure) {
      TString* pdfName_ratio_McdUnfoldedClosureCheck = new TString(pdfTitleBase+"_McdUnfoldedClosureCheck_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcdSplitClosure, textContext, pdfName_ratio_McdUnfoldedClosureCheck, texPtX, texRatioMcdUnfoldedMcd, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }
  }
}
// declare this function here. It is used in the funciton below "Draw_Pt_spectrum_unfolded_parameterVariation_singleDataset"
void DrawRatioWithOffset(TH1D* histList[], int nUnfoldIteration,   const TString& yAxisTitle,const TString& canvasName, int unfoldIterationMax, int step, double yMin, double yMax);


void Draw_Pt_spectrum_unfolded_parameterVariation_singleDataset(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {

  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin)/step) + 1;

  TH1D* H1D_jetPt_measured;
  TH1D* H1D_jetPt_measured_genBinning;
  TH1D* H1D_jetPt_unfolded[nUnfoldIteration];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nUnfoldIteration];
  TH1D* H1D_jetPt_unfolded_mcpComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_measuredComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_ratio_mcp[nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measured[nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nUnfoldIteration];

  bool divideSuccessMcp[nUnfoldIteration];
  bool divideSuccessMeasured[nUnfoldIteration];
  bool divideSuccessMeasuredRefolded[nUnfoldIteration];
  TString partialUniqueSpecifier;


  if (!useFineBinningTest) {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
  }
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);

  int unfoldParameterInput = 0;

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    }
  }

  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){

    // unfoldParameterInput = iUnfoldIteration * step + unfoldIterationMin; 
    unfoldParameterInput = unfoldIterationMax - iUnfoldIteration * step; 
    // unfoldingIterationLegend[iUnfoldIteration] = unfoldParameterInput;

    cout << "((((((((((((()))))))))))))" << endl;
    cout << "Iteration "<< iUnfoldIteration << endl;
    cout << "((((((((((((()))))))))))))" << endl;
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iUnfoldIteration], measuredInput, iDataset, iRadius, unfoldParameterInput, options);

    // comparison with measured
    H1D_jetPt_unfolded_measuredComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measured[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
    divideSuccessMeasured[iUnfoldIteration] = H1D_jetPt_ratio_measured[iUnfoldIteration]->Divide(H1D_jetPt_measured_genBinning);

    // comparison with mcp truth
    H1D_jetPt_unfolded_mcpComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
    divideSuccessMcp[iUnfoldIteration] = H1D_jetPt_ratio_mcp[iUnfoldIteration]->Divide(H1D_jetPt_mcp);


    // comparison with refolded
    Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration], measuredInput, iDataset, iRadius, unfoldParameterInput, options);
    H1D_jetPt_unfolded_refoldedComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
    divideSuccessMeasuredRefolded[iUnfoldIteration] = H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration]->Divide(H1D_jetPt_measured);
  }
  H1D_jetPt_unfolded_measuredComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_mcpComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_measured->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);


  TString unfoldingInfo = (TString)unfoldingMethod+"_"+(TString)unfoldingPrior+"_kmax="+Form("%i", unfoldIterationMax);

  TString unfoldingIterationLegend[nUnfoldIteration+1]; IterationLegend(unfoldingIterationLegend, unfoldIterationMin, unfoldIterationMax, step);

  
  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo);

  TString* yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots || normaliseUnfoldingResultsAtEnd) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  Draw_TH1_Histograms(H1D_jetPt_unfolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");


    // comparison with measured
  // TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured (gen binning)"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
  TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMeasured_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasured[iUnfoldIteration])) {
      divideSuccessMeasured_boolsum = false;
    }
  }
  if (divideSuccessMeasured_boolsum){
    TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_measured, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
  }


    // comparison with mcp truth
  // TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"mcp";
  TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMcp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMcp[iUnfoldIteration])) {
      divideSuccessMcp_boolsum = false;
    }
  }
  if (divideSuccessMcp_boolsum){
    TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMcp");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_mcp, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

  // comparison with refolded
  // TString unfoldedRefoldedCompLegend[2] = {"refolded", "measured"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
  for(int iIteration = 0; iIteration < nUnfoldIteration; iIteration++){
    unfoldingIterationLegend[iIteration] += (TString)" refolded";
  }
  TString* pdfName_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessRefoldedComp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasuredRefolded[iUnfoldIteration])) {
      divideSuccessRefoldedComp_boolsum = false;
    }
  }
  if (divideSuccessRefoldedComp_boolsum){
    TString* pdfName_ratio_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_refoldedComp_zoom = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp_zoom, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////ratio refolded/measured with offset differetn k////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  DrawRatioWithOffset(H1D_jetPt_ratio_measuredRefolded, nUnfoldIteration, "Refolded / Measured", "ratio_RefoldMeasure_different_k_withOffset", unfoldIterationMax, step, 0.5, 1.7);
  DrawRatioWithOffset(H1D_jetPt_ratio_mcp, nUnfoldIteration, "unfolded / mcp", "ratio_Unfolded_mcp_different_k_withOffset", unfoldIterationMax, step, 0.7, 1.4); // last two numbers are y min and max
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

void Draw_Pt_spectrum_unfolded_datasetComparison(int iRadius, int unfoldParameterInput, std::string options) {
  bool splitTestControlMC = true;

  TH1D* H1D_jetPt_measured[nDatasets];
  TH1D* H1D_jetPt_measured_genBinning[nDatasets];
  TH1D* H1D_jetPt_unfolded[nDatasets];
  TH1D* H1D_jetPt_unfolded_ratio_datasets[nDatasets];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nDatasets];
  // TH1D* H1D_jetPt_mcpFolded[nDatasets];
  // TH1D* H1D_jetPt_mcpFolded2[nDatasets];
  // TH1D* H1D_jetPt_mcpFoldedThenUnfolded[nDatasets];
  // TH1D* H1D_jetPt_unfolded_mcpComp[2][nDatasets];
  // TH1D* H1D_jetPt_unfolded_measuredComp[2][nDatasets];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nDatasets];
  // TH1D* H1D_jetPt_unfolded_mcpFoldedComp[2][nDatasets];
  // TH1D* H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[2][nDatasets];
  TH1D* H1D_jetPt_mcp[nDatasets];
  TH1D* H1D_jetPt_mcp_recBinControl[nDatasets];
  TH1D* H1D_jetPt_ratio_mcp[nDatasets];
  TH1D* H1D_jetPt_ratio_measured[nDatasets];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nDatasets];

  TH1D* H1D_jetPt_unfolded_run2Comp[nDatasets+1];
  TH1D* H1D_jetPt_ratio_run2[nDatasets];
  
  // TH1D* H1D_jetPt_ratio_mcpFoldedMcp[nDatasets];
  // TH1D* H1D_jetPt_ratio_mcpFoldedUnfoldedMcp[nDatasets];


  bool divideSuccessDatasets[nDatasets];
  bool divideSuccessMcp[nDatasets];
  bool divideSuccessMeasured[nDatasets];
  bool divideSuccessMeasuredRefolded[nDatasets];
  bool divideSuccessRun2[nDatasets];
  // bool divideSuccessMcpFoldedMcp[nDatasets];
  // bool divideSuccessMcpFoldedUnfoldedMcp[nDatasets];

  TString partialUniqueSpecifier = (TString)"datasetComparison_R="+Form("%.1f",arrayRadius[iRadius]);
  TString datasetNameSpecifier[nDatasets];

  int unfoldParameter[nDatasets];


  // RUN 2 settings
  TH1D* H1D_jetPt_run2_MLPaperFile = new TH1D("H1D_jetPt_run2_MLPaperFile", "H1D_jetPt_run2_MLPaperFile", nBinPtJetsGen_run2[iRadius], ptBinsJetsGen_run2[iRadius]);
  if (comparePbPbWithRun2) {
    TGraph* Graph_jetPt_run2_MLPaperFile;
    double Ncoll;
    if (centralityRange[0] == 00 && centralityRange[1] == 10) {
      // Ncoll = (1780.9+1387.0)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (1956+1722+1521+1346)/4; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else if (centralityRange[0] == 50 && centralityRange[1] == 70) {
      // Ncoll = (103.7+46.1)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (89.8+39.8)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else {
      cout << "comparison with run2: Ncoll hasn't been calculated for this centrality interval" << endl;
    }
    double sigmaNN = 67.6; // value for sqrt(s) = 5.02 TeV https://arxiv.org/abs/1710.07098
    double T_AA = Ncoll / sigmaNN;
    Graph_jetPt_run2_MLPaperFile = ((TGraph*)((TDirectoryFile*)file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObjectAny("Graph1D_y1")); // https://doi.org/10.1016/j.physletb.2023.138412
    // H1D_jetPt_run2_MLPaperFile = (TH1D*)((TH1D*)(file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObject("Graph1D_y1"))->Clone("H1D_jetPt_run2_MLPaperFile");
    int Ngraph = Graph_jetPt_run2_MLPaperFile->GetN();
    for (int i=0; i < Ngraph; ++i) // setting bin contents to the TGraph values
    {
      double x,y;
      Graph_jetPt_run2_MLPaperFile->GetPoint(i, x, y);
      H1D_jetPt_run2_MLPaperFile->Fill(x, y); // uncertainties are of course screwed up
      int iHist = H1D_jetPt_run2_MLPaperFile->GetXaxis()->FindBin(x);
      H1D_jetPt_run2_MLPaperFile->SetBinError(iHist, Graph_jetPt_run2_MLPaperFile->GetErrorY(i));
    }
    H1D_jetPt_run2_MLPaperFile->Scale(T_AA);
    // now the spectre from the file is 1/N d2N/dpTdeta, instead of 1/T_AA 1/N d2N/dpTdeta
    H1D_jetPt_unfolded_run2Comp[nDatasets] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_run2Rebin_Run2"+partialUniqueSpecifier);
  }

  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    // getting inputs to unfolding
    if (!useFineBinningTest) {
      Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp[iDataset], iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_recBinControl[iDataset], iDataset, iRadius, options, splitTestControlMC);
      }
    } else {
      Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp[iDataset], iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp_recBinControl[iDataset], iDataset, iRadius, options, splitTestControlMC);
      }
    }

    TH1D* measuredInput[nDatasets];
    if (!normGenAndMeasByNEvtsForUnfoldingInput) {
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput[iDataset], iDataset, iRadius, options); 
      if (useFineBinningTest) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput[iDataset], iDataset, iRadius, options);
      }
    } else{
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput[iDataset], iDataset, iRadius, options);
      if (useFineBinningTest) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput[iDataset], iDataset, iRadius, options);
      }
    }

    // doing the unfolding
    unfoldParameter[iDataset] = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iDataset], measuredInput[iDataset], iDataset, iRadius, unfoldParameterInput, options).first;
    // TH1D* H1D_jetPt_unfolded2 = (TH1D*)H1D_jetPt_unfolded->Clone(H1D_jetPt_unfolded->GetName()+(TString)"H1D_jetPt_unfolded2");
    H1D_jetPt_unfolded_ratio_datasets[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_unfolded_ratio_datasets"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    divideSuccessDatasets[iDataset] = H1D_jetPt_unfolded_ratio_datasets[iDataset]->Divide(H1D_jetPt_unfolded[0]);
    // Creating to-be-plotted histograms
    datasetNameSpecifier[iDataset] = "_"+DatasetsNames[iDataset]+Form("%.1d",iDataset);


    cout << "comparison with measured" << endl; 
    if (!useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning[iDataset], iDataset, iRadius, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured_genBinning[iDataset], iDataset, iRadius, options);
    }
    H1D_jetPt_ratio_measured[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    // divideSuccessMeasured[iDataset] = H1D_jetPt_ratio_measured[iDataset]->Divide(H1D_jetPt_measured_genBinning[iDataset]);
    divideSuccessMeasured[iDataset] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measured[iDataset], H1D_jetPt_measured_genBinning[iDataset]);

    cout << "comparison with mcp truth" << endl; 
    H1D_jetPt_ratio_mcp[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    divideSuccessMcp[iDataset] = H1D_jetPt_ratio_mcp[iDataset]->Divide(H1D_jetPt_mcp[iDataset]);

    cout << "comparison with refolded" << endl; 
    if (!useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured[iDataset], iDataset, iRadius, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured[iDataset], iDataset, iRadius, options);
    }
    Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iDataset], measuredInput[iDataset], iDataset, iRadius, unfoldParameterInput, options);
    H1D_jetPt_unfolded_refoldedComp[iDataset] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iDataset]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    H1D_jetPt_ratio_measuredRefolded[iDataset] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iDataset]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    // divideSuccessMeasuredRefolded[iDataset] = H1D_jetPt_ratio_measuredRefolded[iDataset]->Divide(H1D_jetPt_measured[iDataset]);
    divideSuccessMeasuredRefolded[iDataset] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measuredRefolded[iDataset], H1D_jetPt_measured[iDataset]);

    cout << "comparison with run2" << endl; 
    if (comparePbPbWithRun2) {
      // // comparison with run2 results rebinned using a fit (errors look way underestimated; tsallis function not great aboe 100+ GeV; where should one eval the function inside a bin? probably not just the center)
      // H1D_jetPt_unfolded_run2Comp_fitRebin[0][iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_unfolded_run2Comp_fitRebin"+partialUniqueSpecifier);
      // double fitPtRange[2] = {ptBinsJetsGen_run2[iRadius][0], ptBinsJetsGen_run2[iRadius][nBinPtJetsGen_run2[iRadius]]};
      // std::pair<TH1D*, TF1*> pairResult = RebinWithTsallisFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange);
      // H1D_jetPt_run2_MLPaperFile_rebinned = pairResult.first;
      // TF1_jetPt_run2_MLPaperFile_fit[0] = pairResult.second;
      // H1D_jetPt_unfolded_run2Comp_fitRebin[1][iDataset] = (TH1D*)H1D_jetPt_run2_MLPaperFile_rebinned->Clone("H1D_jetPt_unfolded_run2_rebinned_fitRebin"+partialUniqueSpecifier);
      // H1D_jetPt_unfolded_run2Comp_fitRebin[2][iDataset] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2_fitRebin"+partialUniqueSpecifier);
      // // H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_unfolded_run2Comp_HannaBossiLauraFile"+partialUniqueSpecifier);


      // H1D_jetPt_ratio_run2_fitRebin[0][iDataset] = (TH1D*)H1D_jetPt_run2_MLPaperFile_rebinned->Clone("H1D_jetPt_ratio_run2_MLPaperFile_fitRebin"+partialUniqueSpecifier);
      // // H1D_jetPt_ratio_run2[1] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_ratio_run2_HannaBossiLauraFile"+partialUniqueSpecifier);
      // divideSuccessRun2_fitRebin[0][iDataset] = H1D_jetPt_ratio_run2_fitRebin[0][iDataset]->Divide(H1D_jetPt_unfolded[iDataset]);
      // // divideSuccessRun2[1] = H1D_jetPt_ratio_run2[1]->Divide(H1D_jetPt_unfolded);


      // comparison with run2 results by rebinning Run3 into Run2 bins
      H1D_jetPt_unfolded_run2Comp[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Rebin(nBinPtJetsGen_run2[iRadius],"H1D_jetPt_unfolded_run2Comp_run2Rebin"+partialUniqueSpecifier+datasetNameSpecifier[iDataset], ptBinsJetsGen_run2[iRadius]);
      int scalingFactorRebin[9] = {2, 2, 2, 2, 2, 3, 3, 1, 1}; //width of run 2 histogram
      for (auto i = 1; i <= H1D_jetPt_unfolded_run2Comp[iDataset]->GetNbinsX(); i++) {
        H1D_jetPt_unfolded_run2Comp[iDataset]->SetBinContent(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[iDataset]->GetBinContent(i));
        H1D_jetPt_unfolded_run2Comp[iDataset]->SetBinError(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[iDataset]->GetBinError(i));
      }
      H1D_jetPt_ratio_run2[iDataset] = (TH1D*)H1D_jetPt_unfolded_run2Comp[iDataset]->Clone("H1D_jetPt_ratio_run2_MLPaperFile"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
      divideSuccessRun2[iDataset] = H1D_jetPt_ratio_run2[iDataset]->Divide(H1D_jetPt_run2_MLPaperFile);
    }

    // if (doClosure_splitMC_mcpFoldedWithFluct) {
    //   cout << "comparison mcp folded with fluctuations vs mcp" << endl; 
    //   Get_Pt_spectrum_mcpFoldedWithFluctuations(H1D_jetPt_mcpFolded, iDataset, iRadius, options);
    //   H1D_jetPt_unfolded_mcpFoldedComp[0] = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_mcpFoldedComp_mcpFolded"+partialUniqueSpecifier);
    //   H1D_jetPt_unfolded_mcpFoldedComp[1] = (TH1D*)H1D_jetPt_mcp_recBinControl->Clone("H1D_jetPt_mcpFoldedComp_mcp"+partialUniqueSpecifier);
    //   H1D_jetPt_ratio_mcpFoldedMcp = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_ratio_mcpFoldedMcp"+partialUniqueSpecifier);
    //   divideSuccessMcpFoldedMcp = H1D_jetPt_ratio_mcpFoldedMcp->Divide(H1D_jetPt_mcp_recBinControl);

    //   // cout << "Integral mcp folded: " << H1D_jetPt_mcpFolded->Integral(1, H1D_jetPt_mcpFolded->GetNbinsX()) << endl;
    //   // cout << "Integral mcp       : " << H1D_jetPt_mcp_recBinControl->Integral(1, H1D_jetPt_mcp_recBinControl->GetNbinsX()) << endl;
    
    //   cout << "comparison mcp folded with fluctuations then unfolded vs mcp" << endl; 
    //   if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    //     Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    //   } else{
    //     Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    //   }  
    //   Get_Pt_spectrum_unfolded(H1D_jetPt_mcpFoldedThenUnfolded, H1D_jetPt_mcpFolded2, iDataset, iRadius, unfoldParameterInput, options+", noKineEff, noPurity, noEff, inputIsMC, inputIsMCPFoldedTest"); // input is mcp with fluctuations smearing: there are no fake jets, and the ptbinrange is the gen one so no kine efficiency
    //   H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[0] = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcpFoldedUnfolded"+partialUniqueSpecifier);
    //   H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcp"+partialUniqueSpecifier);
    //   H1D_jetPt_ratio_mcpFoldedUnfoldedMcp = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_ratio_mcpFoldedUnfoldedMcp"+partialUniqueSpecifier);
    //   divideSuccessMcpFoldedUnfoldedMcp = H1D_jetPt_ratio_mcpFoldedUnfoldedMcp->Divide(H1D_jetPt_mcp);
    // }
  }

  TString unfoldingCode;
  if (useManualRespMatrixSettingMethod){
    unfoldingCode = "myUnfold";
  } else {
    unfoldingCode = "joUnfold";
  }
  TString unfoldingInfo = (TString)unfoldingMethod+"-k="+Form("%i", unfoldParameterInput)+"-"+(TString)unfoldingPrior+"-"+unfoldingCode+"-matrixTransfo"+matrixTransformationOrder;

  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/IterationsDump", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/IterationsDump", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/IterationsDump", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/IterationsDump", &st1) == -1) {
  //     mkdir("pdfFolder/IterationsDump", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/IterationsDump", &st2) == -1) {
  //     mkdir("pngFolder/IterationsDump", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString dummyLegend[1] = {""};

  TString* yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots || normaliseUnfoldingResultsAtEnd) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }

  TString pdfTitleBase = (TString)"IterationsDump/jet_DatasetComp_"+unfoldingInfo;//+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_";
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  TString* pdfName_unfolded = new TString(pdfTitleBase+"_unfoldedSpectrum");
  Draw_TH1_Histograms(H1D_jetPt_unfolded, DatasetsNames, nDatasets, textContext, pdfName_unfolded, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (std::all_of(std::begin(divideSuccessDatasets), std::end(divideSuccessDatasets), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMeasured are true
    TString* pdfName_unfolded_ratio = new TString(pdfTitleBase+"_unfoldedSpectrum_ratio");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_ratio_datasets, DatasetsNames, nDatasets, textContext, pdfName_unfolded_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    // TString* pdfName_ratio_measured_zoom = new TString(pdfTitleBase+"_ratioToMeasured_zoom");
    // Draw_TH1_Histograms(H1D_jetPt_ratio_measured, DatasetsNames, nDatasets, textContext, pdfName_unfolded_ratio, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }

    // comparison with raw measured
  if (std::all_of(std::begin(divideSuccessMeasured), std::end(divideSuccessMeasured), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMeasured are true
    TString* pdfName_ratio_measured = new TString(pdfTitleBase+"_ratioToMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, DatasetsNames, nDatasets, textContext, pdfName_ratio_measured, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_measured_zoom = new TString(pdfTitleBase+"_ratioToMeasured_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, DatasetsNames, nDatasets, textContext, pdfName_ratio_measured_zoom, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }

    // comparison with mcp truth
  if (std::all_of(std::begin(divideSuccessMcp),std::end(divideSuccessMcp), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMcp are true
    TString* pdfName_ratio_mcp = new TString(pdfTitleBase+"_RatioToMcp");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, DatasetsNames, nDatasets, textContext, pdfName_ratio_mcp, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_mcp_zoom = new TString(pdfTitleBase+"_RatioToMcp_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, DatasetsNames, nDatasets, textContext, pdfName_ratio_mcp_zoom, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with refolded
  TString* pdfName_refoldedComp = new TString(pdfTitleBase+"_refoldedSpectrum");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, DatasetsNames, nDatasets, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (std::all_of(std::begin(divideSuccessMeasuredRefolded), std::end(divideSuccessMeasuredRefolded), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMeasuredRefolded are true
    TString* pdfName_ratio_refoldedComp = new TString(pdfTitleBase+"_refoldedSpetrum_ratioToMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, DatasetsNames, nDatasets, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_refoldedComp_zoom = new TString(pdfTitleBase+"_refoldedSpetrum_ratioToMeasured_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, DatasetsNames, nDatasets, textContext, pdfName_ratio_refoldedComp_zoom, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with run2
  TString* pdfName_run2Comp = new TString(pdfTitleBase+"_run2Comp");
  TString DatasetsNamesAppended[nDatasets+1];
  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    DatasetsNamesAppended[iDataset] = DatasetsNames[iDataset];
  }
  DatasetsNamesAppended[nDatasets] = "Run2 ML";

  Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp, DatasetsNamesAppended, nDatasets+1, textContext, pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (std::all_of(std::begin(divideSuccessRun2), std::end(divideSuccessRun2), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessRun2 are true
    TString* pdfName_ratio_run2 = new TString(pdfTitleBase+"_run2Comp_ratio");
    Draw_TH1_Histograms(H1D_jetPt_ratio_run2, DatasetsNames, nDatasets, textContext, pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
  }

  // if (doClosure_splitMC_mcpFoldedWithFluct) {
  //   // comparison mcp folded with fluctuations vs mcp
  //   TString unfoldedMcpFoldedCheckLegend[2] = {"mcp-folded", "mcp"};
  //   TString* pdfName_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp");
  //   Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedComp, unfoldedMcpFoldedCheckLegend, 2, textContext, pdfName_McpFoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  //   if (H1D_jetPt_ratio_mcpFoldedMcp) {
  //     TString* pdfName_ratio_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp_ratio");
  //     Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedMcp, textContext, pdfName_ratio_McpFoldedCheck, texPtX, texRatioMcpFoldedVsMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
  //   }

  //   // comparison mcp folded with fluctuations then unfolded vs mcp
  //   TString unfoldedMcpFoldedUnfoldedCheckLegend[2] = {"mcp-folded unfolded", "mcp"};
  //   TString* pdfName_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck");
  //   Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedUnfoldedComp, unfoldedMcpFoldedUnfoldedCheckLegend, 2, textContext, pdfName_McpFoldedUnfoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  //   if (divideSuccessMcpFoldedUnfoldedMcp) {
  //     TString* pdfName_ratio_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck_ratio");
  //     Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedUnfoldedMcp, textContext, pdfName_ratio_McpFoldedUnfoldedCheck, texPtX, texRatioMcpFoldedUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
  //   }
  // }
}

// void Draw_Pt_spectrum_unfolded_ImprovedStatErrors(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {

//   TH1D* measuredInput;
//   if (!normGenAndMeasByNEvtsForUnfoldingInput) {
//     Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
//     if (useFineBinningTest) {
//       Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
//     }
//   } else{
//     Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
//     if (useFineBinningTest) {
//       Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
//     }
//   }
//   TH1D* H1D_jetPt_unfolded_withImprovedErrors;
//   Get_Pt_spectrum_unfolded_ImprovedStatisticalErrors(H1D_jetPt_unfolded_withImprovedErrors, measuredInput, iDataset, iRadius, unfoldParameterInput, options);

//   TString pdfname = "jet_Pt_spectrum_unfolded_RelativeUncertainty_Dataset"+DatasetsNames[iDataset]+"_R"+Form("%.1f", arrayRadius[iRadius])+"_k"+Form("%i", unfoldParameterInput);
//   // TString textContext = "Unfolded improved errors";
//   TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
//   // Test : 
//   // TCanvas* c_relUnc = new TCanvas(pdfname, pdfname, 800, 800);
//   // H1D_jetPt_unfolded_withImprovedErrors->Draw();
//   // c_relUnc->SetLogy();

//   // error with Draw_TH1_Histogram!!!?
//   // Draw_TH1_Histogram(H1D_jetPt_unfolded_withImprovedErrors, textContext, pdfname, texPtJetRec, texJet_d2Ndptdeta_EventNorm, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

// }

void DrawRatioWithOffset(TH1D* histList[], int nUnfoldIteration, const TString& yAxisTitle,const TString& canvasName, int unfoldIterationMax, int step, double yMin, double yMax){
    // DrawRatioWithOffset(..., -1, -1);
    TString canvasNameFull = canvasName + "_" + unfoldingMethod;
    TCanvas* c = new TCanvas(canvasNameFull, canvasNameFull, 800, 800);
    TMultiGraph* mg = new TMultiGraph();

    int customColors[] = {
        kBlack, kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1
    };
    int nColors = sizeof(customColors) / sizeof(customColors[0]);

    //int nUnfoldIteration = histList.size();  // infer number of iterations from input list

    for (int i = 0; i < nUnfoldIteration; ++i) {
        TH1D* h = histList[i];
        int nBins = h->GetNbinsX();

        std::vector<double> x_vals, y_vals, ex_vals, ey_vals;

        for (int bin = 1; bin <= nBins; ++bin) {
            // Replace ptBinsJetsRec[iRadius] with your own binning logic if needed:
            double binLowEdge = h->GetXaxis()->GetBinLowEdge(bin);
            double binUpEdge  = h->GetXaxis()->GetBinUpEdge(bin);
            double binWidth   = binUpEdge - binLowEdge;
            double offset     = (i - nUnfoldIteration / 2.0) * 0.065 * binWidth;

            double x  = h->GetBinCenter(bin) + offset;
            double y  = h->GetBinContent(bin);
            double ex = 0;
            double ey = h->GetBinError(bin);

            x_vals.push_back(x);
            y_vals.push_back(y);
            ex_vals.push_back(ex);
            ey_vals.push_back(ey);
        }

        TGraphErrors* gr = new TGraphErrors(nBins, &x_vals[0], &y_vals[0], &ex_vals[0], &ey_vals[0]);
        int color = customColors[i % nColors];
        gr->SetMarkerStyle(20 + i);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetMarkerSize(0.8);
        gr->SetTitle(Form("k_{unfold} = %d", unfoldIterationMax-step*i));  // Adjust if you have different logic

        mg->Add(gr, "P");
    }

    mg->Draw("A");
    mg->GetXaxis()->SetLimits(10.0, 140.0);
    if (yMin < yMax) {
    mg->GetYaxis()->SetRangeUser(yMin, yMax);
    }
    mg->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    mg->GetYaxis()->SetTitle(yAxisTitle);

    c->BuildLegend();
    c->Update();

    // Draw horizontal reference line at y=1
    double xmin = mg->GetXaxis()->GetXmin();
    double xmax = mg->GetXaxis()->GetXmax();
    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);
    line->Draw("same");

    c->Update();
    // Auto-save
    c->SaveAs(canvasNameFull + ".pdf");
    c->SaveAs(canvasNameFull + ".png");
}

void MakeRatio(){
    // Open file in UPDATE mode (so we can write result)
    // TFile* JJ_Gap2 = TFile::Open("../20260319_Unf_LHC26a6_637817/output.root", "READ");
     // H1D_jetEfficiency_LHC25b4ab6_R02_Lead5_654885;1
  // H1D_fakeRatioLHC25b4ab6_R02_Lead5_654885;1
  // H1D_Pt_Unfolded_w_MB;1

  //JJ
  // H1D_jetEfficiency_LHC26c5_R02_Lead5_654919;1
  // H1D_fakeRatioLHC26c5_R02_Lead5_654919;1
  // H1D_Pt_Unfolded_w_JJ;1
    TFile* JJ_Gap3 = TFile::Open("output.root", "READ");

    TFile* MB = TFile::Open("../20260417_Unf_656693_usingMB_654885/output.root", "READ");
    
    // JJ_LHC26a6_train615296_Unf_ppref_unfolded //  MB_LHC25b4b5_train533385_Unf_ppref_unfolded
    // Retrieve histograms
    // TH1D* h1 = (TH1D*)JJ_Gap2->Get("JJ_Gap2_LHC26a6_637817_ppref_unfolded");
    TH1D* h2 = (TH1D*)JJ_Gap3->Get("H1D_Pt_Unfolded_w_JJ");
    TH1D* h1 = (TH1D*)MB->Get("H1D_Pt_Unfolded_w_MB");

    // Clone first histogram to store ratio
    // TH1D* hRatio_JJ2_MB = (TH1D*)h1->Clone("hRatio_JJ2_MB");
    TH1D* hRatio_JJ3_MB = (TH1D*)h2->Clone("hRatio_JJ3_MB");
    // hRatio_JJ2_MB->Divide(h1);
    // hRatio_JJ3_MB->Divide(h1);
    hRatio_JJ3_MB->Reset(); // Clear the content, keep the bins

    // 2. Loop through every bin of the new Ratio histogram
    for (int i = 1; i <= hRatio_JJ3_MB->GetNbinsX(); ++i) {
        double binCenter = hRatio_JJ3_MB->GetBinCenter(i);

        // 3. Apply your range constraint
        if (binCenter >= 5.0 && binCenter <= 100.0) {
            
            // Find the corresponding bin index in the other histogram (h1)
            int binH1 = h1->FindBin(binCenter);
            
            double val2 = h2->GetBinContent(i);
            double val1 = h1->GetBinContent(binH1);
            double err2 = h2->GetBinError(i);
            double err1 = h1->GetBinError(binH1);

            if (val1 > 0) {
                double ratio = val2 / val1;
                hRatio_JJ3_MB->SetBinContent(i, ratio);
                
                // Calculate error propagation: (R/V)^2 = (e1/v1)^2 + (e2/v2)^2
                double error = ratio * TMath::Sqrt(TMath::Power(err1/val1, 2) + TMath::Power(err2/val2, 2));
                hRatio_JJ3_MB->SetBinError(i, error);
            }
        } else {
            // Outside the range, set to 0 or leave empty
            hRatio_JJ3_MB->SetBinContent(i, 0);
            hRatio_JJ3_MB->SetBinError(i, 0);
        }
    }

    // TH1D* hRatio_JJ2_JJ3 = (TH1D*)h1->Clone("hRatio_JJ2_JJ3");
    // hRatio_JJ2_JJ3->Divide(h2);
    hRatio_JJ3_MB->Fit("pol0", "Q0");
    TF1* fit = hRatio_JJ3_MB->GetFunction("pol0");
    double value = fit->GetParameter(0);
    double error = fit->GetParError(0);

    TString* pdfName_JJ3_MB = new TString("Ratio_of_Data_649659_Unf_w_JJ_LHC26c5_649618_to_MB_LHC25b4ab6_649683");
    TString* YLabel_JJ3_MB = new TString("Unf with JJ_Gap3 / Unf with MB");
    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, "Ratio_of_Data_649659_Unf_w_JJ_LHC26c5_649618_to_MB_LHC25b4ab6_649683"));
    cout << Form("Fit result: %.3f #pm %.3f", value, error) << endl;
    Draw_TH1_Histogram(hRatio_JJ3_MB, textContext, pdfName_JJ3_MB, texPtJetRec, YLabel_JJ3_MB, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");

    // TString* pdfName = new TString("Ratio of JJ_Gap2 and JJ_Gap3 to MB");
    // TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, "Ratio_of_DataUnfold_w_JJ_Gap2_LHC26a6_637817_and_JJ_Gap3_LHC26c5_637150_to_MB_LHC25b4ab6_637087"));
    // TString* YLabel = new TString("Unf with JJ / Unf with MB");

    // TH1D* H1D_jetPt_Unfolded_ratio[2];
    // H1D_jetPt_Unfolded_ratio[0] = hRatio_JJ2_MB;
    // H1D_jetPt_Unfolded_ratio[1] = hRatio_JJ3_MB;

    // TString LegendNames[2];
    // LegendNames[0] = "JJ_Gap2 / MB";
    // LegendNames[1] = "JJ_Gap3 / MB";
    
    // Draw_TH1_Histograms(H1D_jetPt_Unfolded_ratio, LegendNames, 2, textContext, pdfName, texPtJetRec, YLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");

    // file->Close();

    std::cout << "Ratio successfully created and saved." << std::endl;
}


void Plot_Eff_Same_Convas(){
    // Open file in UPDATE mode (so we can write result)
    // TFile* JJ_Gap2 = TFile::Open("../20260319_Unf_LHC26a6_637817/output.root", "READ");

    TFile* JJ_Gap3 = TFile::Open("output.root", "READ");

    TFile* MB = TFile::Open("../20260417_Unf_656693_usingMB_654885/output.root", "READ");
    
    // JJ_LHC26a6_train615296_Unf_ppref_unfolded //  MB_LHC25b4b5_train533385_Unf_ppref_unfolded
    // Retrieve histograms
    // TH1D* h1 = (TH1D*)JJ_Gap2->Get("JJ_Gap2_LHC26a6_637817_ppref_unfolded");
    TH1D* h2_e = (TH1D*)JJ_Gap3->Get("H1D_jetEfficiency_LHC26c5_R02_Lead5_654919");
    TH1D* h1_e = (TH1D*)MB->Get("H1D_jetEfficiency_LHC25b4ab6_R02_Lead5_654885");

    TH1D* h2_p = (TH1D*)JJ_Gap3->Get("H1D_fakeRatioLHC26c5_R02_Lead5_654919");
    TH1D* h1_p = (TH1D*)MB->Get("H1D_fakeRatioLHC25b4ab6_R02_Lead5_654885");

    TCanvas* c_eff = new TCanvas("Matching_Efficiency_Comparison", "Matching_Efficiency_Comparison", 800, 800);
    h1_e->SetMarkerStyle(20);
    h1_e->SetMarkerColor(kRed);
    h1_e->SetLineColor(kRed);
    h1_e->SetMarkerSize(0.8);
    h1_e->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1_e->GetYaxis()->SetTitle("Matching Efficiency");
    h1_e->Draw("E1");
    h2_e->SetMarkerStyle(21);
    h2_e->SetMarkerColor(kBlue);
    h2_e->SetLineColor(kBlue);
    h2_e->SetMarkerSize(0.8);
    h2_e->Draw("E1 same");
    // --- Build the Legend ---
    // Coordinates: x1, y1, x2, y2 (normalized 0 to 1)
    TLegend* leg = new TLegend(0.15, 0.75, 0.45, 0.88); 
    leg->SetBorderSize(0); // Clean look without a box
    leg->SetFillStyle(0);   // Transparent background
    leg->SetTextSize(0.035);

    // Add entries: "p" means it shows the marker, "l" means it shows the line
    leg->AddEntry(h1_e, "MB ", "pl");
    leg->AddEntry(h2_e, "JJ Gap3", "pl");
    leg->Draw();

    TCanvas* c_pur = new TCanvas("Purity_Comparison", "Purity_Comparison", 800, 800);
    h1_p->SetMarkerStyle(20);
    h1_p->SetMarkerColor(kRed);
    h1_p->SetLineColor(kRed); 
    h1_p->SetMarkerSize(0.8);
    h1_p->SetTitle("Jet Purity Comparison");
    h1_p->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1_p->GetYaxis()->SetTitle("Fake Ratio");
    h1_p->Draw("E1");
    h2_p->SetMarkerStyle(21);
    h2_p->SetMarkerColor(kBlue);
    h2_p->SetLineColor(kBlue);
    h2_p->SetMarkerSize(0.8);
    h2_p->Draw("E1 same");
    
    // --- Build the Legend ---
    TLegend* leg_p = new TLegend(0.15, 0.75, 0.45, 0.88); 
    leg_p->SetBorderSize(0); // Clean look without a box
    leg_p->SetFillStyle(0);   // Transparent background
    leg_p->SetTextSize(0.035);  
    leg_p->AddEntry(h1_p, "MB ", "pl");
    leg_p->AddEntry(h2_p, "JJ Gap3", "pl");
    leg_p->Draw();
}

void Plot_ratio_Unf_Mcp_Same_Convas(){
    // Open file in UPDATE mode (so we can write result)
    // TFile* JJ_Gap2 = TFile::Open("../20260319_Unf_LHC26a6_637817/output.root", "READ");
    TFile* JJ_Gap3 = TFile::Open("output.root", "READ");

    TFile* MB = TFile::Open("../20260409_UnfMB_R02_Lead3_Full650972_649683/output.root", "READ");
    
    // JJ_LHC26a6_train615296_Unf_ppref_unfolded //  MB_LHC25b4b5_train533385_Unf_ppref_unfolded
    // Retrieve histograms
    // TH1D* h1 = (TH1D*)JJ_Gap2->Get("JJ_Gap2_LHC26a6_637817_ppref_unfolded");
    TH1D* h2_e = (TH1D*)JJ_Gap3->Get("H1D_jetPt_ratio_mcp_LHC24ap_pass1_R02_Lead3_650972_using_LHC26c5_R02_Lead3_649618_Unf");
    TH1D* h1_e = (TH1D*)MB->Get("H1D_jetPt_ratio_mcp_LHC24ap_pass1_R02_Lead3_650972_using_LHC25b4ab6_R02_Lead3_649683_Unf");

    TCanvas* c_eff = new TCanvas("Matching_Efficiency_Comparison", "Matching_Efficiency_Comparison", 800, 800);
    h1_e->SetMarkerStyle(20);
    h1_e->SetMarkerColor(kRed);
    h1_e->SetLineColor(kRed);
    h1_e->SetMarkerSize(1);
    h1_e->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1_e->GetYaxis()->SetTitle("unfolded(i) / mcp(i)");
    h1_e->Draw("E1");
    h2_e->SetMarkerStyle(21);
    h2_e->SetMarkerColor(kBlue);
    h2_e->SetLineColor(kBlue);
    h2_e->SetMarkerSize(1);
    h2_e->Draw("E1 same");
    // --- Build the Legend ---
    // Coordinates: x1, y1, x2, y2 (normalized 0 to 1)
    TLegend* leg = new TLegend(0.15, 0.75, 0.45, 0.88); 
    leg->SetBorderSize(0); // Clean look without a box
    leg->SetFillStyle(0);   // Transparent background
    leg->SetTextSize(0.035);

    // Add entries: "p" means it shows the marker, "l" means it shows the line
    leg->AddEntry(h1_e, "Unf w MB ", "pl");
    leg->AddEntry(h2_e, "Unf w JJ Gap3", "pl");
    leg->Draw();

    //Adjust window x from 4 to 140
    h1_e->GetXaxis()->SetRangeUser(4, 140);
    h2_e->GetXaxis()->SetRangeUser(4, 140);
    h1_e->GetYaxis()->SetRangeUser(0.7, 1.6);
    h2_e->GetYaxis()->SetRangeUser(0.7, 1.6);

    //Add ratio line at y=1 from x 4 to 140
    double xmin = h1_e->GetXaxis()->GetXmin();
    double xmax = h1_e->GetXaxis()->GetXmax();
    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);
    line->Draw("same"); 


    
}



void SaveSettingsToText(int iDataset, int iRadius, int unfoldParameter) {
    std::ofstream outFile("settings.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not create settings.txt" << std::endl;
        return;
    }

    if (doClosure_splitMC_mcdUnfoldedVsGen){
        outFile << " Unfolding Closure Test" << std::endl;
    }
    outFile << "##################################################" << std::endl;
    outFile << "###             UNFOLDING SETTINGS             ###" << std::endl;
    outFile << "##################################################" << std::endl;
    outFile << std::endl;

    // --- Dataset Settings ---
    outFile << "--- Dataset Settings ---" << std::endl;
    outFile << "iDataset       : " << iDataset << std::endl;
    for (int d = 0; d < nDatasets; ++d) {
        outFile << "Data  [" << d << "]      : " << Datasets[d] << std::endl;
        outFile << "MC    [" << d << "]      : " << MC_Datasets[d] << std::endl;
        outFile << "Derived[" << d << "]     : " << (derived_data[d] ? "true" : "false") << std::endl;
    }
    outFile << "mcIsWeighted   : " << (mcIsWeighted ? "true" : "false") << std::endl;
    outFile << std::endl;

    // --- Radius Settings ---
    outFile << "--- Radius Settings ---" << std::endl;
    outFile << "iRadius        : " << iRadius << std::endl;
    for (int r = 0; r < nRadius; ++r) {
        outFile << "Radius [" << r << "]     : " << RadiusLegend[r] << std::endl;
    }
    outFile << std::endl;

    // --- Unfolding Settings ---
    outFile << "--- Unfolding Settings ---" << std::endl;
    outFile << "unfoldingPrior  : " << unfoldingPrior << std::endl;
    outFile << "unfoldingMethod : " << unfoldingMethod << " ---- parameter : " << unfoldParameter << std::endl;
    outFile << std::endl;

    // --- Save Reconstruction Bins ---
    outFile << "--- ptBinsJetsRec ---" << std::endl;
    for (int r = 0; r < nRadius; ++r) {
        outFile << "Radius [" << r << "] (" << nBinPtJetsRec[r] << " bins): ";
        for (int i = 0; i <= nBinPtJetsRec[r]; ++i) {
            outFile << ptBinsJetsRec[r][i] << (i == nBinPtJetsRec[r] ? "" : ", ");
        }
        outFile << std::endl;
    }
    outFile << std::endl;
    outFile << "--------------------------------------------------" << std::endl;

    // --- Save Generated Bins ---
    outFile << "--- ptBinsJetsGen ---" << std::endl;
    for (int r = 0; r < nRadius; ++r) {
        outFile << "Radius [" << r << "] (" << nBinPtJetsGen[r] << " bins): ";
        for (int i = 0; i <= nBinPtJetsGen[r]; ++i) {
            outFile << ptBinsJetsGen[iRadius][i] << (i == nBinPtJetsGen[r] ? "" : ", ");
        }
        outFile << std::endl;
    }
    outFile << std::endl;

    outFile.close();
    std::cout << "Settings successfully saved to settings.txt" << std::endl;
}


double Get_zVertex_reconstruction_efficiency(int iDataset, int iRadius){

  // Retrieve histogram
  TH1D* H1D_collisions_zvertex = (TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData + "/h_collisions_zvertex");

  if (!H1D_collisions_zvertex) {
    std::cerr << "Error: z-vertex histogram not found!" << std::endl;
    return -1.;
  }

  H1D_collisions_zvertex = (TH1D*)H1D_collisions_zvertex->Clone(Form("hZvertex_clone_dataset%d_R%d", iDataset, iRadius));

  TCanvas* c = new TCanvas(Form("c_zvtx_%d_R%d", iDataset, iRadius),"Z-vertex fit", 800, 600);
  H1D_collisions_zvertex->SetLineWidth(3);   // thicker line
  H1D_collisions_zvertex->SetLineColor(kBlue);
  H1D_collisions_zvertex->SetTitle("Z-vertex distribution; z (cm); Counts");

  // Optional: also make markers nicer
  H1D_collisions_zvertex->SetMarkerStyle(20);
  H1D_collisions_zvertex->SetMarkerSize(0.8);

  TF1* fGaus = new TF1(Form("fGaus_%d_R%d", iDataset, iRadius), "gaus", -10, 10);
  fGaus->SetLineColorAlpha(kRed, 0.25);
  fGaus->SetLineWidth(15);

  H1D_collisions_zvertex->Fit(fGaus, "R");  // "R" = use fit range

  H1D_collisions_zvertex->Draw("");
  fGaus->Draw("same");

  // c->BuildLegend();

  // --- Extract parameters ---
  double A     = fGaus->GetParameter(0);
  double mean  = fGaus->GetParameter(1);
  double sigma = fGaus->GetParameter(2);

  // --- Integrals ---
  double totalIntegral   = A * sigma * std::sqrt(2 * TMath::Pi());
  double partialIntegral = fGaus->Integral(-10, 10);

  double efficiency = partialIntegral / totalIntegral;

  // Debug print (very useful)
  std::cout << "mean = " << mean
            << ", sigma = " << sigma
            << ", efficiency = " << efficiency << std::endl;

  return efficiency;
}

TH1D* Get_TVX_Eff(int iDataset, int iRadius){
  // Get 2D histogram
  TH2D* h2 = (TH2D*)file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get("jet-cross-section-efficiency/h2_jet_pt_part_eventselection");

  if (!h2) {
    std::cerr << "Error: 2D histogram not found!" << std::endl;
    return nullptr;
  }

  // Find Y bins for kTVX and NColl
  int bin_kTVX = -1;
  int bin_NColl = -1;

  for (int i = 1; i <= h2->GetYaxis()->GetNbins(); i++) {
    TString label = h2->GetYaxis()->GetBinLabel(i);

    if (label.Contains("kTVX"))  bin_kTVX = i;
    if (label.Contains("INEL")) bin_NColl = i;
  }

  if (bin_kTVX < 0 || bin_NColl < 0) {
    std::cerr << "Error: Could not find kTVX or NColl in Y axis labels!" << std::endl;
    return nullptr;
  }
  
  // Project to X (pt) for both selections
  TH1D* h_kTVX = h2->ProjectionX("h_kTVX", bin_kTVX, bin_kTVX);
  TH1D* h_NColl = h2->ProjectionX("h_NColl", bin_NColl, bin_NColl);

  TH1D* h_kTVX_reb = (TH1D*)h_kTVX->Rebin(nBinPtJetsGen[iRadius],Form("TVX_rebinned_dataset%d_R%d", iDataset, iRadius),ptBinsJetsGen[iRadius]);
  TH1D* h_NColl_reb = (TH1D*)h_NColl->Rebin(nBinPtJetsGen[iRadius],Form("INEL_rebinned_dataset%d_R%d", iDataset, iRadius),ptBinsJetsGen[iRadius]);
  h_kTVX_reb->Scale(1.,"width");
  h_NColl_reb->Scale(1.,"width");

  TCanvas* c = new TCanvas("c_TVX_check", "kTVX vs INEL", 800, 600);
  h_kTVX_reb->SetLineColor(kRed);
  h_kTVX_reb->SetLineWidth(2);
  h_kTVX_reb->SetTitle("Jet p_{T} spectra; p_{T}; Counts");
  h_NColl_reb->SetLineColor(kBlue);
  h_NColl_reb->SetLineWidth(2);
  h_kTVX_reb->Draw("hist");
  h_NColl_reb->Draw("hist same");
  TLegend* leg = new TLegend(0.60, 0.70, 0.85, 0.85);
  leg->SetBorderSize(0);      // Remove the black border box
  leg->SetFillStyle(0);       // Make the background transparent
  leg->SetTextFont(42);       // Standard ALICE font
  
  // Add entries: (histogram, "Label", "Option")
  // "l" means draw a line in the legend
  leg->AddEntry(h_NColl_reb, "INEL (Total)", "l");
  leg->AddEntry(h_kTVX_reb, "kTVX (Triggered)", "l");
  
  leg->Draw();
  c->SetLogy();  // VERY useful for spectra

  TH1D* TVX_eff = (TH1D*)h_kTVX_reb->Clone(Form("TVX_Eff_dataset%d_R%d", iDataset, iRadius));
  TVX_eff->Divide(h_kTVX_reb, h_NColl_reb, 1.0, 1.0, "B"); // the probability that an inelastic collision passes the TVX selection, as a function of pT

  TCanvas *cEff = new TCanvas("cEff", "TVX Efficiency", 800, 600);
  TVX_eff->GetYaxis()->SetTitle("#epsilon_{TVX}"); 
  TVX_eff->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV/#it{c})");
  TVX_eff->SetMarkerStyle(20);     // Filled circle
  TVX_eff->SetMarkerSize(1.2);
  TVX_eff->SetMarkerColor(kBlue+1);
  TVX_eff->SetLineColor(kBlue+1);
  TVX_eff->SetMinimum(0.0);        // Set Y-axis to start at 0
  TVX_eff->SetMaximum(1.2);        // Leave some room at the top
  TVX_eff->Draw("E1 P"); 
  cEff->Update();

  TH1D* h_inv = (TH1D*) TVX_eff->Clone("h_inv");
  for (int i = 1; i <= h_inv->GetNbinsX(); ++i) {
      double val = TVX_eff->GetBinContent(i);
      double err = TVX_eff->GetBinError(i);

      if (val > 0) {
          h_inv->SetBinContent(i, 1.0 / val);

          // error propagation: σ(1/x) = σ(x) / x^2
          h_inv->SetBinError(i, err / (val * val));
      } else {
          h_inv->SetBinContent(i, 0);
          h_inv->SetBinError(i, 0);
      }
  }

  TCanvas *cInv = new TCanvas("cInv", "Inverse TVX Efficiency", 800, 600);

  h_inv->GetYaxis()->SetTitle("1 / #epsilon_{TVX}");
  h_inv->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV/#it{c})");

  h_inv->SetMarkerStyle(20);
  h_inv->SetMarkerSize(1.2);
  h_inv->SetMarkerColor(kRed+1);
  h_inv->SetLineColor(kRed+1);

  // Important: range is no longer [0,1]
  h_inv->SetMinimum(0.0);
  h_inv->SetMaximum(5.0); // adjust depending on your values

  h_inv->Draw("E1 P");

  cInv->Update();


  return TVX_eff;
}

struct JetEfficiencies {
  TH1D* eff_zvtxToINEL;  // e_jet^{zvtx->INEL} = INEL / zvtx
  TH1D* eff_CollToTVX;   // e_jet^{Coll->TVX}  = kTVX / splitColl  (used in norm.)
  TH1D* eff_INELToColl;  // e_jet^{INEL->Coll} = splitColl / INEL  (-> +syst)
  TH1D* eff_INELToTVX;  // e_jet^{INEL->TVX} = kTVX / INEL  (-> +syst)
  TH1D* inv_zvtxToINEL;  // 1 / e_jet^{zvtx->INEL}
  TH1D* inv_CollToTVX;   // 1 / e_jet^{Coll->TVX}
  TH1D* inv_INELToColl;  // 1 / e_jet^{INEL->Coll}
  TH1D* inv_INELToTVX;  //  1 / e^{INEL->TVX} 
};

JetEfficiencies Get_MC_Jet_Eff(int iDataset, int iRadius) {

  JetEfficiencies result = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

  // ─── Get 2D histogram ───────────────────────────────────────────────────────
  TH2D* h2 = (TH2D*)file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get("jet-cross-section-efficiency/h2_jet_pt_part_eventselection");
  if (!h2) {
    std::cerr << "Error: 2D histogram not found!" << std::endl;
    return result;
  }

  // ─── Find Y bins for zvtx, kTVX, splitColl, INEL ────────────────────────────
  int bin_kTVX      = -1;
  int bin_splitColl = -1;
  int bin_INEL      = -1;
  int bin_zvtx      = -1;

  for (int i = 1; i <= h2->GetYaxis()->GetNbins(); i++) {
    TString label = h2->GetYaxis()->GetBinLabel(i);
    if (label.Contains("kTVX"))      bin_kTVX      = i;
    if (label.Contains("splitColl")) bin_splitColl = i;
    if (label.Contains("INEL"))      bin_INEL      = i;
    if (label.Contains("zvtx"))      bin_zvtx      = i;
  }

  if (bin_kTVX < 0 || bin_splitColl < 0 || bin_INEL < 0 || bin_zvtx < 0) {
    std::cerr << "Error: Could not find one or more labels in Y axis!" << std::endl;
    std::cerr << "Available labels: ";
    for (int i = 1; i <= h2->GetYaxis()->GetNbins(); i++)
      std::cerr << h2->GetYaxis()->GetBinLabel(i) << "  ";
    std::cerr << std::endl;
    return result;
  }

  // ─── Project to pT for all four selections ───────────────────────────────────
  TH1D* h_kTVX      = h2->ProjectionX("h_kTVX",      bin_kTVX,      bin_kTVX);
  TH1D* h_splitColl = h2->ProjectionX("h_splitColl",  bin_splitColl, bin_splitColl);
  TH1D* h_INEL      = h2->ProjectionX("h_INEL",       bin_INEL,      bin_INEL);
  TH1D* h_zvtx      = h2->ProjectionX("h_zvtx",       bin_zvtx,      bin_zvtx);

  // ─── Rebin ───────────────────────────────────────────────────────────────────
  TH1D* h_kTVX_reb      = (TH1D*)h_kTVX->Rebin(nBinPtJetsGen[iRadius],
                             Form("TVX_rebinned_dataset%d_R%d",       iDataset, iRadius), ptBinsJetsGen[iRadius]);
  TH1D* h_splitColl_reb = (TH1D*)h_splitColl->Rebin(nBinPtJetsGen[iRadius],
                             Form("splitColl_rebinned_dataset%d_R%d", iDataset, iRadius), ptBinsJetsGen[iRadius]);
  TH1D* h_INEL_reb      = (TH1D*)h_INEL->Rebin(nBinPtJetsGen[iRadius],
                             Form("INEL_rebinned_dataset%d_R%d",      iDataset, iRadius), ptBinsJetsGen[iRadius]);
  TH1D* h_zvtx_reb      = (TH1D*)h_zvtx->Rebin(nBinPtJetsGen[iRadius],
                             Form("zvtx_rebinned_dataset%d_R%d",      iDataset, iRadius), ptBinsJetsGen[iRadius]);

  h_kTVX_reb->Scale(1., "width");
  h_splitColl_reb->Scale(1., "width");
  h_INEL_reb->Scale(1., "width");
  h_zvtx_reb->Scale(1., "width");

  // ─── Plot all four spectra ────────────────────────────────────────────────────
  TCanvas* cSpectra = new TCanvas(Form("cSpectra_%d_%d", iDataset, iRadius),
                                  "Jet spectra: zvtx / INEL / splitColl / kTVX", 800, 600);
  h_zvtx_reb->SetLineColor(kViolet+1);   h_zvtx_reb->SetLineWidth(2);
  h_INEL_reb->SetLineColor(kBlue);       h_INEL_reb->SetLineWidth(2);
  h_splitColl_reb->SetLineColor(kGreen+2); h_splitColl_reb->SetLineWidth(2);
  h_kTVX_reb->SetLineColor(kRed);        h_kTVX_reb->SetLineWidth(2);
  h_zvtx_reb->SetTitle("Jet p_{T} spectra; p_{T} (GeV/c); Counts / GeV/c");
  h_zvtx_reb->Draw("hist");
  h_INEL_reb->Draw("hist same");
  h_splitColl_reb->Draw("hist same");
  h_kTVX_reb->Draw("hist same");

  TLegend* leg = new TLegend(0.55, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42);
  leg->AddEntry(h_zvtx_reb,      "zvtx",      "l");
  leg->AddEntry(h_INEL_reb,      "INEL",       "l");
  leg->AddEntry(h_splitColl_reb, "splitColl",  "l");
  leg->AddEntry(h_kTVX_reb,      "kTVX",       "l");
  leg->Draw();
  cSpectra->SetLogy();
  cSpectra->Update();

  // ─── Compute efficiencies ─────────────────────────────────────────────────────
  // e_jet^{zvtx->INEL} = INEL / zvtx
  TH1D* eff_zvtxToINEL = (TH1D*)h_INEL_reb->Clone(
                            Form("eff_zvtxToINEL_dataset%d_R%d", iDataset, iRadius));
  eff_zvtxToINEL->Divide(h_INEL_reb, h_zvtx_reb, 1.0, 1.0, "B");

  // e_jet^{Coll->TVX} = kTVX / splitColl
  TH1D* eff_CollToTVX = (TH1D*)h_kTVX_reb->Clone(
                           Form("eff_CollToTVX_dataset%d_R%d", iDataset, iRadius));
  eff_CollToTVX->Divide(h_kTVX_reb, h_splitColl_reb, 1.0, 1.0, "B");

  // e_jet^{INEL->Coll} = splitColl / INEL  (-> +syst)
  TH1D* eff_INELToColl = (TH1D*)h_splitColl_reb->Clone(
                            Form("eff_INELToColl_dataset%d_R%d", iDataset, iRadius));
  eff_INELToColl->Divide(h_splitColl_reb, h_INEL_reb, 1.0, 1.0, "B");

  // e_jet^{INEL->TVX} = kTVX / INEL 
  TH1D* eff_INELToTVX = (TH1D*)h_kTVX_reb->Clone(
                            Form("eff_INELToTVX_dataset%d_R%d", iDataset, iRadius));
  eff_INELToTVX->Divide(h_kTVX_reb, h_INEL_reb, 1.0, 1.0, "B");

  // ─── Plot all efficiencies ──────────────────────────────────────────────
  TCanvas* cEff = new TCanvas(Form("cEff_%d_%d", iDataset, iRadius),
                               "Jet efficiencies", 800, 600);
  eff_zvtxToINEL->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV/#it{c})");
  eff_zvtxToINEL->GetYaxis()->SetTitle("Efficiency");
  eff_zvtxToINEL->SetMarkerStyle(22); eff_zvtxToINEL->SetMarkerSize(1.2);
  eff_zvtxToINEL->SetMarkerColor(kViolet+1); eff_zvtxToINEL->SetLineColor(kViolet+1);
  eff_zvtxToINEL->SetMinimum(0.0); eff_zvtxToINEL->SetMaximum(1.3);
  eff_zvtxToINEL->SetTitle("Jet efficiencies");
  eff_zvtxToINEL->Draw("E1 P");

  eff_CollToTVX->SetMarkerStyle(20); eff_CollToTVX->SetMarkerSize(1.2);
  eff_CollToTVX->SetMarkerColor(kRed+1); eff_CollToTVX->SetLineColor(kRed+1);
  eff_CollToTVX->Draw("E1 P same");

  eff_INELToColl->SetMarkerStyle(21); eff_INELToColl->SetMarkerSize(1.2);
  eff_INELToColl->SetMarkerColor(kGreen+2); eff_INELToColl->SetLineColor(kGreen+2);
  eff_INELToColl->Draw("E1 P same");

  eff_INELToTVX->SetMarkerStyle(22); eff_INELToTVX->SetMarkerSize(1.2);
  eff_INELToTVX->SetMarkerColor(kBlue+2); eff_INELToTVX->SetLineColor(kBlue+2);
  eff_INELToTVX->Draw("E1 P same");

  TLegend* legEff = new TLegend(0.35, 0.15, 0.88, 0.38);
  legEff->SetBorderSize(0); legEff->SetFillStyle(0); legEff->SetTextFont(42);
  legEff->AddEntry(eff_zvtxToINEL, "#epsilon_{jet}^{zvtx#rightarrowINEL}",         "lp");
  legEff->AddEntry(eff_CollToTVX,  "#epsilon_{jet}^{Coll#rightarrowTVX} ",  "lp");
  legEff->AddEntry(eff_INELToColl, "#epsilon_{jet}^{INEL#rightarrowColl} ", "lp");
  legEff->AddEntry(eff_INELToTVX,  "#epsilon_{jet}^{INEL#rightarrowTVX} ",  "lp");
  legEff->Draw();
  cEff->Update();

  // ─── Compute inverses ─────────────────────────────────────────────────────────
  auto MakeInverse = [](TH1D* h, const char* name) -> TH1D* {
    TH1D* h_inv = (TH1D*)h->Clone(name);
    for (int i = 1; i <= h_inv->GetNbinsX(); ++i) {
      double val = h->GetBinContent(i);
      double err = h->GetBinError(i);
      if (val > 0) {
        h_inv->SetBinContent(i, 1.0 / val);
        h_inv->SetBinError(i, err / (val * val));
      } else {
        h_inv->SetBinContent(i, 0);
        h_inv->SetBinError(i, 0);
      }
    }
    return h_inv;
  };

  TH1D* inv_zvtxToINEL = MakeInverse(eff_zvtxToINEL, Form("inv_zvtxToINEL_dataset%d_R%d",  iDataset, iRadius));
  TH1D* inv_CollToTVX  = MakeInverse(eff_CollToTVX,  Form("inv_CollToTVX_dataset%d_R%d",   iDataset, iRadius));
  TH1D* inv_INELToColl = MakeInverse(eff_INELToColl, Form("inv_INELToColl_dataset%d_R%d",  iDataset, iRadius));
  TH1D* inv_INELToTVX  = MakeInverse(eff_INELToTVX,  Form("inv_INELToTVX_dataset%d_R%d",   iDataset, iRadius));

  // ─── Plot all three inverses ──────────────────────────────────────────────────
  TCanvas* cInv = new TCanvas(Form("cInv_%d_%d", iDataset, iRadius),
                               "Inverse efficiencies", 800, 600);
  inv_zvtxToINEL->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV/#it{c})");
  inv_zvtxToINEL->GetYaxis()->SetTitle("1 / #epsilon");
  inv_zvtxToINEL->SetMarkerStyle(22); inv_zvtxToINEL->SetMarkerSize(1.2);
  inv_zvtxToINEL->SetMarkerColor(kViolet+1); inv_zvtxToINEL->SetLineColor(kViolet+1);
  inv_zvtxToINEL->SetMinimum(0.0); inv_zvtxToINEL->SetMaximum(5.0);
  inv_zvtxToINEL->SetTitle("Inverse jet efficiencies");
  inv_zvtxToINEL->Draw("E1 P");

  inv_CollToTVX->SetMarkerStyle(20); inv_CollToTVX->SetMarkerSize(1.2);
  inv_CollToTVX->SetMarkerColor(kRed+1); inv_CollToTVX->SetLineColor(kRed+1);
  inv_CollToTVX->Draw("E1 P same");

  inv_INELToColl->SetMarkerStyle(21); inv_INELToColl->SetMarkerSize(1.2);
  inv_INELToColl->SetMarkerColor(kGreen+2); inv_INELToColl->SetLineColor(kGreen+2);
  inv_INELToColl->Draw("E1 P same");

  inv_INELToTVX->SetMarkerStyle(22); inv_INELToTVX->SetMarkerSize(1.2);
  inv_INELToTVX->SetMarkerColor(kBlue+2); inv_INELToTVX->SetLineColor(kBlue+2);
  inv_INELToTVX->Draw("E1 P same");

  TLegend* legInv = new TLegend(0.35, 0.65, 0.88, 0.85);
  legInv->SetBorderSize(0); legInv->SetFillStyle(0); legInv->SetTextFont(42);
  legInv->AddEntry(inv_zvtxToINEL, "1/#epsilon_{jet}^{zvtx#rightarrowINEL}", "lp");
  legInv->AddEntry(inv_CollToTVX,  "1/#epsilon_{jet}^{Coll#rightarrowTVX}",  "lp");
  legInv->AddEntry(inv_INELToColl, "1/#epsilon_{jet}^{INEL#rightarrowColl}", "lp");
  legInv->AddEntry(inv_INELToTVX,  "1/#epsilon_{jet}^{INEL#rightarrowTVX}",  "lp");
  legInv->Draw();
  cInv->Update();

  result = {eff_zvtxToINEL, eff_CollToTVX, eff_INELToColl, eff_INELToTVX, inv_zvtxToINEL, inv_CollToTVX, inv_INELToColl, inv_INELToTVX};
  return result;
}

double Get_SBP_Eff(int iDataset, int iRadius){

  TH1D* h_evt = (TH1D*)file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(
    mcIsWeighted ? "jet-cross-section-efficiency/h_mccollisions_eventselection_weighted" : "jet-cross-section-efficiency/h_mccollisions_eventselection");
  int bin_evt_kITSROFBorder = h_evt->GetXaxis()->FindBin("kITSROFBorder");
  int bin_evt_NoSameBunchPileup = h_evt->GetXaxis()->FindBin("NoSameBunchPileup");
  double N_kITSROFBorder = h_evt->GetBinContent(bin_evt_kITSROFBorder);
  double N_NoSameBunchPileup = h_evt->GetBinContent(bin_evt_NoSameBunchPileup);
  double efficiency = N_NoSameBunchPileup / N_kITSROFBorder;

  std::cout << "=========================================="<< std::endl;
  std::cout << "SBP efficiency = " << efficiency << std::endl;
  std::cout << "=========================================="<< std::endl;

  return efficiency; 
}

void Draw_Sigma_spectrum(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  bool splitTestControlMC = true;

  TH1D* H1D_jetPt_unfolded;
  TString partialUniqueSpecifier;
  int unfoldParameter;
  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TH1D* measuredInput;

  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
  unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).first; // 1/N_ev d^2N/dpTdeta where N_ev is bin8

  const double sigma_vdm_mb = 46.46;   // mb
  const double vdmRun3      = 0.045;   // 4.5% relative lumi unc (Run 3)
  const double sigma_O2 =  50.3; //mb line 1685 https://github.com/AliceO2Group/O2Physics/blob/f4ec4e509734228d21fc11a707541c629186c7b5/Common/Tools/EventSelectionModule.h

  TH1* hCounter = (TH1*)file_O2Analysis_list[iDataset]->Get("jet-luminosity-calculator/counter");
  double nEventsSel8 = hCounter->GetBinContent(8);
  double nEventsBin4 = hCounter->GetBinContent(4);
  double nBCTVX_Bin2 = hCounter->GetBinContent(2);
  double nColTVXBin6 = hCounter->GetBinContent(6);
  // cout << "nEventsSel8 = " << nEventsSel8 << ", nEventsBin4 = " << nEventsBin4 << endl;
  // cout << "nColTVXBin6 = " << nColTVXBin6 << ", nBCTVX_Bin2 = " << nBCTVX_Bin2 << endl;
  // cout << "Eff_evt = " << nEventsSel8/ nColTVXBin6 << endl;
  // cout << "Lumi = " << nBCTVX_Bin2/ sigma_vdm_mb << endl;

  auto effs = Get_MC_Jet_Eff(iDataset, iRadius);
  TH1D* eff_INELToTVX  = effs.eff_INELToTVX;   
  double z_vtx_eff = Get_zVertex_reconstruction_efficiency(iDataset, iRadius);
  double SBP_eff = Get_SBP_Eff(iDataset, iRadius);
  double lumi = nEventsBin4 / sigma_vdm_mb ;
  cout << "Lumi = " << lumi << endl;
  const double scalingFactor = nEventsSel8 /(lumi * z_vtx_eff * SBP_eff ); 
  TH1D* H1D_Xsection_woINELToTVXeff = (TH1D*)H1D_jetPt_unfolded->Clone("Xsection_woINELToTVXeff");  
  
  H1D_Xsection_woINELToTVXeff->Scale(scalingFactor);
  TH1D* H1D_Xsection_wINELToTVXeff = (TH1D*)H1D_Xsection_woINELToTVXeff->Clone("H1D_Xsection_wINELToTVXeff");  
  H1D_Xsection_wINELToTVXeff->Divide(eff_INELToTVX);

  // TString* pdfName = new TString("XSextion_jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  // TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
  // TString* texJetXsection_d2Sigmadptdeta = new TString("d^{2}#sigma_{jet}/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
  // Draw_TH1_Histogram(H1D_jetPt_unfolded, textContext, pdfName, texPtX, texJetXsection_d2Sigmadptdeta, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  // double newbins[] = {10., 12, 14, 16, 18, 20., 25, 30., 40., 50., 60., 70., 85., 100., 140.}; 
  // double newbins[] = {10., 20., 30., 40., 50., 60., 70., 80., 100., 140.}; 
  // int n_newbins = sizeof(newbins)/sizeof(newbins[0]) - 1; // = 10
  // TH1D* H1D_Xsection_woINELToTVXeff_reb = ReweightedRebin(H1D_Xsection_woINELToTVXeff, "H1D_jetPt_rebinned", n_newbins, newbins);
  // TH1D* H1D_Xsection_wINELToTVXeff_reb = ReweightedRebin(H1D_Xsection_wINELToTVXeff, "H1D_jetPt_rebinned", n_newbins, newbins);

  // TString* pdfName_Rebinned = new TString("Rebinned_XSextion_jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  // Draw_TH1_Histogram(H1D_Xsection_woINELToTVXeff_reb, textContext, pdfName_Rebinned, texPtX, texJetXsection_d2Sigmadptdeta, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  // if (writeOutputRootFile) {
  //   cout << "######################### FILE CREATED IN PRINCIPLE##################" << endl; 
  //   TFile* outFile = new TFile("Xsection_2ndJune.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  //   TString histoName = Form("H1D_XSection_Pt_Unfolded_w_%s", mcIsWeighted ? "JJ" : "MB");
  //   TString histoName_Rebinned = Form("H1D_XSection_Pt_Unfolded_Rebinned_w_%s", mcIsWeighted ? "JJ" : "MB");
  //   H1D_jetPt_unfolded->Write(histoName);
  //   H1D_jetPt_rebinned->Write(histoName_Rebinned);
  //   outFile->Close();
  //   cout << "######################### HISTO SAVED IN PRINCIPLE ##################" << endl; 
  // }
 
  // ---- PLOTTING ----
  // // ---- Run 2 ---- paper : 
  // TFile* fRun2 = TFile::Open("/Users/tabikh/Documents/PhdWork/MyWork/Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02.root", "READ");
  // TH1D* hRun2      = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1");
  // TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1_e1");
  // TH1D* hRun2_sys  = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1_e2");
  // hRun2->SetDirectory(0); hRun2_sys->SetDirectory(0);

  // // transfer Run-2 stat errors (stored as contents in _e1) into hRun2 bin errors
  // for (int i = 1; i <= hRun2->GetNbinsX(); ++i)
  //   hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));
  // fRun2->Close();

  // ---- Run 2 --- paper : ins2637686
  TFile* fRun2 = TFile::Open("/Users/tabikh/Documents/PhdWork/MyWork/Datasets/Run2_ins2637686/Figure3b_bottom_value_stat_sys.root", "READ");
  TH1D* hRun2      = (TH1D*)fRun2->Get("h_value_R0p2");
  TH1D* hRun2_stat = (TH1D*)fRun2->Get("h_stat_R0p2");
  TH1D* hRun2_sys  = (TH1D*)fRun2->Get("h_sys_R0p2");
  hRun2->SetDirectory(0); hRun2_sys->SetDirectory(0);

  // transfer Run-2 stat errors (stored as contents in _e1) into hRun2 bin errors
  for (int i = 1; i <= hRun2->GetNbinsX(); ++i)
    hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));
  fRun2->Close();

  // ---- Powheg -----
  const char* powhegFile = "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_woPtLeadCut/POWHEG_Uncertainties.root";
  TFile* fPow = TFile::Open(powhegFile);
  TH1F* hPow    = (TH1F*)fPow->Get("Central_Inclusive_R02");
  TH1F* hPowUnc = (TH1F*)fPow->Get("hTotalUnc_Inclusive_R02");   // RELATIVE
  hPow->SetDirectory(0); hPowUnc->SetDirectory(0);

  std::vector<TH1*> hists = {H1D_Xsection_wINELToTVXeff, H1D_Xsection_woINELToTVXeff, hRun2, hRun2_sys, hPow, hPowUnc};
  if (!HarmonizeBinning(hists, 10., 100.)) return;

  H1D_Xsection_wINELToTVXeff  = (TH1D*)hists[0];
  H1D_Xsection_woINELToTVXeff = (TH1D*)hists[1];
  hRun2     = (TH1D*)hists[2];
  hRun2_sys = (TH1D*)hists[3];
  hPow      = (TH1F*)hists[4];
  hPowUnc   = (TH1F*)hists[5];

    // ---- entries ----
  SpecEntry run2;
  run2.h        = hRun2; run2.style = kData; run2.color = kGray+3; run2.marker = 24;
  run2.label    = "ALICE Run 2"; run2.sysKind = kSysAbsolute; run2.hSys = hRun2_sys;
  run2.sysColor = kGray+1; run2.sysLabel = "Sys. unc.";

  SpecEntry pow;
  pow.h = hPow;   pow.style = kTheory;   pow.color = kGreen+2;  pow.lwidth = 3;
  pow.label = "POWHEG+Pythia8 5.36 TeV"; pow.sysKind = kSysRelative; pow.hSys = hPowUnc;
  pow.sysColor = kGreen-6; pow.sysLabel = "POWHEG unc.";

  SpecEntry run3_wINELToTVXeff;
  run3_wINELToTVXeff.h       = H1D_Xsection_wINELToTVXeff;          // with INEL->TVX eff
  run3_wINELToTVXeff.style   = kData;
  run3_wINELToTVXeff.color   = kAzure+2;
  run3_wINELToTVXeff.marker  = 20;               // filled circle
  run3_wINELToTVXeff.label   = "ALICE Run 3 (w #varepsilon_{INEL#rightarrowTVX})";
  run3_wINELToTVXeff.sysKind = kSysNone;         // no systematic

  SpecEntry run3_woINELToTVXeff;
  run3_woINELToTVXeff.h       = H1D_Xsection_woINELToTVXeff;        // without INEL->TVX eff
  run3_woINELToTVXeff.style   = kData;
  run3_woINELToTVXeff.color   = kViolet+1;
  run3_woINELToTVXeff.marker  = 25;              // open square (distinct from run2's 24)
  run3_woINELToTVXeff.label   = "ALICE Run 3 (wo #varepsilon_{INEL#rightarrowTVX})";
  run3_woINELToTVXeff.sysKind = kSysNone;        // no systematic

  DrawSpectraWithRatio(
      { pow, run3_wINELToTVXeff, run3_woINELToTVXeff, run2 },  // reference first
      0,                            // refIndex=0 -> run3_wINELToTVXeff is denominator
      10., 100.,
      "d^{2}#sigma_{jet}/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]",
      "#it{p}_{T} (GeV/#it{c})", "X / POWHEG",
      { "Inclusive jets, anti-#it{k}_{T}, R = 0.2", "|#eta_{jet}| < 0.7" },
      "cFour", 0., 2.5, 0.30,
      -1.,                          // vdmRatioRel (full-width): off
      -1.,        // globalNormRel: off
      "",         // globalNormLabel
      0., 0.,     // yMinUser, yMaxUser: auto
      0.045, 97);     // vdmBoxRel: small 4.5% dark box at unity at pt = 97 GeV/c
}


void Run2vsRun3_comparison(int         iDataset,
                           int         iRadius,
                           int         unfoldParameterInput,
                           std::string options){
  // ---------------------------------------------------------------------------
  // 0.  GLOBAL STYLE
  // ---------------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "X");
  gStyle->SetTickLength(0.02, "Y");
  gStyle->SetEndErrorSize(0);   // no end caps on stat bars

  // ---------------------------------------------------------------------------
  // 1.  RUN-3 CROSS SECTION
  //     Stat errors live in hRun3 bin errors (set by the unfolding step)
  // ---------------------------------------------------------------------------
  TH1D* H1D_jetPt_unfolded = nullptr;
  TH1D* measuredInput      = nullptr;

  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput,iDataset, iRadius, unfoldParameterInput, options);
  int binStart = H1D_jetPt_unfolded->FindBin(10.0);
  int binEnd   = H1D_jetPt_unfolded->FindBin(100.0)-1;

  RewriteHistogramRange(H1D_jetPt_unfolded, binStart, binEnd);

  // vdM luminosity constants
  const double sigma_vdm_mb = 46.46;   // mb
  const double vdmRun3      = 0.045;   // 4.5% relative lumi unc (Run 3)
  const double vdmRun2      = 0.023;   // 2.3% relative lumi unc (Run 2) https://alice-notes.web.cern.ch/system/files/notes/analysis/534/2019-01-16-Analysis_note_pp5TeV_CrossSection_v1.1.pdf
  // Fully uncorrelated (independent vdM campaigns) -> add in quadrature
  const double vdmRatio = TMath::Sqrt(vdmRun3*vdmRun3 + vdmRun2*vdmRun2);

  double z_vtx_eff = Get_zVertex_reconstruction_efficiency(iDataset, iRadius);
  std::cout << "\n>>> z-vertex efficiency = " << z_vtx_eff << "\n\n";

  // TH1D* TVX_eff     = Get_TVX_Eff(iDataset, iRadius);
  TH1* hCounter = (TH1*)file_O2Analysis_list[iDataset]->Get("jet-luminosity-calculator/counter");
  double nEventsSel8 = hCounter->GetBinContent(8);
  double nEventsBin4 = hCounter->GetBinContent(4);
  double nBCTVX_Bin2 = hCounter->GetBinContent(2);
  double nColTVXBin6 = hCounter->GetBinContent(6);
  cout << "nEventsSel8 = " << nEventsSel8 << ", nEventsBin4 = " << nEventsBin4 << endl;
  cout << "nColTVXBin6 = " << nColTVXBin6 << ", nBCTVX_Bin2 = " << nBCTVX_Bin2 << endl;
  cout << "Eff_TVX = " << nEventsSel8/ nColTVXBin6 << endl;
  cout << "Lumi = " << nBCTVX_Bin2/ sigma_vdm_mb << endl;

  // const double scalingFactor = sigma_vdm_mb * nEventsSel8 /(z_vtx_eff * nEventsBin4);
  // const double scalingFactor = sigma_vdm_mb /(z_vtx_eff );
  // const double scalingFactor = sigma_vdm_mb *nColTVXBin6 /(z_vtx_eff * nBCTVX_Bin2);
  const double scalingFactor = sigma_vdm_mb *nColTVXBin6 /(nBCTVX_Bin2); // Nima suggestion

  cout << "scalingFactor = " << scalingFactor << endl;


  TH1D* hRun3 = (TH1D*)H1D_jetPt_unfolded->Clone("hRun3");
  hRun3->Scale(scalingFactor);
  // hRun3->Divide(TVX_eff);
  // hRun3 bin errors = stat uncertainties, scaled + divided consistently

  // Run-3 systematic: hRun3_sys stores RELATIVE uncertainties
  bool  UE_trk_method = true;
  TH1D* hRun3_sys = Get_Total_Systematic(UE_trk_method);  // RELATIVE
  RewriteHistogramRange(hRun3_sys, binStart, binEnd);
  for (int i = 1; i <= hRun3_sys->GetNbinsX(); i++) {
    printf("pT=[%.0f-%.0f]: Sys. unc. = %.4f\n",
          hRun3_sys->GetBinLowEdge(i),
          hRun3_sys->GetBinLowEdge(i) + hRun3_sys->GetBinWidth(i),
          hRun3_sys->GetBinContent(i));
  }


  if (writeOutputRootFile) {
    cout << "######################### FILE CREATED PRINCIPLE##################" << endl; 
    TFile* outFile = new TFile("XSection_pp536_R02_ptLeadcut5.root", "RECREATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
    TString histoName = Form("H1D_XSection_pp536_R02_ptLeadcut5");
    TString histoName_syst = Form("H1D_RelativeSyst_on_XSection");
    hRun3->Write(histoName);
    hRun3_sys->Write(histoName_syst);
    outFile->Close();
    cout << "######################### HISTO SAVED IN PRINCIPLE ##################" << endl; 
  }
  // ------------
  // POWHEG : /Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_ptlead5/POWHEG_Uncertainties.root
  // ------------

  // ---------------------------------------------------------------------------
  // 2.  RUN-2 CROSS SECTION
  //     hRun2_stat has stat errors as bin CONTENTS -> copy to bin errors
  //     hRun2_sys  has systematic errors as ABSOLUTE bin contents
  // ---------------------------------------------------------------------------
  TFile* fRun2 = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02_LeadPtCut5.root", "READ");
  TH1D* hRun2      = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1");
  TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e1");
  TH1D* hRun2_sys  = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e2");

  // *** Transfer Run-2 stat errors: contents of hRun2_stat -> errors of hRun2 ***
  // This must happen before any ratio or TGraphErrors construction
  for (int i = 1; i <= hRun2->GetNbinsX(); i++) {
      hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));
  }
  RewriteHistogramRange(hRun2, binStart, binEnd);
  RewriteHistogramRange(hRun2_sys, binStart, binEnd);

  const int nbins = hRun3->GetNbinsX();

  // ---------------------------------------------------------------------------
  // 3.  SYS-BOX TGRAPHS FOR SPECTRA  (vdM handled separately)
  // ---------------------------------------------------------------------------
  TGraphErrors* gRun3_sys = new TGraphErrors(nbins);
  TGraphErrors* gRun2_sys = new TGraphErrors(nbins);
  // vdM bands: TGraphErrors so the band tracks the curve in log scale
  TGraphErrors* gRun3_vdm = new TGraphErrors(nbins);
  TGraphErrors* gRun2_vdm = new TGraphErrors(nbins);

  for (int i = 1; i <= nbins; i++) {

      // Run 3
      const double x3  = hRun3->GetBinCenter(i);
      const double dx3 = hRun3->GetBinWidth(i) / 2.0;
      const double y3  = hRun3->GetBinContent(i);
      // hRun3_sys is RELATIVE -> convert to absolute
      const double abs_sys3 = hRun3_sys->GetBinContent(i) * y3;

      gRun3_sys->SetPoint(i-1, x3, y3);
      gRun3_sys->SetPointError(i-1, dx3, abs_sys3);

      gRun3_vdm->SetPoint(i-1, x3, y3);
      gRun3_vdm->SetPointError(i-1, dx3, y3 * vdmRun3);

      // Run 2
      const double x2  = hRun2->GetBinCenter(i);
      const double dx2 = hRun2->GetBinWidth(i) / 2.0;
      const double y2  = hRun2->GetBinContent(i);
      // hRun2_sys is already ABSOLUTE
      const double abs_sys2 = hRun2_sys->GetBinContent(i);

      gRun2_sys->SetPoint(i-1, x2, y2);
      gRun2_sys->SetPointError(i-1, dx2, abs_sys2);

      gRun2_vdm->SetPoint(i-1, x2, y2);
      gRun2_vdm->SetPointError(i-1, dx2, y2 * vdmRun2);
  }

  // ---------------------------------------------------------------------------
  // 4.  RATIO + RATIO SYS GRAPH
  // ---------------------------------------------------------------------------
  TH1D* hRatio = (TH1D*)hRun3->Clone("hRatio");
  hRatio->Reset("ICES");

  TGraphErrors* gRatio_sys = new TGraphErrors(nbins);

  for (int i = 1; i <= nbins; i++) {

      const double x  = hRun3->GetBinCenter(i);
      const double dx = hRun3->GetBinWidth(i) / 2.0;

      const double y3 = hRun3->GetBinContent(i);
      const double y2 = hRun2->GetBinContent(i);
      if (y2 <= 0 || y3 <= 0) continue;

      const double r = y3 / y2;
      hRatio->SetBinContent(i, r);

      // Statistical: both histograms now have proper bin errors
      const double e3    = hRun3->GetBinError(i);   // from unfolding
      const double e2    = hRun2->GetBinError(i);   // transferred in step 2
      const double e_rat = r * TMath::Sqrt(TMath::Power(e3/y3, 2) +
                                           TMath::Power(e2/y2, 2));
      hRatio->SetBinError(i, e_rat);

      // Systematic: convert both to relative, then quadrature sum
      const double rel_sys3 = hRun3_sys->GetBinContent(i);          // RELATIVE
      const double rel_sys2 = hRun2_sys->GetBinContent(i) / y2;     // ABSOLUTE / value

      const double rel_sys_ratio = TMath::Sqrt(rel_sys3*rel_sys3 +
                                               rel_sys2*rel_sys2);
      gRatio_sys->SetPoint(i-1, x, r);
      gRatio_sys->SetPointError(i-1, dx, r * rel_sys_ratio);
  }

  // ---------------------------------------------------------------------------
  // 5.  CANVAS   (65% top / 35% bottom)
  // ---------------------------------------------------------------------------
  TCanvas* c = new TCanvas("c_Run2vsRun3",
                            "Run2 vs Run3 Inclusive Jet Cross Section",
                            800, 900);
  c->SetFillColor(0);

  const double split = 0.35;

  TPad* pTop = new TPad("pTop", "pTop", 0, split, 1, 1);
  TPad* pBot = new TPad("pBot", "pBot", 0, 0,     1, split);

  pTop->SetLeftMargin(0.14);  pTop->SetRightMargin(0.04);
  pTop->SetTopMargin(0.05);   pTop->SetBottomMargin(0.01);
  pTop->SetLogy();
  pTop->SetTickx(1);  pTop->SetTicky(1);

  pBot->SetLeftMargin(0.14);  pBot->SetRightMargin(0.04);
  pBot->SetTopMargin(0.01);   pBot->SetBottomMargin(0.30);
  pBot->SetTickx(1);  pBot->SetTicky(1);

  pTop->Draw();
  pBot->Draw();

  // ---------------------------------------------------------------------------
  // 6.  COLORS
  // ---------------------------------------------------------------------------
  const Color_t cRun3     = kAzure+2;
  const Color_t cRun3_sys = kAzure-9;
  const Color_t cRun3_vdm = kAzure-7;

  const Color_t cRun2     = kRed+1;
  const Color_t cRun2_sys = kRed-9;
  const Color_t cRun2_vdm = kRed-10;

  const Color_t cVdmRatio = kOrange+1;

  // ---------------------------------------------------------------------------
  // 7.  TOP PAD  —  spectra
  // ---------------------------------------------------------------------------
  pTop->cd();

  // Shape sys boxes
  gRun3_sys->SetFillColorAlpha(cRun3_sys, 0.50);
  gRun3_sys->SetLineColor(cRun3);
  gRun3_sys->SetMarkerSize(0);

  gRun2_sys->SetFillColorAlpha(cRun2_sys, 0.50);
  gRun2_sys->SetLineColor(cRun2);
  gRun2_sys->SetMarkerSize(0);

  // vdM bands (TGraphErrors tracking the curve in log scale)
  gRun3_vdm->SetFillColorAlpha(cRun3_vdm, 0.22);
  gRun3_vdm->SetLineColorAlpha(cRun3_vdm, 0.55);
  gRun3_vdm->SetMarkerSize(0);

  gRun2_vdm->SetFillColorAlpha(cRun2_vdm, 0.22);
  gRun2_vdm->SetLineColorAlpha(cRun2_vdm, 0.55);
  gRun2_vdm->SetMarkerSize(0);

  // Data point styles
  hRun3->SetMarkerStyle(20);   // filled circle
  hRun3->SetMarkerSize(1.1);
  hRun3->SetMarkerColor(cRun3);
  hRun3->SetLineColor(cRun3);
  hRun3->SetLineWidth(1);

  hRun2->SetMarkerStyle(24);   // open circle
  hRun2->SetMarkerSize(1.1);
  hRun2->SetMarkerColor(cRun2);
  hRun2->SetLineColor(cRun2);
  hRun2->SetLineWidth(1);

  // Axis formatting (suppress x in top pad)
  const double labelSz = 0.055;
  const double titleSz = 0.060;

  hRun3->GetYaxis()->SetTitle(
      "d^{2}#sigma_{jet} / d#it{p}_{T} d#eta  (mb GeV^{-1} #it{c})");
  hRun3->GetYaxis()->SetLabelSize(labelSz);
  hRun3->GetYaxis()->SetTitleSize(titleSz);
  hRun3->GetYaxis()->SetTitleOffset(1.15);
  hRun3->GetXaxis()->SetLabelSize(0);
  hRun3->GetXaxis()->SetTitleSize(0);

  // Draw order: vdM (outermost) -> shape sys -> stat bars -> points on top
  hRun3->Draw("E1 P");
  pTop->Update();
  pTop->GetFrame()->SetX1(10);
  pTop->GetFrame()->SetX2(100);
  // hRun3->GetXaxis()->SetRangeUser(10, 100);
  hRun3->GetXaxis()->SetLimits(10, 100);
  gRun3_vdm->Draw("E2 SAME");
  gRun2_vdm->Draw("E2 SAME");
  gRun3_sys->Draw("E2 SAME");
  gRun2_sys->Draw("E2 SAME");
  hRun3->Draw("E1 P SAME");
  hRun2->Draw("E1 P SAME");

  // Legend
  TLegend* leg = new TLegend(0.52, 0.53, 0.95, 0.93);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetTextFont(42);
  leg->AddEntry(hRun3,     "ALICE Run 3, pp 5.36 TeV",                   "P");
  leg->AddEntry(gRun3_sys, "Sys. unc.",                          "F");
  // leg->AddEntry(gRun3_vdm, Form("Lumi. unc. vdM (%.1f%%)",vdmRun3*100.),"F");
  leg->AddEntry((TObject*)nullptr, "", "");   // spacer
  leg->AddEntry(hRun2,     "ALICE Run 2, pp 5.02 TeV",                   "P");
  leg->AddEntry(gRun2_sys, "Sys. unc.",                          "F");
  // leg->AddEntry(gRun2_vdm, Form("Lumi. unc. vdM (%.1f%%)",vdmRun2*100.),"F");
  leg->Draw();

  // ALICE label
  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.055);
  // tex.DrawLatex(0.17, 0.88, "ALICE");
  tex.SetTextSize(0.046);
  tex.DrawLatex(0.2, 0.88,
      Form("Inclusive jets, anti-#it{k}_{T}, #it{R} = %.1f",
           arrayRadius[iRadius]));
  // tex.DrawLatex(0.17, 0.74, "pp,  #sqrt{#it{s}} = 5.02 TeV");
  tex.DrawLatex(0.2, 0.81, "|#eta_{jet}| < 0.7");

  // ---------------------------------------------------------------------------
  // 8.  BOTTOM PAD  —  ratio
  // ---------------------------------------------------------------------------
  pBot->cd();
 
  const double sf = (1.0 - split) / split;   // font scale for smaller pad
 
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(1.1);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->SetLineWidth(1);
 
  hRatio->GetYaxis()->SetTitle("Run 3 / Run 2");
  hRatio->GetYaxis()->SetRangeUser(0.0, 2.4);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetLabelSize(labelSz * sf);
  hRatio->GetYaxis()->SetTitleSize(titleSz * sf);
  hRatio->GetYaxis()->SetTitleOffset(1.15 / sf);
 
  hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hRatio->GetXaxis()->SetLabelSize(labelSz * sf);
  hRatio->GetXaxis()->SetTitleSize(titleSz * sf);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetXaxis()->SetTickLength(0.05);
 
  // Shape sys boxes
  gRatio_sys->SetFillColorAlpha(kGray, 0.55);
  gRatio_sys->SetLineColor(kGray+2);
  gRatio_sys->SetMarkerSize(0);
 
  // vdM combined band: full-width TBox at unity (linear scale -> TBox correct)
  const double xMin_bot = hRatio->GetXaxis()->GetXmin();
  const double xMax_bot = hRatio->GetXaxis()->GetXmax();
 
  TBox* vdmBand = new TBox(xMin_bot, 1.0 - vdmRatio,
                            xMax_bot, 1.0 + vdmRatio);
  vdmBand->SetFillColorAlpha(cVdmRatio, 0.45);
  vdmBand->SetLineColorAlpha(cVdmRatio, 0.70);
  vdmBand->SetLineWidth(1);
 
  // Draw order: axes -> vdM band -> shape sys -> unity line -> points
  hRatio->Draw("E1 P");
  pBot->Update();
  pBot->GetFrame()->SetX1(10);
  pBot->GetFrame()->SetX2(100);
  hRatio->GetXaxis()->SetRangeUser(10, 100);
  hRatio->GetXaxis()->SetLimits(10, 100);
  vdmBand->Draw("SAME");
  gRatio_sys->Draw("E2 SAME");
 
  TLine unity(xMin_bot, 1.0, xMax_bot, 1.0);
  unity.SetLineColor(kBlack);
  unity.SetLineStyle(2);
  unity.SetLineWidth(1);
  unity.Draw("SAME");
 
  hRatio->Draw("E1 P SAME");
 
  // Legend (ratio pad)
  TH1F hDumSys("hDumSys", "", 1, 0, 1);
  hDumSys.SetFillColorAlpha(kGray, 0.55);
  hDumSys.SetLineColor(kGray+2);
 
  // Use a TBox for the vdM legend entry - matches exactly what is drawn
  // and is always visible regardless of alpha
  TBox* legVdmBox = new TBox(0, 0, 1, 1);
  legVdmBox->SetFillColorAlpha(cVdmRatio, 0.45);  // slightly higher alpha for legibility
  legVdmBox->SetLineColor(cVdmRatio);
  legVdmBox->SetLineWidth(1);
 
  TLegend legR(0.52, 0.76, 0.95, 0.99);
  legR.SetBorderSize(0);
  legR.SetFillStyle(0);
  legR.SetTextSize(0.087);
  legR.SetTextFont(42);
  legR.AddEntry(hRatio,    "Run 3 / Run 2",                              "P");
  legR.AddEntry(&hDumSys,  "Sys. unc. (shape, quadr.)",                  "F");
  legR.AddEntry(legVdmBox,
      Form("Lumi. unc. vdM (%.1f%%, uncorr.)", vdmRatio*100.),          "F");
  legR.Draw();
 
  // Label showing vdM percentage at the right edge of the band
  TLatex texR;
  texR.SetNDC(kFALSE);
  texR.SetTextFont(42);
  texR.SetTextSize(labelSz * sf * 0.85);
  texR.SetTextColor(cVdmRatio);
  texR.DrawLatex(xMax_bot * 0.65,
                 1.0 + vdmRatio + 0.06,
                 Form("#pm%.1f%%", vdmRatio * 100.));
 
  c->Update();

  // // ---------------------------------------------------------------------------
  // // 9.  SAVE
  // // ---------------------------------------------------------------------------
  // TString tag = Form("R%.0f", arrayRadius[iRadius] * 10);
  // c->SaveAs(Form("Run2vsRun3_jetXsec_%s.pdf", tag.Data()));
  // c->SaveAs(Form("Run2vsRun3_jetXsec_%s.png", tag.Data()));

  // std::cout << "\n>>> Saved: Run2vsRun3_jetXsec_" << tag << ".[pdf/png]\n";
}


void Run2vsRun3_comparison_R04(int         iDataset,
                           int         iRadius,
                           int         unfoldParameterInput,
                           std::string options){
  // ---------------------------------------------------------------------------
  // 0.  GLOBAL STYLE
  // ---------------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "X");
  gStyle->SetTickLength(0.02, "Y");
  gStyle->SetEndErrorSize(0);   // no end caps on stat bars

  // ---------------------------------------------------------------------------
  // 1.  RUN-3 CROSS SECTION
  //     Stat errors live in hRun3 bin errors (set by the unfolding step)
  // ---------------------------------------------------------------------------
  TH1D* H1D_jetPt_unfolded = nullptr;
  TH1D* measuredInput      = nullptr;

  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput,iDataset, iRadius, unfoldParameterInput, options);
  int binStart = H1D_jetPt_unfolded->FindBin(10.0);
  int binEnd   = H1D_jetPt_unfolded->FindBin(100.0)-1;

  RewriteHistogramRange(H1D_jetPt_unfolded, binStart, binEnd);

  // vdM luminosity constants
  const double sigma_vdm_mb = 46.46;   // mb
  const double vdmRun3      = 0.045;   // 4.5% relative lumi unc (Run 3)
  const double vdmRun2      = 0.012;   // 1.2% relative lumi unc (Run 2)
  // Fully uncorrelated (independent vdM campaigns) -> add in quadrature
  const double vdmRatio = TMath::Sqrt(vdmRun3*vdmRun3 + vdmRun2*vdmRun2);

  double z_vtx_eff = Get_zVertex_reconstruction_efficiency(iDataset, iRadius);
  std::cout << "\n>>> z-vertex efficiency = " << z_vtx_eff << "\n\n";

  // TH1D* TVX_eff     = Get_TVX_Eff(iDataset, iRadius);
  TH1* hCounter = (TH1*)file_O2Analysis_list[iDataset]->Get("jet-luminosity-calculator/counter");
  double nEventsSel8 = hCounter->GetBinContent(8);
  double nEventsBin4 = hCounter->GetBinContent(4);
  double nBCTVX_Bin2 = hCounter->GetBinContent(2);
  double nColTVXBin6 = hCounter->GetBinContent(6);
  cout << "nEventsSel8 = " << nEventsSel8 << ", nEventsBin4 = " << nEventsBin4 << endl;
  cout << "nColTVXBin6 = " << nColTVXBin6 << ", nBCTVX_Bin2 = " << nBCTVX_Bin2 << endl;
  cout << "Eff_TVX = " << nEventsSel8/ nColTVXBin6 << endl;
  cout << "Lumi = " << nBCTVX_Bin2/ sigma_vdm_mb << endl;

  // const double scalingFactor = sigma_vdm_mb * nEventsSel8 /(z_vtx_eff * nEventsBin4);
  // const double scalingFactor = sigma_vdm_mb /(z_vtx_eff );
  // const double scalingFactor = sigma_vdm_mb *nColTVXBin6 /(z_vtx_eff * nBCTVX_Bin2);
  const double scalingFactor = sigma_vdm_mb *nColTVXBin6 /(nBCTVX_Bin2); // Nima suggestion

  cout << "scalingFactor = " << scalingFactor << endl;


  TH1D* hRun3 = (TH1D*)H1D_jetPt_unfolded->Clone("hRun3");
  hRun3->Scale(scalingFactor);
  // hRun3->Divide(TVX_eff);
  // hRun3 bin errors = stat uncertainties, scaled + divided consistently

  // ---------------------------------------------------------------------------
  // 2.  RUN-2 CROSS SECTION
  //     hRun2_stat has stat errors as bin CONTENTS -> copy to bin errors
  //     hRun2_sys  has systematic errors as ABSOLUTE bin contents
  // ---------------------------------------------------------------------------
  // TFile* fRun2 = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02_LeadPtCut5.root", "READ");
  // TH1D* hRun2      = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y3");
  // TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y3_e1");
  // TH1D* hRun2_sys  = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y3_e2");

  TFile* fRun2 = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02.root", "READ");
  TH1D* hRun2      = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1");
  TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1_e1");

  // *** Transfer Run-2 stat errors: contents of hRun2_stat -> errors of hRun2 ***
  // This must happen before any ratio or TGraphErrors construction
  for (int i = 1; i <= hRun2->GetNbinsX(); i++) {
      hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));
  }
  RewriteHistogramRange(hRun2, binStart, binEnd);

  TCanvas *c = new TCanvas("c", "Run2 vs Run3", 800, 900);

  // Define pads
  TPad *pad1 = new TPad("pad1", "upper pad", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "lower pad", 0, 0.0, 1, 0.3);

  // Pad styling
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);
  pad1->SetTicks();

  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);
  pad2->SetTicks();

  pad1->Draw();
  pad2->Draw();

  // =====================
  // Upper pad
  // =====================
  pad1->cd();

  hRun3->SetLineColor(kBlue+1);
  hRun3->SetMarkerColor(kBlue+1);
  hRun3->SetMarkerStyle(20);
  hRun3->SetMarkerSize(1.0);
  hRun3->SetLineWidth(2);

  hRun2->SetLineColor(kRed+1);
  hRun2->SetMarkerColor(kRed+1);
  hRun2->SetMarkerStyle(24);
  hRun2->SetMarkerSize(1.0);
  hRun2->SetLineWidth(2);

  // Set axis titles
  hRun3->SetTitle("");
  hRun3->GetYaxis()->SetTitle("d^{2}#sigma_{jet} / d#it{p}_{T} d#eta  (mb GeV^{-1} #it{c})");
  hRun3->GetYaxis()->SetTitleSize(0.05);
  hRun3->GetYaxis()->SetLabelSize(0.045);
  hRun3->GetYaxis()->SetTitleOffset(1.2);
  hRun3->GetXaxis()->SetLabelSize(0); // Hide x labels upper pad

  // Dynamic max
  double maxVal = std::max(hRun3->GetMaximum(), hRun2->GetMaximum());
  hRun3->SetMaximum(maxVal * 1.4);

  hRun3->Draw("E");
  hRun2->Draw("E SAME");

  // Legend
  TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hRun3, "Run 3", "lp");
  leg->AddEntry(hRun2, "Run 2", "lp");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.055);
  // tex.DrawLatex(0.17, 0.88, "ALICE");
  tex.SetTextSize(0.046);
  tex.DrawLatex(0.2, 0.85,
      Form("Inclusive jets, anti-#it{k}_{T}, #it{R} = %.1f",
           arrayRadius[iRadius]));
  // tex.DrawLatex(0.17, 0.74, "pp,  #sqrt{#it{s}} = 5.02 TeV");
  tex.DrawLatex(0.2, 0.78, "|#eta_{jet}| < 0.5");

  // Optional log scale
  pad1->SetLogy();

  // =====================
  // Lower pad
  // =====================
  pad2->cd();

  TH1D *hRatio = (TH1D*)hRun3->Clone("hRatio");
  hRatio->Divide(hRun2);

  hRatio->SetTitle("");
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.9);
  hRatio->SetLineWidth(2);

  hRatio->GetYaxis()->SetTitle("Run3 / Run2");
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.10);
  hRatio->GetYaxis()->SetLabelSize(0.08);
  hRatio->GetYaxis()->SetTitleOffset(0.5);

  hRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetTitleOffset(1.0);

  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);

  hRatio->Draw("E");

  // Reference line at ratio = 1
  TLine *line = new TLine(
      hRatio->GetXaxis()->GetXmin(), 1.0,
      hRatio->GetXaxis()->GetXmax(), 1.0
  );
  line->SetLineStyle(2);
  line->SetLineColor(kRed);
  line->Draw("SAME");

  

  c->cd();
  c->Update();
  
}

void Run2vsRun2_comparison_wptLead(int iDataset, int iRadius){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "X");
  gStyle->SetTickLength(0.02, "Y");
  gStyle->SetEndErrorSize(0);   // no end caps on stat bars

  TFile* fRun2wcut = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02_LeadPtCut5.root", "READ");
  TH1D* hRun2wcut      = (TH1D*)fRun2wcut->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1");
  TH1D* hRun2_stat_wcut = (TH1D*)fRun2wcut->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e1");
  TH1D* hRun2_sys_wcut  = (TH1D*)fRun2wcut->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e2");

  TFile* fRun2 = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02.root", "READ");
  TH1D* hRun2      = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1");
  TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets in pp 5.02 TeV/Hist1D_y1_e1");

  int binStart = hRun2->FindBin(10.0);
  int binEnd   = hRun2->FindBin(100.0)-1;

  // *** Transfer Run-2 stat errors: contents of hRun2_stat -> errors of hRun2 ***
  // This must happen before any ratio or TGraphErrors construction
  for (int i = 1; i <= hRun2->GetNbinsX(); i++) {
      hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));
      hRun2wcut->SetBinError(i, hRun2_stat_wcut->GetBinContent(i));
  }
  RewriteHistogramRange(hRun2, binStart, binEnd);
  RewriteHistogramRange(hRun2wcut, binStart, binEnd);

  TCanvas *c = new TCanvas("c", "Run2 w and wo pt cut", 800, 900);

  // Define pads
  TPad *pad1 = new TPad("pad1", "upper pad", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "lower pad", 0, 0.0, 1, 0.3);

  // Pad styling
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);
  pad1->SetTicks();

  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);
  pad2->SetTicks();

  pad1->Draw();
  pad2->Draw();

  // =====================
  // Upper pad
  // =====================
  pad1->cd();

  hRun2wcut->SetLineColor(kBlue+1);
  hRun2wcut->SetMarkerColor(kBlue+1);
  hRun2wcut->SetMarkerStyle(20);
  hRun2wcut->SetMarkerSize(1.0);
  hRun2wcut->SetLineWidth(2);

  hRun2->SetLineColor(kRed+1);
  hRun2->SetMarkerColor(kRed+1);
  hRun2->SetMarkerStyle(24);
  hRun2->SetMarkerSize(1.0);
  hRun2->SetLineWidth(2);

  // Set axis titles
  hRun2wcut->SetTitle("");
  hRun2wcut->GetYaxis()->SetTitle("d^{2}#sigma_{jet} / d#it{p}_{T} d#eta  (mb GeV^{-1} #it{c})");
  hRun2wcut->GetYaxis()->SetTitleSize(0.05);
  hRun2wcut->GetYaxis()->SetLabelSize(0.045);
  hRun2wcut->GetYaxis()->SetTitleOffset(1.2);
  hRun2wcut->GetXaxis()->SetLabelSize(0); // Hide x labels upper pad

  // Dynamic max
  double maxVal = std::max(hRun2wcut->GetMaximum(), hRun2->GetMaximum());
  hRun2wcut->SetMaximum(maxVal * 1.4);

  hRun2wcut->Draw("E");
  hRun2->Draw("E SAME");

  // Legend
  TLegend *leg = new TLegend(0.50, 0.72, 0.88, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);    // Standard Helvetica font
  leg->SetTextSize(0.035); // Slightly larger than default
  leg->AddEntry(hRun2wcut, "pp, #sqrt{s} = 5.02 TeV, p_{T, lead} > 5 GeV/c", "lp");
  leg->AddEntry(hRun2, "pp, #sqrt{s} = 5.02 TeV", "lp");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.055);
  // tex.DrawLatex(0.17, 0.88, "ALICE");
  tex.SetTextSize(0.046);
  tex.DrawLatex(0.2, 0.85,
      Form("Inclusive jets, anti-#it{k}_{T}, #it{R} = %.1f",
           arrayRadius[iRadius]));
  // tex.DrawLatex(0.17, 0.74, "pp,  #sqrt{#it{s}} = 5.02 TeV");
  tex.DrawLatex(0.2, 0.78, "|#eta_{jet}| < 0.5");

  // Optional log scale
  pad1->SetLogy();

  // =====================
  // Lower pad
  // =====================
  pad2->cd();

  TH1D *hRatio = (TH1D*)hRun2wcut->Clone("hRatio");
  // hRatio->Divide(hRun2);
  hRatio->Divide(hRun2wcut, hRun2, 1, 1, "B");

  hRatio->SetTitle("");
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.9);
  hRatio->SetLineWidth(2);

  hRatio->GetYaxis()->SetTitle("Run2 w cut / Run2");
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.10);
  hRatio->GetYaxis()->SetLabelSize(0.08);
  hRatio->GetYaxis()->SetTitleOffset(0.5);

  hRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetTitleOffset(1.0);

  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);

  hRatio->Draw("E");

  // Reference line at ratio = 1
  TLine *line = new TLine(
      hRatio->GetXaxis()->GetXmin(), 1.0,
      hRatio->GetXaxis()->GetXmax(), 1.0
  );
  line->SetLineStyle(2);
  line->SetLineColor(kRed);
  line->Draw("SAME");

  

  c->cd();
  c->Update();
  
}

// void Draw_Sigma_spectrum_comparison(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
//   bool splitTestControlMC = true;

//   TH1D* H1D_jetPt_unfolded;
//   TString partialUniqueSpecifier;
//   int unfoldParameter;
//   partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
//   TH1D* measuredInput;

//   Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
//   unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).first; // 1/N_ev d^2N/dpTdeta where N_ev is bin8

//   const double sigma_vdm_mb = 46.46;   // mb
//   const double vdmRun3      = 0.045;   // 4.5% relative lumi unc (Run 3)
//   const double sigma_O2 =  50.3; //mb line 1685 https://github.com/AliceO2Group/O2Physics/blob/f4ec4e509734228d21fc11a707541c629186c7b5/Common/Tools/EventSelectionModule.h

//   TH1* hCounter = (TH1*)file_O2Analysis_list[iDataset]->Get("jet-luminosity-calculator/counter");
//   double nEventsSel8 = hCounter->GetBinContent(8);
//   double nEventsBin4 = hCounter->GetBinContent(4);
//   double nBCTVX_Bin2 = hCounter->GetBinContent(2);
//   double nColTVXBin6 = hCounter->GetBinContent(6);
//   cout << "nEventsSel8 = " << nEventsSel8 << ", nEventsBin4 = " << nEventsBin4 << endl;
//   cout << "nColTVXBin6 = " << nColTVXBin6 << ", nBCTVX_Bin2 = " << nBCTVX_Bin2 << endl;
//   // cout << "Eff_TVX = " << nEventsSel8/ nColTVXBin6 << endl;
//   // cout << "Lumi = " << nBCTVX_Bin2/ sigma_vdm_mb << endl;
//   // const double scalingFactor = sigma_vdm_mb *nColTVXBin6 /(nBCTVX_Bin2); // Nima suggestion
//   // cout << "scalingFactor = " << scalingFactor << endl;
//   // H1D_jetPt_unfolded->Scale(scalingFactor);

//   TH1D* TVX_eff = Get_TVX_Eff(iDataset, iRadius);
//   double z_vtx_eff = Get_zVertex_reconstruction_efficiency(iDataset, iRadius);
//   double SBP_eff = Get_SBP_Eff(iDataset, iRadius);
//   // double lumi = nBCTVX_Bin2 / sigma_vdm_mb ;
//   double lumi = nEventsBin4 / sigma_vdm_mb ;
//   cout << "Lumi = " << lumi << endl;
//   const double scalingFactor = nEventsSel8 /(lumi * z_vtx_eff * SBP_eff ); 
//   H1D_jetPt_unfolded->Divide(TVX_eff);
//   H1D_jetPt_unfolded->Scale(scalingFactor);


//   TString* pdfName = new TString("XSextion_jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
//   TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
//   TString* texJetXsection_d2Sigmadptdeta = new TString("d^{2}#sigma_{jet}/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
//   Draw_TH1_Histogram(H1D_jetPt_unfolded, textContext, pdfName, texPtX, texJetXsection_d2Sigmadptdeta, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

//   double newbins[] = {10., 12, 14, 16, 18, 20., 25, 30., 40., 50., 60., 70., 85., 100., 140.}; 
//   int n_newbins = sizeof(newbins)/sizeof(newbins[0]) - 1; // = 10

//   TH1D* H1D_jetPt_rebinned = ReweightedRebin(H1D_jetPt_unfolded, "H1D_jetPt_rebinned", n_newbins, newbins);
//   TString* pdfName_Rebinned = new TString("Rebinned_XSextion_jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
//   Draw_TH1_Histogram(H1D_jetPt_rebinned, textContext, pdfName_Rebinned, texPtX, texJetXsection_d2Sigmadptdeta, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");


//   // ---------------------------------------------------------------------------
//   // 2.  RUN-2 CROSS SECTION
//   //     hRun2_stat has stat errors as bin CONTENTS -> copy to bin errors
//   //     hRun2_sys  has systematic errors as ABSOLUTE bin contents
//   // ---------------------------------------------------------------------------
//   TFile* fRun2 = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02_LeadPtCut5.root", "READ");
//   TH1D* hRun2      = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1");
//   TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e1");
//   TH1D* hRun2_sys  = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e2");

//   // *** Transfer Run-2 stat errors: contents of hRun2_stat -> errors of hRun2 ***
//   // This must happen before any ratio or TGraphErrors construction
//   for (int i = 1; i <= hRun2->GetNbinsX(); i++) {
//       hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));
//   }


// }

void Draw_Sigma_spectrum_comparison(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  bool splitTestControlMC = true;

  TH1D* H1D_jetPt_unfolded;
  TString partialUniqueSpecifier;
  int unfoldParameter;
  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  TH1D* measuredInput;

  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
  unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).first;

  const double sigma_vdm_mb = 46.46;   // mb
  const double vdmRun3      = 0.045;   // 4.5% relative lumi unc (Run 3)
  const double sigma_O2     = 50.3;    // mb

  TH1* hCounter = (TH1*)file_O2Analysis_list[iDataset]->Get("jet-luminosity-calculator/counter");
  double nEventsSel8 = hCounter->GetBinContent(8);
  double nEventsBin4 = hCounter->GetBinContent(4);
  double nBCTVX_Bin2 = hCounter->GetBinContent(2);
  double nColTVXBin6 = hCounter->GetBinContent(6);
  cout << "nEventsSel8 = " << nEventsSel8 << ", nEventsBin4 = " << nEventsBin4 << endl;
  cout << "nColTVXBin6 = " << nColTVXBin6 << ", nBCTVX_Bin2 = " << nBCTVX_Bin2 << endl;

  TH1D* TVX_eff   = Get_TVX_Eff(iDataset, iRadius);
  double z_vtx_eff = Get_zVertex_reconstruction_efficiency(iDataset, iRadius);
  double SBP_eff   = Get_SBP_Eff(iDataset, iRadius);

  // Luminosity from MB-trigger counter (bin 4)
  double lumi = nEventsBin4 / sigma_vdm_mb;
  cout << "Lumi = " << lumi << " mb^-1" << endl;

  // Scale unfolded spectrum to cross section: d^2sigma/dpT deta  [mb/(GeV/c)]
  // H1D_jetPt_unfolded is (1/N_ev) d^2N/dpT deta  (N_ev = bin-8 count)
  // d^2sigma/dpT deta = (N_ev / (lumi * eps_zvtx * eps_SBP * eps_TVX(pT))) * (1/N_ev) d^2N/dpT deta
  //                   = (N_ev / (lumi * eps_zvtx * eps_SBP)) * H1D / eps_TVX(pT)
  const double scalingFactor = nEventsSel8 / (lumi * z_vtx_eff * SBP_eff);
  H1D_jetPt_unfolded->Divide(TVX_eff);
  H1D_jetPt_unfolded->Scale(scalingFactor);

  // ---------------------------------------------------------------------------
  // 1.  REBIN RUN-3 TO ANALYSIS BINNING
  // ---------------------------------------------------------------------------
  double newbins[] = {10., 12., 14., 16., 18., 20., 25., 30., 40., 50., 60., 70., 85., 100., 140.};
  int n_newbins = sizeof(newbins)/sizeof(newbins[0]) - 1; // 14 bins

  TH1D* hRun3 = ReweightedRebin(H1D_jetPt_unfolded, "hRun3_rebinned", n_newbins, newbins);
  hRun3->SetDirectory(0);

  // ---------------------------------------------------------------------------
  // 2.  RUN-2 CROSS SECTION
  // ---------------------------------------------------------------------------
  TFile* fRun2 = TFile::Open("../Datasets/Run2_CrossSection/Run2_Jets_in_pp_5.02_LeadPtCut5.root", "READ");
  if (!fRun2 || fRun2->IsZombie()) {
    cerr << "ERROR: Cannot open Run-2 file!" << endl;
    return;
  }
  TH1D* hRun2      = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1");
  TH1D* hRun2_stat = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e1");
  TH1D* hRun2_sys  = (TH1D*)fRun2->Get("Jets with leading p_{T} 5 GeV-c in pp 5.02 TeV/Hist1D_y1_e2");
  if (!hRun2 || !hRun2_stat) {
    cerr << "ERROR: Cannot retrieve Run-2 histograms!" << endl;
    return;
  }
  hRun2->SetDirectory(0);
  hRun2_stat->SetDirectory(0);
  fRun2->Close();

  // Transfer Run-2 statistical errors (stored as bin contents in _e1)
  for (int i = 1; i <= hRun2->GetNbinsX(); ++i)
    hRun2->SetBinError(i, hRun2_stat->GetBinContent(i));

  // ---------------------------------------------------------------------------
  // 3.  FIND COMMON BIN EDGES FOR THE RATIO
  //     Only keep bins whose low AND high edges exist in BOTH histograms
  // ---------------------------------------------------------------------------
  auto GetEdges = [](TH1D* h) -> std::vector<double> {
    std::vector<double> edges;
    for (int i = 1; i <= h->GetNbinsX(); ++i) edges.push_back(h->GetXaxis()->GetBinLowEdge(i));
    edges.push_back(h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
    return edges;
  };

  auto ApproxEqual = [](double a, double b, double tol = 1e-4) -> bool {
    return std::fabs(a - b) < tol * (std::fabs(a) + std::fabs(b) + 1e-10);
  };

  std::vector<double> edgesRun3 = GetEdges(hRun3);
  std::vector<double> edgesRun2 = GetEdges(hRun2);

  // Common edges = edges present in both (within tolerance)
  std::vector<double> commonEdges;
  for (double e3 : edgesRun3)
    for (double e2 : edgesRun2)
      if (ApproxEqual(e3, e2)) { commonEdges.push_back(e3); break; }

  if (commonEdges.size() < 2) {
    cerr << "ERROR: No overlapping bins between Run-2 and Run-3!" << endl;
    return;
  }
  std::sort(commonEdges.begin(), commonEdges.end());

  int nCommon = (int)commonEdges.size() - 1;
  cout << "Common bins for ratio: " << nCommon << " bins, pT in ["
       << commonEdges.front() << ", " << commonEdges.back() << "] GeV/c" << endl;

  // Rebin both to common edges (ReweightedRebin preserves cross section)
  TH1D* hRun3_common = ReweightedRebin(hRun3, "hRun3_common", nCommon, commonEdges.data());
  TH1D* hRun2_common = ReweightedRebin(hRun2, "hRun2_common", nCommon, commonEdges.data());
  hRun3_common->SetDirectory(0);
  hRun2_common->SetDirectory(0);

  // ---------------------------------------------------------------------------
  // 4.  RATIO  Run3 / Run2  (stat errors propagated in quadrature)
  // ---------------------------------------------------------------------------
  TH1D* hRatio = (TH1D*)hRun3_common->Clone("hRatio_Run3overRun2");
  hRatio->SetDirectory(0);
  hRatio->Divide(hRun3_common, hRun2_common, 1., 1., "B"); // "B" = binomial-error option — use "" for uncorr.
  // For uncorrelated stat errors use standard Divide without "B":
  // hRatio->Divide(hRun3_common, hRun2_common);

  // ---------------------------------------------------------------------------
  // 5.  STYLING
  // ---------------------------------------------------------------------------
  // Run-3
  hRun3->SetMarkerStyle(kFullCircle);
  hRun3->SetMarkerColor(kRed+1);
  hRun3->SetLineColor(kRed+1);
  hRun3->SetMarkerSize(1.1);
  // Run-2
  hRun2->SetMarkerStyle(kOpenSquare);
  hRun2->SetMarkerColor(kBlue+1);
  hRun2->SetLineColor(kBlue+1);
  hRun2->SetMarkerSize(1.1);
  // Ratio
  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(1.0);

  // ---------------------------------------------------------------------------
  // 6.  CANVAS WITH TWO PADS
  // ---------------------------------------------------------------------------
  TString canvasName = "cXsection_Run3vsRun2_"+Datasets[iDataset]+Form("_R%.1f", arrayRadius[iRadius]);
  TCanvas* c = new TCanvas(canvasName, canvasName, 800, 900);
  c->SetFillStyle(0);

  // Upper pad  (80 % of height)
  TPad* padUp = new TPad("padUp", "padUp", 0., 0.30, 1., 1.);
  padUp->SetBottomMargin(0.015);
  padUp->SetTopMargin(0.08);
  padUp->SetLeftMargin(0.15);
  padUp->SetRightMargin(0.05);
  padUp->SetLogy();
  padUp->Draw();

  // Lower pad  (30 % of height)
  TPad* padDn = new TPad("padDn", "padDn", 0., 0.00, 1., 0.30);
  padDn->SetTopMargin(0.015);
  padDn->SetBottomMargin(0.35);
  padDn->SetLeftMargin(0.15);
  padDn->SetRightMargin(0.05);
  padDn->Draw();

  // ---- upper pad ----
  padUp->cd();

  // Determine y-axis range from both histograms
  double yMax = std::max(hRun3->GetMaximum(), hRun2->GetMaximum()) * 5.;
  double yMin = 1e-8; // adjust to data

  hRun3->GetYaxis()->SetTitle("d^{2}#sigma_{jet}/d#it{p}_{T}d#it{#eta} [mb (GeV/#it{c})^{-1}]");
  hRun3->GetXaxis()->SetLabelSize(0.);
  hRun3->GetYaxis()->SetTitleSize(0.055);
  hRun3->GetYaxis()->SetTitleOffset(1.2);
  hRun3->GetYaxis()->SetLabelSize(0.05);
  hRun3->SetStats(0);
  hRun3->GetXaxis()->SetRangeUser(commonEdges.front(), commonEdges.back());
  hRun3->SetMinimum(yMin);
  hRun3->SetMaximum(yMax);

  hRun3->Draw("E1 X0");    // stat errors, no x error bars
  hRun2->Draw("E1 X0 SAME");

  // Legend
  TLegend* leg = new TLegend(0.55, 0.65, 0.90, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.048);
  leg->AddEntry(hRun3, Form("Run 3  pp #sqrt{s}=5.02 TeV  R=%.1f", arrayRadius[iRadius]), "lep");
  leg->AddEntry(hRun2, "Run 2  pp #sqrt{s}=5.02 TeV  (ALICE)", "lep");
  leg->Draw();

  // Context label
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.048);
  latex.DrawLatex(0.18, 0.88, "Anti-#it{k}_{T} jets, |#it{#eta}_{jet}| < 0.7 - R");
  latex.DrawLatex(0.18, 0.82, "Stat. errors only");

  // ---- lower pad ----
  padDn->cd();

  // Unity line
  TLine* line1 = new TLine(commonEdges.front(), 1., commonEdges.back(), 1.);
  line1->SetLineStyle(2);
  line1->SetLineColor(kGray+2);

  hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hRatio->GetYaxis()->SetTitle("Run 3 / Run 2");
  hRatio->GetXaxis()->SetTitleSize(0.13);
  hRatio->GetXaxis()->SetLabelSize(0.11);
  hRatio->GetYaxis()->SetTitleSize(0.11);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetTitleOffset(0.55);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetXaxis()->SetRangeUser(commonEdges.front(), commonEdges.back());
  hRatio->SetMinimum(0.0);
  hRatio->SetMaximum(2.5);
  hRatio->SetStats(0);

  hRatio->Draw("E1 X0");
  line1->Draw("SAME");

  // ---- save ----
  // TString pdfOut = "XSextion_Run3vsRun2_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]
  //                  +"_"+Datasets[iDataset]+DatasetsNames[iDataset]
  //                  +Form("_R=%.1f", arrayRadius[iRadius])+".pdf";
  // // c->SaveAs(pdfOut);
  // // cout << "Saved: " << pdfOut << endl;
}


// void Run3vsRun3_comparison_wptLead(int iDataset, int iRadius) {
//   gStyle->SetOptStat(0);
//   gStyle->SetOptTitle(0);
//   gStyle->SetPadTickX(1);
//   gStyle->SetPadTickY(1);
//   gStyle->SetTickLength(0.02, "X");
//   gStyle->SetTickLength(0.02, "Y");
//   gStyle->SetEndErrorSize(0);

//   // ---------------------------------------------------------------------------
//   // 1.  LOAD HISTOGRAMS
//   // ---------------------------------------------------------------------------
//   // Run3 without pT-lead cut
//   TFile* fRun3_woCut = TFile::Open(
//     "/Users/tabikh/Documents/PhdWork/MyWork/20260530_R02_673480/Xsection/Xsection_woptLead.root",
//     "READ");
//   if (!fRun3_woCut || fRun3_woCut->IsZombie()) {
//     cerr << "ERROR: Cannot open Run3 (no cut) file!" << endl; return;
//   }
//   TH1D* hRun3_woCut = (TH1D*)fRun3_woCut->Get("H1D_XSection_Pt_Unfolded_Rebinned_w_MB");
//   if (!hRun3_woCut) {
//     cerr << "ERROR: Cannot retrieve H1D_XSection_Pt_Unfolded_Rebinned_w_MB!" << endl; return;
//   }
//   hRun3_woCut->SetDirectory(0);
//   fRun3_woCut->Close();

//   // Run3 with pT-lead cut > 5 GeV
//   TFile* fRun3_wCut = TFile::Open(
//     "/Users/tabikh/Documents/PhdWork/MyWork/20260504_UnfDD_667972_w_MB_669669/XSection_pp536_R02_ptLeadcut5.root",
//     "READ");
//   if (!fRun3_wCut || fRun3_wCut->IsZombie()) {
//     cerr << "ERROR: Cannot open Run3 (with cut) file!" << endl; return;
//   }
//   TH1D* hRun3_wCut = (TH1D*)fRun3_wCut->Get("H1D_XSection_pp536_R02_ptLeadcut5");
//   if (!hRun3_wCut) {
//     cerr << "ERROR: Cannot retrieve H1D_XSection_pp536_R02_ptLeadcut5!" << endl; return;
//   }
//   hRun3_wCut->SetDirectory(0);
//   fRun3_wCut->Close();

//   // ---------------------------------------------------------------------------
//   // 2.  FIND COMMON BIN EDGES  (numerator = subset of denominator)
//   //     The cut spectrum (wCut) is the numerator — it is a strict subset of the
//   //     no-cut spectrum (woCut), so binomial error propagation is valid.
//   // ---------------------------------------------------------------------------
//   auto GetEdges = [](const TH1D* h) -> std::vector<double> {
//     std::vector<double> edges;
//     for (int i = 1; i <= h->GetNbinsX(); ++i)
//       edges.push_back(h->GetXaxis()->GetBinLowEdge(i));
//     edges.push_back(h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
//     return edges;
//   };

//   auto ApproxEqual = [](double a, double b, double tol = 1e-4) -> bool {
//     return std::fabs(a - b) < tol * (std::fabs(a) + std::fabs(b) + 1e-10);
//   };

//   std::vector<double> edgesWoCut = GetEdges(hRun3_woCut);
//   std::vector<double> edgesWCut  = GetEdges(hRun3_wCut);

//   std::vector<double> commonEdges;
//   for (double ew : edgesWCut)
//     for (double eo : edgesWoCut)
//       if (ApproxEqual(ew, eo)) { commonEdges.push_back(ew); break; }

//   std::sort(commonEdges.begin(), commonEdges.end());
//   // Remove duplicate edges (safety)
//   commonEdges.erase(std::unique(commonEdges.begin(), commonEdges.end(),
//     [&](double a, double b){ return ApproxEqual(a, b); }), commonEdges.end());

//   if (commonEdges.size() < 2) {
//     cerr << "ERROR: No overlapping bin edges between the two Run-3 spectra!" << endl;
//     return;
//   }

//   int nCommon = (int)commonEdges.size() - 1;
//   cout << "Common bins: " << nCommon
//        << "  pT in [" << commonEdges.front() << ", " << commonEdges.back() << "] GeV/c" << endl;

//   // Print edge-by-edge comparison for sanity check
//   cout << "--- Bin edge check ---" << endl;
//   cout << Form("  %-12s  %-12s  %s", "woCut edge", "wCut edge", "match?") << endl;
//   {
//     int j = 0;
//     for (int i = 0; i < (int)edgesWoCut.size() && j < (int)edgesWCut.size(); ) {
//       if (ApproxEqual(edgesWoCut[i], edgesWCut[j])) {
//         cout << Form("  %-12.4f  %-12.4f  OK", edgesWoCut[i], edgesWCut[j]) << endl;
//         ++i; ++j;
//       } else if (edgesWoCut[i] < edgesWCut[j] - 1e-6) {
//         cout << Form("  %-12.4f  %-12s  (woCut only)", edgesWoCut[i], "-") << endl;
//         ++i;
//       } else {
//         cout << Form("  %-12s  %-12.4f  (wCut only)", "-", edgesWCut[j]) << endl;
//         ++j;
//       }
//     }
//   }
//   cout << "----------------------" << endl;

//   // Rebin both spectra to the common grid
//   // (ReweightedRebin must be defined in your framework and preserve the cross section)
//   TH1D* hDenom = ReweightedRebin(hRun3_woCut, "hDenom_common", nCommon, commonEdges.data());
//   TH1D* hNum   = ReweightedRebin(hRun3_wCut,  "hNum_common",   nCommon, commonEdges.data());
//   hDenom->SetDirectory(0);
//   hNum  ->SetDirectory(0);

//   // ---------------------------------------------------------------------------
//   // 3.  RATIO  wCut / woCut  — binomial error propagation
//   //     epsilon_i = N_cut,i / N_nocut,i  (both are cross sections in mb/(GeV/c),
//   //     but since they share the same luminosity and bin widths the ratio cancels
//   //     those factors and we get the bin-by-bin efficiency of the leading-pT cut)
//   //
//   //     Binomial variance:  Var(eps) = eps*(1-eps)/N_trials
//   //     In terms of bin counts (before luminosity scaling):
//   //       sigma_eps = sqrt( eps*(1-eps) / N_denom )
//   //     Here we only have the cross sections, so we propagate via ROOT's "B" option
//   //     which applies the binomial formula using the contents as counts — valid
//   //     because numerator IS a subset of denominator (same events, extra cut applied).
//   // ---------------------------------------------------------------------------
//   TH1D* hRatio = (TH1D*)hNum->Clone("hRatio_wCut_over_woCut");
//   hRatio->SetDirectory(0);
//   hRatio->Divide(hNum, hDenom, 1., 1., "B");  // "B" = binomial error propagation

//   // Sanity print
//   cout << "--- Ratio bin-by-bin ---" << endl;
//   for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
//     double lo  = hRatio->GetXaxis()->GetBinLowEdge(i);
//     double hi  = hRatio->GetXaxis()->GetBinUpEdge(i);
//     double val = hRatio->GetBinContent(i);
//     double err = hRatio->GetBinError(i);
//     cout << Form("  pT [%5.1f, %5.1f]  ratio = %.4f +/- %.4f", lo, hi, val, err) << endl;
//   }
//   cout << "------------------------" << endl;

//   // ---------------------------------------------------------------------------
//   // 4.  STYLING
//   // ---------------------------------------------------------------------------
//   // Run3 without cut  (denominator / reference)
//   hRun3_woCut->SetMarkerStyle(kFullCircle);
//   hRun3_woCut->SetMarkerColor(kBlue+1);
//   hRun3_woCut->SetLineColor(kBlue+1);
//   hRun3_woCut->SetMarkerSize(1.1);

//   // Run3 with pT-lead cut
//   hRun3_wCut->SetMarkerStyle(kFullSquare);
//   hRun3_wCut->SetMarkerColor(kRed+1);
//   hRun3_wCut->SetLineColor(kRed+1);
//   hRun3_wCut->SetMarkerSize(1.1);

//   // Ratio
//   hRatio->SetMarkerStyle(kFullCircle);
//   hRatio->SetMarkerColor(kBlack);
//   hRatio->SetLineColor(kBlack);
//   hRatio->SetMarkerSize(1.0);

//   // ---------------------------------------------------------------------------
//   // 5.  CANVAS + PADS
//   // ---------------------------------------------------------------------------
//   TString canvasName = Form("cRun3_woCut_vs_wCut_R%.1f", arrayRadius[iRadius]);
//   TCanvas* c = new TCanvas(canvasName, canvasName, 800, 900);
//   c->SetFillStyle(0);

//   // Upper pad  (70 % height)
//   TPad* padUp = new TPad("padUp", "padUp", 0., 0.30, 1., 1.);
//   padUp->SetBottomMargin(0.015);
//   padUp->SetTopMargin(0.08);
//   padUp->SetLeftMargin(0.16);
//   padUp->SetRightMargin(0.05);
//   padUp->SetLogy();
//   padUp->SetTickx(1); padUp->SetTicky(1);
//   padUp->Draw();

//   // Lower pad  (30 % height)
//   TPad* padDn = new TPad("padDn", "padDn", 0., 0.00, 1., 0.30);
//   padDn->SetTopMargin(0.015);
//   padDn->SetBottomMargin(0.38);
//   padDn->SetLeftMargin(0.16);
//   padDn->SetRightMargin(0.05);
//   padDn->SetTickx(1); padDn->SetTicky(1);
//   padDn->Draw();

//   // ---------------------------------------------------------------------------
//   // 6.  UPPER PAD — both spectra
//   // ---------------------------------------------------------------------------
//   padUp->cd();

//   double xLo = commonEdges.front();
//   double xHi = commonEdges.back();

//   // Use the no-cut spectrum as the frame (it covers the full common range)
//   hRun3_woCut->GetXaxis()->SetRangeUser(xLo, xHi);
//   hRun3_woCut->GetXaxis()->SetLabelSize(0.);    // hide x labels on upper pad
//   hRun3_woCut->GetYaxis()->SetTitle("d^{2}#sigma_{jet}/d#it{p}_{T}d#it{#eta} mb (GeV/#it{c})^{-1}");
//   hRun3_woCut->GetYaxis()->SetTitleSize(0.055);
//   hRun3_woCut->GetYaxis()->SetTitleOffset(1.30);
//   hRun3_woCut->GetYaxis()->SetLabelSize(0.050);
//   hRun3_woCut->SetStats(0);

//   double yMax = std::max(hRun3_woCut->GetMaximum(), hRun3_wCut->GetMaximum()) * 8.;
//   double yMin = std::min(hRun3_woCut->GetMinimum(0.), hRun3_wCut->GetMinimum(0.)) * 0.1;
//   hRun3_woCut->SetMinimum(yMin > 0. ? yMin : 1e-9);
//   hRun3_woCut->SetMaximum(yMax);

//   hRun3_woCut->Draw("E1 X0");
//   hRun3_wCut ->Draw("E1 X0 SAME");

//   TLegend *leg = new TLegend(0.50, 0.72, 0.88, 0.89);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->SetTextFont(42);    // Standard Helvetica font
//   leg->SetTextSize(0.035); // Slightly larger than default
//   leg->AddEntry(hRun3_woCut, "pp #sqrt{#it{s}} = 5.36 TeV ","lep");
//   leg->AddEntry(hRun3_wCut,
//               "pp #sqrt{#it{s}} = 5.36 TeV, #it{p}_{T}^{lead} > 5 GeV/#it{c}",
//               "lep");
//   leg->Draw();

//   // Info label
//   TLatex ltx;
//   ltx.SetNDC();
//   ltx.SetTextSize(0.048);
//   ltx.DrawLatex(0.18, 0.90, "Anti-#it{k}_{T} charged jets,  |#it{#eta}_{jet}| < 0.7 #minus #it{R}");
//   ltx.SetTextSize(0.042);
//   ltx.DrawLatex(0.18, 0.83, "Stat. errors only");

//   // ---------------------------------------------------------------------------
//   // 7.  LOWER PAD — ratio  wCut / woCut
//   // ---------------------------------------------------------------------------
//   padDn->cd();

//   // Unity and guide lines
//   TLine* lineUnity = new TLine(xLo, 1., xHi, 1.);
//   lineUnity->SetLineStyle(2);
//   lineUnity->SetLineColor(kGray+2);
//   lineUnity->SetLineWidth(1);

//   hRatio->GetXaxis()->SetRangeUser(xLo, xHi);
//   hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
//   hRatio->GetYaxis()->SetTitle("Run3 w cut/Run3");

//   // Scale axes for the compressed pad height
//   const double scaleFactor = (1. - 0.30) / 0.30;   // upper_height / lower_height ~ 2.33
//   hRatio->GetXaxis()->SetTitleSize  (0.055 * scaleFactor);
//   hRatio->GetXaxis()->SetLabelSize  (0.048 * scaleFactor);
//   hRatio->GetXaxis()->SetTitleOffset(0.85);
//   hRatio->GetYaxis()->SetTitleSize  (0.048 * scaleFactor);
//   hRatio->GetYaxis()->SetLabelSize  (0.044 * scaleFactor);
//   hRatio->GetYaxis()->SetTitleOffset(0.48);
//   hRatio->GetYaxis()->SetNdivisions(504);
//   hRatio->GetXaxis()->SetTickLength (0.06);
//   hRatio->GetYaxis()->SetTickLength (0.03);

//   hRatio->SetMinimum(0.0);
//   hRatio->SetMaximum(1.4);
//   hRatio->SetStats(0);

//   hRatio->Draw("E1 X0");
//   lineUnity->Draw("SAME");

//   // Binomial error note
//   TLatex ltxR;
//   ltxR.SetNDC();
//   ltxR.SetTextSize(0.048 * scaleFactor * 0.85);
//   ltxR.DrawLatex(0.18, 0.82, "Binomial stat. errors");

//   // ---------------------------------------------------------------------------
//   // // 8.  SAVE
//   // // ---------------------------------------------------------------------------
//   // TString pdfOut = Form("Run3_wCut_vs_woCut_R%.1f.pdf", arrayRadius[iRadius]);
//   // c->SaveAs(pdfOut);
//   // cout << "Saved: " << pdfOut << endl;
// }

void Run3vsRun3_comparison_wptLead(int iDataset, int iRadius) {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "X");
  gStyle->SetTickLength(0.02, "Y");
  gStyle->SetEndErrorSize(0);

  // ---------------------------------------------------------------------------
  // 1.  LOAD HISTOGRAMS
  // ---------------------------------------------------------------------------
  // Run3 without pT-lead cut
  TFile* fRun3_woCut = TFile::Open(
    "/Users/tabikh/Documents/PhdWork/MyWork/20260530_R02_673480/Xsection/Xsection_woptLead.root",
    "READ");
  if (!fRun3_woCut || fRun3_woCut->IsZombie()) {
    cerr << "ERROR: Cannot open Run3 (no cut) file!" << endl; return;
  }
  TH1D* hRun3_woCut = (TH1D*)fRun3_woCut->Get("H1D_XSection_Pt_Unfolded_Rebinned_w_MB");
  if (!hRun3_woCut) {
    cerr << "ERROR: Cannot retrieve H1D_XSection_Pt_Unfolded_Rebinned_w_MB!" << endl; return;
  }
  hRun3_woCut->SetDirectory(0);
  fRun3_woCut->Close();

  // Run3 with pT-lead cut > 5 GeV
  TFile* fRun3_wCut = TFile::Open(
    "/Users/tabikh/Documents/PhdWork/MyWork/20260504_UnfDD_667972_w_MB_669669/XSection_pp536_R02_ptLeadcut5.root",
    "READ");
  if (!fRun3_wCut || fRun3_wCut->IsZombie()) {
    cerr << "ERROR: Cannot open Run3 (with cut) file!" << endl; return;
  }
  TH1D* hRun3_wCut = (TH1D*)fRun3_wCut->Get("H1D_XSection_pp536_R02_ptLeadcut5");
  if (!hRun3_wCut) {
    cerr << "ERROR: Cannot retrieve H1D_XSection_pp536_R02_ptLeadcut5!" << endl; return;
  }
  hRun3_wCut->SetDirectory(0);
  fRun3_wCut->Close();

  // -------- POWHEG: no pT-lead cut ----------
  TFile* fPow_woCut = TFile::Open(
    "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_woPtLeadCut/POWHEG_Uncertainties.root",
    "READ");
  if (!fPow_woCut || fPow_woCut->IsZombie()) {
    cerr << "ERROR: Cannot open POWHEG (no cut) file!" << endl; return;
  }
  TH1F* hPow_woCut    = (TH1F*)fPow_woCut->Get("Central_Inclusive_R02");
  TH1F* hPow_woCutUnc = (TH1F*)fPow_woCut->Get("hTotalUnc_Inclusive_R02");
  if (!hPow_woCut || !hPow_woCutUnc) {
    cerr << "ERROR: POWHEG (no cut) histograms missing!" << endl; return;
  }
  hPow_woCut    = (TH1F*)hPow_woCut->Clone("hPow_woCut");    hPow_woCut->SetDirectory(0);
  hPow_woCutUnc = (TH1F*)hPow_woCutUnc->Clone("hPow_woCutUnc"); hPow_woCutUnc->SetDirectory(0);
  fPow_woCut->Close();

  // -------- POWHEG: with pT-lead > 5 cut ----------
  TFile* fPow_wCut = TFile::Open(
    "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_ptlead5/POWHEG_Uncertainties.root",
    "READ");
  if (!fPow_wCut || fPow_wCut->IsZombie()) {
    cerr << "ERROR: Cannot open POWHEG (with cut) file!" << endl; return;
  }
  TH1F* hPow_wCut    = (TH1F*)fPow_wCut->Get("Central_Inclusive_R02");
  TH1F* hPow_wCutUnc = (TH1F*)fPow_wCut->Get("hTotalUnc_Inclusive_R02");
  if (!hPow_wCut || !hPow_wCutUnc) {
    cerr << "ERROR: POWHEG (with cut) histograms missing!" << endl; return;
  }
  hPow_wCut    = (TH1F*)hPow_wCut->Clone("hPow_wCut");    hPow_wCut->SetDirectory(0);
  hPow_wCutUnc = (TH1F*)hPow_wCutUnc->Clone("hPow_wCutUnc"); hPow_wCutUnc->SetDirectory(0);
  fPow_wCut->Close();

  // ---------------------------------------------------------------------------
  // 2.  FIND COMMON BIN EDGES  (numerator = subset of denominator)
  // ---------------------------------------------------------------------------
  auto GetEdges = [](const TH1D* h) -> std::vector<double> {
    std::vector<double> edges;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
      edges.push_back(h->GetXaxis()->GetBinLowEdge(i));
    edges.push_back(h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
    return edges;
  };

  auto ApproxEqual = [](double a, double b, double tol = 1e-4) -> bool {
    return std::fabs(a - b) < tol * (std::fabs(a) + std::fabs(b) + 1e-10);
  };

  std::vector<double> edgesWoCut = GetEdges(hRun3_woCut);
  std::vector<double> edgesWCut  = GetEdges(hRun3_wCut);

  std::vector<double> commonEdges;
  for (double ew : edgesWCut)
    for (double eo : edgesWoCut)
      if (ApproxEqual(ew, eo)) { commonEdges.push_back(ew); break; }

  std::sort(commonEdges.begin(), commonEdges.end());
  commonEdges.erase(std::unique(commonEdges.begin(), commonEdges.end(),
    [&](double a, double b){ return ApproxEqual(a, b); }), commonEdges.end());

  if (commonEdges.size() < 2) {
    cerr << "ERROR: No overlapping bin edges between the two Run-3 spectra!" << endl;
    return;
  }

  int nCommon = (int)commonEdges.size() - 1;
  cout << "Common bins: " << nCommon
       << "  pT in [" << commonEdges.front() << ", " << commonEdges.back() << "] GeV/c" << endl;

  // Rebin both spectra to the common grid
  TH1D* hDenom = ReweightedRebin(hRun3_woCut, "hDenom_common", nCommon, commonEdges.data());
  TH1D* hNum   = ReweightedRebin(hRun3_wCut,  "hNum_common",   nCommon, commonEdges.data());
  hDenom->SetDirectory(0);
  hNum  ->SetDirectory(0);

  // ---------------------------------------------------------------------------
  // 3.  RATIO  wCut / woCut  — binomial error propagation
  // ---------------------------------------------------------------------------
  TH1D* hRatio = (TH1D*)hNum->Clone("hRatio_wCut_over_woCut");
  hRatio->SetDirectory(0);
  hRatio->Divide(hNum, hDenom, 1., 1., "B");

  // ---------------------------------------------------------------------------
  // 3b. POWHEG ratio  (wCut / woCut) with uncertainty propagated
  // ---------------------------------------------------------------------------
  // Build a histogram on POWHEG's binning, evaluated bin-center by bin-center.
  // Assume both POWHEG histograms share the same binning (they're produced by
  // the same uncertainty macro). The relative uncertainties are added in
  // quadrature for the ratio (independent variations between the two POWHEG runs).
  TH1F* hPowRatio = (TH1F*)hPow_wCut->Clone("hPow_ratio_wCut_over_woCut");
  hPowRatio->SetDirectory(0);
  hPowRatio->Reset();
  int nPow = hPow_wCut->GetNbinsX();
  TGraphAsymmErrors* gPowRatioBand = new TGraphAsymmErrors();
  int kpr = 0;
  for (int ib = 1; ib <= nPow; ++ib) {
    double x   = hPow_wCut->GetBinCenter(ib);
    double w   = hPow_wCut->GetBinContent(ib);
    int ibWo   = hPow_woCut->FindBin(x);
    double wo  = hPow_woCut->GetBinContent(ibWo);
    if (w <= 0 || wo <= 0) continue;
    double r   = w / wo;
    hPowRatio->SetBinContent(ib, r);
    hPowRatio->SetBinError  (ib, 0.);

    double relW  = hPow_wCutUnc->GetBinContent(ib);                 // relative POWHEG unc (wCut)
    double relWo = hPow_woCutUnc->GetBinContent(hPow_woCutUnc->FindBin(x)); // relative (woCut)
    double relR  = std::sqrt(relW*relW + relWo*relWo);
    double exLow  = x - hPow_wCut->GetBinLowEdge(ib);
    double exHigh = hPow_wCut->GetBinLowEdge(ib) + hPow_wCut->GetBinWidth(ib) - x;
    gPowRatioBand->SetPoint(kpr, x, r);
    gPowRatioBand->SetPointError(kpr, exLow, exHigh, r*relR, r*relR);
    kpr++;
  }

  // ---------------------------------------------------------------------------
  // 4.  STYLING
  // ---------------------------------------------------------------------------
  // Run3 without cut (blue family)
  hRun3_woCut->SetMarkerStyle(kFullCircle);
  hRun3_woCut->SetMarkerColor(kBlue+1);
  hRun3_woCut->SetLineColor(kBlue+1);
  hRun3_woCut->SetMarkerSize(1.1);

  // Run3 with pT-lead cut (red family)
  hRun3_wCut->SetMarkerStyle(kFullSquare);
  hRun3_wCut->SetMarkerColor(kRed+1);
  hRun3_wCut->SetLineColor(kRed+1);
  hRun3_wCut->SetMarkerSize(1.1);

  // POWHEG no-cut (blue line + light blue band)
  hPow_woCut->SetLineColor(kBlue+1);
  hPow_woCut->SetLineWidth(3);
  hPow_woCut->SetLineStyle(1);
  hPow_woCut->SetMarkerSize(0);

  TGraphAsymmErrors* gPowBand_woCut = new TGraphAsymmErrors();
  {
    int k = 0;
    for (int ib = 1; ib <= hPow_woCut->GetNbinsX(); ++ib) {
      double x = hPow_woCut->GetBinCenter(ib);
      double exLow  = x - hPow_woCut->GetBinLowEdge(ib);
      double exHigh = hPow_woCut->GetBinLowEdge(ib) + hPow_woCut->GetBinWidth(ib) - x;
      double c = hPow_woCut->GetBinContent(ib);
      double r = hPow_woCutUnc->GetBinContent(ib);
      if (c <= 0) continue;
      gPowBand_woCut->SetPoint(k, x, c);
      gPowBand_woCut->SetPointError(k, exLow, exHigh, c*r, c*r);
      k++;
    }
  }
  gPowBand_woCut->SetFillColorAlpha(kBlue+1, 0.25);
  gPowBand_woCut->SetLineColor(kBlue+1);
  gPowBand_woCut->SetFillStyle(1001);

  // POWHEG with-cut (red line + light red band)
  hPow_wCut->SetLineColor(kRed+1);
  hPow_wCut->SetLineWidth(3);
  hPow_wCut->SetLineStyle(2);
  hPow_wCut->SetMarkerSize(0);

  TGraphAsymmErrors* gPowBand_wCut = new TGraphAsymmErrors();
  {
    int k = 0;
    for (int ib = 1; ib <= hPow_wCut->GetNbinsX(); ++ib) {
      double x = hPow_wCut->GetBinCenter(ib);
      double exLow  = x - hPow_wCut->GetBinLowEdge(ib);
      double exHigh = hPow_wCut->GetBinLowEdge(ib) + hPow_wCut->GetBinWidth(ib) - x;
      double c = hPow_wCut->GetBinContent(ib);
      double r = hPow_wCutUnc->GetBinContent(ib);
      if (c <= 0) continue;
      gPowBand_wCut->SetPoint(k, x, c);
      gPowBand_wCut->SetPointError(k, exLow, exHigh, c*r, c*r);
      k++;
    }
  }
  gPowBand_wCut->SetFillColorAlpha(kRed+1, 0.25);
  gPowBand_wCut->SetLineColor(kRed+1);
  gPowBand_wCut->SetFillStyle(1001);

  // POWHEG ratio styling
  hPowRatio->SetLineColor(kGreen+2);
  hPowRatio->SetLineWidth(3);
  hPowRatio->SetLineStyle(1);
  hPowRatio->SetMarkerSize(0);
  gPowRatioBand->SetFillColorAlpha(kGreen+2, 0.30);
  gPowRatioBand->SetLineColor(kGreen+2);
  gPowRatioBand->SetFillStyle(1001);

  // Ratio
  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(1.0);

  // ---------------------------------------------------------------------------
  // 5.  CANVAS + PADS
  // ---------------------------------------------------------------------------
  TString canvasName = Form("cRun3_woCut_vs_wCut_R%.1f", arrayRadius[iRadius]);
  TCanvas* c = new TCanvas(canvasName, canvasName, 800, 900);
  c->SetFillStyle(0);

  TPad* padUp = new TPad("padUp", "padUp", 0., 0.30, 1., 1.);
  padUp->SetBottomMargin(0.015);
  padUp->SetTopMargin(0.08);
  padUp->SetLeftMargin(0.16);
  padUp->SetRightMargin(0.05);
  padUp->SetLogy();
  padUp->SetTickx(1); padUp->SetTicky(1);
  padUp->Draw();

  TPad* padDn = new TPad("padDn", "padDn", 0., 0.00, 1., 0.30);
  padDn->SetTopMargin(0.015);
  padDn->SetBottomMargin(0.38);
  padDn->SetLeftMargin(0.16);
  padDn->SetRightMargin(0.05);
  padDn->SetTickx(1); padDn->SetTicky(1);
  padDn->Draw();

  // ---------------------------------------------------------------------------
  // 6.  UPPER PAD — both data spectra + both POWHEG curves
  // ---------------------------------------------------------------------------
  padUp->cd();

  double xLo = commonEdges.front();
  double xHi = commonEdges.back();

  hRun3_woCut->GetXaxis()->SetRangeUser(xLo, xHi);
  hRun3_woCut->GetXaxis()->SetLabelSize(0.);
  hRun3_woCut->GetYaxis()->SetTitle("d^{2}#sigma_{jet}/d#it{p}_{T}d#it{#eta} mb (GeV/#it{c})^{-1}");
  hRun3_woCut->GetYaxis()->SetTitleSize(0.055);
  hRun3_woCut->GetYaxis()->SetTitleOffset(1.30);
  hRun3_woCut->GetYaxis()->SetLabelSize(0.050);
  hRun3_woCut->SetStats(0);

  double yMax = std::max(hRun3_woCut->GetMaximum(), hRun3_wCut->GetMaximum()) * 8.;
  double yMin = std::min(hRun3_woCut->GetMinimum(0.), hRun3_wCut->GetMinimum(0.)) * 0.1;
  hRun3_woCut->SetMinimum(yMin > 0. ? yMin : 1e-9);
  hRun3_woCut->SetMaximum(yMax);

  hRun3_woCut->Draw("E1 X0");
  gPowBand_woCut->Draw("2 SAME");
  gPowBand_wCut ->Draw("2 SAME");
  hPow_woCut->Draw("hist SAME");
  hPow_wCut ->Draw("hist SAME");
  hRun3_woCut->Draw("E1 X0 SAME");   // redraw markers on top of bands
  hRun3_wCut ->Draw("E1 X0 SAME");

  TLegend *leg = new TLegend(0.45, 0.62, 0.88, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.030);
  leg->AddEntry(hRun3_woCut,    "Data pp #sqrt{#it{s}} = 5.36 TeV",                                  "lep");
  leg->AddEntry(hPow_woCut,     "POWHEG+Pythia8, no #it{p}_{T}^{lead} cut",                          "lf");
  leg->AddEntry(gPowBand_woCut, "POWHEG unc. (no cut)",                                              "f");
  leg->AddEntry(hRun3_wCut,     "Data, #it{p}_{T}^{lead} > 5 GeV/#it{c}",                            "lep");
  leg->AddEntry(hPow_wCut,      "POWHEG+Pythia8, #it{p}_{T}^{lead} > 5 GeV/#it{c}",                  "lf");
  leg->AddEntry(gPowBand_wCut,  "POWHEG unc. (cut)",                                                 "f");
  leg->Draw();

  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextSize(0.040);
  ltx.DrawLatex(0.18, 0.90, "Anti-#it{k}_{T} charged jets,  |#it{#eta}_{jet}| < 0.7 #minus #it{R}");
  ltx.SetTextSize(0.036);
  ltx.DrawLatex(0.18, 0.85, "Stat. errors only");

  padUp->RedrawAxis();

  // ---------------------------------------------------------------------------
  // 7.  LOWER PAD — data ratio  wCut / woCut  +  POWHEG ratio
  // ---------------------------------------------------------------------------
  padDn->cd();

  TLine* lineUnity = new TLine(xLo, 1., xHi, 1.);
  lineUnity->SetLineStyle(2);
  lineUnity->SetLineColor(kGray+2);
  lineUnity->SetLineWidth(1);

  hRatio->GetXaxis()->SetRangeUser(xLo, xHi);
  hRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hRatio->GetYaxis()->SetTitle("(w cut) / (no cut)");

  const double scaleFactor = (1. - 0.30) / 0.30;
  hRatio->GetXaxis()->SetTitleSize  (0.055 * scaleFactor);
  hRatio->GetXaxis()->SetLabelSize  (0.048 * scaleFactor);
  hRatio->GetXaxis()->SetTitleOffset(0.85);
  hRatio->GetYaxis()->SetTitleSize  (0.048 * scaleFactor);
  hRatio->GetYaxis()->SetLabelSize  (0.044 * scaleFactor);
  hRatio->GetYaxis()->SetTitleOffset(0.48);
  hRatio->GetYaxis()->SetNdivisions(504);
  hRatio->GetXaxis()->SetTickLength (0.06);
  hRatio->GetYaxis()->SetTickLength (0.03);

  hRatio->SetMinimum(0.0);
  hRatio->SetMaximum(1.4);
  hRatio->SetStats(0);

  hRatio->Draw("E1 X0");
  gPowRatioBand->Draw("2 SAME");
  hPowRatio->Draw("hist SAME");
  hRatio->Draw("E1 X0 SAME");      // data ratio markers on top
  lineUnity->Draw("SAME");

  TLegend* legR = new TLegend(0.50, 0.70, 0.88, 0.95);
  legR->SetBorderSize(0);
  legR->SetFillStyle(0);
  legR->SetTextFont(42);
  legR->SetTextSize(0.048 * scaleFactor * 0.75);
  legR->AddEntry(hRatio,        "Data ratio (binomial stat.)",   "lep");
  legR->AddEntry(hPowRatio,     "POWHEG ratio",                  "l");
  legR->AddEntry(gPowRatioBand, "POWHEG unc. (quad. sum)",       "f");
  legR->Draw();

  padDn->RedrawAxis();

  // ---------------------------------------------------------------------------
  // // 8.  SAVE
  // // ---------------------------------------------------------------------------
  // TString pdfOut = Form("Run3_wCut_vs_woCut_R%.1f.pdf", arrayRadius[iRadius]);
  // c->SaveAs(pdfOut);
  // cout << "Saved: " << pdfOut << endl;
}




// ============================================================================
// PlotDataVsPOWHEG_R04.C
//
// Overlays:
//   - Data: H1D_XSection_Pt_Unfolded_Rebinned_w_MB  (d sigma / dpT d eta)
//   - POWHEG central: Central_Inclusive_R04
//   - POWHEG total uncertainty band: hTotalUnc_Inclusive_R04 (relative, stored as bin content)
//
// Layout:
//   Top pad    : spectra on log scale with POWHEG uncertainty band
//   Bottom pad : Data / POWHEG ratio
// ============================================================================


// void Draw_Data_Vs_POWHEG_R02_L3(){
//     gStyle->SetOptStat(0);
//     gStyle->SetPadTickX(1);
//     gStyle->SetPadTickY(1);

//     // -------------------------------------------------------------------------
//     // 1. Open files
//     // -------------------------------------------------------------------------
//     const char *dataFile   = "/Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/ROOT_file_XSection/output.root";
//     const char *powhegFile = "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_ptLead3/POWHEG_Uncertainties.root";
//     // relative stat : /Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/Tracking_eff/sys_TrackEff.root (has diffrent bining) then add a flat 2.5 % quadratically

//     TFile *fData   = TFile::Open(dataFile);
//     TFile *fPow    = TFile::Open(powhegFile);

//     if (!fData  || fData->IsZombie())  { std::cerr << "ERROR: cannot open data file\n";   return; }
//     if (!fPow   || fPow->IsZombie())   { std::cerr << "ERROR: cannot open POWHEG file\n"; return; }

//     // -------------------------------------------------------------------------
//     // 2. Load histograms
//     // -------------------------------------------------------------------------
//     TH1F *hData    = (TH1F *)fData->Get("H1D_XSection_Pt_Unfolded_Rebinned_w_MB");
//     TH1F *hPow     = (TH1F *)fPow->Get("Central_Inclusive_R02");
//     TH1F *hPowUnc  = (TH1F *)fPow->Get("hTotalUnc_Inclusive_R02");   // relative uncertainty

//     if (!hData)   { std::cerr << "ERROR: H1D_XSection_Pt_Unfolded_Rebinned_w_MB not found\n"; return; }
//     if (!hPow)    { std::cerr << "ERROR: Central_Inclusive_R02 not found\n";                   return; }
//     if (!hPowUnc) { std::cerr << "ERROR: hTotalUnc_Inclusive_R02 not found\n";                 return; }

//     // Clone to avoid modifying originals
//     hData   = (TH1F *)hData->Clone("hData_clone");
//     hPow    = (TH1F *)hPow->Clone("hPow_clone");
//     hPowUnc = (TH1F *)hPowUnc->Clone("hPowUnc_clone");

//     // -------------------------------------------------------------------------
//     // 3. Build POWHEG uncertainty band as TGraphAsymmErrors
//     //    hPowUnc stores relative uncertainty per bin (as computed by CombineQuadrature)
//     // -------------------------------------------------------------------------
//     int nBins = hPow->GetNbinsX();
//     TGraphAsymmErrors *gPowBand = new TGraphAsymmErrors(nBins);
//     for (int ib = 1; ib <= nBins; ib++)
//     {
//         double x      = hPow->GetBinCenter(ib);
//         double exLow  = x - hPow->GetBinLowEdge(ib);
//         double exHigh = hPow->GetBinLowEdge(ib) + hPow->GetBinWidth(ib) - x;
//         double c      = hPow->GetBinContent(ib);
//         double relUnc = hPowUnc->GetBinContent(ib);   // fractional (e.g. 0.15 = 15%)
//         gPowBand->SetPoint(ib-1, x, c);
//         gPowBand->SetPointError(ib-1, exLow, exHigh, c*relUnc, c*relUnc);
//     }
//     gPowBand->SetFillColorAlpha(kAzure+7, 0.35);
//     gPowBand->SetLineColor(kAzure+7);
//     gPowBand->SetFillStyle(1001);

//     // -------------------------------------------------------------------------
//     // 4. Build ratio  Data / POWHEG  (bin-by-bin, matching by pT center)
//     //    Also build POWHEG uncertainty band on ratio (symmetric around 1)
//     // -------------------------------------------------------------------------
//     // Use POWHEG binning as reference for the ratio
//     TH1F *hRatio = (TH1F *)hPow->Clone("hRatio");
//     hRatio->Reset();
//     hRatio->SetTitle("");

//     TGraphAsymmErrors *gRatioBand = new TGraphAsymmErrors(nBins);

//     for (int ibP = 1; ibP <= nBins; ibP++)
//     {
//         double ptCenter = hPow->GetBinCenter(ibP);
//         int ibD = hData->FindBin(ptCenter);

//         double powVal = hPow->GetBinContent(ibP);
//         double datVal = hData->GetBinContent(ibD);
//         double datErr = hData->GetBinError(ibD);

//         if (powVal > 0 && datVal > 0)
//         {
//             hRatio->SetBinContent(ibP, datVal / powVal);
//             hRatio->SetBinError(ibP,   datErr / powVal);
//         }

//         // POWHEG uncertainty band centered at 1
//         double x      = hPow->GetBinCenter(ibP);
//         double exLow  = x - hPow->GetBinLowEdge(ibP);
//         double exHigh = hPow->GetBinLowEdge(ibP) + hPow->GetBinWidth(ibP) - x;
//         double relUnc = hPowUnc->GetBinContent(ibP);
//         gRatioBand->SetPoint(ibP-1, x, 1.0);
//         gRatioBand->SetPointError(ibP-1, exLow, exHigh, relUnc, relUnc);
//     }
//     gRatioBand->SetFillColorAlpha(kAzure+7, 0.35);
//     gRatioBand->SetLineColor(kAzure+7);
//     gRatioBand->SetFillStyle(1001);

//     // -------------------------------------------------------------------------
//     // 5. Styling
//     // -------------------------------------------------------------------------
//     // Data
//     hData->SetMarkerStyle(20);
//     hData->SetMarkerSize(1.1);
//     hData->SetMarkerColor(kBlack);
//     hData->SetLineColor(kBlack);
//     hData->SetLineWidth(2);

//     // POWHEG central line
//     hPow->SetLineColor(kAzure+2);
//     hPow->SetLineWidth(3);
//     hPow->SetMarkerSize(0);

//     // Ratio
//     hRatio->SetMarkerStyle(20);
//     hRatio->SetMarkerSize(1.1);
//     hRatio->SetMarkerColor(kBlack);
//     hRatio->SetLineColor(kBlack);
//     hRatio->SetLineWidth(2);

//     // -------------------------------------------------------------------------
//     // 6. Canvas
//     // -------------------------------------------------------------------------
//     TCanvas *cv = new TCanvas("cDataVsPOWHEG_R02", "Data vs POWHEG R=0.2", 800, 900);
//     cv->cd();

//     // --- Top pad ---
//     TPad *pad1 = new TPad("pad1", "", 0.0, 0.32, 1.0, 1.0);
//     pad1->SetBottomMargin(0.015);
//     pad1->SetLeftMargin(0.14);
//     pad1->SetRightMargin(0.05);
//     pad1->SetLogy();
//     pad1->Draw();
//     pad1->cd();

//     // Determine common axis range
//     double xMin = hPow->GetXaxis()->GetXmin();
//     double xMax = hPow->GetXaxis()->GetXmax();

//     // Draw band first, then lines, then data
//     gPowBand->GetXaxis()->SetLimits(xMin, xMax);
//     gPowBand->Draw("A2");   // "2" = filled band, no axis from TGraph use frame below

//     // Use POWHEG histo to set frame
//     hPow->GetXaxis()->SetLabelSize(0);
//     hPow->GetXaxis()->SetTitleSize(0);
//     hPow->GetYaxis()->SetTitle("d#sigma/dp_{T}d#eta  (mb GeV^{-1}c)");
//     hPow->GetYaxis()->SetTitleSize(0.058);
//     hPow->GetYaxis()->SetTitleOffset(1.10);
//     hPow->GetYaxis()->SetLabelSize(0.052);
//     hPow->Draw("hist same");
//     gPowBand->Draw("2 same");
//     hPow->Draw("hist same");   // redraw line on top of band
//     hData->Draw("E1 same");

//   TLegend *leg = new TLegend(0.45, 0.49, 0.88, 0.71);

//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);
//     leg->SetTextSize(0.03);
//     leg->AddEntry(hData,    "ALICE, pp #sqrt{s}=5.36 TeV", "ep");
//     leg->AddEntry(hPow,     "POWHEG+Pythia8 dijet with CT18NNLO",                  "l");
//     leg->AddEntry(gPowBand, "POWHEG unc. (Scale, PDF, #alpha_{s} variations)",                    "f");
//     leg->Draw();

//     TLatex lat;
//     lat.SetNDC(); lat.SetTextFont(42);
//     lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.83, "Anti-k_{T} jets, R=0.2");
//     lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.75, "|#eta_{jet}| < 0.7,  p_{T,lead} > 3 GeV/c");

//     // Redraw axes on top
//     pad1->RedrawAxis();

//     // --- Bottom pad ---
//     cv->cd();
//     TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.32);
//     pad2->SetTopMargin(0.015);
//     pad2->SetBottomMargin(0.30);
//     pad2->SetLeftMargin(0.14);
//     pad2->SetRightMargin(0.05);
//     pad2->Draw();
//     pad2->cd();

//     hRatio->GetXaxis()->SetTitle("p_{T,jet}  (GeV/c)");
//     hRatio->GetXaxis()->SetTitleSize(0.115);
//     hRatio->GetXaxis()->SetTitleOffset(1.05);
//     hRatio->GetXaxis()->SetLabelSize(0.100);
//     hRatio->GetYaxis()->SetTitle("Data / POWHEG");
//     hRatio->GetYaxis()->SetTitleSize(0.100);
//     hRatio->GetYaxis()->SetTitleOffset(0.58);
//     hRatio->GetYaxis()->SetLabelSize(0.095);
//     hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
//     hRatio->GetYaxis()->SetNdivisions(505);
//     hRatio->SetTitle("");

//     hRatio->Draw("E1");
//     gRatioBand->Draw("2 same");
//     hRatio->Draw("E1 same");   // redraw points on top of band

//     // Unity line
//     TLine *line = new TLine(xMin, 1.0, xMax, 1.0);
//     line->SetLineColor(kAzure+2);
//     line->SetLineWidth(2);
//     line->SetLineStyle(2);
//     line->Draw();

//     pad2->RedrawAxis();

//     // -------------------------------------------------------------------------
//     // 7. Save
//     // -------------------------------------------------------------------------
//     cv->SaveAs("DataVsPOWHEG_R04_Lead3.pdf");
//     cv->SaveAs("DataVsPOWHEG_R04_Lead3.png");
//     std::cout << "\nSaved: DataVsPOWHEG_R04.pdf  and  DataVsPOWHEG_R04.png\n";
// }

void Draw_Data_Vs_POWHEG_R02_L3(){
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // -------------------------------------------------------------------------
    // 1. Open files
    // -------------------------------------------------------------------------
    // const char *dataFile    = "/Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/ROOT_file_XSection/output.root"; // Woeseok method with k 10
    // const char *dataFile    = "/Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/XSection_Nima_Method/Xsection.root"; // Nima methode with k = 10
    // const char *dataFile    = "/Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/ROOT_file_XSection/Xsection_k=13.root"; // wooseok method with k 13
    const char *dataFile    = "/Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/Xsection_option1.root"; // normalisation of 1st June Final
    const char *powhegFile  = "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_ptLead3/POWHEG_Uncertainties.root";
    const char *systFile    = "/Users/tabikh/Documents/PhdWork/MyWork/20260520_R02_Lead3_Xsection/Tracking_eff/sys_TrackEff.root";

    // Flat uncertainty to be added in quadrature (e.g. global normalization, etc.)
    const double kFlatRelSyst = 0.025;   // 2.5 %

    TFile *fData = TFile::Open(dataFile);
    TFile *fPow  = TFile::Open(powhegFile);
    TFile *fSyst = TFile::Open(systFile);

    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: cannot open data file\n";   return; }
    if (!fPow  || fPow->IsZombie())  { std::cerr << "ERROR: cannot open POWHEG file\n"; return; }
    if (!fSyst || fSyst->IsZombie()) { std::cerr << "ERROR: cannot open syst  file\n";  return; }

    // -------------------------------------------------------------------------
    // 2. Load histograms
    // -------------------------------------------------------------------------
    TH1F *hData    = (TH1F *)fData->Get("H1D_XSection_Pt_Unfolded_Rebinned_w_MB");
    TH1F *hPow     = (TH1F *)fPow->Get("Central_Inclusive_R02");
    TH1F *hPowUnc  = (TH1F *)fPow->Get("hTotalUnc_Inclusive_R02");

    TH1 *hRelSyst = (TH1 *)fSyst->Get("Rel_TrackEff_svd");

    if (!hData)    { std::cerr << "ERROR: data histo not found\n";                      return; }
    if (!hPow)     { std::cerr << "ERROR: Central_Inclusive_R02 not found\n";           return; }
    if (!hPowUnc)  { std::cerr << "ERROR: hTotalUnc_Inclusive_R02 not found\n";         return; }
    if (!hRelSyst) { std::cerr << "ERROR: Rel_TrackEff_svd not found in syst file\n";  return; }

    hData    = (TH1F *)hData->Clone("hData_clone");
    hPow     = (TH1F *)hPow->Clone("hPow_clone");
    hPowUnc  = (TH1F *)hPowUnc->Clone("hPowUnc_clone");
    hRelSyst = (TH1  *)hRelSyst->Clone("hRelSyst_clone");
    hRelSyst->SetDirectory(0);

    // Print syst binning so the user can verify
    std::cout << "\n--- Syst (tracking-eff) histogram binning ---\n";
    for (int ib = 1; ib <= hRelSyst->GetNbinsX(); ++ib) {
        printf("  bin %2d: [%6.2f, %6.2f] GeV/c   rel = %.4f\n",
               ib,
               hRelSyst->GetBinLowEdge(ib),
               hRelSyst->GetBinLowEdge(ib) + hRelSyst->GetBinWidth(ib),
               hRelSyst->GetBinContent(ib));
    }
    std::cout << "---------------------------------------------\n\n";

    // Helper: evaluate the syst histogram at any pT (handles different binning).
    // Returns the bin content of the syst histo whose range contains pT;
    // if pT falls outside the syst histo range, uses the closest valid bin.
    auto GetRelSystAt = [&](double pt) -> double {
        int ib = hRelSyst->FindFixBin(pt);
        if (ib < 1)                       ib = 1;
        if (ib > hRelSyst->GetNbinsX())   ib = hRelSyst->GetNbinsX();
        return hRelSyst->GetBinContent(ib);
    };

    // -------------------------------------------------------------------------
    // 3. POWHEG uncertainty band
    // -------------------------------------------------------------------------
    int nBins = hPow->GetNbinsX();
    TGraphAsymmErrors *gPowBand = new TGraphAsymmErrors(nBins);
    for (int ib = 1; ib <= nBins; ib++) {
        double x      = hPow->GetBinCenter(ib);
        double exLow  = x - hPow->GetBinLowEdge(ib);
        double exHigh = hPow->GetBinLowEdge(ib) + hPow->GetBinWidth(ib) - x;
        double c      = hPow->GetBinContent(ib);
        double r      = hPowUnc->GetBinContent(ib);
        gPowBand->SetPoint(ib-1, x, c);
        gPowBand->SetPointError(ib-1, exLow, exHigh, c*r, c*r);
    }
    gPowBand->SetFillColorAlpha(kAzure+7, 0.35);
    gPowBand->SetLineColor(kAzure+7);
    gPowBand->SetFillStyle(1001);

    // -------------------------------------------------------------------------
    // 4. Data SYSTEMATIC band on the spectrum
    //    total relative syst = sqrt( tracking-eff^2  +  flat^2 )
    // -------------------------------------------------------------------------
    int nBinsData = hData->GetNbinsX();
    TGraphAsymmErrors *gDataSyst = new TGraphAsymmErrors();
    int kd = 0;
    for (int ib = 1; ib <= nBinsData; ib++) {
        double x = hData->GetBinCenter(ib);
        double y = hData->GetBinContent(ib);
        if (y <= 0) continue;
        double exLow  = x - hData->GetBinLowEdge(ib);
        double exHigh = hData->GetBinLowEdge(ib) + hData->GetBinWidth(ib) - x;

        double relTrk  = GetRelSystAt(x);
        double relTot  = std::sqrt(relTrk*relTrk + kFlatRelSyst*kFlatRelSyst);

        gDataSyst->SetPoint(kd, x, y);
        gDataSyst->SetPointError(kd, exLow, exHigh, y*relTot, y*relTot);
        kd++;

        printf("pT %6.2f  rel_trk=%.4f  rel_flat=%.4f  rel_tot=%.4f\n",
               x, relTrk, kFlatRelSyst, relTot);
    }
    gDataSyst->SetFillColorAlpha(kGray+1, 0.45);
    gDataSyst->SetLineColor(kGray+2);
    gDataSyst->SetFillStyle(1001);
    gDataSyst->SetMarkerSize(0);

    // -------------------------------------------------------------------------
    // 5. Build ratio  Data / POWHEG  +  POWHEG band at 1  +  Data syst band on ratio
    // -------------------------------------------------------------------------
    TH1F *hRatio = (TH1F *)hPow->Clone("hRatio");
    hRatio->Reset(); hRatio->SetTitle("");

    TGraphAsymmErrors *gRatioBand     = new TGraphAsymmErrors(nBins);   // POWHEG band at 1
    TGraphAsymmErrors *gDataSystRatio = new TGraphAsymmErrors();        // data syst on ratio
    int kr = 0;

    for (int ibP = 1; ibP <= nBins; ibP++) {
        double ptCenter = hPow->GetBinCenter(ibP);
        int    ibD      = hData->FindBin(ptCenter);
        double powVal   = hPow->GetBinContent(ibP);
        double datVal   = hData->GetBinContent(ibD);
        double datErr   = hData->GetBinError(ibD);   // stat only

        double exLow  = ptCenter - hPow->GetBinLowEdge(ibP);
        double exHigh = hPow->GetBinLowEdge(ibP) + hPow->GetBinWidth(ibP) - ptCenter;

        if (powVal > 0 && datVal > 0) {
            double ratio = datVal / powVal;
            hRatio->SetBinContent(ibP, ratio);
            hRatio->SetBinError(ibP,   datErr / powVal);   // stat only

            // Data syst centered on the ratio point (kept separate from stat)
            double relTrk = GetRelSystAt(ptCenter);
            double relTot = std::sqrt(relTrk*relTrk + kFlatRelSyst*kFlatRelSyst);
            gDataSystRatio->SetPoint(kr, ptCenter, ratio);
            gDataSystRatio->SetPointError(kr, exLow, exHigh, ratio*relTot, ratio*relTot);
            kr++;
        }

        // POWHEG band centered at 1
        double relPow = hPowUnc->GetBinContent(ibP);
        gRatioBand->SetPoint(ibP-1, ptCenter, 1.0);
        gRatioBand->SetPointError(ibP-1, exLow, exHigh, relPow, relPow);
    }
    gRatioBand->SetFillColorAlpha(kAzure+7, 0.35);
    gRatioBand->SetLineColor(kAzure+7);
    gRatioBand->SetFillStyle(1001);

    gDataSystRatio->SetFillColorAlpha(kGray+1, 0.45);
    gDataSystRatio->SetLineColor(kGray+2);
    gDataSystRatio->SetFillStyle(1001);
    gDataSystRatio->SetMarkerSize(0);

    // -------------------------------------------------------------------------
    // 6. Styling
    // -------------------------------------------------------------------------
    hData->SetMarkerStyle(20); hData->SetMarkerSize(1.1);
    hData->SetMarkerColor(kBlack); hData->SetLineColor(kBlack); hData->SetLineWidth(2);

    hPow->SetLineColor(kAzure+2); hPow->SetLineWidth(3); hPow->SetMarkerSize(0);

    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.1);
    hRatio->SetMarkerColor(kBlack); hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2);

    // -------------------------------------------------------------------------
    // 7. Canvas
    // -------------------------------------------------------------------------
    TCanvas *cv = new TCanvas("cDataVsPOWHEG_R02", "Data vs POWHEG R=0.2", 800, 900);
    cv->cd();

    // --- Top pad ---
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.32, 1.0, 1.0);
    pad1->SetBottomMargin(0.015); pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);  pad1->SetLogy();
    pad1->Draw(); pad1->cd();

    double xMin = hPow->GetXaxis()->GetXmin();
    double xMax = hPow->GetXaxis()->GetXmax();

    gPowBand->GetXaxis()->SetLimits(xMin, xMax);
    gPowBand->Draw("A2");

    hPow->GetXaxis()->SetLabelSize(0); hPow->GetXaxis()->SetTitleSize(0);
    hPow->GetYaxis()->SetTitle("d#sigma/dp_{T}d#eta  (mb GeV^{-1}c)");
    hPow->GetYaxis()->SetTitleSize(0.058); hPow->GetYaxis()->SetTitleOffset(1.10);
    hPow->GetYaxis()->SetLabelSize(0.052);
    hPow->Draw("hist same");
    gPowBand->Draw("2 same");
    hPow->Draw("hist same");
    gDataSyst->Draw("2 same");          // data syst box
    hData->Draw("E1 same");              // data points + stat err

    TLegend *leg = new TLegend(0.45, 0.45, 0.88, 0.71);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
    leg->AddEntry(hData,     "ALICE, pp #sqrt{s}=5.36 TeV",                            "ep");
    leg->AddEntry(gDataSyst, "Data syst. unc. ",               "f");
    leg->AddEntry(hPow,      "POWHEG+Pythia8 dijet with CT18NNLO",                     "l");
    leg->AddEntry(gPowBand,  "POWHEG unc. (Scale, PDF, #alpha_{s} variations)",        "f");
    leg->Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
    lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.83, "Anti-k_{T} jets, R=0.2");
    lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.75, "|#eta_{jet}| < 0.7,  p_{T,lead} > 3 GeV/c");

    pad1->RedrawAxis();

    // --- Bottom pad ---
    cv->cd();
    TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.32);
    pad2->SetTopMargin(0.015); pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.14); pad2->SetRightMargin(0.05);
    pad2->Draw(); pad2->cd();

    hRatio->GetXaxis()->SetTitle("p_{T,jet}  (GeV/c)");
    hRatio->GetXaxis()->SetTitleSize(0.115); hRatio->GetXaxis()->SetTitleOffset(1.05);
    hRatio->GetXaxis()->SetLabelSize(0.100);
    hRatio->GetYaxis()->SetTitle("Data / POWHEG");
    hRatio->GetYaxis()->SetTitleSize(0.100); hRatio->GetYaxis()->SetTitleOffset(0.58);
    hRatio->GetYaxis()->SetLabelSize(0.095);
    hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->SetTitle("");

    hRatio->Draw("E1");
    gRatioBand->Draw("2 same");          // POWHEG band at 1
    gDataSystRatio->Draw("2 same");      // data syst band on ratio
    hRatio->Draw("E1 same");             // ratio points on top

    TLine *line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(kAzure+2); line->SetLineWidth(2); line->SetLineStyle(2);
    line->Draw();

    pad2->RedrawAxis();

    // // -------------------------------------------------------------------------
    // // 8. Save
    // // -------------------------------------------------------------------------
    // cv->SaveAs("DataVsPOWHEG_R02_Lead3.pdf");
    // cv->SaveAs("DataVsPOWHEG_R02_Lead3.png");
    // std::cout << "\nSaved: DataVsPOWHEG_R02_Lead3.pdf  and  .png\n";
}



void Draw_Data_Vs_POWHEG_R02_woLcut(){
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // -------------------------------------------------------------------------
    // 1. Open files
    // -------------------------------------------------------------------------
    const char *dataFile   = "/Users/tabikh/Documents/PhdWork/MyWork/20260530_R02_673480/Xsection/Xsection_woptLead.root";
    const char *powhegFile = "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_woPtLeadCut/POWHEG_Uncertainties.root";
    // extrapolated : 

    TFile *fData   = TFile::Open(dataFile);
    TFile *fPow    = TFile::Open(powhegFile);

    if (!fData  || fData->IsZombie())  { std::cerr << "ERROR: cannot open data file\n";   return; }
    if (!fPow   || fPow->IsZombie())   { std::cerr << "ERROR: cannot open POWHEG file\n"; return; }

    // -------------------------------------------------------------------------
    // 2. Load histograms
    // -------------------------------------------------------------------------
    TH1F *hData    = (TH1F *)fData->Get("H1D_XSection_Pt_Unfolded_Rebinned_w_MB");
    TH1F *hPow     = (TH1F *)fPow->Get("Central_Inclusive_R02");
    TH1F *hPowUnc  = (TH1F *)fPow->Get("hTotalUnc_Inclusive_R02");   // relative uncertainty

    if (!hData)   { std::cerr << "ERROR: H1D_XSection_Pt_Unfolded_Rebinned_w_MB not found\n"; return; }
    if (!hPow)    { std::cerr << "ERROR: Central_Inclusive_R02 not found\n";                   return; }
    if (!hPowUnc) { std::cerr << "ERROR: hTotalUnc_Inclusive_R02 not found\n";                 return; }

    // Clone to avoid modifying originals
    hData   = (TH1F *)hData->Clone("hData_clone");
    hPow    = (TH1F *)hPow->Clone("hPow_clone");
    hPowUnc = (TH1F *)hPowUnc->Clone("hPowUnc_clone");

    // -------------------------------------------------------------------------
    // 3. Build POWHEG uncertainty band as TGraphAsymmErrors
    //    hPowUnc stores relative uncertainty per bin (as computed by CombineQuadrature)
    // -------------------------------------------------------------------------
    int nBins = hPow->GetNbinsX();
    TGraphAsymmErrors *gPowBand = new TGraphAsymmErrors(nBins);
    for (int ib = 1; ib <= nBins; ib++)
    {
        double x      = hPow->GetBinCenter(ib);
        double exLow  = x - hPow->GetBinLowEdge(ib);
        double exHigh = hPow->GetBinLowEdge(ib) + hPow->GetBinWidth(ib) - x;
        double c      = hPow->GetBinContent(ib);
        double relUnc = hPowUnc->GetBinContent(ib);   // fractional (e.g. 0.15 = 15%)
        gPowBand->SetPoint(ib-1, x, c);
        gPowBand->SetPointError(ib-1, exLow, exHigh, c*relUnc, c*relUnc);
    }
    gPowBand->SetFillColorAlpha(kAzure+7, 0.35);
    gPowBand->SetLineColor(kAzure+7);
    gPowBand->SetFillStyle(1001);

    // -------------------------------------------------------------------------
    // 4. Build ratio  Data / POWHEG  (bin-by-bin, matching by pT center)
    //    Also build POWHEG uncertainty band on ratio (symmetric around 1)
    // -------------------------------------------------------------------------
    // Use POWHEG binning as reference for the ratio
    TH1F *hRatio = (TH1F *)hPow->Clone("hRatio");
    hRatio->Reset();
    hRatio->SetTitle("");

    TGraphAsymmErrors *gRatioBand = new TGraphAsymmErrors(nBins);

    for (int ibP = 1; ibP <= nBins; ibP++)
    {
        double ptCenter = hPow->GetBinCenter(ibP);
        int ibD = hData->FindBin(ptCenter);

        double powVal = hPow->GetBinContent(ibP);
        double datVal = hData->GetBinContent(ibD);
        double datErr = hData->GetBinError(ibD);

        if (powVal > 0 && datVal > 0)
        {
            hRatio->SetBinContent(ibP, datVal / powVal);
            hRatio->SetBinError(ibP,   datErr / powVal);
        }

        // POWHEG uncertainty band centered at 1
        double x      = hPow->GetBinCenter(ibP);
        double exLow  = x - hPow->GetBinLowEdge(ibP);
        double exHigh = hPow->GetBinLowEdge(ibP) + hPow->GetBinWidth(ibP) - x;
        double relUnc = hPowUnc->GetBinContent(ibP);
        gRatioBand->SetPoint(ibP-1, x, 1.0);
        gRatioBand->SetPointError(ibP-1, exLow, exHigh, relUnc, relUnc);
    }
    gRatioBand->SetFillColorAlpha(kAzure+7, 0.35);
    gRatioBand->SetLineColor(kAzure+7);
    gRatioBand->SetFillStyle(1001);

    // -------------------------------------------------------------------------
    // 5. Styling
    // -------------------------------------------------------------------------
    // Data
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.1);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(2);

    // POWHEG central line
    hPow->SetLineColor(kAzure+2);
    hPow->SetLineWidth(3);
    hPow->SetMarkerSize(0);

    // Ratio
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.1);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->SetLineWidth(2);

    // -------------------------------------------------------------------------
    // 6. Canvas
    // -------------------------------------------------------------------------
    TCanvas *cv = new TCanvas("cDataVsPOWHEG_R02", "Data vs POWHEG R=0.2", 800, 900);
    cv->cd();

    // --- Top pad ---
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.32, 1.0, 1.0);
    pad1->SetBottomMargin(0.015);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    // Determine common axis range
    double xMin = hPow->GetXaxis()->GetXmin();
    double xMax = hPow->GetXaxis()->GetXmax();

    // Draw band first, then lines, then data
    gPowBand->GetXaxis()->SetLimits(xMin, xMax);
    gPowBand->Draw("A2");   // "2" = filled band, no axis from TGraph use frame below

    // Use POWHEG histo to set frame
    hPow->GetXaxis()->SetLabelSize(0);
    hPow->GetXaxis()->SetTitleSize(0);
    hPow->GetYaxis()->SetTitle("d#sigma/dp_{T}d#eta  (mb GeV^{-1}c)");
    hPow->GetYaxis()->SetTitleSize(0.058);
    hPow->GetYaxis()->SetTitleOffset(1.10);
    hPow->GetYaxis()->SetLabelSize(0.052);
    hPow->Draw("hist same");
    gPowBand->Draw("2 same");
    hPow->Draw("hist same");   // redraw line on top of band
    hData->Draw("E1 same");

  TLegend *leg = new TLegend(0.45, 0.49, 0.88, 0.71);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(hData,    "ALICE, pp #sqrt{s}=5.36 TeV", "ep");
    leg->AddEntry(hPow,     "POWHEG+Pythia8 dijet with CT18NNLO",                  "l");
    leg->AddEntry(gPowBand, "POWHEG unc. (Scale, PDF, #alpha_{s} variations)",                    "f");
    leg->Draw();

    TLatex lat;
    lat.SetNDC(); lat.SetTextFont(42);
    lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.83, "Anti-k_{T} jets, R=0.2");
    lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.75, "|#eta_{jet}| < 0.7");

    // Redraw axes on top
    pad1->RedrawAxis();

    // --- Bottom pad ---
    cv->cd();
    TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.32);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    hRatio->GetXaxis()->SetTitle("p_{T,jet}  (GeV/c)");
    hRatio->GetXaxis()->SetTitleSize(0.115);
    hRatio->GetXaxis()->SetTitleOffset(1.05);
    hRatio->GetXaxis()->SetLabelSize(0.100);
    hRatio->GetYaxis()->SetTitle("Data / POWHEG");
    hRatio->GetYaxis()->SetTitleSize(0.100);
    hRatio->GetYaxis()->SetTitleOffset(0.58);
    hRatio->GetYaxis()->SetLabelSize(0.095);
    hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->SetTitle("");

    hRatio->Draw("E1");
    gRatioBand->Draw("2 same");
    hRatio->Draw("E1 same");   // redraw points on top of band

    // Unity line
    TLine *line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(kAzure+2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();

    pad2->RedrawAxis();

    // // -------------------------------------------------------------------------
    // // 7. Save
    // // -------------------------------------------------------------------------
    // cv->SaveAs("DataVsPOWHEG_R04_Lead3.pdf");
    // cv->SaveAs("DataVsPOWHEG_R04_Lead3.png");
    // std::cout << "\nSaved: DataVsPOWHEG_R02.pdf  and  DataVsPOWHEG_R04.png\n";
}

void Draw_Data_Vs_POWHEG_R02_PtLcut(){
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // -------------------------------------------------------------------------
    // 1. Open files
    // -------------------------------------------------------------------------
    const char *dataFile   = "/Users/tabikh/Documents/PhdWork/MyWork/20260530_R02_673480/Xsection/Xsection_woptLead.root";
    const char *powhegFile = "/Users/tabikh/Documents/PhdWork/MyWork/Datasets/POWHEG_Run3_woPtLeadCut/POWHEG_Uncertainties.root";

    TFile *fData   = TFile::Open(dataFile);
    TFile *fPow    = TFile::Open(powhegFile);

    if (!fData  || fData->IsZombie())  { std::cerr << "ERROR: cannot open data file\n";   return; }
    if (!fPow   || fPow->IsZombie())   { std::cerr << "ERROR: cannot open POWHEG file\n"; return; }

    // -------------------------------------------------------------------------
    // 2. Load histograms
    // -------------------------------------------------------------------------
    TH1F *hData    = (TH1F *)fData->Get("H1D_XSection_Pt_Unfolded_Rebinned_w_MB");
    TH1F *hPow     = (TH1F *)fPow->Get("Central_Inclusive_R02");
    TH1F *hPowUnc  = (TH1F *)fPow->Get("hTotalUnc_Inclusive_R02");   // relative uncertainty

    if (!hData)   { std::cerr << "ERROR: H1D_XSection_Pt_Unfolded_Rebinned_w_MB not found\n"; return; }
    if (!hPow)    { std::cerr << "ERROR: Central_Inclusive_R02 not found\n";                   return; }
    if (!hPowUnc) { std::cerr << "ERROR: hTotalUnc_Inclusive_R02 not found\n";                 return; }

    // Clone to avoid modifying originals
    hData   = (TH1F *)hData->Clone("hData_clone");
    hPow    = (TH1F *)hPow->Clone("hPow_clone");
    hPowUnc = (TH1F *)hPowUnc->Clone("hPowUnc_clone");

    // -------------------------------------------------------------------------
    // 3. Build POWHEG uncertainty band as TGraphAsymmErrors
    //    hPowUnc stores relative uncertainty per bin (as computed by CombineQuadrature)
    // -------------------------------------------------------------------------
    int nBins = hPow->GetNbinsX();
    TGraphAsymmErrors *gPowBand = new TGraphAsymmErrors(nBins);
    for (int ib = 1; ib <= nBins; ib++)
    {
        double x      = hPow->GetBinCenter(ib);
        double exLow  = x - hPow->GetBinLowEdge(ib);
        double exHigh = hPow->GetBinLowEdge(ib) + hPow->GetBinWidth(ib) - x;
        double c      = hPow->GetBinContent(ib);
        double relUnc = hPowUnc->GetBinContent(ib);   // fractional (e.g. 0.15 = 15%)
        gPowBand->SetPoint(ib-1, x, c);
        gPowBand->SetPointError(ib-1, exLow, exHigh, c*relUnc, c*relUnc);
    }
    gPowBand->SetFillColorAlpha(kAzure+7, 0.35);
    gPowBand->SetLineColor(kAzure+7);
    gPowBand->SetFillStyle(1001);

    // -------------------------------------------------------------------------
    // 4. Build ratio  Data / POWHEG  (bin-by-bin, matching by pT center)
    //    Also build POWHEG uncertainty band on ratio (symmetric around 1)
    // -------------------------------------------------------------------------
    // Use POWHEG binning as reference for the ratio
    TH1F *hRatio = (TH1F *)hPow->Clone("hRatio");
    hRatio->Reset();
    hRatio->SetTitle("");

    TGraphAsymmErrors *gRatioBand = new TGraphAsymmErrors(nBins);

    for (int ibP = 1; ibP <= nBins; ibP++)
    {
        double ptCenter = hPow->GetBinCenter(ibP);
        int ibD = hData->FindBin(ptCenter);

        double powVal = hPow->GetBinContent(ibP);
        double datVal = hData->GetBinContent(ibD);
        double datErr = hData->GetBinError(ibD);

        if (powVal > 0 && datVal > 0)
        {
            hRatio->SetBinContent(ibP, datVal / powVal);
            hRatio->SetBinError(ibP,   datErr / powVal);
        }

        // POWHEG uncertainty band centered at 1
        double x      = hPow->GetBinCenter(ibP);
        double exLow  = x - hPow->GetBinLowEdge(ibP);
        double exHigh = hPow->GetBinLowEdge(ibP) + hPow->GetBinWidth(ibP) - x;
        double relUnc = hPowUnc->GetBinContent(ibP);
        gRatioBand->SetPoint(ibP-1, x, 1.0);
        gRatioBand->SetPointError(ibP-1, exLow, exHigh, relUnc, relUnc);
    }
    gRatioBand->SetFillColorAlpha(kAzure+7, 0.35);
    gRatioBand->SetLineColor(kAzure+7);
    gRatioBand->SetFillStyle(1001);

    // -------------------------------------------------------------------------
    // 5. Styling
    // -------------------------------------------------------------------------
    // Data
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.1);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(2);

    // POWHEG central line
    hPow->SetLineColor(kAzure+2);
    hPow->SetLineWidth(3);
    hPow->SetMarkerSize(0);

    // Ratio
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.1);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->SetLineWidth(2);

    // -------------------------------------------------------------------------
    // 6. Canvas
    // -------------------------------------------------------------------------
    TCanvas *cv = new TCanvas("cDataVsPOWHEG_R02", "Data vs POWHEG R=0.2", 800, 900);
    cv->cd();

    // --- Top pad ---
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.32, 1.0, 1.0);
    pad1->SetBottomMargin(0.015);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    // Determine common axis range
    double xMin = hPow->GetXaxis()->GetXmin();
    double xMax = hPow->GetXaxis()->GetXmax();

    // Draw band first, then lines, then data
    gPowBand->GetXaxis()->SetLimits(xMin, xMax);
    gPowBand->Draw("A2");   // "2" = filled band, no axis from TGraph use frame below

    // Use POWHEG histo to set frame
    hPow->GetXaxis()->SetLabelSize(0);
    hPow->GetXaxis()->SetTitleSize(0);
    hPow->GetYaxis()->SetTitle("d#sigma/dp_{T}d#eta  (mb GeV^{-1}c)");
    hPow->GetYaxis()->SetTitleSize(0.058);
    hPow->GetYaxis()->SetTitleOffset(1.10);
    hPow->GetYaxis()->SetLabelSize(0.052);
    hPow->Draw("hist same");
    gPowBand->Draw("2 same");
    hPow->Draw("hist same");   // redraw line on top of band
    hData->Draw("E1 same");

  TLegend *leg = new TLegend(0.45, 0.49, 0.88, 0.71);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(hData,    "ALICE, pp #sqrt{s}=5.36 TeV", "ep");
    leg->AddEntry(hPow,     "POWHEG+Pythia8 dijet with CT18NNLO",                  "l");
    leg->AddEntry(gPowBand, "POWHEG unc. (Scale, PDF, #alpha_{s} variations)",                    "f");
    leg->Draw();

    TLatex lat;
    lat.SetNDC(); lat.SetTextFont(42);
    lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.83, "Anti-k_{T} jets, R=0.2");
    lat.SetTextSize(0.038); lat.DrawLatex(0.45, 0.75, "|#eta_{jet}| < 0.7");

    // Redraw axes on top
    pad1->RedrawAxis();

    // --- Bottom pad ---
    cv->cd();
    TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.32);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    hRatio->GetXaxis()->SetTitle("p_{T,jet}  (GeV/c)");
    hRatio->GetXaxis()->SetTitleSize(0.115);
    hRatio->GetXaxis()->SetTitleOffset(1.05);
    hRatio->GetXaxis()->SetLabelSize(0.100);
    hRatio->GetYaxis()->SetTitle("Data / POWHEG");
    hRatio->GetYaxis()->SetTitleSize(0.100);
    hRatio->GetYaxis()->SetTitleOffset(0.58);
    hRatio->GetYaxis()->SetLabelSize(0.095);
    hRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->SetTitle("");

    hRatio->Draw("E1");
    gRatioBand->Draw("2 same");
    hRatio->Draw("E1 same");   // redraw points on top of band

    // Unity line
    TLine *line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineColor(kAzure+2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();

    pad2->RedrawAxis();

    // // -------------------------------------------------------------------------
    // // 7. Save
    // // -------------------------------------------------------------------------
    // cv->SaveAs("DataVsPOWHEG_R04_Lead3.pdf");
    // cv->SaveAs("DataVsPOWHEG_R04_Lead3.png");
    // std::cout << "\nSaved: DataVsPOWHEG_R02.pdf  and  DataVsPOWHEG_R04.png\n";
}





// rename refoldedUnfolded as closure test?
// and try and spend 15 min to clean hist names for the spectrum analysis



// WARNING FOR EFFICIENCIES I SHOULD REREAD THIS BELOW!!
// hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

#endif

