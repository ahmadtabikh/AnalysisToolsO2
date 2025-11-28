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
#include <RooUnfold.h>
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

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();



void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_systematics() {
  // Load necessary libraries
  LoadLibs_Systematics();
  // Set the default style
  SetStyle_Systematics();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // gathers the analysis options in a single char[]

  int iDataset = 0;
  int iRadius = 1;


  char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s,%s", mergingPrior, unfoldingPrior);


  const int nUnfoldingMethods = 2;
  char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes"}; // default is the first one in this list
  int unfoldParameterInputList[2] = {7, 5};

  Draw_Systematics_UnfoldMethod(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);
}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs_Systematics() {
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

void SetStyle_Systematics(Bool_t graypalette) {
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


void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {

  TH1D* hTempSystematicUncertainty = new TH1D("hTempSystematicUncertainty", "hTempSystematicUncertainty", nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]);
  TH1D* hTempSystematicUncertainty_PreBarlow = new TH1D("hTempSystematicUncertainty_PreBarlow", "hTempSystematicUncertainty_PreBarlow", nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]);
  hTempSystematicUncertainty->Sumw2();
  hTempSystematicUncertainty_PreBarlow->Sumw2();
  hTempSystematicUncertainty->Reset("M");
  hTempSystematicUncertainty_PreBarlow->Reset("M");
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // return histogram that has the systematics in its contents
  TH1D* H1D_jetPt_unfolded[nUnfoldingMethods];
  TH1D* H1D_jetPt_unfolded_differences[nUnfoldingMethods-1];


  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsBeforeUnfolding) {
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

  char optionsAnalysis_withUnfoldingMethod[100] = "";
  for(int iMethod = 0; iMethod < nUnfoldingMethods; iMethod++){
    snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDataset, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);

    if (iMethod != 0) {
      H1D_jetPt_unfolded_differences[iMethod-1] = (TH1D*)H1D_jetPt_unfolded[iMethod]->Clone("H1D_jetPt_unfolded_differences"+partialUniqueSpecifier);
      H1D_jetPt_unfolded_differences[iMethod-1]->Add(H1D_jetPt_unfolded[0],-1);
    }
    cout << "do I want the absolute value of the difference?" << endl;
  }

  // cout << "Do I apply Barlow condition even though not param variation ? check paper again" << endl; YES, the subset of data thing is only shown for first demonstration, but barlow says it holds true even if that's not the cast

  /////////////////
  // Barlow test //
  /////////////////

  TH1D* H1D_jetPt_unfolded_REF = H1D_jetPt_unfolded[0];
  double SystUncertainty;
  int id_SignalExtractionType_maxDeviation;
  double hSigmaBarlow[nBinPtJetsGen[iRadius]];
  for(int iBinPt = 1; iBinPt <= nBinPtJetsGen[iRadius]; iBinPt++){
    SystUncertainty = 0;
    for(int iMethod = 1; iMethod < nUnfoldingMethods; iMethod++){ // get maximum difference among the nUnfoldingMethods-1 ones, hold value with SystUncertainty, and the id of the method wîth id_SignalExtractionType_maxDeviation
      if (abs(H1D_jetPt_unfolded_differences[iMethod-1]->GetBinContent(iBinPt)) > SystUncertainty) {
        SystUncertainty = abs(H1D_jetPt_unfolded_differences[iMethod-1]->GetBinContent(iBinPt));
        id_SignalExtractionType_maxDeviation = iMethod;
      }
    }

    // Barlow condition for systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = H1D_jetPt_unfolded_REF->GetBinError(iBinPt);
    Double_t StatUncertainty_MaxDeviationCase = H1D_jetPt_unfolded[id_SignalExtractionType_maxDeviation]->GetBinError(iBinPt);
 
    int PtArrayIterator = iBinPt - 1;
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample : sqrt( | sigma1_{unfBayes}^2 - sigma_{unfSVD}^2 | )
 
    hTempSystematicUncertainty_PreBarlow->SetBinContent(iBinPt,SystUncertainty);
    hTempSystematicUncertainty_PreBarlow->SetBinError(iBinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hTempSystematicUncertainty->SetBinContent(iBinPt,SystUncertainty);
      hTempSystematicUncertainty->SetBinError(iBinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hTempSystematicUncertainty->SetBinContent(iBinPt,0.);
    }
  }

  hSystematicUncertainty = (TH1D*)hTempSystematicUncertainty->Clone("hSystematicUncertainty_UnfoldMethod"+partialUniqueSpecifier);
  hSystematicUncertainty_PreBarlow = (TH1D*)hTempSystematicUncertainty_PreBarlow->Clone("hSystematicUncertainty_PreBarlow_UnfoldMethod"+partialUniqueSpecifier);

}



void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {

  TH1D* hSystematicUncertainty;
  TH1D* hSystematicUncertainty_PreBarlow;
  Get_systematics_UnfoldMethod(hSystematicUncertainty, hSystematicUncertainty_PreBarlow, iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, options);

  TH1D* H1D_jetPt_unfolded;
  char optionsAnalysis_withUnfoldingMethod[100] = "";
  snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethod);

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsBeforeUnfolding) {
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
  
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInputList[0], optionsAnalysis_withUnfoldingMethod); //SVD unfolding as reference
  hSystematicUncertainty->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_PreBarlow->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"]_"+unfoldingMethodList[0]+"_kUnfold="+Form("%i", unfoldParameterInputList[0]);

  TString* pdfName = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier);
  TString* pdfName_PreBarlow = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier+"_PreBarlow");

  TString textContext("");

  Draw_TH1_Histogram(hSystematicUncertainty, textContext, pdfName, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSystematicUncertainty_PreBarlow, textContext, pdfName_PreBarlow, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}