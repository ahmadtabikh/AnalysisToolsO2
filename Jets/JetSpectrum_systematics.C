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
#include <vector>
#include <TRandom3.h>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();



void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);
void Draw_Systematics_trackefficiency(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_TrackEfficiency_systematics();


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
  int unfoldParameterInputList[2] = {7, 4}; // first and third for nominal, second and fourth for systematics variation
  Draw_Systematics_UnfoldMethod(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);


  // char optionsAnalysis[100] = "";
  // snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", mergingPrior, unfoldingPrior, unfoldingMethod);
  // int unfoldParameterInputMin = 6;
  // int unfoldParameterInputMax = 8;
  // int unfoldParameterInputStep = 1;
  // Draw_Systematics_parameterVariation(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);


  // Draw_TrackEfficiency_systematics(); // use a root file with all unfolded spectra with methods and track efficiency variations


  // const int nUnfoldingMethodsTrackeff = 4;
  // char* unfoldingMethodListTrackeff[nUnfoldingMethodsTrackeff] = {"Svd", "Svd", "Bayes", "Bayes"}; // default is the first one in this list
  // int unfoldParameterInputListTrackeff[4] = {7, 4, 4, 2}; // first and third for nominal, second and fourth for systematics variation
  // Draw_Systematics_trackefficiency(iDataset, iRadius, unfoldingMethodListTrackeff, unfoldParameterInputListTrackeff, nUnfoldingMethodsTrackeff, optionsAnalysis);


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
  TString partialUniqueSpecifier = Datasets[iDataset]+DatasetsMC[0]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // return histogram that has the systematics in its contents
  TH1D* H1D_jetPt_unfolded[nUnfoldingMethods];
  TH1D* H1D_jetPt_unfolded_differences[nUnfoldingMethods-1];
  TH1D* H1D_jetPt_unfolded_ratio[nUnfoldingMethods-1];


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

    // write in root file the unfolded spectra
    // TFile* outputFile_UnfoldedSpectra = new TFile("TrackEff_Systematics.root","UPDATE");
    // std::string histName = "H1D_jetPt_unfolded_" + 
    //                    std::string(unfoldingMethodList[iMethod]) + 
    //                    "_param" + 
    //                    std::to_string(unfoldParameterInputList[iMethod]) + 
    //                    "_" + 
    //                    std::string(partialUniqueSpecifier.Data());

    // H1D_jetPt_unfolded[iMethod]->Write(histName.c_str());
    // outputFile_UnfoldedSpectra->Close();  

    // --- DIFFERENCE histogram (method_i - method_0)
    if (iMethod != 0) {
      H1D_jetPt_unfolded_differences[iMethod-1] = (TH1D*)H1D_jetPt_unfolded[iMethod]->Clone("H1D_jetPt_unfolded_differences"+partialUniqueSpecifier);
      H1D_jetPt_unfolded_differences[iMethod-1]->Add(H1D_jetPt_unfolded[0],-1);
    
    cout << "do I want the absolute value of the difference?" << endl;

    // // --- RATIO histogram (method_i / method_0)
    //   H1D_jetPt_unfolded_ratio[iMethod-1] = (TH1D*)H1D_jetPt_unfolded[iMethod]->Clone(("H1D_jetPt_unfolded_ratio" + partialUniqueSpecifier).c_str());
    //   H1D_jetPt_unfolded_ratio[iMethod-1]->Divide(H1D_jetPt_unfolded[0]);
    }  
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
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase + StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample : sqrt( | sigma1_{unfBayes}^2 - sigma_{unfSVD}^2 | )
 
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
  hSystematicUncertainty->Scale(100.0);
  hSystematicUncertainty_PreBarlow->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_PreBarlow->Scale(100.0);


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"]_"+unfoldingMethodList[0]+"_kUnfold="+Form("%i", unfoldParameterInputList[0]);

  TString* pdfName = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier);
  TString* pdfName_PreBarlow = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier+"_PreBarlow");

  // TString textContext("sys. unfolding method");
  // TString textContext =
  // "#splitline{sys. unfolding method}"
  //            "{k_{svd} = 7, k_{bayes} = 4}";
  TString textContext = Form(
    "#splitline{sys. unfolding method}"
    "{k_{svd} = %d, k_{bayes} = %d}",
    unfoldParameterInputList[0],
    unfoldParameterInputList[1]
  );

  TString* texRelativeErrPercent = new TString ("relative error (%)");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 25}}};
  Draw_TH1_Histogram(hSystematicUncertainty, textContext, pdfName, texPtJetRecX, texRelativeErrPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSystematicUncertainty_PreBarlow, textContext, pdfName_PreBarlow, texPtJetRecX, texRelativeErrPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {
  cout << "########### Drawing systematics from parameter variation ###############" << endl;
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);

  TH1D* H1D_jetPt_unfolded[nUnfoldIteration];

  TString partialUniqueSpecifier;

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  int unfoldParameterInput = 0;

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

  if (measuredInput == nullptr) {
    cout << "Error: measuredInput histogram is null!" << endl;
    return;
  }
  else {
    cout << "measuredInput histogram successfully retrieved." << endl;
  }

  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    cout << "  entering the for loop "  << endl;
    unfoldParameterInput = unfoldIterationMax - iUnfoldIteration * step; 

    cout << "((((((((((((()))))))))))))" << endl;
    cout << "Iteration "<< iUnfoldIteration << endl;
    cout << "((((((((((((()))))))))))))" << endl;
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iUnfoldIteration], measuredInput, iDataset, iRadius, unfoldParameterInput, options);
  }

  int iNominal = nUnfoldIteration/2; // choose nominal (e.g. central iteration)
  TH1D* hNom = H1D_jetPt_unfolded[iNominal];

  TH1D* hSys_envelope = (TH1D*)hNom->Clone("hSys_envelope");
  hSys_envelope->Reset(); // will store absolute systematic (positive)

  int nBins = hNom->GetNbinsX();
  for (int ib = 1; ib <= nBins; ++ib) {
      double valNom = hNom->GetBinContent(ib);
      double maxAbs = 0.0;
      for (int i = 0; i < nUnfoldIteration; ++i) {
          double val = H1D_jetPt_unfolded[i]->GetBinContent(ib);
          double d = fabs(val - valNom);
          if (d > maxAbs) maxAbs = d;
      }
      hSys_envelope->SetBinContent(ib, maxAbs);
  }

  TH1D* hSys_rel = (TH1D*)hSys_envelope->Clone("hSys_relative");
  hSys_rel->Reset();

  for (int ib = 1; ib <= nBins; ++ib) {
      double absSys = hSys_envelope->GetBinContent(ib);
      double valNom = hNom->GetBinContent(ib);

      double rel = 0.0;
      if (valNom > 0) rel = absSys / valNom;

      hSys_rel->SetBinContent(ib, rel*100.0); // in percent
  }
  TString* pdfName_envolope = new TString("Systematics_deviation_ParameterVariation_"+partialUniqueSpecifier);
  TString* pdfName_relUnc = new TString("Systematics_RelativeUncertainty_ParameterVariation_"+partialUniqueSpecifier);
  TString textContext("SVD unfolding");
  TString* sigma = new TString("#sigma interation variation");
  TString* relativeErrors = new TString("realtive errors (%) ");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 1.6}}};
  // Draw_TH1_Histogram(hSys_envelope, textContext, pdfName_envolope, texPtJetRecX, sigma, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSys_rel, textContext, pdfName_relUnc, texPtJetRecX, relativeErrors, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  
}


void Draw_TrackEfficiency_systematics()
{
    // Open the ROOT file
    TFile *f = TFile::Open("TrackEff_Systematics.root", "READ");
    if (!f || f->IsZombie()) {
        Error("Draw_TrackEfficiency_systematics", "Cannot open file TrackEff_Systematics.root");
        return;
    }

    // Get histograms
    TH1D *H1D_Nominal_SVD = (TH1D*)f->Get("H1D_jetPt_unfolded_Svd_param7_LHC24_ppref_pass1_train380686LHC25b4b5_train533385_R=0.4"); //Svd 7 with nominal track efficiency
    TH1D *H1D_Nominal_Bayes = (TH1D*)f->Get("H1D_jetPt_unfolded_Bayes_param4_LHC24_ppref_pass1_train380686LHC25b4b5_train533385_R=0.4"); //Bayes 4 with nominal track efficiency
    TH1D *H1D_Reduced_SVD = (TH1D*)f->Get("H1D_jetPt_unfolded_Svd_param4_LHC24_ppref_pass1_train380686LHC25b4b5_train533385ReducedTrackEff_R=0.4");  //Svd 4 with reduced track efficiency
    TH1D *H1D_Reduced_Bayes = (TH1D*)f->Get("H1D_jetPt_unfolded_Bayes_param2_LHC24_ppref_pass1_train380686LHC25b4b5_train533385ReducedTrackEff_R=0.4"); //Bayes 2 with reduced track efficiency

    if (!H1D_Nominal_SVD || !H1D_Nominal_Bayes || !H1D_Reduced_SVD || !H1D_Reduced_Bayes) {
        Error("Draw_TrackEfficiency_systematics", "One or more histograms not found");
        f->Close();
        return;
    }

    // TString* pdfName_svdNominal = new TString("Nominal_SVD7");
    // TString* pdfName_bayesNominal = new TString("Nominal_Bayes4");
    // TString textContext("");
    // Draw_TH1_Histogram(H1D_Nominal_SVD, textContext, pdfName_svdNominal, texPtJetRecX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    // Draw_TH1_Histogram(H1D_Nominal_Bayes, textContext, pdfName_bayesNominal, texPtJetRecX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");

    // TString* pdfName_svdReduced = new TString("Reduced_SVD4");
    // TString* pdfName_bayesReduced = new TString("Reduced_Bayes2");
    // Draw_TH1_Histogram(H1D_Reduced_SVD, textContext, pdfName_svdReduced, texPtJetRecX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    // Draw_TH1_Histogram(H1D_Reduced_Bayes, textContext, pdfName_bayesReduced, texPtJetRecX, texJetPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");



    // Create histograms for absolute differences and relative uncertainties
    TH1D *H1D_Delta_SVD = (TH1D*)H1D_Nominal_SVD->Clone("hDiff_NominalReduced_SVD");
    H1D_Delta_SVD->Reset();
    TH1D *H1D_RelativeUncertainty_SVD = (TH1D*)H1D_Nominal_SVD->Clone("H1D_RelativeUncertainty_SVD");
    H1D_RelativeUncertainty_SVD->Reset();

    TH1D *H1D_Delta_Bayes = (TH1D*)H1D_Nominal_Bayes->Clone("H1D_Delta_Bayes");
    H1D_Delta_Bayes->Reset();
    TH1D *H1D_RelativeUncertainty_Bayes = (TH1D*)H1D_Nominal_Bayes->Clone("H1D_RelativeUncertainty_Bayes");
    H1D_RelativeUncertainty_Bayes->Reset();

    // Compute absolute differences bin by bin and compute relative uncertainties
    for (int i = 1; i <= H1D_Nominal_SVD->GetNbinsX(); ++i) {
        H1D_Delta_SVD->SetBinContent(i, abs(H1D_Nominal_SVD->GetBinContent(i) - H1D_Reduced_SVD->GetBinContent(i)));
        H1D_Delta_SVD->SetBinError(i, 0.0);

        if (H1D_Nominal_SVD->GetBinContent(i) != 0) {
            double relUnc = (H1D_Delta_SVD->GetBinContent(i) / H1D_Nominal_SVD->GetBinContent(i)) * 100.0; // in percent
            H1D_RelativeUncertainty_SVD->SetBinContent(i, relUnc);
            H1D_RelativeUncertainty_SVD->SetBinError(i, 0.0);
        } else {
            H1D_RelativeUncertainty_SVD->SetBinContent(i, 0.0);
            H1D_RelativeUncertainty_SVD->SetBinError(i, 0.0);
        }
    }

    for (int i = 1; i <= H1D_Nominal_Bayes->GetNbinsX(); ++i) {
        H1D_Delta_Bayes->SetBinContent(i, abs(H1D_Nominal_Bayes->GetBinContent(i) - H1D_Reduced_Bayes->GetBinContent(i)));
        H1D_Delta_Bayes->SetBinError(i, 0.0);

        if (H1D_Nominal_Bayes->GetBinContent(i) != 0) {
            double relUnc = (H1D_Delta_Bayes->GetBinContent(i) / H1D_Nominal_Bayes->GetBinContent(i)) * 100.0; // in percent
            H1D_RelativeUncertainty_Bayes->SetBinContent(i, relUnc);
            H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
        } else {
            H1D_RelativeUncertainty_Bayes->SetBinContent(i, 0.0);
            H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
        }
    }
    TString textContextSVD("with 3% reduced track efficiency (SVD)");
    TString textContextBayes("with 3% reduced track efficiency (Bayes)");

    TString* pdfName_relUncSvd_logx = new TString("Relative_uncert_Svd_Data_train380686_Nominal_train533385_param7_ReducedTrackEff3perCent_train564527_param4_R=0.4_logx");
    TString* pdfName_relUncBayes_logx = new TString("Relative_uncert_Bayes_Data_train380686_Nominal_train533385_param4_ReducedTrackEff3perCent_train564527_param2_R=0.4_logx");
    TString* pdfName_relUncSvd = new TString("Relative_uncert_Svd_Data_train380686_Nominal_train533385_param7_ReducedTrackEff3perCent_train564527_param4_R=0.4");
    TString* pdfName_relUncBayes = new TString("Relative_uncert_Bayes_Data_train380686_Nominal_train533385_param4_ReducedTrackEff3perCent_train564527_param2_R=0.4");
    std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 80}}};

    Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd_logx, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
    Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes_logx, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
    Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
    Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");

    ///
    // TString* pdfName_Svd = new TString(Form( "Systematics_Trackeff_SVD_R%d_Dataset%d_NominalK=%d_VarK=%d", iRadius, iDataset, unfoldParameterInputList[0], unfoldParameterInputList[1] ));
    // TString* pdfName_Bayes = new TString(Form( "Systematics_Trackeff_Bayes_R%d_Dataset%d_NominalK=%d_VarK=%d", iRadius, iDataset, unfoldParameterInputList[2], unfoldParameterInputList[3]));

    // TString textContext("");
    // TString* AbsoluteDifference = new TString ("Absolute difference");


    // Draw_TH1_Histogram(H1D_jetPt_unfolded_absDiff_01, textContext, pdfName_Svd, texPtJetRecX, AbsoluteDifference, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
    // Draw_TH1_Histogram(H1D_jetPt_unfolded_absDiff_23, textContext, pdfName_Bayes, texPtJetRecX, AbsoluteDifference, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");

}



void Draw_Systematics_trackefficiency(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options){

  TH1D* hTempSystematicUncertainty = new TH1D("hTempSystematicUncertainty", "hTempSystematicUncertainty", nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]);
  hTempSystematicUncertainty->Sumw2();
  hTempSystematicUncertainty->Reset("M");
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // return histogram that has the systematics in its contents
  TH1D* H1D_jetPt_unfolded[nUnfoldingMethods];
  TH1D* H1D_jetPt_unfolded_differences[nUnfoldingMethods-2];
  TH1D* H1D_jetPt_unfolded_ratio[nUnfoldingMethods-1];


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
    if (file_O2Analysis_MCfileForMatrix) {
        file_O2Analysis_MCfileForMatrix->Close();
        delete file_O2Analysis_MCfileForMatrix;
        file_O2Analysis_MCfileForMatrix = nullptr;
    }

    if (iMethod == 0 || iMethod == 2) {
        file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4b5_train588753/AnalysisResults.root");
    } else {
        file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4ab6_train564527_trackLoss3perCent/AnalysisResults.root");
    }
    snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDataset, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);
  }
  

  TH1D* H1D_Nominal_SVD = (TH1D*)H1D_jetPt_unfolded[0]->Clone("H1D_Nominal_SVD");
  TH1D* H1D_Reduced_SVD = (TH1D*)H1D_jetPt_unfolded[1]->Clone("H1D_Reduced_SVD"); 

  TH1D* H1D_Nominal_Bayes = (TH1D*)H1D_jetPt_unfolded[2]->Clone("H1D_Nominal_Bayes"); 
  TH1D* H1D_Reduced_Bayes = (TH1D*)H1D_jetPt_unfolded[3]->Clone("H1D_Reduced_Bayes"); 

  TH1D *H1D_Delta_SVD = (TH1D*)H1D_Nominal_SVD->Clone("hDiff_NominalReduced_SVD");
  H1D_Delta_SVD->Reset();
  TH1D *H1D_RelativeUncertainty_SVD = (TH1D*)H1D_Nominal_SVD->Clone("H1D_RelativeUncertainty_SVD");
  H1D_RelativeUncertainty_SVD->Reset();

  TH1D *H1D_Delta_Bayes = (TH1D*)H1D_Nominal_Bayes->Clone("H1D_Delta_Bayes");
  H1D_Delta_Bayes->Reset();
  TH1D *H1D_RelativeUncertainty_Bayes = (TH1D*)H1D_Nominal_Bayes->Clone("H1D_RelativeUncertainty_Bayes");
  H1D_RelativeUncertainty_Bayes->Reset();

  // Compute absolute differences bin by bin and compute relative uncertainties
  for (int i = 1; i <= H1D_Nominal_SVD->GetNbinsX(); ++i) {
      H1D_Delta_SVD->SetBinContent(i, abs(H1D_Nominal_SVD->GetBinContent(i) - H1D_Reduced_SVD->GetBinContent(i)));
      H1D_Delta_SVD->SetBinError(i, 0.0);

      if (H1D_Nominal_SVD->GetBinContent(i) != 0) {
          double relUnc = (H1D_Delta_SVD->GetBinContent(i) / H1D_Nominal_SVD->GetBinContent(i)) * 100.0; // in percent
          H1D_RelativeUncertainty_SVD->SetBinContent(i, relUnc);
          H1D_RelativeUncertainty_SVD->SetBinError(i, 0.0);
      } else {
          H1D_RelativeUncertainty_SVD->SetBinContent(i, 0.0);
          H1D_RelativeUncertainty_SVD->SetBinError(i, 0.0);
      }
  }

  for (int i = 1; i <= H1D_Nominal_Bayes->GetNbinsX(); ++i) {
      H1D_Delta_Bayes->SetBinContent(i, abs(H1D_Nominal_Bayes->GetBinContent(i) - H1D_Reduced_Bayes->GetBinContent(i)));
      H1D_Delta_Bayes->SetBinError(i, 0.0);

      if (H1D_Nominal_Bayes->GetBinContent(i) != 0) {
          double relUnc = (H1D_Delta_Bayes->GetBinContent(i) / H1D_Nominal_Bayes->GetBinContent(i)) * 100.0; // in percent
          H1D_RelativeUncertainty_Bayes->SetBinContent(i, relUnc);
          H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
      } else {
          H1D_RelativeUncertainty_Bayes->SetBinContent(i, 0.0);
          H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
      }
  }
  TString textContextSVD("with 3% reduced track efficiency (SVD)");
  TString textContextBayes("with 3% reduced track efficiency (Bayes)");

  TString* pdfName_relUncSvd_logx = new TString("Relative_uncert_Svd_Data_train380686_Nominal_train533385_param7_ReducedTrackEff3perCent_train564527_param4_R=0.4_logx");
  TString* pdfName_relUncBayes_logx = new TString("Relative_uncert_Bayes_Data_train380686_Nominal_train533385_param4_ReducedTrackEff3perCent_train564527_param2_R=0.4_logx");
  TString* pdfName_relUncSvd = new TString("Relative_uncert_Svd_Data_train380686_Nominal_train533385_param7_ReducedTrackEff3perCent_train564527_param4_R=0.4");
  TString* pdfName_relUncBayes = new TString("Relative_uncert_Bayes_Data_train380686_Nominal_train533385_param4_ReducedTrackEff3perCent_train564527_param2_R=0.4");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 80}}};

  Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd_logx, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes_logx, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");


  // TString* pdfName_Svd = new TString(Form( "Systematics_Trackeff_SVD_R%d_Dataset%d_NominalK=%d_VarK=%d",
  //       iRadius, iDataset,
  //       unfoldParameterInputList[0],
  //       unfoldParameterInputList[1]
  //   )
  // );


  // TString* pdfName_Bayes = new TString(Form( "Systematics_Trackeff_Bayes_R%d_Dataset%d_NominalK=%d_VarK=%d",
  //     iRadius, iDataset,
  //     unfoldParameterInputList[2],
  //     unfoldParameterInputList[3]
  // ));

  // TString textContext("");
  // TString* AbsoluteDifference = new TString ("Absolute difference");


  // Draw_TH1_Histogram(H1D_jetPt_unfolded_absDiff_01, textContext, pdfName_Svd, texPtJetRecX, AbsoluteDifference, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  // Draw_TH1_Histogram(H1D_jetPt_unfolded_absDiff_23, textContext, pdfName_Bayes, texPtJetRecX, AbsoluteDifference, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");




}



// Assume Get_Pt_spectrum_unfolded(TH1D* out, TH1D* measured, ...)
// is the function you already call and that produces unfolded histograms.
// Not complete yet, i should poperly call the function with the right parameters.
/*
void Draw_SvdBayesRatioWithCovariance(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {
    TH1D* hA_nominal,          // unfolded method 0 (reference)
    TH1D* hB_nominal,          // unfolded method 1
    TH1D* measuredInput,       // measured histogram (same input used by both)
    int nToys = 2000,
    unsigned int randomSeed = 0

    TRandom3 rnd(randomSeed?randomSeed:0);
    snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
      char optionsAnalysis_withUnfoldingMethod[100] = "";

    int nBins = hA_nominal->GetNbinsX();
    // Containers: values[method][bin][toyIndex]
    std::vector<std::vector<double>> valsA(nBins, std::vector<double>(nToys));
    std::vector<std::vector<double>> valsB(nBins, std::vector<double>(nToys));

    // Create temporary histograms for unfolding outputs
    for (int it = 0; it < nToys; ++it) {
        // --- make fluctuated measured histogram (Poisson per bin)
        TH1D* measuredToy = (TH1D*)measuredInput->Clone(Form("measuredToy_%d", it));
        for (int ib = 1; ib <= measuredToy->GetNbinsX(); ++ib) {
            double mu = measuredInput->GetBinContent(ib);
            double fluctuated = rnd.PoissonD(mu);
            measuredToy->SetBinContent(ib, fluctuated);
            // Also copy errors if you want (Poisson: sqrt)
            measuredToy->SetBinError(ib, sqrt(fluctuated));
        }

        // --- unfold with method 0 and method 1
        // Create output histos cloned from nominal shapes (to keep binning)
        TH1D* outA = (TH1D*)hA_nominal->Clone(Form("outA_toy_%d", it));
        TH1D* outB = (TH1D*)hB_nominal->Clone(Form("outB_toy_%d", it));

        // Call your unfolding wrapper for each method.
        // You must adapt the call signature to your function.
        // Example placeholders:
        // Get_Pt_spectrum_unfolded(outA, measuredToy, iDataset, iRadius, unfoldParam0, options0);
        // Get_Pt_spectrum_unfolded(outB, measuredToy, iDataset, iRadius, unfoldParam1, options1);

        // *** ADAPT THESE LINES TO MATCH YOUR FUNCTION CALLS ***
        Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDataset, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);
        // ********************************************************

        // Store bin values
        for (int ib = 0; ib < nBins; ++ib) {
            valsA[ib][it] = outA->GetBinContent(ib+1);
            valsB[ib][it] = outB->GetBinContent(ib+1);
        }

        // cleanup
        delete measuredToy;
        delete outA;
        delete outB;
    }

    // Compute means, variances, covariances
    std::vector<double> meanA(nBins,0.0), meanB(nBins,0.0);
    std::vector<double> varA(nBins,0.0), varB(nBins,0.0), covAB(nBins,0.0);

    for (int ib = 0; ib < nBins; ++ib) {
        // mean
        for (int it = 0; it < nToys; ++it) {
            meanA[ib] += valsA[ib][it];
            meanB[ib] += valsB[ib][it];
        }
        meanA[ib] /= nToys;
        meanB[ib] /= nToys;

        // variances and covariance (unbiased sample variance: divide by N-1)
        for (int it = 0; it < nToys; ++it) {
            double dA = valsA[ib][it] - meanA[ib];
            double dB = valsB[ib][it] - meanB[ib];
            varA[ib] += dA*dA;
            varB[ib] += dB*dB;
            covAB[ib] += dA*dB;
        }
        if (nToys > 1) {
            varA[ib] /= (nToys - 1);
            varB[ib] /= (nToys - 1);
            covAB[ib] /= (nToys - 1);
        }
    }

    // Build ratio histogram with propagated error using covariance
    TH1D* hRatio = (TH1D*)hB_nominal->Clone("hRatio_method1_over_method0");
    hRatio->Reset(); // empty contents/errs
    for (int ib = 0; ib < nBins; ++ib) {
        double A = hB_nominal->GetBinContent(ib+1); // note: using same notation as formula (A/B) adapt if needed
        double B = hA_nominal->GetBinContent(ib+1);

        // guard against zero denominator
        if (B == 0) {
            hRatio->SetBinContent(ib+1, 0);
            hRatio->SetBinError(ib+1, 0);
            continue;
        }
        double R = A / B;

        double varR = (1.0/(B*B)) * varA[ib]
                    + (A*A/(B*B*B*B)) * varB[ib]
                    - 2.0 * A / (B*B*B) * covAB[ib];

        if (varR < 0) varR = 0; // numerical safeguard
        double errR = sqrt(varR);

        hRatio->SetBinContent(ib+1, R);
        hRatio->SetBinError(ib+1, errR);
    }

    // Now you have hRatio with errors that include covariance.
    hRatio->Draw("E1");
}
*/


  // for(int iMethod = 0; iMethod < nUnfoldingMethods; iMethod++){
  //   if (iMethod == 0 || iMethod == 2) {
  //     TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4b5_train533385/AnalysisResults.root");
  //     // oldformatNEventsGenInO2Analysis = true;
  //     snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
  //     Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDataset, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);
  //   } else if (iMethod == 1 || iMethod == 3) {
  //     // TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4ab6_train564527_trackLoss3perCent/AnalysisResults.root");
  //       TFile* file_O2Analysis_MCfileForMatrix = new TFile("../Datasets/LHC25b4b5_train533385/AnalysisResults.root");
  //     // oldformatNEventsGenInO2Analysis = false;

  //     snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
  //     Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDataset, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);  
  //   }
  // }