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
void Draw_Syatematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);


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
  // char optionsAnalysis[100] = "";
  // snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", mergingPrior, unfoldingPrior, unfoldingMethod);
  // int unfoldParameterInputMin = 6;
  // int unfoldParameterInputMax = 8;
  // int unfoldParameterInputStep = 1;
  // Draw_Syatematics_parameterVariation(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);

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
  hSystematicUncertainty->Scale(100.0);
  hSystematicUncertainty_PreBarlow->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_PreBarlow->Scale(100.0);


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"]_"+unfoldingMethodList[0]+"_kUnfold="+Form("%i", unfoldParameterInputList[0]);

  TString* pdfName = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier);
  TString* pdfName_PreBarlow = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier+"_PreBarlow");

  TString textContext("");

  Draw_TH1_Histogram(hSystematicUncertainty, textContext, pdfName, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSystematicUncertainty_PreBarlow, textContext, pdfName_PreBarlow, texPtJetRecX, texSystematicsPercent, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_Syatematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {
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
  // Draw_TH1_Histogram(hSys_envelope, textContext, pdfName_envolope, texPtJetRecX, sigma, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSys_rel, textContext, pdfName_relUnc, texPtJetRecX, relativeErrors, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  
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