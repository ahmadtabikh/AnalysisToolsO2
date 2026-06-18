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
#include "../Utilities/Fits.C"
#include "../Utilities/Fits.h" 

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();


TH1D* ComputeSecondarySystematic(TH1* h, TF1* fitFunc) ;
void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);
void Draw_Systematics_SecondaryContamination(int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Draw_Systematics_TrackEff(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Binning_systematic();
double Get_zVertex_reconstruction_efficiency(int iDataset, int iRadius);
TH1D* Get_TVX_Eff(int iDataset, int iRadius);
TH1D* Get_UE_Systematic_Band();
void Draw_UE_systematic_band();
TH1D* Get_Total_Systematic(bool UE_trk_method );
void RewriteHistogramRange(TH1D*& h, int binStart, int binEnd);


/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_systematics() {
  // Load necessary libraries
  LoadLibs_Systematics();
  // Set the default style
  SetStyle_Systematics();
  cout<<"--- I am here 0----"<<endl; // very bizzar, if i remove it the function track efficiency crashes !!!

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // gathers the analysis options in a single char[]

  int iDataset = 0;
  int iRadius = 0;

  // //// ######################################################### Unf Method Systematics #####################################################
  // char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  // snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);


  // const int nUnfoldingMethods = 2;
  // char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes"}; // default is the first one in this list
  // int unfoldParameterInputList[2] = {13, 9};

  // Draw_Systematics_UnfoldMethod(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);
  // //// #############################################################################################################################################

  // ######################################################### Parameter variation Systematics #####################################################
  char optionsAnalysis[100] = "";
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", unfoldingPrior, unfoldingMethod);
  int unfoldParameterInputMin = 13;
  int unfoldParameterInputMax = 14;
  int unfoldParameterInputStep = 1;
  Draw_Systematics_parameterVariation(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);
  // #############################################################################################################################################



  // ////######################################################### Secondary tracks Systematics #####################################################
  // char optionsAnalysis[100] = "";
  // snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", unfoldingPrior, unfoldingMethod);
  // int unfoldParameterInput = 12;
  // Draw_Systematics_SecondaryContamination(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);
  // ////############################################################################################################################################# 2.2% low pt then 2.5% constant after 30

  // //######################################################### Track efficiency Systematics #####################################################
  // char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  // snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);
  // const int nUnfoldingMethods = 4;
  // char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes", "Svd", "Bayes"}; // first two to be with nominal efficiency, last two with efficiency varied 
  // int unfoldParameterInputList[4] = {10, 9, 10, 9}; // first two to be with nominal efficiency, last two with efficiency varied
  // cout<<"--- I am here ----"<<endl;
  // Draw_Systematics_TrackEff(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);
  // ////#############################################################################################################################################

  // Draw_Binning_systematic();
  // Draw_UE_systematic_band();
  // TH1D* hRun3_sys = Get_Total_Systematic(true);  // RELATIVE
 
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
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
 
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
  
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInputList[0], optionsAnalysis_withUnfoldingMethod);
  hSystematicUncertainty->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  hSystematicUncertainty->Scale(100.0);
  hSystematicUncertainty_PreBarlow->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  TFile* outFile = new TFile("sys_Method.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  hSystematicUncertainty_PreBarlow->Write("Rel_Sys_Method");
  outFile->Close();
  hSystematicUncertainty_PreBarlow->Scale(100.0);


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"]_"+unfoldingMethodList[0]+"_kUnfold="+Form("%i", unfoldParameterInputList[0]);

  TString* pdfName = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier);
  TString* pdfName_PreBarlow = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier+"_PreBarlow");

  // TString textContext("");
  TString textContext = Form(
    "#splitline{sys. unfolding method}"
    "{k_{svd} = %d, k_{bayes} = %d}",
    unfoldParameterInputList[0],
    unfoldParameterInputList[1]
  );

  std::array<std::array<float, 2>, 2> drawnWindow = {{{10, 100}, {0, 10}}};

  TString* texSystematicsPercent = new TString ("relative error (%)");

  Draw_TH1_Histogram(hSystematicUncertainty, textContext, pdfName, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSystematicUncertainty_PreBarlow, textContext, pdfName_PreBarlow, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_Systematics_SecondaryContamination(int iDataset, int iRadius, int unfoldParameterInput, std::string options){
  cout << "########### Drawing systematics from secondary contamination variation ###############" << endl;
  TH1D* H1D_jetPt_unfolded;

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
  
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options); 

  // Define your bins and fit range
  int nBinsX = H1D_jetPt_unfolded->GetNbinsX();
  double* binsX = new double[nBinsX+1];
  for(int i=0; i<=nBinsX; i++) binsX[i] = H1D_jetPt_unfolded->GetBinLowEdge(i+1);

  // double xRangeFit[2] = {5.0, 120.0}; // Fit range in GeV
  double xRangeFit[2] = {10.0, 140.0}; // Fit range in GeV

  // Step 1: Rebin histogram using double Tsallis fit
  std::tuple<TH1D*, TGraphErrors*, TF1*> result = RebinWithDoubleTsallisFit(H1D_jetPt_unfolded, nBinsX, binsX, xRangeFit);

  // Step 2: Extract outputs
  TH1D* hJetPtRebinned = std::get<0>(result);
  TGraphErrors* fitGraph = std::get<1>(result);
  TF1* fitFunctionDrawn = std::get<2>(result);

  TH1D* H1D_FitUnf_ratio = (TH1D*) hJetPtRebinned->Clone("hRatio");
  H1D_FitUnf_ratio->SetTitle("Ratio: Fit / Unfolded; p_{T}^{jet}; Ratio");
  H1D_FitUnf_ratio->Divide(H1D_jetPt_unfolded);

  TString* pdfName_fitUnfRatio = new TString("ratio fit histogram to unfolded spectrum");
  TString textContext("");
  TString* yLabel = new TString("Fit / Unfolded");
  Draw_TH1_Histogram(H1D_FitUnf_ratio, textContext, pdfName_fitUnfRatio, texPtJetRec, yLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "ratioLine");


  // Step 3: Draw original histogram and rebinned fit
  TCanvas* c1 = new TCanvas("c1", "Double Tsallis Fit", 800, 700);
  c1->SetLogy(); // VERY important for spectra
  c1->SetTicks(1,1);
  c1->SetLeftMargin(0.13);
  c1->SetBottomMargin(0.12);
  c1->Modified();
  c1->Update();

  // ── Original histogram ─────────────────────────────
  H1D_jetPt_unfolded->SetMarkerStyle(20);
  H1D_jetPt_unfolded->SetMarkerSize(1.0);
  H1D_jetPt_unfolded->SetMarkerColor(kBlack);
  H1D_jetPt_unfolded->SetLineColor(kBlack);
  H1D_jetPt_unfolded->GetXaxis()->SetRangeUser(10, 140); // example

  H1D_jetPt_unfolded->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  H1D_jetPt_unfolded->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T} d#eta");

  H1D_jetPt_unfolded->GetXaxis()->SetTitleSize(0.045);
  H1D_jetPt_unfolded->GetYaxis()->SetTitleSize(0.045);

  H1D_jetPt_unfolded->Draw("E");

  // ── Fit curve ──────────────────────────────────────
  fitGraph->SetLineColor(kRed+1);
  fitGraph->SetLineWidth(1);
  fitGraph->SetMarkerSize(0);

  // draw band first if needed
  // fitGraph->Draw("3 SAME");

  fitGraph->Draw("L SAME");

  // ── Rebinned histogram ─────────────────────────────
  hJetPtRebinned->SetMarkerStyle(24);
  hJetPtRebinned->SetMarkerSize(1.1);
  hJetPtRebinned->SetMarkerColor(kBlue+1);
  hJetPtRebinned->SetLineColor(kBlue+1);

  hJetPtRebinned->Draw("E SAME");

  // c1->BuildLegend();
  // ================= Legend =================
  TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);   // transparent
  leg->SetTextSize(0.035);

  leg->AddEntry(H1D_jetPt_unfolded, "Unfolded data", "lep");
  leg->AddEntry(fitGraph, "Double Tsallis fit", "l");
  leg->AddEntry(hJetPtRebinned, "Fit sampling (rebinned)", "lep");

  leg->Draw();
  c1->Update();

  TH1D* hSecSys = ComputeSecondarySystematic(H1D_jetPt_unfolded, fitFunctionDrawn);
  
  for (int i = 1; i <= hSecSys->GetNbinsX(); i++) {
    double val = hSecSys->GetBinContent(i);
    hSecSys->SetBinContent(i, val);
  }

  cout << "######################### FILE CREATED IN PRINCIPLE##################" << endl; 
  TFile* outFile_sys = new TFile("sys_Secondary.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  TString histoName = Form("Rel_Sys_secondary");
  hSecSys->Write(histoName);
  outFile_sys->Close();
  cout << "######################### HISTO SAVED IN PRINCIPLE ##################" << endl; 
  
  TCanvas* cSys = new TCanvas("cSys", "Secondary systematic", 800, 600);

  cSys->SetTicks(1,1);

  hSecSys->SetTitle("Secondary track systematic; p_{T} (GeV/c); Relative uncertainty");
  hSecSys->GetYaxis()->SetTitle("Systematic uncertainty (%)");

  hSecSys->SetLineColor(kRed+1);
  hSecSys->SetLineWidth(2);
  hSecSys->Scale(100);

  hSecSys->Draw("HIST");


  
}

TH1D* ComputeSecondarySystematic(TH1* h, TF1* fitFunc) {

    int nBins = h->GetNbinsX();

    TH1D* hSys = (TH1D*)h->Clone("hSecondarySystematic");
    hSys->Reset(); // we only store the systematic

    for (int i = 1; i <= nBins; i++) {

        double x = h->GetBinCenter(i);

        // central value
        double f0 = fitFunc->Eval(x);

        if (f0 <= 0) {
            hSys->SetBinContent(i, 0);
            continue;
        }

        // shifted evaluations (NO refit!)
        double f_up   = fitFunc->Eval(1.005 * x);
        double f_down = fitFunc->Eval(0.995 * x);

        // relative variations
        double delta_up   = (f_up   - f0) / f0;
        double delta_down = (f_down - f0) / f0;

        // take max variation
        double delta_sec = std::max(std::abs(delta_up), std::abs(delta_down));

        hSys->SetBinContent(i, delta_sec);
        hSys->SetBinError(i, 0); // systematic only
    }

    return hSys;
}

void Draw_Systematics_TrackEff(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {
  cout << "########### Drawing systematics from track efficiency variation ###############" << endl;
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // return histogram that has the systematics in its contents
  TH1D* H1D_jetPt_unfolded[nUnfoldingMethods];

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

  char optionsAnalysis_withUnfoldingMethod[100] = "";
  for(int iMethod = 0; iMethod < nUnfoldingMethods; iMethod++){
    snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
    int iDatasetMC = (iMethod == 0 || iMethod == 1) ? 0 : 1; // first two to be with nominal efficiency dataset 0, last two with efficiency varied dataset 1
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDatasetMC, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);
  }

  // Get histograms
  TH1D* H1D_Nominal_SVD = (TH1D*) H1D_jetPt_unfolded[0]->Clone("H1D_Nominal_SVD");
  TH1D* H1D_Nominal_Bayes = (TH1D*) H1D_jetPt_unfolded[1]->Clone("H1D_Nominal_Bayes");
  TH1D* H1D_Reduced_SVD = (TH1D*) H1D_jetPt_unfolded[2]->Clone("H1D_Reduced_SVD");
  TH1D* H1D_Reduced_Bayes = (TH1D*) H1D_jetPt_unfolded[3]->Clone("H1D_Reduced_Bayes");

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
          double relUnc = (H1D_Delta_SVD->GetBinContent(i) / H1D_Nominal_SVD->GetBinContent(i)) ; 
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
          double relUnc = (H1D_Delta_Bayes->GetBinContent(i) / H1D_Nominal_Bayes->GetBinContent(i));
          H1D_RelativeUncertainty_Bayes->SetBinContent(i, relUnc);
          H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
      } else {
          H1D_RelativeUncertainty_Bayes->SetBinContent(i, 0.0);
          H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
      }
  }
  TString textContextSVD("with reduced track efficiency (SVD)");
  TString textContextBayes("with reduced track efficiency (Bayes)");

  TString* pdfName_relUncSvd_logx = new TString("Relative_uncert_Svd_logx");
  TString* pdfName_relUncBayes_logx = new TString("Relative_uncert_Bayes_logx");
  TString* pdfName_relUncSvd = new TString("Relative_uncert_Svd");
  TString* pdfName_relUncBayes = new TString("Relative_uncert_Bayes");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{10, 100}, {0, 0.13}}};
  TString* texSystematicsPercent = new TString ("relative error (%)");

  Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd_logx, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes_logx, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  
  TFile* outFile = new TFile("sys_TrackEff.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  H1D_RelativeUncertainty_SVD->Write("Rel_TrackEff_svd");
  H1D_RelativeUncertainty_Bayes->Write("Rel_TrackEff_bayes");
  outFile->Close();

}


TH1D* Get_Binning_systematic(){
  TFile* Nominal = TFile::Open("../Unfolding_SVD_k12/output.root", "READ");
  TFile* Low_var = TFile::Open("lower_edge/output.root", "READ");
  TFile* End_var = TFile::Open("higher_edge/output.root", "READ");
  
  TH1D* hNominal = (TH1D*)Nominal->Get("H1D_Pt_Unfolded_w_MB");
  TH1D* hLowVar  = (TH1D*)Low_var->Get("H1D_Pt_Unfolded_w_MB");
  TH1D* hEndVar  = (TH1D*)End_var->Get("H1D_Pt_Unfolded_w_MB");

  TH1D* hBinningSys = (TH1D*)hNominal->Clone("hBinningSys");
  hBinningSys->Reset();
  hBinningSys->SetTitle("Systematic Uncertainty: Binning Choice; p_{T} (GeV/c); Rel. Uncertainty");

  double minPtCommon = hLowVar->GetXaxis()->GetBinLowEdge(1);
  double maxPtCommon = hEndVar->GetXaxis()->GetBinUpEdge(hEndVar->GetNbinsX());

  cout << "Evaluating systematic in common range: [" << minPtCommon << ", " << maxPtCommon << "] GeV" << endl;

  for (int i = 1; i <= hNominal->GetNbinsX(); ++i) {
      double binCenter = hNominal->GetBinCenter(i);
      double binLow    = hNominal->GetXaxis()->GetBinLowEdge(i);
      double binUp     = hNominal->GetXaxis()->GetBinUpEdge(i);

      if (binCenter < minPtCommon || binCenter > maxPtCommon) {
          cout << "Skipping Bin " << i << " (Outside common range)" << endl;
          hBinningSys->SetBinContent(i, 0); 
          continue;
      }

      double valNom = hNominal->GetBinContent(i);
      if (valNom <= 0) continue; 

      // 3. Find contents using coordinates (FindBin)
      // FindBin is safe here because we already verified binCenter is within range
      double valLow = hLowVar->GetBinContent(hLowVar->FindBin(binCenter));
      double valEnd = hEndVar->GetBinContent(hEndVar->FindBin(binCenter));

      double diffLow = std::abs(valNom - valLow) / valNom;
      double diffEnd = std::abs(valNom - valEnd) / valNom;

      double maxDiff = std::max(diffLow, diffEnd);

      hBinningSys->SetBinContent(i, maxDiff);
      hBinningSys->SetBinError(i, 0); 
      
      cout << "Bin " << i << " [" << binCenter << " GeV]: Rel Sys = " << maxDiff * 100 << "%" << endl;
  }

  return hBinningSys;

}

void Draw_Binning_systematic(){
  TH1D* hBinningSys = Get_Binning_systematic();
  TFile* outFile = new TFile("sys_Binning.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  hBinningSys->Write("Rel_Sys_Binning");
  outFile->Close();


  hBinningSys->Scale(100);
  TCanvas* c1 = new TCanvas("c1", "Binning Systematic", 800, 700);
  gPad->SetGridy(); // Grid helps see the percentage levels clearly
  gPad->SetLeftMargin(0.15);

  // 1. Style the systematic histogram
  hBinningSys->SetStats(0); // Remove the stats box
  hBinningSys->SetTitle("Systematic Uncertainty: Binning Choice");
  hBinningSys->GetYaxis()->SetTitle("Rel. Uncertainty (|Nom - Var|/Nom)");
  hBinningSys->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hBinningSys->GetYaxis()->SetTitleOffset(1.5);

  // Set Y-axis range from 0 to slightly above your max error (e.g., 20% or 0.2)
  hBinningSys->SetMinimum(0.0);
  hBinningSys->SetMaximum(hBinningSys->GetMaximum() * 1.5); 

  // 2. Draw as a shaded area (The "Band")
  hBinningSys->SetFillColorAlpha(kGray, 0.4); // Light gray transparent band
  hBinningSys->SetFillStyle(1001);
  hBinningSys->Draw("HIST"); // Draw the shaded histogram first

  // 3. Draw a bold step line on top
  hBinningSys->SetLineColor(kBlue+1);
  hBinningSys->SetLineWidth(3);
  hBinningSys->Draw("HIST SAME"); // "HIST" draws the step line without markers

  // 4. Optional: Draw markers if you still want to see the bin centers
  hBinningSys->SetMarkerStyle(20);
  hBinningSys->SetMarkerSize(0.8);
  hBinningSys->SetMarkerColor(kBlack);
  hBinningSys->Draw("P SAME");

}

// double Get_zVertex_reconstruction_efficiency(int iDataset, int iRadius){
//   // Retrieve histogram
//   TH1D* H1D_collisions_zvertex = (TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData + "/h_collisions_zvertex");

//   if (!H1D_collisions_zvertex) {
//     std::cerr << "Error: z-vertex histogram not found!" << std::endl;
//     return -1.;
//   }

//   // Clone to avoid modifying original
//   H1D_collisions_zvertex = (TH1D*)H1D_collisions_zvertex
//     ->Clone(Form("hZvertex_clone_dataset%d_R%d", iDataset, iRadius));

//   // Define Gaussian fit
//   TF1* fGaus = new TF1("fGaus", "gaus", -10, 10);

//   // Fit histogram (quiet mode)
//   H1D_collisions_zvertex->Fit(fGaus, "Q0");

//   // Extract parameters
//   double A     = fGaus->GetParameter(0); // amplitude
//   double mean  = fGaus->GetParameter(1);
//   double sigma = fGaus->GetParameter(2);

//   // Total integral of Gaussian from -inf to +inf
//   double totalIntegral = A * sigma * std::sqrt(2 * TMath::Pi());

//   // Integral between -10 and 10
//   double partialIntegral = fGaus->Integral(-10, 10);

//   // Efficiency
//   double efficiency = partialIntegral / totalIntegral;

//   return efficiency;
// }

// TH1D* Get_TVX_Eff(int iDataset, int iRadius){
//   // Get 2D histogram
//   TH2D* h2 = (TH2D*)file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get("jet-cross-section-efficiency/h2_jet_pt_part_eventselection");

//   if (!h2) {
//     std::cerr << "Error: 2D histogram not found!" << std::endl;
//     return nullptr;
//   }

//   // Find Y bins for kTVX and NColl
//   int bin_kTVX = -1;
//   int bin_NColl = -1;

//   for (int i = 1; i <= h2->GetYaxis()->GetNbins(); i++) {
//     TString label = h2->GetYaxis()->GetBinLabel(i);

//     if (label.Contains("kTVX"))  bin_kTVX = i;
//     if (label.Contains("NColl")) bin_NColl = i;
//   }

//   if (bin_kTVX < 0 || bin_NColl < 0) {
//     std::cerr << "Error: Could not find kTVX or NColl in Y axis labels!" << std::endl;
//     return nullptr;
//   }

//   // Project to X (pt) for both selections
//   TH1D* h_kTVX = h2->ProjectionX("h_kTVX", bin_kTVX, bin_kTVX);
//   TH1D* h_NColl = h2->ProjectionX("h_NColl", bin_NColl, bin_NColl);

//   // Compute ratio
//   TH1D* h_ratio = (TH1D*)h_kTVX->Clone(Form("TVX_Eff_dataset%d_R%d", iDataset, iRadius));
//   h_ratio->Divide(h_kTVX, h_NColl, 1.0, 1.0, "B"); // binomial errors

//   // --- Rebinning ---
//   TH1D* h_rebinned = (TH1D*)h_ratio->Rebin(nBinPtJetsGen[iRadius],Form("TVX_Eff_rebinned_dataset%d_R%d", iDataset, iRadius),ptBinsJetsGen[iRadius]);

//   return h_rebinned;
// }


TH1D* Get_UE_Systematic_Band(){
  TFile* Nominal_RCwoLJ = TFile::Open("../Matrices_Unfolding/output.root", "READ");
  TFile* Var_RC_Track         = TFile::Open("RTrackwoLJ/output.root", "READ");
  TFile* Var_RC         = TFile::Open("RC/output.root", "READ");

  TH1D* hNominal = (TH1D*)Nominal_RCwoLJ->Get("H1D_Pt_Unfolded_w_JJ");
  TH1D* hLowVar  = (TH1D*)Var_RC_Track->Get("H1D_Pt_Unfolded_w_JJ");
  TH1D* h_Rc     = (TH1D*)Var_RC->Get("H1D_Pt_Unfolded_w_JJ");

  if (!hNominal || !hLowVar) {
    std::cerr << "Error: missing histograms!" << std::endl;
    return nullptr;
  }

  TH1D* hUEsys = (TH1D*)hNominal->Clone("hUE_systematic");
  TH1D* hUEsys2 = (TH1D*)hNominal->Clone("hUE_systematic2");
  hUEsys->Reset();
  hUEsys2->Reset();
  int nBins = hNominal->GetNbinsX();
  for (int i = 1; i <= nBins; i++) {
    double N = hNominal->GetBinContent(i);
    if (N == 0) {
      hUEsys->SetBinContent(i, 0);
      continue;
    }
    double v1 = hLowVar->GetBinContent(i);
    double v2 = h_Rc->GetBinContent(i);
    double rel1 = fabs(v1 - N) / N;
    double rel2 = fabs(v2 - N) / N;
    hUEsys->SetBinContent(i, rel1);
    hUEsys2->SetBinContent(i, rel2);
  }
  // hUEsys->Scale(100);
  // hUEsys2->Scale(100);

  // TCanvas* c = new TCanvas("c", "UE systematic ", 800, 600);
  // hUEsys->Draw("E1");
  // hUEsys2->Draw("E1 same");
  TFile* outFile = new TFile("sys_UE.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  hUEsys->Write("Rel_Sys_UE_sub_RCwoLJ_RC_Tracks");
  hUEsys2->Write("Rel_Sys_UE_sub_RCwoLJ_RC");
  outFile->Close();

  return hUEsys2;
}

void Draw_UE_systematic_band(){
  TH1D* hUE = Get_UE_Systematic_Band();
  // hUE->Scale(100);

  if (!hUE) return;

  TCanvas* c1 = new TCanvas("c1", "UE systematic band", 800, 600);

  hUE->Draw("E1");
}

void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {
  cout << "########### Drawing systematics from parameter variation ###############" << endl;
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);

  TH1D* H1D_jetPt_unfolded[nUnfoldIteration];

  TString partialUniqueSpecifier;

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

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

      hSys_rel->SetBinContent(ib, rel); 
  }

  TFile* outFile = new TFile("sys_iteration.root", "UPDATE");   // "UPDATE" Open existing or create if missing / "RECREATE" Always deletes file and creates new one
  hSys_rel->Write("Rel_Sys_iteration");
  outFile->Close();

  hSys_rel->Scale(100);
  TString* pdfName_envolope = new TString("Systematics_deviation_ParameterVariation_"+partialUniqueSpecifier);
  TString* pdfName_relUnc = new TString("Systematics_RelativeUncertainty_ParameterVariation_"+partialUniqueSpecifier);
  TString textContext("SVD unfolding");
  TString* sigma = new TString("#sigma interation variation");
  TString* relativeErrors = new TString("realtive errors (%) ");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{10, 100}, {0, 3}}};
  // Draw_TH1_Histogram(hSys_envelope, textContext, pdfName_envolope, texPtJetRecX, sigma, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSys_rel, textContext, pdfName_relUnc, texPtJetRec, relativeErrors, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  
}

TH1D* Get_Total_Systematic(bool UE_trk_method){

  // =========================
  // LOAD FILES
  // =========================
  TFile* fTrackEff  = TFile::Open("Track_eff_sys/20260523_TrackingEff/sys_TrackEff.root");
  // TFile* fUE        = TFile::Open("Unfolding_UE_sys/sys_UE.root");
  TFile* fMethod    = TFile::Open("Method_sys/sys_Method.root");
  TFile* fBinning   = TFile::Open("Binning_sys/sys_Binning.root");
  TFile* fIteration = TFile::Open("Iteration_sys/sys_iteration.root");
  TFile* fSecondary = TFile::Open("Secondary_sys/sys_Secondary.root");

  // =========================
  // LOAD HISTOS
  // =========================
  TH1D* hTrackEff = (TH1D*) fTrackEff->Get("Rel_TrackEff_svd");
  
  // TH1D* hUE_tracks = (TH1D*) fUE->Get("Rel_Sys_UE_sub_RCwoLJ_RC_Tracks"); 
  // TH1D* hUE_RC     = (TH1D*) fUE->Get("Rel_Sys_UE_sub_RCwoLJ_RC");

  TH1D* hMethod  = (TH1D*) fMethod->Get("Rel_Sys_Method");
  TH1D* hBinning = (TH1D*) fBinning->Get("Rel_Sys_Binning");
  TH1D* hIter    = (TH1D*) fIteration->Get("Rel_Sys_iteration");
  TH1D* hSec    = (TH1D*) fSecondary->Get("Rel_Sys_secondary");

  int binStart = hTrackEff->FindBin(10.0);
  int binEnd   = hTrackEff->FindBin(100.0)-1;
  int nbins = hTrackEff->GetNbinsX();

  RewriteHistogramRange(hTrackEff, binStart, binEnd);
  RewriteHistogramRange(hMethod, binStart, binEnd);
  RewriteHistogramRange(hBinning, binStart, binEnd);
  RewriteHistogramRange(hIter, binStart, binEnd);
  RewriteHistogramRange(hSec, binStart, binEnd);

  // =========================
  // CREATE EXTRA HISTOS
  // =========================

  // // Secondary tracks (propagated)
  // TH1D* hSec = (TH1D*) hTrackEff->Clone("hSecondary");
  // hSec->Reset();

  // Closure test (NOT propagated)
  // TH1D* hClosure = (TH1D*) hTrackEff->Clone("hClosure");
  // hClosure->Reset();

  // Cross section (NOT propagated)
  TH1D* hNorm = (TH1D*) hTrackEff->Clone("hNorm");
  hNorm->Reset();

  for (int i = 1; i <= nbins; i++) {

    // double pt = hTrackEff->GetBinCenter(i);

    // // Secondary: 0.022 → 0.025 up to 30 GeV, then flat
    // double sec;
    // if (pt < 30.0) {
    //   sec = 0.022 + (0.025 - 0.022) * (pt / 30.0);
    // } else {
    //   sec = 0.025;
    // }

    // hSec->SetBinContent(i, sec);

    // // Closure test: constant 0.5%
    // hClosure->SetBinContent(i, 0.005);

    // Cross section: constant 4.5%
    hNorm->SetBinContent(i, 0.045);
  }

  // =========================
  // TOTAL HISTOS
  // =========================
  TH1D* hTot_tracks = (TH1D*) hTrackEff->Clone("hTot_tracks");
  hTot_tracks->Reset();

  // TH1D* hTot_RC = (TH1D*) hTrackEff->Clone("hTot_RC");
  // hTot_RC->Reset();

  // =========================
  // COMPUTE TOTAL
  // =========================
  for (int i = 1; i <= nbins; i++) {

    double t  = hTrackEff->GetBinContent(i);
    double m  = hMethod->GetBinContent(i);
    double b  = hBinning->GetBinContent(i);
    double it = hIter->GetBinContent(i);

    // double ue_tr = hUE_tracks->GetBinContent(i);
    // double ue_rc = hUE_RC->GetBinContent(i);

    double sec = hSec->GetBinContent(i);

    // include secondary in quadrature
    double tot_tr = sqrt(t*t + m*m + b*b + it*it + sec*sec);//+ ue_tr*ue_tr 
    // double tot_rc = sqrt(t*t + m*m + b*b + it*it + ue_rc*ue_rc + sec*sec);

    hTot_tracks->SetBinContent(i, tot_tr);
    // hTot_RC->SetBinContent(i, tot_rc);
  }

  
  // =========================
  // STYLE SETTINGS
  // =========================
  gStyle->SetOptStat(0);

  // =========================
  // SELECT HISTOS ONCE
  // =========================
  // TH1D* hUE  = UE_trk_method ? hUE_tracks : hUE_RC;
  // TH1D* hTot = UE_trk_method ? hTot_tracks : hTot_RC;
  TH1D* hTot =  hTot_tracks;


  // =========================
  // SMALL STYLE HELPER
  // =========================
  auto SetStyle = [](TH1D* h, int color, int width=2, int style=1){
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetLineStyle(style);
    // h->Scale(100);
  };

  // =========================
  // STYLE
  // =========================
  SetStyle(hTrackEff, kBlue+1);
  // SetStyle(hUE,       kGreen+2);
  SetStyle(hMethod,   kMagenta+1);
  SetStyle(hBinning,  kOrange+1);
  SetStyle(hIter,     kCyan+2);

  // Secondary
  SetStyle(hSec, kRed+1);

  // Not propagated
  // SetStyle(hClosure, kGray+2);
  SetStyle(hNorm,    kBlack, 3, 3); // dotted

  // Total (dominant)
  SetStyle(hTot, kBlack, 4);

  // =========================
  // CANVAS + AXIS
  // =========================
  TCanvas* c = new TCanvas("cSys","Systematics",900,700);
  // hTot->GetXaxis()->SetRangeUser(10.0, 100.0);

  hTot->SetTitle("Relative Systematic Uncertainties");
  hTot->GetYaxis()->SetTitle("Relative uncertainty");
  hTot->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV/#it{c})");

  double max = hTot->GetMaximum();
  hTot->SetMaximum(max * 1.1);
  hTot->SetMinimum(0.0);

  // =========================
  // DRAW
  // =========================
  hTot->Draw("HIST");

  hTrackEff->Draw("HIST SAME");
  // hUE->Draw("HIST SAME");
  hMethod->Draw("HIST SAME");
  hBinning->Draw("HIST SAME");
  hIter->Draw("HIST SAME");
  hSec->Draw("HIST SAME");

  // hClosure->Draw("HIST SAME");
  hNorm->Draw("HIST SAME");

  // redraw total on top
  hTot->Draw("HIST SAME");

  // =========================
  // LEGEND
  // =========================
  TLegend* leg = new TLegend(0.55,0.55,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(hTot,"Total","l");
  leg->AddEntry(hTrackEff,"Tracking efficiency","l");

  // if (UE_trk_method)
  //   leg->AddEntry(hUE,"UE (RCwoLJ vs RCtrkswoLJ)","l");
  // else
  //   leg->AddEntry(hUE,"UE (RCwoLJ vs RC)","l");

  leg->AddEntry(hMethod,"Unfolding method","l");
  leg->AddEntry(hBinning,"Binning","l");
  leg->AddEntry(hIter,"Iteration","l");
  leg->AddEntry(hSec,"Secondary particles","l");
  // leg->AddEntry(hClosure,"Closure 0.5% (not in total)","l");
  leg->AddEntry(hNorm,"Normalization 4.5% (not in total)","l");

  leg->Draw();

  // =========================
  // LABEL (optional)
  // =========================
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  // latex.DrawLatex(0.15,0.92,"ALICE simulation");

  c->Update();
  return hTot_tracks;
}

void RewriteHistogramRange(TH1D*& h, int binStart, int binEnd){
    if (!h) {
        std::cerr << "Error: Input histogram is null." << std::endl;
        return;
    }

    int nBins = h->GetNbinsX();

    if (binStart < 1 || binEnd > nBins || binStart > binEnd) {
        std::cerr << "Error: Invalid bin range." << std::endl;
        return;
    }

    int newNBins = binEnd - binStart + 1;

    // Store new bin edges
    std::vector<double> newBins(newNBins + 1);
    for (int i = 0; i <= newNBins; ++i) {
        newBins[i] = h->GetXaxis()->GetBinLowEdge(binStart + i);
    }

    // Create temporary histogram
    TH1D* hTemp = new TH1D(
        Form("%s_temp", h->GetName()),
        h->GetTitle(),
        newNBins,
        newBins.data()
    );

    // Copy contents and errors
    for (int i = 1; i <= newNBins; ++i) {
        int oldBin = binStart + i - 1;
        hTemp->SetBinContent(i, h->GetBinContent(oldBin));
        hTemp->SetBinError(i, h->GetBinError(oldBin));
    }

    // Preserve axis titles
    hTemp->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    hTemp->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

    // Delete old histogram
    delete h;

    // Replace pointer
    h = hTemp;
}