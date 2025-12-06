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
#include <RooUnfold.h> // one should likely do `aliBuild build RooUnfold` then `alienv enter RooUnfold/latest` as alidist roounfold version is usually quite old
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"


#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include "TMath.h"
#include "TRandom3.h"
#include "TProfile.h"

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
#include <stdlib.h>     /* abort, NULL */
#include "TRandom.h"

using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step);

void unfoldSpec(RooUnfoldResponse* response, TH1D* measured, TH1D** hist_unfold, int unfoldParameterInput, int errTreatment, int iRadius, bool doBayes);
void smearResponseTH2D(TH2D* hOrig, TH2D* hSmeared);
void resetErrorsTH2D(TH2D* h);
void reweightHistG(TH1D* histG, TH1D* histOrg, TH1D* histSmear);
void smearTH1(TH1D* hOrig, TH1D* hSmeared);
void GetResponse(RooUnfoldResponse** respJetPt, TH2D** Response_fine_Smeared_out, bool smearResp, int iDataset, int iRadius, bool resetErr );
void FillVariations(TH1D* hist, TProfile* histP);
void GetJetPurity(TH1D** H1D_jetPurity, TH2D* Response_fine_Smeared, TH1D* mcdRecBin, int iRadius);
void GetJetEfficiency(TH1D** H1D_jetEfficiency, TH2D* Response_fine_Smeared, TH2D* Response_fine_Org, int iDataset, int iRadius);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void unfoldSpec_statVariations_Ahmad(){

  int iDataset = 0;
  int iRadius = 1;
  int nRepeats = 100;
  bool doBayes = true;
  bool varySpec = true;
  bool varyResp = false;
  int errTreatment = 1;

  int nIterSpecBayes = 4;
  int nIterSpecSVD   = 7;
  int nIterSpec = doBayes ? nIterSpecBayes : nIterSpecSVD;

  // ##### Original Unfolding #####

  // --- Get Original Response Matrix (roounfold object + fine binned TH2D)
  RooUnfoldResponse* OrgResponse;
  TH2D* TH2D_Response_fine; // fine binning raw event matrix
  bool smearResp = false;
  bool resetErr = false;
  GetResponse(&OrgResponse, &TH2D_Response_fine, smearResp, iDataset, iRadius, resetErr);

  if (!TH2D_Response_fine) {
      cout << "TH2D_Response_fine is NULL!" << endl;
  } else {
      cout << "############# Jet TH2D_Response_fine exist #############" << endl;
  }

  // --- Get Measured Spectrum (rec binning, pre width scaling)
  TH1D* measuredInput = nullptr;
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, ""); //Good 
  
  // --- Correct measured for Purity 
  TH1D* mcdRecBin = nullptr;
  TH1D* purity = nullptr;
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(mcdRecBin, iDataset, iRadius, ""); // Good 
  GetJetPurity(&purity, TH2D_Response_fine, mcdRecBin, iRadius); // Good
  if (!purity) {
      cout << "purity is NULL!" << endl;
  } else {
      cout << "############# Jet purity exist #############" << endl;
  }

  TH1D* measured = (TH1D*) measuredInput->Clone("measured_correctedForPurity");
  measured->Multiply(purity);  // Good 
  if (!measured) {
        cout << "measured is NULL, before the for loop!" << endl;
    } else {
        cout << "############# Jet measured exist #############" << endl;
  }

  // --- Unfold Priginal spectrum
  TH1D* TH1D_Unfolded_original = nullptr;
  unfoldSpec(OrgResponse, measured, &TH1D_Unfolded_original, nIterSpec, errTreatment, iRadius, doBayes); 
  if (!TH1D_Unfolded_original) {
      cout << "TH1D_Unfolded_original is NULL, before the for loop!" << endl;
    } else {
      cout << "############# Jet TH1D_Unfolded_original exist #############" << endl;
  }

  // // --- Kinematic efficiency
  // TH1D* kinematicEfficiency;
  // TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  // TH2D* Response_normYslice = (TH2D*) TH2D_Response_fine->Clone("response normYslice");
  // NormaliseYSlicesToOne(Response_normYslice); 
  // Get_ResponseMatrix_Pt_KinematicEffiency(kinematicEfficiency, Response_normYslice, name_H1D_kinematicEfficiency, iRadius); // Good
  // if (!kinematicEfficiency) {
  //     cout << "kinematicEfficiency is NULL!" << endl;
  // } else {
  //     cout << "############# Jet kinematicEfficiency exist #############" << endl;
  // }

  // --- Matched efficiency
  TH1D* jetEfficiency = nullptr;
  TH2D* TH2D_Response_fine_clone = (TH2D*) TH2D_Response_fine->Clone("TH2D_Response_fine_clone");
  GetJetEfficiency(&jetEfficiency, TH2D_Response_fine, TH2D_Response_fine_clone, iDataset, iRadius);  // Good
  if (!jetEfficiency) {
      cout << "jetEfficiency is NULL!" << endl;
  } else {
      cout << "############# Jet Efficiency exist #############" << endl;
  }
    
  
  TH1D_Unfolded_original->Divide(jetEfficiency);
  NormaliseRawHistToNEvents(TH1D_Unfolded_original, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
  TransformRawHistToYield(TH1D_Unfolded_original); // final unfolded spectrum

  if (!TH1D_Unfolded_original) {
        cout << "TH1D_Unfolded_original is NULL, before the for loop!" << endl;
    } else {
        cout << "############# Jet TH1D_Unfolded_original exist #############" << endl;
  }

  // ---------------------------
  // smear + unfold
  
  TProfile* TProfile_JetPtVar = new TProfile("TProfile_JetPtVar","",nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius],"S");   
    
  for(int rep=0; rep < nRepeats; rep++){
    cout<<" ########################################################################################### " <<endl;
    cout<<" ############################## UNFOLD REPETATION "<<rep << " ############################## " <<endl;
    cout<<" ########################################################################################### " <<endl;

    TH1D* TH1D_Unfolded = nullptr; 
    RooUnfoldResponse* SmearedResponse = nullptr;
    TH2D* Response_fine_Smeared = nullptr;
    TH1D* measuredSmeared = nullptr;
    
    if(varyResp){
      bool smearResp = true;
      GetResponse(&SmearedResponse, &Response_fine_Smeared, smearResp, iDataset, iRadius, resetErr);
    }
    
    if(varySpec){
      measuredSmeared = (TH1D*) measuredInput->Clone("measuredSmeared");
      measuredSmeared->Reset();
      smearTH1(measuredInput,measuredSmeared);
      TH1D* purityS = nullptr;
      GetJetPurity(&purityS, varyResp ? Response_fine_Smeared : TH2D_Response_fine, mcdRecBin, iRadius);
      measuredSmeared->Multiply(purityS);
    }

    // Choose inputs for unfolding (pointers)
    RooUnfoldResponse* response = varyResp ? SmearedResponse : OrgResponse;
    TH1D* spec = varySpec ? measuredSmeared : measured;
    // Debug
    if(!spec) cout << "ERROR: spec is NULL!" << endl;
    unfoldSpec(response, spec, &TH1D_Unfolded, nIterSpec, errTreatment, iRadius, doBayes);

    
    TH1D* jetEfficiency = nullptr;
    GetJetEfficiency(&jetEfficiency, (varyResp ? Response_fine_Smeared : TH2D_Response_fine) , TH2D_Response_fine, iDataset, iRadius);

    TH1D_Unfolded->Divide(jetEfficiency);

    NormaliseRawHistToNEvents(TH1D_Unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    TransformRawHistToYield(TH1D_Unfolded); // final smeared unfolded spectrum


    FillVariations(TH1D_Unfolded, TProfile_JetPtVar); 

    delete TH1D_Unfolded;
    delete Response_fine_Smeared;    
  }

  TProfile* TProfile_MatrixVar = new TProfile("TProfile_MatrixVar","",nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius],"S"); 
  for(int rep=0; rep < nRepeats; rep++){
    cout<<" ################################################################################################### " <<endl;
    cout<<" ############################## MATRIX SMEARED REPETATION "<<rep << " ############################## " <<endl;
    cout<<" ################################################################################################### " <<endl;

    TH1D* TH1D_Unfolded = nullptr; 
    RooUnfoldResponse* SmearedResponse = nullptr;
    TH2D* Response_fine_Smeared = nullptr;
    TH1D* measuredSmeared = nullptr;
    GetResponse(&SmearedResponse, &Response_fine_Smeared, true, iDataset, iRadius, resetErr); // true to smear
    unfoldSpec(SmearedResponse, measured, &TH1D_Unfolded, nIterSpec, errTreatment, iRadius, doBayes); // unfold with smeared response and nominal measured
    TH1D* jetEfficiency = nullptr;
    GetJetEfficiency(&jetEfficiency, Response_fine_Smeared, TH2D_Response_fine, iDataset, iRadius);
    TH1D_Unfolded->Divide(jetEfficiency);
    NormaliseRawHistToNEvents(TH1D_Unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    TransformRawHistToYield(TH1D_Unfolded); // final smeared unfolded spectrum
    FillVariations(TH1D_Unfolded, TProfile_MatrixVar); 

    delete TH1D_Unfolded;
    delete Response_fine_Smeared;    
  }

  // --------------------------------------

  // unfolded spec errors as assigned by RooUnfold
  // TH1D* TH1D_Unfolded_original_error = (TH1D*) TH1D_Unfolded_original->Clone("TH1D_Unfolded_original_error");
  // TH1D_Unfolded_original_error->Reset();
  
  // for(int bin=1; bin<=TH1D_Unfolded_original->GetNbinsX(); bin++){ 
  //   if(TH1D_Unfolded_original->GetBinContent(bin)){
  //     TH1D_Unfolded_original_error->SetBinContent(bin,TH1D_Unfolded_original->GetBinError(bin)/TH1D_Unfolded_original->GetBinContent(bin));
  //     TH1D_Unfolded_original_error->SetBinError(bin,0);
  //   }
  // }

  
  // // error from variations 
  // TH1D* hRelUncert = (TH1D*) TH1D_Unfolded_original->Clone("hRelUncert");
  // hRelUncert->Reset();

  // int nBins = hRelUncert->GetNbinsX();

  // for (int bin = 1; bin <= nBins; bin++) {
  //     double mean = TProfile_JetPtVar->GetBinContent(bin);
  //     double rms  = TProfile_JetPtVar->GetBinError(bin);  // because of option "S"

  //     if (mean > 0) {
  //         hRelUncert->SetBinContent(bin, rms / mean);
  //     } else {
  //         hRelUncert->SetBinContent(bin, 0);
  //     }
  //     hRelUncert->SetBinError(bin, 0);
  // }

  // TH1D* hRelUncertSmearedMatrix = (TH1D*) TH1D_Unfolded_original->Clone("hRelUncert");
  // hRelUncertSmearedMatrix->Reset();

  // int nBins2 = hRelUncertSmearedMatrix->GetNbinsX();

  // for (int bin = 1; bin <= nBins2; bin++) {
  //     double mean = TProfile_MatrixVar->GetBinContent(bin);
  //     double rms  = TProfile_MatrixVar->GetBinError(bin);  // because of option "S"

  //     if (mean > 0) {
  //         hRelUncertSmearedMatrix->SetBinContent(bin, rms / mean);
  //     } else {
  //         hRelUncertSmearedMatrix->SetBinContent(bin, 0);
  //     }
  //     hRelUncertSmearedMatrix->SetBinError(bin, 0);
  // }
  // === 1) Define / clone the histograms ===

  TH1D* TH1D_Unfolded_original_error = (TH1D*) TH1D_Unfolded_original->Clone("TH1D_Unfolded_original_error");
  TH1D_Unfolded_original_error->Reset();

  TH1D* hRelUncert = (TH1D*) TH1D_Unfolded_original->Clone("hRelUncert");
  hRelUncert->Reset();

  TH1D* hRelUncertSmearedMatrix = (TH1D*) TH1D_Unfolded_original->Clone("hRelUncertSmearedMatrix");
  hRelUncertSmearedMatrix->Reset();

  // === 2) Combined single loop for all three calculations ===
  int nBins = TH1D_Unfolded_original->GetNbinsX();
  for (int bin = 1; bin <= nBins; bin++) {
      // --------------------------------------------------------
      // A) Original relative error (statistical)
      // --------------------------------------------------------
      double val = TH1D_Unfolded_original->GetBinContent(bin);
      double err = TH1D_Unfolded_original->GetBinError(bin);

      if (val > 0)
          TH1D_Unfolded_original_error->SetBinContent(bin, err / val);
      else
          TH1D_Unfolded_original_error->SetBinContent(bin, 0);

      TH1D_Unfolded_original_error->SetBinError(bin, 0);

      // --------------------------------------------------------
      // B) Relative uncertainty from spectrum variations (TProfile_JetPtVar)
      // --------------------------------------------------------
      double mean_spec = TProfile_JetPtVar->GetBinContent(bin);   // ⟨x⟩
      double rms_spec  = TProfile_JetPtVar->GetBinError(bin);     // RMS ("S" option)

      if (mean_spec > 0)
          hRelUncert->SetBinContent(bin, rms_spec / mean_spec);
      else
          hRelUncert->SetBinContent(bin, 0);

      hRelUncert->SetBinError(bin, 0);

      // --------------------------------------------------------
      // C) Relative uncertainty from smeared matrix variations
      // --------------------------------------------------------
      double mean_mat = TProfile_MatrixVar->GetBinContent(bin);
      double rms_mat  = TProfile_MatrixVar->GetBinError(bin);

      if (mean_mat > 0)
          hRelUncertSmearedMatrix->SetBinContent(bin, rms_mat / mean_mat);
      else
          hRelUncertSmearedMatrix->SetBinContent(bin, 0);

      hRelUncertSmearedMatrix->SetBinError(bin, 0);
  }


  // ----- Drawing -----
  // TString* pdfName_realtiveError = new TString("realtive Error, SVD unfolding");
  // TString* texPtUnfol = new TString("#it{p}_{T}^{unf} (GeV/#it{c})");
  // TString* relativeErrors = new TString("realtive errors");
  // TString textContext = TString("SVD ");
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{0, 100},{-999, -999}}};
  // std::array<std::array<float, 2>, 2> legendPlacement = {{{0.15, 0.75}, {0.45, 0.88}}}; // {{{x1, y1}, {x2, y2}}}

  // // Determine the third legend dynamically
  // const TString LegendRelativeErrors[3] = {"RooUnfold orginal error", "Smeared measured, Nominal Matrix", "Nominal measured, Smeared Matrix"};

  // TH1D** RelativeErrorsSVD = new TH1D*[3];
  // RelativeErrorsSVD[0] = TH1D_Unfolded_original_error;
  // RelativeErrorsSVD[1] = hRelUncert;
  // RelativeErrorsSVD[2] = hRelUncertSmearedMatrix;
  
  // Draw_TH1_Histograms(RelativeErrorsSVD, LegendRelativeErrors, 3, textContext, pdfName_realtiveError, texPtUnfol, relativeErrors, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacementAuto, "logy");

  // Put histograms in an array
  TH1D* RelativeErrorsSVD[3] = {TH1D_Unfolded_original_error, hRelUncert, hRelUncertSmearedMatrix };

  Color_t colors[3] = { kBlack, kRed, kBlue };

  TString labels[3] = {
      "RooUnfold error",           // later append error treatment
      "Smeared measured, Nominal Matrix",
      "Nominal measured, Smeared Matrix"
  };

  // === Build the full legend name for the first histogram ===
  TString etName;
  switch (errTreatment) {
      case -1: etName = "kCovariance"; break;
      case  1: etName = "kErrors";     break;
      case  2: etName = "kNoError";    break;
      case  3: etName = "kCovToy";     break;
      default: etName = "kErrors";     break;
  }

  labels[0] += " (" + etName + ")";


  // === Canvas ===
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "Relative Errors SVD", 800, 800);
  c->SetLogy();

  // === Loop over histograms to configure and draw ===
  for (int i = 0; i < 3; i++) {
      RelativeErrorsSVD[i]->SetStats(0);
      RelativeErrorsSVD[i]->SetLineColor(colors[i]);
      RelativeErrorsSVD[i]->SetLineWidth(2);

      RelativeErrorsSVD[i]->SetTitle("Relative Uncertainties; p_{T}^{jet} [GeV]; Relative error");
      if (i == 0)
          RelativeErrorsSVD[i]->Draw("HIST");
      else
          RelativeErrorsSVD[i]->Draw("HIST SAME");
  }

  // === Legend ===
  TLegend* leg = new TLegend(0.15, 0.75, 0.6, 0.88);

  for (int i = 0; i < 3; i++)
      leg->AddEntry(RelativeErrorsSVD[i], labels[i], "l");

  leg->Draw();

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

void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step){
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);
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
////////////////////////////////////////////////////////////////////////////////////// functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ------------------------------------------------------------
void unfoldSpec(RooUnfoldResponse* response, TH1D* measured, TH1D** hist_unfold, int unfoldParameterInput, int errTreatment, int iRadius, bool doBayes = kTRUE){
  
  if (!response) {
      cout << "unfoldSpec ERROR: Response Matrix pointer is NULL!" << endl;
      return;
    }

  if (!measured) {
      cout << "unfoldSpec ERROR: Measured histogram pointer is NULL!" << endl;
      return;
    }

  if (!hist_unfold) {
      cout << "unfoldSpec ERROR: HistUnfold pointer is NULL!" << endl;
      return;
    }

  // --- prepare a RooUnfold object (base pointer) ---
  RooUnfold* unfold = nullptr;
  RooUnfoldBayes* unfoldBayes = nullptr;
  RooUnfoldSvd*   unfoldSvd    = nullptr;

  // // unfold spectrum : default Bayes
  // unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);
  // unfold = unfoldBayes; 

  if (!doBayes) {
    int unfoldParameterSvdInitial = 1;
    RooUnfoldSvd* unfoldSvdInit = new RooUnfoldSvd(response, measured, unfoldParameterSvdInitial); 
    unfoldSvdInit->Hreco(); // necessary to have GetD() give a meaningful output
    TSVDUnfold *tsvdUnfold = (TSVDUnfold*)unfoldSvdInit->Impl();

    // Optionally draw D distribution if requested (kept from your original)
    if (tsvdUnfold) {
      TH1D* H1D_D = tsvdUnfold->GetD(); // may be nullptr if not available
      if (H1D_D) {
        TString inputUnfoldingName = "";
        TString* pdfName_regparam = new TString("Svd_regularisationd_distribution");
        TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));
        std::array<std::array<float, 2>, 2> drawnWindowSvdParam = {{{0, 30}, {0.01, 10000}}}; // {{xmin, xmax}, {ymin, ymax}}
        // Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdK, texSvdDvector, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");
      }
    }
    // Now create the real SVD object with the chosen regularisation parameter
    unfoldSvd = new RooUnfoldSvd(response, measured, unfoldParameterInput);
    unfold = unfoldSvd;
    cout<<" SVD unfolding, nIter "<<unfoldParameterInput<<endl;
    
    //delete unfoldSvdInit;
  } else if (doBayes) {
    unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);
    unfold = unfoldBayes;
    cout<<" Bayes unfolding, nIter "<<unfoldParameterInput<<endl;
  }

  if (!unfold) {
        cerr << "unfoldSpec: failed to create RooUnfold object." << endl;
        // clean up whatever was created
        delete unfoldBayes;
        delete unfoldSvd;
        return;
  }

  cout<<" unfold spectrum "<<endl;
    
  // --- Perform unfolding, choose error treatment ---
  RooUnfold::ErrorTreatment et = RooUnfold::kErrors;
  switch (errTreatment) {
    case -1: et = RooUnfold::kCovariance; break;
    case  1: et = RooUnfold::kErrors;     break;
    case  2: et = RooUnfold::kNoError;    break;
    case  3: et = RooUnfold::kCovToy;     break;
    default: et = RooUnfold::kErrors;     break;
  }

  // Hreco returns a TH1* (actually a TH1D* for 1D). RooUnfold owns this histogram? Typically
  // it returns a freshly allocated histogram (check your RooUnfold version). We capture it here.
  TH1D* hReco = (TH1D*) unfold->Hreco(et);    // This is the unfolded histogram
  if (!hReco) {
      cout << "unfoldSpec: Hreco returned nullptr." << endl;
      delete unfold; // delete concrete RooUnfold
      return;
  }

  *hist_unfold = (TH1D*) hReco->Clone("unfolded histogram before efficiencies correction and normalization (bin widht and event)");  // Important: clone it!

  // delete unfold;
  cout<<" ############# unfolding done properly with unfoldSpec function #################"<<endl;
}
// ------------------------------------------------------------

void smearResponseTH2D(TH2D* hOrig, TH2D* hSmeared) {
    int nBinsX = hOrig->GetNbinsX();
    int nBinsY = hOrig->GetNbinsY();

    cout << "smearResponseTH2D, nBinsX: " << nBinsX << endl;
    cout << "smearResponseTH2D, nBinsY: " << nBinsY << endl;

    for (int iY = 1; iY <= nBinsY; iY++) { // should start from 1 to skip underflow?
        if (!(iY % 100)) cout << "iY: " << iY << endl;

        for (int iX = 1; iX <= nBinsX; iX++) { // should start from 1 to skip underflow?
            double resp = hOrig->GetBinContent(iX, iY);
            double err = hOrig->GetBinError(iX, iY);
            double relErr = (resp != 0) ? err / resp : 0.0;

            // Gaussian smearing
            double contSmeared = gRandom->Gaus(resp, err);
            //double errSmeared = relErr * contSmeared;

            if (contSmeared < 0) { // keep original if negative
                contSmeared = resp;
                // errSmeared = err;
            }

            hSmeared->SetBinContent(iX, iY, contSmeared);
            hSmeared->SetBinError(iX, iY, err);
        }
    }
    cout << " ############# smearResponseTH2D done properly #################" << endl;
}
// ------------------------------------------------------------

void resetErrorsTH2D(TH2D* h) {
    int nBinsX = h->GetNbinsX();
    int nBinsY = h->GetNbinsY();

    cout << "resetErrorsTH2D, nBinsX: " << nBinsX << endl;
    cout << "resetErrorsTH2D, nBinsY: " << nBinsY << endl;

    for (int iY = 1; iY <= nBinsY; iY++) {
        for (int iX = 1; iX <= nBinsX; iX++) {
            h->SetBinError(iX, iY, 0.0);
        }
    }
    cout << " ############# resetErrorsTH2D done properly #################" << endl;
}
// ------------------------------------------------------------

void reweightHistG(TH1D* histG, TH1D* histOrg, TH1D* histSmear) {
    int nBins = histG->GetXaxis()->GetNbins();

    for (int bin = 0; bin <= nBins; bin++) { // skip underflow/overflow
        double contG = histG->GetBinContent(bin);
        double errG = histG->GetBinError(bin);
        double contOrg = histOrg->GetBinContent(bin);
        double contSmear = histSmear->GetBinContent(bin);

        if (contOrg != 0) {
            double weight = contSmear / contOrg;
            double contNew = contG * weight;
            double errNew = errG * weight;

            cout << "bin " << bin 
                 << " weight " << weight 
                 << " contG " << contG 
                 << " contNew " << contNew 
                 << " errG " << errG 
                 << " errNew " << errNew << endl;

            histG->SetBinContent(bin, contNew);
            histG->SetBinError(bin, errNew);
        }
    }

    cout << " ############# reweightHistG done properly #################" << endl;
}
// ------------------------------------------------------------

void smearTH1(TH1D* hist, TH1D* hSmeared) {
    hSmeared->Reset();

    int nBins = hist->GetNbinsX();
    for (int binx = 1; binx <= nBins; binx++) { // main bins
        double cont = hist->GetBinContent(binx);
        double err = hist->GetBinError(binx);

        double contSmeared = 0.0;
        double errCorr = 0.0;

        if (cont != 0.0) {
            contSmeared = gRandom->Gaus(cont, err);
            //errCorr = err * contSmeared / cont;
            errCorr = err; // keep original error
        }

        if (contSmeared < 0.0) { // keep original if negative
            contSmeared = cont;
            errCorr = err;
        }

        // cout << "hist: " << hist->GetName() 
        //      << " bin " << binx 
        //      << " x = " << hist->GetXaxis()->GetBinCenter(binx)
        //      << " cont = " << cont 
        //      << " err = " << err 
        //      << " contSmeared = " << contSmeared 
        //      << " errCorr = " << errCorr << endl;

        hSmeared->SetBinContent(binx, contSmeared);
        hSmeared->SetBinError(binx, errCorr);
    }
    cout << " ############# smearTH1 done properly #################" << endl;
}
// ------------------------------------------------------------

void GetResponse(RooUnfoldResponse** respJetPt, TH2D** Response_fine_Smeared_out, bool smearResp, int iDataset, int iRadius, bool resetErr = kFALSE){
  
  TH2D* H2D_jetPtMcdjetPtMcd;

  H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  
  TH2D* Response_fine_Smeared = nullptr;
  TH2D* Response_rebinned = nullptr;

  Response_fine_Smeared = (TH2D*) H2D_jetPtMcdjetPtMcd->Clone("Response_fine_Smeared");

  // --- smear response 
  if(smearResp){ 
    Response_fine_Smeared->Sumw2();
    smearResponseTH2D(H2D_jetPtMcdjetPtMcd, Response_fine_Smeared); // org, smeared (fine binning)
  }

  // fine binning response matrix for checking
  TH2D* Response_fine_Smeared_proba = (TH2D*) Response_fine_Smeared->Clone("Response_fine_Smeared_proba");
  NormaliseYSlicesToOne(Response_fine_Smeared_proba); // convert to probability matrix
  TH2D* H2D_response = (TH2D*)RebinVariableBins2D(Response_fine_Smeared_proba, nBinPtJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], ptBinsJetsFine[iRadius]).Clone("Get_PtResponseMatrix_detectorResponse_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]); // check the binning here

  // constructing the response matrix for unfolding
  Response_rebinned = (TH2D*) H2D_response->Clone("Response_smeared_rebinned");
  TH1D* priorSpectrumWeighting;
  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, ""); // get prior spectrum in fine binning (no width scaling, no event norm)

  WeightMatrixWithPrior(Response_rebinned, priorSpectrumWeighting, doUnfoldingPriorDivision); // Response_rebinned : Fine event weighted matrix
  MergeResponseMatrixBins(Response_rebinned, iDataset, iRadius, "noPriorMerging"); // Response_rebinned : rebined with no bin area scaling

  if(resetErr) resetErrorsTH2D(Response_rebinned);

  RooUnfoldResponse* response = new RooUnfoldResponse(0, 0, Response_rebinned); // create the roounfold matrix object (weighted event matrix and coarse rebinned)

  *respJetPt = response;
  
  *Response_fine_Smeared_out = Response_fine_Smeared;
  
  cout<<" ############# GetResponse done properly with GetResponse function #################"<<endl;
}
// ------------------------------------------------------------

void FillVariations(TH1D* hist, TProfile* histP) {
    // Check for consistency
    if (hist->GetNbinsX() != histP->GetNbinsX() ||
        hist->GetXaxis()->GetXmin() != histP->GetXaxis()->GetXmin() ||
        hist->GetXaxis()->GetXmax() != histP->GetXaxis()->GetXmax()) {
        cout << "FillVariations: discrepancy between hist and histP!" << endl;
        exit(1);
    }

    int nBins = hist->GetNbinsX();
    for (int binx = 1; binx <= nBins; binx++) {
        double cent = hist->GetXaxis()->GetBinCenter(binx);
        double cont = hist->GetBinContent(binx);
        histP->Fill(cent, cont); // chaque Fill dans un for (en bas) ne remplace pas la valeur, il construit une moyenne avec son erreur (\sigma/\sqrt(N)).
        //  Pour variance brute : double sigma = histP->GetBinError(binx) * sqrt(histP->GetBinEntries(binx));
    }

    cout << "############ FillVariations: callet properly #################" << endl;
}
// ------------------------------------------------------------

void GetJetPurity(TH1D** H1D_jetPurity, TH2D* Response_fine_Smeared, TH1D* mcdRecBin, int iRadius){
    cout << "############# GetJetPurity called! ###############" << endl;

    // --- Step 1: Project fine matrix onto reconstructed axis
    TH2D* TH2D_Response_fine_clone = (TH2D*) Response_fine_Smeared->Clone("TH2D_Response_fine_clone");
    TH2D_Response_fine_clone->Sumw2();
    TH1D* H1D_jetPt_smeared_projRec_RecBin = nullptr;
    H1D_jetPt_smeared_projRec_RecBin = (TH1D*)TH2D_Response_fine_clone->ProjectionX("jetPt_mcdMatched_smeared_forRecBin", 1, TH2D_Response_fine_clone->GetNbinsY(), "e");

    // --- Step 2: Rebin the projected histogram
    TH1D* H1D_jetPt_mcdMatched_recBinning = nullptr;
    H1D_jetPt_mcdMatched_recBinning = (TH1D*) H1D_jetPt_smeared_projRec_RecBin->Rebin(nBinPtJetsRec[iRadius], "jetPt_mcdMatched_recBinning_rebinned", ptBinsJetsRec[iRadius] );

    // TCanvas* c5 = new TCanvas("c5", "H1D_jetPt_mcdMatched_recBinning", 800, 600);
    // H1D_jetPt_mcdMatched_recBinning->SetTitle("H1D_jetPt_mcdMatched_recBinning;pT [GeV/c];H1D_jetPt_mcdMatched_recBinning");
    // H1D_jetPt_mcdMatched_recBinning->Draw();
    // for (int bin = 1; bin <= H1D_jetPt_mcdMatched_recBinning->GetNbinsX(); bin++) {
    //     cout << "bin " << bin << endl ;
    //     cout << " mcdMatched content " << H1D_jetPt_mcdMatched_recBinning->GetBinContent(bin) 
    //          << " error " << H1D_jetPt_mcdMatched_recBinning->GetBinError(bin) << endl;
    //     cout << " mcdRecBin content " << mcdRecBin->GetBinContent(bin) 
    //          << " mcdRecBin error " << mcdRecBin->GetBinError(bin) << endl;     
    // }

    // --- Step 3: Clone to create the purity histogram
    *H1D_jetPurity = (TH1D*) H1D_jetPt_mcdMatched_recBinning->Clone("H1D_jetPurity");

    // --- Step 4: Compute purity = matched / measured (binomial errors)
    (*H1D_jetPurity)->Divide(*H1D_jetPurity, mcdRecBin, 1., 1., "b");

    // --- Clean up temporary histograms
    delete H1D_jetPt_smeared_projRec_RecBin;
    delete H1D_jetPt_mcdMatched_recBinning;
    cout << " ############# GetJetPurity done properly with GetJetPurity function #################" << endl;
}
// ------------------------------------------------------------
void GetJetEfficiency(TH1D** H1D_jetEfficiency, TH2D* Response_fine_Smeared, TH2D* Response_fine_Org, int iDataset, int iRadius) {
  TH1D* H1D_jetPt_mcp = nullptr;
  TH1D* H1D_jetPt_mcpMatched = nullptr;

  Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, false, "");
  Get_Pt_spectrum_mcpMatched_genBinning(H1D_jetPt_mcpMatched, iDataset, iRadius, "");

  if (!H1D_jetPt_mcp || !H1D_jetPt_mcpMatched) {
        std::cerr << "GetJetEfficiency ERROR: input histograms are null!" << std::endl;
        *H1D_jetEfficiency = nullptr;
        return;
    }

    // Clone into output histogram
    *H1D_jetEfficiency = (TH1D*) H1D_jetPt_mcpMatched->Clone( "H1D_jetEfficiency" );

    // Compute efficiency = matched / generated
    (*H1D_jetEfficiency)->Divide(*H1D_jetEfficiency, H1D_jetPt_mcp, 1., 1., "b");
  
}
//unused
/*
void GetJetEfficiency(TH1D** H1D_jetEfficiency, TH2D* Response_fine_Smeared, TH2D* Response_fine_Org, int iDataset, int iRadius) {
  cout << "############# GetJetEfficiency called! ###############" << endl;
  if (!Response_fine_Smeared) {
        cout << "GetJetEfficiency ERROR: Response_fine_Smeared is NULL!" << endl;
        *H1D_jetEfficiency = nullptr;
        return;
    }

  // --- Step 1: Project fine Smeared matrix onto rec axis (X) then rebin it as gen binning, that is the numerotr of the efficiency
  TH1D* projGen = Response_fine_Smeared->ProjectionY("projGen", 1, Response_fine_Smeared->GetNbinsY(), "e");
  if (!projGen) {
        cout << "######### GetJetEfficiency ERROR: projGen is NULL! ########" << endl;
        *H1D_jetEfficiency = nullptr;
        return;
    }
  else {
      cout << "projGen successfully created with " << projGen->GetNbinsX() << " bins." << endl;
  }

  TH1D* matchedGen = (TH1D*)projGen->Rebin(nBinPtJetsGen[iRadius], "matchedGen", ptBinsJetsGen[iRadius]); // numerator of the efficiency
  if (!matchedGen) {
        cout << "######### GetJetEfficiency ERROR: matchedGen is NULL! ########" << endl;
        *H1D_jetEfficiency = nullptr;
        return;
    }
  else {
      cout << "matchedGen successfully created with " << matchedGen->GetNbinsX() << " bins." << endl;
  }

  // --- Step 2: Get the mcp matched gen binning spectrum (smeared -> proj Y )

  // TH2D* respRebinned = nullptr;
  // TH1D* projGenSmeared = nullptr;
  // respRebinned = (TH2D*)Response_fine_Smeared->Clone("respRebinned");
  // MergeResponseMatrixBins(respRebinned, iDataset, iRadius, "noPriorMerging");
  // projGenSmeared = respRebinned->ProjectionY("projGenSmeared", 1, respRebinned->GetNbinsX(), "e");

  // // --- Step 3: Get the mcp matched gen binning spectrum (original -> proj Y )

  // TH2D* OrgRespRebinned = nullptr;
  // TH1D* projGenOrg = nullptr;
  // OrgRespRebinned = (TH2D*)Response_fine_Org->Clone("OrgRespRebinned");
  // MergeResponseMatrixBins(OrgRespRebinned, iDataset, iRadius, "noPriorMerging");
  // projGenOrg = OrgRespRebinned->ProjectionY("projGenOrg", 1, OrgRespRebinned->GetNbinsX(), "e");

  TH1D* projGenOrg = Response_fine_Org->ProjectionY("projGenOrg", 1, Response_fine_Org->GetNbinsY(), "e");
  TH1D* matchedGenOrig = (TH1D*)projGenOrg->Rebin(nBinPtJetsGen[iRadius], "matchedGenOrig", ptBinsJetsGen[iRadius]); // numerator of the efficiency

  // --- Step 4: Reweight the original mcp gen binning spectrum "mcpGenBin" with the ratio of (original/smeared) proj gen spectra, weight = projGenOrg/projGenSmeared 

  TH1D* mcpGenBin = nullptr;
  Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(mcpGenBin, iDataset, iRadius, false, ""); //original full mcp gen binning specatrum
  if (!mcpGenBin) {
        cout << "GetJetEfficiency ERROR: mcpGenBin is NULL!" << endl;
        *H1D_jetEfficiency = nullptr;
        return;
    }

  TH1D* mcpGen_new = (TH1D*)mcpGenBin->Clone("mcpGen_new");
  reweightHistG(mcpGen_new, projGenOrg, projGenSmeared); // mcpGen_new is the new denominator of the efficiency

  *H1D_jetEfficiency = (TH1D*)matchedGen->Clone("H1D_jetEfficiency");
  (*H1D_jetEfficiency)->Divide(*H1D_jetEfficiency, mcpGen_new, 1., 1., "b");

  cout << " ############# GetJetEfficiency done properly with GetJetEfficiency function #################" << endl;
  delete respRebinned;
  delete mcpGen_new;
}
  */
// ------------------------------------------------------------


