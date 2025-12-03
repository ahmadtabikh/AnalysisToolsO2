#include "TSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TDatime.h" 

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

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void unfoldSpec_statVariations_Ahmad(){

  int iDataset = 0;
  int iRadius = 1;
  int nRepeats = 100;
  bool doBayes = false;
  bool varySpec = true;
  bool varyResp = false;
  int errTreatment = 0;

  int nIterSpecBayes = 4;
  int nIterSpecSVD   = 7;
  int nIterSpec = doBayes ? nIterSpecBayes : nIterSpecSVD;

  // ##### Original Unfolding #####

  // --- Get Original Response Matrix (roounfold object + fine binned TH2D)
  RooUnfoldResponse* OrgResponse;
  TH2D* TH2D_Response_fine;
  bool smearResp = false;
  bool resetErr = false;
  GetResponse(&OrgResponse, &TH2D_Response_fine, smearResp, resetErr);

  // --- Get Measured Spectrum (rec binning, pre width scaling)
  TH1D* measuredInput;
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
  
  // --- Correct measured for Purity 
  TH1D* mcdRecBin;
  TH1D* purity;
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(mcdRecBin, iDataset, iRadius, options); 
  GetJetPurity(&purity, TH2D_Response_fine, mcdRecBin, iRadius);
  TH1D* measured = (TH1D*) measuredInput->Clone("measured_correctedForPurity");
  measured->Multiply(purity);

  // --- Unfold Priginal spectrum
  TH1D* TH1D_Unfolded_original;
  unfoldSpec(OrgResponse, measured, &TH1D_Unfolded_original, nIterSpec, doBayes, errTreatment);  

  // --- Kinematic efficiency
  TH1D* kinematicEfficiency;
  TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  Get_ResponseMatrix_Pt_KinematicEffiency(kinematicEfficiency, TH2D_Response_fine, name_H1D_kinematicEfficiency, iRadius);

  // --- Matched efficiency
  TH1D* jetEfficiency = nullptr;
  GetJetEfficiency(&jetEfficiency, TH2D_Response_fine, iDataset, iRadius);

  // apply efficiencies to unfolded spectrum
  TH1D_Unfolded_original->Divide(jetEfficiency);
  TH1D_Unfolded_original->Divide(kinematicEfficiency);
  NormaliseRawHistToNEvents(TH1D_Unfolded_original, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
  TransformRawHistToYield(TH1D_Unfolded_original); // final unfolded spectrum

  // ---------------------------
  // smear + unfold

  TProfile* TProfile_JetPtVar = new TProfile("TProfile_JetPtVar","",nBinPtJetsRec[iRadius], ptBinsJetsRec[iRadius],"S");  
    
  for(int rep=0; rep < nRepeats; rep++){

    cout<<"unfold FF rep "<<rep<<endl;

    TH1D* TH1D_Unfolded; 
    
    if(varyResp){
      RooUnfoldResponse* SmearedResponse;
      TH2D* Response_fine_Smeared;
      bool smearResp = true;
      GetResponse(&SmearedResponse, &Response_fine_Smeared, smearResp, resetErr);
    }
    
    if(varySpec){
      TH1D* measuredSmeared = (TH1D*) measuredInput->Clone("measuredSmeared");
      measuredSmeared->Reset();
      smearTH1(measuredInput,measuredSmeared);
      TH1D* purityS;
      GetJetPurity(&purityS, (varyResp ? Response_fine_Smeared : TH2D_Response_fine), mcdRecBin, iRadius);
      measuredSmeared->Multiply(purityS);
    }

    if(varyResp) unfoldSpec(SmearedResponse, measuredSmeared,  &TH1D_Unfolded, nIterSpec, doBayes, errTreatment);
    else         unfoldSpec(OrgResponse,            measured,  &TH1D_Unfolded, nIterSpec, doBayes, errTreatment);    
    
    // ---------------------------
    // Kinematic efficiency from smeared response matrix
    TH1D* kinematicEfficiencyS;
    TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+ "_rep" + Form("%d", rep);
    Get_ResponseMatrix_Pt_KinematicEffiency(kinematicEfficiencyS, (varyResp ? Response_fine_Smeared : TH2D_Response_fine) , name_H1D_kinematicEfficiency, iRadius);

    // matched efficiency
    TH1D* jetEfficiency = nullptr;
    GetJetEfficiency(&jetEfficiency, (varyResp ? Response_fine_Smeared : TH2D_Response_fine) , iDataset, iRadius);

    // apply efficiencies to unfolded spectrum
    TH1D_Unfolded->Divide(jetEfficiency);
    TH1D_Unfolded->Divide(kinematicEfficiencyS);
    NormaliseRawHistToNEvents(TH1D_Unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    TransformRawHistToYield(TH1D_Unfolded); // final smeared unfolded spectrum


    FillVariations(TH1D_Unfolded, TProfile_JetPtVar); 

    delete TH1D_Unfolded;
    delete Response_fine_Smeared;    
  }

  // --------------------------------------

  // unfolded spec errors as assigned by RooUnfold
  TH1D* TH1D_Unfolded_original_error = (TH1D*) TH1D_Unfolded_original->Clone("TH1D_Unfolded_original_error");
  TH1D_Unfolded_original_error->Reset();
  
  for(int bin=1; bin<=TH1D_Unfolded_original->GetNbinsX(); bin++){ 
    if(TH1D_Unfolded_original->GetBinContent(bin)){
      TH1D_Unfolded_original_error->SetBinContent(bin,TH1D_Unfolded_original->GetBinError(bin)/TH1D_Unfolded_original->GetBinContent(bin));
      TH1D_Unfolded_original_error->SetBinError(bin,0);
    }
  }

  // poissonian errors of unfolded spec bin content 
  TH1D* TH1D_Unfolded_original_PoissErr = (TH1D*) TH1D_Unfolded_original->Clone("TH1D_Unfolded_original_PoissErr");
  TH1D_Unfolded_original_PoissErr->Reset();
  
  for(int bin=1; bin<=TH1D_Unfolded_original->GetNbinsX(); bin++){ 
    if(TH1D_Unfolded_original->GetBinContent(bin)){
      double cont = TH1D_Unfolded_original->GetBinContent(bin);
      double poissRelErr = 0;
      if(cont > 0){
	      poissRelErr = 1/TMath::Sqrt(cont);
      }
      TH1D_Unfolded_original_PoissErr->SetBinContent(bin,poissRelErr);
      TH1D_Unfolded_original_PoissErr->SetBinError(bin,0);
    }
  }

  // // error from covariance
  // TH1F* hCovDiagErr = (TH1F*) hSpecUnfolded->Clone("hCovDiagErr");
  // hCovDiagErr->Reset();

  // for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
  //   Double_t cont = hSpecUnfolded->GetBinContent(bin);
  //   if(cont){
  //     Double_t err = TMath::Sqrt(h2Cov->GetBinContent(bin,bin));
  //     hCovDiagErr->SetBinContent(bin,err/cont);
  //     hCovDiagErr->SetBinError(bin,0);
  //   }
  // }

  
  // error from variations 
  TH1D* hRelUncert = (TH1D*) TH1D_Unfolded_original->Clone("hRelUncert");
  hRelUncert->Reset();

  int nBins = hRelUncert->GetNbinsX();

  for (int bin = 1; bin <= nBins; bin++) {
      double mean = TProfile_JetPtVar->GetBinContent(bin);
      double rms  = TProfile_JetPtVar->GetBinError(bin);  // because of option "S"

      if (mean > 0) {
          hRelUncert->SetBinContent(bin, rms / mean);
      } else {
          hRelUncert->SetBinContent(bin, 0);
      }

      hRelUncert->SetBinError(bin, 0);
  }

  
  // ---------------------------------------
  // plot 

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",610,470);
  c1->Divide(2,1);
   
  c1->cd(1);
  gPad->SetLogy();
  hSpecUnfolded->SetTitle("unfolded jet spectrum");
  hSpecUnfolded->SetLineColor(2);
  hSpecUnfolded->GetXaxis()->SetRangeUser(0,100);
  hSpecUnfolded->SetXTitle("p_{T}^{jet}"); 
  hSpecUnfolded->SetYTitle("1/nJets dN/dp_{T}"); 
  hSpecUnfolded->DrawCopy();
  
  // --

  c1->cd(2);
  gPad->SetLogy();
  hJetPtVar->SetTitle(Form("unfolded pseudodata, %d variations, error = spread",nRepeats));
  hJetPtVar->SetLineColor(2);
  hJetPtVar->GetXaxis()->SetRangeUser(0,100);
  hJetPtVar->SetXTitle("p_{T}^{jet}"); 
  hJetPtVar->SetYTitle("1/nJets dN/dp_{T}"); 
  hJetPtVar->DrawCopy();
  
  // --
  
  TCanvas *c2 = new TCanvas("c2","",460,560);
  c2->Divide(1,1);
   
  c2->cd(1);
  TString strTit;
  if(doBayes) strTit = "relative error, Bayes unfolding";
  else        strTit = "relative error, SVD unfolding";
  hVarErr->SetTitle(strTit);
  hVarErr->SetLineColor(2);
  hVarErr->GetXaxis()->SetRangeUser(0,100);
  hVarErr->GetYaxis()->SetRangeUser(0,0.2);
  hVarErr->SetMarkerStyle(20);
  hVarErr->SetMarkerColor(4);
  hVarErr->SetXTitle("p_{T}^{jet} (GeV/c)"); 
  hVarErr->SetYTitle("relative error");
  hVarErr->GetYaxis()->SetTitleOffset(1.2);
  hVarErr->Draw("PM");

  hSpecUnfoldedPoissErr->SetMarkerStyle(24);
  hSpecUnfoldedPoissErr->SetMarkerColor(4);
  hSpecUnfoldedPoissErr->Draw("PM same");

  
  hSpecUnfoldedErr->SetMarkerStyle(20);
  hSpecUnfoldedErr->SetMarkerColor(2);
  hSpecUnfoldedErr->Draw("PM same");
   
  TLegend* leg1 = new TLegend(0.15,0.70,0.61,0.87);
  leg1->SetTextSize(0.02);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(1);
  if(errTreatment == -1)     leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kCovariance, from cov. matrix","P");
  else if(errTreatment == 1) leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kErrors, from diag. elements of cov. matrix","P");
  else if(errTreatment == 2) leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kNoError","P");
  else if(errTreatment == 3) leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kCovToy","P");
  else 	 leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment default ","P");
  leg1->AddEntry(hVarErr,"unfolded pseudodata: spread","P");
  leg1->AddEntry(hSpecUnfoldedPoissErr,"Poisson error unfolded spectrum","P");
  leg1->Draw("");


  TCanvas *c3 = new TCanvas("c3","",460,560);
  c3->Divide(1,1);
   
  c3->cd(1);
  TString strTit2;
  if(doBayes) strTit2 = "relative error cov matrix, Bayes unfolding";
  else        strTit2 = "relative error cov matrix, SVD unfolding";
  hCovDiagErr->SetTitle(strTit2);
  hCovDiagErr->SetLineColor(2);
  hCovDiagErr->GetXaxis()->SetRangeUser(0,100);
  hCovDiagErr->GetYaxis()->SetRangeUser(0,0.2);
  hCovDiagErr->SetMarkerStyle(20);
  hCovDiagErr->SetMarkerColor(4);
  hCovDiagErr->SetXTitle("p_{T}^{jet} (GeV/c)"); 
  hCovDiagErr->SetYTitle("relative error");
  hCovDiagErr->GetYaxis()->SetTitleOffset(1.2);
  hCovDiagErr->Draw("PM");

   
  // --

  TString strApp;
  if(doBayes) strApp.Form("Bayes_statVariations_nRep%d",nRepeats);
  else        strApp.Form("SVD_statVariations_nRep%d",nRepeats);
  if(errTreatment == -1) strApp += "_errkCov";
  if(errTreatment ==  1) strApp += "_errkErr";
  if(errTreatment ==  2) strApp += "_errkNoError";
  if(errTreatment ==  3) strApp += "_errkCovToy";
  if(varySpec) strApp += "_varySpec";
  if(varyResp) strApp += "_varyResp";
  
  // c1->SaveAs(Form("spectra_%s.pdf",strApp.Data()));
  c2->SaveAs(Form("errors_%s.pdf",strApp.Data()));

  // ---- 

  delete hSpecUnfolded;
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
////////////////////////////////////////////////////////////////////////////////////// functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ------------------------------------------------------------
void unfoldSpec(RooUnfoldResponse* response, TH1D* measured, TH1D** hist_unfold, int unfoldParameterInput, bool doBayes = kTRUE, int errTreatment = 0){
  
  if (!response || !measured || !hist_unfold || !h2Cov) {
      cerr << "unfoldSpec: null pointer input." << endl;
      return;
  }

  // --- prepare a RooUnfold object (base pointer) ---
  RooUnfold* unfold = nullptr;
  RooUnfoldBayes* unfoldBayes = nullptr;
  RooUnfoldSvd*   unfoldSvd    = nullptr;

  // unfold spectrum 
  RooUnfoldBayes* unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);

  RooUnfold* unfold = unfoldBayes; // default Bayes

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
        Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdK, texSvdDvector, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");
      }
    }
    // Now create the real SVD object with the chosen regularisation parameter
    unfoldSvd = new RooUnfoldSvd(response, measured, unfoldParameterInput);
    unfold = unfoldSvd;
    cout<<" SVD unfolding, nIter "<<unfoldParameterInput<<endl;
    
    delete unfoldSvdInit;
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
      cerr << "unfoldSpec: Hreco returned nullptr." << endl;
      delete unfold; // delete concrete RooUnfold
      return;
  }

  *hist_unfold = (TH1D*) hReco->Clone();  // Important: clone it!

  delete unfold;
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
                errSmeared = err;
            }

            hSmeared->SetBinContent(iX, iY, contSmeared);
            hSmeared->SetBinError(iX, iY, err);
        }
    }
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

        cout << "hist: " << hist->GetName() 
             << " bin " << binx 
             << " x = " << hist->GetXaxis()->GetBinCenter(binx)
             << " cont = " << cont 
             << " err = " << err 
             << " contSmeared = " << contSmeared 
             << " errCorr = " << errCorr << endl;

        hSmeared->SetBinContent(binx, contSmeared);
        hSmeared->SetBinError(binx, errCorr);
    }
}
// ------------------------------------------------------------

void GetResponse(RooUnfoldResponse** respJetPt, TH2D** Response_fine_Smeared_out, bool smearResp, bool resetErr = kFALSE){
  
  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;                           // detector response as proba fine binning
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;// det + fluct response as proba fine binning

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");
  ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, "mcpPriorUnfolding"); // weight the fine proba combined matrix with the prior (mcp) = fine event matrix 

  TH2D* Response_fine_Smeared = nullptr;
  TH2D* Response_rebinned = nullptr;

  if(smearResp){ //Clones the response matrix, applies smearResponse(), then deletes the temporary.
    // --- smear response 
    Response_fine_Smeared = (TH2D*) H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("Response_fine_Smeared");
    Response_fine_Smeared->Sumw2();
    smearResponseTH2D(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, Response_fine_Smeared); // org, smeared (fine binning)
    
    Response_rebinned = (TH2D*) Response_fine_Smeared->Clone("Response_smeared_rebinned");
    MergeResponseMatrixBins(Response_rebinned, iDataset, iRadius, options); // rebin the smeared response matrix with no bin area scaling
  }
  else{
    Response_fine_Smeared = (TH2D*) H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("Response_org_rebinned");
    Response_rebinned = (TH2D*)Response_fine_Smeared->Clone("Response_org_rebinned");
    MergeResponseMatrixBins(Response_rebinned, iDataset, iRadius, options); // rebin the smeared response matrix with no bin area scaling
  }

  if(resetErr) resetErrorsTH2D(Response_smeared_rebinned);
  
  RooUnfoldResponse* response = new RooUnfoldResponse(0, 0, Response_rebinned); // create the roounfold matrix object

  *respJetPt = response;
  *Response_fine_Smeared_out = Response_fine_Smeared;


  delete H2D_jetPtResponseMatrix_fluctuations;
  delete H2D_jetPtResponseMatrix_detectorResponse;
  delete H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  
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
}
// ------------------------------------------------------------

void GetJetPurity(TH1D** H1D_jetPurity, TH2D* Response_fine_Smeared, TH1D* mcdRecBin, int iRadius){
    // --- Step 1: Project fine matrix onto reconstructed axis
    TH1D* H1D_jetPt_smeared_projRec_RecBin = (TH1D*) Response_fine_Smeared->ProjectionX("jetPt_mcdMatched_smeared_forRecBin", 1, Response_fine_Smeared->GetNbinsY(), "e");

    // --- Step 2: Rebin the projected histogram
    TH1D* H1D_jetPt_mcdMatched_recBinning = (TH1D*) H1D_jetPt_smeared_projRec_RecBin->Rebin(nBinPtJetsRec[iRadius], "jetPt_mcdMatched_recBinning_rebinned", ptBinsJetsRec[iRadius] );

    // --- Step 3: Clone to create the purity histogram
    *H1D_jetPurity = (TH1D*) H1D_jetPt_mcdMatched_recBinning->Clone("H1D_jetPurity");

    // --- Step 4: Compute purity = matched / measured (binomial errors)
    (*H1D_jetPurity)->Divide(*H1D_jetPurity, mcdRecBin, 1., 1., "b");

    // --- Clean up temporary histograms
    delete H1D_jetPt_smeared_projRec_RecBin;
    delete H1D_jetPt_mcdMatched_recBinning;
}
// ------------------------------------------------------------

void GetJetEfficiency(TH1D** H1D_jetEfficiency, TH2D* Response_fine_Smeared, int iDataset, int iRadius) {
    TH1D* projRec = (TH1D*)Response_fine_Smeared->ProjectionX("projRec", 1, Response_fine_Smeared->GetNbinsY(), "e");
    TH1D* matchedGen = (TH1D*)projRec->Rebin(nBinPtJetsGen[iRadius], "matchedGen", ptBinsJetsGen[iRadius]);
    TH2D* respRebinned = (TH2D*)Response_fine_Smeared->Clone("respRebinned");
    MergeResponseMatrixBins(respRebinned, iDataset, iRadius, "");
    TH1D* projGen = (TH1D*)respRebinned->ProjectionY("projGen", 1, respRebinned->GetNbinsX(), "e");

    TH1D* mcpGenBin = nullptr;
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(mcpGenBin, iDataset, iRadius, false, "");
    TH1D* mcpGen_new = (TH1D*)mcpGenBin->Clone("mcpGen_new");
    reweightHistG(mcpGen_new, projGen, projGen);

    *H1D_jetEfficiency = (TH1D*)matchedGen->Clone("H1D_jetEfficiency");
    (*H1D_jetEfficiency)->Divide(*H1D_jetEfficiency, mcpGen_new, 1., 1., "b");

    delete projRec; delete matchedGen; delete respRebinned; delete projGen; delete mcpGenBin; delete mcpGen_new;
}
// ------------------------------------------------------------
