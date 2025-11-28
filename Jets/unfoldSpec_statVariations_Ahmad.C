#if !defined( __CINT__) || defined(__MAKECINT__)

#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include <iostream>
#include "TPaveText.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH3F.h"
#include "THn.h"
#include  "TDatime.h" 

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h" 
#include "RooUnfoldSvd.h" 
#include "TSVDUnfold.h"
#include "RooUnfoldBinByBin.h" 

#include "UnfoldHelper.h"

#endif

using namespace std;



UnfoldHelper *UFH; 


// -------------------------------------------------

void unfoldSpec(RooUnfoldResponse* response, 
	 TH1D* measured, TH1D** hist_unfold, TH2D** h2Cov, Int_t unfoldParameterInput, Bool_t doBayes = kTRUE, Int_t errTreatment = 0){
  
  if (!response || !measured || !hist_unfold || !h2Cov) {
      cerr << "unfoldSpec: null pointer input." << endl;
      return;
  }

  // ----------------
  // --- prepare a RooUnfold object (base pointer) ---
    RooUnfold* unfold = nullptr;
    RooUnfoldBayes* unfoldBayes = nullptr;
    RooUnfoldSvd*   unfoldSvd    = nullptr;

  // unfold spectrum 
  RooUnfoldBayes* unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);

  RooUnfold* unfold = unfoldBayes; // default Bayes
  TH1D* hist_unfold;

  if (options.find("Svd") != std::string::npos) {
    int unfoldParameterSvdInitial = 1;
    RooUnfoldSvd* unfoldSvdInit = new RooUnfoldSvd(response, measured, unfoldParameterSvdInitial); 
    unfoldSvdInit->Hreco(); // necessary to have GetD() give a meaningful output
    TSVDUnfold *tsvdUnfold = (TSVDUnfold*)unfoldSvdInit->Impl();

    // Optionally draw D distribution if requested (kept from your original)
    if (tsvdUnfold) {
      TH1D* H1D_D = tsvdUnfold->GetD(); // may be nullptr if not available
      if (H1D_D) {
        TString inputUnfoldingName = (options.find("inputIsMCPFoldedTest") != std::string::npos) ? "_mcpFoldedTestInput" : "";
        TString* pdfName_regparam = new TString("Svd_regularisationd_distribution_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+(TString)mergingPrior+"_"+(TString)unfoldingPrior+inputUnfoldingName);
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
  } else if (options.find("Bayes") != std::string::npos) {
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
  TH1D* hReco = (TH1D*) unfold->Hreco(et);
  if (!hReco) {
      cerr << "unfoldSpec: Hreco returned nullptr." << endl;
      delete unfold; // delete concrete RooUnfold
      return;
  }

  *hist_unfold = hReco;
  
  delete unfold;
}

// ------------------------------------------------------------------------
//This function takes an original 2D THnSparse histogram (hnResponseOrg) and produces a “smeared” version (hnResponseSmeared) by randomly fluctuating each bin content according to its statistical error.
void smearResponseTH2D(TH2D* hOrig, TH2D* hSmeared) {
    int nBinsX = hOrig->GetNbinsX();
    int nBinsY = hOrig->GetNbinsY();

    cout << "smearResponseTH2D, nBinsX: " << nBinsX << endl;
    cout << "smearResponseTH2D, nBinsY: " << nBinsY << endl;

    for (int iY = 1; iY <= nBinsY; iY++) {
        if (!(iY % 100)) cout << "iY: " << iY << endl;

        for (int iX = 1; iX <= nBinsX; iX++) {
            double resp = hOrig->GetBinContent(iX, iY);
            double err = hOrig->GetBinError(iX, iY);
            double relErr = (resp != 0) ? err / resp : 0.0;

            // Gaussian smearing
            double contSmeared = gRandom->Gaus(resp, err);
            double errSmeared = relErr * contSmeared;

            if (contSmeared < 0) { // keep original if negative
                contSmeared = resp;
                errSmeared = err;
            }

            hSmeared->SetBinContent(iX, iY, contSmeared);
            hSmeared->SetBinError(iX, iY, errSmeared);
        }
    }
}

// ------------------------------------------------------------------------
// resetErrors sets all bin errors of a 2D THnSparse histogram to zero
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

// -------------------------------------------------------------
// This function reweights a histogram histG using the ratio of smeared to original histograms (histSmear / histOrg).
void reweightHistG(TH1D* histG, TH1D* histOrg, TH1D* histSmear) {
    Int_t nBins = histG->GetXaxis()->GetNbins();

    for (Int_t bin = 0; bin <= nBins; bin++) { // skip underflow/overflow
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
            errCorr = err * contSmeared / cont;
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
// ---------------------------------------------------
// Purpose of the function
// Reads response histograms from a ROOT file.
// Optionally smears the response to simulate statistical fluctuations.
// Optionally resets bin errors to zero.
// Builds a RooUnfoldResponse object, which can be used for unfolding.
void GetResponse(RooUnfoldResponse** respJetPt, Int_t nBinsX, Double_t* binsX, Bool_t smearResp, Bool_t resetErr = kFALSE){
  
  TH1D* mcp;
  TH1D* mcd;
  TH1D* mcpGenBin;
  TH1D* mcdRecBin;  
  TH1D* H1D_jetPt_projGen;
  TH1D* H1D_jetPt_smeared_projRec;
  TH1D* H1D_jetPt_smeared_projGen;
  TH1D* h1RespJetPt_smeared_projRec  = 0;
  TH1D* h1RespJetPt_smeared_projGen  = 0;
  TH1D* h1RespJetPtHistG_rew = 0;
  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;                           // detector response as proba fine binning
  TH2D* H2D_jetPtResp_Rebinned;                                             // det + fluct response as evnt large binning
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;// det + fluct response as proba fine binning

  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(mcp, iDataset, iRadius, false, options); 
  Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(mcd, iDataset, iRadius, false, options);

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  
  // ---------------------  
  // rebin before smearing

  Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(mcpGenBin, iDataset, iRadius, false, options); // h1RespJetPtHistG
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(mcdRecBin, iDataset, iRadius, options); // h1RespJetPtHistM
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResp_Rebinned, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // hnRespJetPt

  // -------------------------------
  Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_projGen, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) // mcp matched gen pt spectrum //h1RespJetPt_unsmeared_projGen

  if(smearResp){ //Clones the response matrix, applies smearResponse(), then deletes the temporary.

    // smear response 
    TH2D* Response_temp = (TH2D*) H2D_jetPtResp_Rebinned->Clone("Response_temp");
    smearResponseTH2D(Response_temp, H2D_jetPtResp_Rebinned); // org, smeared : overwrite the original Matrix with the smeared one
    delete Response_temp;

    // project smeared response on rec axis to replace HistM (avoid spurious 'fakes' correction)
    H1D_jetPt_smeared_projRec = (TH1D*)H2D_jetPtResp_Rebinned->ProjectionX("jetPt_mcdMatched_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtResp_Rebinned->GetNbinsY(), "e");

    // project smeared response on gen axis to replace HistG (avoid spurious 'efficiency' correction)
    H1D_jetPt_smeared_projGen = (TH1D*)H2D_jetPtResp_Rebinned->ProjectionY("jetPt_mcpMatched_GenBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtResp_Rebinned->GetNbinsX(), "e");
    
    mcpGenBin_new = (TH1D*) mcpGenBin->Clone("mcpGenBin_new");    
    reweightHistG(mcpGenBin_new, H1D_jetPt_projGen , H1D_jetPt_smeared_projGen);
    //Rec axis projection → replaces HistM for unfolding. Gen axis projection → used to reweight HistG to match smeared response.
  }

  if(resetErr) resetErrorsTH2D(H2D_jetPtResp_Rebinned);
  
  RooUnfoldResponse* response;
  response = new RooUnfoldResponse(0, 0, H2D_jetPtResp_Rebinned);

  *respJetPt = response;

  delete mcp;
  delete mcd;
  delete mcpGenBin;
  delete mcdRecBin;   

  delete H2D_jetPtResponseMatrix_fluctuations;
  delete H2D_jetPtResponseMatrix_detectorResponse;
  delete H2D_jetPtResp_Rebinned;
  delete H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  delete h1RespJetPtHistG;
  delete h1RespJetPtHistM;

  delete H1D_jetPt_smeared_projRec;
  delete H1D_jetPt_smeared_projGen;
  delete mcpGenBin_new;
}
// ------------------------------------------------------------


// -------------------------------------------
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
        histP->Fill(cent, cont);
    }
}

// ---------------------------------------------------

void unfoldSpec_statVariations(Int_t nRepeats = 100, Bool_t doBayes = kTRUE, Bool_t varySpec = kTRUE, Bool_t varyResp = kFALSE, Int_t errTreatment = 0){
  
  // 1D unfolding of spec 
  // test RooUnfold errors 
  // varyFFResp: smear response matrix
  // varyFFSpec: smear input spectrum histo
  // for pure spectrum variations 100 repeats take ~ 2-3 seconds, but response matrix variations very slow (leaky !) 
  
  UFH = new UnfoldHelper(fnameResp);
  
  Double_t binsX[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,38,44,50,58,66,76,86,100,120,150,200};
  Int_t nBinsX     = 29;

  Int_t nIterSpecBayes = 4;
  Int_t nIterSpecSVD   = 11;

  Int_t nIterSpec = doBayes ? nIterSpecBayes : nIterSpecSVD;
  
  // ------------------------------

  TFile f2(strInFileData);

  TList* list = getList(strID_data);

  TH1D* h1JetPtRec_org = (TH1D*) list->FindObject("hJetSpecIncRec");
  
  h1JetPtRec_org->SetDirectory(0);

  f2.Close();
 
  // ---------------------------
  // bias corr

  Bool_t biasCorr = kTRUE;

  if(biasCorr){

    TH1D* fh1JetPtRecMatchBC;
    TH1D* fh1JetPtRecIncBC;

    TH2D* fh2FFZRecMatchBC;
    TH2D* fh2FFZRecIncBC;

    TH1D* hRatioJetPtRecMatch  = UFH->ratioFFJetPtMatchInc(strIDResp, strIDFFRec, &fh1JetPtRecMatchBC, &fh1JetPtRecIncBC, kTRUE);

    h1JetPtRec_org->Multiply(hRatioJetPtRecMatch);
  }

  // ---------------------------
  // rebin
  
  TH1D* h1JetPtRec = (TH1D*) h1JetPtRec_org->Rebin(nBinsX,"h1JetPtRec_reb",binsX);

  // ---------------------------
  // response matrix 

  RooUnfoldResponse* respJetPt = 0;

  Bool_t smearResp = kFALSE;
  Bool_t resetErr  = kFALSE; // don't remember why errors were reset ? (errTreatment == 0) ? kFALSE : kTRUE;

  GetResponse(&respJetPt, nBinsX, binsX, smearResp, resetErr);

  TH1D* hSpecUnfolded  = 0;
  TH2D* h2Cov = 0;
  
  unfoldSpec(respJetPt, h1JetPtRec, &hSpecUnfolded, &h2Cov, nIterSpec, doBayes, errTreatment);

  //return; // TEST!!!
 
  // ---------------------------
  // smear + unfold

  TProfile* tpJetPtVar = new TProfile("tpJetPtVar","",nBinsX,binsX,"S");  
    
  for(Int_t rep=0; rep<nRepeats; rep++){

    cout<<"unfold FF rep "<<rep<<endl;

    TH1D* h1JetPtRec_rep =  new TH1D(*h1JetPtRec);
    h1JetPtRec_rep->SetName(Form("h1JetPtRec_rep%d",rep));
    
    RooUnfoldResponse* respJetPt_rep = 0;
    
    TH1D* hUnfolded_rep = 0; 
    TH2D* h2Cov         = 0; 
    
    if(varyResp){
      Bool_t smearResp = kTRUE;
      GetResponse(&respJetPt_rep, nBinsX, binsX, smearResp);
    }

    // delete respJetPt_rep; // TEST !!!
    // continue; // TEST !!!
    
    if(varySpec){
      h1JetPtRec_rep->Reset();
      smearTH1(h1JetPtRec,h1JetPtRec_rep);
    }

    if(varyResp) unfoldSpec(respJetPt_rep, h1JetPtRec_rep,  &hUnfolded_rep, &h2Cov, nIterSpec, doBayes, errTreatment);
    else         unfoldSpec(respJetPt,     h1JetPtRec_rep,  &hUnfolded_rep, &h2Cov, nIterSpec, doBayes, errTreatment);    
    
    FillVariations(hUnfolded_rep,tpJetPtVar); 

    cout<<" rep "<<rep<<" hUnfolded_rep bin Cont 5 "<<hUnfolded_rep->GetBinContent(5)<<endl;
    
    delete hUnfolded_rep;
    delete respJetPt_rep;
    delete h1JetPtRec_rep;		    
  }

  // --------------------------------------

  // unfolded spec errors as assigned by RooUnfold
  TH1F* hSpecUnfoldedErr = (TH1F*) hSpecUnfolded->Clone("hSpecUnfoldedErr");
  hSpecUnfoldedErr->Reset();
  
  for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
    if(hSpecUnfolded->GetBinContent(bin)){
      hSpecUnfoldedErr->SetBinContent(bin,hSpecUnfolded->GetBinError(bin)/hSpecUnfolded->GetBinContent(bin));
      hSpecUnfoldedErr->SetBinError(bin,0);
    }
  }

  // poissonian errors of unfolded spec bin content 
  TH1F* hSpecUnfoldedPoissErr = (TH1F*) hSpecUnfolded->Clone("hSpecUnfoldedPoissErr");
  hSpecUnfoldedPoissErr->Reset();
  
  for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
    if(hSpecUnfolded->GetBinContent(bin)){
      Double_t cont = hSpecUnfolded->GetBinContent(bin);
      Double_t poissRelErr = 0;
      if(cont > 0){
	poissRelErr = 1/TMath::Sqrt(cont);
      }
      hSpecUnfoldedPoissErr->SetBinContent(bin,poissRelErr);
      hSpecUnfoldedPoissErr->SetBinError(bin,0);
    }
  }

  // poissonian errors of raw spec bin content 
  TH1F* hSpecRawErr = (TH1F*) h1JetPtRec->Clone("hSpecRawErr");
  hSpecRawErr->Reset();
  
  for(Int_t bin=1; bin<=h1JetPtRec->GetNbinsX(); bin++){ 
    if(h1JetPtRec->GetBinContent(bin)){
      hSpecRawErr->SetBinContent(bin,h1JetPtRec->GetBinError(bin)/h1JetPtRec->GetBinContent(bin));
      hSpecRawErr->SetBinError(bin,0);
    }
  }

  // error from covariance
  TH1F* hCovDiagErr = (TH1F*) hSpecUnfolded->Clone("hCovDiagErr");
  hCovDiagErr->Reset();

  for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
    Double_t cont = hSpecUnfolded->GetBinContent(bin);
    if(cont){
      Double_t err = TMath::Sqrt(h2Cov->GetBinContent(bin,bin));
      hCovDiagErr->SetBinContent(bin,err/cont);
      hCovDiagErr->SetBinError(bin,0);
    }
  }

  
  // error from variations 
  TH1F* hVarErr = (TH1F*) hSpecUnfolded->Clone("hVarErr");
  hVarErr->Reset();
  
  for(Int_t bin=1; bin<=hVarErr->GetNbinsX(); bin++){ 

    //     //cout<<" jbin "<<jbin<<" bin "<<bin<<" hVar err "<<hVar->GetBinError(bin)<<endl;
    
    if(tpJetPtVar->GetBinContent(bin)){
      hVarErr->SetBinContent(bin,tpJetPtVar->GetBinError(bin)/tpJetPtVar->GetBinContent(bin));
      hVarErr->SetBinError(bin,0);
    }
  }

  // -----------------
  // scale by binwidth
  hSpecUnfolded->Scale(1,"width");
  TH1D* hJetPtVar = tpJetPtVar->ProjectionX("hJetPtVar"); // option 'width' doesn't work for TProfile, so first project to TH1
  hJetPtVar->Scale(1,"width");
  
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

