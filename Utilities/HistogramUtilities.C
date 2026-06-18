#ifndef HISTOGRAM_UTILITIES_C
#define HISTOGRAM_UTILITIES_C


#include "HistogramUtilities.h"
#include "HistogramPlotting.h"
#include "HistogramPlotting.C"
#include "../Settings/GlobalSettings.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include<array>





std::vector<double> GetTH1Bins(TH1* H1_histo) {
  std::vector<double> bins;
  for(int iBin = 1; iBin <= H1_histo->GetNbinsX(); iBin++){
    bins.push_back(H1_histo->GetBinLowEdge(iBin));
  }
  bins.push_back(H1_histo->GetXaxis()->GetBinUpEdge(H1_histo->GetNbinsX()));
  return bins;
}

std::vector<double> MakeVariableBinning_twoWidths(double xMin, int nLeft, double xMiddle, double xMax, int nRight) {
  if (xMin > xMiddle || xMiddle > xMax){
    cout << "ERROR in MakeVariableBinning_twoWidths(): it should be xMin < xMiddle < xMax but it is not the case" << endl;
  }
  std::vector<double> bins;
  double binWidth;
  binWidth = (xMiddle - xMin) / nLeft;
  for(int i = 0; i < nLeft; i++){
    bins.push_back(xMin + i * binWidth);
  }
  binWidth = (xMax - xMiddle) / nRight;
  for(int i = 0; i < nRight+1; i++){ // using nRight +1 here so that the rightmost edge is added. But due to double prec, might not be exactly xMax; alternative could be to forego the +1 and manually pushback xMax
    bins.push_back(xMiddle + i * binWidth);
  }
  return bins;
}

std::vector<double> MakeConstantSizeBinning(double xMin, double xMax, int nBins) {
  if (xMin > xMax){
    cout << "ERROR in MakeConstantSizeBinning(): it should be xMin < xMax but it is not the case" << endl;
  }
  std::vector<double> bins;
  double binWidth;
  binWidth = (xMax - xMin) / nBins;
  for(int i = 0; i < nBins+1; i++){ // using nRight +1 here so that the rightmost edge is added. But due to double prec, might not be exactly xMax; alternative could be to forego the +1 and manually pushback xMax
    bins.push_back(xMin + i * binWidth);
  }
  return bins;
}

TH1D* TGraphErrorsToTH1D(TGraphErrors* tgraphToConvert) {
  // assumes errors in x axis correspond to binning
  int nBins = tgraphToConvert->GetN();
  std::vector<double> binningVector;
  for(int iPoint = 0; iPoint < nBins; iPoint++){
    binningVector.push_back(tgraphToConvert->GetPointX(iPoint) - tgraphToConvert->GetErrorX(iPoint));
  }
  binningVector.push_back(tgraphToConvert->GetPointX(nBins-1) + tgraphToConvert->GetErrorX(nBins-1)); // right edge of last bin, as only left edge of bins have been filled so far

  double* binning = &binningVector[0];
  TH1D* convertedTH1D = new TH1D("convertedTH1D", "convertedTH1D", nBins, binning);
  for(int iBin = 1; iBin <= nBins; iBin++){
    int iPoint = iBin - 1;
    convertedTH1D->SetBinContent(iBin, tgraphToConvert->GetPointY(iPoint));
    convertedTH1D->SetBinError(iBin, tgraphToConvert->GetErrorY(iPoint)); // for some reason GetErrorY gets the error in X
  }

  return convertedTH1D;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Operations //////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> MakeVariableBinning_logarithmic(TH1D* histogram, int nBinsRough) {
  double xMin = histogram->GetBinCenter(1); // center to avoid 0, which is where most histograms start at
  double xMax = histogram->GetXaxis()->GetXmax();
  std::vector<double> bins;
  double xTemp = xMin;
  double binLogWidth = (log10(xMax) - log10(xMin)) / nBinsRough;
  // int originalHistBin;
  cout << "xMax = " << xMax << endl;
  cout << "xMin = " << xMin << endl;
  cout << "i = " << 0 << endl;
  cout << "rough bin edge:   " << xTemp << endl;
  int originalHistBin = histogram->GetXaxis()->FindBin(xTemp);
  double xTempRefined = histogram->GetXaxis()->GetBinLowEdge(originalHistBin);
  cout << "refined bin edge:   " << xTempRefined << endl;
  bins.push_back(xTempRefined);
  for(int i = 0; i < nBinsRough; i++){
    if (i !=0 && xTempRefined != bins.back()) {
    cout << "i = " << i << endl;
    cout << "refined bin edge: " << xTempRefined << endl;
      bins.push_back(xTempRefined);
    }
    xTemp = xTemp * pow(10, binLogWidth); // at this point, this edge is not necessarily matching any of the initial histogram edges
    originalHistBin = histogram->GetXaxis()->FindBin(xTemp);
    xTempRefined = histogram->GetXaxis()->GetBinLowEdge(originalHistBin);
  }
  bins.push_back(xMax);
  return bins;
}

// looks like it doesn't change my result, so my rebinning is good
// what they do here is renormalise each gen slice to 1; which we do outside of this function
TH2D* RebinVariableBins2D_aliPhysics(TH2D* hRMFine, int nBinsX, int nBinsY, double* binsX, double* binsY, __attribute__ ((unused)) bool useFunctionWeight){ // AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin  in https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx#L769
  TH2D hRM("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);

  // Rebin matrix hRMFine to dimensions of hRM
  // function returns matrix in TH2D format (hRM2) with dimensions from hRM
  TH2D *hRM2 = (TH2D*)hRM.Clone("hRM2");
  hRM2->Reset();

  //First normalize lines of input
  const Int_t nbinsNorm = hRM2->GetNbinsY();

  TArrayF *normVector = new TArrayF(nbinsNorm);

  for(int iy=1; iy<=hRM2->GetNbinsY(); iy++) {
    Double_t yLow = hRM2->GetYaxis()->GetBinLowEdge(iy);
    Double_t yUp = hRM2->GetYaxis()->GetBinUpEdge(iy);
    Int_t binyLowFine = hRMFine->GetYaxis()->FindBin(yLow);
    Int_t binyUpFine = hRMFine->GetYaxis()->FindBin(yUp)-1;
    normVector->SetAt(hRMFine->Integral(1,hRMFine->GetXaxis()->GetNbins(),binyLowFine,binyUpFine),iy-1); // Y and X have been changed as I use transposed matrices compared to aliphysics 
  }

  Double_t content, oldcontent = 0.;
  Int_t ixNew = 0;
  Int_t iyNew = 0;
  Double_t xvalLo, xvalUp, yvalLo, yvalUp;
  Double_t xmin = hRM2->GetXaxis()->GetXmin();
  Double_t ymin = hRM2->GetYaxis()->GetXmin();
  Double_t xmax = hRM2->GetXaxis()->GetXmax();
  Double_t ymax = hRM2->GetYaxis()->GetXmax();
  for(int ix=1; ix<=hRMFine->GetXaxis()->GetNbins(); ix++) {
    xvalLo = hRMFine->GetXaxis()->GetBinLowEdge(ix);
    xvalUp = hRMFine->GetXaxis()->GetBinUpEdge(ix);
    if(xvalLo<xmin || xvalUp>xmax) continue;
    ixNew = hRM2->GetXaxis()->FindBin(hRMFine->GetXaxis()->GetBinCenter(ix));
    for(int iy=1; iy<=hRMFine->GetYaxis()->GetNbins(); iy++) {
      yvalLo = hRMFine->GetYaxis()->GetBinLowEdge(iy);
      yvalUp = hRMFine->GetYaxis()->GetBinUpEdge(iy);
      if(yvalLo<ymin || yvalUp>ymax) continue;
      content = hRMFine->GetBinContent(ix,iy);
      iyNew = hRM2->GetYaxis()->FindBin(hRMFine->GetYaxis()->GetBinCenter(iy));
      oldcontent = hRM2->GetBinContent(ixNew,iyNew);
      Double_t weight = 1.;
      if(normVector->At(iyNew-1)>0.) {
	      weight = 1./normVector->At(iyNew-1);
      }
      hRM2->SetBinContent(ixNew,iyNew,oldcontent+content*weight);
    }
  }
  if(normVector) delete normVector;
  return hRM2;
}


TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, bool debug){
  if (debug == true) {cout << "___________ begin RebinVariableBins2D, nBinsX = " << nBinsX << ", nBinsY = " << nBinsY << endl;}

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);

  int ibinX_low, ibinX_high, ibinY_low, ibinY_high;
  double xLowEdgeRebin, xHighEdgeRebin, yLowEdgeRebin, yHighEdgeRebin;
  double H2D_hist_content, H2D_hist_contentError;
  for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    if (debug == true) {cout << "iBinX = " << iBinX << endl;}
    xLowEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX);
    xHighEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1);
    xLowEdgeRebin < binsX[0] ? ibinX_low = 0 : ibinX_low = H2D_hist->GetXaxis()->FindBin(xLowEdgeRebin); // according to definition of getbinlowedge, this takes into account under/overflows thanks to FindBin() giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    xHighEdgeRebin > binsX[nBinsX] ? ibinX_high = H2D_hist->GetNbinsX()+1 : ibinX_high = H2D_hist->GetXaxis()->FindBin(xHighEdgeRebin) - 1;
    for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      yLowEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY);
      yHighEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1);
      yLowEdgeRebin < binsY[0] ? ibinY_low = 0 : ibinY_low = H2D_hist->GetYaxis()->FindBin(yLowEdgeRebin); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
      yHighEdgeRebin > binsY[nBinsY] ? ibinY_high =  H2D_hist->GetNbinsY()+1 : ibinY_high = H2D_hist->GetYaxis()->FindBin(yHighEdgeRebin) - 1;

      H2D_hist_content = H2D_hist->IntegralAndError(ibinX_low, ibinX_high, ibinY_low, ibinY_high, H2D_hist_contentError);
      H2D_hist_rebinned.SetBinContent(iBinX, iBinY, H2D_hist_content);
      H2D_hist_rebinned.SetBinError(iBinX, iBinY, H2D_hist_contentError);
      if (debug == true) {cout << "ibinX_low = " << ibinX_low << ", ibinX_high = " << ibinX_high << ", ibinY_low = " << ibinY_low << ", ibinY_high = " << ibinY_high << "         --------          H2D_hist_rebinned(" << iBinX <<", "<< iBinY <<") = " << H2D_hist_rebinned.GetBinContent(iBinX, iBinY) << endl;}
    }
  }

  H2D_hist_rebinned.ResetStats(); // setbincontent interacts weirdly with getentries(); resetstats makes it so that getentries gives the sum of bin contents correctly
  if (debug == true) {cout << "RebinVariableBins2D - what of the errors" << endl;}
  if (debug == true) {cout << "H2D_hist nbins X =" << H2D_hist->GetNbinsX() << endl;}
  if (debug == true) {cout << "H2D_hist_rebinned nbins X =" << H2D_hist_rebinned.GetNbinsX() << endl;}

  std::stringstream ss;
  ss << H2D_hist->GetName() << "_RebinVariableBins2D";
  TString histName((TString)ss.str());

  if (debug == true) {cout << "___________ end RebinVariableBins2D" << endl;}

  return H2D_hist_rebinned;
  // why am i using this histName and ss?
}
TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY){
  return RebinVariableBins2D(H2D_hist, nBinsX, nBinsY, binsX, binsY, false);
}

TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum, bool debug){
  // H1D_priorSpectrum is assumed to have same binning as H2D_hist
  if (debug == true) {cout << "___________ begin RebinVariableBins2D, nBinsX = " << nBinsX << ", nBinsY = " << nBinsY << endl;}

  TH2D H2D_hist_rebinned("H2D_hist_rebinned", "H2D_hist_rebinned", nBinsX, binsX, nBinsY, binsY);
  H2D_hist_rebinned.Sumw2();

  TH2D* H2D_hist_temp = (TH2D*)H2D_hist->Clone(H2D_hist->GetName()+(TString)"_temp"); // needs this as its content will be changed for the weighting, we do not want to modify H2D_hist

  // declaration of variables used when calculating content and error of histograms
  double H2D_hist_content;
  double H2D_hist_contentError = 0;
  double H2D_hist_contentErrorA = 0;
  double H2D_hist_contentErrorB = 0;


  // prior weighting
  double priorWeightedSpectrumContent = 0;
  double priorWeightedSpectrumError = 0;
  for(int iBinXFine = 0; iBinXFine <= H2D_hist_temp->GetNbinsX()+1; iBinXFine++){ // 0 and n+1 take underflow and overflow into account
    priorWeightedSpectrumContent = H1D_priorSpectrum->GetBinContent(iBinXFine);
    priorWeightedSpectrumError = H1D_priorSpectrum->GetBinError(iBinXFine);
    for(int iBinYFine = 0; iBinYFine <= H2D_hist_temp->GetNbinsY()+1; iBinYFine++){ // 0 and n+1 take underflow and overflow into account
      H2D_hist_content = H2D_hist_temp->GetBinContent(iBinXFine, iBinYFine);
      H2D_hist_contentError = H2D_hist_temp->GetBinError(iBinXFine, iBinYFine);
      H2D_hist_temp->SetBinContent(iBinXFine, iBinYFine, H2D_hist_content * priorWeightedSpectrumContent);
      H2D_hist_temp->SetBinError(iBinXFine, iBinYFine, sqrt( priorWeightedSpectrumContent*priorWeightedSpectrumContent * H2D_hist_contentError*H2D_hist_contentError + H2D_hist_content*H2D_hist_content * priorWeightedSpectrumError*priorWeightedSpectrumError));
    }
  }

  cout << "##############################################" << endl;
  cout << "######### reweighting with prior #############" << endl;
  cout << "##############################################" << endl;

  int ibinX_low, ibinX_high, ibinY_low, ibinY_high;
  double xLowEdgeRebin, xHighEdgeRebin, yLowEdgeRebin, yHighEdgeRebin;
  double priorBinNorm, priorBinNormError;
  double hist2D_rebinned_content;
  double hist2D_rebinned_error;
  cout << "to think about: do I really give the value 0 if denominator is 0 ? when normalising by 1./priorBinNorm" << endl;
  for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // overflow and underflow are taken care of later in the function (in the truncation part)
    if (debug == true) {cout << "iBinY = " << iBinY << endl;}
    yLowEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY);
    yHighEdgeRebin = H2D_hist_rebinned.GetYaxis()->GetBinLowEdge(iBinY+1);
    yLowEdgeRebin < binsY[0] ? ibinY_low = 0 : ibinY_low = H2D_hist_temp->GetYaxis()->FindBin(yLowEdgeRebin); // according to definition of getbinlowedge, it does work for under/overflows by giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
    yHighEdgeRebin > binsY[nBinsY] ? ibinY_high =  H2D_hist_temp->GetNbinsY()+1 : ibinY_high = H2D_hist_temp->GetYaxis()->FindBin(yHighEdgeRebin) - 1;

    for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // overflow and underflow are taken care of later in the function (in the truncation part)
      xLowEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX);
      xHighEdgeRebin = H2D_hist_rebinned.GetXaxis()->GetBinLowEdge(iBinX+1);
      xLowEdgeRebin < binsX[0] ? ibinX_low = 0 : ibinX_low = H2D_hist_temp->GetXaxis()->FindBin(xLowEdgeRebin); // according to definition of getbinlowedge, this takes into account under/overflows thanks to FindBin() giving the low edge of bin 1 minus the bin width (calculated assuming bins of equal sizes in all the axis); in fact, it works regardless of the value of iBin even for Nbin+2 or more (https://root.cern.ch/doc/master/classTAxis.html#a6d2d17ef00382842f58e88f356516a0d)
      xHighEdgeRebin > binsX[nBinsX] ? ibinX_high =  H2D_hist_temp->GetNbinsX()+1 : ibinX_high = H2D_hist_temp->GetXaxis()->FindBin(xHighEdgeRebin) - 1;

      // dividing by prior integral over the xbin to make the weighting done before a probability
      priorBinNorm = H1D_priorSpectrum->IntegralAndError(ibinX_low, ibinX_high, priorBinNormError);
      for(int iBinYFine = ibinY_low; iBinYFine <= ibinY_high; iBinYFine++){
        for(int iBinXFine = ibinX_low; iBinXFine <= ibinX_high; iBinXFine++){
          H2D_hist_content = H2D_hist_temp->GetBinContent(iBinXFine, iBinYFine);
          H2D_hist_contentError = H2D_hist_temp->GetBinError(iBinXFine, iBinYFine);
          abs(H2D_hist_content) < GLOBAL_epsilon ? H2D_hist_contentErrorA = 0 : H2D_hist_contentErrorA = H2D_hist_contentError*H2D_hist_contentError / (H2D_hist_content*H2D_hist_content);
          abs(priorBinNorm) < GLOBAL_epsilon     ? H2D_hist_contentErrorB = 0 : H2D_hist_contentErrorB = priorBinNormError*priorBinNormError / (priorBinNorm*priorBinNorm);
          abs(priorBinNorm) < GLOBAL_epsilon ? H2D_hist_content = 0 : H2D_hist_content = H2D_hist_content * 1./priorBinNorm; // do I really give the value 0 if denominator is 0 ?
          H2D_hist_temp->SetBinContent(iBinXFine, iBinYFine, H2D_hist_content);
          H2D_hist_temp->SetBinError(iBinXFine, iBinYFine, sqrt(H2D_hist_content*H2D_hist_content * (H2D_hist_contentErrorA + H2D_hist_contentErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
        }
      }
      hist2D_rebinned_content = H2D_hist_temp->IntegralAndError(ibinX_low, ibinX_high, ibinY_low, ibinY_high, hist2D_rebinned_error);
      H2D_hist_rebinned.SetBinContent(iBinX, iBinY, hist2D_rebinned_content);
      H2D_hist_rebinned.SetBinError(iBinX, iBinY, hist2D_rebinned_error);
      if (debug == true) {cout << "ibinX_low = " << ibinX_low << ", ibinX_high = " << ibinX_high << ", ibinY_low = " << ibinY_low << ", ibinY_high = " << ibinY_high << "         --------          H2D_hist_rebinned(" << iBinX <<", "<< iBinY <<") = " << H2D_hist_rebinned.GetBinContent(iBinX, iBinY) << "                      --- priorBinNorm = " << priorBinNorm << "                                (xLowEdgeRebin, xHighEdgeRebin) = (" << xLowEdgeRebin << ", " << xHighEdgeRebin << "), (yLowEdgeRebin, yHighEdgeRebin) = (" << yLowEdgeRebin << ", " << yHighEdgeRebin << ")" << endl;}
    }
  }


  if (debug == true) {
    // il semble se passer quelque chose à 30GeV et 120GeV, limites des recbins;
    // en fait, il y a des démarcations pour chaque bin de recbins; elles sont juste plus visible pour les grand bins puisque le steeply falling spectrum as plus de temps pour diminuer lorsque l'on va de gauche à droite
    // ça me semble ok et normal, pas pathologique

    TString* pdfName_preRebin = new TString("H2D_hist_temp");
    TString* pdfName_preRebin_logz = new TString("H2D_hist_temp_logz");

    TString textContext_preRebin((TString)"");


    Draw_TH2_Histogram(H2D_hist_temp, textContext_preRebin, pdfName_preRebin, texPtJetRec, texPtJetGen, texCollisionDataInfo, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(H2D_hist_temp, textContext_preRebin, pdfName_preRebin_logz, texPtJetRec, texPtJetGen, texCollisionDataInfo, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  std::stringstream ss;
  ss << H2D_hist->GetName() << "_RebinVariableBins2D";
  TString histName((TString)ss.str());

  if (debug == true) {cout << "___________ end RebinVariableBins2D" << endl;}

  return H2D_hist_rebinned;
}
TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum){
  return RebinVariableBins2D_PriorWeightedBinMerging(H2D_hist, nBinsX, nBinsY, binsX, binsY, H1D_priorSpectrum, false);
}



// here if I dont include the under/overflows I have a weird behaviour at low pt gen; but then that looks like what marta actually did; see her non orig response histogram, also https://github.com/alisw/AliPhysics/blob/5e4a7731040f3f2739d98091d5549a44c2ae4e44/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx#L769
void NormaliseYSlicesToOne(TH2D* H2D_hist){ 
  // normalise only what is inside the gen bins, not outside
  // see AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin

  cout << "Normalising y-slices" << endl;
  // Normalisation of the H2D_hist 2D histogram: each pt gen slice is normalised to unity, takes underflow into account if it's present
  double genSliceNorm = 1;
  double genSliceNormError;
  // double dpT;

  double binContent, binError, binErrorA, binErrorB;
  for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
    genSliceNorm = H2D_hist->IntegralAndError(0, H2D_hist->GetNbinsX()+1, iBinY, iBinY, genSliceNormError); // in AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin it goes from 1 to N but  (1,N) in there is the fine matrix min and max bins; while here we already coarsened the matrix, and so that NFine is in the overflow and 1 in the underflow
    for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
      abs(H2D_hist->GetBinContent(iBinX, iBinY)) < GLOBAL_epsilon ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY));
      abs(genSliceNorm) < GLOBAL_epsilon                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm);
      abs(genSliceNorm) < GLOBAL_epsilon ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    }
  }
} 
void NormaliseYSlicesToOneNoUnderOverFlows(TH2D* H2D_hist){ 
  cout << "Normalising y-slices" << endl;
  double genSliceNorm = 1;
  double genSliceNormError;
  double binContent, binError, binErrorA, binErrorB;
  for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){
    genSliceNorm = H2D_hist->IntegralAndError(1, H2D_hist->GetNbinsX(), iBinY, iBinY, genSliceNormError);

    for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){ // 0 and n+1 take underflow and overflow into account
      abs(H2D_hist->GetBinContent(iBinX, iBinY)) < GLOBAL_epsilon ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY));
      abs(genSliceNorm) < GLOBAL_epsilon                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm);
      abs(genSliceNorm) < GLOBAL_epsilon ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    }
  }
} 

void NormaliseXSlicesToOne(TH2D* H2D_hist){ 

    // gustavo
    // nbinsY = H2D_sigmapt_pt[iDataset]->GetNbinsY();
    // for(int j = 1; j <= H2D_sigmapt_pt[iDataset]->GetNbinsX(); j++) {
    //   TH1D* temp1 = h2D->ProjectionY(Form("Bin%d",j),j,j);
    //   float scale1 = temp1 -> Integral(-1,-1);
      
    //   for(int i = 1; i <= nbinsY; i++) {
    //     if ( scale1 > 0 ) {
    //       h2D->SetBinContent(j,i, h2D->GetBinContent(j,i)/scale1);
    //       h2D->SetBinError  (j,i, h2D->GetBinError  (j,i)/scale1);
    //     } else {
    //       h2D->SetBinContent(j,i, 0);
    //       h2D->SetBinError  (j,i, 0);       
    //     }
    //   } // y bin loop 
    // } // x bin loop

  cout << "Normalising x-slices" << endl;
  // Normalisation of the H2D_hist 2D histogram: each pt rec slice is normalised to unity, takes underflow into account if it's present
  double genSliceNorm = 1;
  double genSliceNormError;
  double binContent, binError, binErrorA, binErrorB;
  for(int iBinX = 0; iBinX <= H2D_hist->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    genSliceNorm = H2D_hist->IntegralAndError(iBinX, iBinX, 0, H2D_hist->GetNbinsY()+1, genSliceNormError);
    for(int iBinY = 0; iBinY <= H2D_hist->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      abs(H2D_hist->GetBinContent(iBinX, iBinY))  < GLOBAL_epsilon ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY));
      abs(genSliceNorm)  < GLOBAL_epsilon                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm);
      abs(genSliceNorm)  < GLOBAL_epsilon ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    }
  }
} 
void NormaliseXSlicesToOneNoUnderOverFlows(TH2D* H2D_hist){ 
  // normalise only what is inside the gen bins, not outside
  // see AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin

  cout << "Normalising x-slices" << endl;
  // Normalisation of the H2D_hist 2D histogram: each pt rec slice is normalised to unity, takes underflow into account if it's present
  double genSliceNorm = 1;
  double genSliceNormError;
  double binContent, binError, binErrorA, binErrorB;
  for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){ // 0 and n+1 take underflow and overflow into account
    genSliceNorm = H2D_hist->IntegralAndError(iBinX, iBinX, 1, H2D_hist->GetNbinsY(), genSliceNormError);
    for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){ // 0 and n+1 take underflow and overflow into account
      abs(H2D_hist->GetBinContent(iBinX, iBinY)) < GLOBAL_epsilon ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY));
      abs(genSliceNorm) < GLOBAL_epsilon                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm);
      abs(genSliceNorm) < GLOBAL_epsilon ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    }
  }
} 

void TransformRawResponseToYieldResponse(TH2D* H2D_hist){ 
  // normalise only what is inside the gen bins, not outside

  cout << "TransformRawResponseToYieldResponse: Dividing Response by bin widths" << endl;

  double widthX, widthY;
  for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){
    widthX = H2D_hist->GetXaxis()->GetBinWidth(iBinX);
    for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){
      widthY = H2D_hist->GetYaxis()->GetBinWidth(iBinY);

      H2D_hist->SetBinContent(iBinX, iBinY, H2D_hist->GetBinContent(iBinX, iBinY) *1./widthX *1./widthY);
      H2D_hist->SetBinError(iBinX, iBinY, H2D_hist->GetBinError(iBinX, iBinY) *1./widthX *1./widthY); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    }
  }
} 



TH2* NormalizeResponsMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior) { // aliphysics version; shorter because doesn't take care of uncertainties
  // Normalize such that the Y projection is the prior

  double intPrior = hPrior->Integral();//"width");
  for (Int_t jbin = 1; jbin <= h2RM->GetNbinsY(); jbin++) {
      for (Int_t ibin = 1; ibin <= h2RM->GetNbinsX(); ibin++) {
	double content = h2RM->GetBinContent(ibin,jbin);
	h2RM->SetBinContent(ibin,jbin,hPrior->GetBinContent(jbin)/hPrior->GetBinWidth(jbin)/intPrior*content);
    }
  }

  return h2RM;
}

void WeightMatrixWithPrior(TH2D* H2D_hist, TH1D* priorSpectrum, bool doPriorDivision = false){
  cout << "Weighting y-slices" << endl;

  double binContent, binError, binErrorA, binErrorB;


  // multiplication
  double errorA, errorB;
  for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){ // 0 and n+1 take underflow and overflow into account
      binContent = H2D_hist->GetBinContent(iBinX, iBinY) * priorSpectrum->GetBinContent(iBinY);
      errorA = priorSpectrum->GetBinContent(iBinY)*priorSpectrum->GetBinContent(iBinY) * H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY);
      errorB = H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY) * priorSpectrum->GetBinError(iBinY)*priorSpectrum->GetBinError(iBinY);
      H2D_hist->SetBinContent(iBinX, iBinY, binContent);
      H2D_hist->SetBinError(iBinX, iBinY, sqrt(errorA + errorB)); // sigma(AB)2 = sigma(A)2 B2 + sigma(B)2 A2
    }
  }

  // division
  if (doPriorDivision){ 
    double genSliceNorm, genSliceNormError;
    // double genSliceNorm = priorSpectrum->IntegralAndError(0, -1, genSliceNormError);
    genSliceNorm = priorSpectrum->IntegralAndError(1, priorSpectrum->GetNbinsX(), genSliceNormError); // I think we only prior the 1, N bins

    for(int iBinX = 1; iBinX <= H2D_hist->GetNbinsX(); iBinX++){ // 0 and n+1 take underflow and overflow into account
      for(int iBinY = 1; iBinY <= H2D_hist->GetNbinsY(); iBinY++){ // 0 and n+1 take underflow and overflow into account
        abs(H2D_hist->GetBinContent(iBinX, iBinY)) < GLOBAL_epsilon ? binErrorB = 0 : binErrorB = H2D_hist->GetBinError(iBinX, iBinY)*H2D_hist->GetBinError(iBinX, iBinY) / (H2D_hist->GetBinContent(iBinX, iBinY)*H2D_hist->GetBinContent(iBinX, iBinY)); // sigma(A)2 /A2
        abs(genSliceNorm) < GLOBAL_epsilon                          ? binErrorA = 0 : binErrorA = genSliceNormError*genSliceNormError / (genSliceNorm*genSliceNorm); // sigma(B)2 /B2
        abs(genSliceNorm) < GLOBAL_epsilon ? binContent = 0 : binContent = H2D_hist->GetBinContent(iBinX, iBinY) *1./genSliceNorm; // do I really give the value 0 if denominator is 0 ? 
        H2D_hist->SetBinContent(iBinX, iBinY, binContent);
        H2D_hist->SetBinError(iBinX, iBinY, sqrt(binContent*binContent * (binErrorA + binErrorB)));  // sigma(A/B)2 = (A/B)2 * sigma(A)2 /A2 + sigma(B)2 /B2
      }
    }
  }

  // note: if I can caluclate sigma(priorSpectrum->Integral() - priorSpectrum->GetBinContent(iBinY)) then I can have the sigma calculated at once and have some of it be cancelled; could do integral of ibin1 to ibinY-1, plus ibinY+1 to N ? 
  // that is doable if I really am using the uncertainties on the matrices, which I don't think is actually being used
}

TH2D GetTransposeHistogram(TH2D* inputHist){
  int nBinsX = inputHist->GetNbinsX();
  int nBinsY = inputHist->GetNbinsY();
  std::vector<double> vectorBinsX = GetTH1Bins(inputHist->ProjectionX(inputHist->GetName()+(TString)"projX", 1, nBinsX));
  std::vector<double> vectorBinsY = GetTH1Bins(inputHist->ProjectionY(inputHist->GetName()+(TString)"projY", 1, nBinsY));
  double* binsX = &vectorBinsX[0];
  double* binsY = &vectorBinsY[0];

  TH2D transposedHist("inputHist_transposed", "inputHist_transposed", nBinsY, binsY, nBinsX, binsX);
  transposedHist.Sumw2();
  for(int iBinX = 0; iBinX <= nBinsX+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinY = 0; iBinY <= nBinsY+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      transposedHist.SetBinContent(iBinY, iBinX, inputHist->GetBinContent(iBinX, iBinY));
      transposedHist.SetBinError(iBinY, iBinX, inputHist->GetBinError(iBinX, iBinY));
    }
  }
  return transposedHist;
}

TH2D GetMatrixProductTH2xTH2(TH2D* histA, TH2D* histB){
  int nBinsX_A = histA->GetNbinsX();
  int nBinsY_A = histA->GetNbinsY();
  int nBinsX_B = histB->GetNbinsX();
  int nBinsY_B = histB->GetNbinsY();


  if (nBinsY_A != nBinsX_B) {
    cout << "#########################################################################################" << endl;
    cout << "###################### MATRIX PRODUCT IMPOSSIBLE DUE TO DIMENSIONS ######################" << endl;
    cout << "#########################################################################################" << endl;
    cout << "matrix A times B, with A(" << nBinsY_A << " x " << nBinsX_A << ") and B(" << nBinsY_B << " x " << nBinsX_B << ")" << endl; 
  }

  int nBinsK = nBinsY_A;

  std::vector<double> vectorBinsX_AB = GetTH1Bins(histB->ProjectionX(histB->GetName()+(TString)"projX", 1, nBinsX_B));
  std::vector<double> vectorBinsY_AB = GetTH1Bins(histA->ProjectionY(histA->GetName()+(TString)"projY", 1, nBinsY_A));
  double* binsX_AB = &vectorBinsX_AB[0];
  double* binsY_AB = &vectorBinsY_AB[0];
  int nBinsX_AB = nBinsX_B;
  int nBinsY_AB = nBinsY_A;

  TH2D histAB("histAB", "histAB", nBinsX_B, binsX_AB, nBinsY_A, binsY_AB);
  histAB.Sumw2();
  // Roounfold sees the element of a hist histogram(i,j) as the element R(i,j) of a matrix; thus no transformation has to be done for the indices or axes of the histograms when doing the multiplication

  double productContent_iBinX_iBinY;
  double productError2_iBinX_iBinY;
  for(int iBinX = 0; iBinX <= nBinsX_AB+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    for(int iBinY = 0; iBinY <= nBinsY_AB+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
      productContent_iBinX_iBinY = 0;
      productError2_iBinX_iBinY = 0;
      if (0 < iBinX && iBinX < nBinsX_AB+1 && 0 < iBinY && iBinY < nBinsY_AB+1) { // reason we do this, is because we don't want the under/ofverflows spread inside the matrix when we multiply, which happens when we do ; by separating it, we can still have the under/overflows calculated correctly in the else case, so that we can use them for the kinematic efficiency
        for(int iBinK = 1; iBinK <= nBinsK; iBinK++){ // 1 and n don't take underflow and overflow into account
          productContent_iBinX_iBinY += histA->GetBinContent(iBinX, iBinK) * histB->GetBinContent(iBinK, iBinY); 
          productError2_iBinX_iBinY += pow(histB->GetBinContent(iBinK, iBinY), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histB->GetBinError(iBinK, iBinY), 2); // simple sigma_ab = b2sigma_a2 + a2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated    
          // if (iBinX < 7) {
          //   cout << "A(" << iBinX << ", " << iBinK << ") = " << histA->GetBinContent(iBinX, iBinK) << ", B(" << iBinK << ", " << iBinY << ") = " << histB->GetBinContent(iBinK, iBinY) << endl;
          // }
        }
      } else {
        for(int iBinK = 0; iBinK <= nBinsK+1; iBinK++){ // 0 and n+1 take underflow and overflow into account
          productContent_iBinX_iBinY += histA->GetBinContent(iBinX, iBinK) * histB->GetBinContent(iBinK, iBinY); 
          productError2_iBinX_iBinY += pow(histB->GetBinContent(iBinK, iBinY), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histB->GetBinError(iBinK, iBinY), 2); // simple sigma_ab = b2sigma_a2 + a2sigma_b2 ; that assumes there are no correlations; here it s background fluct from PbPB sim, and detector effects from a pp sim, so we can maybe say theyre not correlated    
        }
      }
      histAB.SetBinContent(iBinX, iBinY, productContent_iBinX_iBinY);
      histAB.SetBinError(iBinX, iBinY, sqrt(productError2_iBinX_iBinY));  
      // if (iBinX < 7) {
      //   cout << ". . . . . . . . . . . . . . . . . . . . . . ." << endl;
      //   cout << "AB(" << iBinX << ", " << iBinY << ") = " << histAB.GetBinContent(iBinX, iBinY) << endl;
      //   cout << "==========================================================" << endl;
      // }
    }
  }
  return histAB;
}

TH1D GetMatrixVectorProductTH2xTH1(TH2D* histA, TH1D* histU){
  // bit different from TH2 x TH2 as a TH1 vector is a (1,N) TH2 and not a (N,1) one
  int nBinsX_A = histA->GetNbinsX();
  int nBinsY_A = histA->GetNbinsY();
  int nBinsX_U = histU->GetNbinsX();

  if (nBinsY_A != nBinsX_U) {
    cout << "#########################################################################################" << endl;
    cout << "###################### MATRIX PRODUCT IMPOSSIBLE DUE TO DIMENSIONS ######################" << endl;
    cout << "#########################################################################################" << endl;
  }

  int nBinsK = nBinsY_A;
  // int nBinsX_AB = nBinsX_B;
  // int nBinsY_AB = nBinsY_A;
  
  std::vector<double> vectorBinsX_AU = GetTH1Bins(histA->ProjectionX(histA->GetName()+(TString)"projY", 1, nBinsY_A));
  double* binsX_AU = &vectorBinsX_AU[0];
  int nBinsX_AU = nBinsX_A;

  TH1D histAU("histAU", "histAU", nBinsX_AU, binsX_AU);
  histAU.Sumw2();

  // Roounfold sees the element of a hist histogram(i,j) as the element R(i,j) of a matrix; thus no transformation has to be done for the indices or axes of the histograms when doing the multiplication

  double productContent_iBinX;
  double productError2_iBinX;
  for(int iBinX = 0; iBinX <= nBinsX_AU+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
    productContent_iBinX = 0;
    productError2_iBinX = 0;
    if (0 < iBinX && iBinX < nBinsX_AU+1) { // reason we do this, is because we don't want the under/ofverflows spread inside the matrix when we multiply, which happens when we do ; by separating it, we can still have the under/overflows calculated correctly in the else case, so that we can use them for the kinematic efficiency
      for(int iBinK = 1; iBinK <= nBinsK; iBinK++){ // 0 and n+1 take underflow and overflow into account
        productContent_iBinX += histA->GetBinContent(iBinX, iBinK) * histU->GetBinContent(iBinK); 
        productError2_iBinX += pow(histU->GetBinContent(iBinK), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histU->GetBinError(iBinK), 2); // simple sigma_ab = b2sigma_a2 + a2sigma_b2 ; that assumes there are no correlations;
        // cout << "-----" << endl;
        // cout << " A(" << iBinX << "," << iBinK << ") = " << histA->GetBinContent(iBinX, iBinK) << ",  U(" << iBinK << ") = " << histU->GetBinContent(iBinK) << endl;
        // cout << "eA(" << iBinX << "," << iBinK << ") = " << histA->GetBinError(iBinX, iBinK) << ", eU(" << iBinK << ") = " << histU->GetBinError(iBinK) << endl;
      }
    } else {
      for(int iBinK = 0; iBinK <= nBinsK+1; iBinK++){ // 0 and n+1 take underflow and overflow into account
        productContent_iBinX += histA->GetBinContent(iBinX, iBinK) * histU->GetBinContent(iBinK); 
        productError2_iBinX += pow(histU->GetBinContent(iBinK), 2)*pow(histA->GetBinError(iBinX, iBinK), 2) + pow(histA->GetBinContent(iBinX, iBinK), 2)*pow(histU->GetBinError(iBinK), 2); // simple sigma_ab = b2sigma_a2 + a2sigma_b2 ; that assumes there are no correlations;
      }
    }
    histAU.SetBinContent(iBinX, productContent_iBinX);
    histAU.SetBinError(iBinX, sqrt(productError2_iBinX));  
  }
  return histAU;
}

bool DivideWithCorrelatedErrors_simpleMax(TH1D* histNumeratorToBeDivided, TH1D* histDenominator) {
  // like ->Divide() this function replaces the first histogram with the ratio
  int nBins = histNumeratorToBeDivided->GetNbinsX();
  std::vector<double> errorsNum, errorsDen;
  for (int iBin = 1; iBin <= nBins; iBin++) {
    errorsNum.push_back(histNumeratorToBeDivided->GetBinError(iBin)/histNumeratorToBeDivided->GetBinContent(iBin));
    errorsDen.push_back(histDenominator->GetBinError(iBin)/histDenominator->GetBinContent(iBin));
  }
  bool divideSuccess = histNumeratorToBeDivided->Divide(histDenominator);
  for (int iBin = 1; iBin <= nBins; iBin++) {
    histNumeratorToBeDivided->SetBinError(iBin, errorsNum[iBin-1] > errorsDen[iBin-1] ? errorsNum[iBin-1] : errorsDen[iBin-1]);
    // cout << "iBin " << iBin << ", eNum = " << errorsNum[iBin-1] << ", eDen = " << errorsDen[iBin-1] << ", divisionError = " << histNumeratorToBeDivided->GetBinError(iBin) << endl;
  }

  return divideSuccess;
}




void Save_PtResponseMatrix(TH2D* &H2D_jetPtResponseMatrix, TString fileName) {
  TFile *newFile = new TFile(fileName+".root", "RECREATE");
  H2D_jetPtResponseMatrix->Write();
  newFile->Close();
}



TH1D* ReweightedRebin(TH1D* h_in, const char* newname, int n_newbins, const double* newbins){
    TH1D* h_out = new TH1D(newname, h_in->GetTitle(),
                            n_newbins, newbins);
    h_out->Sumw2();

    for (int j = 1; j <= n_newbins; j++) {
        double newlo = h_out->GetBinLowEdge(j);
        double newhi = newlo + h_out->GetBinWidth(j);

        double sum_val = 0.0;
        double sum_var = 0.0;  // sum of variances for error propagation
        double sum_width = 0.0;

        for (int i = 1; i <= h_in->GetNbinsX(); i++) {
            double lo = h_in->GetBinLowEdge(i);
            double hi = lo + h_in->GetBinWidth(i);
            double w  = h_in->GetBinWidth(i);

            // check old bin is fully inside new bin
            if (lo >= newlo - 1e-9 && hi <= newhi + 1e-9) {
                // content is dsigma/dpT deta, multiply back by width
                sum_val   += h_in->GetBinContent(i) * w;
                sum_var   += pow(h_in->GetBinError(i) * w, 2);
                sum_width += w;
            }
        }

        double newwidth = h_out->GetBinWidth(j);

        // sanity check: sum_width should equal newwidth
        if (fabs(sum_width - newwidth) > 1e-6)
            std::cerr << "WARNING: bin " << j << " width mismatch: "
                      << sum_width << " vs " << newwidth << std::endl;

        if (newwidth > 0) {
            h_out->SetBinContent(j, sum_val / newwidth);
            h_out->SetBinError(j, sqrt(sum_var) / newwidth);
        }
    }
    return h_out;
}

bool HarmonizeBinning(std::vector<TH1*>& hists, double xLo, double xHi) {
  if (hists.empty()) { std::cerr << "HarmonizeBinning: empty histogram list\n"; return false; }

  auto edges = [](TH1* h){
    std::vector<double> e;
    for (int i=1;i<=h->GetNbinsX();++i) e.push_back(h->GetXaxis()->GetBinLowEdge(i));
    e.push_back(h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
    return e;
  };
  auto approx = [](double a,double b){ return std::fabs(a-b) < 1e-4*(std::fabs(a)+std::fabs(b)+1e-10); };

  // Start from the edges of the first histogram (restricted to [xLo,xHi]),
  // then intersect with every other histogram's edges in turn.
  std::vector<double> common;
  for (double a : edges(hists[0]))
    if (a >= xLo-1e-6 && a <= xHi+1e-6) common.push_back(a);

  for (size_t k = 1; k < hists.size(); ++k) {
    std::vector<double> eK = edges(hists[k]);
    std::vector<double> next;
    for (double a : common)
      for (double b : eK)
        if (approx(a,b)) { next.push_back(a); break; }
    common = std::move(next);
    if (common.size() < 2) {
      std::cerr << "HarmonizeBinning: no overlap in ["<<xLo<<","<<xHi<<"] after intersecting hist "
                 << k << " (" << hists[k]->GetName() << ")\n";
      return false;
    }
  }

  std::sort(common.begin(), common.end());
  common.erase(std::unique(common.begin(),common.end(),approx), common.end());
  if (common.size() < 2) {
    std::cerr << "HarmonizeBinning: no overlap in ["<<xLo<<","<<xHi<<"]\n";
    return false;
  }

  int n = (int)common.size()-1;
  for (auto& h : hists) {
    h = ReweightedRebin((TH1D*)h, Form("%s_h", h->GetName()), n, common.data());
    h->SetDirectory(0);
  }
  return true;
}

#endif
