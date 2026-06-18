#ifndef FITS_C
#define FITS_C

TGraphErrors* GetFunctionTGraphErrorsFromFitResult(double* xRangeFit, TF1* fitFunctionDrawn, TFitResultPtr fitResult, int nPointsGraph = 1000){
  std::vector<double> xAxisGraph= {};
  std::vector<double> yAxisGraph= {};
  std::vector<double> yAxisGraphErrors= {};
  // double* ;

  for(int iPoint = 0; iPoint < nPointsGraph; iPoint++){
    xAxisGraph.push_back(xRangeFit[0]+iPoint*1./nPointsGraph*(xRangeFit[1]-xRangeFit[0]));
    yAxisGraph.push_back(fitFunctionDrawn->Eval(xAxisGraph.back()));
    yAxisGraphErrors.push_back(0);
  }
  double oneSigmaInterval = 0.683;
  fitResult->GetConfidenceIntervals(nPointsGraph, 1, 1, &xAxisGraph[0], &yAxisGraphErrors[0], oneSigmaInterval, false);
  TGraphErrors* fitFunctionTGraphErrors = new TGraphErrors(nPointsGraph, &xAxisGraph[0], &yAxisGraph[0], nullptr, &yAxisGraphErrors[0]);
  return fitFunctionTGraphErrors;
}

TGraphErrors* GetFunctionTGraphErrorsFromCovMatrix(double* xRangeFit, TF1* fitFunctionDrawn, TMatrixDSym* covMatrix, int nPointsGraph = 1000){
  std::vector<double> xAxisGraph= {};
  std::vector<double> yAxisGraph= {};
  std::vector<double> yAxisGraphErrors= {};
  // double* ;

  for(int iPoint = 0; iPoint < nPointsGraph; iPoint++){
    xAxisGraph.push_back(xRangeFit[0]+iPoint*1./nPointsGraph*(xRangeFit[1]-xRangeFit[0]));
    yAxisGraph.push_back(fitFunctionDrawn->Eval(xAxisGraph.back()));
    yAxisGraphErrors.push_back(fitFunctionDrawn->EvalUncertainty(xAxisGraph.back(), covMatrix));
  }
  TGraphErrors* fitFunctionTGraphErrors = new TGraphErrors(nPointsGraph, &xAxisGraph[0], &yAxisGraph[0], nullptr, &yAxisGraphErrors[0]);
  return fitFunctionTGraphErrors;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Fits ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

std::tuple<TF1*, TMatrixDSym, TFitResultPtr> TsallisFit(TH1D* &histogramInput, double* xRangeFit) {
  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *fitFunctionInit;
  TF1 *fitFunctionFinal;
  TF1 *fitFunctionDrawn; // drawn over the full range

  TFitResultPtr fFitResult;

  const int nParameters = 2;
  std::array<double, nParameters> parfitFunctionInit;
  std::array<double, nParameters> parfitFunctionFinal;

  ////////////////////////////////////////////////////////////////////
  //////////////////////////// Fit start /////////////////////////////
  ////////////////////////////////////////////////////////////////////

  // double xHistMax = xRange[1];
  // double xHistMin = xRange[0];
  
  fitFunctionInit = new TF1("fitFunctionInit_", "x*(1+1/([0]*[1])*x)**(-[0])", xRangeFit[0], xRangeFit[1]);
  // fitFunctionInit = new TF1("fitFunctionInit_", "[0]*exp(([1]-x)**[2])", xRangeFit[0], xRangeFit[1]);
  fitFunctionInit->SetParName(0, "n");
  fitFunctionInit->SetParName(1, "T");
  // fitFunctionInit->SetParName(2, "b");
  fitFunctionInit->SetParameters(8, 0.9);

  histogramInput->Fit(fitFunctionInit, "R0QL"); // P: Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the square-root of the bin function value. (WL for weithged likelihood is currently bugged in root, the fit crashes)

  fitFunctionInit->GetParameters(&parfitFunctionInit[0]);

  fitFunctionFinal = new TF1("fitFunctionFinal_", "x*(1+1/([0]*[1])*x)**(-[0])", xRangeFit[0], xRangeFit[1]);
  // fitFunctionFinal = new TF1("fitFunctionFinal_", "[0]*exp(([1]-x)**[2])", xRangeFit[0], xRangeFit[1]);
  fitFunctionInit->SetParName(0, "n");
  fitFunctionInit->SetParName(1, "T");
  // fitFunctionInit->SetParName(2, "b");
  fitFunctionFinal->SetParameters(parfitFunctionInit[0], parfitFunctionInit[1]);
  // fitFunctionFinal->SetParameters(parfitFunctionInit[0], parfitFunctionInit[1], parfitFunctionInit[2]);
  // fitFunctionFinal->SetParLimits(0, 0., 1.1*yHistMax);
  // fitFunctionFinal->SetParLimits(1, -10, 10);
  // fitFunctionFinal->SetParLimits(2, 0.1, 100);

  fFitResult = histogramInput->Fit(fitFunctionFinal, "R0QPS"); // P: Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the square-root of the bin function value. (WL for weithged likelihood is currently bugged in root, the fit crashes)
  // gauss->Draw("same");
  fitFunctionFinal->GetParameters(&parfitFunctionFinal[0]);
  TMatrixDSym covMatrixFit = fFitResult->GetCovarianceMatrix();

  Double_t *pDataSmall = covMatrixFit.GetMatrixArray();
  for (int i = 0; i < nParameters*nParameters; i++) {
    cout << "i = " << i << ", covMatrixFit[i]" << pDataSmall[i] << endl;
  }

  fitFunctionDrawn = new TF1("fitFunctionDrawn_", "x*(1+1/([0]*[1])*x)**(-[0])", xRangeFit[0], xRangeFit[1]);
  // fitFunctionDrawn = new TF1("fitFunctionDrawn_", "[0]*exp(([1]-x)**[2])", xRangeFit[0], xRangeFit[1]);
  fitFunctionDrawn->SetParameters(parfitFunctionFinal[0], parfitFunctionFinal[1]);
  // fitFunctionDrawn->SetParameters(parfitFunctionFinal[0], parfitFunctionFinal[1], parfitFunctionFinal[2]);
  // fitFunctionDrawn->SetParameters(5, 0.9);

  // cout << "init:  n = " << parfitFunctionInit[0] << ", T = " << parfitFunctionInit[1]<< endl;
  // cout << "final: n = " << parfitFunctionFinal[0] << ", T = " << parfitFunctionFinal[1]<< endl;

  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> fitFunctionAndFitParams(fitFunctionDrawn, covMatrixFit, fFitResult);
  return fitFunctionAndFitParams;
}

std::tuple<TF1*, TMatrixDSym, TFitResultPtr> ExponentialFitWithLogTransfo(TH1D* &histogramInput, double* xRangeFit) {
  // Transforming distribution wiht application of ln()
  TH1D* H1D_logInput = (TH1D*)histogramInput->Clone((TString)histogramInput->GetName()+(TString)"logTransfo");
  for (int iBin = 1; iBin < histogramInput->GetNbinsX() +1; iBin++) {
    H1D_logInput->SetBinContent(iBin, std::log(H1D_logInput->GetBinContent(iBin)));
    H1D_logInput->SetBinError(iBin, H1D_logInput->GetBinError(iBin)/H1D_logInput->GetBinContent(iBin)); // df = dx/x if f=ln(x)
  }

  ////////////////////////////////// Fit initialisation //////////////////////////////////
  //Fit tools initialisation
  TF1 *fitFunctionInit;
  TF1 *fitFunctionFinal;
  TF1 *fitFunctionDrawn; // drawn over the full range

  TFitResultPtr fFitResult;

  const int nParameters = 3;
  std::array<double, nParameters> parfitFunctionInit;
  std::array<double, nParameters> parfitFunctionFinal;

  ////////////////////////////////////////////////////////////////////
  //////////////////////////// Fit start /////////////////////////////
  ////////////////////////////////////////////////////////////////////

  // double xHistMax = xRange[1];
  // double xHistMin = xRange[0];
  fitFunctionInit = new TF1("fitFunctionInit_", "[0] + [1]*log(x) + [2]*x*log(x)", xRangeFit[0], xRangeFit[1]); //eventually x**[3] should be replaced with fit to PYTHIA fulljet spectra
  fitFunctionInit->SetParName(0, "C");
  fitFunctionInit->SetParName(1, "expb");
  fitFunctionInit->SetParName(2, "expa");

  int x30 = H1D_logInput->GetXaxis()->FindBin(30.0);
  double y30 = H1D_logInput->GetBinContent(x30);

  fitFunctionInit->SetParameters(std::log(std::max(1e-12, y30)), -5.0, 0.0);

  // M0QR was used by archita, R0QL by me; Archita's better
  H1D_logInput->Fit(fitFunctionInit, "M0QR"); // 

  // parfitFunctionInit = fitFunctionInit->GetParameters();
  fitFunctionInit->GetParameters(&parfitFunctionInit[0]);

  // fitFunctionFinal = new TF1("fitFunctionFinal_", "[0] + [1]*log(x) + [2]*x*log(x)", xRangeFit[0], xRangeFit[1]);
  fitFunctionFinal = new TF1("fitFunctionDrawn_", "exp([0])*pow(x, [1]+[2]*x)", xRangeFit[0], xRangeFit[1]);
  fitFunctionInit->SetParName(0, "C");
  fitFunctionInit->SetParName(1, "expb");
  fitFunctionInit->SetParName(2, "expa");
  fitFunctionFinal->SetParameters(parfitFunctionInit[0], parfitFunctionInit[1], parfitFunctionInit[2]);
  // fitFunctionFinal->SetParLimits(0, 0., 1.1*yHistMax);
  // fitFunctionFinal->SetParLimits(1, -10, 10);
  // fitFunctionFinal->SetParLimits(2, 0.1, 100);

  fFitResult = histogramInput->Fit(fitFunctionFinal, "M0QRS");

  // parfitFunctionFinal = fitFunctionFinal->GetParameters();
  fitFunctionFinal->GetParameters(&parfitFunctionFinal[0]);

  TMatrixDSym covMatrixFit = fFitResult->GetCovarianceMatrix();

  Double_t *pDataSmall = covMatrixFit.GetMatrixArray();
  for (int i = 0; i < nParameters*nParameters; i++) {
    cout << "covMatrixFit[" << i << "] =" << pDataSmall[i] << endl;
  }

  fitFunctionDrawn = new TF1("fitFunctionDrawn_", "exp([0])*pow(x, [1]+[2]*x)", xRangeFit[0], xRangeFit[1]);
  fitFunctionDrawn->SetParName(0, "C");
  fitFunctionDrawn->SetParName(1, "expb");
  fitFunctionDrawn->SetParName(2, "expa");
  fitFunctionDrawn->SetParameters(parfitFunctionFinal[0], parfitFunctionFinal[1], parfitFunctionFinal[2]);
  // fitFunctionDrawn->SetParameters(5, 0.9);

  // cout << "init:  n = " << parfitFunctionInit[0] << ", T = " << parfitFunctionInit[1]<< endl;
  // cout << "final: n = " << parfitFunctionFinal[0] << ", T = " << parfitFunctionFinal[1]<< endl;

  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> fitFunctionAndFitParams(fitFunctionDrawn, covMatrixFit, fFitResult);
  return fitFunctionAndFitParams;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Rebin with Fit ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

std::pair<TH1D*, TGraphErrors*> RebinWithTsallisFit(TH1D* &histogramInput, int nBinsX, double* binsX, double* xRangeFit, TString histName) {
  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tsallisFitFunctionResult = TsallisFit(histogramInput, xRangeFit);
  TF1* fitFunctionDrawn = std::get<0>(tsallisFitFunctionResult);
  TFitResultPtr fitResult = std::get<2>(tsallisFitFunctionResult);
  TGraphErrors* fitFunctionTGraphErrors = GetFunctionTGraphErrorsFromFitResult(xRangeFit, fitFunctionDrawn, fitResult);

  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Rebin of input histogram /////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  TH1D* histogramRebinned = new TH1D(histName+(TString)"_rebinned", histName+(TString)"_rebinned", nBinsX, binsX);
  for(int iBin = 0; iBin < nBinsX; iBin++){
    // histogramRebinned->SetBinContent(iBin, histogramInput->GetBinContent(iBin)); // Getting bin center here not ideal; should try to read and apply "Where to stick your data points: The treatment of measurements within wide bins"
    histogramRebinned->SetBinContent(iBin, fitFunctionDrawn->Eval(histogramRebinned->GetXaxis()->GetBinCenter(iBin))); // Getting bin center here not ideal; should try to read and apply "Where to stick your data points: The treatment of measurements within wide bins"
    // cout << "histogramRebinned(" << iBin << ") = " << histogramRebinned->GetBinContent(iBin) << endl; 
    // histogramRebinned->SetBinError(iBin, fitFunctionDrawn->EvalUncertainty(histogramRebinned->GetXaxis()->GetBinCenter(iBin), nullptr));
    double oneSigmaInterval = 0.683;
    double errorEval[1] = {0};
    double xEval[1] = {(double)histogramRebinned->GetXaxis()->GetBinCenter(iBin)};
    fitResult->GetConfidenceIntervals(1, 1, 1, xEval, errorEval, oneSigmaInterval, false);
    histogramRebinned->SetBinError(iBin, errorEval[0]);
  }

  std::pair<TH1D*, TGraphErrors*> rebinResultAndFitFunction(histogramRebinned, fitFunctionTGraphErrors);
  return rebinResultAndFitFunction;
}


std::pair<TH1D*, TGraphErrors*> RebinWithFit(TH1D* &histogramInput, int nBinsX, double* binsX, double* xRangeFit, TString histName, std::tuple<TF1*, TMatrixDSym, TFitResultPtr> fitFunctionResult) {
  TF1* fitFunctionDrawn = std::get<0>(fitFunctionResult);
  TFitResultPtr fitResult = std::get<2>(fitFunctionResult);
  TGraphErrors* fitFunctionTGraphErrors = GetFunctionTGraphErrorsFromFitResult(xRangeFit, fitFunctionDrawn, fitResult);

  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Rebin of input histogram /////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  TH1D* histogramRebinned = new TH1D(histName+(TString)"_rebinned", histName+(TString)"_rebinned", nBinsX, binsX);
  for(int iBin = 0; iBin < nBinsX; iBin++){
    // histogramRebinned->SetBinContent(iBin, histogramInput->GetBinContent(iBin)); // Getting bin center here not ideal; should try to read and apply "Where to stick your data points: The treatment of measurements within wide bins"
    histogramRebinned->SetBinContent(iBin, fitFunctionDrawn->Eval(histogramRebinned->GetXaxis()->GetBinCenter(iBin))); // Getting bin center here not ideal; should try to read and apply "Where to stick your data points: The treatment of measurements within wide bins"
    // cout << "histogramRebinned(" << iBin << ") = " << histogramRebinned->GetBinContent(iBin) << endl; 
    // histogramRebinned->SetBinError(iBin, fitFunctionDrawn->EvalUncertainty(histogramRebinned->GetXaxis()->GetBinCenter(iBin), nullptr));
    double oneSigmaInterval = 0.683;
    double errorEval[1] = {0};
    double xEval[1] = {(double)histogramRebinned->GetXaxis()->GetBinCenter(iBin)};
    fitResult->GetConfidenceIntervals(1, 1, 1, xEval, errorEval, oneSigmaInterval, false);
    histogramRebinned->SetBinError(iBin, errorEval[0]);
  }

  std::pair<TH1D*, TGraphErrors*> rebinResultAndFitFunction(histogramRebinned, fitFunctionTGraphErrors);
  return rebinResultAndFitFunction;
}



// Fit a histogram with a double Tsallis-like function and return everything needed to propagate uncertainties
std::tuple<TF1*, TMatrixDSym, TFitResultPtr> FitDoubleTsallis(TH1D* &histogramInput, int nBinsX, double* binsX, double* xRangeFit) {
  TF1 *fitFunctionInit;
  TF1 *fitFunctionFinal;
  TF1 *fitFunctionDrawn; // drawn over the full range
  TFitResultPtr fFitResult;

  // double parfitFunctionInit[4];
  // double parfitFunctionFinal[4];
  double parfitFunctionInit[8];
  double parfitFunctionFinal[8];
  const char* doubleTsallis = "([2]+[3]*x)*pow(1 + x/([0]*[1]), -[1]) + ([6]+[7]*x)*pow(1 + x/([4]*[5]), -[5])";
  // const char* doubleTsallis = "([2]+[3]*x)*pow(1 + x/([0]*[1]), -[1])";


  ////////////////////////////////////////////////////////////////////
  //////////////////////////// Fit start /////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  fitFunctionInit = new TF1("fitFunctionInit_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  
  // Set parameter names
  fitFunctionInit->SetParName(0, "p0");
  fitFunctionInit->SetParName(1, "p1");
  fitFunctionInit->SetParName(2, "p2");
  fitFunctionInit->SetParName(3, "p3");
  fitFunctionInit->SetParName(4, "p4");
  fitFunctionInit->SetParName(5, "p5");
  fitFunctionInit->SetParName(6, "p6");
  fitFunctionInit->SetParName(7, "p7");

  fitFunctionInit->SetParameters(0.44,  5.53,   4.58,  0.09,  0.63,  9.48,  3.78, -0,56);  
  //                             p0,   p1,   p2,  p3,  p4,  p5,  p6,  p7

  fitFunctionInit->SetParLimits(0, 0.05, 1.0);
  fitFunctionInit->SetParLimits(1, 4.0, 6.0);
  fitFunctionInit->SetParLimits(2, 3.0, 5.0);
  fitFunctionInit->SetParLimits(3, 0.0, 0.1);
  fitFunctionInit->SetParLimits(4, 0.05, 1.0);
  fitFunctionInit->SetParLimits(5, 9.0, 10.0);
  fitFunctionInit->SetParLimits(6, 2.0, 4.0);
  fitFunctionInit->SetParLimits(7, -0.9, 0.5);

  histogramInput->Fit(fitFunctionInit, "SR0Q"); // R = fit range, Q = quiet, L = likelihood
  fitFunctionInit->GetParameters(&parfitFunctionInit[0]); // Save initial parameters

  fitFunctionFinal = new TF1("fitFunctionFinal_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  
  for(int i=0; i<8; i++) fitFunctionFinal->SetParameter(i, parfitFunctionInit[i]);
  // fitFunctionFinal->SetParLimits(0, 0.05, 1.0);
  // fitFunctionFinal->SetParLimits(1, 4.0, 6.0);
  // fitFunctionFinal->SetParLimits(2, 2.0, 5.0);
  // fitFunctionFinal->SetParLimits(3, -0.5, 0.1);
  // fitFunctionFinal->SetParLimits(4, 0.05, 1.0);
  // fitFunctionFinal->SetParLimits(5, 9.0, 12.0);
  // fitFunctionFinal->SetParLimits(6, 1.0, 4.0);
  // fitFunctionFinal->SetParLimits(7, -0.9, 0.5);
  

  fFitResult = histogramInput->Fit(fitFunctionFinal, "RS");  
  fitFunctionFinal->GetParameters(&parfitFunctionFinal[0]);

  // Check covariance availability
  TMatrixDSym covMatrixFit; // default empty
  if (fFitResult && fFitResult->CovMatrixStatus() == 3) {
      covMatrixFit = fFitResult->GetCovarianceMatrix();
  } else {
      // std::cout << "Warning: Covariance matrix not available!" << std::endl;
      std::cout << "Covariance Status: " << fFitResult->CovMatrixStatus() << std::endl;
  }

  // TMatrixDSym covMatrixFit = fFitResult->GetCovarianceMatrix();

  // Double_t *pDataSmall = covMatrixFit.GetMatrixArray();
  // for (int i = 0; i < 2*2; i++) {
  //   cout << "i = " << i << ", covMatrixFit[i]" << pDataSmall[i] << endl;
  // }

  fitFunctionDrawn = new TF1("fitFunctionDrawn_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  for(int i=0; i<8; i++) fitFunctionDrawn->SetParameter(i, parfitFunctionFinal[i]);

  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> fitFunctionAndFitParams(fitFunctionDrawn, covMatrixFit, fFitResult);
  return fitFunctionAndFitParams;
}

// Use the fitted function to rebin a histogram and propagate fit uncertainties to the new bins
std::tuple<TH1D*, TGraphErrors*, TF1*> RebinWithDoubleTsallisFit(TH1D* &histogramInput, int nBinsX, double* binsX, double* xRangeFit) {
  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tsallisFitFunctionResult = FitDoubleTsallis(histogramInput, nBinsX, binsX, xRangeFit);
  TF1* fitFunctionDrawn = std::get<0>(tsallisFitFunctionResult);
  TFitResultPtr fitResult = std::get<2>(tsallisFitFunctionResult);
  TGraphErrors* fitFunctionTGraphErrors = GetFunctionTGraphErrorsFromFitResult(xRangeFit, fitFunctionDrawn, fitResult);
  
  //////////////////////////// Rebin of input histogram /////////////////////////////

  TH1D* histogramRebinned = new TH1D("Unfolded: fit sampling", "Unfolded: fit sampling", nBinsX, binsX);
  for(int iBin = 0; iBin < nBinsX; iBin++){
    double xCenter = histogramRebinned->GetXaxis()->GetBinCenter(iBin);
    histogramRebinned->SetBinContent(iBin, fitFunctionDrawn->Eval(xCenter)); 
    double oneSigmaInterval = 0.683;
    double errorEval[1] = {0};
    double xEval[1] = {xCenter};
    fitResult->GetConfidenceIntervals(1, 1, 1, xEval, errorEval, oneSigmaInterval, false);
    histogramRebinned->SetBinError(iBin, errorEval[0]);
  }

  std::tuple<TH1D*, TGraphErrors*, TF1*> rebinResultAndFitFunction(histogramRebinned, fitFunctionTGraphErrors, fitFunctionDrawn);
  return rebinResultAndFitFunction;
}

//////////////////////////////////////////////////////////////

#endif