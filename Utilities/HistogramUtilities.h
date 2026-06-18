#ifndef HISTOGRAM_UTILITIES_H
#define HISTOGRAM_UTILITIES_H

#include "TPolyLine.h"


std::vector<double> MakeVariableBinning_twoWidths(double xMin, int nLeft, double xMiddle, double xMax, int nRight);
std::vector<double> MakeConstantSizeBinning(double xMin, double xMax, int nBins);
std::vector<double> GetTH1Bins(TH1 H1_histo);
TH1D* TGraphErrorsToTH1D(TGraphErrors* tgraphToConvert);

TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY);
TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum);
TH2D RebinVariableBins2D(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, bool debug);
TH2D RebinVariableBins2D_PriorWeightedBinMerging(TH2D* H2D_hist, int nBinsX, int nBinsY, double* binsX, double* binsY, TH1D* H1D_priorSpectrum, bool debug);
void NormaliseYSlicesToOne(TH2D* H2D_hist);
void NormaliseXSlicesToOne(TH2D* H2D_hist);
void NormaliseXSlicesToOneNoUnderOverFlows(TH2D* H2D_hist);
void WeightMatrixWithPrior(TH2D* H2D_hist, TH1D* priorSpectrum, bool doPriorDivision);
void TransformRawResponseToYieldResponse(TH2D* H2D_hist);


TH2* NormalizeResponsMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior);
TH2D* RebinVariableBins2D_aliPhysics(TH2D* hRMFine, int nBinsX, int nBinsY, double* binsX, double* binsY, bool useFunctionWeight);

// weird stuff happening here with optional debug thing, maybe remove it completely
TH2D GetTransposeHistogram(TH2D* inputHist);
TH2D GetMatrixProductTH2xTH2(TH2D* histA, TH2D* histB);
TH1D GetMatrixVectorProductTH2xTH1(TH2D* histA, TH1D* histU);

bool DivideWithCorrelatedErrors_simpleMax(TH1D* histNumeratorToBeDivided, TH1D* histDenominator);

void Save_PtResponseMatrix(TH2D* &H2D_jetPtResponseMatrix, TString fileName);

TH1D* ReweightedRebin(TH1D* h_in, const char* newname, int n_newbins, const double* newbins);

bool HarmonizeBinning(std::vector<TH1*>& hists, double xLo, double xHi);

#endif
