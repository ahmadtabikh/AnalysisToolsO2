#ifndef HISTOGRAM_PLOTTING_H
#define HISTOGRAM_PLOTTING_H

float findMinFloat(float* array, int length);
float findMaxFloat(float* array, int length);
bool CreateDirectoryRecursive(std::string const & dirName, std::error_code & err);

TString contextCustomFourFields(TString mainContext, TString secondaryContext, TString tertiaryContext, TString quaternaryContext, std::string options);
TString contextCustomThreeFields(TString mainContext, TString secondaryContext, TString tertiaryContext, std::string options);
TString contextCustomTwoFields(TString mainContext, TString secondaryContext, std::string options);
TString contextCustomOneField(TString mainContext, std::string options);
TString contextEtaRange(float* EtaRange);
TString contextPtRange(float* PtRange);
TString contextJetRadius(float jetRadius);

TString contextDatasetRadiusCompAndVarRange(TString mainContext, int iDataset, float* variableRange, std::string options);
TString contextDatasetCompAndRadiusAndVarRange(TString mainContext, float jetRadius, float* variableRange, std::string options);
TString contextDatasetCompAndRadius(TString mainContext, float jetRadius, std::string options);
TString contextDatasetComp(TString mainContext, std::string options);

void CentralityLegend(TString* centralityLegend, const float** arrayCentralityIntervals, int nCentralityBins);

void Draw_TH1_Histograms_MasterFunction(TH1D** histograms_collection, const TString* legendList_string, TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection);

void Draw_TH1_Histograms_ratioInSameCanvas(TH1D** histograms_collection, const TString* legendList_string,TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection);
void Draw_TH1_Histograms_ratioInSameCanvas(TH1D** histograms_collection, const TString* legendList_string,TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options);

void Draw_TH1_Histograms(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection);
void Draw_TH1_Histograms(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options);
void Draw_TH1_Histogram(TH1D* H1D_Sigma_asFunctionOf_Centrality, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options);

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, const TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options);
void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options, TPolyLine** optionalLine) ;
void Draw_TH2_Histogram(TH2D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options);
void Plot_2D_Ratio(TH2D** histos,std::string options, int nDatasets, int refIndex);
// void customRebin2D(TH2D* h2d, const char* histoNameSuffix);

// Preferred colors and markers
// const int fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7,}; // for syst bands
const int colors[]               = {kRed+1, kBlack, kBlue+1, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2, kGray+1};
const int markers[]              = {kFullSquare, kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kFullCrossX, kFullStar, kFullFourTrianglesX, kFullDoubleDiamond, kFullDiamond, kFullCross, kFullStar, kFullDoubleDiamond};
const int markerstwoByTwoDatasetPairs[]    = {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullTriangleUp, kOpenTriangleUp, kFullTriangleDown, kOpenTriangleDown, kFullDiamond, kOpenDiamond, kFullCross, kOpenCross, kFullCrossX, kOpenCrossX, kFullStar, kOpenStar, kFullFourTrianglesX, kOpenFourTrianglesX, kFullDoubleDiamond, kOpenDoubleDiamond, kFullDiamond, kOpenDiamond, kFullCross, kOpenCross, kFullStar, kOpenStar, kFullDoubleDiamond, kOpenDoubleDiamond};

#endif