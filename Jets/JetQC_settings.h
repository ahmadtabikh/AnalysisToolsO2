
// // Bin edge control
// float GLOBAL_epsilon = 0.00001;

// Analysis settings
const int nJetFinderQaType = 3;
const TString jetFinderQaHistType[nJetFinderQaType] = {"", "_rhoareasubtracted", "_eventwiseconstituentsubtracted"};
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
// const int nRadius = 3;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
// float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// const float areaDisplayMax[nRadius] = {0.5, 1, 1.5};
const int nRadius = 1;
const TString RadiusLegend[nRadius] = {"R = 0.4"};
float arrayRadius[nRadius] = {0.4};
const float areaDisplayMax[nRadius] = {0.5};
// const int nRadius = 9;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
// double arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
// const float areaDisplayMax[nRadius] = {0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.5};
const int nMethodRC = 5;
const TString methodRandomConeHistNames[nMethodRC] = {"", "withoutleadingjet", "randomtrackdirection", "randomtrackdirectionwithoutoneleadingjets", "randomtrackdirectionwithouttwoleadingjets"};
TString methodRandomConeHistLegend[nMethodRC] = {"Random Cones (RC)", "RC w/o leadJet", "RC rand(#eta,#phi)", "rand(#eta,#phi) w/o leadJet", "rand(#eta,#phi) w/o 2leadJet"};

// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;

// Choice of jet QA type (uncorrected jets, background corrected jet (rho area version), event wise constituent subtraction
const int iJetFinderQaType = 0;

// Choice of Random Cone method:
const int iMethodRandomCone = 1; 
// Default window for random cone:
std::array<std::array<float, 2>, 2> drawnWindowRCdefault = {{{-30, 60}, {5E-7, 20}}}; // {{xmin, xmax}, {ymin, ymax}}



// const int nCentralityBins = 6;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 50, 90};
// const int nCentralityBins = 3;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 50, 80};
const int nCentralityBins = 4;
const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 30, 50, 70};
// const int nCentralityBins = 7;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 60, 70, 80};

 
const int nTracksBins = 11;
const float arrayNTracksBinning[nTracksBins+1] = {0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 3000};