#ifndef HISTOGRAM_PLOTTING_C
#define HISTOGRAM_PLOTTING_C


#include "HistogramPlotting.h"
#include "../Settings/AxisTitles.h"

#include "../Settings/GlobalSettings.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <array>
#include <vector>
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <string>
#include <filesystem>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Various Utilities /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float findMinFloat(float* array, int length){
  float min = array[0];
  for (int i = 1; i < length; i++)
      if (array[i] < min)
          min = array[i];
  return min;
}
float findMaxFloat(float* array, int length){
  float max = array[0];
  for (int i = 1; i < length; i++)
      if (array[i] > max)
          max = array[i];
  return max;
}

// Returns:
//   true upon success.
//   false upon failure, and set the std::error_code & err accordingly.
bool CreateDirectoryRecursive(std::string const & dirName, std::error_code & err) //https://stackoverflow.com/questions/71658440/c17-create-directories-automatically-given-a-file-path
{
    err.clear();
    if (!std::filesystem::create_directories(dirName, err))
    {
        if (std::filesystem::exists(dirName))
        {
            // The folder already exists:
            err.clear();
            return true;    
        }
        return false;
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Context /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextCustomFiveFields(TString mainContext, TString secondaryContext, TString tertiaryContext, TString quaternaryContext, TString quinaryContext, __attribute__ ((unused)) std::string options){
  TString texContextFinal;
  // texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{"+tertiaryContext+"}";
  texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{#splitline{"+tertiaryContext+"}{#splitline{"+quaternaryContext+"}{"+quinaryContext+"}}}";
  // texContextFinal = "testtesttestest";
  return texContextFinal;
}

TString contextCustomFourFields(TString mainContext, TString secondaryContext, TString tertiaryContext, TString quaternaryContext, __attribute__ ((unused)) std::string options){
  TString texContextFinal;
  // texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{"+tertiaryContext+"}";
  texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{#splitline{"+tertiaryContext+"}{"+quaternaryContext+"}}";
  // texContextFinal = "testtesttestest";
  return texContextFinal;
}


TString contextCustomThreeFields(TString mainContext, TString secondaryContext, TString tertiaryContext, __attribute__ ((unused)) std::string options){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{"+tertiaryContext+"}";
  // texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{#splitline{2023 QC}{test"+tertiaryContext+"}}";
  // texContextFinal = "testtesttestest";
  return texContextFinal;
}

TString contextCustomTwoFields(TString mainContext, TString secondaryContext, __attribute__ ((unused)) std::string options){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+"}{"+secondaryContext+"}";
  return texContextFinal;
  // return contextCustomThreeFields(mainContext, (TString)"" , secondaryContext, options);

}

TString contextCustomOneField(TString mainContext, __attribute__ ((unused)) std::string options){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+"}{}";
  return texContextFinal;
  // return contextCustomTwoFields(mainContext, (TString)"", options);
}

TString contextPtRange(float* PtRange){
  int lowBoundPrec = PtRange[0] < 10 ? 2 : 0;
  int highBoundPrec = PtRange[1] < 10 ? 2 : 0; 

  std::stringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(lowBoundPrec);
  ss << PtRange[0] << " < #it{p}_{T} < ";
  ss.precision(highBoundPrec);
  ss << PtRange[1];
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
  return textContext;
}

TString contextEtaRange(float* EtaRange){
  std::stringstream ss;
  ss << EtaRange[0] << " < #eta < " << EtaRange[1];
  TString textContext((TString)ss.str());
  return textContext;
}

TString contextCentRange(float* CentRange){
  std::stringstream ss;
    ss << "" << CentRange[0] << "-" << CentRange[1] << "%";
  TString textContext((TString)ss.str());
  return textContext;
}

TString contextJetRadius(float jetRadius){
  std::stringstream ss;
  ss << "#it{R} = " << jetRadius;
  TString textContext((TString)ss.str());
  // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
  return textContext;
}






TString contextDatasetRadiusCompAndVarRange(TString mainContext, int iDataset, float* variableRange, std::string options){
  TString texcontextDatasetRadiusCompAndVarRange;
  if (options.find("pt") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+mainContext+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextPtRange(variableRange)+"}}";
  }
  if (options.find("eta") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetRadiusCompAndVarRange = "#splitline{"+mainContext+" "+DatasetsNames[iDataset]+"}{#splitline{2023 QC}{"+contextEtaRange(variableRange)+"}}";
  }

  return texcontextDatasetRadiusCompAndVarRange;
}

TString contextDatasetCompAndRadiusAndVarRange(TString mainContext, float jetRadius, float* variableRange, std::string options){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  if (options.find("pt") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextPtRange(variableRange)+"}}";
  }
  if (options.find("eta") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextEtaRange(variableRange)+"}}";
  }
  if (options.find("centrality") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
    texcontextDatasetCompAndRadiusAndVarRange = "#splitline{"+mainContext+"}{#splitline{"+contextJetRadius(jetRadius)+"}{"+contextCentRange(variableRange)+"}}";
  }

  return texcontextDatasetCompAndRadiusAndVarRange;
}

TString contextDatasetCompAndRadius(TString mainContext, float jetRadius, __attribute__ ((unused)) std::string options){
  TString texcontextDatasetCompAndRadius;
  texcontextDatasetCompAndRadius = "#splitline{"+mainContext+"}{"+contextJetRadius(jetRadius)+"}";

  return texcontextDatasetCompAndRadius;
}

TString contextDatasetComp(TString mainContext, __attribute__ ((unused)) std::string options){
  TString texcontextDatasetComp;
  texcontextDatasetComp = mainContext;

  return texcontextDatasetComp;
}

void CentralityLegend(TString* centralityLegend, const float arrayCentralityIntervals[][2], int nCentralityBins){
  std::stringstream ss;
  ss.precision(2);
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    ss << "" << arrayCentralityIntervals[iCentralityBin][0] << "-" << arrayCentralityIntervals[iCentralityBin][1] << "%";
    centralityLegend[iCentralityBin] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Histogram Drawing /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_TH1_Histograms_MasterFunction(TH1D** histograms_collection, const TString* legendList_string, TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options, std::vector<TGraphErrors*> optionalFitCollection) {
  // has options:
  // - "autoratio" : if in the options string, the Y range is chosen automatically based on the difference to 1
  // - "zoomToOneLarge" : if in the options string, the Y range is [0,2.2]
  // - "zoomToOneMedium1" : if in the options string, the Y range is [0.6,1.54]
  // - "logy" : if in the options string, then the y axis of the plot is set to a log scale (except for the ratio plot)
  // - "noMarkerFirst", "noMarkerSecond"  : if in the options string, then the first/second histogram of the collection isn't plotted

  int largeCollectionThreshold = 14;
  int collectionSizeColorThreshold = 6;

  // canvas settings
  TCanvas *canvas = new TCanvas ("canvasWithRatio"+*pdfName, "canvasWithRatio"+*pdfName, 800, 800);
  TPad *padMainHist = new TPad("padMainHist","padMainHist",0.0, 0.0, 1.0, 1.0);
  TPad *padRatio    = new TPad("padRatio","padRatio",      0.0, 0.0, 1.0, 0.3); // only used if "ratioInSameCanvas" is in options
  canvas->cd(0);
  padMainHist->Draw();
  padMainHist->SetFillColor(0);
  padMainHist->SetFrameFillStyle(0);

  if (options.find("ratioInSameCanvas") != std::string::npos) {
    gPad->SetBottomMargin(0.01);
    gPad->SetTopMargin(0);
    padMainHist->SetPad(0.0, 0.3, 1.0, 1.0);
    // padMainHist->SetTopMargin(0.01);
    padMainHist->SetBottomMargin(0);
    padMainHist->SetGridx(); // Vertical grid  

    canvas->cd(0);
    padRatio->Draw();
    padRatio->SetTopMargin(0);
    padRatio->SetBottomMargin(0.3);
    padRatio->SetGridx(); // Vertical grid
    padRatio->SetFillColor(0);
    padRatio->SetFrameFillStyle(0);
  }
  padMainHist->cd();

  std::vector<float> minY_collection(collectionSize);
  std::vector<float> maxY_collection(collectionSize);
  std::vector<float> minY_ratio_collection(collectionSize);
  std::vector<float> maxY_ratio_collection(collectionSize);
  std::vector<float> minX_collection(collectionSize);
  std::vector<float> maxX_collection(collectionSize);

  for (int i = 0; i < collectionSize; i++) {
    if (options.find("autoXrange") != std::string::npos) {
      maxX_collection[i] = histograms_collection[i]->FindLastBinAbove(0, 1);
      minX_collection[i] = histograms_collection[i]->FindFirstBinAbove(0,1);
      // cout << "test1.1a" << endl;
    }
    else {
      maxX_collection[i] = histograms_collection[i]->GetXaxis()->GetXmax();
      minX_collection[i] = histograms_collection[i]->GetXaxis()->GetXmin();
      // cout << "test1.1b" << endl;
    }
    maxY_collection[i] = histograms_collection[i]->GetMaximum();
    if (options.find("ratioInSameCanvas") != std::string::npos) {
      maxY_ratio_collection[i] = histograms_collection_ratios[i]->GetMaximum();
    }
    // cout << "test1.2" << endl;
    
    if (options.find("logy") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
      minY_collection[i] = histograms_collection[i]->GetMinimum(1E-10);
      // cout << "test1.3a" << endl;
    }
    else if (options.find("minYnotZero") != std::string::npos) {
      minY_collection[i] = histograms_collection[i]->GetMinimum(GLOBAL_epsilon); // asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
      // cout << "test1.3b" << endl;
    }
    else {
      minY_collection[i] = histograms_collection[i]->GetMinimum();
      // cout << "test1.3c" << endl;
    }

    if (options.find("ratioInSameCanvas") != std::string::npos) {
      if (options.find("ratioLogy") != std::string::npos) { //  || options.find("ratio") != NULL not sure why I had this here
        minY_ratio_collection[i] = histograms_collection_ratios[i]->GetMinimum(1E-10);
        // cout << "test1.3a" << endl;
      }
      else {
        minY_ratio_collection[i] = histograms_collection_ratios[i]->GetMinimum();
        // cout << "test1.3c" << endl;
      }
    }
  }
  // cout << "test2" << endl;

  float yUpMarginScaling, yDownMarginScaling;
  float maxX, minX, maxY, minY;
  float maxY_ratio;
  float minY_ratio;

  yDownMarginScaling = 1;
  maxX = *std::max_element(maxX_collection.begin(), maxX_collection.end());
  minX = *std::min_element(minX_collection.begin(), minX_collection.end());
  maxY = *std::max_element(maxY_collection.begin(), maxY_collection.end());
  minY = *std::min_element(minY_collection.begin(), minY_collection.end());
  maxY > 0 ? yUpMarginScaling = 1.5 : yUpMarginScaling = 0.6;
  maxY_ratio = *std::max_element(maxY_ratio_collection.begin(), maxY_ratio_collection.end());
  minY_ratio = *std::min_element(minY_ratio_collection.begin(), minY_ratio_collection.end());

  if (options.find("logy") != std::string::npos) {
    yUpMarginScaling = 100;
    if (minY < 0) {
    } 
    else if (minY < 1E-10) { // if minY is 0 logY won't like it
      minY = minY+1E-10 ;
    }
    else {
      minY = minY-1E-10 ; // to be sure to see the point?
    }
  } else {
    if (minY < 0 || options.find("minYnotZero") != std::string::npos) {
      minY > 0 ? yDownMarginScaling = 0.95 : yDownMarginScaling = 1.05;
      yUpMarginScaling = 1.1;
    } else {
      minY = 0.;
    }

    if (options.find("autoratio") != std::string::npos) {
      float deltaMax = max(1-minY, maxY-1);
      minY = max((float)0., 1-deltaMax);
      maxY = 1+deltaMax;
      yUpMarginScaling = 1.3;
    }
    if (options.find("zoomToOneLarge") != std::string::npos) {
      minY = 0;
      maxY = 2;
      yUpMarginScaling = 1.1;
    }
    if (options.find("zoomToOneMedium1") != std::string::npos) {
      minY = 0.6;
      maxY = 1.4;
      yUpMarginScaling = 1.1;
    }
    if (options.find("zoomToOneMedium2") != std::string::npos) {
      minY = 0.8;
      maxY = 1.2;
      yUpMarginScaling = 1.1;
    }
    if (options.find("zoomToOneExtra") != std::string::npos) {
      minY = 0.9;
      maxY = 1.05;
      yUpMarginScaling = 1.1;
    }
    if (options.find("zoomToOneExtraExtra") != std::string::npos) {
      minY = 0.93;
      maxY = 1.;
      yUpMarginScaling = 1.1;
    }
    if (options.find("efficiency") != std::string::npos) {
      minY = 0;
      maxY = 1.5;
      yUpMarginScaling = 1.1;
    }

    // choosing y-range for second plot if ratioInSameCanvas is requested
    if (options.find("ratioInSameCanvas") != std::string::npos) {
      float deltaMax = max(1-minY, maxY-1);
      minY_ratio = max((float)0., 1-deltaMax);
      maxY_ratio = 1+deltaMax;
      if (options.find("ratioZoomToOneLarge") != std::string::npos) {
        minY_ratio = 0;
        maxY_ratio = 2;
      }
      if (options.find("ratioZoomToOneMedium1") != std::string::npos) {
        minY_ratio = 0.6;
        maxY_ratio = 1.4;
      }
      if (options.find("ratioZoomToOneMedium2") != std::string::npos) {
        minY_ratio = 0.8;
        maxY_ratio = 1.2;
      }
      if (options.find("ratioZoomToOneExtra") != std::string::npos) {
        minY_ratio = 0.9;
        maxY_ratio = 1.1;
      }
      if (options.find("ratioZoomToOneExtraExtra") != std::string::npos) {
        minY_ratio = 0.93;
        maxY_ratio = 1.07;
      }
    }
  }
    // cout << "test3" << endl;

  if (!std::equal(std::begin(drawnWindow[0]), std::end(drawnWindow[0]), std::begin(drawnWindowAuto[0]), std::end(drawnWindowAuto[0]))) {
    minX = drawnWindow[0][0];
    maxX = drawnWindow[0][1];
  }
  if (!std::equal(std::begin(drawnWindow[1]), std::end(drawnWindow[1]), std::begin(drawnWindowAuto[1]), std::end(drawnWindowAuto[1]))) {
    yUpMarginScaling = 1;
    yDownMarginScaling = 1;
    minY = drawnWindow[1][0];
    maxY = drawnWindow[1][1];
  }

  if (options.find("ratioInSameCanvas") != std::string::npos) {
    if (!std::equal(std::begin(drawnWindowRatio), std::end(drawnWindowRatio), std::begin(drawnWindowAuto[1]), std::end(drawnWindowAuto[1]))) {
      minY = drawnWindowRatio[0];
      maxY = drawnWindowRatio[1];
    }
  }

  TH1 *hFrame = canvas->DrawFrame(minX, yDownMarginScaling*minY, maxX, yUpMarginScaling*maxY);
  TH1 *hFrameRatio;
  if (options.find("ratioInSameCanvas") != std::string::npos) {
    padRatio->cd();
    hFrameRatio = padRatio->DrawFrame(minX, minY_ratio, maxX, maxY_ratio); // REALLY NOT SURE THIS IS WORKING, TO BE CHECKED
    padMainHist->cd();
  }
  if (options.find("logy") != std::string::npos) {
    padMainHist->SetLogy();
  }  
  if (options.find("logx") != std::string::npos) {
    padMainHist->SetLogx();
    if (options.find("ratioInSameCanvas") != std::string::npos) {
      // padRatio->cd();
      padRatio->SetLogx();
      // padMainHist->cd();
    }
  }
    // cout << "test4" << endl;

  hFrame->SetXTitle(texXtitle->Data());
  hFrame->SetYTitle(texYtitle->Data());
  // hFrame->GetYaxis()->SetLabelFont(63); // TESTESTEST;
  // hFrame->GetYaxis()->SetLabelSize(16); //in pixels
  if (options.find("ratioInSameCanvas") != std::string::npos) {
    hFrameRatio->SetXTitle(texXtitle->Data());
    hFrameRatio->SetYTitle(texYtitle->Data());
    hFrameRatio->GetYaxis()->SetLabelFont(63);
    hFrameRatio->GetYaxis()->SetLabelSize(16); //in pixels
  }
  // hFrame->GetYaxis()->SetTitleOffset(2);

  // legend settings
  double xLeftLegend = 0.7;
  double yLowLegend = 0.75;
  double xRightLegend = 0.8;
  double yUpLegend = 0.87;

  double xLeftLegendRatio = 0.7;
  double yLowLegendRatio = 0.75;
  double xRightLegendRatio = 0.8;
  double yUpLegendRatio = 0.87;

  if (collectionSize > largeCollectionThreshold) {
    xLeftLegend = 0.65;
    yLowLegend = 0.6;
    xRightLegend = 0.85;
    yUpLegend = 0.85;

    xLeftLegendRatio = 0.4;
    yLowLegendRatio = 0.4;
    xRightLegendRatio = 0.9;
    yUpLegendRatio = 0.7;
  }
  if (!std::equal(std::begin(legendPlacement[0]), std::end(legendPlacement[0]), std::begin(legendPlacementAuto[0]), std::end(legendPlacementAuto[0]))) {
    xLeftLegend = legendPlacement[0][0];
    yLowLegend = legendPlacement[0][1];
  }
  if (!std::equal(std::begin(legendPlacement[1]), std::end(legendPlacement[1]), std::begin(legendPlacementAuto[1]), std::end(legendPlacementAuto[1]))) {
    xRightLegend = legendPlacement[1][0];
    yUpLegend = legendPlacement[1][1];
  }
  if (!std::equal(std::begin(legendPlacementRatio[0]), std::end(legendPlacementRatio[0]), std::begin(legendPlacementAuto[0]), std::end(legendPlacementAuto[0]))) {
    xLeftLegendRatio = legendPlacementRatio[0][0];
    yLowLegendRatio = legendPlacementRatio[0][1];
  }
  if (!std::equal(std::begin(legendPlacementRatio[1]), std::end(legendPlacementRatio[1]), std::begin(legendPlacementAuto[1]), std::end(legendPlacementAuto[1]))) {
    xRightLegendRatio = legendPlacementRatio[1][0];
    yUpLegendRatio = legendPlacementRatio[1][1];
  }
  TLegend * leg = new TLegend(xLeftLegend, yLowLegend, xRightLegend, yUpLegend);
  TLegend * legRatios = new TLegend(xLeftLegendRatio, yLowLegendRatio, xRightLegendRatio, yUpLegendRatio);

  leg->SetTextSize(gStyle->GetTextSize()*0.7);
  legRatios->SetTextSize(gStyle->GetTextSize()*0.7);
  if (collectionSize >= collectionSizeColorThreshold) { // maybe fine tune that
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
    legRatios->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (options.find("fitCollection") != std::string::npos) {
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
    legRatios->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (collectionSize >= collectionSizeColorThreshold) {
    gStyle->SetPalette(kRainbow); // for the choice of marker's colours; only use this if we have many histograms in the same plot
  }
  // cout << "test5" << endl;

  // draws histograms from collection, and setting the colors, line and marker options
  int nColors = gStyle->GetNumberOfColors();
  int histoPaletteColor;
  for (int i = 0; i < collectionSize; i++) {
    if (options.find("twoByTwoDatasetPairs") == std::string::npos && options.find("twoRatioSetsInOne") == std::string::npos) { // if hasn't found "twoByTwoDatasetPairs" in options
      if (collectionSize >= collectionSizeColorThreshold) {
        histoPaletteColor = (float)nColors / collectionSize * (int)i;
        histograms_collection[i]->SetLineColor(gStyle->GetColorPalette(histoPaletteColor));
        histograms_collection[i]->SetMarkerColor(gStyle->GetColorPalette(histoPaletteColor));
        histograms_collection[i]->Draw("same"); 
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetLineColor(gStyle->GetColorPalette(histoPaletteColor));
          histograms_collection_ratios[i]->SetMarkerColor(gStyle->GetColorPalette(histoPaletteColor));
          histograms_collection_ratios[i]->Draw("same");
          padMainHist->cd();
        }

        if (options.find("histWithLine") != std::string::npos) {
          histograms_collection[i]->Draw("][ Hist same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
          if (options.find("ratioInSameCanvas") != std::string::npos) {
            padRatio->cd();
            histograms_collection_ratios[i]->Draw("][ Hist same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
            padMainHist->cd();
          }
        }
      } else {
        histograms_collection[i]->SetMarkerColor(colors[(int)i]);
        histograms_collection[i]->SetLineColor(colors[(int)i]);
        histograms_collection[i]->Draw("same");
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerColor(colors[(int)i]);
          histograms_collection_ratios[i]->SetLineColor(colors[(int)i]);
          histograms_collection_ratios[i]->Draw("same");
          padMainHist->cd();
        }

        if (options.find("histWithLine") != std::string::npos) {
          histograms_collection[i]->Draw("][ Hist same");
          if (options.find("ratioInSameCanvas") != std::string::npos) {
            padRatio->cd();
            histograms_collection_ratios[i]->Draw("][ Hist same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
            padMainHist->cd();
          }
        }
      }

      histograms_collection[i]->SetMarkerStyle(markers[i]);
      if (options.find("ratioInSameCanvas") != std::string::npos) {
        padRatio->cd();
        histograms_collection_ratios[i]->SetMarkerStyle(markers[i]);
        padMainHist->cd();
      }

      if (collectionSize > largeCollectionThreshold) {
        histograms_collection[i]->SetMarkerStyle(markers[2]);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(markers[2]);
          padMainHist->cd();
        }
      }
      if (i == 0 && options.find("noMarkerFirst") != std::string::npos) { // if i=0 and if option noMarkerFirst is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
        histograms_collection[i]->SetMarkerStyle(1);
        histograms_collection[i]->SetLineColorAlpha(kWhite, 100);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(1);
          histograms_collection_ratios[i]->SetLineColorAlpha(kWhite, 100);
          padMainHist->cd();
        }
      }
      if (i == 1 && options.find("noMarkerSecond") != std::string::npos) { // if i=1 and if option noMarkerSecond is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
        histograms_collection[i]->SetMarkerStyle(1);
        histograms_collection[i]->SetLineColorAlpha(kWhite, 100);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(1);
          histograms_collection_ratios[i]->SetLineColorAlpha(kWhite, 100);
          padMainHist->cd();
        }
      }
    } else { // if has found "twoByTwoDatasetPairs" in options
      if (collectionSize >= 2*collectionSizeColorThreshold) {
        histoPaletteColor = (float)nColors / collectionSize * (int)i/2;
        histograms_collection[i]->SetLineColor(gStyle->GetColorPalette(histoPaletteColor));
        histograms_collection[i]->SetMarkerColor(gStyle->GetColorPalette(histoPaletteColor));
        histograms_collection[i]->Draw("same"); 
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetLineColor(gStyle->GetColorPalette(histoPaletteColor));
          histograms_collection_ratios[i]->SetMarkerColor(gStyle->GetColorPalette(histoPaletteColor));
          histograms_collection_ratios[i]->Draw("same"); 
          padMainHist->cd();
        }

        if (options.find("histWithLine") != std::string::npos) {
          histograms_collection[i]->Draw("][ Hist same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
          if (options.find("ratioInSameCanvas") != std::string::npos) {
            padRatio->cd();
            histograms_collection_ratios[i]->Draw("][ Hist same"); // PMC uses the palette chosen with gStyle->SetPalette() to chose the colours of the markers, PLC for the lines
            padMainHist->cd();
          }
        }
      } else {
        histograms_collection[i]->Draw("same");
        histograms_collection[i]->SetMarkerColor(colors[(int)i/2]);
        histograms_collection[i]->SetLineColor(colors[(int)i/2]);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->Draw("same");
          histograms_collection_ratios[i]->SetMarkerColor(colors[(int)i/2]);
          histograms_collection_ratios[i]->SetLineColor(colors[(int)i/2]);
          padMainHist->cd();
        }

        if (options.find("histWithLine") != std::string::npos) {
          histograms_collection[i]->Draw("][ Hist same");
          if (options.find("ratioInSameCanvas") != std::string::npos) {
            padRatio->cd();
            histograms_collection_ratios[i]->Draw("][ Hist same");
            padMainHist->cd();
          }
        }
      }
      histograms_collection[i]->SetMarkerStyle(markerstwoByTwoDatasetPairs[i]);
      if (options.find("ratioInSameCanvas") != std::string::npos) {
        padRatio->cd();
        histograms_collection_ratios[i]->SetMarkerStyle(markerstwoByTwoDatasetPairs[i]);
        padMainHist->cd();
      }
          
      if (i == 0 && options.find("noMarkerFirst") != std::string::npos) { // if i=0 and if option noMarkerFirst is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
        histograms_collection[i]->SetMarkerStyle(1);
        histograms_collection[i]->SetLineColorAlpha(kWhite, 100);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(1);
          histograms_collection_ratios[i]->SetLineColorAlpha(kWhite, 100);
          padMainHist->cd();
        }
      }
      if (i == 1 && options.find("noMarkerSecond") != std::string::npos) { // if i=1 and if option noMarkerSecond is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
        histograms_collection[i]->SetMarkerStyle(1);
        histograms_collection[i]->SetLineColorAlpha(kWhite, 100);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(1);
          histograms_collection_ratios[i]->SetLineColorAlpha(kWhite, 100);
          padMainHist->cd();
        }
      }
    }
    if (i == 0 && options.find("noMarkerFirst") != std::string::npos) { // if i=1 and if option noMarkerSecond is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
      continue;
    }
    if (i == 1 && options.find("noMarkerSecond") != std::string::npos) { // if i=1 and if option noMarkerSecond is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
      continue;
    }
      leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
    if (options.find("ratioInSameCanvas") != std::string::npos) {
      padRatio->cd();
      legRatios->AddEntry(histograms_collection_ratios[i], legendList_string_ratios[i], "LP");
      padMainHist->cd();
    }
    if (options.find("smallMarkers") != std::string::npos) { // if i=1 and if option noMarkerSecond is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
      histograms_collection[i]->SetMarkerSize(0.3);
      // histograms_collection[i]->SetMarkerStyle(1);
    }

  }
  padMainHist->cd();

  TGaxis *axis = new TGaxis( minX, minY, maxX, maxY, minY, maxY, 510,"");
  // TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  if (options.find("ratioInSameCanvas") != std::string::npos) {
    hFrame->GetYaxis()->SetLabelSize(0.);
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();
  }

  // draws fit functions if requested
  if (options.find("fitCollection") != std::string::npos) {
    for (int i = 0; i < collectionSize; i++) {
      // optionalFitCollection[i]->SetNpx(2000);
      // optionalFitCollection[i]->Draw("E, same");
      optionalFitCollection[i]->SetLineColorAlpha(histograms_collection[i]->GetLineColor(), 0.08);
    }
  }
  if (options.find("fitSingle") != std::string::npos) {
    // optionalFitCollection[0]->SetNpx(2000);
    optionalFitCollection[0]->Draw("E, same");
    optionalFitCollection[0]->SetLineColorAlpha(colors[collectionSize], 0.08);
  }

  if (collectionSize >= 2) {
    leg->Draw("same");
    if (options.find("ratioInSameCanvas") != std::string::npos) {
      padRatio->cd();
      legRatios->Draw("same");
      padMainHist->cd();
    }
  }
  // cout << "test6" << endl;

  if (options.find("ratioLine") != std::string::npos) {
    TLine myline(minX,1,maxX,1);
    myline.SetLineColor(kBlack);
    myline.SetLineWidth(1);
    // myline.SetLineStyle(2);
    myline.DrawLine(minX,1,maxX,1);
		canvas->Modified();
		canvas->Update();
    cout << "minX = " << minX << endl;
  }

  if (options.find("150MevLine") != std::string::npos) {
    float lineEdgesX[4] = {0.150, 0.150};
    float lineEdgesY[4] = {minY, yUpMarginScaling*maxY};
    TPolyLine* Line150Mev = new TPolyLine(2, lineEdgesX, lineEdgesY);
    // cout << "Line150Mev->GetN() = " << Line150Mev->GetN() << endl;
    if (Line150Mev->GetN() > 0) {
      Line150Mev->Draw("");
      Line150Mev->SetLineColor(kBlack);
    }
  }

  if (options.find("datasetXaxisBinLabels") != std::string::npos) {
    // ChangeLabel(labNum, labAngle, labSize, labAlign, labColor, labFont, labText)
    // [in]	labNum	Number of the label to be changed, negative numbers start from the end
    // [in]	labAngle	New angle value
    // [in]	labSize	New size (0 erase the label)
    // [in]	labAlign	New alignment value
    // [in]	labColor	New label color
    // [in]	labFont	New label font
    // [in]	labText	New label text
    int nBins = histograms_collection[0]->GetNbinsX();
    for (int i = 0; i < nBins; i++) {
      hFrame->GetXaxis()->SetBinLabel(hFrame->GetNbinsX()/(2*nBins) * (2*i + 1), DatasetsNames[i]); // DrawFrame creates a histo with 1000 bins; takes bin number as input
      hFrame->GetXaxis()->ChangeLabel(1, 30, -1, -1, -1, -1, -1); // didn't manage to get it to work, but also didnt spend much time; maybe because it takes label number as input?
    }
  }

  // Context drawing 
  TLatex* textInfo = new TLatex();
  textInfo->SetTextSize(0.04);
  textInfo->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram  

  double xTopLeftCornerContext = 0.18;
  double yTopLeftCornerContext = 0.82;
  double deltaYContextsPosition = 0.07;
  if (contextPlacement[0] != contextPlacementAuto[0]) {
    xTopLeftCornerContext = contextPlacement[0];
  }
  if (contextPlacement[1] != contextPlacementAuto[1]) {
    yTopLeftCornerContext = contextPlacement[1];
  }
  textInfo->DrawLatex(xTopLeftCornerContext,yTopLeftCornerContext,texCollisionDataInfo->Data());
  textInfo->DrawLatex(xTopLeftCornerContext,yTopLeftCornerContext - deltaYContextsPosition,Context);



  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/", errEPS);
  canvas->SaveAs("pdfFolder/"+*pdfName+".pdf");
  canvas->SaveAs("pngFolder/"+*pdfName+".png");
  canvas->SaveAs("epsFolder/"+*pdfName+".eps");

  // struct stat st1{};
  // if (stat("pdfFolder/", &st1) == -1) {
  //     mkdir("pdfFolder/", 0700);
  // }
  // canvas->SaveAs("pdfFolder/"+*pdfName+".pdf");

  // struct stat st2{};
  // if (stat("pngFolder/", &st2) == -1) {
  //     mkdir("pngFolder/", 0700);
  // }
  // canvas->SaveAs("pngFolder/"+*pdfName+".png");

  // struct stat st3{};
  // if (stat("epsFolder/", &st3) == -1) {
  //   mkdir("epsFolder/", 0700);
  // }
  // canvas->SaveAs("epsFolder/"+*pdfName+".eps");

  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "histograms_collection[0]->GetBinContent(iCentralityBin) = " << histograms_collection[0]->GetBinContent(iCentralityBin) << endl;
  // }
}

// ratio and distrib in same canvas (ratio just below distrib) Work in progress
void Draw_TH1_Histograms_ratioInSameCanvas(TH1D** histograms_collection, const TString* legendList_string,TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options, std::vector<TGraphErrors*> optionalFitCollection) {
  Draw_TH1_Histograms_MasterFunction(histograms_collection, legendList_string, histograms_collection_ratios, legendList_string_ratios, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, drawnWindowRatio, legendPlacement, legendPlacementRatio, contextPlacement, options+(std::string)"ratioInSameCanvas", optionalFitCollection);
}

// ratio and distrib in same canvas (ratio just below distrib) Work in progress
void Draw_TH1_Histograms_ratioInSameCanvas(TH1D** histograms_collection, const TString* legendList_string,TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms function; it will only be used if 'options' has fit in it
  // TF1* optionalFitCollectionDummy[collectionSize];
  std::vector<TGraphErrors*> optionalFitCollectionDummy(collectionSize);
  Draw_TH1_Histograms_ratioInSameCanvas(histograms_collection, legendList_string, histograms_collection_ratios, legendList_string_ratios, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, drawnWindowRatio, legendPlacement, legendPlacementRatio, contextPlacement, options, optionalFitCollectionDummy);
}

void Draw_TH1_Histograms(TH1D** histograms_collection, const TString* legendList_string, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options, std::vector<TGraphErrors*> optionalFitCollection) {
  TH1D* histograms_collection_ratios_dummy[collectionSize];
  // std::vector<TH1D*> histograms_collection_ratios_dummy(collectionSize);
  const TString* legendList_string_ratios_dummy;
  std::array<float, 2> drawnWindowRatioDummy;
  std::array<std::array<float, 2>, 2> legendPlacementRatioDummy;
  Draw_TH1_Histograms_MasterFunction(histograms_collection, legendList_string, histograms_collection_ratios_dummy, legendList_string_ratios_dummy, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, drawnWindowRatioDummy, legendPlacement, legendPlacementRatioDummy, contextPlacement, options, optionalFitCollection);
}

void Draw_TH1_Histograms(TH1D** histograms_collection, const TString* legendList_string, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms function; it will only be used if 'options' has fit in it
  // TF1* optionalFitCollectionDummy[collectionSize];
  std::vector<TGraphErrors*> optionalFitCollectionDummy(collectionSize);
  Draw_TH1_Histograms(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options, optionalFitCollectionDummy);
}

void Draw_TH1_Histogram(TH1D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options) {
  TH1D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH1_Histograms(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options);
}
void Draw_TH1_Histogram(TH1D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options, TGraphErrors* optionalFit) {
  TH1D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  std::vector<TGraphErrors*> optionalFitCollection{optionalFit};
  Draw_TH1_Histograms(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options, optionalFitCollection);
}



void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options, TPolyLine* optionalLine) {

  double width = collectionSize*900;
  double height = 800;
  auto canvas = new TCanvas("canvas"+*pdfName, "canvas"+*pdfName, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));

  canvas->Divide(collectionSize,1);

  // TH2D* contourHist[collectionSize];
  std::vector<TH2D*> contourHist(collectionSize);

  // draws histograms from collection
  for (int i = 0; i < collectionSize; i++) {
    canvas->cd(i+1);
    histograms_collection[i]->Draw("colz");
    if (th2ContourNumber > 0){
      contourHist[i] = (TH2D*)histograms_collection[i]->Clone(legendList_string[i]);
      contourHist[i]->SetContour(th2ContourNumber, th2Contours);
      contourHist[i]->Draw("cont3 same");
    }
    histograms_collection[i]->SetXTitle(texXtitle->Data());
    histograms_collection[i]->SetYTitle(texYtitle->Data());
    canvas->cd(i+1)->SetRightMargin(0.18); // if the z-axis ever gets hidden, one can play with this
    if (options.find("logz") != std::string::npos) {
      gPad->SetLogz(); // sets log scale for the current pad
    }
    if (options.find("logy") != std::string::npos) {
      gPad->SetLogy();
    }
    if (options.find("logx") != std::string::npos) {
      gPad->SetLogx();
    }
    // leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
  }


  if (std::equal(std::begin(drawnWindow2D[0]), std::end(drawnWindow2D[0]), std::begin(drawnWindow2DAuto[0]), std::end(drawnWindow2DAuto[0]))) { // auto x axis
    if (options.find("autoRangeSame") != std::string::npos) {
      int maxXbin = 0;
      int symBinLimitMin = 9999999;
      for (int i = 0; i < collectionSize; i++) {
        if (maxXbin < histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 1)) {
          maxXbin = histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 1);// (asks for the first/last bin on the x axis (axis number 1) to have strictly more than 1 entry)
        }
      }
      for (int i = 0; i < collectionSize; i++) {
        histograms_collection[i]->GetXaxis()->SetRange(1, maxXbin);
      }
    }
  } else {
    double minX = drawnWindow2D[0][0];
    double maxX = drawnWindow2D[0][1];
    for (int i = 0; i < collectionSize; i++) {
      histograms_collection[i]->GetXaxis()->SetRangeUser(minX, maxX);
    }
  }

  if (std::equal(std::begin(drawnWindow2D[1]), std::end(drawnWindow2D[1]), std::begin(drawnWindow2DAuto[1]), std::end(drawnWindow2DAuto[1]))) { // auto y axis
    if (options.find("autoRangeSame") != std::string::npos) {
      int maxYbin = 0;
      int symBinLimitMin = 9999999;
      for (int i = 0; i < collectionSize; i++) {
        if (maxYbin < histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2)) {
          maxYbin = histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2);// (asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
        }
        if (options.find("autoRangeSameSym") != std::string::npos) {
          int symBinLimit = min(histograms_collection[i]->FindFirstBinAbove(GLOBAL_epsilon, 2), abs(histograms_collection[i]->GetNbinsY() - histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2)));
          if (symBinLimitMin > symBinLimit) {
            symBinLimitMin = symBinLimit;
          }
        }
      }
      for (int i = 0; i < collectionSize; i++) {
        if (options.find("autoRangeSameSym") != std::string::npos) {
          int symBinLimit = min(histograms_collection[i]->FindFirstBinAbove(GLOBAL_epsilon, 2), abs(histograms_collection[i]->GetNbinsY() - histograms_collection[i]->FindLastBinAbove(GLOBAL_epsilon, 2))); //(asks for the first/last bin on the y axis (axis number 2) to have strictly more than 1 entry)
          histograms_collection[i]->GetYaxis()->SetRange(symBinLimitMin, histograms_collection[i]->GetNbinsY() - symBinLimitMin); //getting symmetric window around 0 on Y axis
        }
        else {
          histograms_collection[i]->GetYaxis()->SetRange(1, maxYbin);
        }
      }
    }
  } else {
    double minY = drawnWindow2D[1][0];
    double maxY = drawnWindow2D[1][1];
    for (int i = 0; i < collectionSize; i++) {
      histograms_collection[i]->GetYaxis()->SetRangeUser(minY, maxY);
    }
  }

  if (std::equal(std::begin(drawnWindow2D[2]), std::end(drawnWindow2D[2]), std::begin(drawnWindow2DAuto[2]), std::end(drawnWindow2DAuto[2]))) { // auto z axis
    for (int i = 0; i < collectionSize; i++) {
      histograms_collection[i]->GetZaxis()->SetRangeUser(histograms_collection[i]->GetMinimum(GLOBAL_epsilon), histograms_collection[i]->GetMaximum());
    }
  } else {
    double minZ = drawnWindow2D[2][0];
    double maxZ = drawnWindow2D[2][1];
    for (int i = 0; i < collectionSize; i++) {
      cout << "kjsfkdsjhglkfdjhgkdjhgdkjhgkdfjhgdkjfhkfdjhgkdfjhgkdjhgkdfjhgkdfhjg" << minZ << ", " << maxZ << endl;
      histograms_collection[i]->GetZaxis()->SetRangeUser(minZ, maxZ);
    }
  }

  if (options.find("drawLines") != std::string::npos) {
    // cout << "optionalLine->GetN() = " << optionalLine->GetN() << endl;
    if (optionalLine->GetN() > 0) {
      optionalLine->Draw("");
      optionalLine->SetLineColor(kRed);
    }
  }

  gStyle->SetPalette(kBird); // a better palette than the kRainbow that was used by default; https://root.cern.ch/doc/master/classTColor.html lists it as one of the better palettes for Colour Vision Deficiencies 
  gStyle->SetNumberContours(100);
  
  // // adds some text on the plot
  TLatex* textInfo = new TLatex();
  textInfo->SetTextSize(0.04);
  textInfo->SetNDC(kTRUE); //remove if I want x,y in TLatex to be in the coordinate system of the histogram
  for (int i = 0; i < collectionSize; i++) {
    canvas->cd(i+1);
    textInfo->DrawLatex(0.18,0.82,texCollisionDataInfo->Data());
    textInfo->DrawLatex(0.18,0.75,Context);
    textInfo->DrawLatex(0.18,0.65,legendList_string[i]);
  }



  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/", errEPS);
  canvas->SaveAs("pdfFolder/"+*pdfName+".pdf");
  canvas->SaveAs("pngFolder/"+*pdfName+".png");
  canvas->SaveAs("epsFolder/"+*pdfName+".eps");

  // struct stat st1{};
  // if (stat("pdfFolder/", &st1) == -1) {
  //     mkdir("pdfFolder/", 0700);
  // }
  // canvas->SaveAs("pdfFolder/"+*pdfName+".pdf");

  // struct stat st2{};
  // if (stat("pngFolder/", &st2) == -1) {
  //     mkdir("pngFolder/", 0700);
  // }
  // canvas->SaveAs("pngFolder/"+*pdfName+".png");

  // struct stat st3{};
  // if (stat("epsFolder/", &st3) == -1) {
  //   mkdir("epsFolder/", 0700);
  // }
  // canvas->SaveAs("epsFolder/"+*pdfName+".eps");
}

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, const int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms function; it will only be used if 'options' has fit in it
  TPolyLine* optionalLine;
  Draw_TH2_Histograms(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow2D, th2Contours, th2ContourNumber, options, optionalLine);
}

void Draw_TH2_Histogram(TH2D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options) {
  TH2D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH2_Histograms(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow2D, th2Contours, th2ContourNumber, options);


  // for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
  //   cout << "histogram->GetBinContent(iCentralityBin) = " << histogram->GetBinContent(iCentralityBin) << endl;
  // }
}


//=============================================================================
//  Generic spectrum-comparison plotter with ratio panel.
//
//  Each entry = one curve. Two render styles:
//    kData   -> markers + stat error bars (TH1 "E1 P")
//    kTheory -> solid/dashed line "hist" (e.g. POWHEG central)
//  Optional per-entry systematic band (TGraphAsymmErrors box, "2"/"E2").
//  Optional per-entry vdM/lumi band tracking the curve.
//
//  Ratio panel: ratio[i] = entry[i] / entry[refIndex], for every i != refIndex.
//
//  *** DECORRELATED systematics in the ratio ***
//   - the REFERENCE's own relative systematic is drawn as a band around UNITY
//   - each OTHER entry's own relative systematic is drawn around its ratio point
//   - the two are NOT combined in quadrature: overlap of a ratio-point band with
//     the unity band shows consistency within independent systematics.
//
//  *** Global normalization box ***
//   - pass globalNormRel >= 0 to draw a custom full-width TBox around 1 in the
//     ratio pad of half-height globalNormRel (e.g. 0.05 for a 5% global-norm unc).
//=============================================================================



enum EDrawStyle { kData, kTheory };
enum ESysKind   { kSysNone, kSysRelative, kSysAbsolute }; // how hSys is stored

struct SpecEntry {
  // --- central spectrum ---
  TH1*        h        = nullptr;        // bin contents = values, bin errors = STAT
  EDrawStyle  style    = kData;

  // --- systematic (optional) ---
  ESysKind    sysKind  = kSysNone;
  TH1*        hSys      = nullptr;        // kSysRelative: rel unc per bin
                                          // kSysAbsolute: absolute unc per bin
  // --- lumi / vdM band (optional) ---
  double      vdmRel    = 0.;             // 0 = no vdM band

  // --- cosmetics ---
  int         color    = kBlack;
  int         sysColor  = kGray+1;
  int         marker    = 20;
  double      msize     = 1.1;
  int         lstyle    = 1;             // line style for theory curves
  int         lwidth     = 3;
  TString     label     = "";            // legend (central)
  TString     sysLabel   = "";           // legend (sys band), "" = skip
  TString     legOpt     = "";           // override legend marker option
};

void DrawSpectraWithRatio(
    std::vector<SpecEntry> entries,
    int          refIndex,                 // denominator for the ratio panel
    double       xMin, double xMax,
    TString      yTitle,
    TString      xTitle      = "#it{p}_{T} (GeV/#it{c})",
    TString      ratioTitle   = "Ratio",
    std::vector<TString> labels = {},      // extra TLatex context lines (top pad)
    TString      canvasName     = "cSpecRatio",
    double       ratioMin       = 0.0,
    double       ratioMax       = 2.5,
    double       split          = 0.30,    // bottom-pad fraction
    double       vdmRatioRel    = -1.,     // >=0 -> combined vdM TBox at unity
    double       globalNormRel  = -1.,     // >=0 -> global-norm TBox at unity
    TString      globalNormLabel = "",     // legend text for the global-norm box
    double       yMinUser       = 0.,      // 0 -> auto
    double       yMaxUser       = 0.,      // 0 -> auto
    double       vdmBoxRel      = -1.,     // >=0 -> small localized vdM box at unity
    double       vdmBoxXpos     = 0.)      // pT position for the box (0 -> auto)
{
  if (entries.empty() || refIndex < 0 || refIndex >= (int)entries.size()) {
    std::cerr << "DrawSpectraWithRatio: bad entries / refIndex" << std::endl; return;
  }

  // ---- global style ----
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "X");
  gStyle->SetTickLength(0.02, "Y");
  gStyle->SetEndErrorSize(0);

  //---------------------------------------------------------------------------
  // helper: relative systematic of an entry at pT = xx (value yy)
  //---------------------------------------------------------------------------
  auto relSysAt = [](const SpecEntry& s, double xx, double yy) -> double {
    if (s.sysKind == kSysNone || !s.hSys) return 0.;
    double v = s.hSys->GetBinContent(s.hSys->FindBin(xx));
    return (s.sysKind == kSysRelative) ? v : (yy > 0 ? v / yy : 0.);
  };

  //---------------------------------------------------------------------------
  // helper: build an absolute sys-band (around the spectrum) from an entry
  //---------------------------------------------------------------------------
  auto MakeSysBand = [&](const SpecEntry& e) -> TGraphAsymmErrors* {
    if (e.sysKind == kSysNone || !e.hSys) return nullptr;
    auto* g = new TGraphAsymmErrors();
    int k = 0;
    for (int ib = 1; ib <= e.h->GetNbinsX(); ++ib) {
      double x = e.h->GetBinCenter(ib);
      double y = e.h->GetBinContent(ib);
      if (y <= 0) continue;
      double exLo = x - e.h->GetBinLowEdge(ib);
      double exHi = e.h->GetBinLowEdge(ib) + e.h->GetBinWidth(ib) - x;
      double sAbs = relSysAt(e, x, y) * y;
      g->SetPoint(k, x, y);
      g->SetPointError(k, exLo, exHi, sAbs, sAbs);
      ++k;
    }
    g->SetFillColorAlpha(e.sysColor, 0.45);
    g->SetLineColor(e.sysColor);
    g->SetFillStyle(1001);
    g->SetMarkerSize(0);
    return g;
  };

  //---------------------------------------------------------------------------
  // helper: vdM band tracking a curve
  //---------------------------------------------------------------------------
  auto MakeVdmBand = [](const SpecEntry& e) -> TGraphErrors* {
    if (e.vdmRel <= 0.) return nullptr;
    auto* g = new TGraphErrors(e.h->GetNbinsX());
    for (int ib = 1; ib <= e.h->GetNbinsX(); ++ib) {
      double x = e.h->GetBinCenter(ib);
      double y = e.h->GetBinContent(ib);
      g->SetPoint(ib-1, x, y);
      g->SetPointError(ib-1, e.h->GetBinWidth(ib)/2., y*e.vdmRel);
    }
    g->SetFillColorAlpha(e.color, 0.18);
    g->SetLineColorAlpha(e.color, 0.55);
    g->SetMarkerSize(0);
    return g;
  };

  //---------------------------------------------------------------------------
  // Canvas + pads
  //---------------------------------------------------------------------------
  TCanvas* c = new TCanvas(canvasName, canvasName, 800, 900);
  c->SetFillStyle(0);

  TPad* pTop = new TPad("pTop", "pTop", 0., split, 1., 1.);
  pTop->SetLeftMargin(0.14);  pTop->SetRightMargin(0.05);
  pTop->SetTopMargin(0.05);   pTop->SetBottomMargin(0.015);
  pTop->SetLogy();  pTop->SetTickx(1); pTop->SetTicky(1);
  pTop->Draw();

  TPad* pBot = new TPad("pBot", "pBot", 0., 0., 1., split);
  pBot->SetLeftMargin(0.14);  pBot->SetRightMargin(0.05);
  pBot->SetTopMargin(0.015);  pBot->SetBottomMargin(0.32);
  pBot->SetTickx(1); pBot->SetTicky(1);
  pBot->Draw();

  const double sf = (1. - split) / split;  // font scale for bottom pad

  //---------------------------------------------------------------------------
  // Style entries + auto y-range
  //---------------------------------------------------------------------------
  std::vector<TGraphAsymmErrors*> sysBands(entries.size(), nullptr);
  std::vector<TGraphErrors*>      vdmBands(entries.size(), nullptr);
  double yMax = -1e30, yMin = 1e30;

  for (size_t i = 0; i < entries.size(); ++i) {
    auto& e = entries[i];
    e.h->SetStats(0);
    e.h->SetLineColor(e.color);
    e.h->SetLineWidth(e.lwidth);
    if (e.style == kData) {
      e.h->SetMarkerStyle(e.marker);
      e.h->SetMarkerColor(e.color);
      e.h->SetMarkerSize(e.msize);
    } else {
      e.h->SetLineStyle(e.lstyle);
      e.h->SetMarkerSize(0);
    }
    sysBands[i] = MakeSysBand(e);
    vdmBands[i] = MakeVdmBand(e);

    yMax = std::max(yMax, e.h->GetMaximum());
    double mn = e.h->GetMinimum(0.);
    if (mn > 0) yMin = std::min(yMin, mn);
  }
  if (yMaxUser > 0) yMax = yMaxUser; else yMax *= 8.;
  if (yMinUser > 0) yMin = yMinUser; else yMin = (yMin > 0 ? yMin*0.2 : 1e-9);

  //---------------------------------------------------------------------------
  // TOP PAD
  //---------------------------------------------------------------------------
  pTop->cd();
  TH1* hFrame = entries[refIndex].h;       // reference defines the frame
  hFrame->GetYaxis()->SetTitle(yTitle);
  hFrame->GetYaxis()->SetTitleSize(0.058);
  hFrame->GetYaxis()->SetTitleOffset(1.15);
  hFrame->GetYaxis()->SetLabelSize(0.050);
  hFrame->GetXaxis()->SetLabelSize(0.);
  hFrame->GetXaxis()->SetTitleSize(0.);
  hFrame->GetXaxis()->SetLimits(xMin, xMax);
  hFrame->SetMinimum(yMin);
  hFrame->SetMaximum(yMax);

  hFrame->Draw(entries[refIndex].style == kData ? "E1 P" : "hist");

  // bands first (vdM outermost, then sys), then curves on top
  for (auto* g : vdmBands) if (g) g->Draw("E2 SAME");
  for (auto* g : sysBands) if (g) g->Draw("2 SAME");
  for (auto& e : entries)
    e.h->Draw(e.style == kData ? "E1 P SAME" : "hist SAME");

  // Legend
  TLegend* leg = new TLegend(0.50, 0.55, 0.92, 0.92);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(42);  leg->SetTextSize(0.034);
  for (size_t i = 0; i < entries.size(); ++i) {
    auto& e = entries[i];
    TString opt = e.legOpt.IsNull()
                  ? TString(e.style == kData ? "lep" : "l") : e.legOpt;
    if (!e.label.IsNull()) leg->AddEntry(e.h, e.label, opt);
    if (sysBands[i] && !e.sysLabel.IsNull())
      leg->AddEntry(sysBands[i], e.sysLabel, "f");
  }
  leg->Draw();

  // Context labels
  TLatex tex; tex.SetNDC(); tex.SetTextFont(42); tex.SetTextSize(0.040);
  double y0 = 0.90;
  for (auto& s : labels) { tex.DrawLatex(0.18, y0, s); y0 -= 0.05; }

  pTop->RedrawAxis();

  //---------------------------------------------------------------------------
  // BOTTOM PAD : ratios vs entries[refIndex]
  //---------------------------------------------------------------------------
  pBot->cd();
  TH1* hRef = entries[refIndex].h;

  // --- reference's OWN systematic, drawn as a band around UNITY ---
  TGraphAsymmErrors* gRefSys = nullptr;
  {
    const SpecEntry& r = entries[refIndex];
    if (r.sysKind != kSysNone && r.hSys) {
      gRefSys = new TGraphAsymmErrors();
      int k = 0;
      for (int ib = 1; ib <= r.h->GetNbinsX(); ++ib) {
        double x = r.h->GetBinCenter(ib);
        double y = r.h->GetBinContent(ib);
        if (y <= 0) continue;
        double rel = relSysAt(r, x, y);
        double exLo = x - r.h->GetBinLowEdge(ib);
        double exHi = r.h->GetBinLowEdge(ib) + r.h->GetBinWidth(ib) - x;
        gRefSys->SetPoint(k, x, 1.0);
        gRefSys->SetPointError(k, exLo, exHi, rel, rel);
        ++k;
      }
      gRefSys->SetFillColorAlpha(r.sysColor, 0.45);
      gRefSys->SetLineColor(r.sysColor);
      gRefSys->SetFillStyle(1001);
      gRefSys->SetMarkerSize(0);
    }
  }

  // --- reference's STATISTICAL uncertainty, as bars around UNITY ---
  TH1* hRefStat = (TH1*)hRef->Clone("hRefStat_atUnity");
  hRefStat->SetDirectory(0); hRefStat->Reset("ICES");
  for (int ib = 1; ib <= hRef->GetNbinsX(); ++ib) {
    double y = hRef->GetBinContent(ib);
    double e = hRef->GetBinError(ib);
    if (y <= 0) continue;
    hRefStat->SetBinContent(ib, 1.0);
    hRefStat->SetBinError(ib, e / y);      // relative stat error around 1
  }
  hRefStat->SetMarkerStyle(entries[refIndex].marker);
  hRefStat->SetMarkerColor(entries[refIndex].color);
  hRefStat->SetMarkerSize(entries[refIndex].msize);
  hRefStat->SetLineColor(entries[refIndex].color);

  // --- combined vdM TBox at unity (linear scale) ---
  TBox* vdmBox = nullptr;
//   if (vdmRatioRel >= 0.) {
//     vdmBox = new TBox(xMin, 1.-vdmRatioRel, xMax, 1.+vdmRatioRel);
//     vdmBox->SetFillColorAlpha(kOrange+1, 0.45);
//     vdmBox->SetLineColorAlpha(kOrange+1, 0.70);
//   }

  // --- custom global-normalization TBox at unity ---
  TBox* normBox = nullptr;
//   if (globalNormRel >= 0.) {
//     normBox = new TBox(xMin, 1.-globalNormRel, xMax, 1.+globalNormRel);
//     normBox->SetFillColorAlpha(kGreen+2, 0.25);
//     normBox->SetLineColorAlpha(kGreen+2, 0.60);
//     normBox->SetLineWidth(1);
//   }

  TLine* unity = new TLine(xMin, 1., xMax, 1.);
  unity->SetLineStyle(2); unity->SetLineColor(kBlack);

  TLegend* legR = new TLegend(0.50, 0.74, 0.92, 0.98);
  legR->SetBorderSize(0); legR->SetFillStyle(0);
  legR->SetTextFont(42);  legR->SetTextSize(0.040*sf);

  bool firstRatio = true;
  for (size_t i = 0; i < entries.size(); ++i) {
    if ((int)i == refIndex) continue;
    auto& e = entries[i];

    // ratio histogram (matched by pT center)
    TH1* hR = (TH1*)e.h->Clone(Form("hRatio_%zu", i));
    hR->SetDirectory(0); hR->Reset("ICES");
    TGraphAsymmErrors* gRsys = new TGraphAsymmErrors();
    int kr = 0;

    for (int ib = 1; ib <= e.h->GetNbinsX(); ++ib) {
      double x  = e.h->GetBinCenter(ib);
      double yN = e.h->GetBinContent(ib);
      int    jb = hRef->FindBin(x);
      double yD = hRef->GetBinContent(jb);
      if (yN <= 0 || yD <= 0) continue;
      double r = yN / yD;
      hR->SetBinContent(ib, r);

      // STAT in quadrature
      double eN = e.h->GetBinError(ib), eD = hRef->GetBinError(jb);
      hR->SetBinError(ib, r*TMath::Sqrt((eN>0?eN*eN/(yN*yN):0)+(eD>0?eD*eD/(yD*yD):0)));

      // SYS: this entry's OWN relative sys only (NOT combined with reference)
      double relS = relSysAt(e, x, yN);
      if (relS > 0) {
        double exLo = x - e.h->GetBinLowEdge(ib);
        double exHi = e.h->GetBinLowEdge(ib)+e.h->GetBinWidth(ib)-x;
        gRsys->SetPoint(kr, x, r);
        gRsys->SetPointError(kr, exLo, exHi, r*relS, r*relS);
        ++kr;
      }
    }

    if (e.style == kData) {
      hR->SetMarkerStyle(e.marker); hR->SetMarkerColor(e.color);
      hR->SetMarkerSize(e.msize);   hR->SetLineColor(e.color);
    } else {
      hR->SetLineStyle(e.lstyle); hR->SetLineColor(e.color); hR->SetMarkerSize(0);
    }

    gRsys->SetFillColorAlpha(e.sysColor, 0.45);
    gRsys->SetLineColor(e.sysColor);
    gRsys->SetFillStyle(1001); gRsys->SetMarkerSize(0);

    if (firstRatio) {
      // use hRefStat as the frame so reference stat bars sit at unity
      hRefStat->GetXaxis()->SetLimits(xMin, xMax);
      hRefStat->GetXaxis()->SetTitle(xTitle);
      hRefStat->GetYaxis()->SetTitle(ratioTitle);
      hRefStat->GetYaxis()->CenterTitle();
      hRefStat->GetYaxis()->SetNdivisions(505);
      hRefStat->GetXaxis()->SetTitleSize(0.058*sf); hRefStat->GetXaxis()->SetLabelSize(0.050*sf);
      hRefStat->GetXaxis()->SetTitleOffset(0.95);
      hRefStat->GetYaxis()->SetTitleSize(0.048*sf); hRefStat->GetYaxis()->SetLabelSize(0.044*sf);
      hRefStat->GetYaxis()->SetTitleOffset(0.55);
      hRefStat->GetXaxis()->SetTickLength(0.05);
      hRefStat->SetMinimum(ratioMin); hRefStat->SetMaximum(ratioMax);

      hRefStat->Draw("E1 P");                  // frame + reference stat bars at unity
      if (normBox) normBox->Draw("SAME");
      if (vdmBox)  vdmBox->Draw("SAME");
      if (gRefSys) gRefSys->Draw("2 SAME");    // reference sys band around unity
      hR->Draw(e.style==kData ? "E1 P SAME" : "hist SAME");
    } else {
      hR->Draw(e.style==kData ? "E1 P SAME" : "hist SAME");
    }

    if (gRsys->GetN() > 0) gRsys->Draw("2 SAME");        // this entry's own sys
    hR->Draw(e.style==kData ? "E1 P SAME" : "hist SAME"); // points/line on top

    // ratio-pad legend entries
    TString rOpt = (e.style==kData) ? "lep" : "l";
    legR->AddEntry(hR, e.label, rOpt);
    firstRatio = false;
  }

  unity->Draw("SAME");

  // --- small localized vdM/lumi box at unity ---
  TBox* vdmSmall = nullptr;
  if (vdmBoxRel >= 0.) {
    // box width = ~4% of the x-range, parked near the left edge by default
    double w  = 0.02 * (xMax - xMin);
    double x0 = (vdmBoxXpos > 0.) ? vdmBoxXpos : xMin + 0.02*(xMax - xMin);
    vdmSmall = new TBox(x0, 1.-vdmBoxRel, x0 + w, 1.+vdmBoxRel);
    vdmSmall->SetFillColorAlpha(kRed+2, 0.65);   // small DARK box
    vdmSmall->SetLineColor(kRed+2);
    vdmSmall->SetLineWidth(1);
    vdmSmall->Draw("SAME");
  }

   if (vdmSmall) {
    TLegend* legLumi = new TLegend(0.16, 0.86, 0.42, 0.92);
    legLumi->SetBorderSize(0);
    legLumi->SetFillStyle(0);
    legLumi->SetTextFont(42);
    legLumi->SetTextSize(0.034*sf);     // smaller than the 0.040*sf default
    legLumi->SetMargin(0.18);            // shrink swatch column (default ~0.25)
    legLumi->AddEntry(vdmSmall, Form("Lumi. #pm%.1f%%", vdmBoxRel*100.), "f");
    legLumi->Draw();
  }

  // reference + box legend entries
//   if (gRefSys && !entries[refIndex].sysLabel.IsNull())
//     legR->AddEntry(gRefSys, Form("%s sys. (ref.)", entries[refIndex].label.Data()), "f");
//   if (normBox)
//     legR->AddEntry(normBox,
//         globalNormLabel.IsNull()
//           ? Form("Global norm. (#pm%.1f%%)", globalNormRel*100.)
//           : globalNormLabel, "f");
//   legR->Draw();

  pBot->RedrawAxis();
  c->Update();
}

#endif