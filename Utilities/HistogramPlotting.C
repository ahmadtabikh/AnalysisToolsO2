
#ifndef HISTOGRAM_PLOTTING_C
#define HISTOGRAM_PLOTTING_C


#include "HistogramPlotting.h"
#include "../Settings/GlobalSettings.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <array>
#include "TGaxis.h"
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

TString contextCustomTwoFields(TString mainContext, TString secondaryContext, std::string options){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+" "+secondaryContext+"}{}";
  return texContextFinal;
  // return contextCustomThreeFields(mainContext, (TString)"" , secondaryContext, options);

}

TString contextCustomOneField(TString mainContext, std::string options){
  TString texContextFinal;
  texContextFinal = "#splitline{"+mainContext+"}{}";
  return texContextFinal;
  // return contextCustomTwoFields(mainContext, (TString)"", options);
}

TString contextPtRange(float* PtRange){
  std::stringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(2);
  ss << PtRange[0] << " < #it{p}_{T} < ";
  ss.precision(0);
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
  ss << " R = " << jetRadius;
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


void Draw_TH1_Histograms_MasterFunction(TH1D** histograms_collection, const TString* legendList_string, TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection) {
  // has options:
  // - "autoratio" : if in the options string, the Y range is chosen automatically based on the difference to 1
  // - "zoomToOneLarge" : if in the options string, the Y range is [0,2.2]
  // - "zoomToOneMedium1" : if in the options string, the Y range is [0.6,1.54]
  // - "logy" : if in the options string, then the y axis of the plot is set to a log scale (except for the ratio plot)
  // - "noMarkerFirst" : if in the options string, then the first histogram of the collection isn't plotted

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

  //////////////////////////////////loop over all histograms to determine their plotting ranges//////////////////////////////
  float minY_collection[collectionSize];
  float maxY_collection[collectionSize];
  float minY_ratio_collection[collectionSize];
  float maxY_ratio_collection[collectionSize];
  float minX_collection[collectionSize];
  float maxX_collection[collectionSize];

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
  maxX = findMaxFloat(maxX_collection, collectionSize);
  minX = findMinFloat(minX_collection, collectionSize);
  maxY = findMaxFloat(maxY_collection, collectionSize);
  minY = findMinFloat(minY_collection, collectionSize);
  maxY > 0 ? yUpMarginScaling = 1.4 : yUpMarginScaling = 0.6;
  maxY_ratio = findMaxFloat(maxY_ratio_collection, collectionSize);
  minY_ratio = findMinFloat(minY_ratio_collection, collectionSize);

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
  if (options.find("fit") != std::string::npos) {
    leg->SetTextSize(gStyle->GetTextSize()*0.3);
    legRatios->SetTextSize(gStyle->GetTextSize()*0.3);
  }
  if (collectionSize >= collectionSizeColorThreshold) {
    gStyle->SetPalette(kRainbow); // for the choice of marker's colours; only use this if we have many histograms in the same plot
  }
  // cout << "test5" << endl;

  // draws histograms from collection, and setting the colors
  int nColors = gStyle->GetNumberOfColors();
  int histoPaletteColor;
  for (int i = 0; i < collectionSize; i++) {
    if (options.find("twoByTwoDatasetPairs") == std::string::npos) { // if hasn't found "twoByTwoDatasetPairs" in options
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
      if (i == 0 && options.find("noMarkerFirst") != std::string::npos) { // if i=0 requires and if option noMarkerFirst is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
        histograms_collection[i]->SetMarkerStyle(1);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(1);
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
          
      if (i == 0 && options.find("noMarkerFirst") != std::string::npos) { // if i=0 requires and if option noMarkerFirst is there (!= std::string::npos means it found find it in the elements 0 to npos-1, where npos is the size of the string options)
        histograms_collection[i]->SetMarkerStyle(1);
        if (options.find("ratioInSameCanvas") != std::string::npos) {
          padRatio->cd();
          histograms_collection_ratios[i]->SetMarkerStyle(1);
          padMainHist->cd();
        }
      }
    }
    leg->AddEntry(histograms_collection[i], legendList_string[i], "LP");
    if (options.find("ratioInSameCanvas") != std::string::npos) {
      padRatio->cd();
      legRatios->AddEntry(histograms_collection_ratios[i], legendList_string_ratios[i], "LP");
      padMainHist->cd();
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
  if (options.find("fit") != std::string::npos) {
    for (int i = 0; i < collectionSize; i++) {
      optionalFitCollection[i]->SetNpx(2000);
      optionalFitCollection[i]->Draw("same");
      optionalFitCollection[i]->SetLineColor(histograms_collection[i]->GetLineColor());
    }
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

void Draw_TH1_Histograms_ratioInSameCanvas(TH1D** histograms_collection, const TString* legendList_string,TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection) {
  Draw_TH1_Histograms_MasterFunction(histograms_collection, legendList_string, histograms_collection_ratios, legendList_string_ratios, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, drawnWindowRatio, legendPlacement, legendPlacementRatio, contextPlacement, options+(std::string)"ratioInSameCanvas", optionalFitCollection);
}

void Draw_TH1_Histograms_ratioInSameCanvas(TH1D** histograms_collection, const TString* legendList_string,TH1D** histograms_collection_ratios, const TString* legendList_string_ratios, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<float, 2> drawnWindowRatio, std::array<std::array<float, 2>, 2> legendPlacement, std::array<std::array<float, 2>, 2> legendPlacementRatio, std::array<float, 2> contextPlacement, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms function; it will only be used if 'options' has fit in it
  TF1* optionalFitCollectionDummy[collectionSize];
  Draw_TH1_Histograms_ratioInSameCanvas(histograms_collection, legendList_string, histograms_collection_ratios, legendList_string_ratios, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, drawnWindowRatio, legendPlacement, legendPlacementRatio, contextPlacement, options, optionalFitCollectionDummy);
}

void Draw_TH1_Histograms(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options, TF1** optionalFitCollection) {
  TH1D* histograms_collection_ratios_dummy[collectionSize];
  const TString* legendList_string_ratios_dummy;
  std::array<float, 2> drawnWindowRatioDummy;
  std::array<std::array<float, 2>, 2> legendPlacementRatioDummy;
  Draw_TH1_Histograms_MasterFunction(histograms_collection, legendList_string, histograms_collection_ratios_dummy, legendList_string_ratios_dummy, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, drawnWindowRatioDummy, legendPlacement, legendPlacementRatioDummy, contextPlacement, options, optionalFitCollection);
}

void Draw_TH1_Histograms(TH1D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options) {
  // is here to make optionalFitCollection an actual optional parameter; Draw_TH1_Histograms can be called without, and in that case optionalFitCollection is created empty for use by the actual Draw_TH1_Histograms function; it will only be used if 'options' has fit in it
  TF1* optionalFitCollectionDummy[collectionSize];
  Draw_TH1_Histograms(histograms_collection, legendList_string, collectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options, optionalFitCollectionDummy);
}

void Draw_TH1_Histogram(TH1D* histogram, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 2> drawnWindow, std::array<std::array<float, 2>, 2> legendPlacement, std::array<float, 2> contextPlacement, std::string options) {
  TH1D* singleHistArray[1] = {histogram};
  TString dummyLegend[1] = {(TString)""};
  int dummyCollectionSize = 1;
  Draw_TH1_Histograms(singleHistArray, dummyLegend, dummyCollectionSize, Context, pdfName, texXtitle, texYtitle, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacement, options);
}

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options, TPolyLine* optionalLine) {

  double width = collectionSize*900;
  double height = 800;
  auto canvas = new TCanvas("canvas"+*pdfName, "canvas"+*pdfName, width, height);
  canvas->SetWindowSize(width + (width - canvas->GetWw()), height + (height - canvas->GetWh()));

  canvas->Divide(collectionSize,1);

  TH2D* contourHist[collectionSize];

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

void Draw_TH2_Histograms(TH2D** histograms_collection, const TString* legendList_string, int collectionSize, TString Context, TString* pdfName, TString* &texXtitle, TString* &texYtitle, TString* texCollisionDataInfo, std::array<std::array<float, 2>, 3> drawnWindow2D, double* th2Contours, int th2ContourNumber, std::string options) {
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



// void customRebin2D(TH2D* h2d, const char* histoNameSuffix) {
//     // Get the number of bins along each axis
//     int nBinsX = h2d->GetNbinsX();
//     int nBinsY = h2d->GetNbinsY();
  
//     Double_t newBinsY[] = {0,5,10,15,20,25,30,35,45,60,90,140,200};  // Example y-bin edges: [0, 3, 6, 9]
//     TString histoName = TString("h2dNew_") + histoNameSuffix;  // Name with suffix to avoid conflicts
//     // Create a new 2D histogram with the new binning
//     TH2D* h2dNew = new TH2D("h2dNew", "Custom Rebinning Y", 
//                             nBinsX, h2d->GetXaxis()->GetXmin(), h2d->GetXaxis()->GetXmax(),
//                             sizeof(newBinsY)/sizeof(newBinsY[0])-1, newBinsY);
//     // Loop over the original histogram's bins and rebin the y-axis
//     for (int i = 1; i <= nBinsX; ++i) {
//         for (int j = 1; j <= nBinsY; ++j) {
//             // Get the content of the current bin in the original histogram
//             Double_t binContent = h2d->GetBinContent(i, j);
            
//             // Find the new y-bin based on custom rebinning
//             int newBinY = h2dNew->GetYaxis()->FindBin(h2d->GetYaxis()->GetBinCenter(j));
            
//             // Fill the new histogram with the summed bin content
//             h2dNew->AddBinContent(i, newBinY, binContent);

//         }
//     }
//     delete h2d;  // Delete original histogram if no longer needed
//     h2d = h2dNew;
// }







void Plot_2D_Ratio(TH2D** histos, std::string options, int nDatasets, int refIndex) {
  if (refIndex >= nDatasets) {
      std::cerr << "Error: Reference index out of bounds!" << std::endl;
      return;
  }

  TH2D* reference = histos[refIndex];
  TH2D* ratioHisto[nDatasets];
  TH2D* ratioHistoTranspose[nDatasets];

  for (int iDataset = 0; iDataset < nDatasets; iDataset++) {
      // TCanvas* canvas = new TCanvas(Form("Ratio Jet #eta-#phi run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()) , Form("Ratio Jet #eta-#phi run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()), 800, 600);
      // if (options.find("pteta") != std::string::npos) {
      //   TCanvas* canvas = new TCanvas(Form("Ratio Jet #eta-#pt run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()) , Form("Ratio Jet #eta-#pt run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()), 800, 600);
      // } else if (options.find("ptphi") != std::string::npos){
      //   TCanvas* canvas = new TCanvas(Form("Ratio Jet #pt-#phi run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()) , Form("Ratio Jet #pt-#phi run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()), 800, 600);
      // } else TCanvas* canvas = new TCanvas(Form("Ratio Jet #eta-#phi run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()) , Form("Ratio Jet #eta-#phi run %s/run %s", DatasetsNames[iDataset].Data(), DatasetsNames[refIndex].Data()), 800, 600);

      if (iDataset == refIndex) continue; // Skip the reference itself

      ratioHisto[iDataset] = (TH2D*)histos[iDataset]->Clone(Form("Ratio_%d_to_%d", iDataset, refIndex));
      ratioHisto[iDataset]->SetTitle(Form("Ratio: run %d / run %d", iDataset, refIndex));

      for (int xbin = 1; xbin <= histos[iDataset]->GetNbinsX(); xbin++) {
          for (int ybin = 1; ybin <= histos[iDataset]->GetNbinsY(); ybin++) {
              double refValue = reference->GetBinContent(xbin, ybin);
              double value = histos[iDataset]->GetBinContent(xbin, ybin);
              if (refValue != 0) {
                  ratioHisto[iDataset]->SetBinContent(xbin, ybin, value / refValue);
              } else {
                  ratioHisto[iDataset]->SetBinContent(xbin, ybin, 0); // Avoid division by zero
              }
          }
      }
      ratioHisto[iDataset]->GetZaxis()->SetRangeUser(0.85, 1.15);
      // ratioHisto[iDataset]->Draw("COLZ");
      
  }

  std::string canvasName = "c2";
  std::string canvasTitle = "H2D_jetetaPhi_ratio";

  if (options.find("track") != std::string::npos) {
      canvasName = "c2_tracks";
      canvasTitle = "H2D_tracketaPhi_ratio";
  } else if (options.find("jet") != std::string::npos) {
      canvasName = "c2_jets";
      canvasTitle = "H2D_jetetaPhi_ratio";
  }

  TCanvas* c2 = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 1200, 800);

  //TCanvas* c2 = new TCanvas("c2", "H2D_jetetaPhi_ratio_2", 1200, 800);

  int nCols = ceil(sqrt(nDatasets));
  int nRows = ceil((double)nDatasets / nCols);
  c2->Divide(nCols, nRows, 0.01 , 0.01);

  // Loop over the histograms and draw them
  for (int i = 1; i < nDatasets; i++) {
      c2->cd(i); // Canvas pads are 1-indexed
      if (ratioHisto[i]) {
        if (options.find("jet") != std::string::npos){

          // Assume your original histogram is called ratioHisto[i] (with phi on X and eta on Y)
          int nbinsX = ratioHisto[i]->GetNbinsX();
          int nbinsY = ratioHisto[i]->GetNbinsY();

          // Create a new histogram with axes swapped
          ratioHistoTranspose[i] = new TH2D(Form("h2_transposed_%d", i),
                                            ratioHisto[i]->GetTitle(),
                                            nbinsY, ratioHisto[i]->GetYaxis()->GetXmin(), ratioHisto[i]->GetYaxis()->GetXmax(),
                                            nbinsX, ratioHisto[i]->GetXaxis()->GetXmin(), ratioHisto[i]->GetXaxis()->GetXmax());

          // Copy contents with i<->j swap
          for (int ix = 1; ix <= nbinsX; ++ix) {
              for (int iy = 1; iy <= nbinsY; ++iy) {
                  double content = ratioHisto[i]->GetBinContent(ix, iy);
                  ratioHistoTranspose[i]->SetBinContent(iy, ix, content); // Swap i<->j
              }
          }

          /// Set axis labels
          ratioHistoTranspose[i]->GetXaxis()->SetTitle(ratioHisto[i]->GetYaxis()->GetTitle());
          ratioHistoTranspose[i]->GetYaxis()->SetTitle(ratioHisto[i]->GetXaxis()->GetTitle());

          // Draw the transposed histogram
          ratioHistoTranspose[i]->GetZaxis()->SetRangeUser(0.85, 1.15);
          ratioHistoTranspose[i]->Draw("COLZ");

        } else if (options.find("track") != std::string::npos){
          ratioHisto[i]->Draw("COLZ");
          }
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.5, 0.92, Form("Ratio: %s / %s", DatasetsNames[i].Data(), DatasetsNames[refIndex].Data()));
      } else {
          std::cerr << "Warning: ratioHisto[" << i << "] is null!" << std::endl;
        }
  }

  c2->Update();
   
  // for (int i = 0; i < nDatasets; i += plotsPerCanvas) {
  //     TString canvasName = Form("R%d", canvasIndex);
  //     TString canvasTitle = Form("ratio_jetetaPhi %d", canvasIndex + 1);
  //     TCanvas* R = new TCanvas(canvasName, canvasTitle, 1200, 800);

  //     int nPlotsThisCanvas = std::min(plotsPerCanvas, nDatasets - i);
  //     int nCols = ceil(sqrt(nPlotsThisCanvas));
  //     int nRows = ceil((double)nPlotsThisCanvas / nCols);
  //     R->Divide(nCols, nRows);

  //     for (int j = 0; j < nPlotsThisCanvas; ++j) {
  //         int plotIndex = i + j;
  //         R->cd(j + 1);
  //         ratioHisto[plotIndex]->Draw("COLZ");

  //         TLatex latex;
  //         latex.SetNDC();
  //         latex.SetTextSize(0.04);
  //         latex.DrawLatex(0.1, 0.92, Form("Ratio: %s / %s",
  //                             DatasetsNames[plotIndex].Data(), DatasetsNames[refIndex].Data()));
  //     }

  //     R->Update();
  //     ++canvasIndex;
  // }




}







#endif