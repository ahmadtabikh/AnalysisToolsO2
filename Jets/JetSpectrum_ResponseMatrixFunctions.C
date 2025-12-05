#ifndef JETSPECTRU_RESPONSEMATRIXFUNCTIONS_C
#define JETSPECTRU_RESPONSEMATRIXFUNCTIONS_C

#include "JetSpectrum_ResponseMatrixFunctions.h"

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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// Response matrix functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // matrix product of fluct response times det response; assumes the two are of the same size binning wise
  // Careful: xy of hist and ij of Resp(i,j) are inverted ! hist(j,i) = matrix(i,j) and so if matrix(i,j)=SUM(matrixA(i,k)matrixB(k,j)) then hist(j,i)=SUM(histA(k,i)histB(j,k)), and if we replace j,i by gen,rec we get hist(gen,rec)=SUM(histA(k,rec)histB(gen,k))
  H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning = (TH2D*)GetMatrixProductTH2xTH2(H2D_jetPtResponseMatrix_fluctuations, H2D_jetPtResponseMatrix_detectorResponse).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning"+partialUniqueSpecifier);

  if (drawIntermediateResponseMatrices) {
    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }
    
    TH2D* H2D_jetPtResponseMatrix_fineBinningPreTransforms = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("H2D_jetPtResponseMatrix_fineBinningPreTransforms"+partialUniqueSpecifier);

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName_preRebin = new TString("ResponseMatrices/responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_preRebin_logz = new TString("ResponseMatrices/responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_fineBinningPreTransforms).Clone("responseMatrix_combined_fineBinningPreTransforms"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }
  // cout << "bin(topleft 1,N) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetNbinsY()) << endl;
  // cout << "bin(bottom left 1,1) = " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->GetBinContent(1,1) << endl;
}


void Get_PtResponseMatrix_DetectorAndFluctuationsCombined(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options) {
  // https://github.com/alisw/AliPhysics/blob/master/PWGJE/PWGJE/AliAnaChargedJetResponseMaker.cxx for ann example that works, by marta verveij

  // function should be removed when time allows as it doesn't do anything more than that Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning()

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);
}

void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
   TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
 
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // TH2D* H2D_jetPtResponseMatrix_preReweightWithPrior = H2D_jetPtResponseMatrix->Clone();

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      } else {
        Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      }
    } else {
       if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      } else {
        Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, false, options); 
      }
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
    // for (int i = 1; i < priorSpectrumWeighting->GetNbinsX(); i++)
    // {
    //   cout << "prior(" << i << ")" << priorSpectrumWeighting->GetBinContent(i)<< endl;
    //   cout << "responseIntegralLine(" << i << ")" << H2D_jetPtResponseMatrix->Integral(1,H2D_jetPtResponseMatrix->GetNbinsX(), i, i)<< endl;
    // }
    
  }
  if (options.find("mcdPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } 
    } else {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      }
    } else {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    } else {
      if (matrixTransformationOrder == 0 || matrixTransformationOrder == 3) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options);
      } else {
        Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
      }
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
  }
  // cout << "((((((((((((((((((((((((()))))))))))))))))))))))))" << endl;
  // cout << "pre norm that shouldn't be" << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(1, N, 1, 1)" << H2D_jetPtResponseMatrix->Integral(1, H2D_jetPtResponseMatrix->GetNbinsX(), 1, 1) << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1)" << H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1) << endl;
  // // NormaliseXSlicesToOneNoUnderOverFlows(H2D_jetPtResponseMatrix);
  // // NormaliseYSlicesToOneNoUnderOverFlows(H2D_jetPtResponseMatrix);
  // cout << "post norm that shouldn't be" << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(1, N, 1, 1)" << H2D_jetPtResponseMatrix->Integral(1, H2D_jetPtResponseMatrix->GetNbinsX(), 1, 1) << endl;
  // cout << "H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1)" << H2D_jetPtResponseMatrix->Integral(0, -1, 1, 1) << endl;
  // cout << "((((((((((((((((((((((((()))))))))))))))))))))))))" << endl;


  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_postReweightWithPrior = (TH2D*)H2D_jetPtResponseMatrix->Clone("H2D_jetPtResponseMatrix_postReweightWithPrior"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }
    
    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName_preRebin = new TString("ResponseMatrices/responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_preRebin_logz = new TString("ResponseMatrices/responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_postReweightWithPrior).Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_preRebin_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }
}



void MergeResponseMatrixBins(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  // Merge bins of the combined response matrix with fine binning to get the coarse one
  TH1D* priorSpectrumMerging;
  bool debugBool = false;

  TH2D* H2D_jetPtResponseMatrix_postBinMerge;

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); //take mcp as prior by default
  if (options.find("mcpPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_mcp_fineBinning(priorSpectrumMerging, iDataset, iRadius, false, options); 
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("mcdPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_mcd_fineBinning(priorSpectrumMerging, iDataset, iRadius, options);
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("measuredPriorMerging") != std::string::npos) {
    priorSpectrumMerging->Reset("M");
    Get_Pt_spectrum_bkgCorrected_fineBinning(priorSpectrumMerging, iDataset, iRadius, options);
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_PriorWeightedBinMerging(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], priorSpectrumMerging, debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("noPriorMerging") != std::string::npos) {
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool).Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D_aliPhysics(H2D_jetPtResponseMatrix, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], debugBool)->Clone("MergeResponseMatrixBins"+partialUniqueSpecifier);
  }
  // normalising priorSpectrum with evtNorm doesn't change anything as the weighting does prior_content(i)/prior_integral()
  // dividing priorSpectrum by binwidth doesn't change anything for the same reason

  H2D_jetPtResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_postBinMerge->Clone("RespMatrix_MergeResponseMatrixBins_output_"+partialUniqueSpecifier);
  // H2D_jetPtResponseMatrix_postBinMerge = (TH2D*)RebinVariableBins2D(H2D_jetPtResponseMatrix_preRebin, nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius], ptBinsJetsRec[iRadius], ptBinsJetsGen[iRadius], true).Clone("Get_PtResponseMatrix_DetectorAndFluctuationsCombined"+partialUniqueSpecifier);


  // When looking at combined response matrix before normalisation, large bins in y will look strange, and out of place compared to other bin slices of same size; this is because it potentially merges A LOT of bins together; it'll look a lot better after normalisation:


  // debug
  // for(int iBinY = 0; iBinY <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY()+1; iBinY++){ // 0 and n+1 take underflow and overflow into account
  //   for(int iBinX = 0; iBinX <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX()+1; iBinX++){ // 0 and n+1 take underflow and overflow into account
  //     if (iBinX == 0) {
  //       cout << "iBinX = 0 --> H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(0,"<< iBinY << ") = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(0,iBinY) << endl;
  //       cout << "iBinX = 0 --> H2D_jetPtResponseMatrix_fluctuations->GetBinContent(0,"<< iBinY << ") = " << H2D_jetPtResponseMatrix_fluctuations->GetBinContent(0,iBinY) << endl;
  //     }
  //   }
  // }
  if (drawIntermediateResponseMatrices) {
    // TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postBinMerge = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postBinMerge"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_postBinMerge"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_postBinMerge"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_postBinMerge).Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  // TransformRawResponseToYieldResponse(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined);
  // cout << "Should I normalise the combined matrix?" << endl;
  // cout << "     Marta doesn't do it, but it looks like I need it due to merging (some bins are very large)" << endl;
  // cout << "     it's actually done in AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin" << endl;
}



void NormYSlicesAndScaleRespByWidth(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  if (doYSliceNormToOneCombinedResp) {
    NormaliseYSlicesToOne(H2D_jetPtResponseMatrix); // Marta doesn't do it, but it looks like I need it due to merging (some bins are very large) ; actually marta probably uses it: AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin does it by default inside the rebinning function, and it takes into account the whole Fine range (1fine, Nfine)
  } 
  if (scaleRespByXYWidth) {
    H2D_jetPtResponseMatrix->Scale(1., "width");
  }
  if (scaleRespByYWidth) {
    double binWidthY;
    for(int iBinGen = 1; iBinGen <= H2D_jetPtResponseMatrix->GetNbinsY(); iBinGen++){
      binWidthY = H2D_jetPtResponseMatrix->GetYaxis()->GetBinWidth(iBinGen);
      for(int iBinRec = 1; iBinRec <= H2D_jetPtResponseMatrix->GetNbinsX(); iBinRec++){
        H2D_jetPtResponseMatrix->SetBinContent(iBinRec, iBinGen, 1./binWidthY * H2D_jetPtResponseMatrix->GetBinContent(iBinRec, iBinGen));
        H2D_jetPtResponseMatrix->SetBinError(iBinRec, iBinGen, 1./binWidthY * H2D_jetPtResponseMatrix->GetBinError(iBinRec, iBinGen));
      }
    }
  }

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_postYSliceNormAndWidthNorm = (TH2D*)H2D_jetPtResponseMatrix->Clone("NormYSlicesAndScaleRespByWidth"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_postYSliceNormAndWidthNorm"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_postYSliceNormAndWidthNorm"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_postYSliceNormAndWidthNorm).Clone("responseMatrix_combined_postReweightWithPrior"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }
}

void FinaliseResponseMatrix(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  if (matrixTransformationOrder == 0) {
    ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  } else if (matrixTransformationOrder == 1) {
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  } else if (matrixTransformationOrder == 2) {
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  } else if (matrixTransformationOrder == 3) {
    ReweightResponseMatrixWithPrior(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    NormYSlicesAndScaleRespByWidth(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
    MergeResponseMatrixBins(H2D_jetPtResponseMatrix, iDataset, iRadius, options);
  }
}

void ReweightResponseMatrixWithPrior_fineBinning(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options) {
  
  // before this, all y-slices (ie pt gen slices) have been normalised to 1;means each pt gen slice has a proba of 1
  // withthis function, we give each slice a weight so that they have different normalisation values, corresponding to the prior 

  // prior choice; none by default (flat)
  TH1D* priorSpectrumWeighting;
  if (options.find("mcpPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    } else {
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, false, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("mcdPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("measuredPriorUnfolding") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    WeightMatrixWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting, doUnfoldingPriorDivision);
  }
  if (options.find("testAliPhysics") != std::string::npos) {
    if (!normGenAndMeasByNEvtsBeforeUnfolding) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(priorSpectrumWeighting, iDataset, iRadius, options); 
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(priorSpectrumWeighting, iDataset, iRadius, options); 
    }
    H2D_jetPtResponseMatrix = (TH2D*)NormalizeResponsMatrixYaxisWithPrior(H2D_jetPtResponseMatrix, priorSpectrumWeighting)->Clone(H2D_jetPtResponseMatrix->GetName()+(TString)"_testAliPhysics");
  }
}

void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius) {


  TH2D* H2D_gen_det_geoMatched;
  TH2D* H2D_gen_det_geoMatched_rebinned;
  cout << "Get_PtResponseMatrix_detectorResponse 1" << endl;
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  if(addTwoMatchingMatrixFromDifferentMC){
    TH2D* H2D_jetPtMcdjetPtMcdFirst;
    TH2D* H2D_jetPtMcdjetPtMcdSecond;

    H2D_jetPtMcdjetPtMcdFirst = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_jetPtMcdjetPtMcdSecond = (TH2D*)((TH2D*)file_O2Analysis_MCfileForSecondMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);

    H2D_jetPtMcdjetPtMcdFirst->Sumw2();
    H2D_jetPtMcdjetPtMcdSecond->Sumw2();

    TH2D *H2D_jetPtMcdjetPtMcd = (TH2D*)H2D_jetPtMcdjetPtMcdFirst->Clone("H2D_jetPtMcdjetPtMcd_Sum");
    H2D_jetPtMcdjetPtMcd->Add(H2D_jetPtMcdjetPtMcdSecond);  // This does content addition + error^2 propagation

    H2D_gen_det_geoMatched = (TH2D*)H2D_jetPtMcdjetPtMcd->Clone(partialUniqueSpecifier+"_genrec");

  }else{
    TH3D* H3D_jetRpartPtdetPt;
    TH2D* H2D_jetPtMcdjetPtMcd;

    if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
      H3D_jetRpartPtdetPt = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_PtResponseMatrix_detectorResponse"+partialUniqueSpecifier);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
      H3D_jetRpartPtdetPt->Sumw2();

      int ibinJetRadius = H3D_jetRpartPtdetPt->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
      H3D_jetRpartPtdetPt->GetXaxis()->SetRange(ibinJetRadius, ibinJetRadius);
      // project H3D onto a H2D, option "yz" means y goes on y-axis while z goes on x-axis, and so H2D_gen_det_geoMatched will be (x=mcd, y=mcp)
      H2D_gen_det_geoMatched = (TH2D*)H3D_jetRpartPtdetPt->Project3D(partialUniqueSpecifier+"_genrec_e_yz"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)
    } else if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
      if (doBkgSubtractionInMC) {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);
      } else {
        if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
          H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);
        } else {
          H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"))->Clone("Get_PtResponseMatrix_detectorResponse"+Datasets[iDataset]+DatasetsNames[iDataset]);
        }
      }
      H2D_jetPtMcdjetPtMcd->Sumw2();

      // H2D_gen_det_geoMatched = (TH2D*)GetTransposeHistogram(H2D_jetPtMcdjetPtMcd).Clone(partialUniqueSpecifier+"_genrec");
      H2D_gen_det_geoMatched = (TH2D*)H2D_jetPtMcdjetPtMcd->Clone(partialUniqueSpecifier+"_genrec");
    }
  }

  // keep (gen, gen) for the bins; rec will be introduced in the fluctuation response, and by multiplication will stay in the combined matrix
  TH2D* H2D_response = (TH2D*)RebinVariableBins2D(H2D_gen_det_geoMatched, nBinPtJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], ptBinsJetsFine[iRadius]).Clone("Get_PtResponseMatrix_detectorResponse_rebinned"+partialUniqueSpecifier);

  if (doYSliceNormToOneDetResp) {
    NormaliseYSlicesToOne(H2D_response);
  }
  if (normDetRespByNEvts) {
    if (mcIsWeighted) {
      H2D_response->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      double Nevents = GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC);
      for (int iBinX = 0; iBinX < H2D_response->GetNbinsX(); iBinX++) {
        for (int iBinY = 0; iBinY < H2D_response->GetNbinsY(); iBinY++) {
          H2D_response->SetBinContent(iBinX, iBinY, H2D_response->GetBinContent(iBinX, iBinY) * 1./Nevents);
          H2D_response->SetBinError(iBinX, iBinY, H2D_response->GetBinError(iBinX, iBinY) * 1./Nevents);
        }
      }
      // H2D_response->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
  cout << "Detector response building: errors here should probably be reduced to take into account correlations, as the normalisation factor is built from same matrix" << endl;

  // H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_gen_det_geoMatched_rebinned->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); // should be using the one below; this one is just a test
  H2D_jetPtResponseMatrix_detectorResponse = (TH2D*)H2D_response->Clone("H2D_jetPtResponseMatrix_detectorResponse"+partialUniqueSpecifier); 

  // cout << "............................................................." <<endl;
  // cout << "............................................................." <<endl;
  // for(int iBinX = 1; iBinX <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsX(); iBinX++){ // 0 and n+1 would take underflow and overflow into account, don't want that
  //   for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorResponse->GetNbinsY(); iBinY++){ // 0 and n+1 would take underflow and overflow into account, don't want that
  //     cout << "iBinX = " << iBinX << ", iBinY = " << iBinY << "         --------          detResponseContent = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinContent(iBinX, iBinY) << ", detResponseError = " << H2D_jetPtResponseMatrix_detectorResponse->GetBinError(iBinX, iBinY) << endl;
  //   }
  // }
  // cout << "............................................................." <<endl;
  // cout << "............................................................." <<endl;
}

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius) { 
  // see Hiroki Yokoyama thesis
  // iRadius is for chosing the pT binning
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  if (useFactorisedMatrix == false){
    TH2D H2D_identity = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]);
    for(int iBinRec = 0; iBinRec <= H2D_identity.GetNbinsX()+1; iBinRec++){
      H2D_identity.SetBinContent(iBinRec, iBinRec, 1);
      H2D_identity.SetBinError(iBinRec, iBinRec, 0);
    }
    H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_identity.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier);
  } else {
    cout << "I should check that the average of each ptGen slice is as displaced to the diagonal as the randomCone distrib is; ie should I use GetBinLowEdge or GetBinLowEdge+width for ptGen" << endl;


    TH2D* H2D_fluctuations_centrality;
    TH1D* H1D_fluctuations;

    H2D_fluctuations_centrality = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowBkg+"/h2_centrality_rhorandomcone"+randomConeTypeList[randomConeType]))->Clone("Get_PtResponseMatrix_Fluctuations"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_fluctuations_centrality->Sumw2();


    int ibinCent_low = H2D_fluctuations_centrality->GetXaxis()->FindBin(centralityRange[0]);
    int ibinCent_high = H2D_fluctuations_centrality->GetXaxis()->FindBin(centralityRange[1]);
    H1D_fluctuations = (TH1D*)H2D_fluctuations_centrality->ProjectionY("bkgFluctuationCentrality_highRes_"+partialUniqueSpecifier, ibinCent_low, ibinCent_high, "e");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* ctest = new TCanvas("ctest", "Fluctuations_ProjY_of_bkgFluctuationCentrality", 800, 600);
    H1D_fluctuations->GetXaxis()->SetRange(0, H1D_fluctuations->GetNbinsX() + 1);
    H1D_fluctuations->Draw(); 
    ctest->Update();
    int nBinsFluc = H1D_fluctuations->GetNbinsX();
    cout << "Number of bins = " << nBinsFluc << std::endl;
    cout << "Bin#\tLowEdge\tHighEdge\tContent\tError" << std::endl;
    for (int i = 1; i <= nBinsFluc; i++) {
        double lowEdge  = H1D_fluctuations->GetBinLowEdge(i);
        double highEdge = H1D_fluctuations->GetBinLowEdge(i+1);
        double content  = H1D_fluctuations->GetBinContent(i);
        double error    = H1D_fluctuations->GetBinError(i);
        cout << i << "\t" << lowEdge << "\t" << highEdge << "\t" << content << "\t" << error << std::endl;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName_deltapt = new TString("deltapt_beforeNormalisation");
    TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
    TString* texdeltapt = new TString("#delta#it{p}_{T} (GeV/#it{c})");
    TString* counts_deltapt = new TString("counts ");
    Draw_TH1_Histogram(H1D_fluctuations, textContext, pdfName_deltapt, texdeltapt, counts_deltapt, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    NormaliseRawHistToIntegral(H1D_fluctuations); // normalising fluctuations to 1
    // cout << "Integral H1D_fluctuations: " << H1D_fluctuations->Integral(1, H1D_fluctuations->GetNbinsX()) << endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;
    TString* pdfName = new TString("deltapt_normalizedToIntegral");
    // TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
    TString* norm_deltapt = new TString("pdf");
    Draw_TH1_Histogram(H1D_fluctuations, textContext, pdfName, texdeltapt, norm_deltapt, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    // int underflowBin = -50; // by definition
    int underflowBin = H1D_fluctuations->GetNbinsX() + 1;

    double content = H1D_fluctuations->GetBinContent(underflowBin);
    double error   = H1D_fluctuations->GetBinError(underflowBin);
    double binLowEdge = H1D_fluctuations->GetBinLowEdge(underflowBin);
    double binHighEdge = H1D_fluctuations->GetBinLowEdge(underflowBin + 1);

    cout << "H1D_fluctuations_norm Underflow bin number: " << underflowBin << std::endl;
    cout << "H1D_fluctuations_norm Underflow Bin range: [" << binLowEdge << ", " << binHighEdge << "]" << std::endl;
    cout << "H1D_fluctuations_norm Underflow Bin Content: " << content << " ± " << error << std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    TH2D H2D_response = TH2D("H2D_response_"+partialUniqueSpecifier, "H2D_response_"+partialUniqueSpecifier, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius]); 

    //==================== Build response matrix: shift deltaPt by pT gen along the pT rec axis ====================//
    int ibinZeroFluct= H1D_fluctuations->FindBin(0+GLOBAL_epsilon);
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "ibinZeroFluct = " << ibinZeroFluct << std::endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

    double integralError;
    double ptGen, ptRec_low, ptRec_up;
    int iBin_fluct_low, iBin_fluct_high;
    for(int iBinRec = 0; iBinRec <= H2D_response.GetNbinsX()+1; iBinRec++){
      for(int iBinGen = 0; iBinGen <= H2D_response.GetNbinsY()+1; iBinGen++){
        ptGen = H2D_response.GetYaxis()->GetBinCenter(iBinGen); // was bincenter before but then it'd give 0.5 values of GeV,
        ptRec_low = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec); //lower edge  = (X value) : 0
        ptRec_up = H2D_response.GetXaxis()->GetBinLowEdge(iBinRec+1); //1
        // int iBin_fluct_low = H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon);
        // int iBin_fluct_high = H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen - GLOBAL_epsilon);
        iBin_fluct_low = H1D_fluctuations->GetXaxis()->FindBin(ptRec_low - ptGen + GLOBAL_epsilon); //0-0.5 =-0.5 --> bin for -0.5 = 1
        if (iBinGen == 10 && (iBin_fluct_low >= iBin_fluct_high)) { // checks iBinRec =10 so that the message doesn't appear NbinsX*NBinsY times
          cout << "Get_PtResponseMatrix_Fluctuations: some bins are counted twice in the integral, binning needs to be looked at, right now the fluctuation matrix has too many entries" << endl;
        }
        iBin_fluct_high = H1D_fluctuations->GetXaxis()->FindBin(ptRec_up - ptGen - GLOBAL_epsilon)-1; // 1-0.5 =0.5 --> bin for 0.5 = 2,  //FindBin(x) returns the bin index corresponding to a given x-value in a histogram.
        H2D_response.SetBinContent(iBinRec, iBinGen, H1D_fluctuations->IntegralAndError(iBin_fluct_low, iBin_fluct_high, integralError)); //in (0,0) bin, content is the integration from bin 1 to 2
        H2D_response.SetBinError(iBinRec, iBinGen, integralError); 
        // if (iBinRec == 10) {
        //   cout << "DeltaBin = " << iBin_fluct_high - iBin_fluct_low << ", iBin_low = " << iBin_fluct_low << ", iBin_high = " << iBin_fluct_high << endl;
        // }
      }
    }
    float integralCheck = H2D_response.Integral(1, H2D_response.GetNbinsX(), 10, 10);
    if (!(integralCheck > 0.999)) {
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "!!!!!!Response matrix of fluctuations does not have line integral equal to 1; something is wrong with binning!!!!!!" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    cout << "Integral over X for Y bin 10 = " << integralCheck << std::endl;
    float integralCheck120 = H2D_response.Integral(1, H2D_response.GetNbinsX(), 120, 120);
    cout << "Integral over X for Y bin 120 = " << integralCheck120 << std::endl;



    //========================================= Build response matrix end =========================================//

    H2D_jetPtResponseMatrix_fluctuations = (TH2D*)H2D_response.Clone("H2D_jetPtResponseMatrix_fluctuations"+partialUniqueSpecifier);
  }
}

#endif