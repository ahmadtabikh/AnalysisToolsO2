#ifndef JETSPECTRUM_UNFOLDING_C
#define JETSPECTRUM_UNFOLDING_C

#include "JetSpectrum_Unfolding.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
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
///////////////////////////////////////////////////////////////////////////// RooUnfold Custom Utilities ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int GetSvdBestRegularisationParameter_notYetSatisfying(TSVDUnfold* unfoldTSvd){
  //alternative ideas: average over 3-4 k values, find when ratio of those as k increases becomes 1ish
  TH1D* H1D_D = unfoldTSvd->GetD();
  double minimum = 999999;
  double tempContent;
  int k = -99;
  for (int iBinX = 1; iBinX < H1D_D->GetNbinsX(); iBinX++) {
    tempContent = abs(abs(H1D_D->GetBinContent(iBinX)) - 1); // one wants |k| closest to 1 as possible
    if (minimum > tempContent) {
      minimum = tempContent;
      k = iBinX - 1; // or is it iBinX?
    }
  }
  cout << "k = " << k << ", |H1D-1| minimum = " << minimum << endl;
  return k;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// Spectrum Unfolding functions ////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm()" << endl;};
  //for now makes the assumption gen and rec have the same pT binning

  bool controlMC_useResponseSplit_response = true; // has no effect if controlMC=false; is used with spectraGetters if one should draw from response split
  bool controlMC_useResponseSplit_input = false; // has no effect if controlMC=false; is used with spectraGetters if one should draw from input split

  TString partialUniqueSpecifier = Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+"_R="+Form("%.1f",arrayRadius[iRadius])+"unfoldParameterInput"+Form("%.1d",unfoldParameterInput)+"controlMC"+Form("%.1i",controlMC)+"inputIsGen"+Form("%.1i",inputIsGen);

  TH1D* measured = (TH1D*)measuredInput->Clone("measured_Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);;
  TH1D* mcp; // used if useManualRespMatrixSettingMethod is false, for calculating efficiency and purity
  TH1D* mcd; // used if useManualRespMatrixSettingMethod is false, for calculating efficiency and purity


  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    // Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); now fed as input to function
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response); 
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);

    if (useFineBinningTest) {
      // Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);
    }
  } else{
    // Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);

    if (useFineBinningTest) {
      // Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_response);
    }
  }

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useManualRespMatrixSettingMethod) {
    if (applyFakes && (options.find("noPurity") == std::string::npos)) { // if applyFakes AND if option noPurity has not been found in options; necessary check for the mcp-folded unfolding test as we don't want the fake correction
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options, controlMC);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options, controlMC);
      }
      if (divideSuccessFakes){
        measured->Multiply(H1D_jetFakes);
      } else {
        cout << "################## measured->Multiply(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }

  // TH1D* measured;
  // Get_Pt_spectrum_bkgCorrected_recBinning(measured, iDataset, iRadius, options);
  // TH1D* mcp;
  // Get_Pt_spectrum_mcp_genBinning(mcp, iDataset, iRadius, options);
  // TH1D* mcdMatched;
  // Get_Pt_spectrum_mcdMatched_genBinning(mcdMatched, iDataset, iRadius, options);

  TH1D* H1D_kinematicEfficiency;
  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  
  std::string optionsFluctResp = options.find("useIdentityForFluctResp") != std::string::npos ? (std::string)"useIdentityForFluctResp" : (std::string)"";
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, optionsFluctResp); // gets fluct matrix in fine binning 
  std::string optionsDetResp = options.find("useIdentityForDetResp") != std::string::npos ? (std::string)"useIdentityForDetResp" : (std::string)"";
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, optionsDetResp, controlMC); // gets det matrix in fine binning 
  // compute matrixFluctuations times matrixDetector

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  FinaliseResponseMatrix_priorAndNormYslicesAndMergeBins(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options, controlMC);

  // Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }
    
    TString priorInfo = (TString)unfoldingPrior;

    TString* pdfNamePost = new TString("ResponseMatrices/responseMatrix_combined_postWeighting"+(TString)partialUniqueSpecifier);
    TString* pdfNamePost_logz = new TString("ResponseMatrices/responseMatrix_combined_postWeighting"+(TString)partialUniqueSpecifier+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting).Clone("responseMatrix_combined_postWeighting"+(TString)partialUniqueSpecifier+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting->Clone("responseMatrix_combined_postWeighting"+(TString)partialUniqueSpecifier+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfNamePost, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfNamePost_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "RooUnfoldResponse setting - start" << endl;
  RooUnfoldResponse* response;
  if (useManualRespMatrixSettingMethod) {

    // // based on Marta's work: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JEJetSpectrumUnfolding

    response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measuredInput and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation


    cout << "Get_Pt_spectrum_unfolded(): should I use response->UseOverflow() ? using it gives a ratio unfolded/mcp much higher than without using it" << endl;
    if (useMatrixOverflows) {
      response->UseOverflow(); // is off by default
    }
    // cout << "UseOverflowStatus:" << response->UseOverflowStatus() << endl;

  } else {
    TH2D *Respt;
    Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
    TH1D *mcdMatched = (TH1D*)Respt->ProjectionX("Respt_projX", 1, Respt->GetNbinsY());
    TH1D *mcpMatched = (TH1D*)Respt->ProjectionY("Respt_projY", 1, Respt->GetNbinsX());
    TH1D *fake = (TH1D*)mcd->Clone();
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        fake->Scale(1./GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
      } else {
        fake->Scale(1./GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Fake i = " << i << " content: "<< fake->GetBinContent(i) << " - " << mcdMatched->GetBinContent(i) << " = " << fake->GetBinContent(i) - mcdMatched->GetBinContent(i) << " ------ error of mcd = " << fake->GetBinError(i) << ", error of mcdMatched = " << mcdMatched->GetBinError(i) << endl;
    }
    fake->Add(mcdMatched, -1);
    TH1D *miss = (TH1D*)mcp->Clone(); // for each bin of gen-level jet distribution, how many are not matched to a rec-level jet
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        miss->Scale(1./GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
      } else {
        miss->Scale(1./GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Miss i = " << i << " content: "<< miss->GetBinContent(i) << " - " << mcpMatched->GetBinContent(i) << " = " << miss->GetBinContent(i) - mcpMatched->GetBinContent(i) << " ------ error of mcp = " << miss->GetBinError(i) << ", error of mcpMatched = " << mcpMatched->GetBinError(i) << endl;
    }
    miss->Add(mcpMatched, -1); // for each bin of rec-level jet distribution, how many are not matched to a gen-level jet

    response = new RooUnfoldResponse(mcd, mcp);
    for (auto i = 1; i <= Respt->GetNbinsX(); i++) {
      for (auto j = 1; j <= Respt->GetNbinsY(); j++) { // ptpair
        Double_t bincenx = Respt->GetXaxis()->GetBinCenter(i);
        Double_t binceny = Respt->GetYaxis()->GetBinCenter(j);
        Double_t bincont = Respt->GetBinContent(i, j);
        response->Fill(bincenx, binceny, bincont);
        // cout << "response i,j = " << i << "," << j << " content: "<< bincont << endl;
      }
    }
    for (auto i = 1; i <= miss->GetNbinsX(); i++) { 
      Double_t bincenx = miss->GetXaxis()->GetBinCenter(i);
      Double_t bincont = miss->GetBinContent(i);
      response->Miss(bincenx, bincont);
      // cout << "Miss i = " << i << " content: "<< bincont << endl;
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) {
      Double_t bincenx = fake->GetXaxis()->GetBinCenter(i);
      Double_t bincont = fake->GetBinContent(i);
      response->Fake(bincenx, bincont);
      // cout << "Fake i = " << i << " content: "<< bincont << endl;
    }


    // Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, (TH2D*)response->Hresponse(), partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
    if (useMatrixOverflows) {
      response->UseOverflow(); // is off by default
    }
  }
  cout << "RooUnfoldResponse setting - end" << endl;



  // note about unfolding: if SUM(Mkj over k)=1 then it's also true for the inverse of M (can easily demosntrate it, just write it);
  RooUnfoldBayes* unfoldBayes = new RooUnfoldBayes(response, measured, unfoldParameterInput);
  int unfoldParameterSvdInitial = 1;
  RooUnfoldSvd* unfoldSvdInitialiser = new RooUnfoldSvd(response, measured, unfoldParameterSvdInitial); // instance of RooUnfoldSvd only used to find the best regularisation parameter; the unfolded spectrum returned by it is not retrived
  RooUnfold* unfold = unfoldBayes; // default Bayes
  TH1D* hist_unfold;

  int unfoldParameter;
  if (options.find("Svd") != std::string::npos) {
    unfoldSvdInitialiser->Hreco(); // necessary to have GetD() give a meaningful output
    TSVDUnfold *tsvdUnfold = (TSVDUnfold*)unfoldSvdInitialiser->Impl();
    if (automaticBestSvdParameter) {
      unfoldParameter = GetSvdBestRegularisationParameter_notYetSatisfying(tsvdUnfold);
    } else {
      unfoldParameter = unfoldParameterInput;
    }
    RooUnfoldSvd* unfoldSvd = new RooUnfoldSvd(response, measured, unfoldParameter); // the RooUnfoldSvd instance that is actually used to unfold, with the best regularisation parameter

    // unfoldSvd->SetRegParm(-200);  // doesnt seem to change anything, even if I put 0 here ... annoying; but accoring to mailing list of roounfold it doesn't work well
    hist_unfold = (TH1D*)(unfoldSvd->Hreco());
    unfold = unfoldSvd;

    // plot svd d distribution
    TH1D* H1D_D = tsvdUnfold->GetD();
    TString inputUnfoldingName = (options.find("inputIsMCPFoldedTest") != std::string::npos) ? "_mcpFoldedTestInput" : "";
    TString* pdfName_regparam = new TString("Svd_regularisationd_distribution_"+(TString)partialUniqueSpecifier+"_"+(TString)unfoldingPrior+inputUnfoldingName);
    TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));
    std::array<std::array<float, 2>, 2> drawnWindowSvdParam = {{{0, 30}, {0.01, 10000}}}; // {{xmin, xmax}, {ymin, ymax}}
    Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdK, texSvdDvector_k, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");

  } else if (options.find("Bayes") != std::string::npos) {
    unfoldParameter = unfoldParameterInput;
    hist_unfold = (TH1D*)(unfoldBayes->Hreco());
    unfold = unfoldBayes;
  } else {
    cout << "unfolding technique name is wrong, please change it" << endl;
    exit(EXIT_SUCCESS);
  }


  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  //applying efficiciencies and purity
  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options, controlMC);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options, controlMC);
  }

  if (useManualRespMatrixSettingMethod) {
    if (applyEfficiencies > 0) {
      if (applyEfficiencies == 2 || applyEfficiencies == 3) {
        if (divideSuccessEff){
          if (options.find("noEff") == std::string::npos) { // if option noEff has not been found in options; necessary check for the mcp-folded unfolding test as we don't want the eff correction
            H1D_jetPt_unfolded->Divide(H1D_jetEfficiency);
          }
        } else {
          cout << "################## Get_Pt_JetEfficiency FAILED!!!!! in Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm ##################" << endl;
        }
      }
      if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
        if (options.find("noKineEff") == std::string::npos) { // if option noKineEff has not been found in options; necessary check for the mcp-folded unfolding test as we don't want the kine eff correction
          H1D_jetPt_unfolded->Divide(H1D_kinematicEfficiency);
        }
      }
    }
  }

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding = (TH2D*)unfold->response()->Hresponse()->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }
    
    TString priorInfo = (TString)unfoldingPrior;

    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_postUnfolding"+(TString)partialUniqueSpecifier);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_postUnfolding"+(TString)partialUniqueSpecifier+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding).Clone("responseMatrix_combined_postUnfolding"+(TString)partialUniqueSpecifier+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding->Clone("responseMatrix_combined_postUnfolding"+(TString)partialUniqueSpecifier+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfolded);
  // }
  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm()" << endl;};

  std::pair<int, RooUnfold*> unfoldInfo(unfoldParameter, unfold);
  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_unfolded_preWidthScalingAtEnd()" << endl;};

  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen);
  if (normaliseUnfoldingResultsAtEnd){
    if (!controlMC && options.find("inputIsMC") == std::string::npos) { // if option controlMC is false, and if inputIsMC has not been found in options; necessary check for the mcp-folded unfolding test as we want to normalise by the number of events in the MC file
      NormaliseRawHistToInputN(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]));
      cout<< "normaliseUnfoldingResultsAtEnd: normalising by number of events in data file, which is " << GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]) << endl;
    } else {
      TFile* fileNorm = controlMC ? file_O2Analysis_MCfile_UnfoldingControl_input[iDataset] : file_O2Analysis_MCfile_GeneralResponse[iDataset];
      if (mcIsWeighted) {
        if (!inputIsGen) {
          NormaliseRawHistToInputN(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_weighted(fileNorm, analysisWorkflow_unfoldingControl));
        } else {
          NormaliseRawHistToInputN(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_gen_weighted(fileNorm, analysisWorkflow_unfoldingControl));
        }
      } else {
        if (!inputIsGen) {
          NormaliseRawHistToInputN(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(fileNorm, analysisWorkflow_unfoldingControl));
        } else {
          NormaliseRawHistToInputN(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_gen(fileNorm, analysisWorkflow_unfoldingControl));
        }
      }      
    }
  }

  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_unfolded_preWidthScalingAtEnd()" << endl;};
  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_unfolded()" << endl;};

  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAtEnd(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfolded);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_unfolded, deltaJetEta[iRadius]);

  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_unfolded()" << endl;};

  return unfoldInfo;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Spectrum Refolding tests //////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm()" << endl;};

  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];

  bool controlMC_useResponseSplit_response = true; // has no effect if controlMC=false; is used with spectraGetters if one should draw from response split
  bool controlMC_useResponseSplit_input = false; // has no effect if controlMC=false; is used with spectraGetters if one should draw from input split

  TString partialUniqueSpecifier = Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+"_R="+Form("%.1f",arrayRadius[iRadius])+"unfoldParameterInput"+Form("%.1d",unfoldParameterInput)+"controlMC"+Form("%.1i",controlMC)+"inputIsGen"+Form("%.1i",inputIsGen);

  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen).second;
  // Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: " << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  // }

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "", controlMC);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  FinaliseResponseMatrix_priorAndNormYslicesAndMergeBins(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options, controlMC);


  // Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  

  if (drawIntermediateResponseMatrices) {
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding"+partialUniqueSpecifier);

    struct stat st1{};
    if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
        mkdir("pdfFolder/ResponseMatrices", 0700);
    }
    struct stat st2{};
    if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
        mkdir("pngFolder/ResponseMatrices", 0700);
    }

    TString priorInfo = (TString)unfoldingPrior;

    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_duringRefolding"+(TString)partialUniqueSpecifier);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_duringRefolding"+(TString)partialUniqueSpecifier+"_logz");

    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
    TH2D* MatrixResponse;
    TString* xLabel;
    TString* yLabel;
    if (transposeResponseHistogramsInDrawing) {
      MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding).Clone("responseMatrix_combined_duringRefolding"+(TString)partialUniqueSpecifier+"_"+priorInfo);
      xLabel = texPtJetGen;
      yLabel = texPtJetRec;
    } else {
      MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding->Clone("responseMatrix_combined_duringRefolding"+(TString)partialUniqueSpecifier+"_"+priorInfo);
      xLabel = texPtJetRec;
      yLabel = texPtJetGen;
    }

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: pre efficiency" << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
  //   cout << "  H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << endl;
  //   cout << "e_H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  // }
  

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options, controlMC);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options, controlMC);
  }

  if (applyEfficiencies > 0) {
    if (applyEfficiencies == 2 || applyEfficiencies == 3) {
      if (divideSuccessEff){
        H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
      } else {
        cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
      }
    }
    if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
      H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
    }
  }

  TH2D* refoldingResponseMatrix;
  if (useManualRespMatrixSettingMethod) {
    refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
  } else {
    if (useMatrixOverflows) {
      refoldingResponseMatrix = (TH2D*)unfold->response()->Hresponse();
    } else {
      refoldingResponseMatrix = (TH2D*)unfold->response()->HresponseNoOverflow();
    }
    // refoldingResponseMatrix = (TH2D*)unfold->response()->Hresponse();
    // cout << "refold why use HresponseNoOverflow instead of Hresponse?" << endl;
  }

  if (normaliseRespYSliceForRefold){
    NormaliseYSlicesToOne(refoldingResponseMatrix);
  }
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1((TH2D*)unfold->response()->Hresponse(), H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier); //gives the same
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier); //gives the same
  cout << "--UnfoldedJetPt in RooUnfold Refold Method--" << endl;
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier); //gives the same


  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (applyFakes) { // for useManualRespMatrixSettingMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
    if (!useFineBinningTest){ 
      divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options, controlMC);
    } else {
      divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options, controlMC);
    }
    if (divideSuccessFakes){
      H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
    } else {
      cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
    }
  }

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: post efficiency" << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfoldedThenRefolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfoldedThenRefolded(" << iBinX << ") = " << H1D_jetPt_unfoldedThenRefolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfoldedThenRefolded->GetBinError(iBinX) << endl;
  // }
  

  // cout << "Refolding check of Combined matrix = " << endl; 
  // for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
  //   cout << "combinedMatrix integral of slice iBinY = " << iBinY << " is: " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Integral(0, -1, iBinY, iBinY) << endl;
  // }


  // if (useFineBinningTest) {
  //   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);
  // } else {
  //   H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);
  // }
  

  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  // }

  // cout << "still not giving back the measured used as input to the unfolding; got an issue somewhere" << endl;
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen);

  if (normaliseUnfoldingResultsAtEnd){
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]));
    } else {
      if (mcIsWeighted) {
        if (!inputIsGen) {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
        } else {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
        }
      } else {
        if (!inputIsGen) {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
        } else {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
        }
      }
    }
  }

  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm()" << endl;};
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_dataUnfoldedThenRefolded()" << endl;};

  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEnd(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
  
  TransformYieldToEtaDifferentialYield(H1D_jetPt_unfoldedThenRefolded, deltaJetEta[iRadius]);

  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_dataUnfoldedThenRefolded()" << endl;};
}


void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm()" << endl;};

  // ApplyToTruth function doesn't apply errors on folding

  // Matches exactly with the manual method IF NO PRIOR
  // if I have a non flat prior, then the roounfold method gives me a good closure test, but not the manual method!

  TH1D* H1D_jetPt_unfolded;

  bool controlMC_useResponseSplit_response = true; // has no effect if controlMC=false; is used with spectraGetters if one should draw from response split
  bool controlMC_useResponseSplit_input = false; // has no effect if controlMC=false; is used with spectraGetters if one should draw from input split

  TString partialUniqueSpecifier = Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+"_R="+Form("%.1f",arrayRadius[iRadius])+"unfoldParameterInput"+Form("%.1d",unfoldParameterInput)+"controlMC"+Form("%.1i",controlMC)+"inputIsGen"+Form("%.1i",inputIsGen);

  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen).second;

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "", controlMC);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation to get the kinematic efficiency
  FinaliseResponseMatrix_priorAndNormYslicesAndMergeBins(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options, controlMC);

  // Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options, controlMC);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options, controlMC);
  }

  if (useManualRespMatrixSettingMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyEfficiencies > 0) {
      if (applyEfficiencies == 2 || applyEfficiencies == 3) {
        if (divideSuccessEff){
          H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency);
        } else {
          cout << "################## H1D_jetPt_unfolded->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
        }
      }
      if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
        H1D_jetPt_unfolded->Multiply(H1D_kinematicEfficiency);
      }
    }
  }

  // cout << "Refolding check of Combined matrix = " << endl; 
  // for(int iBinY = 1; iBinY <= H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->GetNbinsY(); iBinY++){
  //   cout << "combinedMatrix integral of slice iBinY = " << iBinY << " is: " << H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Integral(0, -1, iBinY, iBinY) << endl;
  // }


  // RooUnfoldResponse* rooresponse_postUnfold = unfold->response();
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)unfold->response()->ApplyToTruth(H1D_jetPt_unfolded, "H1D_jetPt_unfoldedThenRefolded_withApplyToTruth")->Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);

  // RooUnfoldResponse* response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measured and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)unfold->response()->ApplyToTruth(H1D_jetPt_unfolded, "H1D_jetPt_unfoldedThenRefolded_withApplyToTruth")->Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);


  cout << "--UnfoldedJetPt in RooUnfold Refold Method--" << endl;
  for(int iBinK = 1; iBinK <= H1D_jetPt_unfolded->GetNbinsX(); iBinK++){ // 0 and n+1 take underflow and overflow into account
    cout << "-----" << endl;
    cout << " U(" << iBinK << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinK) << endl;
    cout << "eU(" << iBinK << ") = " << H1D_jetPt_unfolded->GetBinError(iBinK) << endl;
  }

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useManualRespMatrixSettingMethod) { // this condition isn't present in the manual refolding method as in the manual case the response matrix object is just a th2 and not a roounfoldresponse object and doesn't have the efficiencies/fakes encoded in it
    if (applyFakes) {
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options, controlMC);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options, controlMC);
      }
      if (divideSuccessFakes){
        H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes);
      } else {
        cout << "################## H1D_jetPt_unfoldedThenRefolded->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
      }
    }
  }

  //not sure I should normalise --> probably not as H1D_jetPt_unfolded is already normalised in Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm

  // cout << "still not giving back the measured used as input to the unfolding; got an issue somewhere" << endl;
  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm()" << endl;};
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd()" << endl;};

  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen);

  if (normaliseUnfoldingResultsAtEnd){
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]));
    } else {
      if (mcIsWeighted) {
        if (!inputIsGen) {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
        } else {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
        }
      } else {
        if (!inputIsGen) {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
        } else {
          NormaliseRawHistToInputN(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
        }
      }
    }
  }
  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd()" << endl;};
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options, bool controlMC = false, bool inputIsGen = false) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod()" << endl;};

  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options, controlMC, inputIsGen);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_unfoldedThenRefolded, deltaJetEta[iRadius]);

  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod()" << endl;};
}


void Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm()" << endl;};

  TString partialUniqueSpecifier = Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* H1D_jetPt_mcp_control;
  // Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_control, iDataset, iRadius, true, options);
  // Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_control, iDataset, iRadius, true, options);
  // if (useFineBinningTest) {
  bool controlMC = false; //always false here because mcpFolded test doesn't use MC info for the unfolding: we only unfold with fluct matrix from data
  bool controlMC_useResponseSplit_input = false; // doesn't matter because controlMC is false

  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_control, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit_input);
  // }

  // we do not apply efficiencies nor fakes on this folding exercise with fluctuations
  // applying efficiency loss and fake enhancement -> nope
  // bool divideSuccessEff;
  // TH1D* H1D_jetEfficiency;
  // if (!useFineBinningTest){ 
  //   divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  // } else {
  //  divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  // }

  // we do not apply efficiencies nor fakes on this folding exercise with fluctuations
  // if (applyEfficiencies > 0) {
  //   if (applyEfficiencies == 2 || applyEfficiencies == 3) {
  //     if (divideSuccessEff){
  //       H1D_jetPt_mcp_control->Multiply(H1D_jetEfficiency);
  //     } else {
  //       cout << "################## H1D_jetPt_mcp_control->Multiply(H1D_jetEfficiency) failed because Get_Pt_JetEfficiency() FAILED!!!!! ##################" << endl;
  //     }
  //   }
  //   if ((applyEfficiencies == 1 || applyEfficiencies == 3) && !useFineBinningTest) {
      // H1D_jetPt_mcp_control->Multiply(H1D_kinematicEfficiency);
  //   }
  // }


  TH2D* refoldingResponseMatrix;
  Get_PtResponseMatrix_Fluctuations(refoldingResponseMatrix, iDataset, iRadius, "");
  cout << "Integral Line 1 fluct matrix: " << refoldingResponseMatrix->Integral(1, refoldingResponseMatrix->GetNbinsX(), 1, 1) << endl;
  cout << "Integral Line 10 fluct matrix: " << refoldingResponseMatrix->Integral(1, refoldingResponseMatrix->GetNbinsX(), 10, 10) << endl;
  cout << "Integral Line 20 fluct matrix: " << refoldingResponseMatrix->Integral(1, refoldingResponseMatrix->GetNbinsX(), 20, 20) << endl;

  TH1D* H1D_jetPt_mcp_control_folded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_mcp_control).Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);
  // TH1D* H1D_jetPt_mcp_control_folded = (TH1D*)H1D_jetPt_mcp_control->Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);

  if (foldMcpWithFluct_alsoFoldWithDetResp) {
    TH2D* H2D_jetPtResponseMatrix_detectorResponse;
    TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
    bool controlMC = true; // has no effect if controlMC=false; is used with spectraGetters if one should draw from response split
    Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "", controlMC);
    Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, refoldingResponseMatrix, iDataset, iRadius, options, controlMC);
    H1D_jetPt_mcp_control_folded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_mcp_control).Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);
  }

  if (!useFineBinningTest) {
    H1D_jetPt_mcpFolded = (TH1D*)H1D_jetPt_mcp_control_folded->Rebin(nBinPtJetsRec[iRadius],"H1D_jetPt_mcp_control_folded_recBinning_mcpFoldedWithFluctuationsThenUnfolded"+partialUniqueSpecifier, ptBinsJetsRec[iRadius]);
  } else {
    H1D_jetPt_mcpFolded = (TH1D*)H1D_jetPt_mcp_control_folded->Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAtEndAndEvtNorm_Clone"+partialUniqueSpecifier);
  }
  // bool divideSuccessFakes;
  // TH1D* H1D_jetFakes;
  // if (applyFakes) { // for useManualRespMatrixSettingMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
  //   if (!useFineBinningTest){ 
  //     divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
  //   } else {
  //     divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
  //   }
  //   if (divideSuccessFakes){
  //     H1D_jetPt_mcp_control_recBinning->Divide(H1D_jetFakes);
  //   } else {
  //     cout << "################## H1D_jetPt_mcp_control_recBinning->Divide(H1D_jetFakes) failed because Get_Pt_JetFakes() FAILED!!!!! ##################" << endl;
  //   }
  // }

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpFolded);
  }

  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm()" << endl;};
}
void Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd()" << endl;};

  Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded, iDataset, iRadius, options);

  if (normaliseDistribsInComparisonPlots) {
    if (mcIsWeighted) {
      NormaliseRawHistToInputN(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));//always MCfile_generalresponse here because mcpFolded test doesn't use MC info for the unfolding: we only unfold with fluct matrix from data
    } else {
      NormaliseRawHistToInputN(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework_gen(file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));//always MCfile_generalresponse here because mcpFolded test doesn't use MC info for the unfolding: we only unfold with fluct matrix from data
    }
  }
  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd()" << endl;};
}
void Get_Pt_spectrum_mcpFoldedWithFluctuations(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
  if (showFunctionInAndOutLog) {cout << "--- IN  Get_Pt_spectrum_mcpFoldedWithFluctuations()" << endl;};

  Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpFolded);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcpFolded, deltaJetEta[iRadius]);
  if (showFunctionInAndOutLog) {cout << "--- OUT Get_Pt_spectrum_mcpFoldedWithFluctuations()" << endl;};
}

#endif