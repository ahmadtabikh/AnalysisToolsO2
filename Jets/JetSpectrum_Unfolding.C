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


std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  //for now makes the assumption gen and rec have the same pT binning

  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* measured = (TH1D*)measuredInput->Clone("measured_Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);;
  TH1D* mcp;
  TH1D* mcd;

  if (!normGenAndMeasByNEvtsBeforeUnfolding) {
    // Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); now fed as input to function
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      // Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(mcd, iDataset, iRadius, options);
    }
  } else{
    // Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(mcp, iDataset, iRadius, false, options);
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(mcd, iDataset, iRadius, options);

    if (useFineBinningTest) {
      // Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
      Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(mcp, iDataset, iRadius, false, options);
      Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(mcd, iDataset, iRadius, options);
    }
  }

  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (useManualRespMatrixSettingMethod) {
    if (applyFakes && (options.find("noPurity") == std::string::npos)) { // if applyFakes AND if option noPurity has not been found in options; necessary check for the mcp-folded unfolding test as we don't want the fake correction
      if (!useFineBinningTest){ 
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
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
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  
  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);
  // compute matrixFluctuations times matrixDetector

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
  }
  
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

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
    
    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

    TString* pdfNamePost = new TString("ResponseMatrices/responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfNamePost_logz = new TString("ResponseMatrices/responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postWeighting).Clone("responseMatrix_combined_postWeighting"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfNamePost, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfNamePost_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "RooUnfoldResponse setting - start" << endl;
  RooUnfoldResponse* response;
  if (useManualRespMatrixSettingMethod) {

    // // based on Marta's work: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JEJetSpectrumUnfolding

    if (useFineBinningTest) {
      response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning); // measuredInput and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
    } else {
      response = new RooUnfoldResponse(0, 0, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined); // measuredInput and mcp_rebinned are here to take inneficiencies and fakes into account; or is it really what's happening? 'Alternatively, the response matrix can be constructed from a pre-existing TH2D 2-dimensional histogram (with truth and measured distribution TH1D histograms for normalisation).' from https://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html, and I'm already normalising so maybe tehre's no need for more normalisation
    }

    cout << "Get_Pt_spectrum_unfolded(): should I use response->UseOverflow() ? using it gives a ratio unfolded/mcp much higher than without using it" << endl;
    // response->UseOverflow();

  } else {
    TH2D *Respt;
    if (useFineBinningTest){
      Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
    } else {
      Respt = H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
    }
    TH1D *mcdMatched = (TH1D*)Respt->ProjectionX("Respt_projX", 1, Respt->GetNbinsY());
    TH1D *mcpMatched = (TH1D*)Respt->ProjectionY("Respt_projY", 1, Respt->GetNbinsX());
    TH1D *fake = (TH1D*)mcd->Clone();
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        fake->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        fake->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      }
    }
    for (auto i = 1; i <= fake->GetNbinsX(); i++) { 
      cout << "Fake i = " << i << " content: "<< fake->GetBinContent(i) << " - " << mcdMatched->GetBinContent(i) << " = " << fake->GetBinContent(i) - mcdMatched->GetBinContent(i) << " ------ error of mcd = " << fake->GetBinError(i) << ", error of mcdMatched = " << mcdMatched->GetBinError(i) << endl;
    }
    fake->Add(mcdMatched, -1);
    TH1D *miss = (TH1D*)mcp->Clone(); // for each bin of gen-level jet distribution, how many are not matched to a rec-level jet
    if (normDetRespByNEvts) {
      if (mcIsWeighted) {
        miss->Scale(1./GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
      } else {
        miss->Scale(1./GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
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
    // unfoldSvd->SetNtoysSVD(1000);
    hist_unfold = (TH1D*)(unfoldSvd->Hreco());
    // /////////////////////////////////////////////////////////////////////// Print Error /////////////////////////////////////////////////////////////////////
    // TVectorD errVec = unfoldSvd->ErecoV(RooUnfold::kCovToy);
    // TH1D* hErrorToys = (TH1D*)hist_unfold->Clone("hErrorToys");
    // hErrorToys->Reset("ICES"); // clear everything but keep structure
    // for (int i = 1; i <= hist_unfold->GetNbinsX(); ++i) {
    //   double sigma_i = errVec[i-1];
    //   double M_i = hist_unfold->GetBinContent(i);
    //   double relStd = (M_i != 0.0) ? sigma_i / M_i : 0.0;
    //   hErrorToys->SetBinContent(i, relStd);
    //   hErrorToys->SetBinError(i, 0.0); // no error on the stddev itself
    // }

    // TString* pdfName_hErrorToys = new TString("errors on nominal unfolding hErrorToys"+partialUniqueSpecifier);
    // TString* PtUnfNominal = new TString("#it{p}_{T}^{unf, nominal} (GeV/#it{c})");
    // TString* ErecoVNominal = new TString("ErecoV rel errors unfolding nominal ");
    // TString textContext1(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));
    // Draw_TH1_Histogram(hErrorToys, textContext1, pdfName_hErrorToys, PtUnfNominal, ErecoVNominal, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
    // cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
    // cout << "Number of toy experiments used: " << unfoldSvd->GetNtoysSVD() << endl;
    // cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
    // cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
    
    
    unfold = unfoldSvd;

    // plot svd d distribution
    TH1D* H1D_D = tsvdUnfold->GetD();
    TString inputUnfoldingName = (options.find("inputIsMCPFoldedTest") != std::string::npos) ? "_mcpFoldedTestInput" : "";
    TString* pdfName_regparam = new TString("Svd_regularisationd_distribution_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+(TString)mergingPrior+"_"+(TString)unfoldingPrior+inputUnfoldingName);
    TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), ""));
    std::array<std::array<float, 2>, 2> drawnWindowSvdParam = {{{0, 30}, {0.01, 10000}}}; // {{xmin, xmax}, {ymin, ymax}}
    // Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdDvector, texSvdK, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");
    Draw_TH1_Histogram(H1D_D, textContext, pdfName_regparam, texSvdK, texSvdDvector, texCollisionDataInfo, drawnWindowSvdParam, legendPlacementAuto, contextPlacementAuto, "logy");

  } else if (options.find("Bayes") != std::string::npos) {
    unfoldParameter = unfoldParameterInput;
    hist_unfold = (TH1D*)(unfoldBayes->Hreco());
    unfold = unfoldBayes;
  } 


  H1D_jetPt_unfolded = (TH1D*)hist_unfold->Clone("H1D_jetPt_unfolded"+partialUniqueSpecifier);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////draw measured input and unfolded(after eff division) before normalisation by evt and width////////////////////////////////////////////////
  // TH1D** H1D_jetPt_unfolded_measuredComp = new TH1D*[2];
  // Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning_test, iDataset, iRadius, options);
  // H1D_jetPt_unfolded_measuredComp[0] = H1D_jetPt_measured_genBinning;
  // H1D_jetPt_unfolded_measuredComp[1] = H1D_jetPt_unfolded;
  // TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
  // TString* Yaxislabel = new TString("counts");
  // TString Measured_unfolded_Legend[2] = {"unfolded after eff division", "measured"};
  // TString* pdfName_measured_unfolded_comp = new TString("_measured_unfolded_comp_counts");
  // Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, Measured_unfolded_Legend, 2, textContext, pdfName_measured_unfolded_comp, texPtX, Yaxislabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    
    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");


    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_postUnfolding).Clone("responseMatrix_combined_postUnfolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }


  // if (doWidthScalingEarly) {
  //   TransformRawHistToYield(H1D_jetPt_unfolded);
  // }

  std::pair<int, RooUnfold*> unfoldInfo(unfoldParameter, unfold);
  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);
  cout << "Get_Pt_spectrum_unfolded_preWidthScalingAtEnd test 1" << endl;
  if (normaliseDistribsAfterUnfolding){
    if (!controlMC && options.find("inputIsMC") == std::string::npos) { // if option controlMC is false, and if inputIsMC has not been found in options; necessary check for the mcp-folded unfolding test as we want to normalise by the number of events in the MC file
      NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
   cout << "Get_Pt_spectrum_unfolded_preWidthScalingAtEnd test 2" << endl;

  return unfoldInfo;
}
std::pair<int, RooUnfold*> Get_Pt_spectrum_unfolded(TH1D* &H1D_jetPt_unfolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  std::pair<int, RooUnfold*> unfoldInfo = Get_Pt_spectrum_unfolded_preWidthScalingAtEnd(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfolded);
  }
  return unfoldInfo;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Spectrum Refolding tests //////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).second;
  // Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, iDataset, iRadius, unfoldParameterInput, options);

  // cout << "((((((((((((((((()))))))))))))))))" << endl;
  // cout << "REFOLDING TEST: " << endl;
  // for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
  //   cout << "H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << ", error = "<< H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  // }

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation
  if (useFineBinningTest) {
    ReweightResponseMatrixWithPrior_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, iDataset, iRadius, options);
  }

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);
  

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

    TString priorInfo = (TString)(TString)mergingPrior+"-"+(TString)unfoldingPrior;

    TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]);
    TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_logz");

    TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
    TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

    // the matrix natural visualisation is actually the transpose of histograms' visualisations we use
    TH2D* MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_duringRefolding).Clone("responseMatrix_combined_duringRefolding"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);

    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
    Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, texPtJetGenX, texPtJetRecX, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  }

  cout << "((((((((((((((((()))))))))))))))))" << endl;
  cout << "REFOLDING TEST: pre efficiency" << endl;
  for(int iBinX = 0; iBinX <= H1D_jetPt_unfolded->GetNbinsX()+1; iBinX++){
    cout << "  H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinContent(iBinX) << endl;
    cout << "e_H1D_jetPt_unfolded(" << iBinX << ") = " << H1D_jetPt_unfolded->GetBinError(iBinX) << endl;
  }
  

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
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
    if (useFineBinningTest) {
      refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
    } else {
      refoldingResponseMatrix = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("refoldingResponseMatrix"+partialUniqueSpecifier);
    }
  } else {
    refoldingResponseMatrix = (TH2D*)unfold->response()->HresponseNoOverflow();
    // refoldingResponseMatrix = (TH2D*)unfold->response()->Hresponse();
    cout << "refold why use HresponseNoOverflow instead of Hresponse?" << endl;
  }

  if (normaliseRespYSliceForRefold){
    NormaliseYSlicesToOne(refoldingResponseMatrix);
    cout << "###################################################################################################" << endl;
    cout << "###################################################################################################" << endl;
    cout << "###################################################################################################" << endl;
    cout << "################  X axis: " << refoldingResponseMatrix->GetNbinsX() << "################" << endl;
    cout << "################  Y axis: " << refoldingResponseMatrix->GetNbinsY() << "################" << endl;
    cout << "###################################################################################################" << endl;
    cout << "###################################################################################################" << endl;
    cout << "###################################################################################################" << endl;

  }
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1((TH2D*)unfold->response()->Hresponse(), H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier); //gives the same
  // H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier); //gives the same
  cout << "--UnfoldedJetPt in RooUnfold Refold Method--" << endl;
  H1D_jetPt_unfoldedThenRefolded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_unfolded).Clone("Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier); //gives the same


  bool divideSuccessFakes;
  TH1D* H1D_jetFakes;
  if (applyFakes) { // for useManualRespMatrixSettingMethod set to false, it will give a slightly different result as Get_Pt_JetFakes doesn't get fakes from the same histogram as the one used to encode fakes in the response matrix
    if (!useFineBinningTest){ 
      divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
    } else {
      divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
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
void Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);

  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}


void Get_Pt_spectrum_dataUnfoldedThenRefolded(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_preWidthScalingAtEnd(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}


void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  // ApplyToTruth function doesn't apply errors on folding

  // Matches exactly with the manual method IF NO PRIOR
  // if I have a non flat prior, then the roounfold method gives me a good closure test, but not the manual method!

  TH1D* H1D_jetPt_unfolded;
  // TH1D* H1D_jetPt_raw[nRadius];
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);


  RooUnfold* unfold = Get_Pt_spectrum_unfolded_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).second;

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning;
  TH1D* H1D_kinematicEfficiency;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius);
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius);

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_fineBinning(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency, H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, partialUniqueSpecifier, iRadius); // I want the efficiency before the reweighting and normalisation

  Get_PtResponseMatrix_DetectorAndFluctuationsCombined(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  if (!useFineBinningTest){ 
    divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
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
        divideSuccessFakes = Get_Pt_JetFakes(H1D_jetFakes, iDataset, iRadius, options);
      } else {
        divideSuccessFakes = Get_Pt_JetFakes_fineBinning(H1D_jetFakes, iDataset, iRadius, options);
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
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);

  if (normaliseDistribsAfterUnfolding){
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_unfoldedThenRefolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod(TH1D* &H1D_jetPt_unfoldedThenRefolded, TH1D* &measuredInput, int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  Get_Pt_spectrum_dataUnfoldedThenRefolded_RooUnfoldMethod_preWidthScalingAtEnd(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_unfoldedThenRefolded);
  }
}





void Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* H1D_jetPt_mcp_control;
  // Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_control, iDataset, iRadius, true, options);
  // Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_control, iDataset, iRadius, true, options);
  // if (useFineBinningTest) {
    Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_control, iDataset, iRadius, true, options);
  // }


  // applying efficiency loss and fake enhancement
  bool divideSuccessEff;
  TH1D* H1D_jetEfficiency;
  // if (!useFineBinningTest){ 
  //   divideSuccessEff = Get_Pt_JetEfficiency(H1D_jetEfficiency, iDataset, iRadius, options);
  // } else {
    divideSuccessEff = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency, iDataset, iRadius, options);
  // }

  // we do not apply efficiencies on this folding exercise with fluctuations
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
  Get_PtResponseMatrix_Fluctuations(refoldingResponseMatrix, iDataset, iRadius);
  // cout << "Integral Line fluct matrix: " << refoldingResponseMatrix->Integral(1, refoldingResponseMatrix->GetNbinsX(), 10, 10) << endl;

  TH1D* H1D_jetPt_mcp_control_folded = (TH1D*)GetMatrixVectorProductTH2xTH1(refoldingResponseMatrix, H1D_jetPt_mcp_control).Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);
  // TH1D* H1D_jetPt_mcp_control_folded = (TH1D*)H1D_jetPt_mcp_control->Clone("Get_Pt_spectrum_mcpFoldedWithFluctuationsThenUnfolded_preWidthScalingAtEndAndEvtNorm"+partialUniqueSpecifier);

  if (!useFineBinningTest) {
    H1D_jetPt_mcpFolded = (TH1D*)H1D_jetPt_mcp_control_folded->Rebin(nBinPtJetsRec[iRadius],"H1D_jetPt_mcp_control_folded_recBinning_mcpFoldedWithFluctuationsThenUnfolded"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsRec[iRadius]);
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
}
// void Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
//   Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded, iDataset, iRadius, options);

//   if (normaliseDistribsAfterUnfolding){
//     if (!controlMC) {
//       NormaliseRawHistToNEvents(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
//     } else {
//       if (!mcIsWeighted) {
//         NormaliseRawHistToNEvents(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
//       } else {
//         NormaliseRawHistToNEvents(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
//       }
//     }
//   }
// }
// void Get_Pt_spectrum_mcpFoldedWithFluctuations(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
//   Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded, iDataset, iRadius, options);

//   if (doWidthScalingAtEnd) {
//     TransformRawHistToYield(H1D_jetPt_mcpFolded);
//   }
// }
void Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpFolded, GetNEventsSelected_JetFramework_gen( file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
    }
  }
}
void Get_Pt_spectrum_mcpFoldedWithFluctuations(TH1D* &H1D_jetPt_mcpFolded, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpFolded);
  }
}

#endif