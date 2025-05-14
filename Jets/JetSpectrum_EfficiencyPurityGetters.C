#ifndef JETSPECTRUM_EFFICIENCYPURITYGETTERS_C
#define JETSPECTRUM_EFFICIENCYPURITYGETTERS_C

#include "JetSpectrum_EfficiencyPurityGetters.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_Unfolding.h"
#include "./JetSpectrum_Unfolding.C"

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Efficiency functions //////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void  Get_ResponseMatrix_Pt_KinematicEffiency(TH1D* &H1D_kinematicEfficiency, TH2D* H2D_jetPtResponseMatrix_fineBinning, TString name_H1D_kinematicEfficiency, int iRadius){
  // assumes the response matrix has the fine binning, and will get the kinematic efficiency for rec axis binning equals to ptBinsJetsRec
  // cout << "Get_ResponseMatrix_Pt_KinematicEffiency" << endl; 
  
  int ibinRec_min = H2D_jetPtResponseMatrix_fineBinning->GetXaxis()->FindBin(ptBinsJetsRec[iRadius][0]);
  int ibinRec_max = H2D_jetPtResponseMatrix_fineBinning->GetXaxis()->FindBin(ptBinsJetsRec[iRadius][nBinPtJetsRec[iRadius]])-1;

  cout << "ibinRec_min = " << ibinRec_min << ", ibinRec_max = " << ibinRec_max << ", ptmin = " << ptBinsJetsRec[iRadius][0] << ", ptmax = " << ptBinsJetsRec[iRadius][nBinPtJetsRec[iRadius]] << endl;

  TH1D* H1D_kinematicEfficiency_preRebin = H2D_jetPtResponseMatrix_fineBinning->ProjectionY("H1D_kinematicEfficiency_preRebin"+name_H1D_kinematicEfficiency, ibinRec_min, ibinRec_max, "e");
  
  H1D_kinematicEfficiency = (TH1D*)H1D_kinematicEfficiency_preRebin->Rebin(nBinPtJetsGen[iRadius],"H1D_kinematicEfficiency"+name_H1D_kinematicEfficiency+RadiusLegend[iRadius], ptBinsJetsGen[iRadius]);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TString* pdfNameTest = new TString("kinematicEfficiency_TestByWidthNorm");
  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));
    TH1D* H1D_jetKinematicEfficiencyTest;
  H1D_jetKinematicEfficiencyTest = (TH1D*)H1D_kinematicEfficiency->Clone("H1D_jetKinematicEfficiencyTestbeforWidthNorm");
  TransformRawHistToYield(H1D_jetKinematicEfficiencyTest); // errors will be falsly calculated since division by width is different of division by the sum of the pdf inside of the bins that have errors after the sum 
  Draw_TH1_Histogram(H1D_jetKinematicEfficiencyTest, textContext, pdfNameTest, texPtJetGenX, texJetKinematicEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double integralOfResponse_iBinGen, integralOfResponse_iBinGen_error;
  double binContent, binError, binErrorA, binErrorB;
  for(int iBinGen = 1; iBinGen <= nBinPtJetsGen[iRadius]; iBinGen++){
    int ibinGen_low = H2D_jetPtResponseMatrix_fineBinning->GetYaxis()->FindBin(ptBinsJetsGen[iRadius][iBinGen-1]);
    int ibinGen_high = H2D_jetPtResponseMatrix_fineBinning->GetYaxis()->FindBin(ptBinsJetsGen[iRadius][iBinGen])-1;
    integralOfResponse_iBinGen = H2D_jetPtResponseMatrix_fineBinning->IntegralAndError(1, nBinPtJetsFine[iRadius], ibinGen_low, ibinGen_high, integralOfResponse_iBinGen_error);

    H1D_kinematicEfficiency->GetBinContent(iBinGen) == 0 ? binErrorB = 0 : binErrorB = H1D_kinematicEfficiency->GetBinError(iBinGen)*H1D_kinematicEfficiency->GetBinError(iBinGen) / (H1D_kinematicEfficiency->GetBinContent(iBinGen)*H1D_kinematicEfficiency->GetBinContent(iBinGen));
    integralOfResponse_iBinGen == 0                                    ? binErrorA = 0 : binErrorA = integralOfResponse_iBinGen_error*integralOfResponse_iBinGen_error / (integralOfResponse_iBinGen*integralOfResponse_iBinGen);
    integralOfResponse_iBinGen == 0                                    ? binContent = 0 : binContent = H1D_kinematicEfficiency->GetBinContent(iBinGen) *1./integralOfResponse_iBinGen; // do I really give the value 0 if denominator is 0 ? 


    // H1D_kinematicEfficiency->SetBinContent(iBinGen, H1D_kinematicEfficiency->GetBinContent(iBinGen) * 1./H2D_jetPtResponseMatrix_fineBinning->Integral( 0, -1, ibinGen_low, ibinGen_high));
    // cout << "H1D_kinematicEfficiency_numerator(" << iBinGen << ") = " << H1D_kinematicEfficiency->GetBinContent(iBinGen) << ", error= = " << H1D_kinematicEfficiency->GetBinError(iBinGen) << ", ibinGen_low = " << ibinGen_low << ", ibinGen_high = " << ibinGen_high << "Integral divider = " << integralOfResponse_iBinGen << "Integral error = " << integralOfResponse_iBinGen_error << endl;
    H1D_kinematicEfficiency->SetBinContent(iBinGen, H1D_kinematicEfficiency->GetBinContent(iBinGen) * 1./integralOfResponse_iBinGen);
    H1D_kinematicEfficiency->SetBinError(iBinGen, sqrt(binContent*binContent * (binErrorA + binErrorB))); // sigma(A/B)2 / (A/B) = sigma(A)2 /A2 + sigma(B)2 /B2
    // cout << "H1D_kinematicEfficiency(" << iBinGen << ") = " << binContent << ", error= = " << H1D_kinematicEfficiency->GetBinError(iBinGen) << endl;

  }
}




bool  Get_Pt_JetEfficiency(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, false, options);
  Get_Pt_spectrum_mcpMatched_genBinning(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetEfficiency, H1D_jetPt_mcp, 1., 1., "b");
  if (!divideSuccess){
    cout << "################## Get_Pt_JetEfficiency FAILED!!!!! ##################" << endl;
  }
  if (smoothenEfficiency) {
    for(int i = 2; i <= H1D_jetEfficiency->GetNbinsX(); i++){
      if ((H1D_jetEfficiency->GetBinContent(i) - H1D_jetEfficiency->GetBinContent(i-1)) < -0.005) {
        H1D_jetEfficiency->SetBinContent(i, H1D_jetEfficiency->GetBinContent(i-1));
        H1D_jetEfficiency->SetBinError(i, H1D_jetEfficiency->GetBinError(i-1));
      }
    }
  }
  return divideSuccess;
}
bool  Get_Pt_JetEfficiency_fineBinning(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcpMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, false, options);
  Get_Pt_spectrum_mcpMatched_fineBinning(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  H1D_jetEfficiency = (TH1D*)H1D_jetPt_mcpMatched->Clone("H1D_jetEfficiency_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetEfficiency->Divide(H1D_jetEfficiency, H1D_jetPt_mcp, 1., 1., "b");
  if (!divideSuccess){
    cout << "################## Get_Pt_JetEfficiency FAILED!!!!! ##################" << endl;
  }
  if (smoothenEfficiency) {
    for(int i = 2; i <= H1D_jetEfficiency->GetNbinsX(); i++){
      if ((H1D_jetEfficiency->GetBinContent(i) - H1D_jetEfficiency->GetBinContent(i-1)) < -0.005) {
        H1D_jetEfficiency->SetBinContent(i, H1D_jetEfficiency->GetBinContent(i-1));
        H1D_jetEfficiency->SetBinError(i, H1D_jetEfficiency->GetBinError(i-1));
      }
    }
  }
  return divideSuccess;
}

bool Get_Pt_JetFakes(TH1D* &H1D_jetFakes, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcd;
  TH1D* H1D_jetPt_mcdMatched;
  bool divideSuccess;


  Get_Pt_spectrum_mcd_recBinning(H1D_jetPt_mcd, iDataset, iRadius, options);
  Get_Pt_spectrum_mcdMatched_recBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);


  H1D_jetFakes = (TH1D*)H1D_jetPt_mcdMatched->Clone("H1D_jetFakes"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetFakes->Divide(H1D_jetFakes, H1D_jetPt_mcd, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency (purity similar to efficiency)
  if (!divideSuccess){
    cout << "################## Get_Pt_JetFakes FAILED!!!!! ##################" << endl;
  }

  return divideSuccess;
}
bool Get_Pt_JetFakes_fineBinning(TH1D* &H1D_jetFakes, int iDataset, int iRadius, std::string options){
  TH1D* H1D_jetPt_mcd;
  TH1D* H1D_jetPt_mcdMatched;
  bool divideSuccess;

  Get_Pt_spectrum_mcd_fineBinning(H1D_jetPt_mcd, iDataset, iRadius, options);
  Get_Pt_spectrum_mcdMatched_fineBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  H1D_jetFakes = (TH1D*)H1D_jetPt_mcdMatched->Clone("H1D_jetFakes_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius]));
  divideSuccess = H1D_jetFakes->Divide(H1D_jetFakes, H1D_jetPt_mcd, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency (purity similar to efficiency)
  if (!divideSuccess){
    cout << "################## Get_Pt_JetFakes FAILED!!!!! ##################" << endl;
  }
  return divideSuccess;
}

#endif