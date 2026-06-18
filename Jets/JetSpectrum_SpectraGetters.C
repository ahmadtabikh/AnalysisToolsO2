#ifndef JETSPECTRUM_SPECTRAGETTERS_C
#define JETSPECTRUM_SPECTRAGETTERS_C

#include "JetSpectrum_SpectraGetters.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
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
#include "JetSpectrum_inputs.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////// Spectrum getting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }
  TString histogramNameControlMC = "";
  if (doBkgSubtractionInMC) {
    histogramNameControlMC = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramNameControlMC = "h_jet_pt";
  }

  TString controlMcSpecifier = Form("controlMC%i",controlMC);
  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramNameControlMC))->Clone("Get_Pt_spectrum_bkgCorrected_recBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_bkgCorrected_rebinned_recBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier+RadiusLegend[iRadius], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}


void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }
  TString histogramNameControlMC = "";
  if (doBkgSubtractionInMC) {
    histogramNameControlMC = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramNameControlMC = "h_jet_pt";
  }

  TString controlMcSpecifier = Form("controlMC%i",controlMC);
  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramNameControlMC))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_bkgCorrected_rebinned_genBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier+RadiusLegend[iRadius], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }
  TString histogramNameControlMC = "";
  if (doBkgSubtractionInMC) {
    histogramNameControlMC = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramNameControlMC = "h_jet_pt";
  }

  TString controlMcSpecifier = Form("controlMC%i",controlMC);
  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramNameControlMC))->Clone("Get_Pt_spectrum_bkgCorrected_fineBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_bkgCorrected_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier+RadiusLegend[iRadius], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H1D_jetPt_defaultBin->Sumw2();


  if (smoothenMCP) {
    for(int i = 1; i <= H1D_jetPt_defaultBin->GetNbinsX()-1; i++){
      if ((H1D_jetPt_defaultBin->GetBinContent(i+1) - H1D_jetPt_defaultBin->GetBinContent(i)) > 0.01) {
        H1D_jetPt_defaultBin->SetBinContent(i+1, H1D_jetPt_defaultBin->GetBinContent(i));
        H1D_jetPt_defaultBin->SetBinError(i+1, H1D_jetPt_defaultBin->GetBinError(i));
      }
    }
  }

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_genBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H1D_jetPt_defaultBin->Sumw2();

  if (smoothenMCP) {
    for(int i = 1; i <= H1D_jetPt_defaultBin->GetNbinsX()-1; i++){
      if ((H1D_jetPt_defaultBin->GetBinContent(i+1) - H1D_jetPt_defaultBin->GetBinContent(i)) > 0.01) {
        H1D_jetPt_defaultBin->SetBinContent(i+1, H1D_jetPt_defaultBin->GetBinContent(i));
        H1D_jetPt_defaultBin->SetBinError(i+1, H1D_jetPt_defaultBin->GetBinError(i));
      }
    }
  }

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcp_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H1D_jetPt_defaultBin->Sumw2();

  if (smoothenMCP) {
    for(int i = 1; i <= H1D_jetPt_defaultBin->GetNbinsX()-1; i++){
      if ((H1D_jetPt_defaultBin->GetBinContent(i+1) - H1D_jetPt_defaultBin->GetBinContent(i)) > 0.01) {
        H1D_jetPt_defaultBin->SetBinContent(i+1, H1D_jetPt_defaultBin->GetBinContent(i));
        H1D_jetPt_defaultBin->SetBinError(i+1, H1D_jetPt_defaultBin->GetBinError(i));
      }
    }
  }

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcp_rebinned_recBinning"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }  
  
  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_fineBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_fineBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcd_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_recBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_recBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcd_rebinned_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_genBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H1D_jetPt_defaultBin = (TH1D*)((TH1D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_genBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcd_rebinned_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH2D* H2D_jetPtMcdjetPtMcp;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcpMatched_defaultBin;

  TString histogramName = "";
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted";
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo";
      } else {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint";
      }
    }
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcpMatched_genBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcpMatched_genBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H2D_jetPtMcdjetPtMcp->Sumw2();
  H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcp->ProjectionY("jetPt_mcpMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, 1, H2D_jetPtMcdjetPtMcp->GetNbinsX(), "e");

  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcpMatched_genBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH2D* H2D_jetPtMcdjetPtMcp;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcpMatched_defaultBin;

  TString histogramName = "";
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted";
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo";
      } else {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint";
      }
    }
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H2D_jetPtMcdjetPtMcp->Sumw2();
  H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcp->ProjectionY("jetPt_mcpMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, 1, H2D_jetPtMcdjetPtMcp->GetNbinsX(), "e");
  
  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcpMatched_fineBinning_rebinned"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  TString histogramName = "";
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted";
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo";
      } else {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
      }
    }
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_genBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_genBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H2D_jetPtMcdjetPtMcd->Sumw2();
  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcdMatched_genBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  TString histogramName = "";
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted";
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo";
      } else {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
      }
    }
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_recBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_recBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H2D_jetPtMcdjetPtMcd->Sumw2();
  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  
  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcdMatched_recBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  TString histogramName = "";
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted";
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo";
      } else {
        histogramName = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint";
      }
    }
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  TString controlMcSpecifier = Form("controlMC%i_controlMC_useResponseSplit%i",controlMC, controlMC_useResponseSplit);
  if (!controlMC){
    H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_GeneralResponse[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
  } else {
    if (controlMC_useResponseSplit) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_response[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning_unfoldingControl_response"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    } else {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*) file_O2Analysis_MCfile_UnfoldingControl_input[iDataset]->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning_unfoldingControl_input"+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier);
    }
  }
  H2D_jetPtMcdjetPtMcd->Sumw2();
  H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  
  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcdMatched_fineBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset]+Form("%.1d",iDataset)+controlMcSpecifier, ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC = false) {
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options, controlMC);

  if (!controlMC) {
    NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]));
  } else {        
    if (mcIsWeighted) {
      if (!mcpInput_useMcpCollCountForUnfoldingResultNorm) {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    } else {
      if (!mcpInput_useMcpCollCountForUnfoldingResultNorm) {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_recBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(H1D_jetPt, iDataset, iRadius, options, controlMC);
  } else {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options, controlMC);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC = false) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options, controlMC);

  if (!controlMC) {
    NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]));
  } else {
    if (mcIsWeighted) {
      if (!mcpInput_useMcpCollCountForUnfoldingResultNorm) {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    } else {
      if (!mcpInput_useMcpCollCountForUnfoldingResultNorm) {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(H1D_jetPt, iDataset, iRadius, options, controlMC);
  } else {
    Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options, controlMC);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC = false) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options, controlMC);

  if (!controlMC) {
    NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData, collHistIsObsolete[iDataset], derived_data[iDataset]));
  } else {
    if (mcIsWeighted) {
      if (!mcpInput_useMcpCollCountForUnfoldingResultNorm) {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    } else {
      if (!mcpInput_useMcpCollCountForUnfoldingResultNorm) {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(H1D_jetPt, iDataset, iRadius, options, controlMC);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options, controlMC);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt, deltaJetEta[iRadius]);
}


void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcp, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcp, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      }
    }
  }
}
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcd, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      }
    }
  }
}
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcd, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcd, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      }
    }
  }
}
void Get_Pt_spectrum_mcd_genBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcd, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      }
    }
  }
}
void Get_Pt_spectrum_mcdMatched_genBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcdMatched, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      }
    }
  }
}
void Get_Pt_spectrum_mcdMatched_recBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEnd(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcdMatched, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, collHistIsObsolete[iDataset]));
      }
    }
  }
}
void Get_Pt_spectrum_mcdMatched_fineBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcdMatched, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_mcpMatched_genBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcpMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcpMatched, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_mcpMatched_fineBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcpMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcpMatched, deltaJetEta[iRadius]);
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);

  if (mcIsWeighted) {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted(file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  } else {
    if (!controlMC) {
      NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_GeneralResponse[iDataset], analysisWorkflowMC, mcCollHistIsObsolete));
    } else {
      if (controlMC_useResponseSplit) {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_response[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      } else {
        NormaliseRawHistToInputN(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfile_UnfoldingControl_input[iDataset], analysisWorkflow_unfoldingControl, mcCollHistIsObsolete));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC = false, bool controlMC_useResponseSplit = false) {
  if (normaliseDistribsInComparisonPlots) {
    Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEnd(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  } else {
    Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, options, controlMC, controlMC_useResponseSplit);
  }

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }

  TransformYieldToEtaDifferentialYield(H1D_jetPt_mcp, deltaJetEta[iRadius]);
}

#endif