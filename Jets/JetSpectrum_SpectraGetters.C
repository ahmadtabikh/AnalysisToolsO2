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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////// Spectrum getting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_recBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();
  
  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_bkgCorrected_rebinned_recBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}


void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_bkgCorrected_rebinned_genBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInData) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  if (!controlMC) {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflowData+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_bkgCorrected_fineBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
  }
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_bkgCorrected_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  if (!fcontrolMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_genBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
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

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_genBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  if (!fcontrolMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_fineBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
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

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcp_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_part_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt_part";
  }

  if (!fcontrolMC){
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  } else {
    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_ppSimDetectorEffect_unfoldingControl->Get(analysisWorkflow_unfoldingControl+"/"+histogramName))->Clone("Get_Pt_spectrum_mcp_recBinning_unfoldingControl"+Datasets[iDataset]+DatasetsNames[iDataset]);
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

  H1D_jetPt_mcp = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcp_rebinned_recBinning"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcd_rebinned_fineBinning_"+Datasets[iDataset]+DatasetsNames[iDataset]+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

  H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcd_rebinned_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH1D* H1D_jetPt_defaultBin;

  TString histogramName = "";
  if (doBkgSubtractionInMC) {
    histogramName = "h_jet_pt_rhoareasubtracted";
  } else {
    histogramName = "h_jet_pt";
  }

    H1D_jetPt_defaultBin = (TH1D*)((TH1D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/"+histogramName))->Clone("Get_Pt_spectrum_mcd_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
  H1D_jetPt_defaultBin->Sumw2();

  H1D_jetPt_mcd = (TH1D*)H1D_jetPt_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcd_rebinned_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcp;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcpMatched_defaultBin;
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      } else {
        H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      }
    }
    H2D_jetPtMcdjetPtMcp->Sumw2();

    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcp->ProjectionY("jetPt_mcpMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcp->GetNbinsX(), "e");

  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionY("jetPt_mcpMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcpMatched_genBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcp;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcpMatched_defaultBin;

  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      } else {
        H2D_jetPtMcdjetPtMcp = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcpetaconstraint"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      }
    }
    H2D_jetPtMcdjetPtMcp->Sumw2();

    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcp->ProjectionY("jetPt_mcpMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcp->GetNbinsX(), "e");
  
  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcpMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcpMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionY("jetPt_mcpMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcpMatched = (TH1D*)H1D_jetPt_mcpMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcpMatched_fineBinning_rebinned"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;
  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      } else {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      }
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");

  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_genBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_genBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcdMatched_genBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      } else {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      }
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  
  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_recBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_recBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsRec[iRadius],"jetPt_mcdMatched_recBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsRec[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, __attribute__ ((unused)) std::string options) {
  TH2D* H2D_jetPtMcdjetPtMcd;
  TH3D* H3D_jetRjetPtTagjetPtBase;

  TH1D* H1D_jetPt_mcdMatched_defaultBin;

  if (analysisWorkflowMC.Contains("jet-spectra-charged") == true) {
    if (doBkgSubtractionInMC) {
      H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_rhoareasubtracted"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
    } else {
      if (etaCutOnMatchedJetsIsObsoleteVersion == true) {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      } else {
        H2D_jetPtMcdjetPtMcd = (TH2D*)((TH2D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h2_jet_pt_mcd_jet_pt_mcp_matchedgeo_mcdetaconstraint"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);
      }
    }
    H2D_jetPtMcdjetPtMcd->Sumw2();

    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H2D_jetPtMcdjetPtMcd->ProjectionX("jetPt_mcdMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], 1, H2D_jetPtMcdjetPtMcd->GetNbinsY(), "e");
  
  } else if (analysisWorkflowMC.Contains("jet-finder-charged-qa") == true) {
    H3D_jetRjetPtTagjetPtBase = (TH3D*)((TH3D*)file_O2Analysis_MCfileForMatrix->Get(analysisWorkflowMC+"/h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"))->Clone("Get_Pt_spectrum_mcdMatched_fineBinning"+Datasets[iDataset]+DatasetsNames[iDataset]);// base is mcd in jetfinderQA as of 06/2024, thus tag is mcp, and so hist is (x=r, y=mcp, z=mcd)
    H3D_jetRjetPtTagjetPtBase->Sumw2();

    int ibinJetRadius = H3D_jetRjetPtTagjetPtBase->GetXaxis()->FindBin(arrayRadius[iRadius]+GLOBAL_epsilon);
    H1D_jetPt_mcdMatched_defaultBin = (TH1D*)H3D_jetRjetPtTagjetPtBase->ProjectionZ("jetPt_mcdMatched_fineBinning_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ibinJetRadius, ibinJetRadius, 1, H3D_jetRjetPtTagjetPtBase->GetNbinsZ(), "e");
  } else {
    cout << "analysisWorkflowMC name is wrong, need to be corrected" << endl;
    abort();
  }

  H1D_jetPt_mcdMatched = (TH1D*)H1D_jetPt_mcdMatched_defaultBin->Rebin(nBinPtJetsFine[iRadius],"jetPt_mcdMatched_fineBinning_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsFine[iRadius]);

  if (doWidthScalingEarly) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options);


  if (normaliseDistribsBeforeUnfolding) {
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_recBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(H1D_jetPt, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options);


  if (normaliseDistribsBeforeUnfolding) {
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}

void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(H1D_jetPt, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (!controlMC) {
      NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_list[iDataset], analysisWorkflowData));
    } else {
      if (mcIsWeighted) {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt, GetNEventsSelected_JetFramework(file_O2Analysis_ppSimDetectorEffect_unfoldingControl, analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(H1D_jetPt, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt);
  }
}


void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_ppSimDetectorEffect_unfoldingControl[iDataset], analysisWorkflow_unfoldingControl));
      }
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_ppSimDetectorEffect_unfoldingControl[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_ppSimDetectorEffect_unfoldingControl[iDataset], analysisWorkflow_unfoldingControl));
      }
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_ppSimDetectorEffect_unfoldingControl[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcd, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcd_genBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcd_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcd, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcd);
  }
}

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcdMatched_genBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcdMatched_recBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_recBinning_preWidthScalingAtEnd(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework_weighted(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcdMatched, GetNEventsSelected_JetFramework(file_O2Analysis_MCfileForMatrix, analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcdMatched_fineBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcdMatched_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcdMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcdMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcpMatched_genBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEnd(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
    } else {
      NormaliseRawHistToNEvents(H1D_jetPt_mcpMatched, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
    }
  }
}
void Get_Pt_spectrum_mcpMatched_fineBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options) {
  Get_Pt_spectrum_mcpMatched_fineBinning_preWidthScalingAtEnd(H1D_jetPt_mcpMatched, iDataset, iRadius, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcpMatched);
  }
}

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (normaliseDistribsBeforeUnfolding) {
    if (mcIsWeighted) {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen_weighted( file_O2Analysis_ppSimDetectorEffect_unfoldingControl[iDataset], analysisWorkflow_unfoldingControl));
      }
    } else {
      if (!fcontrolMC) {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_MCfileForMatrix[iDataset], analysisWorkflowMC));
      } else {
        NormaliseRawHistToNEvents(H1D_jetPt_mcp, GetNEventsSelected_JetFramework_gen( file_O2Analysis_ppSimDetectorEffect_unfoldingControl[iDataset], analysisWorkflow_unfoldingControl));
      }
    }
  }
}
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, bool fcontrolMC, std::string options) {
  Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEnd(H1D_jetPt_mcp, iDataset, iRadius, fcontrolMC, options);

  if (doWidthScalingAtEnd) {
    TransformRawHistToYield(H1D_jetPt_mcp);
  }
}

#endif