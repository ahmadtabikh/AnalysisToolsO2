#ifndef ANALYSIS_UTILITIES_C
#define ANALYSIS_UTILITIES_C

#include "AnalysisUtilities.h"
#include "../Settings/GlobalSettings.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// Analysis Utilities ////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NormaliseYieldToNEntries(TH1D* histogram) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./histogram->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
  cout << "takes overflow/underflow into account; am I sure that is what I want?" << endl;
}
void NormaliseRawHistToNEntries(TH1D* histogram) {
  histogram->Scale(1./histogram->GetEntries(),"");
  cout << "takes overflow/underflow into account; am I sure that is what I want?" << endl;
}

void NormaliseYieldToIntegral(TH1D* histogram) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./histogram->Integral("width"),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}
void NormaliseRawHistToIntegral(TH1D* histogram) { 
  histogram->Scale(1./histogram->Integral("width"),"");
}


void NormaliseAsYieldToInputN(TH1D* histogram, double inputN) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./inputN,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}
void NormaliseRawHistToInputN(TH1D* histogram, double inputN) { 
  histogram->Scale(1./inputN,""); // If option contains "width" the bin contents and errors are divided by the bin width.
}

void TransformRawHistToYield(TH1D* histogram){
  histogram->Scale(1.,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}

void TransformYieldToEtaDifferentialYield(TH1D* histogram, double deltaEta){
  histogram->Scale(1./deltaEta,"");
}


long int GetNEventsSel8(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries(); //this is only sel8 (no sel8Full for example) and doesn't exclude collisions cut by the vertexZ pos cut
}

// long int GetNEventsSelected_JetFramework(TFile* file_O2Analysis, const char analysisWorkflow[], bool collHistIsObsolete = false) {
//   if (collHistIsObsolete) {
//     return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions"))->GetBinContent(3); //this takes the actual selection AND vertexZ into account;
//   } else {
//     return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions"))->GetBinContent(4); //this takes the actual selection AND vertexZ into account;
//   }
// }
long int GetNEventsSelected_JetFramework(TFile* file_O2Analysis,
                                         const char analysisWorkflow[],
                                         bool collHistIsObsolete = false,
                                         bool isDerived = false) {

  if (isDerived) {
    return ((TH1I*)file_O2Analysis->Get("jet-luminosity-calculator/counter"))->GetBinContent(8);
  }

  TH1I* h = (TH1I*)file_O2Analysis->Get((TString)analysisWorkflow + "/h_collisions");

  if (!h) return -1;

  return collHistIsObsolete ? h->GetBinContent(3) : h->GetBinContent(4);
}

double GetNEventsSelected_JetFramework_weighted(TFile* file_O2Analysis, const char analysisWorkflow[], bool collHistIsObsolete = false) {
  if (collHistIsObsolete) {
    return ((TH1F*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions_weighted"))->GetBinContent(3); //this takes the actual selection AND vertexZ into account;
  } else {
    return ((TH1F*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions_weighted"))->GetBinContent(4); //this takes the actual selection AND vertexZ into account;
  }
}

long int GetNEventsSelected_JetFramework_gen(TFile* file_O2Analysis, const char analysisWorkflow[], bool mcCollHistIsObsolete = false) {
  int ibinFinalSelection;
  if (mcCollHistIsObsolete) {
    ibinFinalSelection = 4;
    return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mcColl_counts"))->GetBinContent(ibinFinalSelection);
  } else {
    ibinFinalSelection = 6;
    return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mccollisions"))->GetBinContent(ibinFinalSelection);
  }
}
double GetNEventsSelected_JetFramework_gen_weighted(TFile* file_O2Analysis, const char analysisWorkflow[], bool mcCollHistIsObsolete = false) {
  int ibinFinalSelection;
  if (mcCollHistIsObsolete) {
    ibinFinalSelection = 4;
    return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mcColl_counts_weight"))->GetBinContent(ibinFinalSelection);
  } else {
    ibinFinalSelection = 6;
    return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mccollisions_weighted"))->GetBinContent(ibinFinalSelection);
  }
}


long int GetNEventsSelected_TrackEffWorkflow(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions"))->GetBinContent(3); //this should be the actual selection AND take vertexZ into account; 
}
double GetNEventsSelected_TrackEffWorkflow_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1F*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions_weighted"))->GetBinContent(3); //this should be the actual selection AND take vertexZ into account; 
}
long int GetNEventsSelected_TrackEffWorkflow_gen(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mccollisions"))->GetBinContent(2); //this should be the actual selection AND take vertexZ into account;
}
double GetNEventsSelected_TrackEffWorkflow_gen_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1F*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mccollisions_weighted"))->GetBinContent(2); //this should be the actual selection AND take vertexZ into account;
}

long int GetNEventsGen(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries();
}

long int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh) { // probably same issue as getneventssel8
  TH1D* H1D_Centrality_FT0C= (TH1D*)file_O2Analysis->Get("centrality-qa/hCentFT0C");
  int iBinCent_low = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityLow + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityHigh - GLOBAL_epsilon);
  return H1D_Centrality_FT0C->Integral(iBinCent_low, iBinCent_high);
}


long int GetNEventsSelectedCentrality_JetFramework(TFile* file_O2Analysis, const char analysisWorkflow[], float centralityLow, float centralityHigh) { // should check it gives the correct number of coll (is posZ taken into account, fullsel8 etc)
  TH2D* H2D_Centrality_FT0C= (TH2D*)file_O2Analysis->Get((TString)analysisWorkflow+"/h2_centrality_collisions");
  TH1D* H1D_Centrality_FT0C= (TH1D*)H2D_Centrality_FT0C->ProjectionX("H1D_Centrality_FT0C", 2, 2);

  int iBinCent_low = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityLow + GLOBAL_epsilon);
  int iBinCent_high = H1D_Centrality_FT0C->GetXaxis()->FindBin(centralityHigh - GLOBAL_epsilon);
  return H1D_Centrality_FT0C->Integral(iBinCent_low, iBinCent_high);
}











#endif
