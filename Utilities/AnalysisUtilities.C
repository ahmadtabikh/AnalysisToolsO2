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


void NormaliseYieldToNEvents(TH1D* histogram, double nEvents) { // also takes care of the transformation of raw count to dCount/dQuantity (like dN/dpT), thanks to option "width"
  histogram->Scale(1./nEvents,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}
void NormaliseRawHistToNEvents(TH1D* histogram, double nEvents) { 
  histogram->Scale(1./nEvents,""); // If option contains "width" the bin contents and errors are divided by the bin width.
}

void TransformRawHistToYield(TH1D* histogram){
  histogram->Scale(1.,"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
}

long int GetNEventsSel8(TFile* file_O2Analysis) {
  return ((TH1I*)file_O2Analysis->Get("event-selection-task/hColCounterAcc"))->GetEntries(); //this is only sel8 (no sel8Full for example) and doesn't exclude collisions cut by the vertexZ pos cut
}

long int GetNEventsSelected_JetFramework(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions"))->GetBinContent(3); //this should be the actual selection AND take vertexZ into account;

}
double GetNEventsSelected_JetFramework_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1F*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_collisions_weighted"))->GetBinContent(3); //this should be the actual selection AND take vertexZ into account; 

}

long int GetNEventsSelected_JetFramework_gen(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  if (oldformatNEventsGenInO2Analysis) {
    return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mcColl_counts"))->GetBinContent(6);
  }
  else {
    return ((TH1I*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mccollisions"))->GetBinContent(6);
  }
}

double GetNEventsSelected_JetFramework_gen_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]) {
  return ((TH1F*)file_O2Analysis->Get((TString)analysisWorkflow+"/h_mcColl_counts_weight"))->GetBinContent(6);
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
