#ifndef ANALYSIS_UTILITIES_H
#define ANALYSIS_UTILITIES_H

// Analysis Utilities
void NormaliseYieldToNEntries(TH1D* histogram);
void NormaliseRawHistToNEntries(TH1D* histogram);
void NormaliseAsYieldToInputN(TH1D* histogram, double nEvents);
void NormaliseRawHistToInputN(TH1D* histogram, double nEvents);
void NormaliseYieldToIntegral(TH1D* histogram);
void NormaliseRawHistToIntegral(TH1D* histogram);
long int GetNEvents(TFile* file_O2Analysis, TString analysisWorkflow);
long int GetNEventsGen(TFile* file_O2Analysis);
long int GetNEventsSelected_JetFramework(TFile* file_O2Analysis, const char analysisWorkflow[], bool collHistIsObsolete, bool isDerived);
double GetNEventsSelected_JetFramework_weighted(TFile* file_O2Analysis, const char analysisWorkflow[], bool collHistIsObsolete);
long int GetNEventsSelected_JetFramework_gen(TFile* file_O2Analysis, const char analysisWorkflow[], bool mcCollHistIsObsolete);
double GetNEventsSelected_JetFramework_gen_weighted(TFile* file_O2Analysis, const char analysisWorkflow[], bool mcCollHistIsObsolete);
long int GetNEventsSelected_TrackEffWorkflow(TFile* file_O2Analysis, const char analysisWorkflow[]);
long int GetNEventsSelected_TrackEffWorkflow_gen(TFile* file_O2Analysis, const char analysisWorkflow[]);
double GetNEventsSelected_TrackEffWorkflow_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]);
double GetNEventsSelected_TrackEffWorkflow_gen_weighted(TFile* file_O2Analysis, const char analysisWorkflow[]);
long int GetNEventsSel8Centrality(TFile* file_O2Analysis, float centralityLow, float centralityHigh);
long int GetNEventsSelectedCentrality_JetFramework(TFile* file_O2Analysis, float centralityLow, float centralityHigh, const char wagonId[]);

void TransformRawHistToYield(TH1D* histogram);
void TransformYieldToEtaDifferentialYield(TH1D* histogram, double deltaEta);

#endif