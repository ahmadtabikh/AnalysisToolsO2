#ifndef AXISTITLES_H
#define AXISTITLES_H

// Commonly used titles

////////////////////////////////////////////
////////////////// x-axis ////////////////// not really a good name for the separation anymore
////////////////////////////////////////////

TString* texRapidity = new TString("y");
TString* texPtX = new TString("#it{p}_{T} (GeV/#it{c})");
TString* texPtLeadX = new TString("#it{p}_{T, lead} (GeV/#it{c})");
TString* texEtaX = new TString("#it{#eta}");
TString* texPhiX = new TString("#it{#phi} (rad)");
TString* texNTracksX = new TString("#it{N}_{tracks in jet}");
TString* texLeadPt = new TString("#it{p}_{T,track}^{leading}");
TString* texJetArea = new TString("#it{A}_{jet}");
TString* texJetNTracks = new TString("#it{N}_{tracks}");
TString* texJetNTrdTracks = new TString("#it{N}_{TRD tracks}");
TString* texJetNTrdTracksRatio = new TString("#frac{#it{N}_{TRD tracks}}{#it{N}_{tracks}}");
TString* texCentrality = new TString("Centrality (%)");
TString* texSelectedMultiplicity = new TString("#it{N}_{tracks}^{selected}");
TString* texLeadJetPt = new TString("#it{p}_{T,lead jet} (GeV/#it{c})");

TString* texPtJetRecX = new TString("#it{p}_{T,jet rec} (GeV/#it{c})");
TString* texPtJetGenX = new TString("#it{p}_{T,jet gen} (GeV/#it{c})");

TString* texPtJetBkgCorrX = new TString("#it{p}_{T,jet rec} - #it{A}_{jet} * #it{#rho} (GeV/#it{c})");
TString* texPtJetBkgFreeX = new TString("#it{p}_{T,jet rec}^{bkg-free} (GeV/#it{c})");
TString* texPtJetRawX = new TString("#it{p}_{T,jet raw} (GeV/#it{c})");

TString* texVx = new TString("part Vx");
TString* texVy = new TString("part Vy");
TString* texVz = new TString("part Vz");

TString* texIUy = new TString("track IU y (cm)");


TString* texSvdK = new TString("k");
TString* texSvdDvector = new TString("d_{k}");


////////////////////////////////////////////
////////////////// y-axis ////////////////// 
////////////////////////////////////////////

TString* texPtDifferentialYield = new TString("1/#it{N}_{ev} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texCount = new TString("count");

TString* texPtMeasured = new TString("#it{p}_{T}^{measured} (GeV/#it{c})");
TString* texPtMC = new TString("#it{p}_{T, gen} (GeV/#it{c})");
TString* texPtRec = new TString("#it{p}_{T, rec} (GeV/#it{c})");
TString* texEtaMC = new TString("#eta_{gen}");
TString* texEtaRec = new TString("#eta_{rec}");
TString* texPhiMC = new TString("#phi_{gen} (rad)");
TString* texPhiRec = new TString("#phi_{rec} (rad)");
TString* texSigmaPt = new TString("#sigma(#it{p}_{T})/#it{p}_{T}");
TString* texSigmaPtMean = new TString("<#sigma(#it{p}_{T})/#it{p}_{T}>");
TString* texSigmaPtMedian = new TString("median(#sigma(#it{p}_{T})/#it{p}_{T})");

TString* texRatio = new TString("ratio");

TString* texRatioDatasets = new TString(*texDatasetsComparisonType+"/("+*texDatasetsComparisonType+" "+DatasetsNames[0]+")"); 
TString* texRatioRun3Run2 = new TString("Run3 / Run2"); 

TString* texRho = new TString("#it{#rho} (GeV/#it{c})");
TString* texBkgFluctuationRandomCone = new TString("#delta #it{p}_{T} = #it{p}_{T,random cone} - #it{A}_{random cone} * #it{#rho} (GeV/c)");
TString* texEntriesNorm_BkgFluctuationYield = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#delta pt}");
TString* texCollNorm_BkgFluctuationYield = new TString("1/#it{N}_{coll} d#it{N}_{coll}/d#it{#delta pt}");
TString* texCollNorm_BkgFluctuationYield_CentWindow = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#delta pt}");
TString* texEntriesNormRho = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#rho}");
TString* texCollNorm_RhoYield = new TString("1/#it{N}_{coll} d#it{N}_{coll}/d#it{#rho}");
TString* texCollNorm_RhoYield_CentWindow = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{#rho}");
TString* texNoNorm_SelMultYield = new TString("d#it{N}_{coll}/d#it{N}_{tracks}^{selected}");
TString* texEntriesNorm_selMultYield = new TString("1/#it{N}_{tracks} d#it{N}_{coll}/d#it{N}_{tracks}^{selected}");
TString* texCollNorm_selMultYield = new TString("1/#it{N}_{coll} d#it{N}_{coll}/d#it{N}_{tracks}^{selected}");
TString* texCollNorm_selMultYield_CentWindow = new TString("1/#it{N}_{coll}^{in cent} d#it{N}_{coll}/d#it{N}_{tracks}^{selected}");


TString* texSigmaBkgFluctFit = new TString("#sigma_{#it{#delta pt}}^{fit}");
TString* texMeanBkgFluctFit = new TString("<#it{#delta pt}>^{fit}");

TString* texRhoMean = new TString("<#it{#rho}>");

TString* texTrackEfficiency = new TString("#epsilon = nTrack_{assoc} / nParticle_{gen}");
TString* texTrackPurity = new TString("p = nTrack_{sign} / (nTrack_{sign}+nTrack_{bkg})");


TString* texJetKinematicEfficiency = new TString("Jet kinematic efficiency");
TString* texJetEfficiency = new TString("#epsilon = nJet_{mcp,matched}^{sel coll} / nJet_{mcp}^{sel coll}");
TString* texJetPurity = new TString("#epsilon = nJet_{signal} / (nJet_{signal}+nJet_{background})");
TString* texFakeRatio = new TString("N^{mcd, matched}_{jets} / N^{mcd}_{jets}");

TString* texPartVxYield_EntriesNorm = new TString("1/#it{N}_{tracks} d#it{N}_{track}/dV_{x}");
TString* texPartVyYield_EntriesNorm = new TString("1/#it{N}_{tracks} d#it{N}_{track}/dV_{y}");
TString* texPartVzYield_EntriesNorm = new TString("1/#it{N}_{tracks} d#it{N}_{track}/dV_{z}");

TString* texPartYYield_EntriesNorm = new TString("1/#it{N}_{tracks} d#it{N}_{track}/dy");

TString* texPtPeak = new TString("#it{p}_{T,jet rec} - #it{A}_{jet} * #it{#rho} (GeV/#it{c}) of peak");
TString* texPtAtThreshold = new TString("#it{p}_{T,jet rec} - #it{A}_{jet} * #it{#rho} (GeV/#it{c}) at 95% threshold");

TString* texSystematicsPercent = new TString ("systematic error (%)");

//////////////////////////
////////// Jets //////////
TString* texNSubJettinessRatio = new TString("#it{#tau_{2}/#tau_{1}}");
TString* texDeltaR = new TString("#it{#Delta R}");
TString* texdN_dsubratio = new TString("1/N^{ev} d#it{N}/d(#it{#tau_{2}/#tau_{1}})");
TString* texdN_dDeltaR = new TString("1/N^{ev} d#it{N}/d(#it{#Delta R})");

TString* texJetPtYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetEtaYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#eta}");
TString* texJetPhiYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{#phi}");
TString* texJetPtYield_EntriesNorm = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetEtaYield_EntriesNorm = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{#eta}");
TString* texJetPhiYield_EntriesNorm = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{#phi}");
TString* texJetPtYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{jet}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texJetEtaYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{jet}/d#it{#eta}");
TString* texJetEtaYield_EntriesNorm_CentWindow = new TString("1/#it{N}_{jet}^{in cent} d#it{N}_{jet}/d#it{#eta}");
TString* texJetPhiYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{jet}/d#it{#phi}");
TString* texJetNormNTracksYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{N}_{tracks}");
TString* texEvtNormAreaYield = new TString("1/#it{N}_{ev} d#it{N}_{jet}/d#it{A}_{jet}");
TString* texJetNormAreaYield = new TString("1/#it{N}_{jet} d#it{N}_{jet}/d#it{A}_{jet}");


TString* texJetRatioEtaComparison = new TString("#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0}");
TString* texJetRatioTRDvsNoTRD = new TString("(#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0})_{onlyTRDtracks} / (#it{N}_{jet}^{ #eta > 0} / #it{N}_{jet}^{ #eta < 0})_{noTRDtracks}");

TString* texNGhost = new TString("#it{N}_{ghost}");
TString* texBinSpacing = new TString("Area steps");
TString* texJetRatioAreaJetVsNGhost = new TString("#it{A}_{jet}^{O2}/ #it{A}_{jet}^{0.005*nGhost}");

TString* texRatioUnfolded = new TString("ratio to unfolded");
TString* texRatioMcpUnfolded = new TString("mcp / unfolded");
TString* texRatioRun2Unfolded = new TString("run2 / run3");
TString* texRatioMeasuredUnfolded = new TString("measured / unfolded");
TString* texRatioRefoldedMeasured = new TString("refolded / measured");
TString* texRatioMcpFoldedVsMcp = new TString("mcp-folded / mcp");
TString* texRatioMcpFoldedUnfoldedMcp = new TString("mcp-folded unfolded / mcp");


TString* texTrackEfficiencyRatioEtaComparison = new TString("#epsilon^{ #eta > 0} / #epsilon^{ #eta < 0}");
TString* texTrackPurityRatioEtaComparison = new TString("purity^{ #eta > 0} / #purity^{ #eta < 0}");

TString* texResMC = new TString("(#it{p}_{T, part} - #it{p}_{T, det}) / #it{p}_{T, part}");


////////////////////////////
////////// Tracks //////////
TString* texTrackPtYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texTrackEtaYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{#eta}");
TString* texTrackPhiYield_EventNorm = new TString("1/#it{N}_{ev} d#it{N}_{track}/d#it{#phi}");
TString* texTrackPtYield_EntriesNorm = new TString("1/#it{N}_{track} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texTrackEtaYield_EntriesNorm = new TString("1/#it{N}_{track} d#it{N}_{track}/d#it{#eta}");
TString* texTrackPhiYield_EntriesNorm = new TString("1/#it{N}_{track} d#it{N}_{track}/d#it{#phi}");
TString* texTrackPtYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");
TString* texTrackEtaYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{track}/d#it{#eta}");
TString* texTrackPhiYield_EventNorm_CentWindow = new TString("1/#it{N}_{ev}^{in cent} d#it{N}_{track}/d#it{#phi}");
TString* texMeanPt = new TString("<#it{p}_{T}>");
TString* texMeanEta = new TString("<#eta>");
TString* texMeanNtracks = new TString("<N_{tracks}>");
TString* texTrackPtYield_JetsNorm = new TString("1/#it{N}_{jets} d#it{N}_{track}/d#it{p}_{T} (GeV/#it{c})^{-1}");


////////////////////////////
////////// Events //////////
TString* texWeightOccupancy = new TString("Weighted Occupancy");




#endif
