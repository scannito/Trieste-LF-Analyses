#include <utility>

#include "Riostream.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

#include "../AnalysisUtils/Parameters.h"
#include "../AnalysisUtils/FitFunctions.h"
#include "../AnalysisUtils/Plot.h"
#include "SignalExtraction.h"

using namespace std;

TH2F* Project2D(THnSparseF* hn, Int_t axistocut, Int_t binlow, Int_t binup, Int_t axistoproj1, Int_t axistoproj2, Option_t* option = "", std::string hname = "") 
{ 
    if (!hn) return 0;
    hn->GetAxis(axistocut)->SetRange(binlow, binup);
    TH2F* h2 = (TH2F*)hn->Projection(axistoproj1, axistoproj2, option);
    h2->SetName(hname.c_str());
    h2->SetDirectory(0);
    return h2;
}

TH2F* Project2D(THnSparseF* hn, Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistoproj1, Int_t axistoproj2, Option_t* option = "", string hname = "") 
{ 
    if (!hn) return 0;
    hn->GetAxis(axistocut1)->SetRange(binlow1, binup1);
    hn->GetAxis(axistocut2)->SetRange(binlow2, binup2);
    TH2F* h2 = (TH2F*)hn->Projection(axistoproj1, axistoproj2, option);
    h2->SetName(hname.c_str());
    h2->SetDirectory(0);
    return h2;
}
TH1F* Project1D(THnSparseF* hn, Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistocut3, Int_t binlow3, Int_t binup3, Int_t axistoproj, Option_t* option = "", string hname = "") 
{ 
    if (!hn) return 0;
    hn->GetAxis(axistocut1)->SetRange(binlow1, binup1);
    hn->GetAxis(axistocut2)->SetRange(binlow2, binup2);
    hn->GetAxis(axistocut3)->SetRange(binlow3, binup3);
    TH1F* h1 = (TH1F*)hn->Projection(axistoproj, option);
    h1->SetName(hname.c_str());
    h1->SetDirectory(0);
    return h1;
}

void fitDataAndMCClosure2D(int mode) 
{
    //**************************************************************************************************************
    // Histograms retrieving and projections depending on the mode: data (old or new) or MC closure test
    //**************************************************************************************************************

    //TFile* file1 = TFile::Open("AnalysisResultsPurMedium.root");
    //TFile* file1 = TFile::Open("AnalysisResultsPurLong.root");
    //TFile* file1 = TFile::Open("AnalysisResultsMediumNoTLV.root");
    //TDirectoryFile* lfphik0shortanalysis1 = (TDirectoryFile*)file1->Get("lf-phik0shortanalysis_id19979");

    //TFile* file1 = TFile::Open("AnalysisResultsPurMediumMod.root");
    //TDirectoryFile* lfphik0shortanalysis1 = (TDirectoryFile*)file1->Get("lf-phik0shortanalysis");

    TFile* file1;
    TDirectoryFile* phik0shortanalysis1;
    TDirectoryFile* PhipurHist;
    TDirectoryFile* eventHist;

    string purPhiHistName;
    array<string, nbin_deltay> purPhiK0SInvMassHistName;
    array<string, nbin_deltay> purPhiPiInvMassHistName;
    string eventHistName;
    string binEventHistName;

    TFile* file2;
    TDirectoryFile* phik0shortanalysis2;
    TDirectoryFile* PhiK0SHist;

    array<string, nbin_deltay> PhiK0SInvMassHistName;

    TDirectoryFile* phik0shortanalysis3;
    TDirectoryFile* PhiPiHist;

    array<string, nbin_deltay> PhiPiInvMassHistName;

    string outPath;
    string outFileName;

    if (mode == 0) {
        file1 = TFile::Open("AnalysisResultsMediumNoTLV.root");
        //file1 = TFile::Open("AnalysisResultsLongNoTLV.root");
        phik0shortanalysis1 = (TDirectoryFile*)file1->Get("lf-phik0shortanalysis_id19979");
        PhipurHist = (TDirectoryFile*)phik0shortanalysis1->Get("PhipurHist");
        eventHist = (TDirectoryFile*)phik0shortanalysis1->Get("eventHist");

        purPhiHistName = "h2PhipurInvMass";
        purPhiK0SInvMassHistName = {"h3PhipurK0SInvMassInc", "h3PhipurK0SInvMassFCut", "h3PhipurK0SInvMassSCut"};
        purPhiPiInvMassHistName = {"h3PhipurPiInvMassInc", "h3PhipurPiInvMassFCut", "h3PhipurPiInvMassSCut"};
        eventHistName = "hEventSelection";
        binEventHistName = "With at least a #phi cand";

        //file2 = TFile::Open("AnalysisResultsMediumNoTLV.root");
        phik0shortanalysis2 = (TDirectoryFile*)file1->Get("lf-phik0shortanalysis_id19977");
        PhiK0SHist = (TDirectoryFile*)phik0shortanalysis2->Get("PhiK0SHist");

        PhiK0SInvMassHistName = {"h4PhiK0SSEInc", "h4PhiK0SSEFCut", "h4PhiK0SSESCut"};

        phik0shortanalysis3 = (TDirectoryFile*)file1->Get("lf-phik0shortanalysis_id19978");
        PhiPiHist = (TDirectoryFile*)phik0shortanalysis3->Get("PhiPionHist");

        PhiPiInvMassHistName = {"h5PhiPiSEInc", "h5PhiPiSEFCut", "h5PhiPiSESCut"};

        outPath = "ResultsDataOld/";
        outFileName = "ResultsDataOld2D.root";

        /*file1 = TFile::Open("AnalysisResultsMedium.root");
        phik0shortanalysis1 = (TDirectoryFile*)file1->Get("lf-phik0shortanalysis_id19979");
        PhipurHist = (TDirectoryFile*)phik0shortanalysis1->Get("PhipurHist");
        eventHist = (TDirectoryFile*)phik0shortanalysis1->Get("eventHist");

        purPhiHistName = "h2PhipurInvMass";
        purPhiK0SInvMassHistName = {"h3PhipurK0SInvMassInclusive", "h3PhipurK0SInvMassFirstCut", "h3PhipurK0SInvMassSecondCut"};
        purPhiPiInvMassHistName = {"h3PhipurPiInvMassInclusive", "h3PhipurPiInvMassFirstCut", "h3PhipurPiInvMassSecondCut"};
        eventHistName = "hEventSelection";
        binEventHistName = "With at least a #phi cand";

        file2 = TFile::Open("AnalysisResultsMedium.root");
        phik0shortanalysis2 = (TDirectoryFile*)file2->Get("lf-phik0shortanalysis_id19977");
        PhiK0SHist = (TDirectoryFile*)phik0shortanalysis2->Get("PhiK0SHist");

        PhiK0SInvMassHistName = {"h4PhiK0SSEInc", "h4PhiK0SSEFCut", "h4PhiK0SSESCut"};

        phik0shortanalysis3 = (TDirectoryFile*)file2->Get("lf-phik0shortanalysis_id19978");
        PhiPiHist = (TDirectoryFile*)phik0shortanalysis3->Get("PhiPionHist");

        PhiPiInvMassHistName = {"h5PhiPiSEInc", "h5PhiPiSEFCut", "h5PhiPiSESCut"};

        outPath = "ResultsDataOld/";
        outFileName = "ResultsDataOld.root";*/
    } else if (mode == 1) {
        //file1 = TFile::Open("AnalysisResultsSmallNewSel.root");
        file1 = TFile::Open("AnalysisResultsMediumNewSel.root");
        phik0shortanalysis1 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id19979");
        PhipurHist = (TDirectoryFile*)phik0shortanalysis1->Get("dataPhiHist");
        eventHist = (TDirectoryFile*)phik0shortanalysis1->Get("dataEventHist");

        purPhiHistName = "h2PhipurInvMass";
        purPhiK0SInvMassHistName = {"h3PhipurK0SInvMassInc", "h3PhipurK0SInvMassFCut", "h3PhipurK0SInvMassSCut"};
        purPhiPiInvMassHistName = {"h3PhipurPiInvMassInc", "h3PhipurPiInvMassFCut", "h3PhipurPiInvMassSCut"};
        eventHistName = "hEventSelection";
        binEventHistName = "With at least a #phi cand";

        phik0shortanalysis2 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id19977");
        PhiK0SHist = (TDirectoryFile*)phik0shortanalysis2->Get("dataPhiK0SHist");

        PhiK0SInvMassHistName = {"h4PhiK0SSEInc", "h4PhiK0SSEFCut", "h4PhiK0SSESCut"};

        phik0shortanalysis3 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id19978");
        PhiPiHist = (TDirectoryFile*)phik0shortanalysis3->Get("dataPhiPionHist");

        PhiPiInvMassHistName = {"h5PhiPiSEInc", "h5PhiPiSEFCut", "h5PhiPiSESCut"};

        outPath = "ResultsDataNew/";
        outFileName = "ResultsDataNew2D.root";
    } else if (mode == 2) {
        //file1 = TFile::Open("AnalysisResultsMCClosure.root");
        //file1 = TFile::Open("AnalysisResultsMCClosureNewSel.root");
        file1 = TFile::Open("AnalysisResultsMCLHC24f3b.root");
        phik0shortanalysis1 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20119");
        PhipurHist = (TDirectoryFile*)phik0shortanalysis1->Get("closureMCPhiHist");
        eventHist = (TDirectoryFile*)phik0shortanalysis1->Get("mcEventHist");

        purPhiHistName = "h2MCPhipurInvMass";
        purPhiK0SInvMassHistName = {"h3MCPhipurK0SInvMassInc", "h3MCPhipurK0SInvMassFCut", "h3MCPhipurK0SInvMassSCut"};
        purPhiPiInvMassHistName = {"h3MCPhipurPiInvMassInc", "h3MCPhipurPiInvMassFCut", "h3MCPhipurPiInvMassSCut"};
        eventHistName = "hRecMCEventSelection";
        binEventHistName = "With at least a #phi";

        //file2 = TFile::Open("AnalysisResultsMCClosure.root");
        phik0shortanalysis2 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20120");
        PhiK0SHist = (TDirectoryFile*)phik0shortanalysis2->Get("closureMCPhiK0SHist");

        PhiK0SInvMassHistName = {"h4ClosureMCPhiK0SSEInc", "h4ClosureMCPhiK0SSEFCut", "h4ClosureMCPhiK0SSESCut"};

        phik0shortanalysis3 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20121");
        PhiPiHist = (TDirectoryFile*)phik0shortanalysis3->Get("closureMCPhiPionHist");

        PhiPiInvMassHistName = {"h5ClosureMCPhiPiSEInc", "h5ClosureMCPhiPiSEFCut", "h5ClosureMCPhiPiSESCut"};

        outPath = "ResultsMCClosure/";
        outFileName = "ResultsMCClosure2D.root";
    } else if (mode == 3) {
        file1 = TFile::Open("AnalysisResultsNewMCClosure.root");
        phik0shortanalysis1 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20119");
        PhipurHist = (TDirectoryFile*)phik0shortanalysis1->Get("closureMCPhiHist");
        eventHist = (TDirectoryFile*)phik0shortanalysis1->Get("mcEventHist");

        purPhiHistName = "h2MCPhipurInvMass";
        purPhiK0SInvMassHistName = {"h3MCPhipurK0SInvMassInc", "h3MCPhipurK0SInvMassFCut", "h3MCPhipurK0SInvMassSCut"};
        purPhiPiInvMassHistName = {"h3MCPhipurPiInvMassInc", "h3MCPhipurPiInvMassFCut", "h3MCPhipurPiInvMassSCut"};
        eventHistName = "hRecMCEventSelection";
        binEventHistName = "With at least a #phi";

        //file2 = TFile::Open("AnalysisResultsMCClosure.root");
        phik0shortanalysis2 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id24153");
        PhiK0SHist = (TDirectoryFile*)phik0shortanalysis2->Get("closureMCPhiK0SHist");

        PhiK0SInvMassHistName = {"h4ClosureMCPhiK0SSEInc", "h4ClosureMCPhiK0SSEFCut", "h4ClosureMCPhiK0SSESCut"};

        phik0shortanalysis3 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id24154");
        PhiPiHist = (TDirectoryFile*)phik0shortanalysis3->Get("closureMCPhiPionHist");

        PhiPiInvMassHistName = {"h5ClosureMCPhiPiSEInc", "h5ClosureMCPhiPiSEFCut", "h5ClosureMCPhiPiSESCut"};

        outPath = "ResultsNewMCClosure/";
        outFileName = "ResultsNewMCClosure2D.root";
    } else if (mode == 4) {
        file1 = TFile::Open("AnalysisResultsNewMCClosure.root");
        phik0shortanalysis1 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20119");
        PhipurHist = (TDirectoryFile*)phik0shortanalysis1->Get("closureMCPhiHist");
        eventHist = (TDirectoryFile*)phik0shortanalysis1->Get("mcEventHist");

        purPhiHistName = "h2MCPhipurInvMass";
        purPhiK0SInvMassHistName = {"h3MCPhipurK0SInvMassInc", "h3MCPhipurK0SInvMassFCut", "h3MCPhipurK0SInvMassSCut"};
        purPhiPiInvMassHistName = {"h3MCPhipurPiInvMassInc", "h3MCPhipurPiInvMassFCut", "h3MCPhipurPiInvMassSCut"};
        eventHistName = "hRecMCEventSelection";
        binEventHistName = "With at least a #phi";

        //file2 = TFile::Open("AnalysisResultsMCClosure.root");
        phik0shortanalysis2 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20120");
        PhiK0SHist = (TDirectoryFile*)phik0shortanalysis2->Get("closureMCPhiK0SHist");

        PhiK0SInvMassHistName = {"h4ClosureMCPhiK0SSEInc", "h4ClosureMCPhiK0SSEFCut", "h4ClosureMCPhiK0SSESCut"};

        phik0shortanalysis3 = (TDirectoryFile*)file1->Get("phik0shortanalysis_id20121");
        PhiPiHist = (TDirectoryFile*)phik0shortanalysis3->Get("closureMCPhiPionHist");

        PhiPiInvMassHistName = {"h5ClosureMCPhiPiSEInc", "h5ClosureMCPhiPiSEFCut", "h5ClosureMCPhiPiSESCut"};

        outPath = "ResultsNewMCClosure/";
        outFileName = "ResultsNewMCClosure22D.root";
    } else return;

    //**************************************************************************************************************

    TH1F* hEventSelection = (TH1F*)eventHist->Get(eventHistName.c_str());
    hEventSelection->SetDirectory(0);
    Int_t binNumber = hEventSelection->GetXaxis()->FindBin(binEventHistName.c_str());
    Double_t nEventsPhi = hEventSelection->GetBinContent(binNumber);

    //**************************************************************************************************************

    TH2F* h2PhipurInvMass = (TH2F*)PhipurHist->Get(purPhiHistName.c_str());
    h2PhipurInvMass->SetDirectory(0);

    TH1F* h1PhipurInvMass[nbin_mult];
    for (int j = 0; j < nbin_mult; j++) {
        h1PhipurInvMass[j] = (TH1F*)h2PhipurInvMass->ProjectionY(Form("Phipur%i",j), j+1, j+1);
    }

    TH1F* h1PhipurInvMassMB = (TH1F*)h2PhipurInvMass->ProjectionY("PhipurMB", 1, nbin_mult);

    //**************************************************************************************************************

    TH3F* h3PhipurK0SInvMass[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        h3PhipurK0SInvMass[i] = (TH3F*)PhipurHist->Get(purPhiK0SInvMassHistName[i].c_str());
        h3PhipurK0SInvMass[i]->SetDirectory(0);
    }

    TH1F* h1PhipurK0SInvMass[nbin_deltay][nbin_mult];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhipurK0SInvMass[i][j] = (TH1F*)h3PhipurK0SInvMass[i]->ProjectionZ(Form("PhipurK0S%i%i", i, j), j+1, j+1, 1, nbin_pTK0S);
        }
    }

    TH1F* h1PhipurK0SInvMassMB[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        h1PhipurK0SInvMassMB[i] = (TH1F*)h3PhipurK0SInvMass[i]->ProjectionZ(Form("PhipurK0SMB%i", i), 1, nbin_mult, 1, nbin_pTK0S);
    }

    //**************************************************************************************************************

    TH3F* h3PhipurPiInvMass[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        h3PhipurPiInvMass[i] = (TH3F*)PhipurHist->Get(purPhiPiInvMassHistName[i].c_str());
        h3PhipurPiInvMass[i]->SetDirectory(0);
    }

    TH1F* h1PhipurPiInvMass[nbin_deltay][nbin_mult];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhipurPiInvMass[i][j] = (TH1F*)h3PhipurPiInvMass[i]->ProjectionZ(Form("PhipurPi%i%i", i, j), j+1, j+1, 1, nbin_pTPi);
        }
    }

    TH1F* h1PhipurPiInvMassMB[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        h1PhipurPiInvMassMB[i] = (TH1F*)h3PhipurPiInvMass[i]->ProjectionZ(Form("PhipurPiMB%i", i), 1, nbin_mult, 1, nbin_pTPi);
    }

    //**************************************************************************************************************

    THnSparseF* h4PhiK0SInvMass[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        h4PhiK0SInvMass[i] = (THnSparseF*)PhiK0SHist->Get(PhiK0SInvMassHistName[i].c_str());
    }

    TH2F* h2PhiK0SInvMass[nbin_deltay][nbin_mult][nbin_pTK0S];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pTK0S; k++) {
                h2PhiK0SInvMass[i][j][k] = Project2D(h4PhiK0SInvMass[i], 0, j+1, j+1, 1, k+1, k+1, 3, 2, "", Form("h2PhiK0SInvMass%i_%i_%i", i, j, k));
            }
        }
    }

    TH2F* h2PhiK0SInvMassMB[nbin_deltay][nbin_pTK0S];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int k = 0; k < nbin_pTK0S; k++) {
            h2PhiK0SInvMassMB[i][k] = Project2D(h4PhiK0SInvMass[i], 0, 1, nbin_mult, 1, k+1, k+1, 3, 2, "", Form("h2PhiK0SInvMassMB%i_%i", i, k));
        }
    }   

    //**************************************************************************************************************

    THnSparseF* h5PhiInvMassPiNSigma[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        h5PhiInvMassPiNSigma[i] = (THnSparseF*)PhiPiHist->Get(PhiPiInvMassHistName[i].c_str());
    }

    TH2F* h2PhiInvMassPiNSigmaTPC[nbin_deltay][nbin_mult][nbin_pTPi];
    TH2F* h2PhiInvMassPiNSigmaTOF[nbin_deltay][nbin_mult][nbin_pTPi];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pTPi; k++) {
                h2PhiInvMassPiNSigmaTPC[i][j][k] = Project2D(h5PhiInvMassPiNSigma[i], 0, j+1, j+1, 1, k+1, k+1, 4, 2, "", Form("h2PhiInvMassPiNSigmaTPC%i_%i_%i", i, j, k));
                h2PhiInvMassPiNSigmaTOF[i][j][k] = Project2D(h5PhiInvMassPiNSigma[i], 0, j+1, j+1, 1, k+1, k+1, 4, 3, "", Form("h2PhiInvMassPiNSigmaTOF%i_%i_%i", i, j, k));
            }
        }
    }

    TH2F* h2PhiInvMassPiNSigmaTPCMB[nbin_deltay][nbin_pTPi];
    TH2F* h2PhiInvMassPiNSigmaTOFMB[nbin_deltay][nbin_pTPi];
    for (int i = 0; i < nbin_deltay; i++) {
        for (int k = 0; k < nbin_pTPi; k++) {
            h2PhiInvMassPiNSigmaTPCMB[i][k] = Project2D(h5PhiInvMassPiNSigma[i], 0, 1, nbin_mult, 1, k+1, k+1, 4, 2, "", Form("h2PhiInvMassPiNSigmaTPCMB%i_%i", i, k));
            h2PhiInvMassPiNSigmaTOFMB[i][k] = Project2D(h5PhiInvMassPiNSigma[i], 0, 1, nbin_mult, 1, k+1, k+1, 4, 3, "", Form("h2PhiInvMassPiNSigmaTOFMB%i_%i", i, k));
        }
    }

    //**************************************************************************************************************

    PhipurHist->Close();
    eventHist->Close();
    PhiK0SHist->Close();
    PhiPiHist->Close();
    phik0shortanalysis1->Close();
    phik0shortanalysis2->Close();
    phik0shortanalysis3->Close();
    
    //****************************************************************************************************
    // Multiplicity-dependent purities
    //****************************************************************************************************

    Double_t purityVoigt[nbin_mult] = {0}, errpurityVoigt[nbin_mult] = {0};
    Double_t purityVoigtK0S[nbin_deltay][nbin_mult] = {0}, errpurityVoigtK0S[nbin_deltay][nbin_mult] = {0};
    Double_t purityVoigtPi[nbin_deltay][nbin_mult] = {0}, errpurityVoigtPi[nbin_deltay][nbin_mult] = {0};

    for (int j = 0; j < nbin_mult; j++) {
        tie(purityVoigt[j], errpurityVoigt[j]) = GetPhiPurityAndError(h1PhipurInvMass[j], Form("Phipur%i", j), mode-1, 0, {j});
        for (int i = 0; i < nbin_deltay; i++) {
            tie(purityVoigtK0S[i][j], errpurityVoigtK0S[i][j]) = GetPhiPurityAndError(h1PhipurK0SInvMass[i][j], Form("PhipurK0S%i%i", i, j), mode-1, 1, {i, j});
            tie(purityVoigtPi[i][j], errpurityVoigtPi[i][j]) = GetPhiPurityAndError(h1PhipurPiInvMass[i][j], Form("PhipurPi%i%i", i, j), mode-1, 2, {i, j});
        }
    }

    //********************************************************************************************

    Double_t purityPhi[nbin_mult] = {0}, errpurityPhi[nbin_mult] = {0};
    Double_t purityPhiK0S[nbin_deltay][nbin_mult] = {0}, errpurityPhiK0S[nbin_deltay][nbin_mult] = {0};
    Double_t purityPhiPi[nbin_deltay][nbin_mult] = {0}, errpurityPhiPi[nbin_deltay][nbin_mult] = {0};

    for (int j = 0; j < nbin_mult; j++) {
        purityPhi[j] = purityVoigt[j];
        errpurityPhi[j] = errpurityVoigt[j];
        
        for (int i = 0; i < nbin_deltay; i++) {
            purityPhiK0S[i][j] = purityVoigtK0S[i][j];
            errpurityPhiK0S[i][j] = errpurityVoigtK0S[i][j];

            purityPhiPi[i][j] = purityVoigtPi[i][j];
            errpurityPhiPi[i][j] = errpurityVoigtPi[i][j];
        }
    }

    TMultiGraph* mgPhipurK0S = new TMultiGraph();
    TMultiGraph* mgPhipurPi = new TMultiGraph();

    TGraphAsymmErrors* PurityPhi = new TGraphAsymmErrors(nbin_mult, mult.data(), purityPhi, errmult.data(), errmult.data(), errpurityPhi, errpurityPhi);
    PlotFeatures(PurityPhi, Markers[0], kViolet, 1, 1, kViolet, 2, 3001, kViolet, 0.4, mgPhipurK0S);
    PlotFeatures(PurityPhi, Markers[0], kViolet, 1, 1, kViolet, 2, 3001, kViolet, 0.4, mgPhipurPi);

    TGraphAsymmErrors* PurityPhiK0S[nbin_deltay];
    TGraphAsymmErrors* PurityPhiPi[nbin_deltay];
    for (int i = 0; i < nbin_deltay; i++) {
        PurityPhiK0S[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), purityPhiK0S[i], errmult.data(), errmult.data(), errpurityPhiK0S[i], errpurityPhiK0S[i]);
        PurityPhiPi[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), purityPhiPi[i], errmult.data(), errmult.data(), errpurityPhiPi[i], errpurityPhiPi[i]);
        if (i == 0) {
            PlotFeatures(PurityPhiK0S[i], Markers[0], kBlack, 1, 1, kBlack, 2, 3001, kBlack, 0.4, mgPhipurK0S);
            PlotFeatures(PurityPhiPi[i], Markers[0], kBlack, 1, 1, kBlack, 2, 3001, kBlack, 0.4, mgPhipurPi);
        } else if (i == 1) {
            PlotFeatures(PurityPhiK0S[i], Markers[0], kGreen+3, 1, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgPhipurK0S);
            PlotFeatures(PurityPhiPi[i], Markers[0], kGreen+3, 1, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgPhipurPi);
        } else if (i == 2) {
            PlotFeatures(PurityPhiK0S[i], Markers[0], kRed+1, 1, 1, kRed+1, 2, 3001, kRed+1, 0.4, mgPhipurK0S);
            PlotFeatures(PurityPhiPi[i], Markers[0], kRed+1, 1, 1, kRed+1, 2, 3001, kRed+1, 0.4, mgPhipurPi);
        }
    }

    //********************************************************************************************

    TCanvas* cPurityPhiK0S = new TCanvas("cPurityPhiK0S", "cPurityPhiK0S", 800, 800);
    cPurityPhiK0S->cd();
    gPad->SetMargin(0.16,0.01,0.13,0.06);
    gStyle->SetOptStat(0);
    mgPhipurK0S->Draw("AP");
    mgPhipurK0S->SetTitle("; #LTd#it{N}_{ch}/d#eta#GT_{|#eta|<0.5} ; S/(S+B)");
    mgPhipurK0S->GetYaxis()->SetTitleOffset(1.4);
    mgPhipurK0S->GetYaxis()->SetTitleSize(0.045);
    mgPhipurK0S->GetYaxis()->SetLabelSize(0.045);
    mgPhipurK0S->GetXaxis()->SetTitleOffset(1.2);
    mgPhipurK0S->GetXaxis()->SetTitleSize(0.045);
    mgPhipurK0S->GetXaxis()->SetLabelSize(0.045);
    //PurityPhiK0S[2]->Draw("AP");

    TLegend* legPurityPhiK0S1 = new TLegend(0.45, 0.85, 0.8, 0.88);
    legPurityPhiK0S1->SetHeader("#bf{This work}");
    legPurityPhiK0S1->SetTextSize(0.05);
    legPurityPhiK0S1->SetLineWidth(0);
    legPurityPhiK0S1->Draw("same");

    TLegend* legPurityPhiK0S2 = new TLegend(0.45, 0.65, 0.8, 0.85);
    legPurityPhiK0S2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legPurityPhiK0S2->AddEntry(PurityPhi, "MB", "p");
    legPurityPhiK0S2->AddEntry(PurityPhiK0S[0], "|#Delta#it{y}| < 1.0 (inclusive)", "p");
    legPurityPhiK0S2->AddEntry(PurityPhiK0S[1], "|#Delta#it{y}| < 0.5", "p");
    legPurityPhiK0S2->AddEntry(PurityPhiK0S[2], "|#Delta#it{y}| < 0.1", "p");
    legPurityPhiK0S2->SetTextSize(0.045);
    legPurityPhiK0S2->SetLineWidth(0);
    legPurityPhiK0S2->Draw("same");

    string outNamePurityPhiK0S = outPath + "purPhiK0S.root";
    cPurityPhiK0S->SaveAs(outNamePurityPhiK0S.c_str());
    outNamePurityPhiK0S = outPath + "purPhiK0S.pdf";
    cPurityPhiK0S->SaveAs(outNamePurityPhiK0S.c_str());

    //********************************************************************************************

    TCanvas* cPurityPhiPi = new TCanvas("cPurityPhiPi", "cPurityPhiPi", 800, 800);
    cPurityPhiPi->cd();
    gPad->SetMargin(0.16,0.01,0.13,0.06);
    gStyle->SetOptStat(0);
    mgPhipurPi->SetTitle("; #LTd#it{N}_{ch}/d#eta#GT_{|#eta|<0.5} ; S/(S+B)");
    mgPhipurPi->Draw("AP");
    mgPhipurPi->GetYaxis()->SetTitleOffset(1.4);
    mgPhipurPi->GetYaxis()->SetTitleSize(0.045);
    mgPhipurPi->GetYaxis()->SetLabelSize(0.045);
    mgPhipurPi->GetXaxis()->SetTitleOffset(1.2);
    mgPhipurPi->GetXaxis()->SetTitleSize(0.045);
    mgPhipurPi->GetXaxis()->SetLabelSize(0.045);

    TLegend* legPurityPhiPi1 = new TLegend(0.45, 0.85, 0.8, 0.88);
    legPurityPhiPi1->SetHeader("#bf{This work}");
    legPurityPhiPi1->SetTextSize(0.05);
    legPurityPhiPi1->SetLineWidth(0);
    legPurityPhiPi1->Draw("same");

    TLegend* legPurityPhiPi2 = new TLegend(0.45, 0.65, 0.8, 0.85);
    legPurityPhiPi2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legPurityPhiPi2->AddEntry(PurityPhi, "MB", "p");
    legPurityPhiPi2->AddEntry(PurityPhiPi[0], "|#Delta#it{y}| < 1.0 (inclusive)", "p");
    legPurityPhiPi2->AddEntry(PurityPhiPi[1], "|#Delta#it{y}| < 0.5", "p");
    legPurityPhiPi2->AddEntry(PurityPhiPi[2], "|#Delta#it{y}| < 0.1", "p");
    legPurityPhiPi2->SetTextSize(0.045);
    legPurityPhiPi2->SetLineWidth(0);
    legPurityPhiPi2->Draw("same");

    string outNamePurityPhiPi = outPath + "purPhiPi.root";
    cPurityPhiPi->SaveAs(outNamePurityPhiPi.c_str());
    outNamePurityPhiPi = outPath + "purPhiPi.pdf";
    cPurityPhiPi->SaveAs(outNamePurityPhiPi.c_str());
    
    //********************************************************************************************
    // Multiplicity-integrated purities
    //********************************************************************************************

    Double_t purityVoigtMB{0}, errpurityVoigtMB{0};
    Double_t purityVoigtK0SMB[nbin_deltay] = {0}, errpurityVoigtK0SMB[nbin_deltay] = {0};
    Double_t purityVoigtPiMB[nbin_deltay] = {0}, errpurityVoigtPiMB[nbin_deltay] = {0};

    tie(purityVoigtMB, errpurityVoigtMB) = GetPhiPurityAndError(h1PhipurInvMassMB, "PhipurMB", mode-1, 0, {});
    for (int i = 0; i < nbin_deltay; i++) {
        tie(purityVoigtK0SMB[i], errpurityVoigtK0SMB[i]) = GetPhiPurityAndError(h1PhipurK0SInvMassMB[i], Form("PhipurK0SMB%i", i), mode-1, 1, {i});
        tie(purityVoigtPiMB[i], errpurityVoigtPiMB[i]) = GetPhiPurityAndError(h1PhipurPiInvMassMB[i], Form("PhipurPiMB%i", i), mode-1, 2, {i});
    }

    //********************************************************************************************

    Double_t purityPhiMB{0}, errpurityPhiMB{0};
    Double_t purityPhiK0SMB[nbin_deltay] = {0}, errpurityPhiK0SMB[nbin_deltay] = {0};
    Double_t purityPhiPiMB[nbin_deltay] = {0}, errpurityPhiPiMB[nbin_deltay] = {0};

    purityPhiMB = purityVoigtMB;
    errpurityPhiMB = errpurityVoigtMB;
    
    for (int i = 0; i < nbin_deltay; i++) {
        purityPhiK0SMB[i] = purityVoigtK0SMB[i];
        errpurityPhiK0SMB[i] = errpurityVoigtK0SMB[i];

        purityPhiPiMB[i] = purityVoigtPiMB[i];
        errpurityPhiPiMB[i] = errpurityVoigtPiMB[i];
    }

    //return;

    //********************************************************************************************
    // Signal extraction
    //********************************************************************************************

    Double_t PhiK0SYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTK0S] = {0}, errPhiK0SYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTK0S] = {0};
    TH1D* h1PhiK0SYield[nbin_deltay][nbin_mult];

    Double_t PhiK0SYieldpTdiffMB[nbin_deltay][nbin_pTK0S] = {0}, errPhiK0SYieldpTdiffMB[nbin_deltay][nbin_pTK0S] = {0};
    TH1D* h1PhiK0SYieldMB[nbin_deltay];
    
    Double_t PhiPiTPCYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTPi] = {0}, errPhiPiTPCYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTPi] = {0};
    TH1D* h1PhiPiTPCYield[nbin_deltay][nbin_mult];

    Double_t PhiPiTOFYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTPi] = {0}, errPhiPiTOFYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTPi] = {0};
    TH1D* h1PhiPiTOFYield[nbin_deltay][nbin_mult];

    Double_t PhiPiYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTPi] = {0}, errPhiPiYieldpTdiff[nbin_deltay][nbin_mult][nbin_pTPi] = {0};
    TH1D* h1PhiPiYield[nbin_deltay][nbin_mult];

    Double_t PhiPiTPCYieldpTdiffMB[nbin_deltay][nbin_pTPi] = {0}, errPhiPiTPCYieldpTdiffMB[nbin_deltay][nbin_pTPi] = {0};
    TH1D* h1PhiPiTPCYieldMB[nbin_deltay];

    Double_t PhiPiTOFYieldpTdiffMB[nbin_deltay][nbin_pTPi] = {0}, errPhiPiTOFYieldpTdiffMB[nbin_deltay][nbin_pTPi] = {0};
    TH1D* h1PhiPiTOFYieldMB[nbin_deltay];

    Double_t PhiPiYieldpTdiffMB[nbin_deltay][nbin_pTPi] = {0}, errPhiPiYieldpTdiffMB[nbin_deltay][nbin_pTPi] = {0};
    TH1D* h1PhiPiYieldMB[nbin_deltay];

    //********************************************************************************************

    string outNameCanvasK0S = outPath + "fitCanvasK0S2D.root";
    TFile* fileCanvasK0S = new TFile(outNameCanvasK0S.c_str(), "RECREATE");

    string outNameCanvasPiTPC = outPath + "fitCanvasPiTPC2D.root";
    TFile* fileCanvasPiTPC = new TFile(outNameCanvasPiTPC.c_str(), "RECREATE");

    string outNameCanvasPiTOF = outPath + "fitCanvasPiTOF2D.root";
    TFile* fileCanvasPiTOF = new TFile(outNameCanvasPiTOF.c_str(), "RECREATE");

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhiK0SYield[i][j] = new TH1D(Form("h1PhiK0SYield%i_%i", i, j), Form("h1PhiK0SYield%i_%i", i, j), nbin_pTK0S, pTK0S_axis.data());
            h1PhiK0SYield[i][j]->SetTitle("; ; 1/N_{ev,#phi} d^{2}N_{K^{0}_{S}}/d#it{y}d#it{p}_{T} [(GeV/#it{c})^{-1}]");

            h1PhiPiTPCYield[i][j] = new TH1D(Form("h1PhiPiTPCYield%i_%i", i, j), Form("h1PhiPiTPCYield%i_%i", i, j), nbin_pTPi, pTPi_axis.data());
            h1PhiPiTPCYield[i][j]->SetTitle("; ; 1/N_{ev,#phi} d^{2}N_{(#pi^{+}+#pi^{#minus})}/d#it{y}d#it{p}_{T} [(GeV/#it{c})^{-1}]");

            h1PhiPiTOFYield[i][j] = new TH1D(Form("h1PhiPiTOFYield%i_%i", i, j), Form("h1PhiPiTOFYield%i_%i", i, j), nbin_pTPi, pTPi_axis.data());
            h1PhiPiTOFYield[i][j]->SetTitle("; ; 1/N_{ev,#phi} d^{2}N_{(#pi^{+}+#pi^{#minus})}/d#it{y}d#it{p}_{T} [(GeV/#it{c})^{-1}]");

            h1PhiPiYield[i][j] = new TH1D(Form("h1PhiPiYield%i_%i", i, j), Form("h1PhiPiYield%i_%i", i, j), nbin_pTPi, pTPi_axis.data());
            h1PhiPiYield[i][j]->SetTitle("; ; 1/N_{ev,#phi} d^{2}N_{(#pi^{+}+#pi^{#minus})}/d#it{y}d#it{p}_{T} [(GeV/#it{c})^{-1}]");

            for (int k = 0; k < nbin_pTK0S; k++) {
                //if (i != 0 || j != 0 || k != 0) continue;
                //if (i != 2 || k != 3) continue;
                //if (i != 0) continue;

                tie(PhiK0SYieldpTdiff[i][j][k], errPhiK0SYieldpTdiff[i][j][k]) = FitPhiK0S(h2PhiK0SInvMass[i][j][k], {i, j, k}, fileCanvasK0S);
                PhiK0SYieldpTdiff[i][j][k] = PhiK0SYieldpTdiff[i][j][k] / deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pTK0S_axis[k+1] - pTK0S_axis[k]) /*/ purityPhiK0S[i][j]*/ / (nEventsPhi * purityPhi[j]);
                errPhiK0SYieldpTdiff[i][j][k] = errPhiK0SYieldpTdiff[i][j][k] / deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pTK0S_axis[k+1] - pTK0S_axis[k]) /*/ purityPhiK0S[i][j]*/ / (nEventsPhi * purityPhi[j]);
                if (mode != 3) {
                    PhiK0SYieldpTdiff[i][j][k] = PhiK0SYieldpTdiff[i][j][k] / purityPhiK0S[i][j];
                    errPhiK0SYieldpTdiff[i][j][k] = errPhiK0SYieldpTdiff[i][j][k] / purityPhiK0S[i][j];
                }

                h1PhiK0SYield[i][j]->SetBinContent(k+1, PhiK0SYieldpTdiff[i][j][k]);
                h1PhiK0SYield[i][j]->SetBinError(k+1, errPhiK0SYieldpTdiff[i][j][k]);
            }

            for (int k = 0; k < nbin_pTPi; k++) {
                //if (i != 2 || j != 6 || k != 6) continue;
                //if (i != 0 || j != 0) continue;
                //if (i != 0) continue;

                tie(PhiPiTPCYieldpTdiff[i][j][k], errPhiPiTPCYieldpTdiff[i][j][k]) = FitPhiPi(h2PhiInvMassPiNSigmaTPC[i][j][k], {i, j, k}, 0, mode-2, fileCanvasPiTPC);
                PhiPiTPCYieldpTdiff[i][j][k] = PhiPiTPCYieldpTdiff[i][j][k] / deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPi[i][j]*/ / (nEventsPhi * purityPhi[j]);
                errPhiPiTPCYieldpTdiff[i][j][k] = errPhiPiTPCYieldpTdiff[i][j][k] / deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPi[i][j]*/ / (nEventsPhi * purityPhi[j]);
                //if (errPhiPiTPCYieldpTdiff[i][j][k] != errPhiPiTPCYieldpTdiff[i][j][k]) errPhiPiTPCYieldpTdiff[i][j][k] = 0.0;
                if (mode != 3) {
                    PhiPiTPCYieldpTdiff[i][j][k] = PhiPiTPCYieldpTdiff[i][j][k] / purityPhiPi[i][j];
                    errPhiPiTPCYieldpTdiff[i][j][k] = errPhiPiTPCYieldpTdiff[i][j][k] / purityPhiPi[i][j];
                }

                h1PhiPiTPCYield[i][j]->SetBinContent(k+1, PhiPiTPCYieldpTdiff[i][j][k]);
                h1PhiPiTPCYield[i][j]->SetBinError(k+1, errPhiPiTPCYieldpTdiff[i][j][k]);
                
                tie(PhiPiTOFYieldpTdiff[i][j][k], errPhiPiTOFYieldpTdiff[i][j][k]) = FitPhiPi(h2PhiInvMassPiNSigmaTOF[i][j][k], {i, j, k}, 1, mode-2, fileCanvasPiTOF);
                PhiPiTOFYieldpTdiff[i][j][k] = PhiPiTOFYieldpTdiff[i][j][k] / deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPi[i][j]*/ / (nEventsPhi * purityPhi[j]);
                errPhiPiTOFYieldpTdiff[i][j][k] = errPhiPiTOFYieldpTdiff[i][j][k] / deltay_axis[i] / ((mult_axis[j+1] - mult_axis[j]) / 100.0) / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPi[i][j]*/ / (nEventsPhi * purityPhi[j]);
                //if (errPhiPiTOFYieldpTdiff[i][j][k] != errPhiPiTOFYieldpTdiff[i][j][k]) errPhiPiTOFYieldpTdiff[i][j][k] = 0.0;
                if (mode != 3) {
                    PhiPiTOFYieldpTdiff[i][j][k] = PhiPiTOFYieldpTdiff[i][j][k] / purityPhiPi[i][j];
                    errPhiPiTOFYieldpTdiff[i][j][k] = errPhiPiTOFYieldpTdiff[i][j][k] / purityPhiPi[i][j];
                }

                h1PhiPiTOFYield[i][j]->SetBinContent(k+1, PhiPiTOFYieldpTdiff[i][j][k]);
                h1PhiPiTOFYield[i][j]->SetBinError(k+1, errPhiPiTOFYieldpTdiff[i][j][k]);

                PhiPiYieldpTdiff[i][j][k] = (PhiPiTPCYieldpTdiff[i][j][k] + PhiPiTOFYieldpTdiff[i][j][k]) / 2.0;
                errPhiPiYieldpTdiff[i][j][k] = TMath::Sqrt(TMath::Power(errPhiPiTPCYieldpTdiff[i][j][k], 2) + TMath::Power(errPhiPiTOFYieldpTdiff[i][j][k], 2)) / 2.0;

                h1PhiPiYield[i][j]->SetBinContent(k+1, PhiPiYieldpTdiff[i][j][k]);
                h1PhiPiYield[i][j]->SetBinError(k+1, errPhiPiYieldpTdiff[i][j][k]);
            }
        }
    }

    //********************************************************************************************

    for (int i = 0; i < nbin_deltay; i++) {
        h1PhiK0SYieldMB[i] = new TH1D(Form("h1PhiK0SYieldMB%i", i), Form("h1PhiK0SYieldMB%i", i), nbin_pTK0S, pTK0S_axis.data());

        h1PhiPiTPCYieldMB[i] = new TH1D(Form("h1PhiPiTPCYieldMB%i", i), Form("h1PhiPiTPCYieldMB%i", i), nbin_pTPi, pTPi_axis.data());
        h1PhiPiTOFYieldMB[i] = new TH1D(Form("h1PhiPiTOFYieldMB%i", i), Form("h1PhiPiTOFYieldMB%i", i), nbin_pTPi, pTPi_axis.data());
        h1PhiPiYieldMB[i] = new TH1D(Form("h1PhiPiYieldMB%i", i), Form("h1PhiPiYieldMB%i", i), nbin_pTPi, pTPi_axis.data());

        for (int k = 0; k < nbin_pTK0S; k++) {

            //if (i != 0 || k != 1) continue;
            //if (i != 2) continue;
            //continue;

            tie(PhiK0SYieldpTdiffMB[i][k], errPhiK0SYieldpTdiffMB[i][k]) = FitPhiK0S(h2PhiK0SInvMassMB[i][k], {i, k}, fileCanvasK0S);
            PhiK0SYieldpTdiffMB[i][k] = PhiK0SYieldpTdiffMB[i][k] / deltay_axis[i] / (pTK0S_axis[k+1] - pTK0S_axis[k]) /*/ purityPhiK0SMB[i]*/ / (nEventsPhi * purityPhiMB);
            errPhiK0SYieldpTdiffMB[i][k] = errPhiK0SYieldpTdiffMB[i][k] / deltay_axis[i] / (pTK0S_axis[k+1] - pTK0S_axis[k]) /*/ purityPhiK0SMB[i]*/ / (nEventsPhi * purityPhiMB);
            if (mode != 3) {
                PhiK0SYieldpTdiffMB[i][k] = PhiK0SYieldpTdiffMB[i][k] / purityPhiK0SMB[i];
                errPhiK0SYieldpTdiffMB[i][k] = errPhiK0SYieldpTdiffMB[i][k] / purityPhiK0SMB[i];
            }

            h1PhiK0SYieldMB[i]->SetBinContent(k+1, PhiK0SYieldpTdiffMB[i][k]);
            h1PhiK0SYieldMB[i]->SetBinError(k+1, errPhiK0SYieldpTdiffMB[i][k]);
        }

        for (int k = 0; k < nbin_pTPi; k++) {
            //if (i != 0 || k != 0) continue;
            //if (i != 0) continue;
            //continue;

            tie(PhiPiTPCYieldpTdiffMB[i][k], errPhiPiTPCYieldpTdiffMB[i][k]) = FitPhiPi(h2PhiInvMassPiNSigmaTPCMB[i][k], {i, k}, 0, mode-2, fileCanvasPiTPC);
            PhiPiTPCYieldpTdiffMB[i][k] = PhiPiTPCYieldpTdiffMB[i][k] / deltay_axis[i] / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPiMB[i]*/ / (nEventsPhi * purityPhiMB);
            errPhiPiTPCYieldpTdiffMB[i][k] = errPhiPiTPCYieldpTdiffMB[i][k] / deltay_axis[i] / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPiMB[i]*/ / (nEventsPhi * purityPhiMB);
            //if (errPhiPiTPCYieldpTdiffMB[i][k] != errPhiPiTPCYieldpTdiffMB[i][k]) errPhiPiTPCYieldpTdiffMB[i][k] = 0.0;
            if (mode != 3) {
                PhiPiTPCYieldpTdiffMB[i][k] = PhiPiTPCYieldpTdiffMB[i][k] / purityPhiPiMB[i];
                errPhiPiTPCYieldpTdiffMB[i][k] = errPhiPiTPCYieldpTdiffMB[i][k] / purityPhiPiMB[i];
            }

            h1PhiPiTPCYieldMB[i]->SetBinContent(k+1, PhiPiTPCYieldpTdiffMB[i][k]);
            h1PhiPiTPCYieldMB[i]->SetBinError(k+1, errPhiPiTPCYieldpTdiffMB[i][k]);

            tie(PhiPiTOFYieldpTdiffMB[i][k], errPhiPiTOFYieldpTdiffMB[i][k]) = FitPhiPi(h2PhiInvMassPiNSigmaTOFMB[i][k], {i, k}, 1, mode-2, fileCanvasPiTOF);
            PhiPiTOFYieldpTdiffMB[i][k] = PhiPiTOFYieldpTdiffMB[i][k] / deltay_axis[i] / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPiMB[i]*/ / (nEventsPhi * purityPhiMB);
            errPhiPiTOFYieldpTdiffMB[i][k] = errPhiPiTOFYieldpTdiffMB[i][k] / deltay_axis[i] / (pTPi_axis[k+1] - pTPi_axis[k]) /*/ purityPhiPiMB[i]*/ / (nEventsPhi * purityPhiMB);
            //if (errPhiPiTOFYieldpTdiffMB[i][k] != errPhiPiTOFYieldpTdiffMB[i][k]) errPhiPiTOFYieldpTdiffMB[i][k] = 0.0;
            if (mode != 3) {
                PhiPiTOFYieldpTdiffMB[i][k] = PhiPiTOFYieldpTdiffMB[i][k] / purityPhiPiMB[i];
                errPhiPiTOFYieldpTdiffMB[i][k] = errPhiPiTOFYieldpTdiffMB[i][k] / purityPhiPiMB[i];
            }

            h1PhiPiTOFYieldMB[i]->SetBinContent(k+1, PhiPiTOFYieldpTdiffMB[i][k]);
            h1PhiPiTOFYieldMB[i]->SetBinError(k+1, errPhiPiTOFYieldpTdiffMB[i][k]);

            PhiPiYieldpTdiffMB[i][k] = (PhiPiTPCYieldpTdiffMB[i][k] + PhiPiTOFYieldpTdiffMB[i][k]) / 2.0;
            errPhiPiYieldpTdiffMB[i][k] = TMath::Sqrt(TMath::Power(errPhiPiTPCYieldpTdiffMB[i][k], 2) + TMath::Power(errPhiPiTOFYieldpTdiffMB[i][k], 2)) / 2.0;

            h1PhiPiYieldMB[i]->SetBinContent(k+1, PhiPiYieldpTdiffMB[i][k]);
            h1PhiPiYieldMB[i]->SetBinError(k+1, errPhiPiYieldpTdiffMB[i][k]);
        }
    }

    //fileCanvasK0S->Close();
    //fileCanvasPiTPC->Close();
    //fileCanvasPiTOF->Close();

    //********************************************************************************************

    array<TCanvas*, nbin_deltay> cPhiK0SYield = PlotHistograms(h1PhiK0SYield, h1PhiK0SYieldMB, outPath, "rawSpectrumK0SDY%i2D");
    array<TCanvas*, nbin_deltay> cPhiPiTPCYield = PlotHistograms(h1PhiPiTPCYield, h1PhiPiTPCYieldMB, outPath, "rawSpectrumPiTPCDY%i2D");
    array<TCanvas*, nbin_deltay> cPhiPiTOFYield = PlotHistograms(h1PhiPiTOFYield, h1PhiPiTOFYieldMB, outPath, "rawSpectrumPiTOFDY%i2D");

    //********************************************************************************************

    Double_t PhiK0SYieldpTint[nbin_deltay][nbin_mult] = {0}, errPhiK0SYieldpTint[nbin_deltay][nbin_mult] = {0};

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pTK0S; k++) {
                PhiK0SYieldpTint[i][j] += PhiK0SYieldpTdiff[i][j][k] * (pTK0S_axis[k+1] - pTK0S_axis[k]);
                errPhiK0SYieldpTint[i][j] += TMath::Power(errPhiK0SYieldpTdiff[i][j][k] * (pTK0S_axis[k+1] - pTK0S_axis[k]), 2);
            }
            errPhiK0SYieldpTint[i][j] = TMath::Sqrt(errPhiK0SYieldpTint[i][j]);
        }
    }

    TMultiGraph* mgK0S = new TMultiGraph();
    TGraphAsymmErrors* K0S[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        K0S[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), PhiK0SYieldpTint[i], errmult.data(), errmult.data(), errPhiK0SYieldpTint[i], errPhiK0SYieldpTint[i]);
        if (i == 0) {
            PlotFeatures(K0S[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgK0S);
        } else if (i == 1) {
            PlotFeatures(K0S[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgK0S);
        } else if (i == 2) {
            PlotFeatures(K0S[i], 33, kRed+1, 2, 1, kRed+1, 2, 3001, kRed+1, 0.4, mgK0S);
        }
    }

    TCanvas* cK0S = new TCanvas("cK0S", "cK0S", 1000, 800);
    cK0S->cd();
    //gPad->SetLogy();
    gStyle->SetOptStat(0);
    //gPad->SetMargin(0.16,0.01,0.13,0.06)
    mgK0S->SetTitle("; #LTdN_{ch}/d#eta#GT_{|#eta|<0.5} ; 1/N_{ev,#phi} dN_{K^{0}_{S}}/d#it{y}");
    mgK0S->Draw("AP");

    TLegend* legK0S1 = new TLegend(0.15, 0.82, 0.35, 0.85);
    legK0S1->SetHeader("#bf{This work}");
    legK0S1->SetTextSize(0.05);
    legK0S1->SetLineWidth(0);
    legK0S1->Draw("same");

    TLegend* legK0S2 = new TLegend(0.15, 0.62, 0.35, 0.82);
    legK0S2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legK0S2->AddEntry(K0S[0], "Inclusive", "p");
    legK0S2->AddEntry(K0S[1], "|#it{#Deltay}| < 0.5", "p");
    legK0S2->AddEntry(K0S[2], "|#it{#Deltay}| < 0.1", "p");
    legK0S2->SetTextSize(0.05);
    legK0S2->SetLineWidth(0);
    legK0S2->Draw("same");

    /*string outNameRawK0SYield = outPath + "rawK0SYield.root";
    cK0S->SaveAs(outNameRawK0SYield.c_str());
    outNameRawK0SYield = outPath + "rawK0SYield.pdf";
    cK0S->SaveAs(outNameRawK0SYield.c_str());*/

    //********************************************************************************************

    Double_t PhiPiYieldpTint[nbin_deltay][nbin_mult] = {0}, errPhiPiYieldpTint[nbin_deltay][nbin_mult] = {0};

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            for (int k = 0; k < nbin_pTPi; k++) {
                PhiPiYieldpTint[i][j] += PhiPiYieldpTdiff[i][j][k] * (pTPi_axis[k+1] - pTPi_axis[k]);
                errPhiPiYieldpTint[i][j] += TMath::Power(errPhiPiYieldpTdiff[i][j][k] * (pTPi_axis[k+1] - pTPi_axis[k]), 2);
            }
            errPhiPiYieldpTint[i][j] = TMath::Sqrt(errPhiPiYieldpTint[i][j]);
        }
    }

    TMultiGraph* mgPi = new TMultiGraph();
    TGraphAsymmErrors* Pi[nbin_deltay];

    for (int i = 0; i < nbin_deltay; i++) {
        Pi[i] = new TGraphAsymmErrors(nbin_mult, mult.data(), PhiPiYieldpTint[i], errmult.data(), errmult.data(), errPhiPiYieldpTint[i], errPhiPiYieldpTint[i]);
        if (i == 0) {
            PlotFeatures(Pi[i], 33, kBlack, 2, 1, kBlack, 2, 3001, kBlack, 0.4, mgPi);
        } else if (i == 1) {
            PlotFeatures(Pi[i], 33, kGreen+3, 2, 1, kGreen+3, 2, 3001, kGreen+3, 0.4, mgPi);
        } else if (i == 2) {
            PlotFeatures(Pi[i], 33, kRed+1, 2, 1, kRed+1, 2, 3001, kRed+1, 0.4, mgPi);
        }
    }

    TCanvas* cPi = new TCanvas("cPi", "cPi", 1000, 800);
    cPi->cd();
    //gPad->SetLogy();
    gStyle->SetOptStat(0);
    mgPi->SetTitle("; #LTdN_{ch}/d#eta#GT_{|#eta|<0.5} ; #frac{1}{N_{ev,#phi}} (#pi^{+}+#pi^{#minus})");
    mgPi->Draw("AP");

    TLegend* legPi1 = new TLegend(0.15, 0.82, 0.35, 0.85);
    legPi1->SetHeader("#bf{This work}");
    legPi1->SetTextSize(0.05);
    legPi1->SetLineWidth(0);
    legPi1->Draw("same");

    TLegend* legPi2 = new TLegend(0.15, 0.62, 0.35, 0.82);
    legPi2->SetHeader("pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5");
    legPi2->AddEntry(Pi[0], "Inclusive", "p");
    legPi2->AddEntry(Pi[1], "|#it{#Deltay}| < 0.5", "p");
    legPi2->AddEntry(Pi[2], "|#it{#Deltay}| < 0.1", "p");
    legPi2->SetTextSize(0.05);
    legPi2->SetLineWidth(0);
    legPi2->Draw("same");

    /*string outNameRawPiYield = outPath + "rawPiYield.root";
    cPi->SaveAs(outNameRawPiYield.c_str());
    outNameRawPiYield = outPath + "rawPiYield.pdf";
    cPi->SaveAs(outNameRawPiYield.c_str());*/

    //********************************************************************************************

    TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");
    outFile->cd();

    for (int i = 0; i < nbin_deltay; i++) {
        for (int j = 0; j < nbin_mult; j++) {
            h1PhiK0SYield[i][j]->Write();
            h1PhiPiTPCYield[i][j]->Write();
            h1PhiPiTOFYield[i][j]->Write();
        }
        h1PhiK0SYieldMB[i]->Write();
        h1PhiPiTPCYieldMB[i]->Write();
        h1PhiPiTOFYieldMB[i]->Write();
    }

    outFile->Close();

    return;
}
