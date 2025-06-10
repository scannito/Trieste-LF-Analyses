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

#include "AnalysisUtils/Parameters.h"

#include "ObjectHandler.h"
#include "LFInvMassFitter.h"

int main(int argc, char* argv[])
{
    LFInvMassFitter PhiK0SFitter("../AnalysisSte/data/PhiK0SAnalysisProjections.root", nbin_pT::K0S, "h2PhiK0SInvMass");

    LFInvMassFitter PhiPiTPCFitter("../AnalysisSte/data/PhiPiTPCAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC");
    LFInvMassFitter PhiPiTOFFitter("../AnalysisSte/data/PhiPiTOFAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF");
}
