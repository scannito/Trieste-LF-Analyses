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
//#include "AnalysisUtils/Plot.h"

#include "LFInvMassFitter.h"

int main(int argc, char* argv[])
{
    std::cout << "Starting Phi K0S Analysis..." << std::endl;

    LFInvMassFitter PhiK0SFitter("K0S", "../AnalysisSte/data/PhiK0SAnalysisProjections.root", nbin_pT::K0S, "h2PhiK0SInvMass");
    //PhiK0SFitter.ExportYields("../AnalysisSte/data/PhiK0SYields.root", nbin_pT::K0S, pT_axis::K0S, "h1PhiK0SYield", 0);

    //LFInvMassFitter PhiPiTPCFitter("Pi", "../AnalysisSte/data/PhiPiTPCAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC");
    //PhiPiTPCFitter.ExportYields("../AnalysisSte/data/PhiPiTPCYields.root", nbin_pT::Pi, pT_axis::Pi, "h1PhiPiTPCYield", 0);

    //LFInvMassFitter PhiPiTOFFitter("Pi", "../AnalysisSte/data/PhiPiTOFAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF");
    //PhiPiTOFFitter.ExportYields("../AnalysisSte/data/PhiPiTOFYields.root", nbin_pT::Pi, pT_axis::Pi, "h1PhiPiTOFYield", 1);
}
