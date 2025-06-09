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
    std::cout << "Starting Phi K0S Analysis..." << std::endl;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.json>" << std::endl;
        return 1;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "analysisDir", "eventHistDir", "eventHistName", "binEventHistName", 
                                             "PhiAssocDir", "PhiAssocInvMassHistName", "outputPath", "outputFile"};

    ObjectHandler objectHandlerPhiK0S(argv[1], requiredKeys);
    objectHandlerPhiK0S.CheckValidMembers();
    std::array<std::array<std::vector<TH2*>, nbin_mult>, nbin_deltay> setHisto2D = objectHandlerPhiK0S.GetSetHisto2D(nbin_pT::K0S, "h2PhiK0SInvMass", {4, 3});
    std::array<std::vector<TH2*>, nbin_deltay> setHistoMultInt2D = objectHandlerPhiK0S.GetSetHistoMultInt2D(nbin_pT::K0S, "h2PhiK0SInvMassMB", {4, 3});
    std::string outPath = objectHandlerPhiK0S.GetOutPath();
    std::string outFileName = objectHandlerPhiK0S.GetOutFileName();
    LFInvMassFitter PhiK0SFitter(objectHandlerPhiK0S.GetSetHisto2D(nbin_pT::K0S, "h2PhiK0SInvMass", {4, 3}), objectHandlerPhiK0S.GetSetHistoMultInt2D(nbin_pT::K0S, "h2PhiK0SInvMassMB", {4, 3}),
                                 objectHandlerPhiK0S.GetOutPath(), objectHandlerPhiK0S.GetOutFileName());

    ObjectHandler objectHandlerPhiPi(argv[2], requiredKeys);
    LFInvMassFitter PhiPiTPCFitter(objectHandlerPhiPi.GetSetHisto2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC", {5, 3}), objectHandlerPhiPi.GetSetHistoMultInt2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPCMB", {5, 3}),
                                   objectHandlerPhiPi.GetOutPath(), objectHandlerPhiPi.GetOutFileName());
    LFInvMassFitter PhiPiTOFFitter(objectHandlerPhiPi.GetSetHisto2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF", {5, 4}), objectHandlerPhiPi.GetSetHistoMultInt2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOFMB", {5, 4}),
                                   objectHandlerPhiPi.GetOutPath(), objectHandlerPhiPi.GetOutFileName());
}
