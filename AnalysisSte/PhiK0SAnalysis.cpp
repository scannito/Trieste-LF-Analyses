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
    LFInvMassFitter PhiK0SFitter(objectHandlerPhiK0S.GetSetHisto2D(nbin_pT::K0S, "h2PhiK0SInvMass"), objectHandlerPhiK0S.GetSetHistoMultInt2D(nbin_pT::K0S, "h2PhiK0SInvMassMB"),
                                 objectHandlerPhiK0S.GetOutPath(), objectHandlerPhiK0S.GetOutFileName());

    ObjectHandler objectHandlerPhiPi(argv[2], requiredKeys);
    LFInvMassFitter PhiPiTPCFitter(objectHandlerPhiPi.GetSetHisto2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC"), objectHandlerPhiPi.GetSetHistoMultInt2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPCMB"),
                                   objectHandlerPhiPi.GetOutPath(), objectHandlerPhiPi.GetOutFileName());
    LFInvMassFitter PhiPiTOFFitter(objectHandlerPhiPi.GetSetHisto2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF"), objectHandlerPhiPi.GetSetHistoMultInt2D(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOFMB"),
                                   objectHandlerPhiPi.GetOutPath(), objectHandlerPhiPi.GetOutFileName());
}
