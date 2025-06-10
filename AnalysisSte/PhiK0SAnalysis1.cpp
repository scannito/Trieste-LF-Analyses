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
    objectHandlerPhiK0S.ExportProjections("../AnalysisSte/data/PhiK0SAnalysisProjections.root", nbin_pT::K0S, "h2PhiK0SInvMass", {4, 3});

    ObjectHandler objectHandlerPhiPi(argv[2], requiredKeys);
    objectHandlerPhiPi.ExportProjections("../AnalysisSte/data/PhiPiTPCAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC", {5, 3});
    objectHandlerPhiPi.ExportProjections("../AnalysisSte/data/PhiPiTOFAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF", {5, 4});
}
