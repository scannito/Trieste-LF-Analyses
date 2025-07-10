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

#include "THnSparseProjector.h"

int main(int argc, char* argv[])
{
    std::cout << "Starting Phi K0S Analysis..." << std::endl;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.json>" << std::endl;
        return 1;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "objectPath", "outputFile"};

    THnSparseProjector THnSparseProjectorPhiK0S(argv[1], requiredKeys);
    THnSparseProjectorPhiK0S.ExportProjections(nbin_pT::K0S, "h2PhiK0SInvMass", {4, 3});

    THnSparseProjector THnSparseProjectorPhiPi(argv[2], requiredKeys);
    THnSparseProjectorPhiPi.ExportProjections(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC", {5, 3});
    THnSparseProjectorPhiPi.ExportProjections(nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF", {5, 4});
}
