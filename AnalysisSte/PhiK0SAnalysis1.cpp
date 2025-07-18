#include "Riostream.h"
#include "TApplication.h"
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

#include "THnSparseProjector.h"
#include "LFInvMassFitter.h"
#include "EfficiencyHandler.h"

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/Plot.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/AxesUtils.h"

// Main logics declaration
int runTHnSparseProjector(int argc, char* argv[]);
int runLFInvMassFitter(int argc, char* argv[]);
int runEfficiencyHandler(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    //TApplication app("PhiK0SAnalysis", &argc, argv);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mode> [args...]\n"
                  << "Available modes:\n"
                  << "  projector <input1.json> <input2.json>\n"
                  << "  fitter\n";
        return 1;
    }

    std::string mode = argv[1];

    std::cout << "Executing step " << mode << " of Phi-K0S rapidity correlation analysis\n";

    if (mode == "projector") {
        return runTHnSparseProjector(argc - 1, &argv[1]);  // Shift args so that argv[0] = "projector/fitter/effiency"
    } else if (mode == "fitter") {
        return runLFInvMassFitter(argc - 1, &argv[1]);
    } else if (mode == "efficiency") {
        return runEfficiencyHandler(argc - 1, &argv[1]);
    } else {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 1;
    }

    //app.Run();
}

int runTHnSparseProjector(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: projector <inputK0S.json> <inputPi.json>\n";
        return 1;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "objectPath", "outputFile"};

    THnSparseProjector THnSparseProjectorPhiK0S(argv[1], requiredKeys);
    std::vector<AxisCut> slicingK0S = { {0, 1, nbin_deltay}, {1, 1, nbin_mult}, {2, -1, nbin_pT::K0S}  };
    THnSparseProjectorPhiK0S.ExportProjections(nbin_pT::K0S, slicingK0S, "h2PhiK0SInvMass", {4, 3});

    THnSparseProjector THnSparseProjectorPhiPi(argv[2], requiredKeys);
    std::vector<AxisCut> slicingPion = { {0, 1, nbin_deltay}, {1, 1, nbin_mult}, {2, -1, nbin_pT::Pi} };
    THnSparseProjectorPhiPi.ExportProjections(nbin_pT::Pi, slicingPion, "h2PhiInvMassPiNSigmaTPC", {5, 3});
    THnSparseProjectorPhiPi.ExportProjections(nbin_pT::Pi, slicingPion, "h2PhiInvMassPiNSigmaTOF", {5, 4});

    return 0;
}

int runLFInvMassFitter(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: fitter <inputFile.json>\n";
        return 1;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "objectsPath", "outputFile"};

    std::vector<AxisCut> slicingK0S = { {0, 1, nbin_deltay}, {1, 1, nbin_mult} };
    LFInvMassFitter PhiK0SFitter("Phi-K0S", argv[1], requiredKeys, slicingK0S, nbin_pT::K0S, "h2PhiK0SInvMass");
    PhiK0SFitter.ExportYields(nbin_pT::K0S, pT_axis::K0S, "h1PhiK0SYield", 0);

    std::vector<AxisCut> slicingPion = {  };
    LFInvMassFitter PhiPiTPCFitter("Phi-Pion", argv[2], requiredKeys, slicingPion, nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC");
    PhiPiTPCFitter.ExportYields(nbin_pT::Pi, pT_axis::Pi, "h1PhiPiTPCYield", 0);
    LFInvMassFitter PhiPiTOFFitter("Phi-Pion", argv[2], requiredKeys, slicingPion, nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF");
    PhiPiTOFFitter.ExportYields(nbin_pT::Pi, pT_axis::Pi, "h1PhiPiTOFYield", 1);
    
    return 0;
}

int runEfficiencyHandler(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: efficiency <inputPhi.json> <inputK0S.json> <inputPi.json>\n";
        return 1;
    }

    std::vector<std::string> requiredKeys1 = {"inputFile", "recoPath", "genPath", "genAssocRecoPath", "outputFile"};

    EfficiencyHandler PhiEfficiencyHandler("Phi", argv[1], requiredKeys1);
    //PhiEfficiencyHandler.ExportCorrections();
    PhiEfficiencyHandler.ExportCorrectionsForCCDB();

    EfficiencyHandler K0SEfficiencyHandler("K0S", argv[2], requiredKeys1);
    //K0SEfficiencyHandler.ExportCorrections();
    K0SEfficiencyHandler.ExportCorrectionsForCCDB();

    std::vector<std::string> requiredKeys2 = {"inputFile", "recoPath", "recoPath2", "genPath", "genAssocRecoPath", "outputFile"};

    EfficiencyHandler PiEfficiencyHandler("Pion", argv[3], requiredKeys2);
    //PiEfficiencyHandler.ExportCorrections();
    PiEfficiencyHandler.ExportCorrectionsForCCDB();

    return 0;
}
