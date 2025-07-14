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
    std::vector<AxisCut> slicingK0S = {  };
    THnSparseProjectorPhiK0S.ExportProjections(nbin_pT::K0S, slicingK0S, "h2PhiK0SInvMass", {4, 3});

    THnSparseProjector THnSparseProjectorPhiPi(argv[2], requiredKeys);
    std::vector<AxisCut> slicingPi = {  };
    THnSparseProjectorPhiPi.ExportProjections(nbin_pT::Pi, slicingPi, "h2PhiInvMassPiNSigmaTPC", {5, 3});
    THnSparseProjectorPhiPi.ExportProjections(nbin_pT::Pi, slicingPi, "h2PhiInvMassPiNSigmaTOF", {5, 4});

    return 0;
}

int runLFInvMassFitter(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: fitter <inputFile.json>\n";
        return 1;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "objectsPath", "outputFile"};

    std::vector<AxisCut> slicingK0S = {  };
    LFInvMassFitter PhiK0SFitter("K0S", "../../AnalysisSte/data/PhiK0SAnalysisProjections.root", requiredKeys, slicingK0S, nbin_pT::K0S, "h2PhiK0SInvMass");
    PhiK0SFitter.ExportYields("../../AnalysisSte/data/PhiK0SYields.root", nbin_pT::K0S, pT_axis::K0S, "h1PhiK0SYield", 0);

    //LFInvMassFitter PhiPiTPCFitter("Pi", "../../AnalysisSte/data/PhiPiTPCAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTPC");
    //PhiPiTPCFitter.ExportYields("../../AnalysisSte/data/PhiPiTPCYields.root", nbin_pT::Pi, pT_axis::Pi, "h1PhiPiTPCYield", 0);

    //LFInvMassFitter PhiPiTOFFitter("Pi", "../../AnalysisSte/data/PhiPiTOFAnalysisProjections.root", nbin_pT::Pi, "h2PhiInvMassPiNSigmaTOF");
    //PhiPiTOFFitter.ExportYields("../../AnalysisSte/data/PhiPiTOFYields.root", nbin_pT::Pi, pT_axis::Pi, "h1PhiPiTOFYield", 1);
    
    return 0;
}

int runEfficiencyHandler(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: efficiency <inputPhi.json> <inputK0S.json> <inputPi.json>\n";
        return 1;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "recoPath", "genPath", "genAssocRecoPath", "outputFile"};

    EfficiencyHandler PhiEfficiencyHandler(argv[1], requiredKeys);
    PhiEfficiencyHandler.ExportCorrections();
    //PhiEfficiencyHandler.ExportCorrectionsForCCDB();

    //EfficiencyHandler K0SEfficiencyHandler(argv[2], requiredKeys);
    //K0SEfficiencyHandler.ExportCorrections();
    //K0SEfficiencyHandler.ExportCorrectionsForCCDB();

    //EfficiencyHandler PiEfficiencyHandler(argv[3], requiredKeys);
    //PiEfficiencyHandler.ExportCorrections();
    //PiEfficiencyHandler.ExportCorrectionsForCCDB();

    return 0;
}
