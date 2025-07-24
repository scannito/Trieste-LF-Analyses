#include "Riostream.h"
#include "TApplication.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
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

#include <array>
#include <string>
#include <vector>

#include "THnSparseProjector.h"
#include "LFInvMassFitter.h"
#include "EfficiencyHandler.h"

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/Plot.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/AxesUtils.h"

// Main logics declaration
void runTHnSparseProjector(int argc, char* argv[]);
void runLFInvMassFitter(int argc, char* argv[]);
void runEfficiencyHandler(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    int argcCopy = argc;
    char** argvCopy = new char*[argcCopy];
    for (int i = 0; i < argcCopy; ++i) {
        argvCopy[i] = argv[i];
    }

    TApplication app("PhiK0SAnalysis", &argc, argv);

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
        runTHnSparseProjector(argc - 1, &argv[1]);  // Shift args so that argv[0] = "projector/fitter/effiency"
    } else if (mode == "fitter") {
        runLFInvMassFitter(argc - 1, &argv[1]);
    } else if (mode == "efficiency") {
        runEfficiencyHandler(argcCopy - 1, &argvCopy[1]);
    } else {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 1;
    }
    app.Run();
}

void runTHnSparseProjector(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: projector <inputK0S.json> <inputPi.json>\n";
        return;
    }

    std::vector<std::string> requiredKeys = {"inputFile", "objectPath", "outputFile"};

    std::vector<AxisCut> slicingK0S = { {0, 1, nbin_deltay}, {1, 1, nbin_mult}, {2, -1, nbin_pT::K0S}  };
    THnSparseProjector THnSparseProjectorPhiK0S(argv[1], requiredKeys);
    THnSparseProjectorPhiK0S.ExportProjections(nbin_pT::K0S, slicingK0S, "h2PhiK0SInvMass", {4, 3});

    std::vector<AxisCut> slicingPion = { {0, 1, nbin_deltay}, {1, 1, nbin_mult}, {2, -1, nbin_pT::Pi} };
    THnSparseProjector THnSparseProjectorPhiPiTPC(argv[2], requiredKeys);
    THnSparseProjectorPhiPiTPC.ExportProjections(nbin_pT::Pi, slicingPion, "h2PhiInvMassPiNSigmaTPC", {4, 3});
    THnSparseProjector THnSparseProjectorPhiPiTOF(argv[3], requiredKeys);
    THnSparseProjectorPhiPiTOF.ExportProjections(nbin_pT::Pi, slicingPion, "h2PhiInvMassPiNSigmaTOF", {4, 3});
}

void runLFInvMassFitter(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: fitter <inputFile.json>\n";
        return;
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
}

void runEfficiencyHandler(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: efficiency <inputPhi.json> <inputK0S.json> <inputPi.json>\n";
        return;
    }

    std::vector<std::string> requiredKeys1 = {"inputFile", "recoPath", "genPath", "genAssocRecoPath", "outputFile"};

    EfficiencyHandler PhiEfficiencyHandler("Phi", argv[1], requiredKeys1);
    //PhiEfficiencyHandler.ExportCorrectionsForCCDB();

    EfficiencyHandler K0SEfficiencyHandler("K0S", argv[2], requiredKeys1);
    //K0SEfficiencyHandler.ExportCorrectionsForCCDB();

    std::vector<std::string> requiredKeys2 = {"inputFile", "recoPath", "recoPath2", "genPath", "genAssocRecoPath", "outputFile"};

    EfficiencyHandler PiEfficiencyHandler("Pion", argv[3], requiredKeys2);
    //PiEfficiencyHandler.ExportCorrectionsForCCDB();

    std::vector<double> rebinnedpTAxisPhi{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                                          1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                          2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                          3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                                          4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.4, 5.8, 6.2, 6.6, 
                                          7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0};
    std::vector<double> rebinnedpTAxisK0S{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                                          1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                          2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                          3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                                          4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
                                          5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 
                                          7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 10.0};

    std::array<TH1*, nbin_mult> h1EfficiencyPhi;
    std::array<TH1*, nbin_mult> h1SignalLossPhi;
    std::array<TH1*, nbin_mult> h1EffXSigLossPhi;
    std::array<TH1*, nbin_mult> h1EfficiencyK0S;
    std::array<TH1*, nbin_mult> h1SignalLossK0S;
    std::array<TH1*, nbin_mult> h1EffXSigLossK0S;
    std::array<TH1*, nbin_mult> h1Efficiency1Pi;
    std::array<TH1*, nbin_mult> h1Efficiency2Pi;
    std::array<TH1*, nbin_mult> h1SignalLossPi;
    std::array<TH1*, nbin_mult> h1EffXSigLossPi;
    std::array<TH1*, nbin_mult> h1CombinedEfficiencyPi;
    std::array<TH1*, nbin_mult> h1CombEffXSigLossPi;

    for (int i = 0; i < nbin_mult; ++i) {
        /*h1EfficiencyPhi[i] = PhiEfficiencyHandler.GetEfficiencySpectrum(i);
        h1SignalLossPhi[i] = PhiEfficiencyHandler.GetSignalLossSpectrum(i);
        h1EffXSigLossPhi[i] = PhiEfficiencyHandler.GetEffXSigLossSpectrum(i);
        h1EfficiencyK0S[i] = K0SEfficiencyHandler.GetEfficiencySpectrum(i);
        h1SignalLossK0S[i] = K0SEfficiencyHandler.GetSignalLossSpectrum(i);
        h1EffXSigLossK0S[i] = K0SEfficiencyHandler.GetEffXSigLossSpectrum(i);*/

        bool rebin = true; // Set to true to use rebinned pT axis

        h1EfficiencyPhi[i] = PhiEfficiencyHandler.GetEfficiencySpectrum(i, rebin, rebinnedpTAxisPhi);
        h1SignalLossPhi[i] = PhiEfficiencyHandler.GetSignalLossSpectrum(i, rebin, rebinnedpTAxisPhi);
        h1EffXSigLossPhi[i] = PhiEfficiencyHandler.GetEffXSigLossSpectrum(i, rebin, rebinnedpTAxisPhi);
        h1EfficiencyK0S[i] = K0SEfficiencyHandler.GetEfficiencySpectrum(i, rebin, rebinnedpTAxisK0S);
        h1SignalLossK0S[i] = K0SEfficiencyHandler.GetSignalLossSpectrum(i, rebin, rebinnedpTAxisK0S);
        h1EffXSigLossK0S[i] = K0SEfficiencyHandler.GetEffXSigLossSpectrum(i, rebin, rebinnedpTAxisK0S);

        h1Efficiency1Pi[i] = PiEfficiencyHandler.GetEfficiencySpectrum(i);
        h1Efficiency2Pi[i] = PiEfficiencyHandler.GetEfficiencySpectrum2(i);
        h1SignalLossPi[i] = PiEfficiencyHandler.GetSignalLossSpectrum(i);
        h1EffXSigLossPi[i] = PiEfficiencyHandler.GetEffXSigLossSpectrum(i);
        h1CombinedEfficiencyPi[i] = PiEfficiencyHandler.GetCombinedEfficiencySpectrum(i);
        h1CombEffXSigLossPi[i] = PiEfficiencyHandler.GetCombEffXSigLossSpectrum(i);
    }

    std::vector<TCanvas*> cEfficiencyPhi = PlotSpectra({ h1EfficiencyPhi });
    std::vector<TCanvas*> cSignalLossPhi = PlotSpectra({ h1SignalLossPhi });
    std::vector<TCanvas*> cCombEfficiencyPhi = PlotSpectra({ h1EffXSigLossPhi });
    std::vector<TCanvas*> cEfficiencyK0S = PlotSpectra({ h1EfficiencyK0S });
    std::vector<TCanvas*> cSignalLossK0S = PlotSpectra({ h1SignalLossK0S });
    std::vector<TCanvas*> cCombEfficiencyK0S = PlotSpectra({ h1EffXSigLossK0S });
    std::vector<TCanvas*> cEfficiency1Pi = PlotSpectra({ h1Efficiency1Pi });
    std::vector<TCanvas*> cEfficiency2Pi = PlotSpectra({ h1Efficiency2Pi });
    std::vector<TCanvas*> cSignalLossPi = PlotSpectra({ h1SignalLossPi });
    std::vector<TCanvas*> cCombEfficiencyPi = PlotSpectra({ h1EffXSigLossPi });
    std::vector<TCanvas*> cCombinedEfficiencyPi = PlotSpectra({ h1CombinedEfficiencyPi });
    std::vector<TCanvas*> cCombEfficiency2Pi = PlotSpectra({ h1CombEffXSigLossPi });
}
