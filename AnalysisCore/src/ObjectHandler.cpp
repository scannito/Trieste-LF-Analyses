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

#include <utility>
#include <vector>
#include <string>

#include "AnalysisUtils/Parameters.h"
#include "AnalysisUtils/FitFunctions.h"
#include "AnalysisUtils/Plot.h"
#include "AnalysisUtils/Projections.h"

#include "AnalysisCore/include/ObjectHandler.h"

ObjectHandler::ObjectHandler() : fName("", ""), fFile(nullptr) {}
ObjectHandler::~ObjectHandler() {
    if (fFile) {
        fFile->Close();
        delete fFile;
    }
}
