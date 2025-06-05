#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

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

    // Read JSON from file
    const char* filename = argv[1];
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for reading." << std::endl;
        return 1;
    }

    std::string jsonStr((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    // Parse JSON
    rapidjson::Document document;
    document.Parse(jsonStr.data());

    if (document.HasParseError()) {
        std::cerr << "Error parsing JSON" << std::endl;
        return 1;
    }

    // Convert JSON to std::map
    std::map<std::string, std::string> meta;
    for (auto itr = document.MemberBegin(); itr != document.MemberEnd(); ++itr) {
        meta[itr->name.GetString()] = itr->value.GetString();
    }

    std::vector<std::string> requiredKeys = {"inputFile", "analysisDir", "eventHistDir", "eventHistName", "binEventHistName", 
                                             "PhiK0SDir", "PhiK0SInvMassHistName", "PhiPiDir", "PhiPiInvMassHistName",
                                             "outputPath", "outputFile"};

    bool allKeysPresent = true;
    for (const auto& key : requiredKeys) {
        if (meta.find(key) == meta.end()) {
            std::cerr << "Missing required key: " << key << std::endl;
            allKeysPresent = false;
        }
    }

    if (!allKeysPresent) {
        std::cerr << "Please provide all required keys in the JSON file" << std::endl;
        return 1;
    }

    TFile* file1 = TFile::Open(meta["inputFile"].data());
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error opening input file: " << meta["inputFile"] << std::endl;
        return 1;
    }
    TDirectoryFile* phik0shortanalysis = (TDirectoryFile*)file1->Get(meta["analysisDir"].data());
    if (!phik0shortanalysis) {
        std::cerr << "Error retrieving analysis directory: " << meta["analysisDir"] << std::endl;
        return 1;
    }

    TDirectoryFile* eventHist = (TDirectoryFile*)phik0shortanalysis->Get(meta["eventHistDir"].data());
    if (!eventHist) {
        std::cerr << "Error retrieving event histogram directory: " << meta["eventHistDir"] << std::endl;
        return 1;
    }
    std::string eventHistName = meta["eventHistName"];
    std::string binEventHistName = meta["binEventHistName"];

    TDirectoryFile* PhiK0SHist = (TDirectoryFile*)phik0shortanalysis->Get(meta["PhiK0SDir"].data());
    if (!PhiK0SHist) {
        std::cerr << "Error retrieving Phi K0S histogram directory: " << meta["PhiK0SDir"] << std::endl;
        return 1;
    }
    std::string PhiK0SInvMassHistName = meta["PhiK0SInvMassHistName"];

    TDirectoryFile* PhiPiHist = (TDirectoryFile*)phik0shortanalysis->Get(meta["PhiPiDir"].data());
    if (!PhiPiHist) {
        std::cerr << "Error retrieving Phi Pi histogram directory: " << meta["PhiPiDir"] << std::endl;
        return 1;
    }
    std::string PhiPiInvMassHistName = meta["PhiPiInvMassHistName"];

    std::string outPath = meta["outputPath"];
    std::string outFileName = meta["outputFile"];

    LFInvMassFitter PhiK0SFitter;

}
