#include "Riostream.h"

#include <rapidjson/document.h>

#include "JSONReader.h"

JSONReader::JSONReader(const char* filename, const std::vector<std::string>& requiredKeys)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for reading." << std::endl;
        return;
    }

    std::string jsonStr((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    // Parse JSON
    rapidjson::Document document;
    document.Parse(jsonStr.data());

    if (document.HasParseError()) {
        std::cerr << "Error parsing JSON" << std::endl;
        return;
    }

    // Convert JSON to std::map
    for (auto itr = document.MemberBegin(); itr != document.MemberEnd(); ++itr) {
        mMeta[itr->name.GetString()] = itr->value.GetString();
    }

    bool allKeysPresent = true;
    for (const auto& key : requiredKeys) {
        if (mMeta.find(key) == mMeta.end()) {
            std::cerr << "Missing required key: " << key << std::endl;
            allKeysPresent = false;
        }
    }

    if (!allKeysPresent) {
        std::cerr << "Please provide all required keys in the JSON file" << std::endl;
    }
}
