#pragma once

#include <map>
#include <string>

class JSONReader 
{
public:
    JSONReader() = default;
    JSONReader(const char* filename, const std::vector<std::string>& requiredKeys);
    ~JSONReader() = default;

    std::map<std::string, std::string> GetMeta() const { return mMeta; }

private:
    std::map<std::string, std::string> mMeta;
};
