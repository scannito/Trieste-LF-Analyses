#pragma once

#include <vector>

struct AxisCut {
    int axis;
    int binLow;
    int binHigh;

    AxisCut(int a, int l, int h) : axis(a), binLow(l), binHigh(h) {}
};

// Generate all cut combinations, according to each desidered upper and lower bin
inline std::vector<std::vector<AxisCut>> ExpandAxisCuts(const std::vector<AxisCut>& cuts) {
    std::vector<std::vector<AxisCut>> result;

    auto recurse = [&](auto&& self, size_t i, std::vector<AxisCut> current) -> void {
        if (i == cuts.size()) {
            result.push_back(std::move(current));
            return;
        }

        const AxisCut& cut = cuts[i];
        for (int b = cut.binLow; b <= cut.binHigh; ++b) {
            auto tmp = current;
            tmp.emplace_back(cut.axis, b, b);
            self(self, i + 1, std::move(tmp));
        }
    };

    recurse(recurse, 0, {});
    return result;
}

struct VectorHash {
    std::size_t operator()(const std::vector<int>& vec) const {
        std::size_t seed = vec.size();
        for (auto& i : vec) {
            seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};