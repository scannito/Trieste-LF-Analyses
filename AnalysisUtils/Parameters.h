#include "TROOT.h"

#include <vector>
#include <array>
#include <string>

const std::vector<Int_t> Colors = {634, 628, 807, 797, kOrange - 4, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
const std::vector<Int_t> ColorsFinal = {kBlack, kBlue, kGreen+3, 797, kRed+1};
const std::vector<Int_t> FullMarkers  = {20, 21, 33, 34, 29, 41, 47, 43};
const std::vector<Int_t> EmptyMarkers = {53, 56, 57, 58, 64, 67, 54, 65};
const std::vector<Int_t> Markers = {20, 21, 33, 34, 29, 41, 47, 43, 53, 56, 57, 58, 64, 67, 54, 65};

constexpr Int_t nbin_deltay = 5, nbin_deltay_red = 3, nbin_mult = 10, nbin_pTK0S = 9, nbin_pTPi = 11, nbin_massPhi = 13, nbin_pTPhi = 7;

const std::array<Double_t, nbin_deltay> deltay_axis = {1.0, 0.8, 0.5, 0.3, 0.1};
const std::array<Double_t, nbin_mult+1> mult_axis = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
const std::array<Double_t, nbin_pTK0S+1> pTK0S_axis = {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0};
const std::array<Double_t, nbin_pTPi+1> pTPi_axis = {0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0};

const std::array<Double_t, nbin_mult> mult = {20.3, 17.1, 14.6, 12.9, 11.6, 10.1, 8.5, 7.2, 5.7, 3.8};
const std::array<Double_t, nbin_mult> errmult = {0.5, 0.4, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.1, 0.1};

constexpr Double_t topPadHeight = 0.7;
constexpr Double_t bottomPadHeight = 1.0 - topPadHeight;
constexpr Double_t scaleFactor = bottomPadHeight / topPadHeight;

constexpr Double_t lowmPhiPur = 1.0095, upmPhiPur = 1.029;

const std::array<std::string, nbin_deltay> yCut = {"Inclusive", "|#Delta#it{y}| < 0.5", "|#Delta#it{y}| < 0.1"};
