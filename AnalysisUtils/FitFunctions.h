#pragma once

#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"

#include <vector>
#include <array>
#include <string>

inline Double_t VoigtPoly(Double_t* x, Double_t* par)
{
    return par[0] * TMath::Voigt(x[0] - par[2], par[3], par[1], 4) + (par[4] + par[5] * x[0] + par[6] * TMath::Sqrt(x[0] - 0.987));
}

inline Double_t Voigt(Double_t* x, Double_t* par)
{
    return par[0] * TMath::Voigt(x[0] - par[2], par[3], par[1], 4);
}

inline Double_t Voigt1(Double_t* x, Double_t* par)
{
    return TMath::Voigt(x[0] - par[1], par[2], par[0], 4);
}

inline Double_t PolySqrt(Double_t* x, Double_t* par)
{
    return par[0] + par[1] * x[0] + par[2] * TMath::Sqrt(x[0] - 0.987); 
}

inline Double_t DoubleSidedCrystalBallPoly(Double_t *x, Double_t *par)
{ 
    Double_t alpha_l = par[0]; 
    Double_t alpha_h = par[1]; 
    Double_t n_l = par[2]; 
    Double_t n_h = par[3]; 
    Double_t mean = par[4]; 
    Double_t sigma = par[5];
    Double_t N = par[6];

    Double_t t = (x[0] - mean) / sigma;
    Double_t fact1TLowerAlphaL = alpha_l / n_l;
    Double_t fact2TLowerAlphaL = (n_l / alpha_l) - alpha_l - t;
    Double_t fact1THigherAlphaH = alpha_h / n_h;
    Double_t fact2THigherAlphaH = (n_h / alpha_h) - alpha_h + t;
    Double_t result;
    
    if (-alpha_l <= t && t <= alpha_h) {
        result = exp(-0.5 * TMath::Power(t, 2));
    } else if (t < -alpha_l) {
        result = exp(-0.5 * TMath::Power(alpha_l, 2)) * TMath::Power(fact1TLowerAlphaL * fact2TLowerAlphaL, -n_l);
    } else if (t > alpha_h) {
        result = exp(-0.5 * TMath::Power(alpha_h, 2)) * TMath::Power(fact1THigherAlphaH * fact2THigherAlphaH, -n_h);
    }

    return N * result + (par[7] + par[8] * x[0]);
}

inline Double_t DoubleSidedCrystalBall(Double_t *x, Double_t *par)
{ 
    Double_t alpha_l = par[0]; 
    Double_t alpha_h = par[1]; 
    Double_t n_l = par[2]; 
    Double_t n_h = par[3]; 
    Double_t mean = par[4]; 
    Double_t sigma = par[5];
    Double_t N = par[6];

    Double_t t = (x[0] - mean) / sigma;
    Double_t fact1TLowerAlphaL = alpha_l / n_l;
    Double_t fact2TLowerAlphaL = (n_l / alpha_l) - alpha_l - t;
    Double_t fact1THigherAlphaH = alpha_h / n_h;
    Double_t fact2THigherAlphaH = (n_h / alpha_h) - alpha_h + t;
    Double_t result;
    
    if (-alpha_l <= t && t <= alpha_h) {
        result = exp(-0.5 * TMath::Power(t, 2));
    } else if (t < -alpha_l) {
        result = exp(-0.5 * TMath::Power(alpha_l, 2)) * TMath::Power(fact1TLowerAlphaL * fact2TLowerAlphaL, -n_l);
    } else if (t > alpha_h) {
        result = exp(-0.5 * TMath::Power(alpha_h, 2)) * TMath::Power(fact1THigherAlphaH * fact2THigherAlphaH, -n_h);
    }

    return N * result;
}

inline Double_t DoubleSidedCrystalBall1(Double_t *x, Double_t *par)
{ 
    Double_t alpha_l = par[0]; 
    Double_t alpha_h = par[1]; 
    Double_t n_l = par[2]; 
    Double_t n_h = par[3]; 
    Double_t mean = par[4]; 
    Double_t sigma = par[5];

    Double_t t = (x[0] - mean) / sigma;
    Double_t fact1TLowerAlphaL = alpha_l / n_l;
    Double_t fact2TLowerAlphaL = (n_l / alpha_l) - alpha_l - t;
    Double_t fact1THigherAlphaH = alpha_h / n_h;
    Double_t fact2THigherAlphaH = (n_h / alpha_h) - alpha_h + t;
    Double_t result;
    
    if (-alpha_l <= t && t <= alpha_h) {
        result = exp(-0.5 * TMath::Power(t, 2));
    } else if (t < -alpha_l) {
        result = exp(-0.5 * TMath::Power(alpha_l, 2)) * TMath::Power(fact1TLowerAlphaL * fact2TLowerAlphaL, -n_l);
    } else if (t > alpha_h) {
        result = exp(-0.5 * TMath::Power(alpha_h, 2)) * TMath::Power(fact1THigherAlphaH * fact2THigherAlphaH, -n_h);
    }

    return result;
}

inline Double_t Poly1(Double_t* x, Double_t* par)
{
    return par[0] + par[1] * x[0];
}

inline Double_t PhiInvMassK0SNSigmadEdx(Double_t* x, Double_t* par) 
{
    Double_t K0Ssig = DoubleSidedCrystalBall1(&x[0], &par[0]);
    Double_t K0Sbkg = Poly1(&x[0], &par[6]);

    Double_t Phisig = Voigt1(&x[1], &par[8]);
    Double_t Phibkg = PolySqrt(&x[1], &par[11]);

    return par[14] * Phisig * K0Ssig + par[15] * Phisig * K0Sbkg + par[16] * Phibkg * K0Ssig + par[17] * Phibkg * K0Sbkg;
}

inline Double_t PhiInvMassK0SNSigmadEdxSig(Double_t* x, Double_t* par) 
{
    Double_t K0Ssig = DoubleSidedCrystalBall1(&x[0], &par[0]);
    Double_t Phisig = Voigt1(&x[1], &par[6]);

    return par[9] * Phisig * K0Ssig;
}
