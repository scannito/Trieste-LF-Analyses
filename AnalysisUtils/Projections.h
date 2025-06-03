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

#include <string>

TH2F* Project2D(THnSparseF* hn, Int_t axistocut, Int_t binlow, Int_t binup, Int_t axistoproj1, Int_t axistoproj2, Option_t* option = "", std::string hname = "") 
{ 
    if (!hn) return 0;
    hn->GetAxis(axistocut)->SetRange(binlow, binup);
    TH2F* h2 = (TH2F*)hn->Projection(axistoproj1, axistoproj2, option);
    h2->SetName(hname.c_str());
    h2->SetDirectory(0);
    return h2;
}

TH2F* Project2D(THnSparseF* hn, Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistoproj1, Int_t axistoproj2, Option_t* option = "", std::string hname = "") 
{ 
    if (!hn) return 0;
    hn->GetAxis(axistocut1)->SetRange(binlow1, binup1);
    hn->GetAxis(axistocut2)->SetRange(binlow2, binup2);
    TH2F* h2 = (TH2F*)hn->Projection(axistoproj1, axistoproj2, option);
    h2->SetName(hname.c_str());
    h2->SetDirectory(0);
    return h2;
}
TH1F* Project1D(THnSparseF* hn, Int_t axistocut1, Int_t binlow1, Int_t binup1, Int_t axistocut2, Int_t binlow2, Int_t binup2, Int_t axistocut3, Int_t binlow3, Int_t binup3, Int_t axistoproj, Option_t* option = "", std::string hname = "") 
{ 
    if (!hn) return 0;
    hn->GetAxis(axistocut1)->SetRange(binlow1, binup1);
    hn->GetAxis(axistocut2)->SetRange(binlow2, binup2);
    hn->GetAxis(axistocut3)->SetRange(binlow3, binup3);
    TH1F* h1 = (TH1F*)hn->Projection(axistoproj, option);
    h1->SetName(hname.c_str());
    h1->SetDirectory(0);
    return h1;
}