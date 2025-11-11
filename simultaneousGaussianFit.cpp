// simultaneousGaussianFit.cpp
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TMinuit.h"
#include <iostream>

TH1F *h1 = nullptr;
TH1F *h2 = nullptr;

// Global pointer for TMinuit
TMinuit *gMinuit = nullptr;

// Number of parameters: mu, sigma, A1, A2, b10, b11, b20, b21
const int NPAR = 8;

// Chi2 function
void chi2Func(int &npar, double *grad, double &fval, double *par, int flag) {
    double mu    = par[0];
    double sigma = par[1];
    double A1    = par[2];
    double A2    = par[3];
    double b10   = par[4];
    double b11   = par[5];
    double b20   = par[6];
    double b21   = par[7];

    double chi2 = 0.0;

    // Histogram 1
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double x = h1->GetBinCenter(i);
        double y = h1->GetBinContent(i);
        double err = h1->GetBinError(i);
        if (err == 0) err = 1.0; // avoid division by zero
        double model = b10 + b11*x + A1*TMath::Gaus(x, mu, sigma, true);
        chi2 += (y - model)*(y - model)/(err*err);
    }

    // Histogram 2
    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        double x = h2->GetBinCenter(i);
        double y = h2->GetBinContent(i);
        double err = h2->GetBinError(i);
        if (err == 0) err = 1.0;
        double model = b20 + b21*x + A2*TMath::Gaus(x, mu, sigma, true);
        chi2 += (y - model)*(y - model)/(err*err);
    }

    fval = chi2;
}

int main() {
    gStyle->SetOptStat(0);
    
    TFile *f = TFile::Open("experiments.root");
    if (!f) {
        std::cerr << "Could not open experiments.root!" << std::endl;
        return 1;
    }

    h1 = (TH1F*)f->Get("hexp1");
    h2 = (TH1F*)f->Get("hexp2");
    if (!h1 || !h2) {
        std::cerr << "Histograms not found!" << std::endl;
        return 1;
    }

    // Create TMinuit object
    gMinuit = new TMinuit(NPAR);
    gMinuit->SetFCN(chi2Func);

    // Set initial parameter guesses: mu, sigma, A1, A2, b10, b11, b20, b21
    double vstart[NPAR] = {5.0, 1.0, 100.0, 100.0, 10.0, 0.0, 10.0, 0.0};
    double step[NPAR]   = {0.1, 0.1, 10.0, 10.0, 1.0, 0.1, 1.0, 0.1};

    gMinuit->DefineParameter(0, "mu",    vstart[0], step[0], 0, 0);
    gMinuit->DefineParameter(1, "sigma", vstart[1], step[1], 0, 0);
    gMinuit->DefineParameter(2, "A1",    vstart[2], step[2], 0, 0);
    gMinuit->DefineParameter(3, "A2",    vstart[3], step[3], 0, 0);
    gMinuit->DefineParameter(4, "b10",   vstart[4], step[4], 0, 0);
    gMinuit->DefineParameter(5, "b11",   vstart[5], step[5], 0, 0);
    gMinuit->DefineParameter(6, "b20",   vstart[6], step[6], 0, 0);
    gMinuit->DefineParameter(7, "b21",   vstart[7], step[7], 0, 0);

    // Minimize chi2
    int ierflg = 0;
    gMinuit->Migrad();

    // Get results
    double parVal[NPAR], parErr[NPAR];
    for (int i = 0; i < NPAR; i++) gMinuit->GetParameter(i, parVal[i], parErr[i]);

    std::cout << "===== Simultaneous Fit Results =====" << std::endl;
    std::cout << "mu     = " << parVal[0] << " ± " << parErr[0] << std::endl;
    std::cout << "sigma  = " << parVal[1] << " ± " << parErr[1] << std::endl;
    std::cout << "A1     = " << parVal[2] << " ± " << parErr[2] << std::endl;
    std::cout << "A2     = " << parVal[3] << " ± " << parErr[3] << std::endl;
    std::cout << "b10,b11= " << parVal[4] << "," << parVal[5] << std::endl;
    std::cout << "b20,b21= " << parVal[6] << "," << parVal[7] << std::endl;

    // Create fit functions for plotting
    TF1 *fit1 = new TF1("fit1", "[0] + [1]*x + [2]*TMath::Gaus(x,[3],[4],1)", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
    fit1->SetParameters(parVal[4], parVal[5], parVal[2], parVal[0], parVal[1]);

    TF1 *fit2 = new TF1("fit2", "[0] + [1]*x + [2]*TMath::Gaus(x,[3],[4],1)", h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
    fit2->SetParameters(parVal[6], parVal[7], parVal[3], parVal[0], parVal[1]);

    // Draw histograms and fits
    TCanvas *c = new TCanvas("c", "Simultaneous Fit", 1200, 500);
    c->Divide(2,1);

    c->cd(1);
    h1->SetLineColor(kBlue);
    h1->Draw("E");
    fit1->SetLineColor(kRed);
    fit1->Draw("same");

    c->cd(2);
    h2->SetLineColor(kBlue);
    h2->Draw("E");
    fit2->SetLineColor(kRed);
    fit2->Draw("same");

    c->SaveAs("ex2.pdf");

    return 0;
}

