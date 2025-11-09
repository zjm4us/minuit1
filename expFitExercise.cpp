#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TROOT.h>

using namespace std;

// --- Global pointers for TMinuit (needed for FCN) ---
TH1* gHist = nullptr;

// Exponential fit function for TMinuit
void expFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) {
    fval = 0.0;
    for (int i = 1; i <= gHist->GetNbinsX(); ++i) {
        double x = gHist->GetBinCenter(i);
        double y = gHist->GetBinContent(i);
        double yfit = par[0] * exp(par[1] * x);
        double err = (y > 0) ? sqrt(y) : 1.0; // simple Poisson errors
        fval += ((y - yfit) / err) * ((y - yfit) / err);
    }
}

int main() {
    gROOT->SetBatch(); // run without GUI windows

    // --- Open ROOT file ---
    TFile *file = TFile::Open("/home/zjm4us/compphys/minuit1/distros.root");
    if (!file || file->IsZombie()) {
        cerr << "Error opening ROOT file!" << endl;
        return 1;
    }
    cout << "Successfully opened ROOT file: " << file->GetName() << endl;
    file->ls();

    // --- Part A: Fit dist1 ---
    TH1F *dist1 = (TH1F*)file->Get("dist1");
    if (!dist1) {
        cerr << "dist1 histogram not found!" << endl;
        return 1;
    }

    gHist = dist1;
    TMinuit minuitA(2);
    minuitA.SetFCN(expFCN);
    minuitA.DefineParameter(0, "p0", 1.0, 0.1, 0, 0);
    minuitA.DefineParameter(1, "p1", -0.1, 0.01, 0, 0);
    minuitA.Migrad();

    double p0, p1, err0, err1;
    minuitA.GetParameter(0, p0, err0);
    minuitA.GetParameter(1, p1, err1);
    cout << "Fit results for dist1:\n";
    cout << "p0 = " << p0 << " ± " << err0 << "\n";
    cout << "p1 = " << p1 << " ± " << err1 << "\n";

    // Draw and save Part A plot
    TCanvas *c1 = new TCanvas("c1", "Fit dist1", 800, 600);
    dist1->Draw();
    TF1 *fit1 = new TF1("fit1", "[0]*exp([1]*x)", dist1->GetXaxis()->GetXmin(), dist1->GetXaxis()->GetXmax());
    fit1->SetParameters(p0, p1);
    dist1->Fit(fit1, "R");
    c1->SaveAs("dist1_fit.pdf");

    // --- Part B: Fit dist2 projections ---
    TH2F *dist2 = (TH2F*)file->Get("dist2");
    if (!dist2) {
        cerr << "dist2 histogram not found!" << endl;
        return 1;
    }

    // Projection X
    TH1D *projX_d = dist2->ProjectionX("projX");
    TH1F *projX = (TH1F*)projX_d; // cast is okay for plotting
    gHist = projX;
    TMinuit minuitX(2);
    minuitX.SetFCN(expFCN);
    minuitX.DefineParameter(0, "p0", 1.0, 0.1, 0, 0);
    minuitX.DefineParameter(1, "p1", -0.1, 0.01, 0, 0);
    minuitX.Migrad();
    minuitX.GetParameter(0, p0, err0);
    minuitX.GetParameter(1, p1, err1);
    cout << "Fit results for dist2 Projection X:\n";
    cout << "p0 = " << p0 << " ± " << err0 << "\n";
    cout << "p1 = " << p1 << " ± " << err1 << "\n";

    TCanvas *c2 = new TCanvas("c2", "Projection X", 800, 600);
    projX->Draw();
    TF1 *fitX = new TF1("fitX", "[0]*exp([1]*x)", projX->GetXaxis()->GetXmin(), projX->GetXaxis()->GetXmax());
    fitX->SetParameters(p0, p1);
    projX->Fit(fitX, "R");
    c2->SaveAs("dist2_projX_fit.pdf");

    // Projection Y
    TH1D *projY_d = dist2->ProjectionY("projY");
    TH1F *projY = (TH1F*)projY_d;
    gHist = projY;
    TMinuit minuitY(2);
    minuitY.SetFCN(expFCN);
    minuitY.DefineParameter(0, "p0", 1.0, 0.1, 0, 0);
    minuitY.DefineParameter(1, "p1", 0.05, 0.001, 0, 0);
    minuitY.Migrad();
    minuitY.GetParameter(0, p0, err0);
    minuitY.GetParameter(1, p1, err1);
    cout << "Fit results for dist2 Projection Y:\n";
    cout << "p0 = " << p0 << " ± " << err0 << "\n";
    cout << "p1 = " << p1 << " ± " << err1 << "\n";

    TCanvas *c3 = new TCanvas("c3", "Projection Y", 800, 600);
    projY->Draw();
    TF1 *fitY = new TF1("fitY", "[0]*exp([1]*x)", projY->GetXaxis()->GetXmin(), projY->GetXaxis()->GetXmax());
    fitY->SetParameters(p0, p1);
    projY->Fit(fitY, "R");
    c3->SaveAs("dist2_projY_fit.pdf");

    cout << "\nAll plots saved. Combine them into ex1.pdf for submission." << endl;

    file->Close();
    return 0;
}

