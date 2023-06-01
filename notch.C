
#include "TF1.h"
#include "TMath.h"

Double_t fitfunc(Double_t *x, Double_t *par) {
  Double_t xx = x[0];
  Double_t z = 2 *TMath::Pi() * xx;
  Double_t A = 1 - z * z * par[0] * par[1];
  Double_t a = 1 - z * z * par[2] * par[3];
  Double_t B = z * par[0] * par[4];
  Double_t b = z * par[2] * par[5];
  Double_t r = par[6] + 50; 
  Double_t C = z * par[1] * a + par[4] * b + r * B * a + r * b * A + z * par[3] * A + par[5] * B;
  Double_t D = par[4] * a - z * par[1] * b + r * A * a - r * b * B + par[5] * A - z * par[3] * B;
  Double_t E = A * a - B * b;
  Double_t F = A * b + B * a;
  Double_t val =  5 * r * TMath::Sqrt((E * E + F * F) / (C * C + D * D));
  return val;
}

Double_t fitphase(Double_t *x, Double_t *par) {
  Double_t xx = x[0];
  Double_t z = 2 * TMath::Pi() * xx;
  Double_t A = 1 - z * z * par[0] * par[1];
  Double_t a = 1 - z * z * par[2] * par[3];
  Double_t B = z * par[0] * par[4];
  Double_t b = z * par[2] * par[5];
  Double_t r = par[6];
  Double_t C = z * par[1] * a + par[4] * b + r * B * a + r * b * A +
               z * par[3] * A + par[5] * B;
  Double_t D = par[4] * a - z * par[1] * b + r * A * a - r * b * B +
               par[5] * A - z * par[3] * B;
  Double_t E = A * a - B * b;
  Double_t F = A * b + B * a;
  Double_t val = 30 + 4.5e-04 * xx +  180 * TMath::ATan((-C * E + D * F) / (D * E + F * C)) / TMath::Pi();
  return val;
}

  void notch() {
    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    auto multi = new TMultiGraph("mg", "mg");

    c1->Divide(1);
    c2->Divide(3);

    TGraph *graphbf4 = new TGraph("bf4.dat", "%lg %lg ");
    TGraph *graphbf5 = new TGraph("bf5.dat", "%lg %lg ");
    TGraph *graphbf6 = new TGraph("bf6.dat", "%lg %lg ");
    TH1F *ondaqua = new TH1F("h1", "Errore onda quadra", 30, 2.53, 2.54);
    TGraph *fase2194 = new TGraph("fase4(bella).dat", "%lg %lg");
    TGraph *phaserr = new TGraph("phaserr.dat", "%lg %lg");
    TGraph *fase10k = new TGraph("fase5(10k).dat", "%lg %lg");
    TGraph *fase680 = new TGraph("fase6(680).dat", "%lg %lg");
    TF1 *fit1 = new TF1("linear", "fitfunc", 300, 3000, 7);
    TF1 *fit2 = new TF1("linear", "fitfunc", 300, 3000, 7);
    TF1 *fit3 = new TF1("linear", "fitfunc", 300, 3000, 7);
    TF1 *fit4 = new TF1("linear", "fitfunc", 300, 14000, 7);
    TF1 *fit5 = new TF1("linear", "fitfunc", 300, 14000, 7);
    TF1 *fit6 = new TF1("linear", "fitfunc", 300, 13000, 7);
    TF1 *fn1 = new TF1("fi", "[A]+[B]*x", 1100.0, 10000.0);
    TF1 *line = new TF1("linear", " 30 + 6.96467e-05*x", 300.0, 14000.0);
    TF1 *phasefit10 = new TF1("linear", " fitphase", 400.0, 14000.0, 7);
    TF1 *phasefit22 = new TF1("linear", " fitphase", 400.0, 14000.0, 7);
    TF1 *phasefit68 = new TF1("linear", " fitphase", 400.0, 14000.0, 7);

    multi->Add(graphbf4);
    multi->Add(graphbf5);
    multi->Add(graphbf6);

    fit1->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 2194 );
    fit2->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 680 );
    fit3->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 9974 );
    fit4->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 2194 );
    fit5->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 680 );
    fit6->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 9974 );
    phasefit10->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 9974);
    phasefit22->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 2194);
    phasefit68->SetParameters(0.000000158, 0.04786, 0.0000004657, 0.000479, 128., 2.25, 680);
    ifstream in;
    in.open("ondaq.dat");
    float x;
    float y;
    while (1) {
      in >> x >> y;
      if (!in.good()) {
        break;
      }
      ondaqua->Fill(y);
    }
    in.close();
    c1->cd(1);
    graphbf4->GetXaxis()->SetRangeUser(0, 14000);
    graphbf4->GetYaxis()->SetRangeUser(1., 5.5);

    graphbf4->SetMarkerStyle(8);
    graphbf4->SetLineWidth(4);
    graphbf4->SetLineColor(kBlue);
    graphbf4->SetMarkerSize(0.6);
    graphbf4->SetMarkerColor(kBlue);

    graphbf5->SetMarkerStyle(8);
    graphbf5->SetLineWidth(4);
    graphbf5->SetLineColor(kViolet);
    graphbf5->SetMarkerSize(0.6);
    graphbf5->SetMarkerColor(kViolet);

    graphbf6->SetMarkerStyle(8);
    graphbf6->SetLineWidth(4);
    graphbf6->SetLineColor(kOrange);
    graphbf6->SetMarkerSize(0.6);
    graphbf6->SetMarkerColor(kOrange);

    fit1->SetLineColor(kCyan);
    fit4->SetLineColor(kCyan);
    fit2->SetLineColor(41);
    fit5->SetLineColor(41);
    fit3->SetLineColor(8);
    fit6->SetLineColor(8);
    fit1->SetLineWidth(4);
    fit2->SetLineWidth(4);
    fit3->SetLineWidth(4);
    fit4->SetLineWidth(4);
    fit5->SetLineWidth(4);
    fit6->SetLineWidth(4);

   // graphbf4->Fit(fit1, "R");
    graphbf4->Fit(fit4, "R");
   // graphbf5->Fit(fit2, "R");
    graphbf5->Fit(fit5, "R");
   // graphbf6->Fit(fit3, "R");
    graphbf6->Fit(fit6, "R");


    multi->GetXaxis()->SetLabelSize(0.03);
    multi->GetXaxis()->SetNdivisions(28, 10 , 0 , kTRUE);
    gPad->SetGrid();
    c1->BuildLegend();
    gStyle->SetOptFit(1100);
    multi->Draw("ALP");

    c2->cd(1);
    gPad->SetGrid();
    fase680->Draw("AL");
    fase680->Fit(phasefit68);
    //phasefit68->Draw("SAME");
    c2->cd(2);
    gPad->SetGrid();
    fase2194->Draw("AL");
    fase2194->Fit(phasefit22);
    //phasefit22->Draw("SAME");
    c2->cd(3);
    gPad->SetGrid();
    fase10k->Draw("AL");
    //line->Draw("SAME"); 
    fase10k->Fit(phasefit10);
    //phasefit10->Draw("SAME");
  }