#include <bits/stdc++.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>
#include <TRint.h>
using namespace std;

int main(int argc, char** argv){
	TRint rootapp("app",&argc,argv);
	TCanvas *c1 = new TCanvas();
	c1->Divide(2,2);

	// Part 1: formulation
	c1->cd(1);
	TLatex *l_form = new TLatex();
	l_form->DrawLatex(0.1,0.9,"Energy loss model:");
	l_form->DrawLatex(0.1,0.8,"E_{loss}(d,E_{#alpha}) #approx #left[ #frac{#partial E_{loss}}{#partial d} (E_{#alpha}) #right] (d + d_{0})");
	l_form->DrawLatex(0.1,0.65,"#frac{#partial E_{loss}}{#partial d} (E_{#alpha}) #approx #frac{C_{0}}{E_{#alpha}} + C_{1}");
	l_form->DrawLatex(0.2,0.5,"E_{loss}: energy loss in the SUS303 [eV]");
	l_form->DrawLatex(0.2,0.4,"d: depth [#AA]");
	l_form->DrawLatex(0.2,0.3,"E_{#alpha}: #alpha ray energy [MeV]");

	// Part 2: depth dependence
	c1->cd(2);
	double depth[4] = {5.,10.,15.,20.}; // AA
	double depthe[4] = {0.,0.,0.,0.};
	double eloss_5MeV[4] = {251.442,417.111,593.732,775.864}; // eV (mean)
	double elose_5MeV[4] = {134.677,182.910,246.030,311.039}; // eV (stdev)
	TGraphErrors *eloss5 = new TGraphErrors(4,depth,eloss_5MeV,depthe,elose_5MeV);
	eloss5->SetTitle("E_{#alpha} = 5 MeV");
	eloss5->SetLineColor(2);
	eloss5->SetMarkerColor(2);
	eloss5->Draw();
	TF1 *feloss5 = new TF1("feloss5","[0]*(x+[1])",4.,21.);
	feloss5->SetParameters(1.,1.);
	feloss5->SetLineColor(2);
	feloss5->SetTitle("Linear fit (5 MeV)");
	eloss5->Fit("feloss5","B");
	double dedd5 = feloss5->GetParameter(0);
	double dedd5e = feloss5->GetParError(0);
	double d05 = feloss5->GetParameter(1);
	double d05e = feloss5->GetParError(1);

	double eloss_6MeV[4] = {225.178,372.781,535.411,688.817}; // eV (mean)
	double elose_6MeV[4] = {120.941,165.764,228.367,283.768}; // eV (stdev)
	TGraphErrors *eloss6 = new TGraphErrors(4,depth,eloss_6MeV,depthe,elose_6MeV);
	eloss6->SetTitle("E_{#alpha} = 6 MeV");
	eloss6->SetLineColor(3);
	eloss6->SetMarkerColor(3);
	eloss6->Draw();
	TF1 *feloss6 = new TF1("feloss6","[0]*(x+[1])",4.,21.);
	feloss6->SetParameters(1.,1.);
	feloss6->SetLineColor(3);
	feloss6->SetTitle("Linear fit (6 MeV)");
	eloss6->Fit("feloss6","B");
	double dedd6 = feloss6->GetParameter(0);
	double dedd6e = feloss6->GetParError(0);
	double d06 = feloss6->GetParameter(1);
	double d06e = feloss6->GetParError(1);

	double eloss_7MeV[4] = {204.638,335.169,481.610,625.967}; // eV (mean)
	double elose_7MeV[4] = {109.773,151.335,203.452,262.193}; // eV (stdev)
	TGraphErrors *eloss7 = new TGraphErrors(4,depth,eloss_7MeV,depthe,elose_7MeV);
	eloss7->SetTitle("E_{#alpha} = 7 MeV");
	eloss7->SetLineColor(4);
	eloss7->SetMarkerColor(4);
	eloss7->Draw();
	TF1 *feloss7 = new TF1("feloss7","[0]*(x+[1])",4.,21.);
	feloss7->SetParameters(0,1.);
	feloss7->SetLineColor(4);
	feloss7->SetTitle("Linear fit (7 MeV)");
	eloss7->Fit("feloss7","B");
	double dedd7 = feloss7->GetParameter(0);
	double dedd7e = feloss7->GetParError(0);
	double d07 = feloss7->GetParameter(1);
	double d07e = feloss7->GetParError(1);

	TMultiGraph *mg_eloss = new TMultiGraph();
	mg_eloss->SetTitle("Depth dependence of energy loss of #alpha particles in SUS303;Depth d [#AA];Energy loss E_{loss} [eV]");
	mg_eloss->Add(eloss5);
	mg_eloss->Add(eloss6);
	mg_eloss->Add(eloss7);
	mg_eloss->Draw("AP*");
	feloss5->Draw("SAME");
	feloss6->Draw("SAME");
	feloss7->Draw("SAME");
	c1->cd(2)->BuildLegend();
	double dedd = (dedd5 + dedd6 + dedd7) / 3.0;
	double dedde = TMath::Sqrt(dedd5e*dedd5e + dedd6e*dedd6e + dedd7e*dedd7e) / 3.0;
	double d0 = (d05 + d06 + d07) / 3.0;
	double d0e = TMath::Sqrt(d05e*d05e + d06e*d06e + d07e*d07e) / 3.0;

	// Part 3: energy dependence
	c1->cd(3);
	double ealpha[3] = {5.,6.,7.}; // MeV
	double ealphae[3] = {0.,0.,0.}; // MeV
	double elossae[3] = {dedd5,dedd6,dedd7};
	double elossaee[3] = {dedd5e,dedd6e,dedd7e};
	TGraphErrors *elossedep = new TGraphErrors(3,ealpha,elossae,ealphae,elossaee);
	elossedep->SetTitle("E_{#alpha} dependence of (#partialE_{loss}/#partiald);#alpha particle energy E_{#alpha} [MeV];#frac{#partial E_{loss}}{#partial d} [eV/#AA]");
	elossedep->Draw("AP*");
	TF1 *felossedep = new TF1("felossedep","[0]/x + [1]");
	felossedep->SetParameters(1.,1.);
	elossedep->Fit("felossedep");
	double deddcoef = felossedep->GetParameter(0);
	double deddcoefe = felossedep->GetParError(0);
	double deddconst = felossedep->GetParameter(1);
	double deddconste = felossedep->GetParError(1);

	
	// Part 4: result
	c1->cd(4);
	TLatex *l_eloss = new TLatex();
	l_eloss->DrawLatex(0.1,0.9,Form("d_{0}|_{E_{#alpha} = 5 MeV} = %3.3f #pm %3.3f #AA",d05,d05e));
	l_eloss->DrawLatex(0.1,0.8,Form("d_{0}|_{E_{#alpha} = 6 MeV} = %3.3f #pm %3.3f #AA",d06,d06e));
	l_eloss->DrawLatex(0.1,0.7,Form("d_{0}|_{E_{#alpha} = 7 MeV} = %3.3f #pm %3.3f #AA",d07,d07e));
	l_eloss->DrawLatex(0.2,0.6,Form("#rightarrow (AVERAGE) d_{0} = %3.3f #pm %3.3f #AA",d0,d0e));
	l_eloss->DrawLatex(0.1,0.5,Form("C_{0} = %3.3f #pm %3.3f eV MeV/#AA",deddcoef,deddcoefe));
	l_eloss->DrawLatex(0.1,0.4,Form("C_{1} = %3.3f #pm %3.3f eV/#AA",deddconst,deddconste));
	l_eloss->DrawLatex(0.1,0.3,"Ignoring fit parameter errors,");
	l_eloss->DrawLatex(0.1,0.2,Form("E_{loss} = #left[ #frac{%3.1f}{E_{#alpha}} + %3.1f #right] (d + %3.1f)",deddcoef,deddconst,d0));
	l_eloss->DrawLatex(0.4,0.1,Form("#pm #sqrt{ #left[ #frac{%3.1f}{E_{#alpha}^{2}} + %3.1f #right](#Deltad)^{2} + #frac{%3.1f}{E_{#alpha}^{4}}(d^{2} + %3.1f)(#DeltaE_{#alpha})^{2} }",deddcoef*deddcoef,deddconst*deddconst,deddcoef*deddcoef,d0*d0));


	rootapp.Run();
	c1->Update();
	c1->Modified();

	return 0;
}
