#include <iostream>
#include <TRandom3.h>
#include <bits/stdc++.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <THStack.h>
#include <TF1.h>
#include <TLatex.h>
#include <TRint.h>
using namespace std;

double eloss(double E, double d){
  double d_0 = 2.236; // AA
  double C_0 = 118.663; // eV MeV / AA
  double C_1 = 10.966; // eV / AA
  return ( (C_0/E) + C_1 ) * (d + d_0); // eV
}

double eloss_stdv(double E, double d){
  double d_0p = 22.907; // AA
  double C_2 = 0.247; // eV / AA^2
  double C_3 = 0.052; // eV / AA^2
  return ( C_2 - C_3*TMath::Log(E) ) * (d + d_0p) * (d + d_0p); // eV
}

int main(int argc, char** argv){
  TRint rootapp("app",&argc,argv);
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);

  // RANGE file
  c1->cd(1);

  fstream file;
  string word, filename;

  filename = "./RANGE_3D_Fr_in_SUS303_5keV.txt";
  int V_A = 5000;

  file.open(filename.c_str());

  TH1D *spectrum = new TH1D("spectrum",Form("Stopping Range of %d eV Fr ions Injected into SUS303",V_A),7000,0.,70.);
  spectrum->GetXaxis()->SetTitle("Range (#AA)");
  spectrum->GetYaxis()->SetTitle("Ions / 1 pm");

  int Nion = 0;

  int dataflag = 0;
  int datacol = 0;
  int newion = 0; // 0: record data, 1: has been recorded

  double prevDepth = 0.0;

  double mean = 0.0;
  double stdev = 0.0;

  while (file >> word){
//		cout << word << endl;
    try{
      if (word == "0000001") dataflag = 1;
      if (dataflag == 1){
	++datacol;
//        cout << "datacol = " << datacol << ", ";
        if (datacol%4 == 1){
	  Nion = stoi(word);
	}else if (datacol%4 == 2){
	  double range = stod(word); // AA
	  spectrum->Fill(range);
	  mean += range;
	  stdev += range*range;
	}else if (datacol%4 == 3){
	  double latY = stod(word); // AA
        }else{
	  double latZ = stod(word); // AA
	}
      }
    }
    catch (const invalid_argument& e){
      cout << "EOF" << endl;
    }
  }

  spectrum->Draw();



  c1->cd(2);

  int Npar = 1000000;
  double E_alpha = 6.545; // MeV; 210-Fr main branch
  double E_alpha_stdv = 0.005; // MeV; 210-Fr main branch
  double SSD_stdv = 0.02; // MeV; ORTEC SSD

  double range_ll = 6.4;
  double range_ul = 6.7;

  int color_before = 38;
  int color_after = 46;
  TRandom3 *rnd = new TRandom3();

  // Default energy spectrum: Fr piled up on FC surface
  TH1D *def_spect = new TH1D("def_spect","Adsorbed on surface",int((range_ul-range_ll)*1000),range_ll,range_ul);
  for (int i=0; i<Npar; ++i){
    def_spect->Fill(rnd->Gaus(E_alpha, TMath::Sqrt(E_alpha_stdv*E_alpha_stdv + SSD_stdv*SSD_stdv)));
  }
  def_spect->SetLineColor(color_before);

  // Degraded energy spectrum; Fr penetrated into FC
  TH1D *deg_spect = new TH1D("deg_spect",Form("V_{A} = %d V extraction",V_A),int((range_ul-range_ll)*1000),range_ll,range_ul);
  for (int i=0; i<Npar; ++i){
    double orig = rnd->Gaus(E_alpha, TMath::Sqrt(E_alpha_stdv*E_alpha_stdv + SSD_stdv*SSD_stdv));
    double degr = spectrum->GetRandom();
    double energy = orig - rnd->Gaus(TMath::Power(10.,-6)*eloss(orig,degr), TMath::Power(10.,-6)*eloss_stdv(orig,degr));
    deg_spect->Fill(energy);
  }
  deg_spect->SetLineColor(color_after);

  THStack *energy_spectrum = new THStack();
  energy_spectrum->SetTitle(Form("Energy Spectrum of 6.545 MeV #alpha ray from {}^{210}Fr;Energy (MeV);Counts (/%3.3f MeV)",def_spect->GetXaxis()->GetBinWidth(0)));
  energy_spectrum->Add(deg_spect);
  energy_spectrum->Add(def_spect);
  energy_spectrum->Draw("nostack");

  c1->cd(2)->BuildLegend();

  rootapp.Run();

  c1->Update();
  c1->Modified();

  return 0;

}
