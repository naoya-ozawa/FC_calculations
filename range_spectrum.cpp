#include <iostream>
#include <bits/stdc++.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TRint.h>
using namespace std;

int main(int argc, char** argv){
  TRint rootapp("app",&argc,argv);
  TCanvas *c1 = new TCanvas();
//	c1->Divide(2,1);

  fstream file;
  string word, filename;

  filename = "./RANGE_3D_Fr_in_SUS303_5keV.txt";
//  filename = "./RANGE_3D_Fr_in_SUS303_1keV.txt";
//  filename = "./RANGE_3D_Fr_in_SUS303_50eV.txt";
//  filename = "./RANGE_3D_Fr_in_SUS303_10eV.txt";
//  filename = "./RANGE_3D_Fr_in_SUS303_500eV.txt";
//  filename = "./RANGE_3D_Fr_in_SUS303_100eV.txt";
//  filename = "./EXYZ_alpha_in_SUS303_5MeV.txt";
//  filename = "./EXYZ_alpha_in_SUS303_6MeV.txt";
//  filename = "./EXYZ_alpha_in_SUS303_7MeV.txt";

  int energy = 5000;

  file.open(filename.c_str());

  TH1D *spectrum = new TH1D("spectrum",Form("Stopping Range of %d eV Fr ions Injected into SUS303",energy),5000,0.,50.);
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

  cout << Nion << " ions analyzed" << endl;

  spectrum->Draw();

  mean /= double(Nion);
  cout << "Mean Range = " << mean << " AA" << endl;
  stdev /= double(Nion);
  stdev = TMath::Sqrt(stdev - (mean*mean));
  cout << "StDev = " << stdev << " AA" << endl;

  TF1 *f_gaus = new TF1("f_gaus","gaus(0)",0.,40.);
  f_gaus->SetParameter(0,1.);
  f_gaus->SetParameter(1,mean);
  f_gaus->FixParameter(1,mean);
  f_gaus->SetParameter(2,stdev);
  f_gaus->FixParameter(2,stdev);
  spectrum->Fit("f_gaus","B");

  TLatex *l = new TLatex();
  l->SetTextSize(0.05);
  l->DrawLatex(3,10.,Form("Stopping Range = %5.4f #pm %5.4f #AA",mean,stdev));

  rootapp.Run();

  c1->Update();
  c1->Modified();

  return 0;

}
