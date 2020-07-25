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
	
	fstream file;
	string word, filename;

	filename = "./TRANSMIT_alpha_in_SUS303_5MeV_5A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_5MeV_10A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_5MeV_15A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_5MeV_20A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_6MeV_5A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_6MeV_10A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_6MeV_15A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_6MeV_20A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_7MeV_5A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_7MeV_10A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_7MeV_15A.txt";
//	filename = "./TRANSMIT_alpha_in_SUS303_7MeV_20A.txt";

	file.open(filename.c_str());

	double Einit = 5.; // MeV
	int Tsus = 5; // AA

	TH1D *spectrum = new TH1D("spectrum",Form("Energy loss of %3.f MeV #alpha transmitted through %d #AA SUS303",Einit,Tsus),100,0.,1000.);
	spectrum->GetXaxis()->SetTitle("Energy (eV)");
	spectrum->GetYaxis()->SetTitle("Ions / 10 eV");

	int dataflag = 0;
	int datacol = 0;

	double mean = 0.0;
	double stdev = 0.0;

	while (file >> word){
		try{
			if (word == "T") dataflag = 1;
			if (dataflag == 1){
				++datacol;
				if (datacol%10 == 1){
					string head = word;
					if (head != "T") dataflag = 0;
				}else if (datacol%10 == 2){
					int Nion = stoi(word);
				}else if (datacol%10 == 3){
					int A = stoi(word);
				}else if (datacol%10 == 4){
					// cout << word << endl;
					double eloss = Einit*TMath::Power(10.,6) - stod(word);
					spectrum->Fill(eloss);
					// cout << eloss << endl;
					mean += eloss;
					stdev += eloss*eloss;
				}else if (datacol%10 == 5){
					double depth = stod(word);
				}else if (datacol%10 == 6){
					double LposY = stod(word);
				}else if (datacol%10 == 7){
					double LposZ = stod(word);
				}else if (datacol%10 == 8){
					double cosx = stod(word);
				}else if (datacol%10 == 9){
					double cosy = stod(word);
				}else{
					double cosz = stod(word);
				}
			}
		}
		catch (const invalid_argument& e){
			cout << "EOF" << endl;
		}
	}
	int Nions = datacol/10;
	cout << Nions << " ions analyzed" << endl;
	mean /= double(Nions);
	stdev /= double(Nions);
	stdev = TMath::Sqrt(stdev - (mean*mean));
	cout << "Energy loss: " << mean << " +- " << stdev << " eV" << endl;

	spectrum->Draw();
	TF1 *normdist = new TF1("normdist","gaus(0)",0.,2000.);
	normdist->SetParameter(0,1.); // coef
	normdist->SetParameter(1,mean); // mean
	normdist->FixParameter(1,mean); // mean
	normdist->SetParameter(2,stdev); // stdev
	normdist->FixParameter(2,stdev); // stdev
	spectrum->Fit("normdist","B");
	TLatex *l_spect = new TLatex();
	l_spect->SetTextAlign(12);
	l_spect->SetTextSize(0.05);
	l_spect->DrawLatex(0.,100.,Form("E_{loss} = %3.3f #pm %3.3f eV",mean,stdev));
	
	// The emittance can be calculated using cosx = v_x/v_tot, etc.
	c1->Update();
	c1->Modified();
	rootapp.Run();
	return 0;
}
