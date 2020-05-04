#include "TMath.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TSpectrum.h"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

void fitWaveformCurr(){
  
  TFile *file  = TFile::Open("WaveformTrigger.root");
  int event[11] = {1470807, 1470825, 1471161, 1475775, 1479273,1479586,1479973,1483363,1485864,1486305,1511133};
  
  //TCanvas *canTr = new TCanvas("canTr","canTr");
  //canTr->Divide(3,3);
  
  for (int i = 0; i < 10; i++){
    TH1D *hwfSm = (TH1D*)file->Get(Form("hs%d",event[i]));
    
    TCanvas *canTr = new TCanvas(Form("canTr%d",i),Form("canTr%d",i));    
    //canTr->Divide(2,1);
    canTr->cd();
    //gStyle->SetOptStat(0);
    hwfSm->GetXaxis()->SetRangeUser(84000,160384);
    hwfSm->Draw();
    canTr->Update();
    
    TSpectrum *s1 = new TSpectrum(10);
    s1->Search(hwfSm, 2.5, "", 0.30);
    
    int nfound1 = s1->GetNPeaks();
    double *xpeaksb1 = s1->GetPositionX();
    TH1 *hb1 = s1->Background(hwfSm, 20, "same");
    
    cout << "PeakFound = " << nfound1 << " PeakPosition = " << xpeaksb1[0] << " " << xpeaksb1[1] << endl;
    
  }
}
