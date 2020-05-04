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
//#include "DataLoader.h"
//#include "FileMap.h"

//#include "MGTEvent.hh"
//#include "MGTWaveform.hh"

using namespace std;
using namespace gada;

void BiPo_analysis(){
  gROOT->SetBatch();
  const int nChn = 40;
  int firstRun = 53;
  int lastRun = 93;
  char *mapping_file = Form("mapping_phaseII.txt");
  
  //------- DETECTOR MAPPING -------
  struct detector {
    string name;
    int channel;
    int sim;
    string volume;
    int geoID;
    int strin;
    int position;
    string bkgmodel_position;
    int backtobackid;
    int upsidedown;
    string type;
  };
  
  vector<detector> mapping;
  
  std::ifstream file;
  file.open(mapping_file);
  
  while(!file.eof()){
    detector temp;
    file >> temp.name >> temp.channel >> temp.sim >> temp.volume >> temp.geoID >> temp.strin >> temp.position >> temp.bkgmodel_position >> temp.backtobackid >> temp.upsidedown >> temp.type;
    mapping.push_back(temp);
  }
  
  //---------- BUILD MASTERCHAIN ------------
  
  TStopwatch timer_Build_masterChain;
  TStopwatch timer_DoLoops;
  
  std::string dataPath = "/nfs/gerda5/gerda-data/blind/v07.00/gen";
  std::string metadataPath = "/nfs/gerda5/gerda-data/blind/v07.00/meta/data-sets/phy/";
  
  //path cuts.json file /nfs/gerda5/gerda-data/blind/v04.00/meta/config/_aux/
  
  // Initialize map of files
  FileMap myMap;
  // Set root of the data file system
  myMap.SetRootDir(dataPath);
  // Autogenerate list of files
  for (int run = firstRun; run <= lastRun; run++)
    myMap.BuildFromListOfKeys(metadataPath + Form("run00%d-phy-analysis.txt",run));
  
  timer_Build_masterChain.Start();	
  // initialize data loader
  DataLoader l;
  l.AddFileMap(&myMap);
  l.BuildTier1(1,0,0);
  l.BuildTier2(1,0,0);
  l.BuildTier3();
  l.BuildTier4();
  //l.BuildChain("tier2_aux","tier2","tier2_aux");
  
  TChain* masterChain = l.GetUniqueMasterChain();
  
  timer_Build_masterChain.Stop();
  cout << " " << endl;
  cout << " Time for master chain (tier 1+2+3+4) building: RealTime=" <<  timer_Build_masterChain.RealTime() << "s. CpuTime=" <<  timer_Build_masterChain.CpuTime() << "s. Number of events in the Chain=" << masterChain->GetEntries() << endl;
  cout << " " << endl;
  
  //tier1
  MGTEvent *event = 0;
  masterChain -> SetBranchAddress ("tier1_ged.event", &event);
  //tier2
  vector<int>       *GEMDFADC_channelID = new vector<int>;
  vector<int>       *GEMDFADC_eventType = new vector<int>;
  vector<ULong64_t> *GEMDFADC_timestamp = new vector<ULong64_t>;
  vector<UInt_t>    *GEMDFADC_decimalTimestamp = new vector<UInt_t>;
  vector<int>       *GEMDFTTrigger_triggerNumber = new vector<int>;
  vector<int>       *GEMDTrigger_triggerNumber = new vector<int>;
  vector<int>       *GEMDFADC_waveformTag = new vector<int>;
  vector<int>       *GEMDQuality_isProcessed = new vector<int>;
  vector<double> *GEMDFTTrigger_trigger = new vector<double>;
  vector<double> *GEMDFTTrigger_trigger1 = new vector<double>;
  vector<double> *GEMDFTTrigger_trigger2 = new vector<double>;
  vector<double> *GEMDFTTrigger_trigger3 = new vector<double>;
  vector<double> *GEMDFTTrigger_triggerThr = new vector<double>;
  vector<double> *GEMDFTTrigger_baseline = new vector<double>;
  vector<double> *GEMDBaseline_fitExpCoefficient = new vector<double>;
  vector<double> *GEMDEnergyGauss_energy = new vector<double>;
  vector<double> *GEMDEnergyGauss_maxAmpTime  = new vector<double>;
  vector<double> *GEMDEnergyGauss_energyRevPol  = new vector<double>;
  vector<double> *GEMDBaseline_baselineSigma = new vector<double>;
  vector<double> *GEMDBaseline_baseline = new vector<double>;
  vector<double> *GEMDRiseTimeHF_risetime  = new vector<double>;
  vector<double> *GEMDDecayTailFit_exponentialConstant = new vector<double>;
  vector<double> *GEMDDecayTailFit_offset = new vector<double>;
  vector<double> *GEMDDecayTailFit_RMS = new vector<double>;
  masterChain -> SetBranchAddress ("GEMDFADC_channelID",  &GEMDFADC_channelID);
  masterChain -> SetBranchAddress ("GEMDFADC_eventType",  &GEMDFADC_eventType);
  masterChain -> SetBranchAddress ("GEMDFADC_timestamp",  &GEMDFADC_timestamp);
  masterChain -> SetBranchAddress ("GEMDFADC_decimalTimestamp",  &GEMDFADC_decimalTimestamp);
  masterChain -> SetBranchAddress ("GEMDFTTrigger_triggerNumber",  &GEMDFTTrigger_triggerNumber);
  masterChain -> SetBranchAddress ("GEMDFADC_waveformTag",  &GEMDFADC_waveformTag);
  masterChain -> SetBranchAddress ("GEMDQuality_isProcessed",  &GEMDQuality_isProcessed);
  masterChain -> SetBranchAddress ("GEMDFTTrigger_trigger",  &GEMDFTTrigger_trigger);
  masterChain -> SetBranchAddress ("GEMDFTTrigger_trigger1",  &GEMDFTTrigger_trigger1);			     
  masterChain -> SetBranchAddress ("GEMDFTTrigger_trigger2",  &GEMDFTTrigger_trigger2);
  masterChain -> SetBranchAddress ("GEMDFTTrigger_trigger3",  &GEMDFTTrigger_trigger3);
  masterChain -> SetBranchAddress ("GEMDFTTrigger_triggerThr",  &GEMDFTTrigger_triggerThr);
  masterChain -> SetBranchAddress ("GEMDFTTrigger_baseline",  &GEMDFTTrigger_baseline);
  masterChain -> SetBranchAddress ("GEMDBaseline_fitExpCoefficient",  &GEMDBaseline_fitExpCoefficient);
  masterChain -> SetBranchAddress ("GEMDEnergyGauss_energy",  &GEMDEnergyGauss_energy);
  masterChain -> SetBranchAddress ("GEMDEnergyGauss_maxAmpTime",&GEMDEnergyGauss_maxAmpTime);
  masterChain -> SetBranchAddress ("GEMDEnergyGauss_energyRevPol",&GEMDEnergyGauss_energyRevPol);
  masterChain -> SetBranchAddress ("GEMDBaseline_baselineSigma",&GEMDBaseline_baselineSigma);
  masterChain -> SetBranchAddress ("GEMDBaseline_baseline",&GEMDBaseline_baseline);
  masterChain -> SetBranchAddress ("GEMDRiseTimeHF_risetime",&GEMDRiseTimeHF_risetime);
  masterChain -> SetBranchAddress ("GEMDRiseTimeHF_risetime",&GEMDRiseTimeHF_risetime);
  masterChain -> SetBranchAddress ("GEMDDecayTailFit_exponentialConstant",&GEMDDecayTailFit_exponentialConstant);
  masterChain -> SetBranchAddress ("GEMDDecayTailFit_offset",&GEMDDecayTailFit_offset);
  masterChain -> SetBranchAddress ("GEMDDecayTailFit_RMS",&GEMDDecayTailFit_RMS);
  
  //tier3
  vector<int> *failedFlag = 0;
  vector<int> *firedFlag = 0;
  vector<int> *failedFlag_isPhysical = 0;
  vector<int> *isOverflow = 0;
  vector<double> *rawEnergy = new vector<double>;
    
  Int_t isTP;
  Int_t multiplicity;
  Int_t eventChannelNumber;
  Int_t eventNumber;
  Int_t isVetoedInTime;
  ULong64_t timestamp;
  UInt_t decimalTimestamp;
  masterChain -> SetBranchAddress ("firedFlag",  &firedFlag);
  masterChain -> SetBranchAddress ("failedFlag",  &failedFlag);
  masterChain -> SetBranchAddress ("failedFlag_isPhysical",  &failedFlag_isPhysical);
  masterChain -> SetBranchAddress ("isTP",  &isTP);
  masterChain -> SetBranchAddress ("isOverflow",  &isOverflow);
  masterChain -> SetBranchAddress ("multiplicity",  &multiplicity);
  masterChain -> SetBranchAddress ("timestamp",  &timestamp);
  masterChain -> SetBranchAddress ("decimalTimestamp",  &decimalTimestamp);
  masterChain -> SetBranchAddress ("eventChannelNumber",  &eventChannelNumber);
  masterChain -> SetBranchAddress ("eventNumber",  &eventNumber);	 
  masterChain -> SetBranchAddress ("rawEnergy",  &rawEnergy);	 
  //tier4
  Int_t isLArVetoed;
  Int_t isMuVetoed;
  masterChain -> SetBranchAddress ("isLArVetoed",  &isLArVetoed);
  masterChain -> SetBranchAddress ("isMuVetoed",  &isMuVetoed); 
  
  //vector<int>    *psdFlag_ANN_alpha = new vector<int>;   	 
  //vector<bool>   *psdIsEval_ANN_alpha = new vector<bool>;
  //vector<int>    *psdFlag_ANN_mse = new vector<int>;
  //vector<bool>   *psdIsEval_ANN_mse = new vector<bool>;
  vector<int>    *psdFlag_risetime1090 = new vector<int>;
  vector<bool>   *psdIsEval_risetime1090 = new vector<bool>;
  vector<int>    *psdIsEval = new vector<int>;
  vector<int>    *isPSDVetoed = new vector<int>;		 
  vector<double> *AoEclassifier = new vector<double>;	 
  vector<int>    *isAoEvetoed = new vector<int>;
  //masterChain -> SetBranchAddress ("psdFlag_ANN_alpha", &psdFlag_ANN_alpha); 
  //masterChain -> SetBranchAddress ("psdIsEval_ANN_alpha", &psdIsEval_ANN_alpha); 
  //masterChain -> SetBranchAddress ("psdFlag_ANN_mse", &psdFlag_ANN_mse); 
  //masterChain -> SetBranchAddress ("psdIsEval_ANN_mse", &psdIsEval_ANN_mse); 
  masterChain -> SetBranchAddress ("psdFlag_risetime1090", &psdFlag_risetime1090); 
  masterChain -> SetBranchAddress ("psdIsEval_risetime1090", &psdIsEval_risetime1090); 
  masterChain -> SetBranchAddress ("psdIsEval", &psdIsEval); 
  masterChain -> SetBranchAddress ("isPSDVetoed", &isPSDVetoed); 
  masterChain -> SetBranchAddress ("AoEclassifier", &AoEclassifier); 
  masterChain -> SetBranchAddress ("isAoEvetoed", &isAoEvetoed); 
  
  //Codice per vedere le waveform--> scommentare l.BuildTier1(1,0,0);
  //TCanvas  can = (TCanvas) gROOT->GetListOfCanvases()->FindObject("EventDisplay");
  //TCanvas  *can = new TCanvas("EventDisplay","EventDisplay");
  //can->Divide(6,7);
  //if (!can) can = new TCanvas ("EventDisplay", "EventDisplay");
  /*if (!options.Contains("same")) {
      can->Clear();
      can->Divide(4, 4);
  }
  gStyle->SetOptStat(kFALSE);*/

  TFile *out = new TFile("Plot_BiPo.root","RECREATE"); 
  TFile *outWf = new TFile("Waveform.root","RECREATE"); 
  TFile *outTr = new TFile("WaveformTrigger.root","RECREATE"); 
    
  int ch = 0;
  vector<double> DeltaT;
  vector<int> ch_t; 
  vector<unsigned long long> time_selEv;
  time_t tst;
  
  TCanvas *c0 = new TCanvas("FailedFlag","FailedFlag");
  c0->Divide(2,3);
  TH1D *h1_raw= new TH1D("h1_raw"," h1_raw",5,-2.,3.);
  TH1D *h2_raw= new TH1D(" h2_raw"," h2_raw",2000,0.,2000.);
  //TH1D *h2_raw_fF1 = new TH1D(" h2_raw_fF1"," h2_raw_fF1",2000,0.,2000.);
  TH1D *h1= new TH1D("h1"," h1",5,-2.,3.);
  TH1D *h2= new TH1D(" h2"," h2",2000,0.,2000.);
  TH1D *h1_sel= new TH1D("h1_SelEv"," h1_SelEv",5,-2.,3.);
  TH1D *h2_sel= new TH1D(" h2_SelEv"," h2_SelEv",2000,0.,2000.);
  
  vector<int> *fF=0;
  vector<int> *fF_isP=0;

  TH1D *henergy = new TH1D("henergy","Energy event n.1; energy [keV]; counts", 9000, 0 , 9000);
  TH1D *henergy2 = new TH1D("henergy2","Energy event n.2; energy [keV]; counts", 9000, 0 , 9000);
  TH2D *henergy_corr = new TH2D("henergy_corr","Energy of events; energy event n.1 [keV]; energy event n.2 [keV]", 9000, 0 , 9000, 9000, 0 , 9000);
  
  TCanvas *c1 = new TCanvas("TimeStamp","TimeStamp");
  c1->Divide(1,2);
  TH1D * h_t = new TH1D("TimeStamp_PhaseII","TimeStamp_PhaseII",100,1450000000,1520000000);
  TH1D * h_SelEv_t = new TH1D("SelEv_t","SelEv_t",100,1450000000,1520000000);
  
  TCanvas *c1_bis = new TCanvas("deltaT","deltaT");
  c1_bis->Divide(1,3);
  
  TH1D * h_BaselineSigma = new TH1D("h_BaselineSigma","h_BaselineSigma", 300, 0, 10000);
  TH1D * h_BaselineSigmaTr = new TH1D("h_BaselineSigmaTr","h_BaselineSigmaTr", 300, 0, 10000);
  h_BaselineSigma->SetTitle("Baseline Sigma");
  h_BaselineSigmaTr->SetTitle("Baseline Sigma - Trigger == 2");
  
  TCanvas *c11 = new TCanvas("BaselineSigma","BaselineSigma");
  c11->Divide(2,1);
  
  TH1D * h_Delta_t_sameCh= new TH1D("h_delta_t_sameCh","h_delta_t_sameCh", 300, 0, 2000);
  h_Delta_t_sameCh->SetTitle("Time Diff - Events in Same Detector");
  h_Delta_t_sameCh->GetXaxis()->SetTitle("time [#mus]");
  h_Delta_t_sameCh->GetYaxis()->SetTitle("counts");
  
  TH1D * h_Delta_t_otherCh= new TH1D("h_delta_t_otherCh","h_delta_t_otherCh", 300, 0, 2000);
  h_Delta_t_otherCh->SetTitle("Time Diff - Events in Different Detectors");
  h_Delta_t_otherCh->GetXaxis()->SetTitle("time [#mus]");
  h_Delta_t_otherCh->GetYaxis()->SetTitle("counts");
  
  TH1D * h_Delta_t= new TH1D("delta_t","delta_t", 300, 0, 2000); //
  h_Delta_t->SetTitle("Time Diff - All Events");
  h_Delta_t->GetXaxis()->SetTitle("time [#mus]");
  h_Delta_t->GetYaxis()->SetTitle("counts");
  
  TH1D * h_Delta_t_nTr= new TH1D("delta_t_nTr","delta_t_nTr", 300, 0, 100); 
  h_Delta_t_nTr->GetXaxis()->SetTitle("time [#mus]");
  h_Delta_t_nTr->GetYaxis()->SetTitle("counts");
    
  TCanvas *c2 = new TCanvas("trigger","trigger");
  c2->Divide(3,2);
  
  TH1D * h_triggerNo = new TH1D("TriggerNumber","TriggerNumber",25,0.,25.);
  TH1D * h_trigger = new TH1D("Trigger","Trigger",100,0.,160000.);
  TH1D * h_trigger1 = new TH1D("Trigger1","Trigger1",100,0.,160000.);
  TH1D * h_trigger2 = new TH1D("Trigger2","Trigger2",100,0.,160000.);
  TH1D * h_DeltaTrigger = new TH1D("DeltaTrigger","DeltaTrigger",100,0.,1000.);
  TH1D * h_DeltaTrigger21 = new TH1D("DeltaTrigger21","DeltaTrigger21",100,0.,100000.);
  
  TCanvas *c21 = new TCanvas("tailfit","tailfit");
  c21->Divide(3,2);
  
  TH1D *h_TailFit_exp = new TH1D("h_TailFit_exp","TailFit Exponential Constant",1000,1.,1000000.);
  TH1D *h_TailFit_off = new TH1D("h_TailFit_off","TailFit Offset",1000,0.1,10.);
  TH1D *h_TailFit_RMS = new TH1D("h_TailFit_RMS","TailFit RMS",1000,0.1,10.);
  TH1D *h_TailFit_expTr = new TH1D("h_TailFit_expTr","TailFit Exponential Constant - Trigger = 2",1000,1.,1000000.);
  TH1D *h_TailFit_offTr = new TH1D("h_TailFit_offTr","TailFit Offset - Trigger = 2",1000,0.1,10.);
  TH1D *h_TailFit_RMSTr = new TH1D("h_TailFit_RMSTr","TailFit RMS - Trigger = 2",1000,0.1,10.);
  
  TCanvas *c22 = new TCanvas("TriggTailFit","TriggTailFit");
  c22->Divide(3,2);
  
  TH2D *htt_corr1 = new TH2D("htt_corr1","Correlation DeltaTrigger vs TailFit Exp Const; Delta Trigger; Tail Fit Exp Const", 100,0,100000,1000,1.,1000000);
  TH2D *htt_corr2 = new TH2D("htt_corr2","Correlation DeltaTrigger vs TailFit Offset; Delta Trigger; Tail Fit Offset", 100,0.,100000,1000,0.1,10);
  TH2D *htt_corr3 = new TH2D("htt_corr3","Correlation DeltaTrigger vs TailFit RMS; Delta Trigger; Tail Fit Offset", 100,0.,100000,1000,0.1,10);
  TH2D *htt_corr4 = new TH2D("htt_corr4","Correlation Trigger2 vs TailFit Exp Const; Trigger2; Tail Fit Exp Const", 100,0,160000.,1000,1,1000000);
  TH2D *htt_corr5 = new TH2D("htt_corr5","Correlation Trigger2 vs TailFit Offset; Trigger2; Tail Fit Offset",100,0,160000,1000,0.1,10);
  TH2D *htt_corr6 = new TH2D("htt_corr6","Correlation Trigger2 vs TailFit RMS; Trigger2; Tail Fit Offset",100,0,160000,1000,0.1,10);
  htt_corr1->SetMarkerStyle(7);
  htt_corr1->SetMarkerColor(kBlue+2);
  htt_corr2->SetMarkerStyle(7);
  htt_corr2->SetMarkerColor(kBlue+2);
  htt_corr3->SetMarkerStyle(7);
  htt_corr3->SetMarkerColor(kBlue+2);
  htt_corr4->SetMarkerStyle(7);
  htt_corr4->SetMarkerColor(kBlue+2);
  htt_corr5->SetMarkerStyle(7);
  htt_corr5->SetMarkerColor(kBlue+2);
  htt_corr6->SetMarkerStyle(7);
  htt_corr6->SetMarkerColor(kBlue+2);
  
  TCanvas *c23 = new TCanvas("CorrEneTail","CorrEneTail");
  c23->Divide(2,2);
  
  TH2D *henet_corr1 = new TH2D("henet_corr1","Correlation Energy vs TailFit Exp Const; Energy; Tail Fit Exp Const", 100,0,9000.,1000,1,1000000);
  TH2D *henet_corr2 = new TH2D("henet_corr2","Correlation Energy vs TailFit Offset; Energy; Tail Fit Offset",9000,1,9000,1000,0.1,10);
  TH2D *henet_corr3 = new TH2D("henet_corr3","Correlation Energy vs TailFit RMS; Energy; Tail Fit Offset",9000,1,9000,1000,0.1,10);
  TH1D *henet = new TH1D("henet","Energy events; energy [keV]; counts", 9000, 0 , 9000);
  
  TCanvas *c3 = new TCanvas("FailedFlag2","FailedFlag2");
  c3->Divide(4,2);
  
  TH1D * h_fF_nTr2             = new TH1D("h_fF_nTr2","h_fF_nTr2",5,-2.,3.);
  TH1D * h_fF_isP_nTr2         = new TH1D("h_fF_isP_nTr2","h_fF_isP_nTr2",2000,0.,2000.);
  TH1D * h_mulitplicity_nTr2   = new TH1D("h_mulitplicity_nTr2","h_mulitplicity_nTr2",nChn,0,nChn-1);
  TH1D * h_isTP_nTr2           = new TH1D("h_isTP_nTr2","h_isTP_nTr2",4,-1,3);
  TH1D * h_fF_NOfitExpCoef     = new TH1D("h_fF_NOfitExpCoef","h_fF_NOfitExpCoef",5,-2.,3.);
  TH1D * h_fF_isP_NOfitExpCoef = new TH1D("h_fF_isP_NOfitExpCoef","h_fF_isP_NOfitExpCoef",2000,0.,2000.);
  TH1D * h_mulitplicity_NOfitExpCoef   = new TH1D("h_mulitplicity_NOfitExpCoef","h_mulitplicity_NOfitExpCoef",nChn,0,nChn-1);
  TH1D * h_isTP_NOfitExpCoef           = new TH1D("h_isTP_NOfitExpCoef","h_isTP_NOfitExpCoef",4,-1,3);
  
  TCanvas *c4 = new TCanvas("energy","energy");
  c4->Divide(1,3);
  
  ofstream file_typeI("BiPo-events-typeI.txt");
  ofstream file_typeI_diff("BiPo-events-typeI-diff.txt");
  ofstream file_typeII("BiPo-events-typeII.txt");
  
  timer_DoLoops.Start();
  
  long double time_ref = 0.;
  long double time_ref2 = 0.;
  long double time_j = 0.;
  long double time_j2 = 0.;
  long double firstTimestamp = 0.;
  int cont = 0; 
  int cont_coppiaAll = 0; 
  int cont_coppiaChn = 0; 
  int cont_nTr2 = 0;
  bool first = true;
  
  //loop sulle entries della masterChain		
  //for (int i = 0; i < masterChain->GetEntries(); i++) {
  
  for (int i = 0; i < masterChain->GetEntries(); i++) {
    masterChain->GetEntry(i);
    
    if( i % 1000 == 0 ) cout << " ...processing event " << i << " of " <<  masterChain->GetEntries() << endl;
    h_t->Fill( timestamp + 1e-9*decimalTimestamp );
    
    //START loop for cross-check on failedFlags (to be removed)
    for ( int ch = 0; ch < eventChannelNumber; ch++ ) {  
      h1_raw->Fill(failedFlag->at(ch));
      h2_raw->Fill(failedFlag_isPhysical->at(ch));
      //if(failedFlag->at(ch)==1) h2_raw_fF1->Fill(failedFlag_isPhysical->at(ch));
      
      h_triggerNo->Fill(GEMDFTTrigger_triggerNumber->at(ch));
      h_trigger->Fill(GEMDFTTrigger_trigger->at(ch));
      h_trigger1->Fill(GEMDFTTrigger_trigger1->at(ch));
      h_trigger2->Fill(GEMDFTTrigger_trigger2->at(ch));
      
      h_TailFit_exp->Fill(GEMDDecayTailFit_exponentialConstant->at(ch));
      h_TailFit_off->Fill(GEMDDecayTailFit_offset->at(ch));
      h_TailFit_RMS->Fill(GEMDDecayTailFit_RMS->at(ch));
      
      h_BaselineSigma->Fill(GEMDBaseline_baselineSigma->at(ch));
      
      if( GEMDFTTrigger_triggerNumber->at(ch) == 2 ) {
	//h_DeltaTrigger->Fill( GEMDFTTrigger_trigger1->at(ch)-GEMDFTTrigger_trigger->at(ch) ); 
	h_DeltaTrigger21->Fill( GEMDFTTrigger_trigger2->at(ch)-GEMDFTTrigger_trigger1->at(ch) ); 
	h_fF_nTr2->Fill(failedFlag->at(ch)); 
	h_fF_isP_nTr2->Fill(failedFlag_isPhysical->at(ch)); 
	h_mulitplicity_nTr2->Fill(multiplicity);
	h_isTP_nTr2->Fill(isTP);
		
      }
      
      if( TMath::Abs( GEMDBaseline_fitExpCoefficient->at(ch) )>1000 ) {
	h_fF_NOfitExpCoef->Fill(failedFlag->at(ch)); 
	h_fF_isP_NOfitExpCoef->Fill(failedFlag_isPhysical->at(ch)); 
	h_mulitplicity_NOfitExpCoef->Fill(multiplicity);
	h_isTP_NOfitExpCoef->Fill(isTP);
      }
    }  //END loop for cross-check on failedFlags (to be removed)
    
    int cont_ch = 0;
    int cont_ch_nTr = 0;
    bool nTr = false;
    int isFiredDet1[nChn] = {0};
    int isFiredDet2[nChn] = {0};
    int firedDet1 = -99;
    int firedDet2 = -99;
    int firedDet_nTr = -99;
    double ene1[nChn] = {0.};
    double ene2[nChn] = {0.};
    int multi1[nChn] = {0};
    int multi2[nChn] = {0};
    
    TH1D *hwf1[nChn];
    
    if( isMuVetoed == 0 && isTP == 0 ){// && isLArVetoed==0 ){
      //cout << "GOOD event: " << i << " MuVeto=" <<  isMuVetoed << " Multiplicity= " << multiplicity << " isTP=" << isTP << " LarVetoed=" << isLArVetoed << endl; 
      int LArVeto1 = isLArVetoed;
      //int psdFlag_ANN_alpha1 = 0;       
      //bool psdIsEval_ANN_alpha1 = 0;
      //int psdFlag_ANN_mse1 = 0;	 
      //bool psdIsEval_ANN_mse1 = 0;	 
      int psdFlag_risetime10901= 0;	 
      bool psdIsEval_risetime10901 = 0; 
      int psdIsEval1 = 0;		       
      int isPSDVetoed1 = 0;		 
      double AoEclassifier1 = 0;	       
      int isAoEvetoed1 = 0;                   
      
      for (int ch = 0; ch < eventChannelNumber; ch++) {
	h1->Fill(failedFlag->at(ch));
	h2->Fill(failedFlag_isPhysical->at(ch));
	
	//cout << "	channel= " << ch  <<  " failedFlag= " << failedFlag->at(ch) << " failedFlag_isPhysical= " << failedFlag_isPhysical->at(ch) << " firedFlag= " << firedFlag->at(ch) << " rawEnergy = " << rawEnergy->at(ch) << " GEMDEnergyGauss_energy= " << GEMDEnergyGauss_energy->at(ch)  << endl;
	
	//------ BiPo Events of Type 2 -----
	
	if( failedFlag_isPhysical->at(ch) == 32 && GEMDFTTrigger_triggerNumber->at(ch) == 2 && GEMDFTTrigger_trigger2->at(ch) > 82000 && rawEnergy->at(ch) > 200 ) {
	  
	  
	  MGTWaveform *wftr = event->GetWaveformID(ch);
	  TH1D *hwftr = wftr->GimmeUniqueHist();
	  hwftr->SetTitle( Form("failedFlag = %d",failedFlag_isPhysical->at(ch) ) );
	  double trigg1 = GEMDFTTrigger_trigger1->at(ch);
	  double trigg2 = GEMDFTTrigger_trigger2->at(ch);
	  double tailFitExp = GEMDDecayTailFit_exponentialConstant->at(ch);
	  double tailFitRMS = GEMDDecayTailFit_RMS->at(ch);
	  

	  TH1D *hwftr_d = new TH1D("hwftr_d","hwftr_d", 4095, 40, 163840);
	  TH1D *hwfSm = new TH1D(Form("hs%d",i),Form("hs%d",i), 4095, 40, 163840);
	  
	  for( int k = 1; k < 4096; k++ )
	    hwftr_d->SetBinContent(k+1, hwftr->GetBinContent(k)-hwftr->GetBinContent(k+1));
	  
	  for( int k = 1; k < 4096; k++ ){
	    double contSum = 0.;
	    for (int l = 0; l < 20; l ++ )
	      contSum += hwftr_d->GetBinContent(k+l);
	    hwfSm->SetBinContent(k, contSum);
	  }
	  
	  TLegend *leg = new TLegend(0.7, 0.7, 0.99, 0.99);
	  leg->AddEntry( (TObject*)0, Form("Trigger1 = %g ns", GEMDFTTrigger_trigger1->at(ch)),"");
	  leg->AddEntry( (TObject*)0, Form("Trigger2 = %g ns", GEMDFTTrigger_trigger2->at(ch)),"");
	  leg->AddEntry( (TObject*)0, Form("TailFit-exp = %g ns", GEMDDecayTailFit_exponentialConstant->at(ch)), "");
	  leg->AddEntry( (TObject*)0, Form("TailFit-RMS = %g", GEMDDecayTailFit_RMS->at(ch)),"");
	  leg->AddEntry( (TObject*)0, Form("BaselineSigma = %g", GEMDBaseline_baselineSigma->at(ch)),"");
	  
	  TCanvas *canTr = new TCanvas(Form("wf%d",i),Form("wf%d",i));
	  canTr->Divide(2,1);
	  canTr->cd(1);
	  gStyle->SetOptStat(0);
	  hwftr->DrawCopy();
	  leg->Draw();
	  canTr->cd(2);
	  gStyle->SetOptStat(0);
	  //hwftr_d->Draw();
	  hwfSm->GetXaxis()->SetRangeUser(84000,160384);
	  hwfSm->Draw();
	  
	  TSpectrum *s = new TSpectrum(10);
	  s->Search(hwfSm, 2.5, "", 0.30);
	  //cout << "Found peaks: " << nfound << endl;
	  
	  int nfound = s->GetNPeaks();
	  double *xpeaksb = s->GetPositionX();
	  TH1 *hb = s->Background(hwftr_d, 20, "same");
	  
	  if ( nfound == 1  && GEMDDecayTailFit_exponentialConstant->at(ch) > 0 ){
	    cout << "BiPo Event of type II in ch " << ch << endl;
	    cout << "event " << i << " energy " << rawEnergy->at(ch) << " BaselineSigma " << GEMDBaseline_baselineSigma->at(ch) << " Trigger2Position " << GEMDFTTrigger_trigger2->at(ch) << " PeakPosition: " << xpeaksb[0] << " " << xpeaksb[1] << endl;
	    
	    file_typeII << i << " " << timestamp << " " << ch << " " << rawEnergy->at(ch) << " " << GEMDFTTrigger_trigger->at(ch) << " " << GEMDFTTrigger_trigger2->at(ch) << " " << " " << xpeaksb[0] << endl;
	    
	    h_BaselineSigmaTr->Fill(GEMDBaseline_baselineSigma->at(ch));
	    h_TailFit_expTr->Fill(GEMDDecayTailFit_exponentialConstant->at(ch));
	    h_TailFit_offTr->Fill(GEMDDecayTailFit_offset->at(ch));
	    h_TailFit_RMSTr->Fill(GEMDDecayTailFit_RMS->at(ch));
	    
	    htt_corr1->Fill(GEMDFTTrigger_trigger2->at(ch)-GEMDFTTrigger_trigger1->at(ch), GEMDDecayTailFit_exponentialConstant->at(ch));
	    htt_corr2->Fill(GEMDFTTrigger_trigger2->at(ch)-GEMDFTTrigger_trigger1->at(ch), GEMDDecayTailFit_offset->at(ch));
	    htt_corr3->Fill(GEMDFTTrigger_trigger2->at(ch)-GEMDFTTrigger_trigger1->at(ch), GEMDDecayTailFit_RMS->at(ch));
	    htt_corr4->Fill(GEMDFTTrigger_trigger2->at(ch), GEMDDecayTailFit_exponentialConstant->at(ch));
	    htt_corr5->Fill(GEMDFTTrigger_trigger2->at(ch), GEMDDecayTailFit_offset->at(ch));
	    htt_corr6->Fill(GEMDFTTrigger_trigger2->at(ch), GEMDDecayTailFit_RMS->at(ch));
	    henet->Fill(rawEnergy->at(ch));
	    henet_corr1->Fill(rawEnergy->at(ch), GEMDDecayTailFit_exponentialConstant->at(ch));
	    henet_corr2->Fill(rawEnergy->at(ch), GEMDDecayTailFit_offset->at(ch));
	    henet_corr3->Fill(rawEnergy->at(ch), GEMDDecayTailFit_RMS->at(ch));
	    
	    canTr->Update();
	    outTr->cd();
	    canTr->Write();
	    hwfSm->Write();
	  }
	  delete hwfSm;
	  delete hwftr_d;
	  delete canTr;
	}
	
	//------ END: BiPo Events of Type 2 -----
	
	if( rawEnergy->at(ch) > 0. && rawEnergy->at(ch) < 10000 && failedFlag->at(ch) == 0 ){
	  ene1[ch] = rawEnergy->at(ch);
	  multi1[ch] = multiplicity;
	  firedDet1 = ch;
	  isFiredDet1[ch] = 1;
	  cont_ch++;
	  
	  //psdFlag_ANN_alpha1 =   	  psdFlag_ANN_alpha->at(ch);   	 
	  //psdIsEval_ANN_alpha1 =	  psdIsEval_ANN_alpha->at(ch);	 
	  //psdFlag_ANN_mse1 =	  psdFlag_ANN_mse->at(ch);	 
	  //psdIsEval_ANN_mse1 =	  psdIsEval_ANN_mse->at(ch);	 
	  psdFlag_risetime10901=	  psdFlag_risetime1090->at(ch);	 
	  psdIsEval_risetime10901 =  psdIsEval_risetime1090->at(ch); 
	  psdIsEval1 =		  psdIsEval->at(ch);		 
	  isPSDVetoed1 =		  isPSDVetoed->at(ch);		 
	  AoEclassifier1 =	  AoEclassifier->at(ch);		 
	  isAoEvetoed1 =              isAoEvetoed->at(ch);            	  
	}
	MGTWaveform *wf1 = event->GetWaveformID(ch);
	hwf1[ch] = wf1->GimmeUniqueHist();
	//hwf1->SetTitle( Form("Flag = %d",failedFlag_isPhysical->at(ch)));
      }
      //check on the multplicity
      if(cont_ch == 0) {
	continue;
      }
      
      //h1_sel->Fill(failedFlag->at(firedDet));
      //h2_sel->Fill(failedFlag_isPhysical->at(firedDet));
      
      time_ref=0;
      time_ref2=0;
      
      //ch_t.push_back(firedDet);
      time_selEv.push_back(timestamp+1e-9*decimalTimestamp);		
      time_ref= (long double) timestamp;
      time_ref2= (long double) 1.0e-9*decimalTimestamp;
      
      h_SelEv_t->Fill(timestamp+1e-9*decimalTimestamp);
      
      //----- Star loop on 2nd event -----
      
      for (int j = i+1; j < masterChain->GetEntries(); j++) {
	//cout << "i " << i << " j " << j << endl;
	masterChain->GetEntry(j);
	
	time_j= (long double) timestamp;
	time_j2= (long double) 1.0e-9*decimalTimestamp;
	long double timeDiff = time_j+time_j2-(time_ref+time_ref2);
	
	if( timeDiff < 0.00164 ) { //cut on timeDiff between 2 events, deltaT < 10*164um
	  if( isMuVetoed == 0  && isTP == 0){// && isLArVetoed==0 ){
	    cout << i << " timeDiff = " << timeDiff << endl;
	    cont_coppiaAll++;
	    DeltaT.push_back(timeDiff);
	    
	    TCanvas *can1 = new TCanvas(Form("wf%d",i),Form("wf%d",i));
	    can1->Divide(6,7);
	    
	    TCanvas *can2 = new TCanvas(Form("wf%d",j),Form("wf%d",j));
	    can2->Divide(6,7);
	    
	    //TCanvas  *can = new TCanvas(Form("wf%d",j),Form("wf%d",j));
	    //can->Divide(6,7);
	    
	    int LArVeto2 = isLArVetoed;	    
	    int cont_ch2 = 0;
	    bool coppia = false;
	    for (int ch = 0; ch < eventChannelNumber; ch++){
	      if( rawEnergy->at(ch) > 300. && rawEnergy->at(ch) < 10000 && ( failedFlag->at(ch) == 0 || failedFlag_isPhysical->at(ch) == 8 )){
		ene2[ch] = rawEnergy->at(ch);
		multi2[ch] = multiplicity;

		h_Delta_t->Fill( timeDiff*1e6 );
		cont_ch2++;
		firedDet2 = ch;
		isFiredDet2[ch] = 1;
		
		coppia = true;
	      }
	      
	      //----- Plot of the waveforms -----
	      
	      can1->cd(ch+1);
	      hwf1[ch]->DrawCopy();
	      //hwf1->Write();
	      can1->Update();
	      
	      MGTWaveform* wf2 = event->GetWaveformID(ch);
	      TH1D *hwf2 = wf2->GimmeUniqueHist();
	      hwf2->SetTitle( Form("Flag = %d", failedFlag_isPhysical->at(ch)));
	      can2->cd(ch+1);
	      hwf2->DrawCopy();
	      //hwf2->Write();
	      can2->Update();
	    }
	    
	    //---- Save Plot ----
	    if (coppia){
	      //plotWf(masterChain, outWf, i, j);
	      outWf->cd();
	      can1->Write();
	      can2->Write();
	    }
	    
	    //for (int ch = 0; ch < eventChannelNumber; ch++){
	    if ( firedDet1 == firedDet2 ){//BiPo Events type1 in the same detector
	      //if ( firedDet[ch] == 1 && firedDet2[ch] == 1 ){
	      int ch = firedDet1;
	      cont_coppiaChn++;
	      h_Delta_t_sameCh->Fill( timeDiff*1e6 );
	      henergy->Fill(ene1[ch]);
	      henergy2->Fill(ene2[ch]);
	      henergy_corr->Fill(ene1[ch],ene2[ch]);
	      
	      cout << "BiPo Event of type I in ch " << ch << endl;
	      cout << "multiplicity = " << multi1[ch] << " " << multi2[ch] << " energy = " << ene1[ch] << " " << ene2[ch] << ", failedFlag2 = " << failedFlag_isPhysical->at(ch) << " isLarVetoed1 = " << LArVeto1 << " isLarVetoed2 = " << LArVeto2 << " timeDiff = " << timeDiff << endl;
	      
	      file_typeI << i << " " << timestamp <<  " " << ch << " " << " " << multi1[ch] << " " << multi2[ch] << " " << ene1[ch] << " " << ene2[ch] << " " << failedFlag_isPhysical->at(ch) << " " << LArVeto1 << " " << LArVeto2 << " " << psdIsEval1 << " " << psdIsEval->at(ch) << " " << isPSDVetoed1 << " " << isPSDVetoed->at(ch) << " " << timeDiff << endl;
	      
	      cout << "PSD flags event " << i << " and " << j << ":" << endl;
	      //cout <<  psdFlag_ANN_alpha1    << " " << psdFlag_ANN_alpha->at(ch) << endl;   	 
	      //cout <<  psdIsEval_ANN_alpha1  << " " << 	  psdIsEval_ANN_alpha->at(ch) << endl;	 
	      //cout <<  psdFlag_ANN_mse1      << " " << 	  psdFlag_ANN_mse->at(ch) << endl;	 
	      //cout <<  psdIsEval_ANN_mse1    << " " << 	  psdIsEval_ANN_mse->at(ch) << endl;	 
	      cout <<  psdFlag_risetime10901 << " " << 	  psdFlag_risetime1090->at(ch) << endl;	 
	      cout <<  psdIsEval_risetime10901 << " " <<   psdIsEval_risetime1090->at(ch) << endl; 
	      cout <<  psdIsEval1 << " " << 		  psdIsEval->at(ch) << endl;		 
	      cout <<  isPSDVetoed1 << " " << 		  isPSDVetoed->at(ch) << endl;		 
	      cout <<  AoEclassifier1 << " " << 	  AoEclassifier->at(ch) << endl;		 
	      cout <<  isAoEvetoed1 << " " <<               isAoEvetoed->at(ch) << endl;            
	    }
	    else {//Event in different detectors 
	      h_Delta_t_otherCh->Fill( timeDiff*1e6 );
	      
	      cout << firedDet1 << " " << firedDet2 << endl;
	      
	      //----- Mapping check ------
	      int stringDiff = fabs( mapping[firedDet1].strin - mapping[firedDet2].strin );
	      int posDiff = fabs( mapping[firedDet1].position - mapping[firedDet2].position );
	      if ( stringDiff == 0 ){
		if ( posDiff == 1 ){
		  cout << "Events in same string and near detectors" << endl; 
		  cout << firedDet1 << " " << firedDet2 << " multiplicity = " << multi1[firedDet1] << " " << multi2[firedDet2] << " energy = " << ene1[firedDet1] << " " << ene2[firedDet2] << ", failedFlag2 = " << failedFlag_isPhysical->at(firedDet2) << " isLarVetoed1 = " << LArVeto1 << " isLarVetoed2 = " << LArVeto2 << " timeDiff = " << timeDiff << endl;
		  
		  file_typeI_diff << i << " " << timestamp << " " << firedDet1 << " " << firedDet2 << " " << multi1[firedDet1] << " " << multi2[firedDet2] << " " << ene1[firedDet1] << " " << ene2[firedDet2] << " 0 " << failedFlag_isPhysical->at(firedDet2) << " " << LArVeto1 << " " << LArVeto2 << " " << timeDiff << endl;
	     	}
	      }
	      else if ( stringDiff == 1 || mapping[firedDet1].strin == 7 || mapping[firedDet2].strin == 7 || ( mapping[firedDet1].strin == 1 && mapping[firedDet2].strin == 6 ) || ( mapping[firedDet1].strin == 6 && mapping[firedDet2].strin == 1 ) ) {
		if ( mapping[firedDet1].bkgmodel_position == mapping[firedDet2].bkgmodel_position ) {
		  cout << "Events in strings " << mapping[firedDet1].strin << " and " << mapping[firedDet2].strin << " with near detectors" << endl; 
		  cout << firedDet1 << " " << firedDet2 << " multiplicity = " << multi1[firedDet1] << " " << multi2[firedDet2] << " energy = " << ene1[firedDet1] << " " << ene2[firedDet2] << ", failedFlag2 = " << failedFlag_isPhysical->at(firedDet2) << " isLarVetoed1 = " << LArVeto1 << " isLarVetoed2 = " << LArVeto2 << " timeDiff = " << timeDiff << endl;
		  
		  file_typeI_diff << i << " " << timestamp << " " << firedDet1 << " " << firedDet2 << " " << multi1[firedDet1] << " " << multi2[firedDet2] << " " << ene1[firedDet1] << " " << ene2[firedDet2] << " 0 " << failedFlag_isPhysical->at(firedDet2) << " " << LArVeto1 << " " << LArVeto2 << " " << timeDiff << endl;
		}
	      }
	      else {
		cout << "Events in chn " << firedDet1 << " and " << firedDet2 << ", in strings " << mapping[firedDet1].strin << " and " << mapping[firedDet2].strin << " -> NOT a BiPo EVENT" << endl;
	      }
	    }
	    
	    //if (cont_ch2 != 1) break;
	    delete can1;
	    delete can2;
	  }//first check on isMuVetoed and isTP
	}//check on time diff
	else break;
      }//second loop on entries
    } //first check on isMuVetoed and isTP
  } //first loop on entries
  timer_DoLoops.Stop();
  file_typeI.close();
  file_typeI_diff.close();
  file_typeII.close();
  
  cout << " "  << endl;
  cout << " Time for master Chain entries loops: real time = " << timer_DoLoops.RealTime() << "s. CpuTime =  " << timer_DoLoops.CpuTime() << "s" << endl;
  cout << " "  << endl;
  cout << " Triggered Events                               " << masterChain->GetEntries() <<  endl;
  cout << " Events whithin 10*164 us                       " << cont_coppiaAll << endl;
  cout << " Events whithin 10*164 us in the same detector  " << cont_coppiaChn << endl;
  //cout << " Pile-up with nTrigger = 2                      " << cont_nTr2 << endl; 
  //cout << " Total no. of BiPo-like events                  " << cont_coppia + cont_nTr2 << endl;
  cout << " "  << endl;
  
  out->cd();
  c0->cd(1); gPad->SetLogy(); h1_raw->Draw("hist");
  c0->cd(2); gPad->SetLogy(); h2_raw->Draw("hist");
  c0->cd(3); gPad->SetLogy(); h1->Draw("hist");
  c0->cd(4); gPad->SetLogy(); h2->Draw("hist");
  c0->cd(5); gPad->SetLogy(); h1_sel->Draw("hist");
  c0->cd(6); gPad->SetLogy(); h2_sel->Draw("hist");
  c0->Update();
  c0->Write();
  
  h_t->GetXaxis()->SetTimeOffset(0,"gmt");
  h_t->GetXaxis()->SetTimeDisplay(1);	
  h_t->GetXaxis()->SetTimeFormat("%d-%m-%Y");
  h_t->GetXaxis()->SetNdivisions(-510);
  
  h_SelEv_t->GetXaxis()->SetTimeOffset(0,"gmt");
  h_SelEv_t->GetXaxis()->SetTimeDisplay(1);	
  h_SelEv_t->GetXaxis()->SetTimeFormat("%d-%m-%Y");
  h_SelEv_t->GetXaxis()->SetNdivisions(-510);
  
  c1->cd(1); h_t->Draw();	
  c1->cd(2); h_SelEv_t->Draw();
  c1->Update();
  c1->Write();
  
  c1_bis->cd(1); h_Delta_t_sameCh->Draw();
  c1_bis->cd(2); h_Delta_t_otherCh->Draw();
  c1_bis->cd(3); h_Delta_t->Draw();
  c1_bis->Update();
  c1_bis->Write();
  c1_bis->Print("Plots/deltaT.pdf");
  
  c11->cd(1); gPad->SetLogy(); h_BaselineSigma->Draw();
  c11->cd(2); gPad->SetLogy(); h_BaselineSigmaTr->Draw();
  c11->Update();
  c11->Write();
  
  c2->cd(1); gPad->SetLogy(); h_triggerNo->Draw("hist");
  c2->cd(2); h_DeltaTrigger21->Draw("hist");
  c2->cd(3); gPad->SetLogy(); h_trigger->Draw("hist");
  c2->cd(4); gPad->SetLogy(); h_trigger1->Draw("hist");
  c2->cd(5); gPad->SetLogy(); h_trigger2->Draw("hist");
  c2->cd(6); h_DeltaTrigger->Draw("hist");
  c2->Update();
  c2->Write();
  c2->Print("Plots/trigger.pdf");
  
  c21->cd(1); gPad->SetLogy(); h_TailFit_exp->Draw("hist");
  c21->cd(2); gPad->SetLogy(); h_TailFit_off->Draw("hist");
  c21->cd(3); gPad->SetLogy(); h_TailFit_RMS->Draw("hist");
  c21->cd(4); gPad->SetLogy(); h_TailFit_expTr->Draw("hist");
  c21->cd(5); gPad->SetLogy(); h_TailFit_offTr->Draw("hist");
  c21->cd(6); gPad->SetLogy(); h_TailFit_RMSTr->Draw("hist");
  c21->Update();
  c21->Write();
  c21->Print("Plots/tailfit.pdf");
    
  c22->cd(1); htt_corr1->Draw("COLZ");
  c22->cd(2); htt_corr2->Draw("COLZ");
  c22->cd(3); htt_corr3->Draw("COLZ");
  c22->cd(4); htt_corr4->Draw("COLZ");
  c22->cd(5); htt_corr5->Draw("COLZ");
  c22->cd(6); htt_corr6->Draw("COLZ");
  c22->Update();
  c22->Write();
  c22->Print("Plots/trigger-tailfit.pdf");
  
  c23->cd(1); henet->Draw();
  c23->cd(2); henet_corr1->Draw("COLZ");
  c23->cd(3); henet_corr2->Draw("COLZ");
  c23->cd(4); henet_corr3->Draw("COLZ");
  c23->Update();
  c23->Write();
  c23->Print("Plots/energy-tailfit.pdf");
    
  c3->cd(1);  gPad->SetLogy(); h_fF_nTr2->Draw();
  c3->cd(2);  gPad->SetLogy(); h_fF_isP_nTr2->Draw();
  c3->cd(3);  gPad->SetLogy(); h_mulitplicity_nTr2->Draw();
  c3->cd(4);  gPad->SetLogy(); h_isTP_nTr2->Draw();
  c3->cd(5);  gPad->SetLogy(); h_fF_NOfitExpCoef->Draw();
  c3->cd(6);  gPad->SetLogy(); h_fF_isP_NOfitExpCoef->Draw();
  c3->cd(7);  gPad->SetLogy(); h_mulitplicity_NOfitExpCoef->Draw();
  c3->cd(8);  gPad->SetLogy(); h_isTP_NOfitExpCoef->Draw();
  c3->Update();
  c3->Write();
  
  c4->cd(1);  henergy->Draw();
  c4->cd(2);  henergy2->Draw();
  c4->cd(3);  henergy_corr->Draw();
  c4->Update();
  c4->Write();
  c4->Print("Plots/energy.pdf");
  
  out->Close();	
  outWf->Close();	
  outTr->Close();	
  
  delete GEMDFADC_channelID;
  delete GEMDFADC_eventType;
  delete GEMDFADC_timestamp;
  delete GEMDFADC_decimalTimestamp;
  delete GEMDFTTrigger_triggerNumber;
  delete GEMDFADC_waveformTag;
  delete GEMDQuality_isProcessed;
  delete GEMDFTTrigger_trigger;
  delete GEMDFTTrigger_trigger1;
  delete GEMDFTTrigger_trigger2;
  delete GEMDFTTrigger_trigger3;
  delete GEMDBaseline_fitExpCoefficient;
  delete GEMDEnergyGauss_energy;
  delete GEMDEnergyGauss_maxAmpTime;
  delete GEMDEnergyGauss_energyRevPol;
  delete GEMDBaseline_baselineSigma;
  delete GEMDBaseline_baseline;
  delete GEMDRiseTimeHF_risetime;

  delete out;
  return;
}


