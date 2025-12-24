
#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"

#include "StMaker.h"
#include "StChain.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"

#include "macros/loadSharedHFLibraries.C"

#include <iostream>
#include <ctime>
#include <cstdio>

#include "StPicoLambdaAnaMaker/StPicoLambdaAnaMaker.h"

using namespace std;

#else
class StChain;
#endif

StChain *chain;

//triggerSetup = 0 - MB, 1 - JP2
void runPicoLambdaAnaMaker(const Char_t *inputFile="test.list", const Char_t *outputFile="outputBaseName",  
    const unsigned int makerMode = 0 /*kAnalyze*/,
    const Char_t *badRunListFileName = "picoList_bad_MB.list", const Char_t *treeName = "picoHFtree",
    const Char_t *productionBasePath = "/star/data100/reco/AuAu_200_production_2016/ReversedFullField/P16ij/2016",
    const unsigned int decayChannel = 0 /* kChannel0 */, string SL_version = "SL20c",
    const int triggerSetup = 0)
{ 
  // -- Check STAR Library. Please set SL_version to the original star library used in the production 
  //    from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
  //string SL_version = "SL20c"; //also SL21d for Run17 is available with smaller production sample 08/21/2022
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
    cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
    exit(1);
  }

#ifdef __CINT__
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif


  chain = new StChain();

  // ========================================================================================
  //makerMode    = StPicoHFMaker::kAnalyze;
  // ========================================================================================

  cout << "Maker Mode    " << makerMode << endl;
  cout << "TreeName      " << treeName << endl; 
  cout << "Decay Channel " << decayChannel << endl; 

  TString sInputFile(inputFile);
  TString sInputListHF("");  
  TString sProductionBasePath(productionBasePath);
  TString sTreeName(treeName);

  if (makerMode == StPicoHFMaker::kAnalyze || makerMode == StPicoHFMaker::kQA)
  {
    if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root"))
    {
      cout << "No input list or picoDst root file provided! Exiting..." << endl;
      exit(1);
    }
  }
  else if (makerMode == StPicoHFMaker::kWrite)
  {
    if (!sInputFile.Contains("picoDst.root"))
    {
      cout << "No input picoDst root file provided! Exiting..." << endl;
      exit(1);
    }
  }
  else if (makerMode == StPicoHFMaker::kRead)
  {
    if (!sInputFile.Contains(".list"))
    {
      cout << "No input list provided! Exiting..." << endl;
      exit(1);
    }

    // -- prepare filelist for picoDst from trees
    sInputListHF = sInputFile;
    sInputFile = "tmpPico.list";

    TString command = "sed 's|" + sTreeName + ".root|picoDst.root|g' " + sInputListHF + " > " + sInputFile;
    cout << "COMMAND : " << command << endl; 
    gSystem->Exec(command.Data());

    command = "sed -i 's|^.*" + sTreeName + "|" + sProductionBasePath + "|g' " + sInputFile; // + " > " + sInputFile;
    cout << "COMMAND : " << command << endl; 
    gSystem->Exec(command.Data());
  }
  else
  {
    cout << "Unknown makerMode! Exiting..." << endl;
    exit(1);
  }

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, sInputFile, "picoDstMaker"); //for local testing only

  StPicoLambdaAnaMaker* picoLambdaAnaMaker = new StPicoLambdaAnaMaker("picoLambdaAnaMaker", picoDstMaker, outputFile, sInputListHF);
  picoLambdaAnaMaker->setMakerMode(makerMode);
  picoLambdaAnaMaker->setDecayChannel(StPicoLambdaAnaMaker::kChannel1);//not needed?
  picoLambdaAnaMaker->setTreeName(treeName);


  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
  picoLambdaAnaMaker->setHFBaseCuts(hfCuts);

  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis

  // -- File name of bad run list
  hfCuts->setBadRunListFileName(badRunListFileName);

  // -- ADD USER CUTS HERE ----------------------------
  hfCuts->setCutVzMax(60.); //original 30
  hfCuts->setCutVzVpdVzMax(1e6); //original 3 - now open

  //MB
  if(triggerSetup == 0)
  {
    //Run12 pp200GeV triggers
    hfCuts->addTriggerId(370001); //VPDMB
    hfCuts->addTriggerId(370011); //VPDMB-nosmd
    
    
    //Run15 pp200GeV triggers
     
    //st_ssdmb triggers
    hfCuts->addTriggerId(470011); //VPDMB-5-ssd
    hfCuts->addTriggerId(470021); //VPDMB-5-ssd
    hfCuts->addTriggerId(480001); //VPDMB-5-ssd    
    hfCuts->addTriggerId(490001); //VPDMB-5-ssd
    hfCuts->addTriggerId(500001); //VPDMB-5-ssd
    hfCuts->addTriggerId(510009); //VPDMB-5-ssd
    

    //Run17 pp510GeV triggers
    hfCuts->addTriggerId(570001);    // VPDMB-30 (st_physics)  
  }
  //JP
  else if(triggerSetup == 1)
  {
    //Run12 pp200GeV triggers
    //hfCuts->addTriggerId(370611); // JP1
    hfCuts->addTriggerId(370621); // JP2
    
    //Run12 pp200GeV triggers we don't want
    hfCuts->addBadTriggerId(370001); //VPDMB
    hfCuts->addBadTriggerId(370011); //VPDMB-nosmd
    

    //Run17 pp510GeV triggers
    hfCuts->addTriggerId(570401); // JP2  

    //Run17 pp510GeV triggers we don't want
    hfCuts->addBadTriggerId(570001); // VPDMB-30 (st_physics)
  
  }
  else
  {
    cout<<"Wrong trigger setup in runPicoLambdaAnaMaker!"<<endl;
    return;  
  }

  hfCuts->setCutNHitsFitMin(20); //for analysis (TTree)
  hfCuts->setCutNHitsFitMinHist(20); //for histograms

  hfCuts->setCutRequireHFT(false);
  hfCuts->setCutRequireISTorSST(false); //just for Run15!!!

  hfCuts->setCutNHitsFitnHitsMax(0.52); 

  // ---------------------------------------------------

  picoLambdaAnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay); //to setup secondary vertices as StHFPair


  hfCuts->setCutEta(1.);
  hfCuts->setCutPtMin(0.15); //global min. pT cut

  hfCuts->setCutDcaMin(0.3,StHFCuts::kPion); //original 0.3, loose 0.1
  //hfCuts->setCutDcaMin(0.01,StHFCuts::kKaon); 
  hfCuts->setCutDcaMin(0.1,StHFCuts::kProton); //original 0.1, lose 0.05

  //-----------lambda selection cuts----------------------------
  float dcaDaughters12Max;
  float decayLengthMin, decayLengthMax;
  float cosThetaMin, massMin, massMax;

  dcaDaughters12Max = 1.; //original 1, 0.5 for Run15

  decayLengthMin = 2.; //original 2, loose 0.5
  decayLengthMax = 25.;

  cosThetaMin = 0.996; //original 0.996, loose 0.98 or 0.99

  massMin = 0.9;
  massMax = 1.3; //default 1.3, 1.15 for Run15 to reduce running time

  hfCuts->setCutSecondaryPair(dcaDaughters12Max, decayLengthMin, decayLengthMax, cosThetaMin, massMin, massMax);
  //______________________________________________________________________________________________________________

  //-----------K0s selection cuts----------------------------
  float dcaDaughtersMax_K0s;
  float decayLengthMin_K0s, decayLengthMax_K0s;
  float cosThetaMin_K0s, massMin_K0s, massMax_K0s;
 
  //K0s cuts
  dcaDaughtersMax_K0s = 1.; //original 1., 0.5 for Run15

  decayLengthMin_K0s = 0.5; //original 0.5, loose 0.1
  decayLengthMax_K0s = 25.;

  cosThetaMin_K0s = 0.996; //original 0.996, loose 0.98 or 0.99

  massMin_K0s = 0.45;
  massMax_K0s = 0.55;

  hfCuts->setCutSecondaryPair_2(dcaDaughtersMax_K0s, decayLengthMin_K0s, decayLengthMax_K0s, cosThetaMin_K0s, massMin_K0s, massMax_K0s);

  //hfCuts->setCutSecondaryPairDcaToPvMax_2(0.5); //for K0s


  //Single track pt
  hfCuts->setCutPtRange(0.15,20.0,StHFCuts::kPion);
  hfCuts->setCutPtRange(0.15,20.0,StHFCuts::kKaon);
  hfCuts->setCutPtRange(0.15,20.0,StHFCuts::kProton);
  //___________________________________________________________________________________________________________

  
  //TPC setters
  hfCuts->setCutTPCNSigmaPion(3.0);
  hfCuts->setCutTPCNSigmaKaon(2.0); 
  hfCuts->setCutTPCNSigmaProton(2.0);

  //for TOF matching histograms
  hfCuts->setCutTPCNSigmaHadronHist(1.0, 1); //1 = pion
  hfCuts->setCutTPCNSigmaHadronHist(1.0, 2); //2 = kaon
  hfCuts->setCutTPCNSigmaHadronHist(1.0, 3); //3 = proton

  //TOF setters, need to set pt range as well
  hfCuts->setCutRequireStrictTOF(false); // setter for strict TOF for pions (or second particle in pair in general)

  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kProton);
  hfCuts->setCutPtotRangeHybridTOF(0.15,20.0,StHFCuts::kProton); 

  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kPion);
  hfCuts->setCutPtotRangeHybridTOF(0.15,20.0,StHFCuts::kPion);

  // ========================================================================================

  chain->Init();
  cout << "chain->Init();" << endl;
  int nEvents = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << nEvents << endl;

  for (Int_t i=0; i<nEvents; i++)
  {
    if(i%10000==0) cout << "Working on eventNumber " << i << endl;
    //if(i%500==0) cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if(iret) 
    { 
      cout << "Bad return code!" << iret << endl;
      break;
    }

    //total++;
  }

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();

  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  cout << "****************************************** " << endl;

  delete chain;

  // -- clean up if in read mode
  if (makerMode == StPicoHFMaker::kRead)
    gSystem->Exec(Form("rm -f %s", sInputFile.Data()));
}

