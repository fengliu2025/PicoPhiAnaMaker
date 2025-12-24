#include "StPicoPhiAnaMaker.h"
//#include "StPicoHFMaker/StHFCuts.h"
#include <iostream>
#include <cmath>

ClassImp(StPicoPhiAnaMaker)

using namespace std;

// _________________________________________________________
StPicoPhiAnaMaker::StPicoPhiAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,
               char const* inputHFListHFtree = "") :
  StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
  mDecayChannel(kChannel1), mOutFileBaseName(outputBaseFileName){


  // constructor
}

// _________________________________________________________
StPicoPhiAnaMaker::~StPicoPhiAnaMaker() {
  // destructor
}

// _________________________________________________________
int StPicoPhiAnaMaker::InitHF() {
  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written

  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());

  //QA histograms and TOF matching histograms
  //mOutList->Add(new TH1F("h_time_per_event","h_time_per_event", 2000., 0., 20.));

  if(isMakerMode() == StPicoHFMaker::kWrite)
  {

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    

    //create TTree to store Kaon candidates
    ntp_Kaon = new TTree("ntp_Kaon", "Kaon TTree"); //create TTree


    //---Set TTree branches------------------------------------------------------------------------------

    //event
    //ntp_kaon->Branch("runId", &runId, "runId/I");            //Int_t runId
    ntp_kaon->Branch("eventId", &eventId, "eventId/I");       //Int_t eventId
    ntp_Kaon->Branch("Vz", &Vz, "Vz/F"); //VzVzVPDmax
    //ntp_Kaon->Branch("VzVzVPDmax", &VzVzVPDmax, "VzVzVPDmax/F"); //VzVzVPDmax
    
    ntp_Kaon->Branch("mNTrigs", &mNTrigs, "mNTrigs/I");
    ntp_Kaon->Branch("mTrigId", mTrigId, "mTrigId[mNTrigs]/I");

    //high pt particles
    //ntp_Kaon->Branch("mNTrks", &mNTrks, "mNTrks/I");
    //ntp_Kaon->Branch("high_pt", high_pt, "high_pt[mNTrks]/F");  
    //ntp_Kaon->Branch("high_phi", high_phi, "high_phi[mNTrks]/F");
    //ntp_Kaon->Branch("high_eta", high_eta, "high_eta[mNTrks]/F");

    /*
    ntp_Kaon->Branch("lead_pt", &lead_pt, "lead_pt/F");               //Float_t p1_pt
    ntp_Kaon->Branch("lead_phi", &lead_phi, "lead_phi/F");             //Float_t p1_phi
    ntp_Kaon->Branch("lead_eta", &lead_eta, "lead_eta/F");             //Float_t p1_eta
    
    ntp_Kaon->Branch("sublead_pt", &sublead_pt, "sublead_pt/F");               //Float_t p1_pt
    ntp_Kaon->Branch("sublead_phi", &sublead_phi, "sublead_phi/F");             //Float_t p1_phi
    ntp_Kaon->Branch("sublead_eta", &sublead_eta, "sublead_eta/F");             //Float_t p1_eta
    */
    ntp_Kaon->Branch("NKaon", &NKaon, "NKaon/I");
    
    //Kaon
    ntp_Kaon->Branch("kaon_InEventID", kaon_InEventID, "kaon_InEventID[NKaon]/I");               //Float_t kaon_InEventID
    ntp_Kaon->Branch("kaon_pt", kaon_pt, "kaon_pt[NKaon]/F");                                    //Float_t kaon_pt
    ntp_Kaon->Branch("kaon_phi", kaon_phi, "kaon_phi[NKaon]/F");                                 //Float_t kaon_phi
    ntp_Kaon->Branch("kaon_eta", kaon_eta, "kaon_eta[NKaon]/F");                                 //Float_t kaon_eta
    ntp_Kaon->Branch("kaon_dca", kaon_dca, "kaon_dca[NKaon]/F");                                 //Float_t kaon_dca
    ntp_Kaon->Branch("kaon_ch", kaon_ch, "kaon_ch[NKaon]/I");                                    //Float_t kaon_ch
    ntp_Kaon->Branch("kaon_hasTOFinfo", kaon_hasTOFinfo, "kaon_hasTOFinfo[NKaon]/I");             //Float_t kaon_hasTOFinfo
    ntp_Kaon->Branch("kaon_dedx", kaon_dedx, "kaon_dedx[NKaon]/F");
    ntp_Kaon->Branch("kaon_beta", kaon_beta, "kaon_beta[NKaon]/F");

    
  }


  if(isMakerMode() == StPicoHFMaker::kAnalyze){
    h1D_PhiCandidates_mass_spectrum = new TH1D("h1D_PhiCandidates_mass_spectrum","h1D_PhiCandidates_mass_spectrum",1000,0.5,1.5);

  }

  if(isMakerMode() == StPicoHFMaker::kQA)
  {    
      
    mOutList->Add(new TH1F("h_piTPC","h_piTPC",100,0,10));
    mOutList->Add(new TH1F("h_kTPC","h_kTPC",100,0,10));
    mOutList->Add(new TH1F("h_pTPC","h_pTPC",100,0,10));

    mOutList->Add(new TH1F("h_piTPC_prim","h_piTPC_prim",100,0,10));
    mOutList->Add(new TH1F("h_kTPC_prim","h_kTPC_prim",100,0,10));
    mOutList->Add(new TH1F("h_pTPC_prim","h_pTPC_prim",100,0,10));

    mOutList->Add(new TH2F("h_piTOF","h_piTOF",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_kTOF","h_kTOF",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_pTOF","h_pTOF",100,0,10, 2, 0, 2));

    mOutList->Add(new TH2F("h_piTOF_20","h_piTOF_20",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_kTOF_20","h_kTOF_20",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_pTOF_20","h_pTOF_20",100,0,10, 2, 0, 2));
    
    mOutList->Add(new TH2F("h_piTOF_BEMC_match","h_piTOF_BEMC_match",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_kTOF_BEMC_match","h_kTOF_BEMC_match",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_pTOF_BEMC_match","h_pTOF_BEMC_match",100,0,10, 2, 0, 2));

    mOutList->Add(new TH2F("h_piTOF_1sig","h_piTOF_1sig",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_kTOF_1sig","h_kTOF_1sig",100,0,10, 2, 0, 2));
    mOutList->Add(new TH2F("h_pTOF_1sig","h_pTOF_1sig",100,0,10, 2, 0, 2));

    mOutList->Add(new TH2F("h_piTOFbeta","h_piTOFbeta",500,0,10, 500, 0, 5));
    mOutList->Add(new TH2F("h_kTOFbeta","h_kTOFbeta",500,0,10, 500, 0, 5));
    mOutList->Add(new TH2F("h_pTOFbeta","h_pTOFbeta",500,0,10, 500, 0, 5));
    
    
    mOutList->Add(new TH1F("h_piHFT","h_piHFT",100,0,10));
    mOutList->Add(new TH1F("h_kHFT","h_kHFT",100,0,10));
    mOutList->Add(new TH1F("h_pHFT","h_pHFT",100,0,10));
    
    mOutList->Add(new TH1F("h_piHFT_TOF","h_piHFT_TOF",100,0,10));
    mOutList->Add(new TH1F("h_kHFT_TOF","h_kHFT_TOF",100,0,10));
    mOutList->Add(new TH1F("h_pHFT_TOF","h_pHFT_TOF",100,0,10));
    
    

  }


  mRunNumber = 0;
  return kStOK;
}

// _________________________________________________________
void StPicoPhiAnaMaker::ClearHF(Option_t *opt="") {
  return;
}

// _________________________________________________________
int StPicoPhiAnaMaker::FinishHF() 
{
   if( isMakerMode() == StPicoHFMaker::kWrite )
   {    
     ntp_Kaon->Write(); //for candidates
   }

   if (isMakerMode() == StPicoHFMaker::kAnalyze){
    h1D_PhiCandidates_mass_spectrum->Write();
   }


  return kStOK;
}

// _________________________________________________________
int StPicoPhiAnaMaker::MakeHF() {

  //create and analyze D+- candidates

  //createCandidates() makes triplets of K and pi and checks that they pass cuts
  //analyzeCandidates() saves relevant information about the triplets to TTree

  //std::clock_t start1 = std::clock();//kvapil

  if (isMakerMode() == StPicoHFMaker::kWrite) {
    fillCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kRead) {
    // -- the reading back of the perviously written trees happens in the background
    //analyzeCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kAnalyze)
  {
    
    analyzeCandidates();
  

    /*
      
     */ //end comment TOF histograms

  }// end if (isMakerMode() == StPicoHFMaker::kAnalyze)
  
  else if (isMakerMode() == StPicoHFMaker::kQA)
  {
    createQA();
  }


/*
    double duration = (double) (std::clock() - start1) / (double) CLOCKS_PER_SEC;
  TH1F *h_time_per_event = static_cast<TH1F*>(mOutList->FindObject("h_time_per_event"));
  h_time_per_event->Fill(duration);
*/
  return kStOK;
}

//_____________________________________________________________
int StPicoPhiAnaMaker::fillCandidates(){
  // fill the Kaon candidates into the ntp_Kaon tree 
  cout<< "start fillCandidates" << std::endl;
  
  //fill some event information
  StPicoEvent *picoEvent = mPicoDst->event();
  runId = picoEvent->runId();
  eventId = picoEvent->eventId();
  Vz = mPrimVtx.z();
  VzVzVPDmax = fabs( mPrimVtx.z() - picoEvent->vzVpd());

  //fill the trigger 
  mNTrigs = 0 ; 
  unsigned int mStPhysics_TriggerIDs[5]={910001, 910003, 910013, 910802, 910804};

  for(int itrig=0; itrig<5;itrig++){
    if( picoEvent->isTrigger(mStPhysics_TriggerIDs[itrig]) ){
      mTrigId[mNTrigs] = mStPhysics_TriggerIDs[itrig];
      mNTrigs++;
    }
  }



  //fill the number of Kaon candidates in this event
  NKaon = mIdxPicoKaons.size();
  //fill the Kaon candicate information
  for( unsigned short idxKaon =0; idxKaon < mIdxPicoKaons.size(); ++idxKaon){
    StPicoTrack const *kaon =  mPicoDst->track(mIdxPicoKaons[idxKaon]);

    //check if the kaon candiate has TOF information 
    float kaonBetaBase = -1; 
    kaonBetaBase = mHFCuts->getTofBetaBase(kaon,mPicoDst->event()->bField());
    kaon_hasTOFinfo[idxKaon]= 0;
    if(!isnan(kaonBetaBase) && kaonBetaBase >0){
      kaon_hasTOFinfo[idxKaon] = 1;
    }

    //fill the beta of kaon candidates 
    kaon_beta[idxKaon]      = -999.0;

    Int_t bTofPidTraitsIndex =  kaon->bTofPidTraitsIndex();

    if(bTofPidTraitsIndex>=0){
      StPicoBTofPidTraits *btofPidTraits = mPicoDst->botPidTraits(bTofPidTraitsIndex);
      if(btofPidTraits->btofMatchFlag()){
        kaon_beta[idxKaon]  = btofPidTraits->btofBeta();
      }
    }
    


    kaon_InEventID[idxKaon] =  mIdxPicoKaons[idxKaon];
    kaon_phi[idxKaon]       =  kaon->gMom().Phi() ;
    kaon_eta[idxKaon]       =  kaon->gMom().Eta() ; 
    kaon_pt[idxKaon]        =  kaon->gMom().Pt()  ; 
    kaon_dca[idxKaon]       =  kaon->dDCAxy(mPrimVtx.x(),mPrimVtx.y()); 
    kaon_ch[idxKaon]        =  kaon->charge();   
    kaon_dedx[idxKaon]      =  kaon->dEdx();
    


  }//end idxKaon loop

  //if the number of Kaon Candiadates is non-zero, fill all branches into the tree 

  if(NKaon>0) ntp_Kaon->Fill();

  return kStOK;
}
//_____________________________________________________________


int StPicoPhiAnaMaker::analyzeCandidates(){
  //loop all Kaon Candadited, and pair those with opposite signs of charges. 
  for(int idxKaon = 0 ; idxKaon < mIdxPicoKaons.size(); ++ idxKaon){
      StPicoTrack const* kaon_1 = mPicoDst->track(mIdxPicoKaons[idxKaon]);
      if( fabs( kaon_1->charge() ) !=1 ) continue ;

      for(int jdxKaon = idxKaon+1 ; jdxKaon < mIdxPicoKaons.size(); ++ jdxKaon){
          StPicoTrack const* kaon_2 = mPicoDst->track(mIdxPicoKaons[jdxKaon]);
          if( fabs( kaon_1->charge() + kaon_2->charge() ) != 0 ) continue;

          TLorentzVector LorzVec_kaon_1;
          TLorenrzVector LorzVec_kaon_2;

          LorzVec_kaon_1.SetPtEtaPhiM( kaon_1->gMom().Pt(),
                                       kaon_1->gMom().Eta(),
                                       kaon_1->gMom().Phi(),
                                       mHFCuts->getHypotheticalMass(StHFCuts::kKaon)
                                     ); 
          LorzVec_kaon_2.SetPtEtaPhiM( kaon_1->gMom().Pt(),
                                       kaon_1->gMom().Eta(),
                                       kaon_1->gMom().Phi(),
                                       mHFCuts->getHypotheticalMass(StHFCuts::kKaon)                          
                                     );
          //pair 
          TLorentzVector LorzVec_phi = LorzVec_kaon_1 + LorzVec_kaon_2;
          h1D_PhiCandidates_mass_spectrum->Fill(LorzVec_phi.M());


      }

  }

  return kStOK;
}







int StPicoPhiAnaMaker::createQA()
{
 
  //QA and TOF matching histograms
  TH1F *h_piTPC = static_cast<TH1F*>(mOutList->FindObject("h_piTPC"));
  TH1F *h_kTPC = static_cast<TH1F*>(mOutList->FindObject("h_kTPC"));
  TH1F *h_pTPC = static_cast<TH1F*>(mOutList->FindObject("h_pTPC"));

  TH1F *h_piTPC_prim = static_cast<TH1F*>(mOutList->FindObject("h_piTPC_prim"));
  TH1F *h_kTPC_prim = static_cast<TH1F*>(mOutList->FindObject("h_kTPC_prim"));
  TH1F *h_pTPC_prim = static_cast<TH1F*>(mOutList->FindObject("h_pTPC_prim"));

  TH2F *h_piTOF = static_cast<TH2F*>(mOutList->FindObject("h_piTOF"));
  TH2F *h_kTOF = static_cast<TH2F*>(mOutList->FindObject("h_kTOF"));
  TH2F *h_pTOF = static_cast<TH2F*>(mOutList->FindObject("h_pTOF"));

  TH2F *h_piTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_20"));
  TH2F *h_kTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_20"));
  TH2F *h_pTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_20"));
  
  TH2F *h_piTOF_BEMC_match = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_BEMC_match"));
  TH2F *h_kTOF_BEMC_match = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_BEMC_match"));
  TH2F *h_pTOF_BEMC_match = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_BEMC_match"));

  TH2F *h_piTOF_1sig = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_1sig"));
  TH2F *h_kTOF_1sig = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_1sig"));
  TH2F *h_pTOF_1sig = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_1sig"));

  TH2F *h_piTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_piTOFbeta"));
  TH2F *h_kTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_kTOFbeta"));
  TH2F *h_pTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_pTOFbeta"));
  
  
  TH1F *h_piHFT =  static_cast<TH1F*>(mOutList->FindObject("h_piHFT"));
  TH1F *h_kHFT =  static_cast<TH1F*>(mOutList->FindObject("h_kHFT"));
  TH1F *h_pHFT =  static_cast<TH1F*>(mOutList->FindObject("h_pHFT"));
  
  TH1F *h_piHFT_TOF =  static_cast<TH1F*>(mOutList->FindObject("h_piHFT_TOF"));
  TH1F *h_kHFT_TOF =  static_cast<TH1F*>(mOutList->FindObject("h_kHFT_TOF"));
  TH1F *h_pHFT_TOF =  static_cast<TH1F*>(mOutList->FindObject("h_pHFT_TOF"));



  //primary vertex
  TVector3 pVtx = mPicoDst->event()->primaryVertex();

 

  for (unsigned short idxPion = 0; idxPion < mIdxPicoPions.size(); ++idxPion)
  {
    StPicoTrack const *pion = mPicoDst->track(mIdxPicoPions[idxPion]);
   
   
    
    //calculate beta of track from its momentum
    float piBeta = mHFCuts->getTofBetaBase(pion, mPicoDst->event()->bField()); //SL16j, Vanek

    //get beta of track from TOF, if TOF information is available          
    Int_t piTofAvailable = 0;
    if(!isnan(piBeta) && piBeta > 0)
    {
       piTofAvailable = 1;
       float tofPion = fabs(1. / piBeta - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion->gMom().Mag()*pion->gMom().Mag())));
       h_piTOFbeta->Fill(pion->gPt(),tofPion);
    }

    //fill QA and TOF matching histograms
    h_piTOF->Fill(pion->gPt(),piTofAvailable);
    
    if (mHFCuts->hasGoodNHitsFitMinHist(pion))
    {
      h_piTPC->Fill(pion->gPt());

      if(pion->isPrimary())
      {
        h_piTPC_prim->Fill(pion->pPt());
      }

      h_piTOF_20->Fill(pion->gPt(),piTofAvailable);
      
      if(pion->isBemcTrack())
      {
        h_piTOF_BEMC_match->Fill(pion->gPt(), piTofAvailable);     
      
      }
      
      if(pion->hasIstHit() || pion->hasSstHit())
      {
        h_piHFT->Fill(pion->gPt());
        
        if(piTofAvailable)
        {
          h_piHFT_TOF->Fill(pion->gPt());
        }
      }

    }

    if (mHFCuts->hasGoodNSigmaHist(pion, 1)) h_piTOF_1sig->Fill(pion->gPt(),piTofAvailable); //hasGoodNSigmaHist(pion, 1) => check nSigma of pion
      

  }
  
  
  for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
  {
    StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
   
   
    
    float kBeta = mHFCuts->getTofBetaBase(kaon, mPicoDst->event()->bField()); //SL16j, Vanek
    Int_t kTofAvailable = 0;
    if(!isnan(kBeta) && kBeta > 0)
    {
      kTofAvailable = 1;
      float tofKaon = fabs(1. / kBeta - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(kaon->gMom().Mag()*kaon->gMom().Mag())));
      h_pTOFbeta->Fill(kaon->gPt(),tofKaon);
    }

    h_pTOF->Fill(kaon->gPt(),kTofAvailable);

    if (mHFCuts->hasGoodNHitsFitMinHist(kaon))
    {
      h_pTPC->Fill(kaon->gPt());

      if(kaon->isPrimary())
      {
        h_pTPC_prim->Fill(kaon->pPt());
      }
    
      h_pTOF_20->Fill(kaon->gPt(),kTofAvailable);
      
      if(kaon->isBemcTrack())
      {
        h_kTOF_BEMC_match->Fill(kaon->gPt(), kTofAvailable);      
      }
      
      if(kaon->hasIstHit() || kaon->hasSstHit())
      {
        h_kHFT->Fill(kaon->gPt());
        
        if(kTofAvailable)
        {
          h_kHFT_TOF->Fill(kaon->gPt());
        }
      }
      
    }
        
    if (fabs(mHFCuts->hasGoodNSigmaHist(kaon, 3))) h_pTOF_1sig->Fill(kaon->gPt(),kTofAvailable); //hasGoodNSigmaHist(kaon, 3) => check nSigma of kaon
       

  }
  


  for (unsigned short idxProton = 0; idxProton < mIdxPicoProtons.size(); ++idxProton)
  {
    StPicoTrack const *proton = mPicoDst->track(mIdxPicoProtons[idxProton]);
   
   
    
    float pBeta = mHFCuts->getTofBetaBase(proton, mPicoDst->event()->bField()); //SL16j, Vanek
    Int_t pTofAvailable = 0;
    if(!isnan(pBeta) && pBeta > 0)
    {
      pTofAvailable = 1;
      float tofProton = fabs(1. / pBeta - sqrt(1+M_PROTON*M_PROTON/(proton->gMom().Mag()*proton->gMom().Mag())));
      h_pTOFbeta->Fill(proton->gPt(),tofProton);
    }

    h_pTOF->Fill(proton->gPt(),pTofAvailable);

    if (mHFCuts->hasGoodNHitsFitMinHist(proton))
    {
      h_pTPC->Fill(proton->gPt());

      if(proton->isPrimary())
      {
        h_pTPC_prim->Fill(proton->pPt());
      }
    
      h_pTOF_20->Fill(proton->gPt(),pTofAvailable);
      
      if(proton->isBemcTrack())
      {
        h_pTOF_BEMC_match->Fill(proton->gPt(), pTofAvailable);      
      }
      
      if(proton->hasIstHit() || proton->hasSstHit())
      {
        h_pHFT->Fill(proton->gPt());
        
        if(pTofAvailable)
        {
          h_pHFT_TOF->Fill(proton->gPt());
        }
      }
      
    }
        
    if (fabs(mHFCuts->hasGoodNSigmaHist(proton, 3))) h_pTOF_1sig->Fill(proton->gPt(),pTofAvailable); //hasGoodNSigmaHist(proton, 3) => check nSigma of proton
       

  }


       
  return 0;
}




// _________________________________________________________
bool StPicoPhiAnaMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
  // -- good hadron
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

// _________________________________________________________
bool StPicoPhiAnaMaker::isPion(StPicoTrack const * const trk) const {
  // -- good pion

   TVector3 t = trk->gMom();
   //if ( !(trk->hasIstHit() || trk->hasSstHit()) ) return false; //for testing - match pions to IST or SST for Run15
   if (!mHFCuts->hasGoodEta(t)) return false; 
   if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kPion) ) return false;
   if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
   return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoPhiAnaMaker::isKaon(StPicoTrack const * const trk) const {
  // -- good kaon

  TVector3 t = trk->gMom();
  if (!mHFCuts->hasGoodEta(t)) return false;
  if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kKaon) ) return false;
  if (!mHFCuts->cutMaxDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
}

// _________________________________________________________
bool StPicoPhiAnaMaker::isProton(StPicoTrack const * const trk) const {
  // -- good proton
  TVector3 t = trk->gMom();
  //if ( !(trk->hasIstHit() || trk->hasSstHit()) ) return false; //for testing - match pions to IST or SST for Run15
  if (!mHFCuts->hasGoodEta(t)) return false; 
  if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kProton) ) return false;
  if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kProton)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}

double StPicoPhiAnaMaker::DCA(StPicoTrack const * const trk, TVector3 const & vtx) const {
  // -- particle DCA

  return ((trk->origin() - vtx).Mag());
}




//-----------------------------------------------------------------------------
