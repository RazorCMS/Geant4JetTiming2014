#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#endif

#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "Rtypes.h"
#include "TInterpreter.h"
#include "TClonesArray.h"
#include<vector>
#include <TRandom3.h>                

// define structures to read in ntuple
#include "CMSAna/DataTree/interface/Types.h"
#include "CMSAna/DataTree/interface/TEventInfo.hh"
#include "CMSAna/DataTree/interface/TGenParticle.hh"
#include "CMSAna/DataTree/interface/TGenJet.hh"
#include "CMSAna/DataTree/interface/TPhoton.hh"
#include "CMSAna/DataTree/interface/TMuon.hh"
#include "CMSAna/DataTree/interface/TElectron.hh"
#include "CMSAna/DataTree/interface/TJet.hh"
#include "CMSAna/DataTree/interface/TPFCandidate.hh"

#include "CMSAna/Utils/CommonTools.hh"

// FastJet Stuff
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

using namespace std;

void processSignalPlusPileup() {
  
  bool processAll =  true;
  TRandom3 *MyRandom = new TRandom3( 9980 );
  const uint mixNPU = 140;

  //*****************************************************************************************
  // Histograms
  //*****************************************************************************************
  TH1F *jetEcalTime = new TH1F("jetEcalTime","; ECAL GEANT hit time; Number of Events", 10000, 0,50);
  TH1F *jetEcalEnergyWeightedTime = new TH1F("jetEcalEnergyWeightedTime","; ECAL GEANT hit energy-weighted time; Number of Events", 10000, 0,50);
  TH2F *jetEcalEnergyWeightedTimeEtaPhi = new TH2F("jetEcalEnergyWeightedTimeEtaPhi","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalEnergyWeightedTimeIEtaIPhi = new TH2F("jetEcalEnergyWeightedTimeIEtaIPhi","; i#eta; i#phi; Number of Events", 171, -85,86,360,1,361);
  TH2F *jetEcalEnergyWeightedTimeIXIY = new TH2F("jetEcalEnergyWeightedTimeIXIY","; ix; iy; Number of Events", 100, 1,101,100,1,101);
  TH2F *jetEcalEnergyIEtaIPhi = new TH2F("jetEcalEnergyIEtaIPhi","; i#eta; i#phi; Number of Events", 171, -85,86,360,1,361);
  TH2F *jetEcalTimeEnergyByCrystal = new TH2F("jetEcalTimeEnergyByCrystal", "; time; energy; Number of Events", 100, 5, 50, 100, 0, 300);
  TH2F *jetEcalEnergyIXIY = new TH2F("jetEcalEnergyIXIY","; ix; iy; Number of Events", 100, 1,101,100,1,101);
  TH2F *jetEcalCountIEtaIPhi = new TH2F("jetEcalCountIEtaIPhi","; i#eta; i#phi; Number of Events", 171, -85,86,360,1,361);
  TH2F *jetEcalCountIXIY = new TH2F("jetEcalCountIXIY","; ix; iy; Number of Events", 100, 1,101,100,1,101);
  TH2F *jetEcalTimeEtaPhi = new TH2F("jetEcalTimeEtaPhi","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalEnergyEtaPhi = new TH2F("jetEcalEnergyEtaPhi","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetParticleEnergyEtaPhi = new TH2F("jetParticleEnergyEtaPhi","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalCountEtaPhi = new TH2F("jetEcalCountEtaPhi","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalTimeEnergy = new TH2F("jetEcalTimeEnergy", "; time; energy; Number of Events", 10000, 0, 50, 10000, 0, 20);
  TH3F *jetEcal_EtaPhiTime = new TH3F("jetEcal_EtaPhiTime","; #eta; #phi; time; Number of Events", 50, -0.5,0.5,50,-0.5,0.5, 10000, 0, 50);
  TH3F *jetEcal_IEtaIPhiTime = new TH3F("jetEcal_IEtaIPhiTime","; i#eta; i#phi; time; Number of Events", 171, -85,86,360,1,361, 10000, 0, 50);
  TH3F *jetEcal_IXIYTime = new TH3F("jetEcal_IXIYTime","; ix; iy; time; Number of Events", 100, 1,101,100,1,101, 10000, 0, 50);

  TH1F *jetEcalTimePU = new TH1F("jetEcalTimePU","; ECAL GEANT hit time; Number of Events", 10000, 0,50);
  TH1F *jetEcalEnergyWeightedTimePU = new TH1F("jetEcalEnergyWeightedTimePU","; ECAL GEANT hit energy-weighted time; Number of Events", 10000, 0,50);
  TH2F *jetEcalEnergyWeightedTimeEtaPhiPU = new TH2F("jetEcalEnergyWeightedTimeEtaPhiPU","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalEnergyWeightedTimeIEtaIPhiPU = new TH2F("jetEcalEnergyWeightedTimeIEtaIPhiPU","; i#eta; i#phi; Number of Events", 171, -85,86,360,1,361);
  TH2F *jetEcalEnergyWeightedTimeIXIYPU = new TH2F("jetEcalEnergyWeightedTimeIXIYPU","; ix; iy; Number of Events", 100, 1,101,100,1,101);
  TH2F *jetEcalEnergyIEtaIPhiPU = new TH2F("jetEcalEnergyIEtaIPhiPU","; i#eta; i#phi; Number of Events", 171, -85,86,360,1,361);
  TH2F *jetEcalTimeEnergyByCrystalPU = new TH2F("jetEcalTimeEnergyByCrystalPU", "; time; energy; Number of Events", 100, 5, 50, 100, 0, 300);
  TH2F *jetEcalEnergyIXIYPU = new TH2F("jetEcalEnergyIXIYPU","; ix; iy; Number of Events", 100, 1,101,100,1,101);
  TH2F *jetEcalCountIEtaIPhiPU = new TH2F("jetEcalCountIEtaIPhiPU","; i#eta; i#phi; Number of Events", 171, -85,86,360,1,361);
  TH2F *jetEcalCountIXIYPU = new TH2F("jetEcalCountIXIYPU","; ix; iy; Number of Events", 100, 1,101,100,1,101);
  TH2F *jetEcalTimeEtaPhiPU = new TH2F("jetEcalTimeEtaPhiPU","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalEnergyEtaPhiPU = new TH2F("jetEcalEnergyEtaPhiPU","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetParticleEnergyEtaPhiPU = new TH2F("jetParticleEnergyEtaPhiPU","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalCountEtaPhiPU = new TH2F("jetEcalCountEtaPhiPU","; #eta; #phi; Number of Events", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalTimeEnergyPU = new TH2F("jetEcalTimeEnergyPU", "; time; energy; Number of Events", 10000, 0, 50, 10000, 0,  20);
  TH3F *jetEcal_EtaPhiTimePU = new TH3F("jetEcal_EtaPhiTimePU","; #eta; #phi; time; Number of Events", 50, -0.5,0.5,50,-0.5,0.5, 10000, 0, 50);
  TH3F *jetEcal_IEtaIPhiTimePU = new TH3F("jetEcal_IEtaIPhiTimePU","; #eta; #phi; time; Number of Events", 171, -85, 86, 360, 1, 361, 10000, 0, 50);
  TH3F *jetEcal_IXIYTimePU = new TH3F("jetEcal_IXIYTimePU","; ix; iy; time; Number of Events", 100, 1,101,100,1,101, 10000, 0, 50);

  TH1F *numHitsAbove4GeV = new TH1F("numHitsAbove4GeV", "; num rechits; Number of Events", 50, 0, 50);
  TH1F *numHitsAbove4GeVPU = new TH1F("numHitsAbove4GeVPU", "; num rechits; Number of Events", 50, 0, 50);
  
    TH2F *jetMeanTimeFirstTimeCumulative = new TH2F("jetMeanTimeVsFirstTimeCumulative", "; Mean time; First hit time; Number of events", 50, 8, 16, 50, 8, 16);
    TH2F *jetMeanTimeFirstTimeCumulativePU = new TH2F("jetMeanTimeVsFirstTimeCumulativePU", "; Mean time; First hit time; Number of events", 50, 8, 16, 50, 8, 16);
  //*****************************************************************************************
  // Set up input for analysis tree
  //*****************************************************************************************
  TFile *infile=0;  
  TFile *infilePU=0;
  TTree *eventTree=0;
  TTree *eventTreePU=0;

  infile = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/Timing/Simulation/VBFHgg/BACONNtuple_VBFHgg.root","READ");
  infilePU = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/Timing/Simulation/MinBias/BACONNtuple_MinBias_2.root","READ");

  // Data structures to store info from TTrees
  //infile->cd();
  cmsana::TEventInfo *info  = new cmsana::TEventInfo();
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");

  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("GenParticle", &genparticleArr); TBranch *genparticleBr = eventTree->GetBranch("GenParticle");
  eventTree->SetBranchAddress("GenJet", &genjetArr); TBranch *genjetBr = eventTree->GetBranch("GenJet");

  //infilePU->cd();
  cmsana::TEventInfo *infoPU  = new cmsana::TEventInfo();
  TClonesArray *genparticleArrPU = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArrPU = new TClonesArray("cmsana::TGenJet");

  eventTreePU = (TTree*)infilePU->Get("Events"); assert(eventTreePU);
  eventTreePU->SetBranchAddress("Info",     &infoPU);        TBranch *infoBrPU     = eventTreePU->GetBranch("Info");
  eventTreePU->SetBranchAddress("GenParticle", &genparticleArrPU); TBranch *genparticleBrPU = eventTreePU->GetBranch("GenParticle");
  eventTreePU->SetBranchAddress("GenJet", &genjetArrPU); TBranch *genjetBrPU = eventTreePU->GetBranch("GenJet");
  //*****************************************************************************************
  // Set up geant hits tree
  //*****************************************************************************************
  TFile *fSignal = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/Timing/Simulation/VBFHgg/g4simhitsEcal.root");
  TTree *g4tSignal = (TTree*)fSignal->Get("G4SIM");
  g4tSignal->SetMakeClass(1);

  // set up the Signal tree variables
  int numCrystsB= 0;
  int numCrystsE = 0;
  int ieta[61200];
  int iphi[61200];
  int ix[14648];
  int iy[14648];
  int iz[14648];
  std::vector<std::vector<float> > *eg4EBSignal = 0;
  std::vector<std::vector<float> > *tg4EBSignal = 0;
  std::vector<std::vector<float> > *prexg4EBSignal = 0;
  std::vector<std::vector<float> > *preyg4EBSignal = 0;
  std::vector<std::vector<float> > *prezg4EBSignal = 0;
  std::vector<std::vector<float> > *eg4EESignal = 0;
  std::vector<std::vector<float> > *tg4EESignal = 0;
  std::vector<std::vector<float> > *prexg4EESignal = 0;
  std::vector<std::vector<float> > *preyg4EESignal = 0;
  std::vector<std::vector<float> > *prezg4EESignal = 0;

  g4tSignal->SetBranchAddress("eg4EB",&eg4EBSignal);
  g4tSignal->SetBranchAddress("tg4EB",&tg4EBSignal);
  g4tSignal->SetBranchAddress("prexg4EB",&prexg4EBSignal);
  g4tSignal->SetBranchAddress("preyg4EB",&preyg4EBSignal);
  g4tSignal->SetBranchAddress("prezg4EB",&prezg4EBSignal);
  g4tSignal->SetBranchAddress("eg4EE",&eg4EESignal);
  g4tSignal->SetBranchAddress("tg4EE",&tg4EESignal);
  g4tSignal->SetBranchAddress("prexg4EE",&prexg4EESignal);
  g4tSignal->SetBranchAddress("preyg4EE",&preyg4EESignal);
  g4tSignal->SetBranchAddress("prezg4EE",&prezg4EESignal);
  g4tSignal->SetBranchAddress("ietag4EB", ieta);
  g4tSignal->SetBranchAddress("iphig4EB", iphi);
  g4tSignal->SetBranchAddress("ixg4EE", ix);
  g4tSignal->SetBranchAddress("iyg4EE", iy);
  g4tSignal->SetBranchAddress("izg4EE", iz);
  g4tSignal->SetBranchAddress("ng4EB", &numCrystsB);
  g4tSignal->SetBranchAddress("ng4EE", &numCrystsE);

  //*****************************************************************************************
  // Set up geant hits tree
  //*****************************************************************************************
  TFile *fPU = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/Timing/Simulation/MinBias/g4simhitsEcal_2.root");
  TTree *g4tPU = (TTree*)fPU->Get("G4SIM");
  g4tPU->SetMakeClass(1);

  // set up the PU tree variables
  int numCrystsBPU= 0;
  int numCrystsEPU = 0;
  int ietaPU[61200];
  int iphiPU[61200];
  int ixPU[14648];
  int iyPU[14648];
  int izPU[14648];
  std::vector<std::vector<float> > *eg4EBPU = 0;
  std::vector<std::vector<float> > *tg4EBPU = 0;
  std::vector<std::vector<float> > *prexg4EBPU = 0;
  std::vector<std::vector<float> > *preyg4EBPU = 0;
  std::vector<std::vector<float> > *prezg4EBPU = 0;
  std::vector<std::vector<float> > *eg4EEPU = 0;
  std::vector<std::vector<float> > *tg4EEPU = 0;
  std::vector<std::vector<float> > *prexg4EEPU = 0;
  std::vector<std::vector<float> > *preyg4EEPU = 0;
  std::vector<std::vector<float> > *prezg4EEPU = 0;

  g4tPU->SetBranchAddress("eg4EB",&eg4EBPU);
  g4tPU->SetBranchAddress("tg4EB",&tg4EBPU);
  g4tPU->SetBranchAddress("prexg4EB",&prexg4EBPU);
  g4tPU->SetBranchAddress("preyg4EB",&preyg4EBPU);
  g4tPU->SetBranchAddress("prezg4EB",&prezg4EBPU);
  g4tPU->SetBranchAddress("eg4EE",&eg4EEPU);
  g4tPU->SetBranchAddress("tg4EE",&tg4EEPU);
  g4tPU->SetBranchAddress("prexg4EE",&prexg4EEPU);
  g4tPU->SetBranchAddress("preyg4EE",&preyg4EEPU);
  g4tPU->SetBranchAddress("prezg4EE",&prezg4EEPU);
  g4tPU->SetBranchAddress("ietag4EB", ietaPU);
  g4tPU->SetBranchAddress("iphig4EB", iphiPU);
  g4tPU->SetBranchAddress("ixg4EE", ixPU);
  g4tPU->SetBranchAddress("iyg4EE", iyPU);
  g4tPU->SetBranchAddress("izg4EE", izPU);
  g4tPU->SetBranchAddress("ng4EB", &numCrystsBPU);
  g4tPU->SetBranchAddress("ng4EE", &numCrystsEPU);

  bool jetdone = false;

  uint eventToPlot = MyRandom->Integer(eventTree->GetEntries()); 
  cout << "Will plot a jet in event " << eventToPlot << endl;
  vector< vector<double> > jetEcalEnergyNoPU;
  vector< vector<double> > jetEcalEnergyPU;
  vector< vector<double> > jetEcalEta;
  //*****************************************************************************************
  // Do pileup Event Mixing
  //*****************************************************************************************
  //********************************************************
  // Events
  //********************************************************
  //for(UInt_t ientry=0; ientry < 2; ientry++) { 
  infile->cd(); //all these cd()'s are for debugging the segfaults I keep getting
  //(they shouldn't be necessary)
  for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) { 
    
    cout << "\n\nEvent " << ientry << "\n";

    //****************************************************
    //Generate list of minbias events to mix
    //****************************************************
    infilePU->cd();
    vector<uint> pileupEventsToMix;
    for (uint p=0; p < mixNPU; ++p) {
      uint eventToMix = MyRandom->Integer(eventTreePU->GetEntries());
      bool alreadyUsed = false;
      for (uint q=0; q<pileupEventsToMix.size(); ++q) {
        if (eventToMix == pileupEventsToMix[q]) alreadyUsed = true;
      }
      if (!alreadyUsed) pileupEventsToMix.push_back(eventToMix);
      //cout << "Event to Mix: " << pileupEventsToMix[p] << "\n";
    }

    vector<PseudoJet> VisibleParticles;
    VisibleParticles.clear();
    
    infile->cd();
    genparticleArr->Clear(); 
    genparticleBr->GetEntry(ientry);
    infoBr->GetEntry(ientry);

    double vertexZ = info->genVertexZ;
    cout << "vertexZ = " << vertexZ << "\n";

    vector<const cmsana::TGenParticle *> genPhotons;

    for(Int_t j=0; j<genparticleArr->GetEntries(); j++) {
      const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[j]);
      
      if (! (gen->status == 1)) continue;
      if ( abs(gen->pdgid) == 12 || abs(gen->pdgid) == 14 || abs(gen->pdgid) == 16) continue;
      
      if ( gen->pdgid == 22 && gen->motherPdgID == 25 ) {
        genPhotons.push_back(gen);
      }

      cmsana::FourVectorM genpartv;
      genpartv.SetPt(gen->pt);
      genpartv.SetEta(gen->eta);
      genpartv.SetPhi(gen->phi);
      genpartv.SetM(gen->mass);
      
      VisibleParticles.push_back(PseudoJet(genpartv.Px(),genpartv.Py(),genpartv.Pz(),genpartv.E()));
    }
  
    JetDefinition DefinitionNoPU(antikt_algorithm, 0.5);
    ClusterSequence SequenceNoPU(VisibleParticles, DefinitionNoPU);
    vector<PseudoJet> JetsNoPU = SequenceNoPU.inclusive_jets(0.5);

    //********************************************************
    // Print out original Jets with no pileup mixing
    //********************************************************
    cout << "GenJetsNoPU " << JetsNoPU.size() << "\n";
    fSignal->cd();
    vector<double> eventJetEcalEnergy;
    vector<double> eventJetEcalEta;
    for (int j=0; j< JetsNoPU.size()  ; ++j) {

      bool isPhoton = false;
      for (uint p=0; p<genPhotons.size(); ++p) {
        if (cmsana::deltaR(genPhotons[p]->eta, genPhotons[p]->phi,JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.3) isPhoton = true;
      }
      if (isPhoton) continue;
      if (JetsNoPU[j].pt() > 30 && fabs(JetsNoPU[j].eta()) < 3) {
        cout << "GenJet " << j << " : " << JetsNoPU[j].pt() << " " << JetsNoPU[j].eta() << " " << JetsNoPU[j].phi() << "\n";
        
        if(processAll){
	//get the GEANT hits for the jet
        fSignal->cd();
        g4tSignal->GetEntry(ientry);
          double jetEcalEnergySum = 0;
          cout << " : Finding Ecal Hits. \n";
          
          //look in barrel
          int hitsAbove4GeV = 0;
          for(Int_t l=0; l<eg4EBSignal->size(); l++) { //loop over crystals
              double earliestTime = 999;
              double totEnergy = 0;
              double weightedTime = 0;
            for(int m=0; m<eg4EBSignal->at(l).size(); m++) //loop over Geant hits within crystal
            {
              TVector3 vEBSignal;
              vEBSignal.SetXYZ((*prexg4EBSignal)[l][m], (*preyg4EBSignal)[l][m], (*prezg4EBSignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEBSignal.Eta(),vEBSignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalEnergySum += (*eg4EBSignal)[l][m];
                totEnergy += (*eg4EBSignal)[l][m];
                weightedTime += (*eg4EBSignal)[l][m]*(*tg4EBSignal)[l][m];
                if((*tg4EBSignal)[l][m] < earliestTime) earliestTime = (*tg4EBSignal)[l][m];
              }
            }
            weightedTime /= totEnergy;
            jetMeanTimeFirstTimeCumulative->Fill(weightedTime, earliestTime);
            if(totEnergy > 4000) hitsAbove4GeV++;
          }
            
          //look in endcap
          for(Int_t l=0; l<eg4EESignal->size(); l++) {
              double earliestTime = 999;
              double totEnergy = 0;
              double weightedTime = 0;
            for(int m=0; m<eg4EESignal->at(l).size(); m++)	    
            {
              TVector3 vEESignal;
              vEESignal.SetXYZ((*prexg4EESignal)[l][m], (*preyg4EESignal)[l][m], (*prezg4EESignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEESignal.Eta(),vEESignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalEnergySum += (*eg4EESignal)[l][m];
                totEnergy += (*eg4EESignal)[l][m];
                weightedTime += (*eg4EESignal)[l][m]*(*tg4EESignal)[l][m];
              }                
            }
            weightedTime /= totEnergy;
            jetMeanTimeFirstTimeCumulative->Fill(weightedTime, earliestTime);
            if(totEnergy > 4000) hitsAbove4GeV++;
          }
          numHitsAbove4GeV->Fill(hitsAbove4GeV);
          cout << "Hits above 4 GeV: " << hitsAbove4GeV << endl;

          cout << "jetEcalEnergySum [MeV]: " << jetEcalEnergySum << " , measured jet energy [GeV]: " << JetsNoPU[j].E() << " jet eta : " << JetsNoPU[j].eta() << "\n";
	  eventJetEcalEnergy.push_back(jetEcalEnergySum);
          eventJetEcalEta.push_back(fabs(JetsNoPU[j].eta()));
        }//endif(processAll)

        //Pick one jet to use for histograms
        if (ientry == eventToPlot && !jetdone && JetsNoPU[j].pt() > 30) 
        {
            infile->cd();
            cout << "Processing the 'special' jet (before PU added)" << endl;
          for(Int_t g=0; g<genparticleArr->GetEntries(); g++) {
            const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[g]);
            if (! (gen->status == 1)) continue;
            if ( abs(gen->pdgid) == 12 || abs(gen->pdgid) == 14 || abs(gen->pdgid) == 16) continue;

            if (cmsana::deltaR(gen->eta,gen->phi, JetsNoPU[j].eta(),  JetsNoPU[j].phi()) < 0.5) {
              jetParticleEnergyEtaPhi->Fill(gen->eta - JetsNoPU[j].eta(),cmsana::deltaPhi(gen->phi,JetsNoPU[j].phi()),gen->pt*cosh(gen->eta));
              cout << "gen particle: " << gen->pdgid << " " << gen->eta - JetsNoPU[j].eta() << " " << cmsana::deltaPhi(gen->phi,JetsNoPU[j].phi()) << " " << gen->pt*cosh(gen->eta) << "\n";

            }
          }

            fSignal->cd();
            g4tSignal->GetEntry(ientry);
          //look in barrel
          for(Int_t l=0; l<eg4EBSignal->size(); l++) {
            for(int m=0; m<eg4EBSignal->at(l).size(); m++)	    
            {
              TVector3 vEBSignal;
              vEBSignal.SetXYZ((*prexg4EBSignal)[l][m], (*preyg4EBSignal)[l][m], (*prezg4EBSignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEBSignal.Eta(),vEBSignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalTime->Fill((*tg4EBSignal)[l][m]);
                jetEcalEnergyWeightedTime->Fill((*tg4EBSignal)[l][m],(*eg4EBSignal)[l][m]);

                jetEcalEnergyWeightedTimeEtaPhi->Fill(vEBSignal.Eta() - JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*tg4EBSignal)[l][m] * (*eg4EBSignal)[l][m]);
                jetEcalEnergyWeightedTimeIEtaIPhi->Fill(ieta[l], iphi[l],(*tg4EBSignal)[l][m] * (*eg4EBSignal)[l][m]);
                jetEcalEnergyIEtaIPhi->Fill(ieta[l], iphi[l],(*eg4EBSignal)[l][m]);
                jetEcalCountIEtaIPhi->Fill(ieta[l], iphi[l]);
                jetEcalTimeEnergy->Fill((*tg4EBSignal)[l][m], (*eg4EBSignal)[l][m]);
                jetEcalTimeEtaPhi->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*tg4EBSignal)[l][m]);
                jetEcalEnergyEtaPhi->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*eg4EBSignal)[l][m]);
                jetEcalCountEtaPhi->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()));
                jetEcal_EtaPhiTime->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*tg4EBSignal)[l][m],(*eg4EBSignal)[l][m]);
                jetEcal_IEtaIPhiTime->Fill(ieta[l], iphi[l],(*tg4EBSignal)[l][m],(*eg4EBSignal)[l][m]);
                //cout << "Yeah " << ieta[l] << " " << iphi[l] << " " << (*tg4EBSignal)[l][m] << " " << (*eg4EBSignal)[l][m] << endl;
              }
            }
          }
            
          //look in endcap
          // run over endcap
          for(Int_t l=0; l<eg4EESignal->size(); l++) {
            for(int m=0; m<eg4EESignal->at(l).size(); m++)	    
            {
              TVector3 vEESignal;
              vEESignal.SetXYZ((*prexg4EESignal)[l][m], (*preyg4EESignal)[l][m], (*prezg4EESignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEESignal.Eta(),vEESignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalTime->Fill((*tg4EESignal)[l][m]);
                jetEcalEnergyWeightedTime->Fill((*tg4EESignal)[l][m],(*eg4EESignal)[l][m]);

                jetEcalEnergyWeightedTimeEtaPhi->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*tg4EESignal)[l][m] * (*eg4EESignal)[l][m]);
                jetEcalEnergyWeightedTimeIXIY->Fill(ix[l], iy[l], (*tg4EESignal)[l][m] * (*eg4EESignal)[l][m]);
                jetEcalEnergyIXIY->Fill(ix[l], iy[l],(*eg4EESignal)[l][m]);
                jetEcalCountIXIY->Fill(ix[l], iy[l]);
                jetEcalTimeEnergy->Fill((*tg4EESignal)[l][m], (*eg4EESignal)[l][m]);
                jetEcalTimeEtaPhi->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*tg4EESignal)[l][m]);
                jetEcalEnergyEtaPhi->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*eg4EESignal)[l][m]);
                jetEcalCountEtaPhi->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()));
                jetEcal_EtaPhiTime->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*tg4EESignal)[l][m],(*eg4EESignal)[l][m]);
                jetEcal_IXIYTime->Fill(ix[l], iy[l],(*tg4EESignal)[l][m],(*eg4EESignal)[l][m]);
              }                
            }
          }

          jetdone = true;

        } //end if(special jet)

      }//end loop over jets
    }
    if(processAll){
        jetEcalEnergyNoPU.push_back(eventJetEcalEnergy);
        jetEcalEta.push_back(eventJetEcalEta);
    }
    //********************************************************
    // Mix Events
    //********************************************************
    infilePU->cd();
    for(UInt_t pEvent=0; pEvent < pileupEventsToMix.size(); pEvent++) {

      genparticleArrPU->Clear(); 
      genparticleBrPU->GetEntry(pileupEventsToMix[pEvent]);

      for(Int_t j=0; j<genparticleArrPU->GetEntries(); j++) {
        const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArrPU)[j]);

        if (! (gen->status == 1)) continue;
        if ( abs(gen->pdgid) == 12 || abs(gen->pdgid) == 14 || abs(gen->pdgid) == 16) continue;
      
        cmsana::FourVectorM genpartv;
        genpartv.SetPt(gen->pt);
        genpartv.SetEta(gen->eta);
        genpartv.SetPhi(gen->phi);
        genpartv.SetM(gen->mass);
      
        VisibleParticles.push_back(PseudoJet(genpartv.Px(),genpartv.Py(),genpartv.Pz(),genpartv.E()) );
      }
    }

    //Loop over the original signal jets again, now with PU particles added
    cout << "Added Pileup!" << endl;
    vector<double> eventJetEcalEnergyWithPU;
    jetdone = false;
    for (int j=0; j< JetsNoPU.size()  ; ++j) {

      bool isPhoton = false;
      for (uint p=0; p<genPhotons.size(); ++p) {
        if (cmsana::deltaR(genPhotons[p]->eta, genPhotons[p]->phi,JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.3) isPhoton = true;
      }
      if (isPhoton) continue;

      if (JetsNoPU[j].pt() > 30 && fabs(JetsNoPU[j].eta()) < 3) {
        cout << "GenJet " << j << " : " << JetsNoPU[j].pt() << " " << JetsNoPU[j].eta() << " " << JetsNoPU[j].phi() << "\n";
      
        if(processAll){
        double jetEcalEnergySum = 0;

        cout << " : Finding Ecal Hits. \n";
        fSignal->cd();
        g4tSignal->GetEntry(ientry);
        //look in barrel -- signal
        for(Int_t l=0; l<eg4EBSignal->size(); l++) {
            for(int m=0; m<eg4EBSignal->at(l).size(); m++)	    
            {
              TVector3 vEBSignal;
              vEBSignal.SetXYZ((*prexg4EBSignal)[l][m], (*preyg4EBSignal)[l][m], (*prezg4EBSignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEBSignal.Eta(),vEBSignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalEnergySum += (*eg4EBSignal)[l][m];
              }
            }
        }
            
        //look in endcap -- signal
        for(Int_t l=0; l<eg4EESignal->size(); l++) {
            for(int m=0; m<eg4EESignal->at(l).size(); m++)	    
            {
              TVector3 vEESignal;
              vEESignal.SetXYZ((*prexg4EESignal)[l][m], (*preyg4EESignal)[l][m], (*prezg4EESignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEESignal.Eta(),vEESignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalEnergySum += (*eg4EESignal)[l][m];
              }                
            }
        }
        
        fPU->cd();
        for(int pEvent = 0; pEvent < pileupEventsToMix.size(); pEvent++){
            g4tPU->GetEntry(pileupEventsToMix[pEvent]);
            //look in barrel -- pileup
            for(Int_t l=0; l<eg4EBPU->size(); l++) {
                double totEnergy = 0;
                for(int m=0; m<eg4EBPU->at(l).size(); m++)	    
                {
                  TVector3 vEBPU;
                  vEBPU.SetXYZ((*prexg4EBPU)[l][m], (*preyg4EBPU)[l][m], (*prezg4EBPU)[l][m] - vertexZ);
                  if (cmsana::deltaR(vEBPU.Eta(),vEBPU.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                    jetEcalEnergySum += (*eg4EBPU)[l][m];
                    totEnergy += (*eg4EBPU)[l][m];
                  }
                }
            }
                
            //look in endcap -- pileup
            for(Int_t l=0; l<eg4EEPU->size(); l++) {
                double totEnergy = 0;
                for(int m=0; m<eg4EEPU->at(l).size(); m++)	    
                {
                  TVector3 vEEPU;
                  vEEPU.SetXYZ((*prexg4EEPU)[l][m], (*preyg4EEPU)[l][m], (*prezg4EEPU)[l][m] - vertexZ);
                  if (cmsana::deltaR(vEEPU.Eta(),vEEPU.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                      jetEcalEnergySum += (*eg4EEPU)[l][m];
                      totEnergy += (*eg4EEPU)[l][m];
                  }                
                }
           }
        }

          cout << "jetEcalEnergySum [MeV]: " << jetEcalEnergySum << " , measured jet energy [GeV]: " << JetsNoPU[j].E() << " jet eta : " << JetsNoPU[j].eta() << "\n";
	  eventJetEcalEnergyWithPU.push_back(jetEcalEnergySum);
          double sumBeforePU = jetEcalEnergyNoPU[ientry][eventJetEcalEnergyWithPU.size()-1];
          cout << "Energy before PU: " << sumBeforePU << "; Energy after PU: " << jetEcalEnergySum << "; Difference = " << jetEcalEnergySum - sumBeforePU << endl;
        } //endif(processAll)

        //Pick one jet to use for histograms
        if (ientry == eventToPlot && !jetdone && JetsNoPU[j].pt() > 30) 
        {
            cout << "Processing the 'special' jet (with PU added)" << endl;
            infile->cd();
         /* for(Int_t g=0; g<genparticleArr->GetEntries(); g++) {
            const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArr)[g]);
            if (! (gen->status == 1)) continue;
            if ( abs(gen->pdgid) == 12 || abs(gen->pdgid) == 14 || abs(gen->pdgid) == 16) continue;

            if (cmsana::deltaR(gen->eta,gen->phi, JetsNoPU[j].eta(),  JetsNoPU[j].phi()) < 0.5) {
              jetParticleEnergyEtaPhiPU->Fill(gen->eta - JetsNoPU[j].eta(),cmsana::deltaPhi(gen->phi,JetsNoPU[j].phi()),gen->pt*cosh(gen->eta));
              cout << "gen particle: " << gen->pdgid << " " << gen->eta - JetsNoPU[j].eta() << " " << cmsana::deltaPhi(gen->phi,JetsNoPU[j].phi()) << " " << gen->pt*cosh(gen->eta) << "\n";

            }
          }
*/
            fSignal->cd();
            g4tSignal->GetEntry(ientry);
          //look in barrel
          for(Int_t l=0; l<eg4EBSignal->size(); l++) {
            for(int m=0; m<eg4EBSignal->at(l).size(); m++)	    
            {
              TVector3 vEBSignal;
              vEBSignal.SetXYZ((*prexg4EBSignal)[l][m], (*preyg4EBSignal)[l][m], (*prezg4EBSignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEBSignal.Eta(),vEBSignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalTimePU->Fill((*tg4EBSignal)[l][m]);
                jetEcalEnergyWeightedTimePU->Fill((*tg4EBSignal)[l][m],(*eg4EBSignal)[l][m]);

                jetEcalEnergyWeightedTimeEtaPhiPU->Fill(vEBSignal.Eta() - JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*tg4EBSignal)[l][m] * (*eg4EBSignal)[l][m]);
                jetEcalEnergyWeightedTimeIEtaIPhiPU->Fill(ieta[l], iphi[l],(*tg4EBSignal)[l][m] * (*eg4EBSignal)[l][m]);
                jetEcalEnergyIEtaIPhiPU->Fill(ieta[l], iphi[l],(*eg4EBSignal)[l][m]);
                jetEcalCountIEtaIPhiPU->Fill(ieta[l], iphi[l]);
                jetEcalTimeEnergyPU->Fill((*tg4EBSignal)[l][m], (*eg4EBSignal)[l][m]);
                jetEcalTimeEtaPhiPU->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*tg4EBSignal)[l][m]);
                jetEcalEnergyEtaPhiPU->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*eg4EBSignal)[l][m]);
                jetEcalCountEtaPhiPU->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()));
                jetEcal_EtaPhiTimePU->Fill(vEBSignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBSignal.Phi(),JetsNoPU[j].phi()),(*tg4EBSignal)[l][m], (*eg4EBSignal)[l][m]);
                jetEcal_IEtaIPhiTimePU->Fill(ieta[l], iphi[l],(*tg4EBSignal)[l][m],(*eg4EBSignal)[l][m]);
              }
            }
          }
            
          //look in endcap
          // run over endcap
          for(Int_t l=0; l<eg4EESignal->size(); l++) {
            for(int m=0; m<eg4EESignal->at(l).size(); m++)	    
            {
              TVector3 vEESignal;
              vEESignal.SetXYZ((*prexg4EESignal)[l][m], (*preyg4EESignal)[l][m], (*prezg4EESignal)[l][m] - vertexZ);
              if (cmsana::deltaR(vEESignal.Eta(),vEESignal.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                jetEcalTimePU->Fill((*tg4EESignal)[l][m]);
                jetEcalEnergyWeightedTimePU->Fill((*tg4EESignal)[l][m],(*eg4EESignal)[l][m]);

                jetEcalEnergyWeightedTimeEtaPhiPU->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*tg4EESignal)[l][m] * (*eg4EESignal)[l][m]);
                jetEcalEnergyWeightedTimeIXIYPU->Fill(ix[l], iy[l], (*tg4EESignal)[l][m] * (*eg4EESignal)[l][m]);
                jetEcalEnergyIXIYPU->Fill(ix[l], iy[l],(*eg4EESignal)[l][m]);
                jetEcalCountIXIYPU->Fill(ix[l], iy[l]);
                jetEcalTimeEnergyPU->Fill((*tg4EESignal)[l][m], (*eg4EESignal)[l][m]);
                jetEcalTimeEtaPhiPU->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*tg4EESignal)[l][m]);
                jetEcalEnergyEtaPhiPU->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*eg4EESignal)[l][m]);
                jetEcalCountEtaPhiPU->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()));
                jetEcal_EtaPhiTimePU->Fill(vEESignal.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEESignal.Phi(),JetsNoPU[j].phi()),(*tg4EESignal)[l][m],(*eg4EESignal)[l][m]);
                jetEcal_IXIYTimePU->Fill(ix[l], iy[l],(*tg4EESignal)[l][m], (*eg4EESignal)[l][m]);
              }                
            }
          }

          //include pileup hits
          fPU->cd();
          for(UInt_t iEvent=0; iEvent < pileupEventsToMix.size(); iEvent++) {
            g4tPU->GetEntry(pileupEventsToMix[iEvent]);
            /*  
            genparticleArrPU->Clear(); 
            genparticleBrPU->GetEntry(pileupEventsToMix[iEvent]);
            for(Int_t g=0; g<genparticleArrPU->GetEntries(); g++) {
              const cmsana::TGenParticle *gen = (cmsana::TGenParticle*)((*genparticleArrPU)[g]);
              if (! (gen->status == 1)) continue;
              if ( abs(gen->pdgid) == 12 || abs(gen->pdgid) == 14 || abs(gen->pdgid) == 16) continue;
                
              if (cmsana::deltaR(gen->eta,gen->phi, JetsNoPU[j].eta(),  JetsNoPU[j].phi()) < 0.5) {
                jetParticleEnergyEtaPhiPU->Fill(gen->eta - JetsNoPU[j].eta(),cmsana::deltaPhi(gen->phi,JetsNoPU[j].phi()),gen->pt*cosh(gen->eta));                  
              }
            }
              */
            //look in barrel
            for(Int_t l=0; l<eg4EBPU->size(); l++) {
              for(int m=0; m<eg4EBPU->at(l).size(); m++)	    
              {
                TVector3 vEBPU;
                vEBPU.SetXYZ((*prexg4EBPU)[l][m], (*preyg4EBPU)[l][m], (*prezg4EBPU)[l][m]);
                if (cmsana::deltaR(vEBPU.Eta(),vEBPU.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                  jetEcalTimePU->Fill((*tg4EBPU)[l][m]);
                  jetEcalEnergyWeightedTimePU->Fill((*tg4EBPU)[l][m],(*eg4EBPU)[l][m]);

                  jetEcalEnergyWeightedTimeEtaPhiPU->Fill(vEBPU.Eta() - JetsNoPU[j].eta(),cmsana::deltaPhi(vEBPU.Phi(),JetsNoPU[j].phi()),(*tg4EBPU)[l][m] * (*eg4EBPU)[l][m]);
                  jetEcalEnergyWeightedTimeIEtaIPhiPU->Fill(ietaPU[l], iphiPU[l],(*tg4EBPU)[l][m] * (*eg4EBPU)[l][m]);
                  jetEcalEnergyIEtaIPhiPU->Fill(ietaPU[l], iphiPU[l],(*eg4EBPU)[l][m]);
                  jetEcalCountIEtaIPhiPU->Fill(ietaPU[l], iphiPU[l]);
                  jetEcalTimeEnergyPU->Fill((*tg4EBPU)[l][m], (*eg4EBPU)[l][m]);
                  jetEcalTimeEtaPhiPU->Fill(vEBPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBPU.Phi(),JetsNoPU[j].phi()),(*tg4EBPU)[l][m]);
                  jetEcalEnergyEtaPhiPU->Fill(vEBPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBPU.Phi(),JetsNoPU[j].phi()),(*eg4EBPU)[l][m]);
                  jetEcalCountEtaPhiPU->Fill(vEBPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBPU.Phi(),JetsNoPU[j].phi()));
                  jetEcal_EtaPhiTimePU->Fill(vEBPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEBPU.Phi(),JetsNoPU[j].phi()),(*tg4EBPU)[l][m], (*eg4EBPU)[l][m]);
                  jetEcal_IEtaIPhiTimePU->Fill(ietaPU[l], iphiPU[l],(*tg4EBPU)[l][m], (*eg4EBPU)[l][m]);
                }
              }
            }
            
            //look in endcap
            // run over endcap
            for(Int_t l=0; l<eg4EEPU->size(); l++) {
              for(int m=0; m<eg4EEPU->at(l).size(); m++)	    
              {
                TVector3 vEEPU;
                vEEPU.SetXYZ((*prexg4EEPU)[l][m], (*preyg4EEPU)[l][m], (*prezg4EEPU)[l][m]);
                if (cmsana::deltaR(vEEPU.Eta(),vEEPU.Phi(),JetsNoPU[j].eta(),JetsNoPU[j].phi()) < 0.5) {
                  jetEcalTimePU->Fill((*tg4EEPU)[l][m]);
                  jetEcalEnergyWeightedTimePU->Fill((*tg4EEPU)[l][m],(*eg4EEPU)[l][m]);

                  jetEcalEnergyWeightedTimeEtaPhiPU->Fill(vEEPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEEPU.Phi(),JetsNoPU[j].phi()),(*tg4EEPU)[l][m] * (*eg4EEPU)[l][m]);
                  jetEcalEnergyWeightedTimeIXIYPU->Fill(ixPU[l], iyPU[l], (*tg4EEPU)[l][m] * (*eg4EEPU)[l][m]);
                  jetEcalEnergyIXIYPU->Fill(ixPU[l], iyPU[l],(*eg4EEPU)[l][m]);
                  jetEcalCountIXIYPU->Fill(ixPU[l], iyPU[l]);
                  jetEcalTimeEnergyPU->Fill((*tg4EEPU)[l][m], (*eg4EEPU)[l][m]);
                  jetEcalTimeEtaPhiPU->Fill(vEEPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEEPU.Phi(),JetsNoPU[j].phi()),(*tg4EEPU)[l][m]);
                  jetEcalEnergyEtaPhiPU->Fill(vEEPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEEPU.Phi(),JetsNoPU[j].phi()),(*eg4EEPU)[l][m]);                    
                  jetEcalCountEtaPhiPU->Fill(vEEPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEEPU.Phi(),JetsNoPU[j].phi()));
                  jetEcal_EtaPhiTimePU->Fill(vEEPU.Eta()- JetsNoPU[j].eta(),cmsana::deltaPhi(vEEPU.Phi(),JetsNoPU[j].phi()),(*tg4EEPU)[l][m],(*eg4EEPU)[l][m]);
                  jetEcal_IXIYTimePU->Fill(ixPU[l], iyPU[l],(*tg4EEPU)[l][m], (*eg4EEPU)[l][m]);
                }
              }
            }
          } // loop over mixed pileup events
          jetdone = true;

        } //end if(special jet)

      } //end if(jet passes PT cut)
    }
    jetEcalEnergyPU.push_back(eventJetEcalEnergyWithPU);

    //Now recluster the jets, including pileup
    JetDefinition Definition(antikt_algorithm, 0.5);
    ClusterSequence Sequence(VisibleParticles, Definition);
    vector<PseudoJet> Jets = Sequence.inclusive_jets(0.5);
  
    int jetcount = 0;

    cout << "Particles " << VisibleParticles.size() << "\n";
    cout << "GenJets " << Jets.size() << "\n";
    for (int j=0; j< Jets.size()  ; ++j) {

      bool isPhoton = false;
      for (uint p=0; p<genPhotons.size(); ++p) {
        if (cmsana::deltaR(genPhotons[p]->eta, genPhotons[p]->phi,Jets[j].eta(),Jets[j].phi()) < 0.3) isPhoton = true;
      }
      if (isPhoton) continue;

      if (Jets[j].pt() > 30 && fabs(Jets[j].eta()) < 3) {
        //cout << "GenJet " << j << " : " << Jets[j].pt() << " " << Jets[j].eta() << " " << Jets[j].phi() << "\n";
      
	jetcount++;
      }   
    } // loop over jets
    cout << "Number of GenJets after PU: " << jetcount << endl;

  } // loop over events

  TFile *file = new TFile("jetTimingStudy_signalpluspileup.root", "recreate");
  file->cd();
  //Make a histogram of the change in jet energy due to pileup
  TH2D* jetEnergyChangeDueToPileup = new TH2D("jetEnergyChangeDueToPileup", "; #eta; energy ; Change in jet energy due to PU", 60, 0.0, 3.0, 100, 0, 200);
  if(processAll){
  for(int j = 0; j < jetEcalEnergyNoPU.size(); j++){
      for(int k = 0; k < jetEcalEnergyNoPU[j].size(); k++){
          jetEnergyChangeDueToPileup->Fill(jetEcalEta[j][k], (jetEcalEnergyPU[j][k] - jetEcalEnergyNoPU[j][k])/1000.0);
      }
  }
  } //endif(processAll)
  //*****************************************************************************************
  //Rescale jet time histograms
  //*****************************************************************************************

  for (UInt_t b=0; int(b)<jetEcalEnergyWeightedTimeEtaPhi->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalEnergyWeightedTimeEtaPhi->GetYaxis()->GetNbins()+2; ++c) {
      jetEcalEnergyWeightedTimeEtaPhi->SetBinContent(b,c,jetEcalEnergyWeightedTimeEtaPhi->GetBinContent(b,c)/(jetEcalEnergyEtaPhi->GetBinContent(b,c)*jetEcalCountEtaPhi->GetBinContent(b,c)));
      jetEcalTimeEtaPhi->SetBinContent(b,c,jetEcalTimeEtaPhi->GetBinContent(b,c)/jetEcalCountEtaPhi->GetBinContent(b,c));
    }
  }
  for (UInt_t b=0; int(b)<jetEcalEnergyWeightedTimeEtaPhiPU->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalEnergyWeightedTimeEtaPhiPU->GetYaxis()->GetNbins()+2; ++c) {
      jetEcalEnergyWeightedTimeEtaPhiPU->SetBinContent(b,c,jetEcalEnergyWeightedTimeEtaPhiPU->GetBinContent(b,c)/(jetEcalEnergyEtaPhiPU->GetBinContent(b,c)*jetEcalCountEtaPhiPU->GetBinContent(b,c)));
      jetEcalTimeEtaPhiPU->SetBinContent(b,c,jetEcalTimeEtaPhiPU->GetBinContent(b,c)/jetEcalCountEtaPhiPU->GetBinContent(b,c));
    }
  }
//for ieta/iphi
  for (UInt_t b=0; int(b)<jetEcalEnergyWeightedTimeIEtaIPhi->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalEnergyWeightedTimeIEtaIPhi->GetYaxis()->GetNbins()+2; ++c) {
      jetEcalEnergyWeightedTimeIEtaIPhi->SetBinContent(b,c,jetEcalEnergyWeightedTimeIEtaIPhi->GetBinContent(b,c)/(jetEcalEnergyIEtaIPhi->GetBinContent(b,c)*jetEcalCountIEtaIPhi->GetBinContent(b,c)));
    }
  }
//for ix/iy
  for (UInt_t b=0; int(b)<jetEcalEnergyWeightedTimeIXIY->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalEnergyWeightedTimeIXIY->GetYaxis()->GetNbins()+2; ++c) {
      jetEcalEnergyWeightedTimeIXIY->SetBinContent(b,c,jetEcalEnergyWeightedTimeIXIY->GetBinContent(b,c)/(jetEcalEnergyIXIY->GetBinContent(b,c)*jetEcalCountIXIY->GetBinContent(b,c)));
    }
  }


//for PU
  for (UInt_t b=0; int(b)<jetEcalEnergyWeightedTimeIEtaIPhiPU->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalEnergyWeightedTimeIEtaIPhiPU->GetYaxis()->GetNbins()+2; ++c) {
      jetEcalEnergyWeightedTimeIEtaIPhiPU->SetBinContent(b,c,jetEcalEnergyWeightedTimeIEtaIPhiPU->GetBinContent(b,c)/(jetEcalEnergyIEtaIPhiPU->GetBinContent(b,c)*jetEcalCountIEtaIPhiPU->GetBinContent(b,c)));
    }
  }
  for (UInt_t b=0; int(b)<jetEcalEnergyWeightedTimeIXIYPU->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalEnergyWeightedTimeIXIYPU->GetYaxis()->GetNbins()+2; ++c) {
      jetEcalEnergyWeightedTimeIXIYPU->SetBinContent(b,c,jetEcalEnergyWeightedTimeIXIYPU->GetBinContent(b,c)/(jetEcalEnergyIXIYPU->GetBinContent(b,c)*jetEcalCountIXIYPU->GetBinContent(b,c)));
    }
  }

  //*****************************************************************************************
  //Histogram the first hit time and the mean hit time for the 'selected' jet
  //*****************************************************************************************
  TH2F *jetEcalFirstArrivalTimeEtaPhi = new TH2F("jetEcalFirstArrivalTimeEtaPhi","; #eta; #phi; First Arrival Time [ns]", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalMeanTimeEtaPhi = new TH2F("jetEcalMeanTimeEtaPhi","; #eta; #phi; Mean Time [ns]", 50, -0.5,0.5,50,-0.5,0.5);
  TH1F *jetEcalMeanTime = new TH1F("jetEcalMeanTime", "Mean time [ns]", 50, 8, 16);
  TH1F *jetEcalFirstTime = new TH1F("jetEcalFirstTime", "First hit time [ns]", 50, 8, 16);
  TH2F *jetEcalMeanVsFirstTime = new TH2F("jetEcalMeanVsFirstTime", "Mean time vs First hit time [ns]", 50, 8, 16, 50, 8, 16);
  for (UInt_t b=0; int(b)<jetEcalFirstArrivalTimeEtaPhi->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalFirstArrivalTimeEtaPhi->GetYaxis()->GetNbins()+2; ++c) {

      double firstArrivalTime = 0;
      double energySum = 0;
      double energyWeightedTime = 0;
      for (UInt_t d=0; int(d)<jetEcal_EtaPhiTime->GetZaxis()->GetNbins()+2; ++d) {
        if (jetEcal_EtaPhiTime->GetBinContent(b,c,d) > 0) {
          firstArrivalTime = jetEcal_EtaPhiTime->GetZaxis()->GetBinCenter(d);
          break;
        }
      }

      for (UInt_t d=0; int(d)<jetEcal_EtaPhiTime->GetZaxis()->GetNbins()+2; ++d) {
        energySum += jetEcal_EtaPhiTime->GetBinContent(b,c,d);
        energyWeightedTime += jetEcal_EtaPhiTime->GetZaxis()->GetBinCenter(d) * jetEcal_EtaPhiTime->GetBinContent(b,c,d);
      }

      if (energySum > 50) {
        jetEcalFirstArrivalTimeEtaPhi->SetBinContent(b,c,firstArrivalTime);
        jetEcalFirstTime->Fill(firstArrivalTime);
        jetEcalMeanTimeEtaPhi->SetBinContent(b,c,energyWeightedTime/energySum);
        jetEcalMeanTime->Fill(energyWeightedTime/energySum);
        jetEcalMeanVsFirstTime->Fill(energyWeightedTime/energySum, firstArrivalTime);
      } else {
        jetEcalFirstArrivalTimeEtaPhi->SetBinContent(b,c,0);
        jetEcalMeanTimeEtaPhi->SetBinContent(b,c,0);
      }
    }
  }

  TH2F *jetEcalFirstArrivalTimeEtaPhiPU = new TH2F("jetEcalFirstArrivalTimeEtaPhiPU","; #eta; #phi; First Arrival Time [ns]", 50, -0.5,0.5,50,-0.5,0.5);
  TH2F *jetEcalMeanTimeEtaPhiPU = new TH2F("jetEcalMeanTimeEtaPhiPU","; #eta; #phi; Mean Time [ns]", 50, -0.5,0.5,50,-0.5,0.5);
  TH1F *jetEcalMeanTimePU = new TH1F("jetEcalMeanTimePU", "Mean time [ns]", 50, 8, 16);
  TH1F *jetEcalFirstTimePU = new TH1F("jetEcalFirstTimePU", "First hit time [ns]", 50, 8, 16);
  TH2F *jetEcalMeanVsFirstTimePU = new TH2F("jetEcalMeanVsFirstTimePU", "Mean time vs First hit time [ns]", 50, 8, 16, 50, 8, 16);
  for (UInt_t b=0; int(b)<jetEcalFirstArrivalTimeEtaPhiPU->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalFirstArrivalTimeEtaPhiPU->GetYaxis()->GetNbins()+2; ++c) {

      double firstArrivalTime = 0;
      double energySum = 0;
      double energyWeightedTime = 0;
      for (UInt_t d=0; int(d)<jetEcal_EtaPhiTimePU->GetZaxis()->GetNbins()+2; ++d) {
        if (jetEcal_EtaPhiTimePU->GetBinContent(b,c,d) > 0) {
          firstArrivalTime = jetEcal_EtaPhiTimePU->GetZaxis()->GetBinCenter(d);
          break;
        }
      }

      for (UInt_t d=0; int(d)<jetEcal_EtaPhiTimePU->GetZaxis()->GetNbins()+2; ++d) {
        energySum += jetEcal_EtaPhiTimePU->GetBinContent(b,c,d);
        energyWeightedTime += jetEcal_EtaPhiTimePU->GetZaxis()->GetBinCenter(d) * jetEcal_EtaPhiTimePU->GetBinContent(b,c,d);
      }

      if (energySum > 50) {
        jetEcalFirstArrivalTimeEtaPhiPU->SetBinContent(b,c,firstArrivalTime);
        jetEcalFirstTimePU->Fill(firstArrivalTime);
        jetEcalMeanTimeEtaPhiPU->SetBinContent(b,c,energyWeightedTime/energySum);
        jetEcalMeanTimePU->Fill(energyWeightedTime/energySum);
        jetEcalMeanVsFirstTimePU->Fill(energyWeightedTime/energySum, firstArrivalTime);
      } else {
        jetEcalFirstArrivalTimeEtaPhiPU->SetBinContent(b,c,0);
        jetEcalMeanTimeEtaPhiPU->SetBinContent(b,c,0);
      }
    }
  }

  //Now do it on a crystal-by-crystal basis
  TH2F *jetEcalFirstArrivalTimeIEtaIPhi = new TH2F("jetEcalFirstArrivalTimeIEtaIPhi","; i#eta; i#phi; First Arrival Time [ns]", 171, -85,86,360,1,361);
  TH2F *jetEcalMeanTimeIEtaIPhi = new TH2F("jetEcalMeanTimeIEtaIPhi","; i#eta; i#phi; First Arrival Time [ns]", 171, -85,86,360,1,361);
  TH2F *jetEcalFirstArrivalTimeIXIY = new TH2F("jetEcalFirstArrivalTimeIXIY","; ix; iy; First Arrival Time [ns]", 100, 1, 101, 100, 1, 101);
  TH2F *jetEcalMeanTimeIXIY = new TH2F("jetEcalMeanTimeIXIY","; ix; iy; First Arrival Time [ns]", 100, 1, 101, 100, 1, 101);
  TH1F *jetEcalMeanTimeCry = new TH1F("jetEcalMeanTimeCry", "Mean time [ns]", 50, 8, 16);
  TH1F *jetEcalFirstTimeCry = new TH1F("jetEcalFirstTimeCry", "First hit time [ns]", 50, 8, 16);
  TH2F *jetEcalMeanVsFirstTimeCry = new TH2F("jetEcalMeanVsFirstTimeCry", "Mean time vs First hit time [ns]", 50, 8, 16, 50, 8, 16);
  //barrel
  for (UInt_t b=0; int(b)<jetEcalFirstArrivalTimeIEtaIPhi->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalFirstArrivalTimeIEtaIPhi->GetYaxis()->GetNbins()+2; ++c) {

      double firstArrivalTime = 0;
      double energySum = 0;
      double energyWeightedTime = 0;
      for (UInt_t d=0; int(d)<jetEcal_IEtaIPhiTime->GetZaxis()->GetNbins()+2; ++d) {
        if (jetEcal_IEtaIPhiTime->GetBinContent(b,c,d) > 0) {
          firstArrivalTime = jetEcal_IEtaIPhiTime->GetZaxis()->GetBinCenter(d);
          break;
        }
      }

      for (UInt_t d=0; int(d)<jetEcal_IEtaIPhiTime->GetZaxis()->GetNbins()+2; ++d) {
        energySum += jetEcal_IEtaIPhiTime->GetBinContent(b,c,d);
        energyWeightedTime += jetEcal_IEtaIPhiTime->GetZaxis()->GetBinCenter(d) * jetEcal_IEtaIPhiTime->GetBinContent(b,c,d);
      }

      if (energySum > 0) jetEcalTimeEnergyByCrystal->Fill(energyWeightedTime/energySum, energySum);
      if (energySum > 50) {
        jetEcalFirstArrivalTimeIEtaIPhi->SetBinContent(b,c,firstArrivalTime);
        cout << "Setting bin " << b << " " << c << " to " << firstArrivalTime << endl;
        jetEcalFirstTimeCry->Fill(firstArrivalTime);
        jetEcalMeanTimeIEtaIPhi->SetBinContent(b,c,energyWeightedTime/energySum);
        jetEcalMeanTimeCry->Fill(energyWeightedTime/energySum);
        jetEcalMeanVsFirstTimeCry->Fill(energyWeightedTime/energySum, firstArrivalTime);
      } else {
        jetEcalFirstArrivalTimeIEtaIPhi->SetBinContent(b,c,0);
        jetEcalMeanTimeIEtaIPhi->SetBinContent(b,c,0);
      }
    }
  }
  //endcaps
  for (UInt_t b=0; int(b)<jetEcalFirstArrivalTimeIXIY->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalFirstArrivalTimeIXIY->GetYaxis()->GetNbins()+2; ++c) {

      double firstArrivalTime = 0;
      double energySum = 0;
      double energyWeightedTime = 0;
      for (UInt_t d=0; int(d)<jetEcal_IXIYTime->GetZaxis()->GetNbins()+2; ++d) {
        if (jetEcal_IXIYTime->GetBinContent(b,c,d) > 0) {
          firstArrivalTime = jetEcal_IXIYTime->GetZaxis()->GetBinCenter(d);
          break;
        }
      }

      for (UInt_t d=0; int(d)<jetEcal_IXIYTime->GetZaxis()->GetNbins()+2; ++d) {
        energySum += jetEcal_IXIYTime->GetBinContent(b,c,d);
        energyWeightedTime += jetEcal_IXIYTime->GetZaxis()->GetBinCenter(d) * jetEcal_IXIYTime->GetBinContent(b,c,d);
      }

      if (energySum > 0) jetEcalTimeEnergyByCrystal->Fill(energyWeightedTime/energySum, energySum);
      if (energySum > 50) {
        jetEcalFirstArrivalTimeIXIY->SetBinContent(b,c,firstArrivalTime);
        jetEcalFirstTimeCry->Fill(firstArrivalTime);
        jetEcalMeanTimeIXIY->SetBinContent(b,c,energyWeightedTime/energySum);
        jetEcalMeanTimeCry->Fill(energyWeightedTime/energySum);
        jetEcalMeanVsFirstTimeCry->Fill(energyWeightedTime/energySum, firstArrivalTime);
      } else {
        jetEcalFirstArrivalTimeIXIY->SetBinContent(b,c,0);
        jetEcalMeanTimeIXIY->SetBinContent(b,c,0);
      }
    }
  }

  //for signal plus pileup
  TH2F *jetEcalFirstArrivalTimeIEtaIPhiPU = new TH2F("jetEcalFirstArrivalTimeIEtaIPhiPU","; i#eta; i#phi; First Arrival Time [ns]", 171, -85,86,360,1,361);
  TH2F *jetEcalMeanTimeIEtaIPhiPU = new TH2F("jetEcalMeanTimeIEtaIPhiPU","; i#eta; i#phi; First Arrival Time [ns]", 171, -85,86,360,1,361);
  TH2F *jetEcalFirstArrivalTimeIXIYPU = new TH2F("jetEcalFirstArrivalTimeIXIYPU","; ix; iy; First Arrival Time [ns]", 100, 1, 101, 100, 1, 101);
  TH2F *jetEcalMeanTimeIXIYPU = new TH2F("jetEcalMeanTimeIXIYPU","; ix; iy; First Arrival Time [ns]", 100, 1, 101, 100, 1, 101);
  TH1F *jetEcalMeanTimeCryPU = new TH1F("jetEcalMeanTimeCryPU", "Mean time [ns]", 50, 8, 16);
  TH1F *jetEcalFirstTimeCryPU = new TH1F("jetEcalFirstTimeCryPU", "First hit time [ns]", 50, 8, 16);
  TH2F *jetEcalMeanVsFirstTimeCryPU = new TH2F("jetEcalMeanVsFirstTimeCryPU", "Mean time vs First hit time [ns]", 50, 8, 16, 50, 8, 16);
  for (UInt_t b=0; int(b)<jetEcalFirstArrivalTimeIEtaIPhiPU->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalFirstArrivalTimeIEtaIPhiPU->GetYaxis()->GetNbins()+2; ++c) {

      double firstArrivalTime = 0;
      double energySum = 0;
      double energyWeightedTime = 0;
      for (UInt_t d=0; int(d)<jetEcal_IEtaIPhiTimePU->GetZaxis()->GetNbins()+2; ++d) {
        if (jetEcal_IEtaIPhiTimePU->GetBinContent(b,c,d) > 0) {
          firstArrivalTime = jetEcal_IEtaIPhiTimePU->GetZaxis()->GetBinCenter(d);
          break;
        }
      }

      for (UInt_t d=0; int(d)<jetEcal_IEtaIPhiTimePU->GetZaxis()->GetNbins()+2; ++d) {
        energySum += jetEcal_IEtaIPhiTimePU->GetBinContent(b,c,d);
        energyWeightedTime += jetEcal_IEtaIPhiTimePU->GetZaxis()->GetBinCenter(d) * jetEcal_IEtaIPhiTimePU->GetBinContent(b,c,d);
      }

      if (energySum > 0) jetEcalTimeEnergyByCrystalPU->Fill(energyWeightedTime/energySum, energySum);
      if (energySum > 50) {
        jetEcalFirstArrivalTimeIEtaIPhiPU->SetBinContent(b,c,firstArrivalTime);
        jetEcalFirstTimeCryPU->Fill(firstArrivalTime);
        jetEcalMeanTimeIEtaIPhiPU->SetBinContent(b,c,energyWeightedTime/energySum);
        jetEcalMeanTimeCryPU->Fill(energyWeightedTime/energySum);
        jetEcalMeanVsFirstTimeCryPU->Fill(energyWeightedTime/energySum, firstArrivalTime);
      } else {
        jetEcalFirstArrivalTimeIEtaIPhiPU->SetBinContent(b,c,0);
        jetEcalMeanTimeIEtaIPhiPU->SetBinContent(b,c,0);
      }
    }
  }
  //endcaps
  for (UInt_t b=0; int(b)<jetEcalFirstArrivalTimeIXIYPU->GetXaxis()->GetNbins()+2; ++b) {
    for (UInt_t c=0; int(c)<jetEcalFirstArrivalTimeIXIYPU->GetYaxis()->GetNbins()+2; ++c) {

      double firstArrivalTime = 0;
      double energySum = 0;
      double energyWeightedTime = 0;
      for (UInt_t d=0; int(d)<jetEcal_IXIYTimePU->GetZaxis()->GetNbins()+2; ++d) {
        if (jetEcal_IXIYTimePU->GetBinContent(b,c,d) > 0) {
          firstArrivalTime = jetEcal_IXIYTimePU->GetZaxis()->GetBinCenter(d);
          break;
        }
      }

      for (UInt_t d=0; int(d)<jetEcal_IXIYTimePU->GetZaxis()->GetNbins()+2; ++d) {
        energySum += jetEcal_IXIYTimePU->GetBinContent(b,c,d);
        energyWeightedTime += jetEcal_IXIYTimePU->GetZaxis()->GetBinCenter(d) * jetEcal_IXIYTimePU->GetBinContent(b,c,d);
      }

      if (energySum > 0) jetEcalTimeEnergyByCrystalPU->Fill(energyWeightedTime/energySum, energySum);
      if (energySum > 50) {
        jetEcalFirstArrivalTimeIXIYPU->SetBinContent(b,c,firstArrivalTime);
        jetEcalFirstTimeCryPU->Fill(firstArrivalTime);
        jetEcalMeanTimeIXIYPU->SetBinContent(b,c,energyWeightedTime/energySum);
        jetEcalMeanTimeCryPU->Fill(energyWeightedTime/energySum);
        jetEcalMeanVsFirstTimeCryPU->Fill(energyWeightedTime/energySum, firstArrivalTime);
      } else {
        jetEcalFirstArrivalTimeIXIYPU->SetBinContent(b,c,0);
        jetEcalMeanTimeIXIYPU->SetBinContent(b,c,0);
      }
    }
  }

  TCanvas *cv = new TCanvas ("cv","cv",800,600);

  cv->SetRightMargin(0.20);
  jetEcalFirstArrivalTimeEtaPhi->SetMinimum( 8);
  jetEcalFirstArrivalTimeEtaPhi->SetMaximum(16);
  jetEcalFirstArrivalTimeEtaPhi->Draw("colz");
  cv->SaveAs("jetEcalFirstArrivalTimeEtaPhi_signalonly.gif");

  cv->SetRightMargin(0.20);
  jetEcalMeanTimeEtaPhi->SetMinimum( 8);
  jetEcalMeanTimeEtaPhi->SetMaximum(16);
  jetEcalMeanTimeEtaPhi->Draw("colz");
  cv->SaveAs("jetEcalMeanTimeEtaPhi_signalonly.gif");

  cv->SetRightMargin(0.20);
  jetEcalFirstArrivalTimeEtaPhiPU->SetMinimum( 8);
  jetEcalFirstArrivalTimeEtaPhiPU->SetMaximum(16);
  jetEcalFirstArrivalTimeEtaPhiPU->Draw("colz");
  cv->SaveAs("jetEcalFirstArrivalTimeEtaPhi_signalpluspileup.gif");

  cv->SetRightMargin(0.20);
  jetEcalMeanTimeEtaPhiPU->SetMinimum( 8);
  jetEcalMeanTimeEtaPhiPU->SetMaximum(16);
  jetEcalMeanTimeEtaPhiPU->Draw("colz");
  cv->SaveAs("jetEcalMeanTimeEtaPhi_signalpluspileup.gif");

  delete cv;
  file->WriteTObject(jetEcalTime, jetEcalTime->GetName(), "Write");
  delete jetEcalTime;
  file->WriteTObject(jetEcalEnergyWeightedTime, jetEcalEnergyWeightedTime->GetName(), "Write");
  delete jetEcalEnergyWeightedTime;
  file->WriteTObject(jetEcalEnergyWeightedTimeEtaPhi,jetEcalEnergyWeightedTimeEtaPhi->GetName(), "Write");
  delete jetEcalEnergyWeightedTimeEtaPhi;
  file->WriteTObject(jetEcalEnergyWeightedTimeIEtaIPhi,jetEcalEnergyWeightedTimeIEtaIPhi->GetName(), "Write");
  delete jetEcalEnergyWeightedTimeIEtaIPhi;
  file->WriteTObject(jetEcalEnergyWeightedTimeIXIY,jetEcalEnergyWeightedTimeIXIY->GetName(), "Write");
  delete jetEcalEnergyWeightedTimeIXIY;
  file->WriteTObject(jetEcalEnergyIEtaIPhi,jetEcalEnergyIEtaIPhi->GetName(), "Write");
  delete jetEcalEnergyIEtaIPhi;
  file->WriteTObject(jetEcalEnergyIXIY,jetEcalEnergyIXIY->GetName(), "Write");
  delete jetEcalEnergyIXIY;
  file->WriteTObject(jetEcalCountIEtaIPhi,jetEcalCountIEtaIPhi->GetName(), "Write");
  delete jetEcalCountIEtaIPhi;
  file->WriteTObject(jetEcalCountIXIY,jetEcalCountIXIY->GetName(), "Write");
  delete jetEcalCountIXIY;
  file->WriteTObject(jetEcalTimeEtaPhi,jetEcalTimeEtaPhi->GetName(), "Write");
  delete jetEcalTimeEtaPhi;
  file->WriteTObject(jetEcalEnergyEtaPhi,jetEcalEnergyEtaPhi->GetName(), "Write");
  delete jetEcalEnergyEtaPhi;
  file->WriteTObject(jetParticleEnergyEtaPhi,jetParticleEnergyEtaPhi->GetName(), "Write");
  delete jetParticleEnergyEtaPhi;
  file->WriteTObject(jetEcalTimeEnergy, jetEcalTimeEnergy->GetName(), "Write");
  delete jetEcalTimeEnergy;
  file->WriteTObject(jetEcal_EtaPhiTime,jetEcal_EtaPhiTime->GetName(),"Write");
  delete jetEcal_EtaPhiTime;
//  file->WriteTObject(jetEcal_IEtaIPhiTime,jetEcal_IEtaIPhiTime->GetName(),"Write");
  delete jetEcal_IEtaIPhiTime;
//  file->WriteTObject(jetEcal_IXIYTime,jetEcal_IXIYTime->GetName(),"Write");
  delete jetEcal_IXIYTime;
  file->WriteTObject(jetEcalFirstArrivalTimeEtaPhi,jetEcalFirstArrivalTimeEtaPhi->GetName(),"Write");
  delete jetEcalFirstArrivalTimeEtaPhi;
  file->WriteTObject(jetEcalFirstArrivalTimeIEtaIPhi,jetEcalFirstArrivalTimeIEtaIPhi->GetName(),"Write");
  delete jetEcalFirstArrivalTimeIEtaIPhi;
  file->WriteTObject(jetEcalFirstArrivalTimeIXIY,jetEcalFirstArrivalTimeIXIY->GetName(),"Write");
  delete jetEcalFirstArrivalTimeIXIY;
  file->WriteTObject(jetEcalMeanTimeEtaPhi,jetEcalMeanTimeEtaPhi->GetName(),"Write");
  delete jetEcalMeanTimeEtaPhi;
  file->WriteTObject(jetEcalMeanTimeIEtaIPhi,jetEcalMeanTimeIEtaIPhi->GetName(),"Write");
  delete jetEcalMeanTimeIEtaIPhi;
  file->WriteTObject(jetEcalMeanTimeIXIY,jetEcalMeanTimeIXIY->GetName(),"Write");
  delete jetEcalMeanTimeIXIY;
  file->WriteTObject(jetEcalMeanTime, jetEcalMeanTime->GetName(), "Write");
  delete jetEcalCountEtaPhi;
  delete jetEcalMeanTime;
  file->WriteTObject(jetEcalMeanTimeCry, jetEcalMeanTimeCry->GetName(), "Write");
  delete jetEcalMeanTimeCry;
  file->WriteTObject(jetEcalFirstTime, jetEcalFirstTime->GetName(), "Write");
  delete jetEcalFirstTime;
  file->WriteTObject(jetEcalFirstTimeCry, jetEcalFirstTimeCry->GetName(), "Write");
  delete jetEcalFirstTimeCry;
  file->WriteTObject(jetEcalMeanVsFirstTime, jetEcalMeanVsFirstTime->GetName(), "Write");
  delete jetEcalMeanVsFirstTime;
  file->WriteTObject(jetEcalMeanVsFirstTimeCry, jetEcalMeanVsFirstTimeCry->GetName(), "Write");
  delete jetEcalMeanVsFirstTimeCry;
  file->WriteTObject(jetMeanTimeFirstTimeCumulative, jetMeanTimeFirstTimeCumulative->GetName(), "Write");
  delete jetMeanTimeFirstTimeCumulative;
  file->WriteTObject(jetEcalTimeEnergyByCrystal, jetEcalTimeEnergyByCrystal->GetName(), "Write");
  delete jetEcalTimeEnergyByCrystal;
  file->WriteTObject(jetEnergyChangeDueToPileup, jetEnergyChangeDueToPileup->GetName(), "Write");
  delete jetEnergyChangeDueToPileup;
  file->WriteTObject(jetEcalTimePU, jetEcalTimePU->GetName(), "Write");
  delete jetEcalTimePU;
  file->WriteTObject(jetEcalEnergyWeightedTimePU, jetEcalEnergyWeightedTimePU->GetName(), "Write");
  delete jetEcalEnergyWeightedTimePU;
  file->WriteTObject(jetEcalEnergyWeightedTimeEtaPhiPU,jetEcalEnergyWeightedTimeEtaPhiPU->GetName(), "Write");
  delete jetEcalEnergyWeightedTimeEtaPhiPU;
  file->WriteTObject(jetEcalEnergyWeightedTimeIXIYPU,jetEcalEnergyWeightedTimeIXIYPU->GetName(), "Write");
  delete jetEcalEnergyWeightedTimeIXIYPU;
  file->WriteTObject(jetEcalEnergyWeightedTimeIEtaIPhiPU,jetEcalEnergyWeightedTimeIEtaIPhiPU->GetName(), "Write");
  delete jetEcalEnergyWeightedTimeIEtaIPhiPU;
  file->WriteTObject(jetEcalEnergyIEtaIPhiPU,jetEcalEnergyIEtaIPhiPU->GetName(), "Write");
  delete jetEcalEnergyIEtaIPhiPU;
  file->WriteTObject(jetEcalEnergyIXIYPU,jetEcalEnergyIXIYPU->GetName(), "Write");
  delete jetEcalEnergyIXIYPU;
  file->WriteTObject(jetEcalCountIEtaIPhiPU,jetEcalCountIEtaIPhiPU->GetName(), "Write");
  delete jetEcalCountIEtaIPhiPU;
  file->WriteTObject(jetEcalCountIXIYPU,jetEcalCountIXIYPU->GetName(), "Write");
  delete jetEcalCountIXIYPU;
  file->WriteTObject(jetEcalTimeEtaPhiPU,jetEcalTimeEtaPhiPU->GetName(), "Write");
  delete jetEcalTimeEtaPhiPU;
  file->WriteTObject(jetEcalEnergyEtaPhiPU,jetEcalEnergyEtaPhiPU->GetName(), "Write");
  delete jetEcalEnergyEtaPhiPU;
  file->WriteTObject(jetParticleEnergyEtaPhiPU,jetParticleEnergyEtaPhiPU->GetName(), "Write");
  delete jetParticleEnergyEtaPhiPU;
  file->WriteTObject(jetEcalTimeEnergyPU, jetEcalTimeEnergyPU->GetName(), "Write");
  delete jetEcalTimeEnergyPU;
  file->WriteTObject(jetEcal_EtaPhiTimePU,jetEcal_EtaPhiTimePU->GetName(),"Write");
  delete jetEcal_EtaPhiTimePU;
//  file->WriteTObject(jetEcal_IEtaIPhiTimePU,jetEcal_IEtaIPhiTimePU->GetName(),"Write");
  delete jetEcal_IEtaIPhiTimePU;
//  file->WriteTObject(jetEcal_IXIYTimePU,jetEcal_IXIYTimePU->GetName(),"Write");
  delete jetEcal_IXIYTimePU;
  file->WriteTObject(jetEcalFirstArrivalTimeEtaPhiPU,jetEcalFirstArrivalTimeEtaPhiPU->GetName(),"Write");
  delete jetEcalFirstArrivalTimeEtaPhiPU;
  file->WriteTObject(jetEcalFirstArrivalTimeIEtaIPhiPU,jetEcalFirstArrivalTimeIEtaIPhiPU->GetName(),"Write");
  delete jetEcalFirstArrivalTimeIEtaIPhiPU;
  file->WriteTObject(jetEcalFirstArrivalTimeIXIYPU,jetEcalFirstArrivalTimeIXIYPU->GetName(),"Write");
  delete jetEcalFirstArrivalTimeIXIYPU;
  file->WriteTObject(jetEcalMeanTimeEtaPhiPU,jetEcalMeanTimeEtaPhiPU->GetName(),"Write");
  delete jetEcalMeanTimeEtaPhiPU;
  file->WriteTObject(jetEcalMeanTimeIEtaIPhiPU,jetEcalMeanTimeIEtaIPhiPU->GetName(),"Write");
  delete jetEcalMeanTimeIEtaIPhiPU;
  file->WriteTObject(jetEcalMeanTimeIXIYPU,jetEcalMeanTimeIXIYPU->GetName(),"Write");
  delete jetEcalMeanTimeIXIYPU;
  file->WriteTObject(jetEcalMeanTimePU, jetEcalMeanTimePU->GetName(), "Write");
  delete jetEcalMeanTimePU;
  file->WriteTObject(jetEcalMeanTimeCryPU, jetEcalMeanTimeCryPU->GetName(), "Write");
  delete jetEcalMeanTimeCryPU;
  file->WriteTObject(jetEcalFirstTimePU, jetEcalFirstTimePU->GetName(), "Write");
  delete jetEcalFirstTimePU;
  file->WriteTObject(jetEcalFirstTimeCryPU, jetEcalFirstTimeCryPU->GetName(), "Write");
  delete jetEcalFirstTimeCryPU;
  delete jetEcalCountEtaPhiPU;
  file->WriteTObject(jetEcalMeanVsFirstTimePU, jetEcalMeanVsFirstTimePU->GetName(), "Write");
  delete jetEcalMeanVsFirstTimePU;
  file->WriteTObject(jetEcalMeanVsFirstTimeCryPU, jetEcalMeanVsFirstTimeCryPU->GetName(), "Write");
  delete jetEcalMeanVsFirstTimeCryPU;
  file->WriteTObject(jetMeanTimeFirstTimeCumulativePU, jetMeanTimeFirstTimeCumulativePU->GetName(), "Write");
  delete jetMeanTimeFirstTimeCumulativePU;
  file->WriteTObject(jetEcalTimeEnergyByCrystalPU, jetEcalTimeEnergyByCrystalPU->GetName(), "Write");
  delete jetEcalTimeEnergyByCrystalPU;
  file->WriteTObject(numHitsAbove4GeV, numHitsAbove4GeV->GetName(), "Write");
  delete numHitsAbove4GeV;
  file->WriteTObject(numHitsAbove4GeVPU, numHitsAbove4GeVPU->GetName(), "Write");
  delete numHitsAbove4GeVPU;
  file->Close();
  delete file;       
  infile->Close();
  delete infile;
  infilePU->Close();
  delete infilePU;
  fSignal->Close();
  delete fSignal;
  fPU->Close();
  delete fPU;
  delete MyRandom;
  delete info;
  delete infoPU;
  delete genparticleArr;
  delete genparticleArrPU;
  delete genjetArr;
  delete genjetArrPU;
  return;


}

void JetEnergyResolutionStudy() {
  
    gROOT->SetBatch();
    processSignalPlusPileup();
    return;
}
