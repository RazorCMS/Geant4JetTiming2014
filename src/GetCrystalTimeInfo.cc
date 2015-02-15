//#ifdef __MAKECINT__
//#pragma link C++ class vector<vector<float> >+;
//#pragma link C++ class vector<vector<int> >+;
//#endif

#include <iostream>
#include <fstream>
#include <map>
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
#include <TChain.h>
#include <sstream>

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
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
using namespace fastjet;

using namespace std;
const double weightExp = 2;
class CrystTimeInfo{
    
    private:
    TChain *baconChain;
    TChain *geantChain;
    string baconFileList;
    string g4FileList;
    int numFiles;

    public:
    CrystTimeInfo(string, string);
    ~CrystTimeInfo() {delete baconChain; delete geantChain;}
    void GetChains();
    void GetCrystalTimeInfo();
};

CrystTimeInfo::CrystTimeInfo(string baconList, string g4List){
    baconFileList = baconList;
    g4FileList = g4List;
    numFiles = 0;
    GetChains();
}

void CrystTimeInfo::GetChains(){

    cout << "Chaining event files..." << endl;
    baconChain = new TChain("Events");
    ifstream fin(baconFileList.c_str());
    string tmp;
    while(getline(fin, tmp)){
        baconChain->Add(tmp.c_str());
        cout << "   Chaining " << tmp << endl;
        numFiles++;
    }
    fin.close();

    cout << "Chaining Geant hit files..." << endl;
    geantChain = new TChain("G4SIM");
    ifstream fin2(g4FileList.c_str());
    while(getline(fin2, tmp)){
            geantChain->Add(tmp.c_str());
            cout << "   Chaining " << tmp << endl;
    }       
    fin2.close();
}

void CrystTimeInfo::GetCrystalTimeInfo() {
  //gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  const uint mixNPU = 140;
  bool mixPU = false;

 //*****************************************************************************************
  // Set up input for analysis tree
  //*****************************************************************************************

  // Data structures to store info from TTrees
  cmsana::TEventInfo *info  = new cmsana::TEventInfo();
  TClonesArray *genparticleArr = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArr = new TClonesArray("cmsana::TGenJet");

  cout << "Setting event chain branches..." << endl;
  baconChain->SetBranchAddress("Info",     &info);        TBranch *infoBr     = baconChain->GetBranch("Info");
  baconChain->SetBranchAddress("GenParticle", &genparticleArr); TBranch *genparticleBr = baconChain->GetBranch("GenParticle");
  baconChain->SetBranchAddress("GenJet", &genjetArr); TBranch *genjetBr = baconChain->GetBranch("GenJet");

  /*infilePU->cd();
  cmsana::TEventInfo *infoPU  = new cmsana::TEventInfo();
  TClonesArray *genparticleArrPU = new TClonesArray("cmsana::TGenParticle");
  TClonesArray *genjetArrPU = new TClonesArray("cmsana::TGenJet");
*/
/*  baconChainPU = (TTree*)infilePU->Get("Events"); assert(baconChainPU);
  baconChainPU->SetBranchAddress("Info",     &infoPU);        TBranch *infoBrPU     = baconChainPU->GetBranch("Info");
  baconChainPU->SetBranchAddress("GenParticle", &genparticleArrPU); TBranch *genparticleBrPU = baconChainPU->GetBranch("GenParticle");
  baconChainPU->SetBranchAddress("GenJet", &genjetArrPU); TBranch *genjetBrPU = baconChainPU->GetBranch("GenJet");
  */

  //*****************************************************************************************
  // Set up geant hits tree
  //*****************************************************************************************
  //geantChain->SetMakeClass(1);

  // set up the Signal tree variables
  int numCrystsB= 0;
  int numCrystsE = 0;
  int ieta[61200];
  int iphi[61200];
  int ix[14648];
  int iy[14648];
  int iz[14648];
  float tming4EB[61200];
  float xtming4EB[61200];
  float ytming4EB[61200];
  float ztming4EB[61200];
  float tming4EE[14648];
  float xtming4EE[14648];
  float ytming4EE[14648];
  float ztming4EE[14648];

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

  geantChain->SetBranchAddress("eg4EB",&eg4EBSignal);
  geantChain->SetBranchAddress("tg4EB",&tg4EBSignal);
  geantChain->SetBranchAddress("prexg4EB",&prexg4EBSignal);
  geantChain->SetBranchAddress("preyg4EB",&preyg4EBSignal);
  geantChain->SetBranchAddress("prezg4EB",&prezg4EBSignal);
  geantChain->SetBranchAddress("eg4EE",&eg4EESignal);
  geantChain->SetBranchAddress("tg4EE",&tg4EESignal);
  geantChain->SetBranchAddress("prexg4EE",&prexg4EESignal);
  geantChain->SetBranchAddress("preyg4EE",&preyg4EESignal);
  geantChain->SetBranchAddress("prezg4EE",&prezg4EESignal);
  geantChain->SetBranchAddress("ietag4EB", ieta);
  geantChain->SetBranchAddress("iphig4EB", iphi);
  geantChain->SetBranchAddress("ixg4EE", ix);
  geantChain->SetBranchAddress("iyg4EE", iy);
  geantChain->SetBranchAddress("izg4EE", iz);
  geantChain->SetBranchAddress("tming4EB", tming4EB);
  geantChain->SetBranchAddress("xtming4EB", xtming4EB);
  geantChain->SetBranchAddress("ytming4EB", ytming4EB);
  geantChain->SetBranchAddress("ztming4EB", ztming4EB);
  geantChain->SetBranchAddress("tming4EE", tming4EE);
  geantChain->SetBranchAddress("xtming4EE", xtming4EE);
  geantChain->SetBranchAddress("ytming4EE", ytming4EE);
  geantChain->SetBranchAddress("ztming4EE", ztming4EE);
  geantChain->SetBranchAddress("ng4EB", &numCrystsB);
  geantChain->SetBranchAddress("ng4EE", &numCrystsE);
  cout << "Done setting chain branches.  Real-life chains don't have branches, but whatever." << endl;

  //*****************************************************************************************
  // Set up geant hits tree
  //*****************************************************************************************
/*  TFile *fPU = new TFile("/afs/cern.ch/user/s/sixie/work/public/Phase2Upgrade/Timing/Simulation/MinBias/g4simhitsEcal_2.root");
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
*/
  double crystEThreshold = 1.0; //minimum energy to store a crystal in the output tree

  //output file
  cout << "Creating output tree..." << endl;
  
  baconFileList.erase(baconFileList.length() - 5, baconFileList.length());
  string fileListName = baconFileList.erase(0, baconFileList.length()- 7);
  TFile *fnew = new TFile(Form("/home/djanders/CMSSW_5_3_9/src/CrystalTimeInfoFiles/Hgg/out/CrystalTimeInfo%s.root", fileListName.c_str()), "recreate");
  TTree *ecalTree = new TTree("ecalTree", "Jet ECAL info, including timing");

  //tree variables
  double vertexZ = 0;
  vector<double> *jEcalE;
  vector<double> *jEta;
  vector<double> *jPhi;
  vector<double> *jPt;
  vector< vector<int> > *crystIEta;
  vector< vector<int> > *crystIPhi;
  vector< vector<double> > *crystEB;
  vector< vector<double> > *crystTB;
  vector< vector<double> > *crystE3x3B;
  vector< vector<double> > *crystT3x3B;
  vector< vector<double> > *crystE5x5B;
  vector< vector<double> > *crystT5x5B;
  vector< vector<int> > *crystIX;
  vector< vector<int> > *crystIY;
  vector< vector<int> > *crystIZ;
  vector< vector<double> > *crystEE;
  vector< vector<double> > *crystTE;
  vector< vector<double> > *crystE3x3E;
  vector< vector<double> > *crystT3x3E;
  vector< vector<double> > *crystE5x5E;
  vector< vector<double> > *crystT5x5E;

  vector< vector<double> > *rhoB;
  vector< vector<double> > *rhoE;
  vector< vector<double> > *hitZB;
  vector< vector<double> > *hitZE;

  //new variables
  vector<double> *photEta;
  vector<double> *photPhi;
  vector<double> *photPt;
  vector< vector<int> > *photIEta;
  vector< vector<int> > *photIPhi;
  vector< vector<double> > *photEB;
  vector< vector<double> > *photE3x3B;
  vector< vector<double> > *photTB;
  vector< vector<double> > *photT3x3B;
  vector< vector<int> > *photIX;
  vector< vector<int> > *photIY;
  vector< vector<int> > *photIZ;
  vector< vector<double> > *photEE;
  vector< vector<double> > *photE3x3E;
  vector< vector<double> > *photTE;
  vector< vector<double> > *photT3x3E;
  vector< vector<double> > *photRhoB;
  vector< vector<double> > *photRhoE;
  vector< vector<double> > *photZB;
  vector< vector<double> > *photZE;
  vector< vector<double> > *photFirstTimeB;
  vector< vector<double> > *photFirstTimeE;
  vector< vector<double> > *photFirstRhoB;
  vector< vector<double> > *photFirstRhoE;
  vector< vector<double> > *photFirstZB;
  vector< vector<double> > *photFirstZE;

  vector< vector<double> > *firstTimeB;
  vector< vector<double> > *firstTimeE;
  vector< vector<double> > *firstHitRhoB;
  vector< vector<double> > *firstHitZB;
  vector< vector<double> > *firstHitRhoE;
  vector< vector<double> > *firstHitZE;
  
  jEcalE = new std::vector<double>; jEcalE->clear();
  jEta = new std::vector<double>; jEta->clear();
  jPhi = new std::vector<double>; jPhi->clear();
  jPt = new std::vector<double>; jPt->clear();

  crystEB = new std::vector<std::vector<double> >; crystEB->clear();
  crystTB = new std::vector<std::vector<double> >; crystTB->clear();
  crystE3x3B = new std::vector<std::vector<double> >; crystE3x3B->clear();
  crystT3x3B = new std::vector<std::vector<double> >; crystT3x3B->clear();
  crystE5x5B = new std::vector<std::vector<double> >; crystE5x5B->clear();
  crystT5x5B = new std::vector<std::vector<double> >; crystT5x5B->clear();
  crystEE = new std::vector<std::vector<double> >; crystEE->clear();
  crystTE = new std::vector<std::vector<double> >; crystTE->clear();
  crystE3x3E = new std::vector<std::vector<double> >; crystE3x3E->clear();
  crystT3x3E = new std::vector<std::vector<double> >; crystT3x3E->clear();
  crystE5x5E = new std::vector<std::vector<double> >; crystE5x5E->clear();
  crystT5x5E = new std::vector<std::vector<double> >; crystT5x5E->clear();
  
  rhoB = new std::vector<std::vector<double> >; rhoB->clear();
  rhoE = new std::vector<std::vector<double> >; rhoE->clear();
  hitZB = new std::vector<std::vector<double> >; hitZB->clear();
  hitZE = new std::vector<std::vector<double> >; hitZE->clear();
  
  crystIEta = new std::vector<std::vector<int> >; crystIEta->clear();
  crystIPhi = new std::vector<std::vector<int> >; crystIPhi->clear();
  crystIX = new std::vector<std::vector<int> >; crystIX->clear();
  crystIY = new std::vector<std::vector<int> >; crystIY->clear();
  crystIZ = new std::vector<std::vector<int> >; crystIZ->clear();

  firstTimeB = new std::vector< vector<double> >; firstTimeB->clear();
  firstTimeE = new std::vector< vector<double> >; firstTimeE->clear();
  firstHitRhoB = new std::vector< vector<double> >; firstHitRhoB->clear();
  firstHitZB = new std::vector< vector<double> >; firstHitZB->clear();
  firstHitRhoE = new std::vector< vector<double> >; firstHitRhoE->clear();
  firstHitZE = new std::vector< vector<double> >; firstHitZE->clear();

  photEta = new std::vector<double>; photEta->clear();
  photPhi = new std::vector<double>; photPhi->clear();
  photPt = new std::vector<double>; photPt->clear();
  photIEta = new vector< vector<int> >; photIEta->clear();
  photIPhi = new vector< vector<int> >; photIPhi->clear();
  photIX = new vector< vector<int> >; photIX->clear();
  photIY = new vector< vector<int> >; photIY->clear();
  photIZ = new vector< vector<int> >; photIZ->clear();
  photEB = new vector< vector<double> >; photEB->clear();
  photE3x3B = new vector< vector<double> >; photE3x3B->clear();
  photTB = new vector< vector<double> >; photTB->clear();
  photT3x3B = new vector< vector<double> >; photT3x3B->clear();
  photEE = new vector< vector<double> >; photEE->clear();
  photE3x3E = new vector< vector<double> >; photE3x3E->clear();
  photTE = new vector< vector<double> >; photTE->clear();
  photT3x3E = new vector< vector<double> >; photT3x3E->clear();
  photRhoB = new vector< vector<double> >; photRhoB->clear();
  photRhoE = new vector< vector<double> >; photRhoE->clear();
  photZB = new vector< vector<double> >; photZB->clear();
  photZE = new vector< vector<double> >; photZE->clear();
  photFirstTimeB = new vector< vector<double> >; photFirstTimeB->clear();
  photFirstTimeE = new vector< vector<double> >; photFirstTimeE->clear();
  photFirstRhoB = new vector< vector<double> >; photFirstRhoB->clear();
  photFirstRhoE = new vector< vector<double> >; photFirstRhoE->clear();
  photFirstZB = new vector< vector<double> >; photFirstZB->clear();
  photFirstZE = new vector< vector<double> >; photFirstZE->clear();

  ecalTree->Branch("vertexZ", &vertexZ, "vertexZ/D");
  ecalTree->Branch("jEcalE", "std::vector<double>", &jEcalE);
  ecalTree->Branch("jEta", "std::vector<double>", &jEta);
  ecalTree->Branch("jPhi", "std::vector<double>", &jPhi);
  ecalTree->Branch("jPt", "std::vector<double>", &jPt);

  ecalTree->Branch("crystIEta", "std::vector< std::vector<int> >", &crystIEta);
  ecalTree->Branch("crystIPhi", "std::vector< std::vector<int> >", &crystIPhi);
  ecalTree->Branch("crystEB", "std::vector<std::vector<double> >", &crystEB);
  ecalTree->Branch("crystTB", "std::vector<std::vector<double> >", &crystTB);
  ecalTree->Branch("crystE3x3B", "std::vector<std::vector<double> >", &crystE3x3B);
  ecalTree->Branch("crystT3x3B", "std::vector<std::vector<double> >", &crystT3x3B);
  ecalTree->Branch("crystE5x5B", "std::vector<std::vector<double> >", &crystE5x5B);
  ecalTree->Branch("crystT5x5B", "std::vector<std::vector<double> >", &crystT5x5B);

  ecalTree->Branch("crystIX", "std::vector< std::vector<int> >", &crystIX);
  ecalTree->Branch("crystIY", "std::vector< std::vector<int> >", &crystIY);
  ecalTree->Branch("crystIZ", "std::vector< std::vector<int> >", &crystIZ);
  ecalTree->Branch("crystEE", "std::vector<std::vector<double> >", &crystEE);
  ecalTree->Branch("crystTE", "std::vector<std::vector<double> >", &crystTE);
  ecalTree->Branch("crystE3x3E", "std::vector<std::vector<double> >", &crystE3x3E);
  ecalTree->Branch("crystT3x3E", "std::vector<std::vector<double> >", &crystT3x3E);
  ecalTree->Branch("crystE5x5E", "std::vector<std::vector<double> >", &crystE5x5E);
  ecalTree->Branch("crystT5x5E", "std::vector<std::vector<double> >", &crystT5x5E);
  
  ecalTree->Branch("rhoB", "std::vector<std::vector<double> >", &rhoB);
  ecalTree->Branch("rhoE", "std::vector<std::vector<double> >", &rhoE);
  ecalTree->Branch("hitZB", "std::vector<std::vector<double> >", &hitZB);
  ecalTree->Branch("hitZE", "std::vector<std::vector<double> >", &hitZE);  

  ecalTree->Branch("firstTimeB", "std::vector< std::vector<double> >", &firstTimeB);
  ecalTree->Branch("firstTimeE", "std::vector< std::vector<double> >", &firstTimeE);
  ecalTree->Branch("firstHitRhoB", "std::vector< std::vector<double> >", &firstHitRhoB);
  ecalTree->Branch("firstHitRhoE", "std::vector< std::vector<double> >", &firstHitRhoE);
  ecalTree->Branch("firstHitZB", "std::vector< std::vector<double> >", &firstHitZB);
  ecalTree->Branch("firstHitZE", "std::vector< std::vector<double> >", &firstHitZE);

  ecalTree->Branch("photEta", "std::vector<double>", &photEta);
  ecalTree->Branch("photPhi", "std::vector<double>", &photPhi);
  ecalTree->Branch("photPt", "std::vector<double>", &photPt);
  ecalTree->Branch("photIEta", "std::vector< std::vector<int> >", &photIEta);
  ecalTree->Branch("photIPhi", "std::vector< std::vector<int> >", &photIPhi);
  ecalTree->Branch("photIX", "std::vector< std::vector<int> >", &photIX);
  ecalTree->Branch("photIY", "std::vector< std::vector<int> >", &photIY);
  ecalTree->Branch("photIZ", "std::vector< std::vector<int> >", &photIZ);
  ecalTree->Branch("photEB", "std::vector< std::vector<double> >", &photEB);
  ecalTree->Branch("photE3x3B", "std::vector< std::vector<double> >", &photE3x3B);
  ecalTree->Branch("photTB", "std::vector< std::vector<double> >", &photTB);
  ecalTree->Branch("photT3x3B", "std::vector< std::vector<double> >", &photT3x3B);
  ecalTree->Branch("photEE", "std::vector< std::vector<double> >", &photEE);
  ecalTree->Branch("photE3x3E", "std::vector< std::vector<double> >", &photE3x3E);
  ecalTree->Branch("photTE", "std::vector< std::vector<double> >", &photTE);
  ecalTree->Branch("photT3x3E", "std::vector< std::vector<double> >", &photT3x3E);
  ecalTree->Branch("photRhoB", "std::vector< std::vector<double> >", &photRhoB);
  ecalTree->Branch("photRhoE", "std::vector< std::vector<double> >", &photRhoE);
  ecalTree->Branch("photZB", "std::vector< std::vector<double> >", &photZB);
  ecalTree->Branch("photZE", "std::vector< std::vector<double> >", &photZE);
  ecalTree->Branch("photFirstTimeB", "std::vector< std::vector<double> >", &photFirstTimeB);
  ecalTree->Branch("photFirstTimeE", "std::vector< std::vector<double> >", &photFirstTimeE);
  ecalTree->Branch("photFirstRhoB", "std::vector< std::vector<double> >", &photFirstRhoB);
  ecalTree->Branch("photFirstRhoE", "std::vector< std::vector<double> >", &photFirstRhoE);
  ecalTree->Branch("photFirstZB", "std::vector< std::vector<double> >", &photFirstZB);
  ecalTree->Branch("photFirstZE", "std::vector< std::vector<double> >", &photFirstZE);

  //LOOP OVER EVENTS
  cout << "Looping over events..." << endl;
  //numFiles*10 is the total number of events
  for(UInt_t ientry=0; ientry < numFiles*10; ientry++) { 
    
    //clear tree variables
    firstTimeB->clear();
    firstTimeE->clear();
    firstHitRhoB->clear();
    firstHitRhoE->clear();
    firstHitZB->clear();
    firstHitZE->clear();

    photEta->clear();
    photPhi->clear();
    photPt->clear();
    photIEta->clear();
    photIPhi->clear();
    photIX->clear();
    photIY->clear();
    photIZ->clear();
    photEB->clear();
    photE3x3B->clear();
    photTB->clear();
    photT3x3B->clear();
    photEE->clear();
    photE3x3E->clear();
    photTE->clear();
    photT3x3E->clear();
    photRhoB->clear();
    photRhoE->clear();
    photZB->clear();
    photZE->clear();
    photFirstTimeB->clear();
    photFirstTimeE->clear();
    photFirstRhoB->clear();
    photFirstRhoE->clear();
    photFirstZB->clear();
    photFirstZE->clear();

    jEcalE->clear();
    jEta->clear();
    jPhi->clear();
    jPt->clear();

    crystIEta->clear();
    crystIPhi->clear();
    crystEB->clear();
    crystTB->clear();
    crystE3x3B->clear();
    crystT3x3B->clear();
    crystE5x5B->clear();
    crystT5x5B->clear();
    
    crystIX->clear();
    crystIY->clear();
    crystIZ->clear();
    crystEE->clear();
    crystTE->clear();
    crystE3x3E->clear();
    crystT3x3E->clear();
    crystE5x5E->clear();
    crystT5x5E->clear();

    rhoB->clear();
    rhoE->clear();
    hitZB->clear();
    hitZE->clear();

    cout << "\n\nEvent " << ientry << "\n";

    //****************************************************
    //Generate list of minbias events to mix
    //****************************************************
/*    infilePU->cd();
    vector<uint> pileupEventsToMix;
    for (uint p=0; p < mixNPU; ++p) {
      uint eventToMix = MyRandom->Integer(baconChainPU->GetEntries());
      bool alreadyUsed = false;
      for (uint q=0; q<pileupEventsToMix.size(); ++q) {
        if (eventToMix == pileupEventsToMix[q]) alreadyUsed = true;
      }
      if (!alreadyUsed) pileupEventsToMix.push_back(eventToMix);
    }
*/
    vector<PseudoJet> VisibleParticles;
    VisibleParticles.clear();
    
    genparticleArr->Clear(); 
    baconChain->GetEntry(ientry);
    //genparticleBr->GetEntry(ientry);
    //infoBr->GetEntry(ientry);

    vertexZ = info->genVertexZ;
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
  
    //********************************************************
    // Mix Events
    //********************************************************
/*    if(mixPU){
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
    }
*/
    JetDefinition Definition(antikt_algorithm, 0.5);
    ClusterSequence Sequence(VisibleParticles, Definition);
    vector<PseudoJet> Jets = Sequence.inclusive_jets(0.5);
    //sort by decreasing pT
    Jets = fastjet::sorted_by_pt(Jets);
    //********************************************************
    // Print out Jets 
    //********************************************************
    cout << "GenJets " << Jets.size() << "\n";
    geantChain->GetEntry(ientry);
    for (int j=0; j< Jets.size()  ; ++j) {

      bool isPhoton = false;
      int photIndex = -1;
      for (uint p=0; p<genPhotons.size(); ++p) {
        if (cmsana::deltaR(genPhotons[p]->eta, genPhotons[p]->phi,Jets[j].eta(),Jets[j].phi()) < 0.3){
            isPhoton = true;
            photIndex = p;
        }
      }

      if (isPhoton && fabs(genPhotons[photIndex]->eta) < 3.0 ){
        
	//get the GEANT hits for the photon
        double photEcalEnergySum = 0;
        vector<double> thisPhotFirstTB;
        vector<double> thisPhotFirstRhoB;
        vector<double> thisPhotFirstZB;
     
          //store the crystal information in a map:
          //the key pair is ieta,iphi, and the pair of doubles is <energy, energy-weighted time>
          map<pair<int,int>, pair<double, double> > barrelMap;
          //second map to store position info -- the two doubles are rho and z
          map<pair<int,int>, pair<double, double> > barrelPosMap; 
          //initialize -- just check which crystals have hits in the jet
          for(Int_t l = 0; l < eg4EBSignal->size(); l++) {
              for(int m = 0; m < eg4EBSignal->at(l).size(); m++){  
                TVector3 vEBSignal;
                vEBSignal.SetXYZ((*prexg4EBSignal)[l][m], (*preyg4EBSignal)[l][m], (*prezg4EBSignal)[l][m] - vertexZ);
                if (cmsana::deltaR(vEBSignal.Eta(),vEBSignal.Phi(),genPhotons[photIndex]->eta,genPhotons[photIndex]->phi) < 0.3 && (*eg4EBSignal)[l][m] > 0) {
                    barrelMap[make_pair(ieta[l], iphi[l])] = make_pair(0.0, 0.0);
                    break; //found a hit in crystal -- move on to next crystal
                }
              }
          }

          //check ECAL hits in barrel
          for(Int_t l=0; l<eg4EBSignal->size(); l++) { //loop over crystals
              if(barrelMap.count(make_pair(ieta[l], iphi[l])) > 0) { //if has hits in jet
                double crystEnergy=0;
                double totalWeight = 0;
                double crystWeightedTime=0;
                double crystWeightedRho =0;
                double crystWeightedZ = 0;
                for(int m=0; m<eg4EBSignal->at(l).size(); m++) { //loop over Geant hits within crystal -- add up all energy and get mean time
                    crystEnergy += (*eg4EBSignal)[l][m];
                    double thisWeight = pow((*eg4EBSignal)[l][m], weightExp);
                    totalWeight += thisWeight;
                    crystWeightedTime += thisWeight*(*tg4EBSignal)[l][m];
                    crystWeightedRho += thisWeight*sqrt(pow((*prexg4EBSignal)[l][m], 2) + pow((*preyg4EBSignal)[l][m], 2));
                    crystWeightedZ += thisWeight*(*prezg4EBSignal)[l][m];
                }
                crystWeightedTime /= totalWeight;
                crystWeightedRho /= totalWeight;
                crystWeightedZ /= totalWeight;
                //update the maps
                barrelMap[make_pair(ieta[l], iphi[l])] = make_pair(crystEnergy/1000, crystWeightedTime);
                barrelPosMap[make_pair(ieta[l], iphi[l])] = make_pair(crystWeightedRho, crystWeightedZ);
                //store first hit time and position for this crystal
                if(crystEnergy/1000 > crystEThreshold){
                    thisPhotFirstTB.push_back(tming4EB[l]);
                    thisPhotFirstRhoB.push_back(sqrt(xtming4EB[l]*xtming4EB[l] + ytming4EB[l]*ytming4EB[l]));
                    thisPhotFirstZB.push_back(ztming4EB[l]);
                }
              }   
          }
        
          photFirstTimeB->push_back(thisPhotFirstTB);
          photFirstRhoB->push_back(thisPhotFirstRhoB);
          photFirstZB->push_back(thisPhotFirstZB);
          
          //Repeat for the endcap
          int thisPhotIZ = 0;
          map<pair<int,int>, pair<double,double> > endcapMap;
          map<pair<int,int>, pair<double,double> > endcapPosMap;
          for(Int_t l = 0; l < eg4EESignal->size(); l++){
            for(int m = 0; m < eg4EESignal->at(l).size(); m++){  
                TVector3 vEESignal;
                vEESignal.SetXYZ((*prexg4EESignal)[l][m], (*preyg4EESignal)[l][m], (*prezg4EESignal)[l][m] - vertexZ);
                if (cmsana::deltaR(vEESignal.Eta(),vEESignal.Phi(),genPhotons[photIndex]->eta,genPhotons[photIndex]->phi) < 0.3 && (*eg4EESignal)[l][m] > 0) {
                    if(thisPhotIZ == 0) {
                        thisPhotIZ = iz[l]; //(all endcap hits in the jet should come from the same endcap)
                    }
                    endcapMap[make_pair(ix[l], iy[l])] = make_pair(0.0, 0.0);
                    break; //found a hit in crystal -- move on to next crystal
                }
              }
          }
          
          vector<double> thisPhotFirstTE;
          vector<double> thisPhotFirstRhoE;
          vector<double> thisPhotFirstZE;

          //check ECAL hits in endcap
          for(Int_t l=0; l<eg4EESignal->size(); l++) { //loop over crystals
              if(endcapMap.count(make_pair(ix[l], iy[l])) > 0 && iz[l] == thisPhotIZ) { //if has hits in jet
                double crystEnergy=0;
                double totalWeight = 0;
                double crystWeightedTime=0;
                double crystWeightedRho = 0;
                double crystWeightedZ = 0;
                for(int m=0; m<eg4EESignal->at(l).size(); m++) { //loop over Geant hits within crystal -- add up all energy and get mean time
                    crystEnergy += (*eg4EESignal)[l][m];
                    double thisWeight = pow((*eg4EESignal)[l][m], weightExp);
                    totalWeight += thisWeight;
                    crystWeightedTime += thisWeight*(*tg4EESignal)[l][m];
                    crystWeightedRho += sqrt(pow((*prexg4EESignal)[l][m],2) + pow((*preyg4EESignal)[l][m], 2))*thisWeight;
                    crystWeightedZ += (*prezg4EESignal)[l][m]*thisWeight;
                }
                crystWeightedTime /= totalWeight;
                crystWeightedRho /= totalWeight;
                crystWeightedZ /= totalWeight;
                //update the map
                endcapMap[make_pair(ix[l], iy[l])] = make_pair(crystEnergy/1000, crystWeightedTime);
                endcapPosMap[make_pair(ix[l], iy[l])] = make_pair(crystWeightedRho, crystWeightedZ);
                if(crystEnergy/1000 > crystEThreshold){
                    thisPhotFirstTE.push_back(tming4EE[l]);
                    thisPhotFirstRhoE.push_back(sqrt(xtming4EE[l]*xtming4EE[l] + ytming4EE[l]*ytming4EE[l]));
                    thisPhotFirstZE.push_back(ztming4EE[l]);
                }
              }   
          }

          photFirstTimeE->push_back(thisPhotFirstTE);
          photFirstRhoE->push_back(thisPhotFirstRhoE);
          photFirstZE->push_back(thisPhotFirstZE);

          //TODO: add hits for pileup if mixPU = true

          //Fill tree variables
          photEta->push_back(genPhotons[photIndex]->eta);
          photPhi->push_back(genPhotons[photIndex]->phi);
          photPt->push_back(genPhotons[photIndex]->pt);

          //Get crystal energy & time
          
          vector<int> thisPhotIEtaB;
          vector<int> thisPhotIPhiB;
          vector<double> thisPhotEB;
          vector<double> thisPhotTB;
          vector<double> thisPhotE3x3B;
          vector<double> thisPhotT3x3B;
          vector<double> thisPhotRhoB;
          vector<double> thisPhotZB;

          //barrel
          for(map<pair<int,int>,pair<double,double> >::iterator iter = barrelMap.begin(); iter != barrelMap.end(); ++iter){
              double theE = (iter->second).first;
              if(theE > crystEThreshold){ //store cluster info in tree
                int theIEta = (iter->first).first;
                int theIPhi = (iter->first).second;
                cout << "Photons -- Storing info for barrel crystal with energy " << theE << " at (" << theIEta << ", " << theIPhi << ") " << endl;
                double theT = (iter->second).second;

                //seed crystal
                thisPhotEB.push_back(theE);
                thisPhotTB.push_back(theT);
                thisPhotIEtaB.push_back(theIEta);
                thisPhotIPhiB.push_back(theIPhi);
                pair<double,double> rhoZPair = barrelPosMap[make_pair(theIEta, theIPhi)];
                thisPhotRhoB.push_back(rhoZPair.first);
                thisPhotZB.push_back(rhoZPair.second);
                cout << "Rho: " << rhoZPair.first << "; Z: " << rhoZPair.second << endl;

                //3x3 array
                double e33 = 0; //energy
                double t33 = 0; //mean time
                for(int e = theIEta - 1; e <= theIEta + 1; e++){
                    for(int p = theIPhi - 1; p <= theIPhi + 1; p++){
                        if(barrelMap.count(make_pair(e, p)) > 0){ //if has hits
                            pair<double,double> thisCryst = barrelMap[make_pair(e,p)];
                            e33 += thisCryst.first;
                            t33 += thisCryst.second*thisCryst.first;
                        }
                    }
                }
                t33/=e33;
                thisPhotE3x3B.push_back(e33);
                thisPhotT3x3B.push_back(t33);

              }
          } 
//TODO: order entries by seed crystal energy
          photIEta->push_back(thisPhotIEtaB);
          photIPhi->push_back(thisPhotIPhiB);
          photEB->push_back(thisPhotEB);
          photTB->push_back(thisPhotTB);
          photE3x3B->push_back(thisPhotE3x3B);
          photT3x3B->push_back(thisPhotT3x3B);
          photRhoB->push_back(thisPhotRhoB);
          photZB->push_back(thisPhotZB);
 
          vector<int> thisPhotIXE;
          vector<int> thisPhotIYE;
          vector<int> thisPhotIZE;
          vector<double> thisPhotEE;
          vector<double> thisPhotTE;
          vector<double> thisPhotE3x3E;
          vector<double> thisPhotT3x3E;
          vector<double> thisPhotRhoE;
          vector<double> thisPhotZE;

          //endcap
          for(map<pair<int,int> , pair<double,double> >::iterator iter = endcapMap.begin(); iter != endcapMap.end(); ++iter){
              double theE = (iter->second).first;
              if(theE > crystEThreshold){ //store cluster info in tree
                int theIX = (iter->first).first;
                int theIY = (iter->first).second;
                cout << "Photons -- Storing info for endcap crystal with energy " << theE << " at (" << theIX << ", " << theIY << ", " << thisPhotIZ << ")" << endl;
                double theT = (iter->second).second;

                //seed crystal
                thisPhotEE.push_back(theE);
                thisPhotTE.push_back(theT);
                thisPhotIXE.push_back(theIX);
                thisPhotIYE.push_back(theIY);
                thisPhotIZE.push_back(thisPhotIZ);

                pair<double,double> rhoZPair = endcapPosMap[make_pair(theIX, theIY)];
                thisPhotRhoE.push_back(rhoZPair.first);
                thisPhotZE.push_back(rhoZPair.second);
                cout << "Rho: " << rhoZPair.first << "; Z: " << rhoZPair.second << endl;

                //3x3 array
                double e33 = 0; //energy
                double t33 = 0; //mean time
                for(int x = theIX - 1; x <= theIX + 1; x++){
                    for(int y = theIY - 1; y <= theIY + 1; y++){
                        if(endcapMap.count(make_pair(x, y)) > 0){ //if has hits
                            pair<double,double> thisCryst = endcapMap[make_pair(x,y)];
                            e33 += thisCryst.first;
                            t33 += thisCryst.second*thisCryst.first;
                        }
                    }
                }
                t33/=e33;
                thisPhotE3x3E.push_back(e33);
                thisPhotT3x3E.push_back(t33);
              }
          } 

          photIX->push_back(thisPhotIXE);
          photIY->push_back(thisPhotIYE);
          photIZ->push_back(thisPhotIZE);
          photEE->push_back(thisPhotEE);
          photTE->push_back(thisPhotTE);
          photE3x3E->push_back(thisPhotE3x3E);
          photT3x3E->push_back(thisPhotT3x3E);
          photRhoE->push_back(thisPhotRhoE);
          photZE->push_back(thisPhotZE);

          continue;
      }

      if (Jets[j].pt() > 30 && fabs(Jets[j].eta()) < 3) {
          //analyze selected jets
        cout << "GenJet " << j << " : " << Jets[j].pt() << " " << Jets[j].eta() << " " << Jets[j].phi() << "\n";
        
	//get the GEANT hits for the jet
        double jetEcalEnergySum = 0;
        vector<double> thisJetFirstTB;
        vector<double> thisJetFirstRhoB;
        vector<double> thisJetFirstZB;
     
          //store the crystal information in a map:
          //the key pair is ieta,iphi, and the pair of doubles is <energy, energy-weighted time>
          map<pair<int,int>, pair<double, double> > barrelMap;
          //second map to store position info -- the two doubles are rho and z
          map<pair<int,int>, pair<double, double> > barrelPosMap; 
          //initialize -- just check which crystals have hits in the jet
          for(Int_t l = 0; l < eg4EBSignal->size(); l++) {
              for(int m = 0; m < eg4EBSignal->at(l).size(); m++){  
                TVector3 vEBSignal;
                vEBSignal.SetXYZ((*prexg4EBSignal)[l][m], (*preyg4EBSignal)[l][m], (*prezg4EBSignal)[l][m] - vertexZ);
                if (cmsana::deltaR(vEBSignal.Eta(),vEBSignal.Phi(),Jets[j].eta(),Jets[j].phi()) < 0.5 && (*eg4EBSignal)[l][m] > 0) {
                    barrelMap[make_pair(ieta[l], iphi[l])] = make_pair(0.0, 0.0);
                    break; //found a hit in crystal -- move on to next crystal
                }
              }
          }

          //check ECAL hits in barrel
          for(Int_t l=0; l<eg4EBSignal->size(); l++) { //loop over crystals
              if(barrelMap.count(make_pair(ieta[l], iphi[l])) > 0) { //if has hits in jet
                double crystEnergy=0;
                double totalWeight=0;
                double crystWeightedTime=0;
                double crystWeightedRho =0;
                double crystWeightedZ = 0;
                for(int m=0; m<eg4EBSignal->at(l).size(); m++) { //loop over Geant hits within crystal -- add up all energy and get mean time
                    jetEcalEnergySum += (*eg4EBSignal)[l][m];
                    crystEnergy += (*eg4EBSignal)[l][m];
                    double thisWeight = pow((*eg4EBSignal)[l][m], weightExp);
                    totalWeight += thisWeight;
                    crystWeightedTime += thisWeight*(*tg4EBSignal)[l][m];
                    crystWeightedRho += thisWeight*sqrt(pow((*prexg4EBSignal)[l][m], 2) + pow((*preyg4EBSignal)[l][m], 2));
                    crystWeightedZ += thisWeight*(*prezg4EBSignal)[l][m];
                }
                crystWeightedTime /= totalWeight;
                crystWeightedRho /= totalWeight;
                crystWeightedZ /= totalWeight;
                //update the maps
                barrelMap[make_pair(ieta[l], iphi[l])] = make_pair(crystEnergy/1000, crystWeightedTime);
                barrelPosMap[make_pair(ieta[l], iphi[l])] = make_pair(crystWeightedRho, crystWeightedZ);
                //store first hit time and position for this crystal
                if(crystEnergy/1000 > crystEThreshold){
                    thisJetFirstTB.push_back(tming4EB[l]);
                    cout << "Got " << tming4EB[l] << endl;
                    thisJetFirstRhoB.push_back(sqrt(xtming4EB[l]*xtming4EB[l] + ytming4EB[l]*ytming4EB[l]));
                    thisJetFirstZB.push_back(ztming4EB[l]);
                }
              }   
          }

          firstTimeB->push_back(thisJetFirstTB);
          firstHitRhoB->push_back(thisJetFirstRhoB);
          firstHitZB->push_back(thisJetFirstZB);
          
          //Repeat for the endcap
          int thisJetIZ = 0;
          map<pair<int,int>, pair<double,double> > endcapMap;
          map<pair<int,int>, pair<double,double> > endcapPosMap;
          for(Int_t l = 0; l < eg4EESignal->size(); l++){
            for(int m = 0; m < eg4EESignal->at(l).size(); m++){  
                TVector3 vEESignal;
                vEESignal.SetXYZ((*prexg4EESignal)[l][m], (*preyg4EESignal)[l][m], (*prezg4EESignal)[l][m] - vertexZ);
                if (cmsana::deltaR(vEESignal.Eta(),vEESignal.Phi(),Jets[j].eta(),Jets[j].phi()) < 0.5 && (*eg4EESignal)[l][m] > 0) {
                    if(thisJetIZ == 0) {
                        thisJetIZ = iz[l]; //(all endcap hits in the jet should come from the same endcap)
                    }
                    endcapMap[make_pair(ix[l], iy[l])] = make_pair(0.0, 0.0);
                    break; //found a hit in crystal -- move on to next crystal
                }
              }
          }
          
          vector<double> thisJetFirstTE;
          vector<double> thisJetFirstRhoE;
          vector<double> thisJetFirstZE;

          //check ECAL hits in endcap
          for(Int_t l=0; l<eg4EESignal->size(); l++) { //loop over crystals
              if(endcapMap.count(make_pair(ix[l], iy[l])) > 0 && iz[l] == thisJetIZ) { //if has hits in jet
                double crystEnergy=0;
                double crystWeightedTime=0;
                double crystWeightedRho = 0;
                double crystWeightedZ = 0;
                double totalWeight = 0;
                for(int m=0; m<eg4EESignal->at(l).size(); m++) { //loop over Geant hits within crystal -- add up all energy and get mean time
                    jetEcalEnergySum += (*eg4EESignal)[l][m];
                    crystEnergy += (*eg4EESignal)[l][m];
                    double thisWeight = pow((*eg4EESignal)[l][m], weightExp);
                    totalWeight += thisWeight;
                    crystWeightedTime += thisWeight*(*tg4EESignal)[l][m];
                    crystWeightedRho += sqrt(pow((*prexg4EESignal)[l][m],2) + pow((*preyg4EESignal)[l][m], 2))*thisWeight;
                    crystWeightedZ += (*prezg4EESignal)[l][m]*thisWeight;
                }
                crystWeightedTime /= totalWeight;
                crystWeightedRho /= totalWeight;
                crystWeightedZ /= totalWeight;
                //update the map
                endcapMap[make_pair(ix[l], iy[l])] = make_pair(crystEnergy/1000, crystWeightedTime);
                endcapPosMap[make_pair(ix[l], iy[l])] = make_pair(crystWeightedRho, crystWeightedZ);

               if(crystEnergy/1000 > crystEThreshold){
                    thisJetFirstTE.push_back(tming4EE[l]);
                    thisJetFirstRhoE.push_back(sqrt(xtming4EE[l]*xtming4EE[l] + ytming4EE[l]*ytming4EE[l]));
                    thisJetFirstZE.push_back(ztming4EE[l]);
               }
              }   
          }

          firstTimeE->push_back(thisJetFirstTE);
          firstHitRhoE->push_back(thisJetFirstRhoE);
          firstHitZE->push_back(thisJetFirstZE);

          //TODO: add hits for pileup if mixPU = true

          jetEcalEnergySum /= 1000;
          cout << "jetEcalEnergySum [GeV]: " << jetEcalEnergySum << " , measured jet energy [GeV]: " << Jets[j].E() << " jet eta : " << Jets[j].eta() << "\n";
	  
          //Fill tree variables
          jEcalE->push_back(jetEcalEnergySum);
          jEta->push_back(Jets[j].eta());
          jPhi->push_back(Jets[j].phi());
          jPt->push_back(Jets[j].pt());

          //Get crystal energy & time
          
          vector<int> thisJetIEtaB;
          vector<int> thisJetIPhiB;
          vector<double> thisJetEB;
          vector<double> thisJetTB;
          vector<double> thisJetE3x3B;
          vector<double> thisJetT3x3B;
          vector<double> thisJetE5x5B;
          vector<double> thisJetT5x5B;
          vector<double> thisJetRhoB;
          vector<double> thisJetZB;

          //barrel
          for(map<pair<int,int>,pair<double,double> >::iterator iter = barrelMap.begin(); iter != barrelMap.end(); ++iter){
              double theE = (iter->second).first;
              if(theE > crystEThreshold){ //store cluster info in tree
                int theIEta = (iter->first).first;
                int theIPhi = (iter->first).second;
                cout << "Storing info for barrel crystal with energy " << theE << " at (" << theIEta << ", " << theIPhi << ") " << endl;
                double theT = (iter->second).second;

                //seed crystal
                thisJetEB.push_back(theE);
                thisJetTB.push_back(theT);
                thisJetIEtaB.push_back(theIEta);
                thisJetIPhiB.push_back(theIPhi);
                pair<double,double> rhoZPair = barrelPosMap[make_pair(theIEta, theIPhi)];
                thisJetRhoB.push_back(rhoZPair.first);
                thisJetZB.push_back(rhoZPair.second);
                cout << "Rho: " << rhoZPair.first << "; Z: " << rhoZPair.second << endl;

                //3x3 array
                double e33 = 0; //energy
                double t33 = 0; //mean time
                for(int e = theIEta - 1; e <= theIEta + 1; e++){
                    for(int p = theIPhi - 1; p <= theIPhi + 1; p++){
                        if(barrelMap.count(make_pair(e, p)) > 0){ //if has hits
                            pair<double,double> thisCryst = barrelMap[make_pair(e,p)];
                            e33 += thisCryst.first;
                            t33 += thisCryst.second*thisCryst.first;
                        }
                    }
                }
                t33/=e33;
                thisJetE3x3B.push_back(e33);
                thisJetT3x3B.push_back(t33);

                //5x5 array
                double e55 = 0; //energy
                double t55 = 0; //mean time
                for(int e = theIEta - 2; e <= theIEta + 2; e++){
                    for(int p = theIPhi - 2; p <= theIPhi + 2; p++){
                        if(barrelMap.count(make_pair(e,p)) > 0){
                            pair<double,double> thisCryst = barrelMap[make_pair(e,p)];
                            e55 += thisCryst.first;
                            t55 += thisCryst.second*thisCryst.first;
                        }
                    }
                }
                t55/=e55;
                thisJetE5x5B.push_back(e55);
                thisJetT5x5B.push_back(t55);
                cout << "Energies: " << theE << " " << e33 << " " << e55 << endl;
                cout << "Times: " << theT << " " << t33 << " " << t55 << endl;
              }
          } 
//TODO: order entries by seed crystal energy
          crystIEta->push_back(thisJetIEtaB);
          crystIPhi->push_back(thisJetIPhiB);
          crystEB->push_back(thisJetEB);
          crystTB->push_back(thisJetTB);
          crystE3x3B->push_back(thisJetE3x3B);
          crystT3x3B->push_back(thisJetT3x3B);
          crystE5x5B->push_back(thisJetE5x5B);
          crystT5x5B->push_back(thisJetT5x5B);
          rhoB->push_back(thisJetRhoB);
          hitZB->push_back(thisJetZB);
 
          vector<int> thisJetIXE;
          vector<int> thisJetIYE;
          vector<int> thisJetIZE;
          vector<double> thisJetEE;
          vector<double> thisJetTE;
          vector<double> thisJetE3x3E;
          vector<double> thisJetT3x3E;
          vector<double> thisJetE5x5E;
          vector<double> thisJetT5x5E;
          vector<double> thisJetRhoE;
          vector<double> thisJetZE;

          //endcap
          for(map<pair<int,int> , pair<double,double> >::iterator iter = endcapMap.begin(); iter != endcapMap.end(); ++iter){
              double theE = (iter->second).first;
              if(theE > crystEThreshold){ //store cluster info in tree
                int theIX = (iter->first).first;
                int theIY = (iter->first).second;
                cout << "Storing info for endcap crystal with energy " << theE << " at (" << theIX << ", " << theIY << ", " << thisJetIZ << ")" << endl;
                double theT = (iter->second).second;

                //seed crystal
                thisJetEE.push_back(theE);
                thisJetTE.push_back(theT);
                thisJetIXE.push_back(theIX);
                thisJetIYE.push_back(theIY);
                thisJetIZE.push_back(thisJetIZ);

                pair<double,double> rhoZPair = endcapPosMap[make_pair(theIX, theIY)];
                thisJetRhoE.push_back(rhoZPair.first);
                thisJetZE.push_back(rhoZPair.second);
                cout << "Rho: " << rhoZPair.first << "; Z: " << rhoZPair.second << endl;

                //3x3 array
                double e33 = 0; //energy
                double t33 = 0; //mean time
                for(int x = theIX - 1; x <= theIX + 1; x++){
                    for(int y = theIY - 1; y <= theIY + 1; y++){
                        if(endcapMap.count(make_pair(x, y)) > 0){ //if has hits
                            pair<double,double> thisCryst = endcapMap[make_pair(x,y)];
                            e33 += thisCryst.first;
                            t33 += thisCryst.second*thisCryst.first;
                        }
                    }
                }
                t33/=e33;
                thisJetE3x3E.push_back(e33);
                thisJetT3x3E.push_back(t33);

                //5x5 array
                double e55 = 0; //energy
                double t55 = 0; //mean time
                for(int x = theIX - 2; x <= theIX + 2; x++){
                    for(int y = theIY - 2; y <= theIY + 2; y++){
                        if(endcapMap.count(make_pair(x,y)) > 0){
                            pair<double,double> thisCryst = endcapMap[make_pair(x,y)];
                            e55 += thisCryst.first;
                            t55 += thisCryst.second*thisCryst.first;
                        }
                    }
                }
                t55/=e55;
                thisJetE5x5E.push_back(e55);
                thisJetT5x5E.push_back(t55);
                cout << "Energies: " << theE << " " << e33 << " " << e55 << endl;
                cout << "Times: " << theT << " " << t33 << " " << t55 << endl;
              }
          } 

          crystIX->push_back(thisJetIXE);
          crystIY->push_back(thisJetIYE);
          crystIZ->push_back(thisJetIZE);
          crystEE->push_back(thisJetEE);
          crystTE->push_back(thisJetTE);
          crystE3x3E->push_back(thisJetE3x3E);
          crystT3x3E->push_back(thisJetT3x3E);
          crystE5x5E->push_back(thisJetE5x5E);
          crystT5x5E->push_back(thisJetT5x5E);
          rhoE->push_back(thisJetRhoE);
          hitZE->push_back(thisJetZE);

      }//end this jet
    } //end loop over jets

    ecalTree->Fill();
  }//end of event loop

  fnew->cd();
  cout << "Writing the tree..." << endl;
  ecalTree->Write();
  fnew->Close();
  delete fnew;

  delete info;
  delete genparticleArr;
  delete genjetArr;

    return;

}

# ifndef __CINT__
int main(int argc, char* argv[]){
        gROOT->SetBatch();
        if(argc < 3){
            cerr << "usage: GetCrystalTimeInfo eventFileList geant4FileList" << endl;
            return -1;
        }

        string eventFileList;
        stringstream str1(argv[1]);
        str1 >> eventFileList;

        string geant4FileList;
        stringstream str2(argv[2]);
        str2 >> geant4FileList;

        CrystTimeInfo* cInfo = new CrystTimeInfo(eventFileList, geant4FileList);
        cInfo->GetCrystalTimeInfo();
        delete cInfo;
        return 0;
}
#endif
