//#ifdef __MAKECINT__
//#pragma link C++ class vector<vector<float> >+;
//#pragma link C++ class vector<vector<int> >+;
//#endif

#include <iostream>
#include <map>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "Rtypes.h"
#include "TInterpreter.h"
#include "TClonesArray.h"
#include <vector>
#include <fstream>
#include "RooRandom.h"

#include "CMSAna/Utils/CommonTools.hh"

// FastJet Stuff
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

using namespace std;

int minNumHits = 1; //minimum number of energetic hits required in each jet
int minNumPhotHits = 1;
int hitEThreshold = 30; //in GeV
int photEThreshold = 4;
const double speedoflight = 299.79245800; //in mm/ns
const double jPtCut = 30;
const double jEtaCut = 3.0;
const double photPtCut = 30;
const double photEtaCut = 3.0;
double timeSigma = 0;

//don't show this code to a computer scientist.

//dummy class for making the averaged imaginary time vertex function
class  AveragedTF1Object {
public:
    AveragedTF1Object(vector<TF1*> funcs, vector<double> barrelEnergy, vector<double> endcapEnergy, double thresh){
        theFuncs = funcs;
        eBarrel = barrelEnergy;
        eEndcap = endcapEnergy;
        eThreshold = thresh;
    }

    ~AveragedTF1Object(){
        for(int i = 0; i < theFuncs.size(); i++){
            delete theFuncs[i];
        }
    }

    double operator() (double *x, double *p) { //computes the average of all the TF1's
 
        double funcWeightedSum = 0;
        double funcTotalWeight = 0;
        bool isNaN = true; //check to see if there are any non-NaNs
        int funcIndex = 0;
        for(int i = 0; i < eBarrel.size(); i++){
            if(eBarrel[i] > eThreshold){
                if(!isnan(theFuncs[funcIndex]->Eval(x[0]))){
                    funcTotalWeight += sqrt(eBarrel[i]);
                    funcWeightedSum += sqrt(eBarrel[i])*theFuncs[funcIndex]->Eval(x[0]);
                    isNaN = false;
                }
                funcIndex++;
            }
        }
        for(int i = 0; i < eEndcap.size(); i++){
            if(eEndcap[i] > eThreshold){
                if(!isnan(theFuncs[funcIndex]->Eval(x[0]))){
                    funcTotalWeight += sqrt(eEndcap[i]);
                    funcWeightedSum += sqrt(eEndcap[i])*theFuncs[funcIndex]->Eval(x[0]);
                    isNaN = false;
                }
                funcIndex++;
            }
        }
        if(isNaN) return nan("");
        double avgFuncVal = funcWeightedSum / funcTotalWeight;
        return avgFuncVal;
     
    }

    vector<TF1*> theFuncs;
    vector<double> eBarrel;
    vector<double> eEndcap;
    double eThreshold;
};

//dummy function to find the intersection between two TF1 objects
//(why is ROOT so terrible?)
TF1 *avgVertFuncJet1, *avgVertFuncJet2;
double finter(double *x, double*par) {
   return TMath::Abs(avgVertFuncJet1->EvalPar(x,par) - avgVertFuncJet2->EvalPar(x,par));
}

//gets the number of hits in a jet having energy higher than some threshold
int numHitsAboveThreshold(vector<double> hitEnergiesB, vector<double> hitEnergiesE, int threshold){
    int numHits = 0;
    for(int i = 0; i < hitEnergiesB.size(); i++){
        if(hitEnergiesB[i] >= threshold) numHits++;
    }
    for(int i = 0; i < hitEnergiesE.size(); i++){
        if(hitEnergiesE[i] >= threshold) numHits++;
    }
    return numHits;
}

//gets the z position on the beamline corresponding to the time of the rechit
/*double tToZ(double t, double hitTime, double hitRho, double hitZ){
    double travelTime = hitTime - t;
    double sol1 = 9999999;
    double sol2 = 9999999;
    //if(speedoflight*speedoflight*travelTime*travelTime > hitRho*hitRho){
        sol1 = hitZ - sqrt(speedoflight*speedoflight*travelTime*travelTime - hitRho*hitRho);
        sol2 = hitZ + sqrt(speedoflight*speedoflight*travelTime*travelTime - hitRho*hitRho);
    //}
    if(fabs(sol1) < fabs(sol2)){
        return sol1;
    }
    else{
        return sol2;
    }
}
*/
double tToZ(Double_t *x, Double_t *pars){ //pars[0] = hitTime, pars[1] = hitRho, pars[2] = hitZ
    double travelTime = pars[0] - x[0];
    double sol1 = 9999999;
    double sol2 = 9999999;
    sol1 = pars[2] - sqrt(speedoflight*speedoflight*travelTime*travelTime - pars[1]*pars[1]);
    sol2 = pars[2] + sqrt(speedoflight*speedoflight*travelTime*travelTime - pars[1]*pars[1]);
    if(fabs(sol1) < fabs(sol2)){
        return sol1;
    }
    else{
        return sol2;
    }
}

double timeSmear(double time, double sigma){
    if(sigma == 0) return time;
    else{
        return time + RooRandom::gaussian()*sigma;
    }
}

//Returns a vector of TF1's corresponding to the imaginary vertex functions for the hits in the crystal
vector<TF1*> makeImVertFuncs(vector<double> crystTB, vector<double> crystEB, vector<double> rhoB, vector<double> hitZB, vector<double> crystTE, vector<double> crystEE, vector<double> rhoE, vector<double> hitZE, int threshold){
    vector<TF1*> theFuncs;
    //loop over barrel
    for(int crystI = 0; crystI < crystTB.size(); crystI++){
        if(crystEB[crystI] > threshold){
            TF1 *f = new TF1("f",tToZ,-7, 7, 3);
            double smearedTime = timeSmear(crystTB[crystI], timeSigma);
            f->SetParameters(smearedTime, rhoB[crystI], hitZB[crystI]);
            theFuncs.push_back(f);
        }
    }
    //loop over endcap
    for(int crystI = 0; crystI < crystTE.size(); crystI++){
        if(crystEE[crystI] > threshold){
            TF1 *f = new TF1("f", tToZ, -7, 7, 3);
            double smearedTime = timeSmear(crystTE[crystI], timeSigma);
            f->SetParameters(smearedTime, rhoE[crystI], hitZE[crystI]);
            theFuncs.push_back(f);
        }
    }
    return theFuncs;
}

void VertexEcalJets(int option) {
  //gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  //input file
  //TFile *infile = new TFile("CrystalTimeInfoTestAll.root");
  TFile *infile = new TFile("CrystalTimeInfoWeightESquared.root");
  TTree *ecalTree = (TTree*)infile->Get("ecalTree");

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

  jEcalE = new vector<double>; jEcalE->clear();
  jEta = new vector<double>; jEta->clear();
  jPhi = new vector<double>; jPhi->clear();
  jPt = new vector<double>; jPt->clear();

  crystIEta = new vector< vector<int> >; crystIEta->clear();
  crystIPhi = new vector< vector<int> >; crystIPhi->clear();
  crystIX = new vector< vector<int> >; crystIX->clear();
  crystIY = new vector< vector<int> >; crystIY->clear();
  crystIZ = new vector< vector<int> >; crystIZ->clear();

  crystEB = new vector< vector<double> >; crystEB->clear();
  crystTB = new vector< vector<double> >; crystTB->clear();
  crystE3x3B = new vector< vector<double> >; crystE3x3B->clear();
  crystT3x3B = new vector< vector<double> >; crystT3x3B->clear();
  crystE5x5B = new vector< vector<double> >; crystE5x5B->clear();
  crystT5x5B = new vector< vector<double> >; crystT5x5B->clear();
  crystEE = new vector< vector<double> >; crystEE->clear();
  crystTE = new vector< vector<double> >; crystTE->clear();
  crystE3x3E = new vector< vector<double> >; crystE3x3E->clear();
  crystT3x3E = new vector< vector<double> >; crystT3x3E->clear();
  crystE5x5E = new vector< vector<double> >; crystE5x5E->clear();
  crystT5x5E = new vector< vector<double> >; crystT5x5E->clear();

  rhoB = new vector< vector<double> >; rhoB->clear();
  rhoE = new vector< vector<double> >; rhoE->clear();
  hitZB = new vector< vector<double> >; hitZB->clear();
  hitZE = new vector< vector<double> >; hitZE->clear();

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

  ecalTree->SetBranchAddress("vertexZ", &vertexZ);

  ecalTree->SetBranchAddress("jEcalE", &jEcalE);
  ecalTree->SetBranchAddress("jEta", &jEta);
  ecalTree->SetBranchAddress("jPhi", &jPhi);
  ecalTree->SetBranchAddress("jPt", &jPt);

  ecalTree->SetBranchAddress("crystIEta", &crystIEta);
  ecalTree->SetBranchAddress("crystIPhi", &crystIPhi);
  ecalTree->SetBranchAddress("crystEB", &crystEB);
  ecalTree->SetBranchAddress("crystTB", &crystTB);
  ecalTree->SetBranchAddress("crystE3x3B", &crystE3x3B);
  ecalTree->SetBranchAddress("crystT3x3B", &crystT3x3B);
  ecalTree->SetBranchAddress("crystE5x5B", &crystE5x5B);
  ecalTree->SetBranchAddress("crystT5x5B", &crystT5x5B);

  ecalTree->SetBranchAddress("crystIX", &crystIX);
  ecalTree->SetBranchAddress("crystIY", &crystIY);
  ecalTree->SetBranchAddress("crystIZ", &crystIZ);
  ecalTree->SetBranchAddress("crystEE", &crystEE);
  ecalTree->SetBranchAddress("crystTE", &crystTE);
  ecalTree->SetBranchAddress("crystE3x3E", &crystE3x3E);
  ecalTree->SetBranchAddress("crystT3x3E", &crystT3x3E);
  ecalTree->SetBranchAddress("crystE5x5E", &crystE5x5E);
  ecalTree->SetBranchAddress("crystT5x5E", &crystT5x5E);

  ecalTree->SetBranchAddress("rhoB", &rhoB);
  ecalTree->SetBranchAddress("rhoE", &rhoE);
  ecalTree->SetBranchAddress("hitZB", &hitZB);
  ecalTree->SetBranchAddress("hitZE", &hitZE);

  ecalTree->SetBranchAddress("firstTimeB", &firstTimeB);
  ecalTree->SetBranchAddress("firstTimeE", &firstTimeE);
  ecalTree->SetBranchAddress("firstHitRhoB", &firstHitRhoB);
  ecalTree->SetBranchAddress("firstHitRhoE", &firstHitRhoE);
  ecalTree->SetBranchAddress("firstHitZB", &firstHitZB);
  ecalTree->SetBranchAddress("firstHitZE", &firstHitZE);

  ecalTree->SetBranchAddress("photEta", &photEta);
  ecalTree->SetBranchAddress("photPhi", &photPhi);
  ecalTree->SetBranchAddress("photPt", &photPt);
  ecalTree->SetBranchAddress("photIEta", &photIEta);
  ecalTree->SetBranchAddress("photIPhi", &photIPhi);
  ecalTree->SetBranchAddress("photIX", &photIX);
  ecalTree->SetBranchAddress("photIY", &photIY);
  ecalTree->SetBranchAddress("photIZ", &photIZ);
  ecalTree->SetBranchAddress("photEB", &photEB);
  ecalTree->SetBranchAddress("photE3x3B", &photE3x3B);
  ecalTree->SetBranchAddress("photTB", &photTB);
  ecalTree->SetBranchAddress("photT3x3B", &photT3x3B);
  ecalTree->SetBranchAddress("photEE", &photEE);
  ecalTree->SetBranchAddress("photE3x3E", &photE3x3E);
  ecalTree->SetBranchAddress("photTE", &photTE);
  ecalTree->SetBranchAddress("photT3x3E", &photT3x3E);
  ecalTree->SetBranchAddress("photRhoB", &photRhoB);
  ecalTree->SetBranchAddress("photRhoE", &photRhoE);
  ecalTree->SetBranchAddress("photZB", &photZB);
  ecalTree->SetBranchAddress("photZE", &photZE);
  ecalTree->SetBranchAddress("photFirstTimeB", &photFirstTimeB);
  ecalTree->SetBranchAddress("photFirstTimeE", &photFirstTimeE);
  ecalTree->SetBranchAddress("photFirstRhoB", &photFirstRhoB);
  ecalTree->SetBranchAddress("photFirstRhoE", &photFirstRhoE);
  ecalTree->SetBranchAddress("photFirstZB", &photFirstZB);
  ecalTree->SetBranchAddress("photFirstZE", &photFirstZE);

  string fOutName = Form("JetVertexInfoMinHits%dMinEnergy%dOption%dSmear%3f.root", minNumHits, hitEThreshold, option, timeSigma);
  TFile fOut(fOutName.c_str(), "recreate");
  TH1D *vertexHist = new TH1D("vertexHist", "; vertex z (timing) - vertex z (tracker) (mm); Number of events; Distance between reconstructed vertex and tracker vertex", 2000, -250, 250);
  TH2D *vertexHist2D = new TH2D("vertexHist2D", ";vertex z (tracker) (mm); vertex z (timing) (mm); Number of events; Vertexing correlation plot", 250, -250, 250, 250, -250, 250);
  TH2D *meanJetTimeEta = new TH2D("meanJetTimeEta", ";jet #eta; mean time of jet hits; Number of events; Mean time of jet hits", 60, -3, 3, 100, -5, 5);
  TH2D *meanCrystTimeIEta = new TH2D("meanCrystTimeIEta", ";crystal i#eta; mean time of hits; Number of events; Mean time of crystal hits", 171, -85, 86, 200, -20, 20);
  //TH3D *meanCrystTimeIXIY = new TH3D("meanCrystTimeIXIY", ";ix; iy; mean time of hits; Mean time of crystal hits", 100, 1, 101, 100, 1, 101, 100, -20, 20);
  TH2D *meanCrystTimeIXIY = new TH2D("meanCrystTimeIXIY", ";#sqrt{ix^{2} + iy^{2}}; mean time of hits; Number of Events; Mean time of crystal hits in endcap", 145, 0, 145, 200, -20, 20);

  TH1D *etaOfVertexedJets = new TH1D("etaOfVertexedJets", ";#eta; Number of Events", 30, -3, 3);
  TH1D *ptOfVertexedJets = new TH1D("ptOfVertexedJets", ";p_{T}; Number of Events", 50, 30, 200);
  TH2D *etaECrysts = new TH2D("etaECrysts", ";#eta; E; Number of Events", 60, -3, 3, 60, 30, 200);
  TH2D* vertexHist2DBarrel = new TH2D("vertexHist2DBarrel", ";vertex z (tracker) (mm); vertex z (timing) (mm); Number of Events; Vertexing correlation plot", 50, -250, 250, 50, -250, 250);
  TH2D* vertexHist2DEndcap = new TH2D("vertexHist2DEndcap", ";vertex z (tracker) (mm); vertex z (timing) (mm); Number of Events; Vertexing correlation plot", 50, -250, 250, 50, -250, 250);


  TH1D *vertexHistPhot = new TH1D("vertexHistPhot", "; vertex z (timing) - vertex z (tracker) (mm); Number of events; Distance between reconstructed vertex and tracker vertex", 1000, -250, 250);
  TH2D *vertexHistPhot2D = new TH2D("vertexHistPhot2D", ";vertex z (tracker) (mm); vertex z (timing) (mm); Number of events; Vertexing correlation plot", 250, -250, 250, 250, -250, 250);
  TH2D *jetVertexVsPhotVertex = new TH2D("jetVertexVsPhotVertex", ";vertex z (jets) (mm); vertex z (photons) (mm); Number of Events", 250, -250, 250, 250, -250, 250);
  TH2D *meanTimeVsFirstTime = new TH2D("meanTimeVsFirstTime", ";mean time;first hit time;Number of Events", 100, 3, 20, 100, 3, 15);
  TH2D *meanTimeVsFirstTimeBarrel = new TH2D("meanTimeVsFirstTimeBarrel", ";mean time;first hit time;Number of Events", 100, 3, 20, 100, 3, 15);
  TH2D *meanTimeVsFirstTimePhot = new TH2D("meanTimeVsFirstTimePhot", ";mean time;first hit time;Number of Events", 100, 3, 20, 100, 3, 15);
  TH2D *meanTimeMinusFirstTimeVsEta = new TH2D("meanTimeMinusFirstTimeVsEta", ";ieta;mean time - first hit time;Number of Events", 171, -85, 86, 100, -5, 5);
  TH2D *meanRhoMinusFirstRhoVsEta = new TH2D("meanRhoMinusFirstRhoVsEta", ";ieta;mean #rho - first hit #rho;Number of Events", 171, -85, 86, 100, -500, 500);
  TH2D *meanZMinusFirstZVsEta = new TH2D("meanZMinusFirstZVsEta", ";ieta;mean Z - first hit Z;Number of Events", 171, -85, 86, 100, -500, 500);
  TH2D *meanTimeMinusFirstTimeEE = new TH2D("meanTimeMinusFirstTimeEE", ";#sqrt{ix^{2}+iy^{2}};mean time - first hit time;Number of Events", 145,0,145, 100, -5, 5);
  TH2D *meanRhoMinusFirstRhoEE = new TH2D("meanRhoMinusFirstRhoEE", ";#sqrt{ix^{2}+iy^{2}};mean #rho - first hit #rho;Number of Events", 145,0,145, 100, -500, 500);
  TH2D *meanZMinusFirstZEE = new TH2D("meanZMinusFirstZEE", ";#sqrt{ix^{2}+iy^{2}};mean Z - first hit Z;Number of Events", 145,0,145, 100, -500, 500);

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  int eventsDone = 0;
  bool jetdone = false;
  //LOOP OVER EVENTS
  double vertexFromJets = -999;
  double vertexFromPhots = -999;
  int jetsInOppositeHalves = 0;
  int jetsHaveEnoughHits = 0;
  int jetsPassCuts = 0;
  int jetsVertexed = 0;
  int photsInOppositeHalves = 0;
  int photsHaveEnoughHits = 0;
  int photsPassCuts = 0;
  int photsVertexed = 0;
  int enoughJets = 0;
  int enoughPhots = 0;
  for(UInt_t ientry=0; ientry < ecalTree->GetEntries(); ientry++) { 
    
    jEcalE->clear();
    jEta->clear();
    jPhi->clear();
    jPt->clear();

    crystIEta->clear();
    crystIPhi->clear();
    crystIX->clear();
    crystIY->clear();
    crystIZ->clear();

    crystEB->clear();
    crystTB->clear();
    crystE3x3B->clear();
    crystT3x3B->clear();
    crystE5x5B->clear();
    crystT5x5B->clear();
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

    ecalTree->GetEntry(ientry);
    //cout << "\n\nEvent " << ientry << "\n";

    bool vertexedJets = false;
    bool vertexedPhots = false;
    //cout << "jEcalE->size() = " << jEcalE->size() << endl;
    //if(jEcalE->size() < 1) cout << "No jets." << endl;
    //else if(jEcalE->size() == 1) cout << "Only one jet." << endl;
    if(jEcalE->size() >= 2){ //if we have two jets, try to vertex them
        //Note that jets are sorted in decreasing order by pT, and we will just try to vertex the two most energetic jets for now.  
        enoughJets++;
        bool vertexable = true;

        //Get mean time of jets
        //cout << "Getting jet times..." << endl;
        for(int i = 0; i < crystTB->size(); i++){
            double meanJetTime = 0;
            double totalJetE = 0;
            double meanJetRho = 0;
            double meanJetZ = 0;
            for(int cr = 0; cr < crystTB->at(i).size(); cr++){
                if((*crystEB)[i][cr] >= hitEThreshold){
                    totalJetE += (*crystEB)[i][cr];
                    meanJetTime += (*crystEB)[i][cr]*(*crystTB)[i][cr];
                    meanJetRho += (*crystEB)[i][cr]*(*rhoB)[i][cr];
                    meanJetZ += (*crystEB)[i][cr]*(*hitZB)[i][cr];
                    meanCrystTimeIEta->Fill((*crystIEta)[i][cr], (*crystTB)[i][cr] - sqrt((*hitZB)[i][cr]*(*hitZB)[i][cr] + (*rhoB)[i][cr]*(*rhoB)[i][cr])/speedoflight, (*crystEB)[i][cr]);
                    meanTimeVsFirstTime->Fill((*crystT3x3B)[i][cr], (*firstTimeB)[i][cr]);
                    meanTimeVsFirstTimeBarrel->Fill((*crystT3x3B)[i][cr], (*firstTimeB)[i][cr]);
                    meanTimeMinusFirstTimeVsEta->Fill((*crystIEta)[i][cr], (*crystT3x3B)[i][cr]-(*firstTimeB)[i][cr]);
                    meanRhoMinusFirstRhoVsEta->Fill((*crystIEta)[i][cr],(*rhoB)[i][cr]-(*firstHitRhoB)[i][cr]);
                    meanZMinusFirstZVsEta->Fill((*crystIEta)[i][cr], (*hitZB)[i][cr]-(*firstHitZB)[i][cr]);
                }
            }
            for(int cr = 0; cr < crystTE->at(i).size(); cr++){
                if((*crystEE)[i][cr] >= hitEThreshold){
                    totalJetE += (*crystEE)[i][cr];
                    meanJetTime += (*crystEE)[i][cr]*(*crystTE)[i][cr];
                    meanJetRho += (*crystEE)[i][cr]*(*rhoE)[i][cr];
                    meanJetZ += (*crystEE)[i][cr]*(*hitZE)[i][cr];
                    //meanCrystTimeIXIY->Fill((*crystIX)[i][c], (*crystIY)[i][c], (*crystTE)[i][c] - sqrt((*hitZE)[i][c]*(*hitZE)[i][c] + (*rhoE)[i][c]*(*rhoE)[i][c])/speedoflight, (*crystEE)[i][c]);
                    meanCrystTimeIXIY->Fill(sqrt((*crystIX)[i][cr]*(*crystIX)[i][cr] + (*crystIY)[i][cr]*(*crystIY)[i][cr]), (*crystTE)[i][cr] - sqrt((*hitZE)[i][cr]*(*hitZE)[i][cr] + (*rhoE)[i][cr]*(*rhoE)[i][cr])/speedoflight, (*crystEE)[i][cr]);
                    meanTimeVsFirstTime->Fill((*crystT3x3E)[i][cr], (*firstTimeE)[i][cr]);
                    meanTimeMinusFirstTimeEE->Fill(sqrt((*crystIX)[i][cr]*(*crystIX)[i][cr] + (*crystIY)[i][cr]*(*crystIY)[i][cr]),(*crystT3x3E)[i][cr]-(*firstTimeE)[i][cr]);
                    meanRhoMinusFirstRhoEE->Fill(sqrt((*crystIX)[i][cr]*(*crystIX)[i][cr] + (*crystIY)[i][cr]*(*crystIY)[i][cr]),(*rhoE)[i][cr]-(*firstHitRhoE)[i][cr]);
                    meanZMinusFirstZEE->Fill(sqrt((*crystIX)[i][cr]*(*crystIX)[i][cr] + (*crystIY)[i][cr]*(*crystIY)[i][cr]),(*hitZE)[i][cr]-(*firstHitZE)[i][cr]);
                }
            }
            meanJetTime /= totalJetE;
            meanJetRho /= totalJetE;
            //meanJetZ /= totalJetE;
            meanJetZ = meanJetRho*sinh((*jEta)[i]);
            meanJetTimeEta->Fill((*jEta)[i], meanJetTime - sqrt(meanJetZ*meanJetZ + meanJetRho*meanJetRho)/speedoflight/*, (*jEcalE)[i]*/);
        }

        //Check if they are in opposite halves of the ECAL
        if(((*jEta)[0] >= 0 && (*jEta)[1] >= 0) || ((*jEta)[0] <= 0 && (*jEta)[1] <= 0)){
            //cout << "Jets are in the same half of the ECAL." << endl;
            vertexable = false;
        }
        else if(vertexable) jetsInOppositeHalves++;
        //Check if they pass pT and eta cuts
        if((*jPt)[0] < jPtCut || (*jPt)[1] < jPtCut || fabs((*jEta)[0]) > jEtaCut || fabs((*jEta)[1]) > jEtaCut){
            //cout << "A jet failed a cut!" << endl;
            vertexable = false;
        }
        else if(vertexable) jetsPassCuts++;

        //Check if the jet has enough hits 
        int numHitsJet1 = numHitsAboveThreshold((*crystEB)[0], (*crystEE)[0], hitEThreshold);
        int numHitsJet2 = numHitsAboveThreshold((*crystEB)[1], (*crystEE)[1], hitEThreshold);
        if(numHitsJet1 < minNumHits || numHitsJet2 < minNumHits){
            //cout << "At least one jet has not enough hits." << endl;
            vertexable = false;
        }
        else if(vertexable) jetsHaveEnoughHits++;

    if(vertexable){
        //cout << "Vertexing jets..." << endl;
        //cout << "Number of hits: " << numHitsJet1 << " " << numHitsJet2 << endl;
        
        //Make tVertex function for each hit -- use the 3x3 energies and times
        //cout << "   Making imaginary vertex functions..." << endl;
        vector<TF1*> imVertFuncsJet1;
        vector<TF1*> imVertFuncsJet2;

        if(option == 0){ //vertex jets using mean position and time in each crystal
            imVertFuncsJet1 = makeImVertFuncs((*crystT3x3B)[0], (*crystEB)[0], (*rhoB)[0], (*hitZB)[0], (*crystT3x3E)[0], (*crystEE)[0], (*rhoE)[0], (*hitZE)[0], hitEThreshold);
            imVertFuncsJet2 = makeImVertFuncs((*crystT3x3B)[1], (*crystEB)[1], (*rhoB)[1], (*hitZB)[1], (*crystT3x3E)[1], (*crystEE)[1], (*rhoE)[1], (*hitZE)[1], hitEThreshold);
        }
        else if(option == 1){ //vertex jets using first hit time and position in each crystal
            imVertFuncsJet1 = makeImVertFuncs((*firstTimeB)[0], (*crystEB)[0], (*firstHitRhoB)[0], (*firstHitZB)[0], (*firstTimeE)[0], (*crystEE)[0], (*firstHitRhoE)[0], (*firstHitZE)[0], hitEThreshold);
            imVertFuncsJet2 = makeImVertFuncs((*firstTimeB)[1], (*crystEB)[1], (*firstHitRhoB)[1], (*firstHitZB)[1], (*firstTimeE)[1], (*crystEE)[1], (*firstHitRhoE)[1], (*firstHitZE)[1], hitEThreshold);
        }
        else if(option == 2){ //use first hit time, mean position
            imVertFuncsJet1 = makeImVertFuncs((*firstTimeB)[0], (*crystEB)[0], (*rhoB)[0], (*hitZB)[0], (*firstTimeE)[0], (*crystEE)[0], (*rhoE)[0], (*hitZE)[0], hitEThreshold);
            imVertFuncsJet2 = makeImVertFuncs((*firstTimeB)[1], (*crystEB)[1], (*rhoB)[1], (*hitZB)[1], (*firstTimeE)[1], (*crystEE)[1], (*rhoE)[1], (*hitZE)[1], hitEThreshold);
        }
        else if(option == 3){//use mean time and first hit position
            imVertFuncsJet1 = makeImVertFuncs((*crystT3x3B)[0], (*crystEB)[0], (*firstHitRhoB)[0], (*firstHitZB)[0], (*crystT3x3E)[0], (*crystEE)[0], (*firstHitRhoE)[0], (*firstHitZE)[0], hitEThreshold);
            imVertFuncsJet2 = makeImVertFuncs((*crystT3x3B)[1], (*crystEB)[1], (*firstHitRhoB)[1], (*firstHitZB)[1], (*crystT3x3E)[1], (*crystEE)[1], (*firstHitRhoE)[1], (*firstHitZE)[1], hitEThreshold);
        }
        else if(option == 4){ //use mean time in single crystal
            imVertFuncsJet1 = makeImVertFuncs((*crystTB)[0], (*crystEB)[0], (*rhoB)[0], (*hitZB)[0], (*crystTE)[0], (*crystEE)[0], (*rhoE)[0], (*hitZE)[0], hitEThreshold);
            imVertFuncsJet2 = makeImVertFuncs((*crystTB)[1], (*crystEB)[1], (*rhoB)[1], (*hitZB)[1], (*crystTE)[1], (*crystEE)[1], (*rhoE)[1], (*hitZE)[1], hitEThreshold);
        }   
         else if(option == 5){ //use mean time in single crystal, and first hit position
            imVertFuncsJet1 = makeImVertFuncs((*crystTB)[0], (*crystEB)[0], (*firstHitRhoB)[0], (*firstHitZB)[0], (*crystTE)[0], (*crystEE)[0], (*firstHitRhoE)[0], (*firstHitZE)[0], hitEThreshold);
            imVertFuncsJet2 = makeImVertFuncs((*crystTB)[1], (*crystEB)[1], (*firstHitRhoB)[1], (*firstHitZB)[1], (*crystTE)[1], (*crystEE)[1], (*firstHitRhoE)[1], (*firstHitZE)[1], hitEThreshold);
        }   
        else{
            cout << "ERROR: invalid run option.  Use 0-5." << endl;
            break;
        }

        //Apply averaging -- sqrt(hit energy)
        //cout << "   Averaging hits..." << endl;
        AveragedTF1Object *avObj1 = new AveragedTF1Object(imVertFuncsJet1, (*crystEB)[0], (*crystEE)[0], hitEThreshold);
        avgVertFuncJet1 = new TF1("avgVertFuncJet1", avObj1,-7, 7, 0, "AveragedTF1Object"); 
        AveragedTF1Object *avObj2 = new AveragedTF1Object(imVertFuncsJet2, (*crystEB)[1], (*crystEE)[1], hitEThreshold);
        avgVertFuncJet2 = new TF1("avgVertFuncJet2", avObj2,-7, 7, 0, "AveragedTF1Object"); 
        
        //Print vertex functions for individual hits, and overlay the average
        if(!jetdone){
            
            avgVertFuncJet1->SetLineColor(kRed);
            avgVertFuncJet2->SetLineColor(kRed);
            avgVertFuncJet1->SetLineWidth(3);
            avgVertFuncJet2->SetLineWidth(3);
            avgVertFuncJet1->Draw();
            avgVertFuncJet2->Draw("same");
            for(int i = 0; i < imVertFuncsJet1.size(); i++){
                imVertFuncsJet1[i]->SetLineColor(kBlack);
                imVertFuncsJet1[i]->Draw("same");
            }
            for(int i = 0; i < imVertFuncsJet2.size(); i++){
                imVertFuncsJet2[i]->SetLineColor(kBlack);
                imVertFuncsJet2[i]->Draw("same");
            }
            avgVertFuncJet1->Draw("same");
            avgVertFuncJet2->Draw("same");
            fOut.cd();
            //fOut.WriteTObject(c, c->GetName(), "Write");
            c->Print("vertexJetDiagramAll.pdf");

            jetdone = true;
        }

        //TODO: Beam spot correction
        
        //Intersect the two averaged imaginary vertex functions to estimate the vertex position
        //cout << "   Getting estimated vertex position..." << endl;
        TF1 *fint = new TF1("fint",finter,-7, 7, 0);
        if(fint->GetMinimum(-7, 7) > 0.001){
            //cout << "Jets could not be vertexed! (tVertex functions don't intersect)" <<endl;
            vertexable = false;
        }
        double xint = fint->GetMinimumX(-7, 7);
        double vertexGuess = avgVertFuncJet1->Eval(xint);
        if(vertexable){
            vertexHist->Fill(vertexGuess - vertexZ*10);
            vertexHist2D->Fill(vertexZ*10, vertexGuess);
            cout << "Vertexing successful.  Estimated vertex Z: " << vertexGuess << ".  Tracker vertex Z: " << vertexZ*10 << endl;
            vertexedJets = true;
            jetsVertexed++;
            vertexFromJets = vertexGuess;
            etaOfVertexedJets->Fill((*jEta)[0]);
            ptOfVertexedJets->Fill((*jPt)[0]);
            if(fabs((*jEta)[0]) < 1.479){ //barrel
                vertexHist2DBarrel->Fill(vertexZ*10, vertexGuess);
            }
            else{
                vertexHist2DEndcap->Fill(vertexZ*10, vertexGuess);
            }
            for(int i = 0; i < (*crystEB)[0].size(); i++){
                etaECrysts->Fill(asinh((*hitZB)[0][i]/(*rhoB)[0][i]), (*crystEB)[0][i]);
            }
            for(int i = 0; i < (*crystEE)[0].size(); i++){
                etaECrysts->Fill(asinh((*hitZE)[0][i]/(*rhoE)[0][i]), (*crystEE)[0][i]);
            }
            eventsDone++;
        }

        delete fint;
        delete avObj1;
        delete avObj2;
        delete avgVertFuncJet1;
        delete avgVertFuncJet2;
    }
    }//done vertexing jets

    //Now vertex photons
    //if(photPt->size() < 1) cout << "No photons." << endl;
    //else if(photPt->size() == 1) //cout << "Only one photon." << endl;
    if(photPt->size() >= 2){ //if we have two jets, try to vertex them
        enoughPhots++;

        //Get photon info
        for(int i = 0; i < photTB->size(); i++){
            for(int cr = 0; cr < photTB->at(i).size(); cr++){
                if((*photEB)[i][cr] >= hitEThreshold){
                    meanTimeVsFirstTimePhot->Fill((*photT3x3B)[i][cr], (*photFirstTimeB)[i][cr]);
                }
            }
        }

        for(int i = 0; i < photTE->size(); i++){
            for(int cr = 0; cr < photTE->at(i).size(); cr++){
                if((*photEE)[i][cr] >= hitEThreshold){
                    meanTimeVsFirstTimePhot->Fill((*photT3x3E)[i][cr], (*photFirstTimeE)[i][cr]);
                }
            }
        }

        //Check if they are in opposite halves of the ECAL
        if(((*photEta)[0] >= 0 && (*photEta)[1] >= 0) || ((*photEta)[0] <= 0 && (*photEta)[1] <= 0)){
            //cout << "Photons are in the same half of the ECAL." << endl;
            continue;
        }
        else photsInOppositeHalves++;
        //Check if they pass pT and eta cuts
        if((*photPt)[0] < photPtCut || (*photPt)[1] < photPtCut || fabs((*photEta)[0]) > photEtaCut || fabs((*photEta)[1]) > photEtaCut){
            //cout << "A photon failed a cut!" << endl;
            continue;
        }
        else photsPassCuts++;

        //Check if the photon has enough hits 
        int numHitsPhot1 = numHitsAboveThreshold((*photEB)[0], (*photEE)[0], photEThreshold);
        int numHitsPhot2 = numHitsAboveThreshold((*photEB)[1], (*photEE)[1], photEThreshold);
        if(numHitsPhot1 < minNumPhotHits || numHitsPhot2 < minNumPhotHits){
            //cout << "At least one photon has not enough hits." << endl;
            continue;
        }
        else photsHaveEnoughHits++;

        //cout << "Vertexing photons..." << endl;
        //cout << "Number of hits: " << numHitsPhot1 << " " << numHitsPhot2 << endl;
        
        //Make tVertex function for each hit -- use the 3x3 energies and times
        //cout << "   Making tVertex functions..." << endl;
        vector<TF1*> imVertFuncsPhot1;
        vector<TF1*> imVertFuncsPhot2;

        //Find the most energetic hit in the R=0.3 cone around the gen photon
        int bestPhotHitIndex1 = -1;
        double bestPhotEnergy1 = 0;
        bool phot1IsBarrel = true;
        for(int i = 0; i < (*photEB)[0].size(); i++){
            if((*photEB)[0][i] > bestPhotEnergy1){
                bestPhotEnergy1 = (*photEB)[0][i];
                bestPhotHitIndex1 = i;
            }
        }
        for(int i = 0; i < (*photEE)[0].size(); i++){
            if((*photEE)[0][i] > bestPhotEnergy1){
                bestPhotEnergy1 = (*photEE)[0][i];
                bestPhotHitIndex1 = i;
                phot1IsBarrel = false; //the most energetic hit is in the endcap
            }
        }

        int bestPhotHitIndex2 = -1;
        double bestPhotEnergy2 = 0;
        bool phot2IsBarrel = true;
        for(int i = 0; i < (*photEB)[1].size(); i++){
            if((*photEB)[1][i] > bestPhotEnergy2){
                bestPhotEnergy2 = (*photEB)[1][i];
                bestPhotHitIndex2 = i;
            }
        }
        for(int i = 0; i < (*photEE)[1].size(); i++){
            if((*photEE)[1][i] > bestPhotEnergy2){
                bestPhotEnergy2 = (*photEE)[1][i];
                bestPhotHitIndex2 = i;
                phot2IsBarrel = false;
            }
        }

        //kludge -- give the vertex function maker only a single hit for each photon
        vector<double> phot1EB;
        vector<double> phot1TB;
        vector<double> phot1RhoB;
        vector<double> phot1ZB;
        vector<double> phot1EE;
        vector<double> phot1TE;
        vector<double> phot1RhoE;
        vector<double> phot1ZE;
        vector<double> phot2EB;
        vector<double> phot2TB;
        vector<double> phot2RhoB;
        vector<double> phot2ZB;
        vector<double> phot2EE;
        vector<double> phot2TE;
        vector<double> phot2RhoE;
        vector<double> phot2ZE;
        if(option == 0){
            if(phot1IsBarrel){
                phot1EB.push_back((*photEB)[0][bestPhotHitIndex1]);
                phot1TB.push_back((*photT3x3B)[0][bestPhotHitIndex1]);
                phot1RhoB.push_back((*photRhoB)[0][bestPhotHitIndex1]);
                phot1ZB.push_back((*photZB)[0][bestPhotHitIndex1]);
            }
            else{
                phot1EE.push_back((*photEE)[0][bestPhotHitIndex1]);
                phot1TE.push_back((*photT3x3E)[0][bestPhotHitIndex1]);
                phot1RhoE.push_back((*photRhoE)[0][bestPhotHitIndex1]);
                phot1ZE.push_back((*photZE)[0][bestPhotHitIndex1]);
            }

            if(phot2IsBarrel){
                phot2EB.push_back((*photEB)[1][bestPhotHitIndex2]);
                phot2TB.push_back((*photT3x3B)[1][bestPhotHitIndex2]);
                phot2RhoB.push_back((*photRhoB)[1][bestPhotHitIndex2]);
                phot2ZB.push_back((*photZB)[1][bestPhotHitIndex2]);
            }
            else{
                phot2EE.push_back((*photEE)[1][bestPhotHitIndex2]);
                phot2TE.push_back((*photT3x3E)[1][bestPhotHitIndex2]);
                phot2RhoE.push_back((*photRhoE)[1][bestPhotHitIndex2]);
                phot2ZE.push_back((*photZE)[1][bestPhotHitIndex2]);
            }
        }

        else if(option == 1){
            if(phot1IsBarrel){
                phot1EB.push_back((*photEB)[0][bestPhotHitIndex1]);
                phot1TB.push_back((*photFirstTimeB)[0][bestPhotHitIndex1]);
                phot1RhoB.push_back((*photFirstRhoB)[0][bestPhotHitIndex1]);
                phot1ZB.push_back((*photFirstZB)[0][bestPhotHitIndex1]);
            }
            else{
                phot1EE.push_back((*photEE)[0][bestPhotHitIndex1]);
                phot1TE.push_back((*photFirstTimeE)[0][bestPhotHitIndex1]);
                phot1RhoE.push_back((*photFirstRhoE)[0][bestPhotHitIndex1]);
                phot1ZE.push_back((*photFirstZE)[0][bestPhotHitIndex1]);
            }

            if(phot2IsBarrel){
                phot2EB.push_back((*photEB)[1][bestPhotHitIndex2]);
                phot2TB.push_back((*photFirstTimeB)[1][bestPhotHitIndex2]);
                phot2RhoB.push_back((*photFirstRhoB)[1][bestPhotHitIndex2]);
                phot2ZB.push_back((*photFirstZB)[1][bestPhotHitIndex2]);
            }
            else{
                phot2EE.push_back((*photEE)[1][bestPhotHitIndex2]);
                phot2TE.push_back((*photFirstTimeE)[1][bestPhotHitIndex2]);
                phot2RhoE.push_back((*photFirstRhoE)[1][bestPhotHitIndex2]);
                phot2ZE.push_back((*photFirstZE)[1][bestPhotHitIndex2]);
            }
        }

        else if(option == 2){ //first hit time, mean position
            if(phot1IsBarrel){
                phot1EB.push_back((*photEB)[0][bestPhotHitIndex1]);
                phot1TB.push_back((*photFirstTimeB)[0][bestPhotHitIndex1]);
                phot1RhoB.push_back((*photRhoB)[0][bestPhotHitIndex1]);
                phot1ZB.push_back((*photZB)[0][bestPhotHitIndex1]);
            }
            else{
                phot1EE.push_back((*photEE)[0][bestPhotHitIndex1]);
                phot1TE.push_back((*photFirstTimeE)[0][bestPhotHitIndex1]);
                phot1RhoE.push_back((*photRhoE)[0][bestPhotHitIndex1]);
                phot1ZE.push_back((*photZE)[0][bestPhotHitIndex1]);
            }

            if(phot2IsBarrel){
                phot2EB.push_back((*photEB)[1][bestPhotHitIndex2]);
                phot2TB.push_back((*photFirstTimeB)[1][bestPhotHitIndex2]);
                phot2RhoB.push_back((*photRhoB)[1][bestPhotHitIndex2]);
                phot2ZB.push_back((*photZB)[1][bestPhotHitIndex2]);
            }
            else{
                phot2EE.push_back((*photEE)[1][bestPhotHitIndex2]);
                phot2TE.push_back((*photFirstTimeE)[1][bestPhotHitIndex2]);
                phot2RhoE.push_back((*photRhoE)[1][bestPhotHitIndex2]);
                phot2ZE.push_back((*photZE)[1][bestPhotHitIndex2]);
            }
        }

        else if(option == 3){ //mean time, first hit position
            if(phot1IsBarrel){
                phot1EB.push_back((*photEB)[0][bestPhotHitIndex1]);
                phot1TB.push_back((*photT3x3B)[0][bestPhotHitIndex1]);
                phot1RhoB.push_back((*photFirstRhoB)[0][bestPhotHitIndex1]);
                phot1ZB.push_back((*photFirstZB)[0][bestPhotHitIndex1]);
            }
            else{
                phot1EE.push_back((*photEE)[0][bestPhotHitIndex1]);
                phot1TE.push_back((*photT3x3E)[0][bestPhotHitIndex1]);
                phot1RhoE.push_back((*photFirstRhoE)[0][bestPhotHitIndex1]);
                phot1ZE.push_back((*photFirstZE)[0][bestPhotHitIndex1]);
            }

            if(phot2IsBarrel){
                phot2EB.push_back((*photEB)[1][bestPhotHitIndex2]);
                phot2TB.push_back((*photT3x3B)[1][bestPhotHitIndex2]);
                phot2RhoB.push_back((*photFirstRhoB)[1][bestPhotHitIndex2]);
                phot2ZB.push_back((*photFirstZB)[1][bestPhotHitIndex2]);
            }
            else{
                phot2EE.push_back((*photEE)[1][bestPhotHitIndex2]);
                phot2TE.push_back((*photT3x3E)[1][bestPhotHitIndex2]);
                phot2RhoE.push_back((*photFirstRhoE)[1][bestPhotHitIndex2]);
                phot2ZE.push_back((*photFirstZE)[1][bestPhotHitIndex2]);
            }
        }

        if(option == 4){
            if(phot1IsBarrel){
                phot1EB.push_back((*photEB)[0][bestPhotHitIndex1]);
                phot1TB.push_back((*photTB)[0][bestPhotHitIndex1]);
                phot1RhoB.push_back((*photRhoB)[0][bestPhotHitIndex1]);
                phot1ZB.push_back((*photZB)[0][bestPhotHitIndex1]);
            }
            else{
                phot1EE.push_back((*photEE)[0][bestPhotHitIndex1]);
                phot1TE.push_back((*photTE)[0][bestPhotHitIndex1]);
                phot1RhoE.push_back((*photRhoE)[0][bestPhotHitIndex1]);
                phot1ZE.push_back((*photZE)[0][bestPhotHitIndex1]);
            }

            if(phot2IsBarrel){
                phot2EB.push_back((*photEB)[1][bestPhotHitIndex2]);
                phot2TB.push_back((*photTB)[1][bestPhotHitIndex2]);
                phot2RhoB.push_back((*photRhoB)[1][bestPhotHitIndex2]);
                phot2ZB.push_back((*photZB)[1][bestPhotHitIndex2]);
            }
            else{
                phot2EE.push_back((*photEE)[1][bestPhotHitIndex2]);
                phot2TE.push_back((*photTE)[1][bestPhotHitIndex2]);
                phot2RhoE.push_back((*photRhoE)[1][bestPhotHitIndex2]);
                phot2ZE.push_back((*photZE)[1][bestPhotHitIndex2]);
            }
        }

        else if(option == 5){ //single crystal mean time, first hit position
            if(phot1IsBarrel){
                phot1EB.push_back((*photEB)[0][bestPhotHitIndex1]);
                phot1TB.push_back((*photTB)[0][bestPhotHitIndex1]);
                phot1RhoB.push_back((*photFirstRhoB)[0][bestPhotHitIndex1]);
                phot1ZB.push_back((*photFirstZB)[0][bestPhotHitIndex1]);
            }
            else{
                phot1EE.push_back((*photEE)[0][bestPhotHitIndex1]);
                phot1TE.push_back((*photTE)[0][bestPhotHitIndex1]);
                phot1RhoE.push_back((*photFirstRhoE)[0][bestPhotHitIndex1]);
                phot1ZE.push_back((*photFirstZE)[0][bestPhotHitIndex1]);
            }

            if(phot2IsBarrel){
                phot2EB.push_back((*photEB)[1][bestPhotHitIndex2]);
                phot2TB.push_back((*photTB)[1][bestPhotHitIndex2]);
                phot2RhoB.push_back((*photFirstRhoB)[1][bestPhotHitIndex2]);
                phot2ZB.push_back((*photFirstZB)[1][bestPhotHitIndex2]);
            }
            else{
                phot2EE.push_back((*photEE)[1][bestPhotHitIndex2]);
                phot2TE.push_back((*photTE)[1][bestPhotHitIndex2]);
                phot2RhoE.push_back((*photFirstRhoE)[1][bestPhotHitIndex2]);
                phot2ZE.push_back((*photFirstZE)[1][bestPhotHitIndex2]);
            }
        }


        imVertFuncsPhot1 = makeImVertFuncs(phot1TB, phot1EB, phot1RhoB, phot1ZB, phot1TE, phot1EE, phot1RhoE, phot1ZE, photEThreshold);
        imVertFuncsPhot2 = makeImVertFuncs(phot2TB, phot2EB, phot2RhoB, phot2ZB, phot2TE, phot2EE, phot2RhoE, phot2ZE, photEThreshold);

        //Apply averaging -- sqrt(hit energy)
        //cout << "   Averaging hits..." << endl;
        //AveragedTF1Object *avObj1 = new AveragedTF1Object(imVertFuncsPhot1, (*photEB)[0], (*photEE)[0], photEThreshold);
        AveragedTF1Object *avObj1 = new AveragedTF1Object(imVertFuncsPhot1, phot1EB, phot1EE, photEThreshold);
        avgVertFuncJet1 = new TF1("avgVertFuncJet1", avObj1,-7, 7, 0, "AveragedTF1Object"); 
        AveragedTF1Object *avObj2 = new AveragedTF1Object(imVertFuncsPhot2, phot2EB, phot2EE, photEThreshold);
        //AveragedTF1Object *avObj2 = new AveragedTF1Object(imVertFuncsPhot2, (*photEB)[1], (*photEE)[1], photEThreshold);
        avgVertFuncJet2 = new TF1("avgVertFuncJet2", avObj2,-7, 7, 0, "AveragedTF1Object"); 
        
        //Intersect the two averaged imaginary vertex functions to estimate the vertex position
        //cout << "   Getting estimated vertex position..." << endl;
        TF1 *fint = new TF1("fint",finter,-7, 7, 0);
        if(fint->GetMinimum(-7, 7) > 0.001){
            //cout << "Photons could not be vertexed! (tVertex functions don't intersect)" <<endl;
            continue;
        }
        else photsVertexed++;
        double xint = fint->GetMinimumX(-7, 7);
        double vertexGuess = avgVertFuncJet1->Eval(xint);
        vertexHistPhot->Fill(vertexGuess - vertexZ*10);
        vertexHistPhot2D->Fill(vertexZ*10, vertexGuess);
        //cout << "Photon Vertexing successful.  Estimated vertex Z: " << vertexGuess << ".  Tracker vertex Z: " << vertexZ*10 << endl;
        vertexedPhots = true;
        vertexFromPhots = vertexGuess;

        delete fint;
        delete avObj1;
        delete avObj2;
        delete avgVertFuncJet1;
        delete avgVertFuncJet2;
    }

    if(vertexedJets && vertexedPhots) jetVertexVsPhotVertex->Fill(vertexFromJets, vertexFromPhots);

  }
    
    delete c;
    fOut.cd();
    vertexHist->Write();

    //output resolution estimate to file
    if(vertexHist->GetEntries() > 0){
        vertexHist->Fit("gaus");
        TF1 *gausFit = vertexHist->GetFunction("gaus");
        ofstream resOut("jetVertexRes.txt", std::ofstream::app);
        resOut << minNumHits << " " << hitEThreshold << " " << timeSigma << " " << vertexHist->GetRMS() << " " << gausFit->GetParameter(2) << " " << eventsDone << endl;
    }

    vertexHist2D->Write();
    vertexHistPhot->Write();
    vertexHistPhot2D->Write();
    meanJetTimeEta->Write();
    meanCrystTimeIEta->Write();
    meanCrystTimeIXIY->Write();
    vertexHist2DBarrel->Write();
    vertexHist2DEndcap->Write();
    etaOfVertexedJets->Write();
    ptOfVertexedJets->Write();
    jetVertexVsPhotVertex->Write();
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    vertexHist->Draw();
    c2->Print("vertexHist.pdf");
    vertexHist2D->Draw();
    c2->Print("vertexHist2D.pdf");
    meanJetTimeEta->Draw();
    c2->Print("vertexMeanJetTime.pdf");
    meanCrystTimeIEta->Draw();
    c2->Print("vertexMeanCrystTimeB.pdf");
    meanCrystTimeIXIY->Draw();
    c2->Print("vertexMeanCrystTimeE.pdf");
    meanTimeVsFirstTime->Write();
    meanTimeVsFirstTimeBarrel->Write();
    meanTimeVsFirstTimePhot->Write();
    meanTimeMinusFirstTimeVsEta->Write();
    meanRhoMinusFirstRhoVsEta->Write();
    meanZMinusFirstZVsEta->Write();
    meanTimeMinusFirstTimeEE->Write();
    meanRhoMinusFirstRhoEE->Write();
    meanZMinusFirstZEE->Write();
    cout << "Total number of events: " << ecalTree->GetEntries() << endl;
    cout << "Events with two jets: " << enoughJets << endl;
    cout << "Events with jets in opposite halves of ECAL: " << jetsInOppositeHalves << endl;
    cout << "Events where jets pass cuts: " << jetsPassCuts << endl;
    cout << "Events where jets have enough hits: " << jetsHaveEnoughHits << endl;
    cout << "Events where jets are vertexed: " << jetsVertexed << endl;
    cout << endl << "Events with two photons: " << enoughPhots << endl;
    cout << "Events with photons in opposite halves of ECAL: " << photsInOppositeHalves << endl;
    cout << "Events where photons pass cuts: " << photsPassCuts << endl;
    cout << "Events where photons have enough hits: " << photsHaveEnoughHits << endl;
    cout << "Events where photons are vertexed: " << photsVertexed << endl;
    delete c2;
    fOut.Close();
    infile->Close();    
    return;

}

# ifndef __CINT__
int main(int argc, char *argv[]){
    if(argc < 5){
        cerr << "usage: VertexEcalJets minHits hitThresh timeSigma option <option = 0 for mean time jets, 1 for first hit time jets>" << endl;
        return -1;
    }

    gROOT->SetBatch();
    minNumHits = atoi(argv[1]);
    hitEThreshold = atoi(argv[2]);
    timeSigma = atof(argv[3]);
    int choice = atoi(argv[4]);
    VertexEcalJets(choice);
    return 0;
}
#endif
