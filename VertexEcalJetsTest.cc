#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<int> >+;
#endif

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
#include <TRandom3.h>                
#include <algorithm>

#include "CMSAna/Utils/CommonTools.hh"

// FastJet Stuff
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

using namespace std;

const int minNumHits = 6; //minimum number of energetic hits required in each jet
const double hitEThreshold = 4; //in GeV
const double speedoflight = 299.79245800; //in mm/ns
const double jPtCut = 80;

//don't show this code to a computer scientist.

//dummy class for making the averaged imaginary time vertex function
class  AveragedTF1Object {
public:
    AveragedTF1Object(vector<TF1*> funcs, vector<double> barrelEnergy, vector<double> endcapEnergy){
        theFuncs = funcs;
        eBarrel = barrelEnergy;
        eEndcap = endcapEnergy;
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
        for(int i = 0; i < eBarrel.size(); i++){
            if(!isnan(theFuncs[i]->Eval(x[0]))){
                funcTotalWeight += sqrt(eBarrel[i]);
                funcWeightedSum += sqrt(eBarrel[i])*theFuncs[i]->Eval(x[0]);
                isNaN = false;
            }
        }
        for(int i = 0; i < eEndcap.size(); i++){
            if(!isnan(theFuncs[i+eBarrel.size()]->Eval(x[0]))){
                funcTotalWeight += sqrt(eEndcap[i]);
                funcWeightedSum += sqrt(eEndcap[i])*theFuncs[i + eBarrel.size()]->Eval(x[0]);
                isNaN = false;
            }
        }
        if(isNaN) return nan("");
        double avgFuncVal = funcWeightedSum / funcTotalWeight;
        return avgFuncVal;
     
    }

    vector<TF1*> theFuncs;
    vector<double> eBarrel;
    vector<double> eEndcap;
};

//dummy function to find the intersection between two TF1 objects
//(why is ROOT so terrible?)
TF1 *avgVertFuncJet1, *avgVertFuncJet2;
double finter(double *x, double*par) {
   return TMath::Abs(avgVertFuncJet1->EvalPar(x,par) - avgVertFuncJet2->EvalPar(x,par));
}

//gets the number of hits in a jet having energy higher than some threshold
int numHitsAboveThreshold(vector<double> hitEnergiesB, vector<double> hitEnergiesE){
    int numHits = 0;
    for(int i = 0; i < hitEnergiesB.size(); i++){
        if(hitEnergiesB[i] >= hitEThreshold) numHits++;
    }
    for(int i = 0; i < hitEnergiesE.size(); i++){
        if(hitEnergiesE[i] >= hitEThreshold) numHits++;
    }
    return numHits;
}

//gets the z position on the beamline corresponding to the time of the rechit
double tToZ(double t, double hitTime, double hitRho, double hitZ){
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

//Returns a vector of TF1's corresponding to the imaginary vertex functions for the hits in the crystal
vector<TF1*> makeImVertFuncs(vector<double> crystTB, vector<double> rhoB, vector<double> hitZB, vector<double> crystEB, vector<double> crystTE, vector<double> rhoE, vector<double> hitZE, vector<double> crystEE){
    vector<TF1*> theFuncs;
    /*//loop over barrel
    for(int crystI = 0; crystI < crystTB.size(); crystI++){
        TF1 *f = new TF1("f","tToZ(x,[0],[1],[2])",-7 , 7);
        f->SetParameters(crystTB[crystI], rhoB[crystI], hitZB[crystI]);
        theFuncs.push_back(f);
    }
    //loop over endcap
    for(int crystI = 0; crystI < crystTE.size(); crystI++){
        TF1 *f = new TF1("f", "tToZ(x, [0], [1], [2])", -7 ,  7);
        f->SetParameters(crystTE[crystI], rhoE[crystI], hitZE[crystI]);
        theFuncs.push_back(f);
    }*/
    //get most energetic hit
    TF1 *f = new TF1("f", "tToZ(x, [0], [1], [2])", -7, 7);
    double bestT = 0;
    double bestRho = 0;
    double bestZ = 0;
    double bestEB = 0;
    double bestIndexB = -1;
    double bestEE = 0;
    double bestIndexE = -1;
    for(int i = 0; i < crystEB.size(); i++){
        if(crystEB[i] > bestEB){
            bestEB = crystEB[i];
            bestIndexB = i;
        }
    }
    for(int i = 0; i < crystEE.size(); i++){
        if(crystEE[i] > bestEE){
            bestEE = crystEE[i];
            bestIndexE = i;
        }

    cout << *std::max_element(crystEB.begin(), crystEB.end()) << endl;
    cout << *largestEB << " " << *largestEE << endl;
   if(*largestEB > *largestEE){
        int bestIndex = largestEB - crystEB.begin();
        bestT = crystTB[bestIndex];
        bestRho = rhoB[bestIndex];
        bestZ = hitZB[bestIndex];
    }
    else{
        int bestIndex = largestEE - crystEE.begin();
        bestT = crystTE[bestIndex];
        bestRho = rhoE[bestIndex];
        bestZ = hitZE[bestIndex];
    }
    f->SetParameters(bestT, bestRho, bestZ);

    theFuncs.push_back(f);
        
    cout << "didn't crasy there" << endl;
    return theFuncs;
}

void VertexEcalJetsTest() {
  //gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  //input file
  TFile *infile = new TFile("CrystalTimeInfo.root");
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

  
  TFile fOut("dummy.root", "recreate");
  TH1D *vertexHist = new TH1D("vertexHist", "; vertex z (timing) - vertex z (tracker); Number of events; Distance between reconstructed vertex and tracker vertex", 50, -500, 500);
  TH2D *vertexHist2D = new TH2D("vertexHist2D", ";vertex z (tracker); vertex z (timing); Number of events; Vertexing correlation plot", 50, -500, 500, 50, -500, 500);
  TH2D *meanJetTimeEta = new TH2D("meanJetTimeEta", ";jet #eta; mean time of jet hits; Number of events; Mean time of jet hits", 60, -3, 3, 100, -5, 5);
  TH2D *meanCrystTimeIEta = new TH2D("meanCrystTimeIEta", ";crystal i#eta; mean time of hits; Number of events; Mean time of crystal hits", 171, -85, 86, 200, -20, 20);
  //TH3D *meanCrystTimeIXIY = new TH3D("meanCrystTimeIXIY", ";ix; iy; mean time of hits; Mean time of crystal hits", 100, 1, 101, 100, 1, 101, 100, -20, 20);
  TH2D *meanCrystTimeIXIY = new TH2D("meanCrystTimeIXIY", ";#sqrt{ix^{2} + iy^{2}}; mean time of hits; Number of Events; Mean time of crystal hits in endcap", 145, 0, 145, 200, -20, 20);

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  int eventsDone = 0;
  bool jetdone = false;
  //LOOP OVER EVENTS
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
      
    ecalTree->GetEntry(ientry);
    cout << "\n\nEvent " << ientry << "\n";

    if(jEcalE->size() < 1) cout << "No jets." << endl;
    else if(jEcalE->size() == 1) cout << "Only one jet." << endl;
    else if(jEcalE->size() >= 2){ //if we have two jets, try to vertex them
        //Note that jets are sorted in decreasing order by pT, and we will just try to vertex the two most energetic jets for now.  

        //Get mean time of jets
        cout << "Getting jet times..." << endl;
        for(int i = 0; i < crystTB->size(); i++){
            double meanJetTime = 0;
            double totalJetE = 0;
            double meanJetRho = 0;
            double meanJetZ = 0;
            for(int c = 0; c < crystTB->at(i).size(); c++){
                if((*crystEB)[i][c] >= hitEThreshold){
                    totalJetE += (*crystEB)[i][c];
                    meanJetTime += (*crystEB)[i][c]*(*crystTB)[i][c];
                    meanJetRho += (*crystEB)[i][c]*(*rhoB)[i][c];
                    meanJetZ += (*crystEB)[i][c]*(*hitZB)[i][c];
                    meanCrystTimeIEta->Fill((*crystIEta)[i][c], (*crystTB)[i][c] - sqrt((*hitZB)[i][c]*(*hitZB)[i][c] + (*rhoB)[i][c]*(*rhoB)[i][c])/speedoflight, (*crystEB)[i][c]);
                }
            }
            for(int c = 0; c < crystTE->at(i).size(); c++){
                if((*crystEE)[i][c] >= hitEThreshold){
                    totalJetE += (*crystEE)[i][c];
                    meanJetTime += (*crystEE)[i][c]*(*crystTE)[i][c];
                    meanJetRho += (*crystEE)[i][c]*(*rhoE)[i][c];
                    meanJetZ += (*crystEE)[i][c]*(*hitZE)[i][c];
                    //meanCrystTimeIXIY->Fill((*crystIX)[i][c], (*crystIY)[i][c], (*crystTE)[i][c] - sqrt((*hitZE)[i][c]*(*hitZE)[i][c] + (*rhoE)[i][c]*(*rhoE)[i][c])/speedoflight, (*crystEE)[i][c]);
                    meanCrystTimeIXIY->Fill(sqrt((*crystIX)[i][c]*(*crystIX)[i][c] + (*crystIY)[i][c]*(*crystIY)[i][c]), (*crystTE)[i][c] - sqrt((*hitZE)[i][c]*(*hitZE)[i][c] + (*rhoE)[i][c]*(*rhoE)[i][c])/speedoflight, (*crystEE)[i][c]);
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
            cout << "Jets are in the same half of the ECAL." << endl;
            continue;
        }

        /*//Check if the jet has enough hits (decide what is enough)
        int numHitsJet1 = numHitsAboveThreshold((*crystEB)[0], (*crystEE)[0]);
        int numHitsJet2 = numHitsAboveThreshold((*crystEB)[1], (*crystEE)[1]);
        if(numHitsJet1 < minNumHits || numHitsJet2 < minNumHits){
            cout << "At least one jet has not enough hits." << endl;
            continue;
        }*/

        cout << "Vertexing jets..." << endl;
        //TODO: Calibrate times so that times in each crystal have mean zero for particles coming from (0,0,0)
        
        //Make imaginary vertex function for each hit -- use the 3x3 energies and times
        cout << "   Making imaginary vertex functions..." << endl;
        vector<TF1*> imVertFuncsJet1;
        vector<TF1*> imVertFuncsJet2;
        imVertFuncsJet1 = makeImVertFuncs((*crystT3x3B)[0], (*rhoB)[0], (*hitZB)[0], (*crystEB)[0], (*crystT3x3E)[0], (*rhoE)[0], (*hitZE)[0], (*crystEE)[0]);
        imVertFuncsJet2 = makeImVertFuncs((*crystT3x3B)[1], (*rhoB)[1], (*hitZB)[1], (*crystEB)[1], (*crystT3x3E)[1], (*rhoE)[1], (*hitZE)[1], (*crystEE)[1]);
        cout << imVertFuncsJet1.size() << " " << imVertFuncsJet2.size() << endl;

        //Apply averaging -- sqrt(hit energy)
        cout << "   Averaging hits..." << endl;
        AveragedTF1Object *avObj1 = new AveragedTF1Object(imVertFuncsJet1, (*crystEB)[0], (*crystEE)[0]);
        avgVertFuncJet1 = new TF1("avgVertFuncJet1", avObj1,-7, 7, 0, "AveragedTF1Object"); 
        AveragedTF1Object *avObj2 = new AveragedTF1Object(imVertFuncsJet2, (*crystEB)[1], (*crystEE)[1]);
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

            //draw averaged vertex functions only
            /*TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
            avgVertFuncJet1->Draw();
            avgVertFuncJet2->Draw("same");
            fOut.cd();
            fOut.WriteTObject(c1, c1->GetName(), "Write");
            c1->Print("vertexJetDiagram.pdf");
            delete c1;*/
            //jetdone = true;
        }

        //TODO: Sign correction (what is this exactly)
        //
        //TODO: Slope correction (afterwards)
        //
        //TODO: Beam spot correction
        
        //Intersect the two averaged imaginary vertex functions to estimate the vertex position
        cout << "   Getting estimated vertex position..." << endl;
        TF1 *fint = new TF1("fint",finter,0,10,0);
        if(fabs(fint->GetMinimum()) > 0.001){
            cout << "Jets could not be vertexed! (tVertex functions don't intersect)" <<endl;
            continue;
        }
        double xint = fint->GetMinimumX();
        double vertexGuess = avgVertFuncJet1->Eval(xint);
        vertexHist->Fill(vertexGuess - vertexZ*10);
        vertexHist2D->Fill(vertexZ*10, vertexGuess);
        cout << "Vertexing successful.  Estimated vertex Z: " << vertexGuess << ".  Tracker vertex Z: " << vertexZ*10 << endl;
        eventsDone++;

        delete fint;
        delete avObj1;
        delete avObj2;
        delete avgVertFuncJet1;
        delete avgVertFuncJet2;
    }
  
  }
    
    delete c;
    fOut.cd();
    vertexHist->Write();
    vertexHist2D->Write();
    meanJetTimeEta->Write();
    meanCrystTimeIEta->Write();
    meanCrystTimeIXIY->Write();
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
    delete c2;
    fOut.Close();
    infile->Close();    
    return;

}

# ifndef __CINT__
int main(){
        gROOT->SetBatch();
        VertexEcalJetsTest();
        return 0;
}
#endif
