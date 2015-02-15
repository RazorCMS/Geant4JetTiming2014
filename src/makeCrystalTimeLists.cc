#include <iostream>
#include <fstream>
#include <map>
#include <sstream>

using namespace std;

//Opens the lists of ROOT files and adds them to the "good" list if their index appears in both lists.  
int main(){
    ifstream bIn("baconList.txt");
    string tmp;
    int index;
    map<int, string> fileMap;
    cout << "Reading ntuple files..." << endl;
    int dups = 0;
    while(getline(bIn, tmp)){
        string tmp2 = tmp;
        tmp2.erase(0, 66);
        tmp2.erase(tmp2.length() - 11);
        istringstream(tmp2) >> index;
        if(index == 1138 || index == 1435) continue; //bad files
        if(fileMap.count(index) > 0){
            cout << "File " << index << " is a duplicate!" << endl;
            dups++;
            fileMap.erase(index); //don't risk it
        }
        else fileMap[index] = tmp;
    }
    bIn.close();

    cout << "Reading geant files..." << endl;
    ifstream gIn("geantList.txt");
    ofstream bOut("goodBaconList.txt");
    ofstream gOut("goodGeantList.txt");

    int missed = 0;
    int dups2 = 0;
    int prevIndex = 0;
    while(getline(gIn, tmp)){
        prevIndex = index;
        string tmp2 = tmp;
        tmp2.erase(0, 61);
        tmp2.erase(tmp2.length() - 11);
        istringstream(tmp2) >> index;
        if(index == prevIndex){
            cout << "File " << index << " is a duplicate!" << endl;
            dups2++;
            continue;
        }
        if(fileMap.count(index) > 0){
            bOut << fileMap[index] << endl;
            gOut << tmp << endl;
            fileMap[index] = "checked!";
        }
        else{
            cout << "Missed file " << index << "!" << endl;
            missed++;
        }
    }
    cout << "Checking for missed events in geant files..." << endl;
    int missed2 = 0;
    for(map<int, string>::iterator it = fileMap.begin(); it != fileMap.end(); ++it){
        if(it->second != "checked!"){
            cout << "Missed file " << it->first << "!" << endl;
            missed2++;
        }
    }
        gIn.close();
        bOut.close();
        gOut.close();
        cout << "Missed: in event files: " << missed << endl;
        cout << "Missed in geant files: " << missed2 << endl;
        cout << "Duplicate event files: " << dups << endl;
        cout << "Duplicate geant files: " << dups2 << endl;

        return 0;
}
