// $Id: rootlogon.C,v 1.1 2013/06/09 13:30:16 sixie Exp $

{
  
  if (gSystem->Getenv("CMSSW_VERSION")) {
    
  //  TString rfitpath("/cvmfs/cms.cern.ch/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
  //  path += rfitpath;
      path += " -I/home/djanders/CMSSW_5_3_9/src/fastjet-install/include";
    gSystem->SetIncludePath(path.Data());
 
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
  }

  gSystem->AddIncludePath("-I/home/djanders/CMSSW_5_3_9/src/fastjet-install/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+"/src/");
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_RELEASE_BASE"))+"/src/");
  gInterpreter->AddIncludePath(TString("/home/djanders/CMSSW_5_3_9/src/fastjet-install/include/"));

  gSystem->Load("libCMSAnaDataTree.so");
  gSystem->Load("libCMSAnaUtils.so");
  gSystem->Load("libCMSAnaJetEnergyCorrections.so");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  gSystem->Load("/home/djanders/CMSSW_5_3_9/src/fastjet-install/lib/libfastjet.so");
  gSystem->Load("/home/djanders/CMSSW_5_3_9/src/fastjet-install/lib/libsiscone.so");
  gSystem->Load("/home/djanders/CMSSW_5_3_9/src/fastjet-install/lib/libsiscone_spherical.so");
  gSystem->Load("/home/djanders/CMSSW_5_3_9/src/fastjet-install/lib/libfastjettools.so");
  gSystem->Load("/home/djanders/CMSSW_5_3_9/src/fastjet-install/lib/libfastjetplugins.so");


}
