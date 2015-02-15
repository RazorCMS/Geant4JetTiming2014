#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
process = "CrystalTimeInfoFiles"
dataset = "Hgg"

inputlist =  "list_BACON_FINAL.txt"
inputlist2 = "list_g4Hits_FINAL.txt"

#Get grid proxy
mydir = "/home/djanders/CMSSW_5_3_9/src/"
myproxy = "x509up_u2619"
os.system("voms-proxy-init -voms cms")
if os.path.isfile("/tmp/"+myproxy):
    print "Proxy file found in /tmp."
else:
    print "Couldn't find the correct proxy in /tmp."
    exit()
os.system("cp /tmp/"+myproxy+" "+mydir)

output = dataset
ijobmax = 949
queue = "all.q@compute-3-10.local,all.q@compute-3-12.local,all.q@compute-3-2.local,all.q@compute-3-3.local,all.q@compute-3-4.local,all.q@compute-3-5.local,all.q@compute-3-9.local"
submitScript = open("submitAll.sh", 'w')
submitScript.write('#!/bin/sh\n')

#there are 2849 files total

################################################
os.system("mkdir -p "+process+"/"+output)
os.system("mkdir -p "+process+"/"+output+"/log/")
os.system("mkdir -p "+process+"/"+output+"/input/")
os.system("mkdir -p "+process+"/"+output+"/src/")
os.system("mkdir -p "+process+"/"+output+"/out/")
#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
input2 = open(inputlist2)
######################################

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = pwd+"/"+process+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    inputfilename2 = pwd+"/"+process+"/"+output+"/input/input2_"+str(ijob)+".list"
    inputfile2 = open(inputfilename2, 'w')
    # if it is a normal job get filesperjob lines
    if ijob != (ijobmax-1):
        for line in range(filesperjob):
            ntpfile = input.readline() 
            inputfile.write(ntpfile)
            ntpfile2 = input2.readline()
            inputfile2.write(ntpfile2)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        ntpfile2 = input2.readline()
        while ntpfile != '':
            inputfile.write(ntpfile)
            inputfile2.write(ntpfile2)
            ntpfile = input.readline()
            ntpfile2 = input2.readline()
            continue
    inputfile.close()
    inputfile2.close()

    # prepare the script to run
    outputname = process+"/"+output+"/src/submit_"+str(ijob)+".sge"
    print "OUTPUTFILE: ", outputname
    basedir = pwd+"/"+process+"/"+output+"/log/";
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/sh\n')
    outputfile.write('#$ -S /bin/sh\n')
    outputfile.write('#$ -V\n')
    outputfile.write('#$ -q '+queue+'\n')
    outputfile.write('#$ -o '+basedir+str(ijob)+'.out -e '+basedir+str(ijob)+'.err\n\n')
    outputfile.write('echo halloo\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc462\n')
    outputfile.write('export PATH=/cms/sw/bin:/cms/swslc5/bin:${PATH}\n')
    outputfile.write('export CMS_PATH=/cms/sw\n')
    outputfile.write('cp '+mydir+myproxy+' /tmp/\n')
    outputfile.write('cd /home/djanders/CMSSW_5_3_9/src/\n')
    outputfile.write('source /home/djanders/cmsset_default.sh\n')
    outputfile.write('cmsenv\n')
    outputfile.write('g++ -g GetCrystalTimeInfo.cc /home/djanders/CMSSW_5_3_9/src/Dict.o -o GetCrystalTimeInfo'+str(ijob)+' `/home/djanders/CMSSW_5_3_9/src/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --cflags` `root-config --glibs` -lRooFitCore -lRooFit -lMinuit -lFoam -lRooStats -I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/fastjet/3.0.1-cms2/include -I$ROOFITSYS/include/ -L$ROOFITSYS/lib/ -I$CMSSW_BASE/src/ -I$CMSSW_RELEASE_BASE/src/ -L/home/djanders/CMSSW_5_3_9/lib/slc5_amd64_gcc462/ -lCMSAnaDataTree -lGenVector\n')
    outputfile.write('./GetCrystalTimeInfo'+str(ijob)+' '+inputfilename+' '+inputfilename2+"\n")
    outputfile.write('rm GetCrystalTimeInfo'+str(ijob)+'\n')
    outputfile.write('rm /tmp/'+myproxy+'\n')
    outputfile.close
    submitScript.write("sleep 1; qsub "+pwd+"/"+outputname+"\n")
    ijob = ijob+1
    continue
