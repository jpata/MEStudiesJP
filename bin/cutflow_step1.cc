#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMatrixT.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTH/TTHNtupleAnalyzer/interface/HypoEnums.hh"
#include "TTH/TTHNtupleAnalyzer/interface/tth_tree.hh"

using namespace std;

int main(int argc, const char* argv[])
{
    gROOT->SetBatch(true);
    
    gSystem->Load("libFWCoreFWLite");
    gSystem->Load("libDataFormatsFWLite");
    
    AutoLibraryLoader::enable();
    
    PythonProcessDesc builder(argv[1]);
    const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->
    getParameter<edm::ParameterSet>("fwliteInput");
    
    //list of input samples
    const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples");
    
    //limits with [first, last] events to process, indexed by full list of samples
    const std::vector<int> ev_limits = in.getParameter<std::vector<int>>("evLimits");
    
    std::map<std::string, TH1*> histmap;
    
    TStopwatch sw;
    
    long n_total_entries = 0;
    double tottime = 0.0;
    
    
    TFile* outfile = new TFile("outfile.root", "RECREATE");
    outfile->cd();

    for(auto& sample : samples ) {
        
        //LFN of sample to read
        const string sample_fn = sample.getParameter<string>("fileName");
        
        //nickname of sample (must be unique)
        const string sample_nick = sample.getParameter<string>("nickName");
        
        //sample type, which may affect the meaning/contents of the TTrees
        //const SampleType sample_type = static_cast<SampleType>(sample.getParameter<int>("type"));
        
        TFile* tf = new TFile(sample_fn.c_str());
        if (tf==0 || tf->IsZombie()) {
            std::cerr << "ERROR: could not open file " << sample_fn << " " << tf << std::endl;
            throw std::exception();
        }
        
        //create ME TTree and branch variables
        TTHTree t((TTree*)(tf->Get("tthNtupleAnalyzer/events")));
        std::cout << sample_fn << " " << sample_nick << " entries " << t.tree->GetEntries() << std::endl;
        
        t.set_branch_addresses();
        
        //count number of bytes read
        long nbytes = 0;
        for (int i = 0 ;i < t.tree->GetEntries(); i++) {
            n_total_entries += 1;
            
            if(n_total_entries % 1000000 == 0) {
                sw.Stop();
                float t = sw.CpuTime();
                tottime += t;
                std::cout << n_total_entries << " " << t << " "
                << tottime << std::endl;
                sw.Start();
            }
            
            //check if we are within limits
            if ((n_total_entries<ev_limits[0]) ||
                (ev_limits[1]>ev_limits[0] && n_total_entries > ev_limits[1])
                ) {
                std::cout << "stopping loop with " << n_total_entries << " "
                << i << " in " << sample_nick << std::endl;
                break;
            }
            
            //zero all branch variables
            t.loop_initialize();
            nbytes += t.tree->GetEntry(i);
            
            int nlep = t.n__lep;
            for(int nl=0; nl < nlep; nl++) {
                float pt = t.lep__pt[nl];
                float iso = t.lep__rel_iso[nl];
            }
        }
        
        std::cout << "read " << nbytes/1024/1024 << " Mb" << std::endl;
        std::cout << "processed " << n_total_entries << " in "
        << tottime << " seconds" << std::endl ;
        tf->Close();
        
    }
    
    for(auto& kv : histmap) {
        const std::string& k = kv.first;
        TH1* h = kv.second;
        h->SetDirectory(outfile);
        cout << "saving " << k << " " << h->Integral() << endl;
    }
    
    outfile->Write();
    outfile->Close();
    
    
}