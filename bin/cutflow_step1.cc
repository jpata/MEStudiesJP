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
#include "TH2D.h"
#include "TH3D.h"
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

template <class T>
T* add_hist_1d(std::map<std::string, TH1*>& histmap, string hname, int b1, int b2) {
    histmap[hname] = new T(hname.c_str(), hname.c_str(), b2-b1, b1, b2);
    histmap[hname]->Sumw2();
    return (T*)histmap[hname];
}

template <class T>
T* add_hist_1d(std::map<std::string, TH1*>& histmap, string hname, double b1, double b2, int nb) {
    histmap[hname] = new T(hname.c_str(), hname.c_str(), nb, b1, b2);
    histmap[hname]->Sumw2();
    return (T*)histmap[hname];
}

template <class T>
T* add_hist_2d(std::map<std::string, TH1*>& histmap, string hname, double b11, double b21, int nb1, double b12, double b22, int nb2) {
    histmap[hname] = new T(hname.c_str(), hname.c_str(), nb1, b11, b21, nb2, b12, b22);
    histmap[hname]->Sumw2();
    return (T*)histmap[hname];
}

template <class T>
T* add_hist_3d(std::map<std::string, TH1*>& histmap, string hname, double b11, double b21, int nb1, double b12, double b22, int nb2, double b13, double b23, int nb3) {
    histmap[hname] = new T(hname.c_str(), hname.c_str(), nb1, b11, b21, nb2, b12, b22, nb3, b13, b23);
    histmap[hname]->Sumw2();
    return (T*)histmap[hname];
}

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
    
    
    TFile* outfile = new TFile("outfile_step1.root", "RECREATE");

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
        
        const std::string pf(sample_nick + "_");
        outfile->cd();
        TH2D* h_lep_pt_iso_mu = add_hist_2d<TH2D>(histmap, pf+"lep_pt_iso_mu", 0, 600, 60, 0, 5, 20);
        TH2D* h_lep_pt_iso_ele = add_hist_2d<TH2D>(histmap, pf+"lep_pt_iso_ele", 0, 600, 60, 0, 5, 20);
        
        TH3D* h_lep_pt_iso_npv_mu = add_hist_3d<TH3D>(histmap, pf+"lep_pt_iso_npv_mu", 0, 600, 60, 0, 5, 20, 0, 50, 50);
        TH3D* h_lep_pt_iso_npv_ele = add_hist_3d<TH3D>(histmap, pf+"lep_pt_iso_npv_ele", 0, 600, 60, 0, 5, 20, 0, 50, 50);
        
        TH3D* h_lep_pt_iso2_npv_mu = add_hist_3d<TH3D>(histmap, pf+"lep_pt_iso2_npv_mu", 0, 600, 60, 0, 0.2, 100, 0, 50, 50);
        TH3D* h_lep_pt_iso2_npv_ele = add_hist_3d<TH3D>(histmap, pf+"lep_pt_iso2_npv_ele", 0, 600, 60, 0, 0.2, 100, 0, 50, 50);
        
        //create ME TTree and branch variables
        TTHTree t((TTree*)(tf->Get("tthNtupleAnalyzer/events")));
        std::cout << sample_fn << " " << sample_nick << " entries " << t.tree->GetEntries() << std::endl;
       	 
        t.set_branch_addresses();
       	t.tree->SetBranchStatus("*", false); 
       	t.tree->SetBranchStatus("n__lep", true); 
       	t.tree->SetBranchStatus("lep__pt", true); 
       	t.tree->SetBranchStatus("lep__iso", true); 
       	t.tree->SetBranchStatus("lep__rel_iso", true); 
       	t.tree->SetBranchStatus("lep__id", true); 
       	t.tree->SetBranchStatus("n__pv", true); 
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
            
            float npv = t.n__pv;
            
            for(int nl=0; nl < nlep; nl++) {
                
                float pt = t.lep__pt[nl];
                int id = t.lep__id[nl];
                float iso = t.lep__rel_iso[nl];
                
                if (std::abs(id) == 13) {
                    h_lep_pt_iso_mu->Fill(pt, iso, 1.0);
                    h_lep_pt_iso_npv_mu->Fill(pt, iso, npv, 1.0);
                    h_lep_pt_iso2_npv_mu->Fill(pt, iso, npv, 1.0);
                } else if (std::abs(id) == 11) {
                    h_lep_pt_iso_ele->Fill(pt, iso, 1.0);
                    h_lep_pt_iso_npv_ele->Fill(pt, iso, npv, 1.0);
                    h_lep_pt_iso2_npv_ele->Fill(pt, iso, npv, 1.0);
                }
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
