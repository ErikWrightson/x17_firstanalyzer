/**
 * Header file for MassFinder class that is used for finding the invariant mass and other distributions before and after each cut for X17.
 *
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 07.13.2026
 * @creation 07.13.2026
 */

#ifndef MassFinder_H
#define MassFinder_H

//ROOT Includes that may be handy to have.
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TColor.h>
#include <TString.h>
//Needed for reading in vector types from root files.
#include <TInterpreter.h>
#include <TPrincipal.h>
#include <TVectorD.h>
#include <PxPyPzE4D.h>
#include "ROOT/TThreadedObject.hxx"
#include "Physics_Utils.h"
#include "Utils.h"

//Generally useful includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <iomanip>

using namespace std;
using namespace Physics_Utils;

class MassFinder{

    public:
        TChain* chain;


        static constexpr Int_t CUT_NUM = 11;
        static constexpr TString XCUT_NAME[CUT_NAME] = {"None", "3 or More", "Timing", "Fiducial", "Cluster E", "E Sum", "Coplanarity", "X_E", "1 GEM Match", "2 GEM Match", "Vertex Z"}

        static constexpr Int_t  MAX_CLUSTERS = 400;      //Maximum number of clusters.
        static constexpr Int_t  MAX_GEMS = 4;            //Maximum number of GEMs.

        static constexpr Int_t MAX_VTP_CLUSTERS = 100;    //Maximum Number of clusters found by the VTP for the online clustering trigger.

        static constexpr Double_t HC_XMIN = -600.0;
        static constexpr Double_t HC_XMAX = 600.0;
        static constexpr Double_t HC_XBINS = 1200.0;

        static constexpr Double_t HC_YMIN = -600.0;
        static constexpr Double_t HC_YMAX = 600.0;
        static constexpr Double_t HC_YBINS = 1200.0;

        static constexpr Double_t MIN_TH = 0.0;
        static constexpr Double_t MAX_TH = 5.0;
        static constexpr Double_t TH_BINS = 50.0;

        static constexpr Double_t MIN_INVM = 0.0;
        static constexpr Double_t MAX_INVM = 70.0;
        Double_t INVM_BINS = (MAX_INVM - MIN_INVM)*3.0;

        static constexpr Double_t MIN_PT = 0.0;
        static constexpr Double_t MAX_PT = 40.0;
        Double_t PT_BINS = (MAX_PT - MIN_PT)*10.0;

        static constexpr Double_t MIN_PHI = 150.0;
        static constexpr Double_t MAX_PHI = 210.0;
        Double_t PHI_BINS = (MAX_PHI - MIN_PHI)*10.0;

        static constexpr Double_t MIN_TIME = 150.0;
        static constexpr Double_t MAX_TIME = 210.0;
        Double_t TIME_BINS = (MAX_TIME - MIN_TIME);

        static constexpr Double_t MIN_VZ = -2000.0;
        static constexpr Double_t MAXVZ = 7000.0;
        Double_t VZ_BINS = (MAX_VZ - MIN_VZ)/10.0;

        //Locations to store the event data
        UInt_t eventNum;
        UChar_t trigType;
        UInt_t trigger_bits;
        Long64_t time;
        vector<unsigned int> sspRawBuf;
        vector<unsigned int>* sspRawPtr;
        PRadTrigger trig;

        //Vectors containing the vtp bank information from the clustering triggers.
        vector<unsigned int> vtp_roc_tagsBuf;
        vector<unsigned int>* vtp_roc_tagsPtr;

        vector<unsigned int> vtp_nwordsBuf;
        vector<unsigned int>* vtp_nwordsPtr;
        
        vector<unsigned int> vtp_wordsBuf;
        vector<unsigned int>* vtp_wordsPtr;

        //Decoded VTP cluster information
        Int_t vtp_cl_n;
        UShort_t vtp_cl_time[MAX_VTP_CLUSTERS];
        UShort_t vtp_cl_E[MAX_VTP_CLUSTERS];
        UShort_t vtp_cl_center[MAX_VTP_CLUSTERS];
        UChar_t vtp_cl_nblocks[MAX_VTP_CLUSTERS];

        //---------------------------------------------------------------
        //----RECON Tree Branch Variables----
        //---------------------------------------------------------------
        Int_t nClust;
        Float_t cl_x[MAX_CLUSTERS];
        Float_t cl_y[MAX_CLUSTERS];
        Float_t cl_z[MAX_CLUSTERS];
        Float_t cl_E[MAX_CLUSTERS];
        UChar_t cl_nblocks[MAX_CLUSTERS];
        UShort_t cl_center[MAX_CLUSTERS];
        Float_t cl_time[MAX_CLUSTERS];
        UInt_t cl_flag[MAX_CLUSTERS];

        //GEM Matching Information - Useful for event selection for trigger function on specific event types.
        UInt_t match_flag[MAX_CLUSTERS];
        Float_t mgx[MAX_CLUSTERS][MAX_GEMS];
        Float_t mgy[MAX_CLUSTERS][MAX_GEMS];
        Float_t mgz[MAX_CLUSTERS][MAX_GEMS];

        //GEM Quick two layer match
        Float_t matchGEMx[MAX_CLUSTERS][2];
        Float_t matchGEMy[MAX_CLUSTERS][2];
        Float_t matchGEMz[MAX_CLUSTERS][2];

        MassFinder(TChain* c);
        void search_events_electrons();
        void save_histos(TString rootFile);

    private:

        Int_t runNum;
        Double_t Ebeam;
        Long64_t entries;

        //Relevant values of detected particles.
        Double_t theta[3];
        Double_t phi[3];
        DirVector dir[3];
        PxPyPzE4D particle_e[3];
        vector<Double_t> partEnergies;
        
        //Relevant values of potential X particles. 
        PxPyPzE4D X[3];
        Double_t X_th[3];
        Double_t X_phi[3];

        //X histograms
        TH2F* h_X_HC_XY[CUT_NUM];
        TH2F* h_X_E_theta[CUT_NUM];
        TH1F* h_X_invM[CUT_NUM];
        TH1F* h_X_pt[CUT_NUM];
        TH2F* h_X_pxVpy[CUT_NUM];
        TH1F* h_X_vZ[CUT_NUM][4];
        TH1F* h_X_diffPhi[CUT_NUM];
        TH1F* h_X_timing[CUT_NUM];

        void setup_X_histos();
        void fillCutHistos(Int_t c);

};

#endif