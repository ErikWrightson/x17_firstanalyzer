/**
 * Header file for MassFinder class that is used for finding the invariant mass and other distributions before and after each cut for X17.
 *
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 07.22.2026
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
#include <Math/GenVector/PxPyPzE4D.h>
#include <Math/Vector4D.h>
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
using ROOT::Math::PxPyPzE4D;

class MassFinder{

    public:
        TChain* chain;

        static constexpr Double_t sigma_E = 5.0; // the sigma level of the energy cuts.

        static constexpr Int_t XCUT_NUM = 11;
        static const TString XCUT_NAME[XCUT_NUM];

        static constexpr Int_t MCUT_NUM = 8;
        static const TString MCUT_NAME[MCUT_NUM];

        static constexpr Int_t XGGCUT_NUM = 12;
        static const TString XGGCUT_NAME[XGGCUT_NUM];


        static constexpr Int_t  MAX_CLUSTERS = 10000;      //Maximum number of clusters.
        static constexpr Int_t  MAX_GEMS = 4;            //Maximum number of GEMs.

        static constexpr Int_t MAX_VTP_CLUSTERS = 100;    //Maximum Number of clusters found by the VTP for the online clustering trigger.

        static constexpr Double_t HC_XMIN = -400.0;
        static constexpr Double_t HC_XMAX = 400.0;
        static constexpr Double_t HC_XBINS = (HC_XMAX - HC_XMIN)/4.0;

        static constexpr Double_t HC_YMIN = -400.0;
        static constexpr Double_t HC_YMAX = 400.0;
        static constexpr Double_t HC_YBINS = (HC_YMAX - HC_YMIN)/4.0;

        static constexpr Double_t MIN_TH = 0.0;
        static constexpr Double_t MAX_TH = 5.0;
        static constexpr Double_t TH_BINS = 50.0;

        static constexpr Double_t MIN_INVM = 0.0;
        static constexpr Double_t MAX_INVM = 60.0;
        Double_t INVM_BINS = (MAX_INVM - MIN_INVM)*2.0;

        static constexpr Double_t MIN_PT = 0.0;
        static constexpr Double_t MAX_PT = 40.0;
        Double_t PT_BINS = (MAX_PT - MIN_PT)*5.0;

        static constexpr Double_t MIN_PHI = 0.0;
        static constexpr Double_t MAX_PHI = 360.0;
        Double_t PHI_BINS = (MAX_PHI - MIN_PHI)*10.0;

        static constexpr Double_t MIN_TIME = 0.0;
        static constexpr Double_t MAX_TIME = 32.0;
        Double_t TIME_BINS = (MAX_TIME - MIN_TIME);

        static constexpr Double_t MIN_VZ = -2000.0;
        static constexpr Double_t MAX_VZ = 9000.0;
        Double_t VZ_BINS = (MAX_VZ - MIN_VZ)/10.0;

        //Locations to store the event data
        UInt_t eventNum;
        UChar_t trigType;
        UInt_t trigger_bits;
        Long64_t time;
        vector<unsigned int> sspRawBuf;
        vector<unsigned int>* sspRawPtr;

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

        MassFinder(TChain* c, Int_t rN);
        void search_events();
        void save_histos(TString rootFile);
        void delete_histos();

    private:

        Int_t runNum;
        Double_t Ebeam;
        Long64_t entries;

        //Relevant values of detected particles.
        Double_t theta[3];
        Double_t phi[3];
        DirVector dir[3];
        LorentzV particle_e[3];
        LorentzV particle_g[3];
        Int_t layer[3];
        vector<Double_t> partEnergies;
        Double_t numMatches[3];
        Double_t vZ[3];
        Utils::Point close[3];
        bool oneGEM[3];
        bool twoGEM[3];
        Int_t cut;
        Int_t cut_gg;
        Double_t maxT;
        Double_t minT;
        Double_t dT;
        
        //Relevant values of potential X particles. 
        LorentzV X[3];
        Double_t X_th[3];
        Double_t X_phi[3];
        Double_t phi_dif[3];
        vector<Int_t> X_eli;

        //Relevant values for the sum of all 3 candidate particles if all are e+ or e-
        Double_t sumRes;
        LorentzV sum;

        //Relevant values of potential X particles with gamma gamma decay
        LorentzV Xgg[3];
        Double_t Xgg_th[3];
        Double_t Xgg_phi[3];
        Double_t phi_dif_gg[3];
        bool Xgg_match[3];
        vector<Int_t> Xgg_eli;

        //Relevant values for the sum of all 3 candidate particles with gamma gamma decay
        Double_t sumResGG[3];
        LorentzV sumGG[3];

        

        LorentzV Mol;
        Double_t MsumRes;
        Double_t dT_Mol;
        Double_t phi_dif_Mol;
        Double_t exp_E[2];

        //X histograms
        TH2F* h_X_HC_XY[XCUT_NUM];
        TH2F* h_X_E_theta[XCUT_NUM];
        TH1F* h_X_invM[XCUT_NUM];
        TH1F* h_X_Sum_pt[XCUT_NUM];
        TH2F* h_X_Sum_pxVpy[XCUT_NUM];
        TH1F* h_X_vZ[XCUT_NUM][4];
        TH1F* h_X_diffPhi[XCUT_NUM];
        TH1F* h_X_timing[XCUT_NUM];
        TH1F* h_X_sumE[XCUT_NUM];
        TH1F* h_X_minE[XCUT_NUM];
        TH1F* h_X_medE[XCUT_NUM];
        TH1F* h_X_maxE[XCUT_NUM];

        //X histograms
        TH2F* h_Xgg_HC_XY[XGGCUT_NUM];
        TH2F* h_Xgg_E_theta[XGGCUT_NUM];
        TH1F* h_Xgg_invM[XGGCUT_NUM];
        TH1F* h_Xgg_Sum_pt[XGGCUT_NUM];
        TH2F* h_Xgg_Sum_pxVpy[XGGCUT_NUM];
        TH1F* h_Xgg_vZ[XGGCUT_NUM][4];
        TH1F* h_Xgg_diffPhi[XGGCUT_NUM];
        TH1F* h_Xgg_timing[XGGCUT_NUM];
        TH1F* h_Xgg_sumE[XGGCUT_NUM];
        TH1F* h_Xgg_minE[XGGCUT_NUM];
        TH1F* h_Xgg_medE[XGGCUT_NUM];
        TH1F* h_Xgg_maxE[XGGCUT_NUM];

        //Moller histograms
        TH2F* h_M_HC_XY[MCUT_NUM];
        TH2F* h_M_E_theta[MCUT_NUM];
        TH1F* h_M_invM[MCUT_NUM];
        TH1F* h_M_Sum_pt[MCUT_NUM];
        TH2F* h_M_Sum_pxVpy[MCUT_NUM];
        TH1F* h_M_vZ[MCUT_NUM];
        TH1F* h_M_diffPhi[MCUT_NUM];
        TH1F* h_M_timing[MCUT_NUM];
        TH1F* h_M_sumE[MCUT_NUM];

        //Private X functions expecting a e+e- decay
        void setup_X_histos();
        void searchXEvent_electron(Int_t j, Int_t k, Int_t m);
        void fillCutHistos_electron(Int_t c, Int_t j, Int_t k, Int_t m);

        //Private X functions expecting a gamma-gamma decay
        void setup_Xgg_histos();
        void searchXEvent_photon(Int_t j, Int_t k, Int_t m);
        void fillCutHistos_photon(Int_t c, Int_t j, Int_t k, Int_t m);

        //Private Moller functions
        void setup_Moller_histos();
        void searchMollerEvent(Int_t j, Int_t k);
        void fillMollerCutHistos(Int_t c, Int_t j, Int_t k);

        void findMaxAndMinTime(Int_t j, Int_t k, Int_t m);
        Utils::Point findVertZ(Int_t j);
        void fillIndInfo(Int_t i, Int_t j);

};

#endif