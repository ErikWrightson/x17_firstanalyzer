/**
 * Header file for Physics_Utils class which provides useful physics utilities functions to this project.
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 07.13.2026
 * @creation 07.09.2026
 */

#ifndef Physics_Utils_H
#define Physics_Utils_H

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

//Generally useful includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <iomanip>

using namespace std;

namespace Physics_Utils{

    public:

        static constexpr Double_t m_e = 0.511;
        static constexpr Double_t rad2Deg = 180/TMath::Pi(); //Conversion from radians to degrees

        struct DirVector{
            Double_t d_x;
            Double_t d_y;
            Double_t d_z;
        };

        static Double_t EnergyRes(Double_t E);

        static DirVector makeDirVector(Double_t x, Double_t y, Double_t z);
        static DirVector makeDirVector2Point(Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2);
        static PxPyPzE4D make4D_EMomVector_Electron(Double_t E, DirVector p);
        static PxPyPzE4D make4D_EMomVector_Photon(Double_t E, DirVector p);
        static Int_t hycal_layer(Float_t x, Float_t y);
        static Double_t findTheta(Double_t x, Double_t y, Double_t z);
        static Double_t findPhi(Double_t x, Double_t y);
        static Double_t combined_EnergyRes(Double_t E1, Double_t E2, Double_t E3);

};

#endif