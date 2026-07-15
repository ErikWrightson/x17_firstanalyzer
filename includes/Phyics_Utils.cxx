/**
 * Contains Utility functions that can be used independently.
 * 
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 07.13.2026
 * @creation 07.09.2026
 */
#include "Physics_Utils.h"

namespace Physics_Utils{
    /**
     * Find the energy resolution from HyCal at the given energy in MeV
     *
     * @param E - current particle energy in MeV
     *
     * @return - the energy resolution in MeV
     */
    static Double_t EnergyRes(Double_t E){
        return (0.033/TMath::Sqrt(E/1000.0))*E;
    }

    /**
     * Makes a direction vector assuming that this particle came from (0,0,0).
     * 
     * @param x - the x position of the particle
     * @param y - the y position of the particle
     * @param z - the z position of the particle
     *
     * @return - a direction vector struct that holds the length one direction vector for this particle.
     */
    static DirVector makeDirVector(Double_t x, Double_t y, Double_t z){
        Double_t mag = TMath::Sqrt(x*x + y*y + z*z);

        DirVector p;
        p.d_x = x/mag;
        p.d_y = y/mag;
	    p.d_z = z/mag;

        return p;
    }

    /**
     * Makes a direction vector given two points.
     *
     * @param x1 - the first x point
     * @param y1 - the first y point
     * @param z1 - the first z point
     *
     * @param x2 - the second x point
     * @param y2 - the second y point
     * @param z2 - the second z point
     *
     * @return - a normalized direction vector 
     */
    static DirVector makeDirVector2Point(Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2){
    
        Double_t z = z2 - z1;

        Int_t c = 1;
        if(z < 0){
            c = -1;
            z = z*c;
        }

        x = (x2 - x1)*c;
        y = (y2 - y1)*c;

        Double_t mag = TMath::Sqrt(x*x + y*y + z*z);

        DirVector p;
        p.d_x = x/mag;
        p.d_y = y/mag;
	    p.d_z = z/mag;

        return p;
    }

    /**
     * Makes a 4D Energy Momentum Vector assuming that the particle is an electron.
     *
     * @param E - energy of this particle
     * @param p - the direction vector for this particle
     *
     * @return - a 4D Energy-Momentum Vector using the PxPyPzE4D class from ROOT
     */
    static PxPyPzE4D<Double_t> make4D_EMomVector_Electron(Double_t E, DirVector p){
        Double_t P_mag = TMath::Sqrt(E*E - m_e*m_e);

        PxPyPzE4D<Double_t> v(P_mag*p.d_x, P_mag*p.d_y, P_mag*p.d_z, E);

        return v;
    }

    /**
     * Makes a 4D Energy Momentum Vector assuming that the particle is a photon.
     *
     * @param E - energy of this particle
     * @param p - the direction vector for this particle
     *
     * @return - a 4D Energy-Momentum Vector using the PxPyPzE4D class from ROOT
     */
    static PxPyPzE4D<Double_t> make4D_EMomVector_Photon(Double_t E, DirVector p){
        Double_t P_mag = TMath::Sqrt(E*E);

        PxPyPzE4D<Double_t> v(P_mag*p.d_x, P_mag*p.d_y, P_mag*p.d_z, E);

        return v;
    }

    /** --------------------------
     * @author Tyler Hague <tjhague@jlab.org>
     * @author Erik Wrightson <wrightso@jlab.org> - modified for more accuracy and the first layer cut
     *
     * Determine the hycal layer of the hit
     * Always at least cut to require that it's greater than 0
     * 0 means it went through the center (should never have data for that)
     * or that it didn't hit the PbWO4
     * Only
     * @param x - X Location in Hycal
     * @param y - Y location in Hycal
     *
     * @return an integer correcpsonding to the layer of HyCal struck
     * -------------------------
     */
    static Int_t hycal_layer(Float_t x, Float_t y){
	
	    if(TMath::Abs(x)>20.77*17.0 || TMath::Abs(y)>20.75*17 || (TMath::Abs(x)<(20.77+20.77/2.0) && TMath::Abs(y)<(20.75+20.75/2.0))){
		    //Event is out of the crystals or through the center
		    return 0;
	    }
	
	    Int_t row = TMath::Floor(TMath::Abs(x)/20.77);
	    Int_t column = TMath::Floor(TMath::Abs(y)/20.75);
	
	    Int_t layer = TMath::Max(row, column);
	    return layer;
    }

    /**
     * Finds the theta of the given particle.
     * 
     * @param x - the x coordinate
     * @param y - the y coordinate
     * @param z - the z coordinate
     *
     * @return - the theta of this particle
     */
    static Double_t findTheta(Double_t x, Double_t y, Double_t z){
        return TMath::ATan2(TMath::Sqrt(x*x+y*y),z)*rad2Deg;
    }

    /**
     * Finds the phi given the x and y of the particle provided.
     *
     * @param x - the x-direction used to find phi
     * @param y - the y-direction used to find phi
     *
     * @return - the phi of the particle
     */
    static Double_t findPhi(Double_t x, Double_t y){
        Double_t phi = TMath::ATan2(y, x);
		//Shift from (-Pi < phi_X < Pi) branch to (0 <= x_phi < 2Pi) branch for ease of visualization.
		if(phi<0){
			phi+=2*PI;
		}
		return phi*rad2Deg;
    }

    /**
     * Finds the Energy resolution of multiple particles with their resolutions added in quadrature.
     *
     * @param E - a vector of all of the energies as measured by HyCal.
     *
     * @return - the energy resolution of all of the particles added in quadrature
     */
    static Double_t combined_EnergyRes(vector<Double_t> E){
        Double_t res = 0;
        
        for(UInt_t i = 0; i < E.size(); i++){
            Double_t temp = EnergyRes(E.at(i));
            res = res + temp*temp;
        }

        return TMath::Sqrt(res);
    }

    /**
     * Calculates the expected energy of an electron in Moller (e-e) Scattering at a given angle and beam energy.
     * 
     * @param theta - the theta position of the electron on the calorimeter.
     *
     * @return - the expected energy if this is a Moller event
     */
    static Double_t ee_ExpectedE(Double_t theta, Double_t EBeam){
    
        Double_t cosTheta = TMath::Cos(theta);
        Double_t cosTheta_2 = cosTheta*cosTheta; 

        Double_t num = m_e*(EBeam + m_e + (EBeam-m_e)*cosTheta_2);
        Double_t denom  = EBeam + m_e - (EBeam-m_e)*cosTheta_2;

        return num/denom;
    }
}