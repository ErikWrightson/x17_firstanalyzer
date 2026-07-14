/**
 * Contains functions for a class object that find the invariant mass before and after each cut for X17.
 * 
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 07.13.2026
 * @creation 07.13.2026
 */
#include "MassFinder.h"

MassFinder::MassFinder(TChain* c, Int_t rN){
    runNum = rN;
    Ebeam = 2239.51;

    chain = (TChain*) c;

    chain->SetMakeClass(1);

    //Get the amount of entries from each file to limit looping through them.
	entries = chain->GetEntries();

    //General Event Data
    chain->SetBranchAddress("event_num",       &eventNum);        //Event Number for this event so it can be matched between raw and reconstructed files.
    chain->SetBranchAddress("trigger_type",    &trigType);        //Trigger type for this event during the reconstruction window.
    chain->SetBranchAddress("trigger_bits",    &trigger_bits);    //Trigger bit word of triggers that happened in the reconstruction window.
    chain->SetBranchAddress("timestamp",       &time);            //Overall timestamp for this event.
    chain->SetBranchAddress("ssp_raw",         &sspRawPtr);       //Raw SSP trigger decision information to be further decoded.
    chain->SetBranchAddress("vtp_roc_tags",    &vtp_roc_tagsPtr); //ROC tag vector for the VTP to be further decoded.
    chain->SetBranchAddress("vtp_nwords",      &vtp_nwordsPtr);   //Number of words in each entry of the VTP words vector to be further decoded.
    chain->SetBranchAddress("vtp_words",       &vtp_wordsPtr);    //Vector of words from the VTP to be further decoded.

    //HyCal Information
    chain->SetBranchAddress("n_clusters", &nClust);
    chain->SetBranchAddress("cl_x",       cl_x);       //Cluster x position
    chain->SetBranchAddress("cl_y",       cl_y);       //Cluster y position
    chain->SetBranchAddress("cl_z",       cl_z);       //Cluster z position
    chain->SetBranchAddress("cl_energy",  cl_E);       //Cluster energy
    chain->SetBranchAddress("cl_nblocks", cl_nblocks); //Number of blocks in the cluster
    chain->SetBranchAddress("cl_center",  cl_center);  //center module id for this cluster
    chain->SetBranchAddress("cl_time",    cl_time);    //time of the seed module for this cluster
    chain->SetBranchAddress("cl_flag",    cl_flag);    //Cluster flags

    chain->SetBranchAddress("vtp_cl_n",      &vtp_cl_n);      //Number of clusters found by the vtp for each event.
    chain->SetBranchAddress("vtp_cl_time",   vtp_cl_time);    //Timing for each cluster found by the vtp.
    chain->SetBranchAddress("vtp_cl_energy", vtp_cl_E);       //Cluster energy for each cluster found by the vtp.
    chain->SetBranchAddress("vtp_cl_center", vtp_cl_center);  //Module ID of the seed for each cluster found by the vtp.
    chain->SetBranchAddress("vtp_cl_blocks", vtp_cl_nblocks); //Number of blocks in each cluster found by the vtp.

    chain->SetBranchAddress("matchFlag", match_flag); //Matching Flag bit 0 for GEM0, bit 1 for GEM1 etc.
    chain->SetBranchAddress("mHit_gx", matchGEMx);    //The x-coordinate of the best match found on each GEM plane.
    chain->SetBranchAddress("mHit_gy", matchGEMy);    //The y-coordinate of the best match found on each GEM plane.
    chain->SetBranchAddress("mHit_gz", matchGEMz);    //The z-coordinate of the best match found on each GEM plane.

    chain->SetBranchAddress("matchGEMx",  mgx);       //All x-coordinate matches found for on the GEMs.
    chain->SetBranchAddress("matchGEMy",  mgy);       //All y-coordinate matches found for on the GEMs.
    chain->SetBranchAddress("matchGEMz",  mgz);       //All z-coordinate matches found for on the GEMs.

    setup_X_histos();

}

/**
 * Sets up the histograms for the X particle analysis.
 */
void MassFinder::setup_X_histos(){
    TString m[4] = {"All", "Min_{E}", "Med_{E}", "Max_{E}"};
    for(Int_t i = 0; i < CUT_NUM; i++){

        TString xy_Name  = TString::Format("h_X_HC_XY_cut%d", i);
        TString xy_Title = TString::Format("HyCal XY Position Run %b - Cut: %s;x [mm]; y [mm]", runNum, CUT_NAME[i]);
        h_X_HC_XY[i] = new TH2F(xy_Name, xy_Title, HC_XBINS, HC_XMIN, HC_XMAX, HC_YBINS, HC_YMIN, HC_YMAX);

        TString th_Name  = TString::Format("h_X_E_theta_cut%d", i);
        TString th_Title = TString::Format("Energy v. #theta Run %b - Cut: %s;#theta [#circ]; Energy [MeV]", runNum, CUT_NAME[i]);
        h_X_E_theta[i] = new TH2F(th_Name, th_Title, TH_BINS, MIN_TH, MAX_TH, Ebeam+200.0, 0.0 ,Ebeam+200.0);

        TString invM_Name = TString::Format("h_X_invM_cut%d",i);
        TString invM_Title = TString::Format("Invariant Mass of Potential Particle 1/3 MeV per bin Run %d - Cut: %s;Mass [MeV/c^{2}]; Count", runNum, CUT_NAME[i]);
        h_X_invM[i] = new TH1F(invM_Name, invM_Title, INVM_BINS, MIN_INVM, MAX_INVM);

        TString pt_Name = TString::Format("h_X_pt_cut%d", i);
        TString pt_Title = TString::Format("p_{t} Run %d - Cut: %s;p_{t} [MeV/c];Count", runNum, CUT_NAME[i]);
        h_X_pt[i] = TH1F(pt_Name, pt_Title, PT_BINS, MIN_PT, MAX_PT);

        TString pxVpy_Name = TString::Format("h_X_pxVpy_cut%d", i);
        TString pxVpy_Title = TString::Format("p_{y} v. p_{x} Run %d - Cut: %s", runNum, CUT_NAME[i]);
        h_X_pxVpy[i] = TH2F(pxVpy_Name, pxVpy_Title, PT_BINS, MIN_PT, MAX_PT, PT_BINS, MIN_PT, MAX_PT);

        TString diffPhi_Name = TString::Format("h_X_diffPhi_cut%d", i);
        TString diffPhi_Title = TString::Format("#Delta#phi of particle candidate and e' Run %d - Cut: %s", runNum, CUT_NAME[i]);
        h_X_diffPhi[i] = new TH1F(diffPhi_Name, diffPhi_Title, PHI_BINS, MIN_PHI, MAX_PHI);

        TString timing_Name = TString::Format("h_X_timing_cut%d", i);
        TString timing_Title = TString::Format("#Deltat (t_{cl_{max}} - t_{cl_{min}}) Run %d - Cut: %s", runNum, CUT_NAME[i]);
        h_X_timing[i] = TH1F(timing_Name, timing_Title, TIME_BINS, MIN_TIME, MAX_TIME);

        for(Int_t j = 0; j < 4; j++){

            TString vZ_Name = TString::Format("h_X_vZ_%s_cut%d",m[j], i);
            TString vZ_Title = TString::Format("V_{z} %s Run %d - Cut: %s", m[j], runNum, CUT_NAME[i]);
            h_X_vZ[i][j] = new TH1F(vZ_Name, vZ_Title, VZ_BINS, MIN_VZ, MAX_VZ);
        }
    }
}

/**
 * Save the histograms from this object to a ROOT file.
 *
 * @param rootFile - the full path of where to save these histograms to a ROOT file.
 */
void save_histos(TString rootFile){

    TObjArray* arr = new TObjArray(0,0);

    for(Int_t i = 0; i < CUT_NUM; i++){
        (*arr).Add(h_X_HC_XY[i]);
        (*arr).Add(h_X_E_theta[i]);
        (*arr).Add(h_X_invM[i]);
        (*arr).Add(h_X_pt[i]);
        (*arr).Add(h_X_pxVpy[i]);
        (*arr).Add(h_X_diffPhi[i]);
        (*arr).Add(h_X_timing[i]);

        for(Int_t j = 0; j < 4; j++){
            (*arr).Add(h_X_vZ[i][j]);
        }
    }

    TFile file1(rootFile,"RECREATE");
    (*arr).Write();
	file1.Close();
}

void search_events_electrons(){
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        for(Int_t j = 0; j < nClust; j++){

            //Fill the angles, direction vector, and 4D energy-momentum vector for this particle
            theta[0] = findTheta(cl_x[j], cl_y[j], cl_z[j]);
            phi[0] = findPhi(cl_x[j], cl_y[j]);
            dir[0] = makeDirVector(cl_x[j], cl_y[j], cl_z[j]);
            particle_e[0] = make4D_EMomVector_Electron(cl_E[j], dir[0]);
            partEnergies.push_back(cl_E[j]);

            //Fills the histograms.
            h_X_HC_XY[0]->Fill(cl_x[j], cl_y[j]);
            h_X_E_theta[0]->Fill(theta[0], cl_E[j]);
            h_X_pt[0]->Fill(particle_e.Pt());
            h_X_pxVpy[0]->Fill(particle_e.Px());
            
            for(Int_t k = j + 1; k < nClust; k++){

                //Fill the angles, direction vector, and 4D energy-momentum vector for this particle
                theta[1] = findTheta(cl_x[k], cl_y[k], cl_z[k]);
                phi[1] = findPhi(cl_x[k], cl_y[k]);
                dir[1] = makeDirVector(cl_x[k], cl_y[k], cl_z[k]);
                particle_e[1] = make4D_EMomVector_Electron(cl_E[k], dir[1]);
                partEnergies.push_back(cl_E[k]);

                for(Int_t m = k + 1; m < nClust; m++){

                    //Fill the angles, direction vector, and 4D energy-momentum vector for this particle
                    theta[2] = findTheta(cl_x[m], cl_y[m], cl_z[m]);
                    phi[2] = findPhi(cl_x[m], cl_y[m]);
                    dir[2] = makeDirVector(cl_x[m], cl_y[m], cl_z[m]);
                    particle_e[2] = make4D_EMomVector_Electron(cl_E[m], dir[2]);
                    partEnergies.push_back(cl_E[m]);

                    
                    Double_t sumRes = combined_EnergyRes(partEnergies);

                    X[0] = particle_e[0] + particle_e[1]; //Particle 0 and 1
                    X[1] = particle_e[0] + particle_e[2]; //Particle 0 and 2
                    X[2] = particle_e[1] + particle_e[2]; //Particle 0 and 3
                    for(Int_t n = 0; n < 3; n++){
                        X_th[n] = findTheta(X[n].Px(), X[n].Py(), X[n].Pz());
                        X_phi[n] = findPhi(X[n].Px(), X[n].Py());
                    }
                    
                }
            }
        }
    }
}

void fillCutHistos(Int_t c){
    
}