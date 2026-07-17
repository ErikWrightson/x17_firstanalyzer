/**
 * Contains functions for a class object that find the invariant mass before and after each cut for X17.
 * 
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 07.13.2026
 * @creation 07.13.2026
 */
#include "MassFinder.h"

const TString MassFinder::MCUT_NAME[MCUT_NUM] = {"2 or More", "Timing", "Fiducial", "Kinematic (E_{exp})", "Coplanarity", "Elasticity", "1 GEM Match", "2 GEM Match", "Vertex Z"};
const TString MassFinder::XCUT_NAME[XCUT_NUM] = {"None", "3 or More", "Timing", "Fiducial", "Cluster E", "E Sum", "Coplanarity", "X_E", "1 GEM Match", "2 GEM Match", "Vertex Z"};

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

    chain->SetBranchStatus("cl_dt_rf", 0);

    setup_X_histos();
    setup_Moller_histos();

    partEnergies = {0,0,0};

}

/**
 * Sets up the histograms for the X particle analysis.
 */
void MassFinder::setup_X_histos(){
    TString m[4] = {"All", "Min_E", "Med_E", "Max_E"};
    for(Int_t i = 0; i < XCUT_NUM; i++){

        TString xy_Name  = TString::Format("h_X_HC_XY_cut%d", i);
        TString xy_Title = TString::Format("HyCal XY Position Run %d - Cut: %s;x [mm]; y [mm]", runNum, XCUT_NAME[i].Data());
        h_X_HC_XY[i] = new TH2F(xy_Name, xy_Title, HC_XBINS, HC_XMIN, HC_XMAX, HC_YBINS, HC_YMIN, HC_YMAX);

        TString th_Name  = TString::Format("h_X_E_theta_cut%d", i);
        TString th_Title = TString::Format("Energy v. #theta Run %d - Cut: %s;#theta [#circ]; Energy [MeV]", runNum, XCUT_NAME[i].Data());
        h_X_E_theta[i] = new TH2F(th_Name, th_Title, TH_BINS, MIN_TH, MAX_TH, Ebeam+200.0, 0.0 ,Ebeam+200.0);

        TString invM_Name = TString::Format("h_X_invM_cut%d",i);
        TString invM_Title = TString::Format("Invariant Mass of Potential Particle 0.5 MeV per bin Run %d - Cut: %s;Mass [MeV/c^{2}]; Count", runNum, XCUT_NAME[i].Data());
        h_X_invM[i] = new TH1F(invM_Name, invM_Title, INVM_BINS, MIN_INVM, MAX_INVM);

        TString sum_pt_Name = TString::Format("h_X_Sum_pt_cut%d", i);
        TString sum_pt_Title = TString::Format("3 particle #Sigmap_{t} Run %d - Cut: %s;p_{t} [MeV/c];Count", runNum, XCUT_NAME[i].Data());
        h_X_Sum_pt[i] = new TH1F(sum_pt_Name, sum_pt_Title, PT_BINS, MIN_PT, MAX_PT);

        TString sum_pxVpy_Name = TString::Format("h_X_pxVpy_cut%d", i);
        TString sum_pxVpy_Title = TString::Format("3 particle #Sigmap_{y} v. #Sigmap_{x} Run %d - Cut: %s;p_{x} [MeV/c];p_{y} [MeV/c]", runNum, XCUT_NAME[i].Data());
        h_X_Sum_pxVpy[i] = new TH2F(sum_pxVpy_Name, sum_pxVpy_Title, 2*PT_BINS, -1.0*MAX_PT, MAX_PT, 2*PT_BINS, -1.0*MAX_PT, MAX_PT);

        TString diffPhi_Name = TString::Format("h_X_diffPhi_cut%d", i);
        TString diffPhi_Title = TString::Format("#Delta#phi of particle candidate and e' Run %d - Cut: %s;#Delta#phi [#circ];Count/1#circ", runNum, XCUT_NAME[i].Data());
        h_X_diffPhi[i] = new TH1F(diffPhi_Name, diffPhi_Title, PHI_BINS, MIN_PHI, MAX_PHI);

        TString timing_Name = TString::Format("h_X_timing_cut%d", i);
        TString timing_Title = TString::Format("#Deltat (t_{cl_{max}} - t_{cl_{min}}) Run %d - Cut: %s;#Deltat [ns]; Count/ns", runNum, XCUT_NAME[i].Data());
        h_X_timing[i] = new TH1F(timing_Name, timing_Title, TIME_BINS, MIN_TIME, MAX_TIME);

        TString sum_E_Name = TString::Format("h_X_sumE_cut%d", i);
        TString sum_E_Title = TString::Format("3 particle E_{Sum} Run %d - Cut %s;E_{sum} [MeV]; Count/MeV", runNum, XCUT_NAME[i].Data());
        h_X_sumE[i] = new TH1F(sum_E_Name, sum_E_Title, Ebeam+500.0, 0.0 ,Ebeam+500.0);

        TString min_E_Name = TString::Format("h_X_minE_cut%d", i);
        TString min_E_Title = TString::Format("3 particle Minimum Cluster Energy E_{Min} Run %d - Cut %s;E_{min} [MeV];Count/MeV", runNum, XCUT_NAME[i].Data());
        h_X_minE[i] = new TH1F(min_E_Name, min_E_Title, Ebeam/2.0, 0.0 ,Ebeam/2.0);

        TString med_E_Name = TString::Format("h_X_medE_cut%d", i);
        TString med_E_Title = TString::Format("3 particle Median Cluster Energy E_{Med} Run %d - Cut %s;E_{med} [MeV];Count/MeV", runNum, XCUT_NAME[i].Data());
        h_X_medE[i] = new TH1F(med_E_Name, med_E_Title, Ebeam, 0.0 ,Ebeam);

        TString max_E_Name = TString::Format("h_X_maxE_cut%d", i);
        TString max_E_Title = TString::Format("3 particle Maximum Cluster Energy E_{Max} Run %d - Cut %s;E_{max} [MeV];Count/MeV", runNum, XCUT_NAME[i].Data());
        h_X_maxE[i] = new TH1F(max_E_Name, max_E_Title, Ebeam, 0.0 ,Ebeam);

        for(Int_t j = 0; j < 4; j++){

            TString vZ_Name = TString::Format("h_X_vZ_%s_cut%d",m[j].Data(), i);
            TString vZ_Title = TString::Format("V_{z} %s Run %d - Cut: %s", m[j].Data(), runNum, XCUT_NAME[i].Data());
            h_X_vZ[i][j] = new TH1F(vZ_Name, vZ_Title, VZ_BINS, MIN_VZ, MAX_VZ);
        }
    }
    cout<<"Set up X particle histograms." << endl;
}

/**
 * Save the histograms from this object to a ROOT file.
 *
 * @param rootFile - the full path of where to save these histograms to a ROOT file.
 */
void MassFinder::save_histos(TString rootFile){
    TObjArray* arr = new TObjArray(0,0);

    for(Int_t i = 0; i < XCUT_NUM; i++){
        (*arr).Add(h_X_HC_XY[i]);
        (*arr).Add(h_X_E_theta[i]);
        (*arr).Add(h_X_invM[i]);
        (*arr).Add(h_X_Sum_pt[i]);
        (*arr).Add(h_X_Sum_pxVpy[i]);
        (*arr).Add(h_X_diffPhi[i]);
        (*arr).Add(h_X_timing[i]);
        (*arr).Add(h_X_sumE[i]);
        (*arr).Add(h_X_minE[i]);
        (*arr).Add(h_X_medE[i]);
        (*arr).Add(h_X_maxE[i]);

        for(Int_t j = 0; j < 4; j++){
            (*arr).Add(h_X_vZ[i][j]);
        }
    }

    for(Int_t j = 0; j < MCUT_NUM; j++){
        (*arr).Add(h_M_HC_XY[j]);
        (*arr).Add(h_M_E_theta[j]);
        (*arr).Add(h_M_invM[j]);
        (*arr).Add(h_M_Sum_pt[j]);
        (*arr).Add(h_M_Sum_pxVpy[j]);
        (*arr).Add(h_M_diffPhi[j]);
        (*arr).Add(h_M_timing[j]);
        (*arr).Add(h_M_sumE[j]);
        (*arr).Add(h_M_vZ[j]);
    }

    TFile file1(rootFile,"RECREATE");
    (*arr).Write();
	file1.Close();
}

/**
 * @brief Searches through each event assuming all particles are electrons and fills histograms as such.
 */
void MassFinder::search_events_electrons(){
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        if((i+1)%10000 == 0 || (entries-i+1)<10000){
            cout<<"\rEvents searched: " << i+1 << "/" << entries << flush;
            if(i+1 == entries){
                cout<<endl;
            }
        }

        if(nClust > 10000){continue;} // This deals with odd events where a ridiculous amount of clusters are found.

        for(Int_t j = 0; j < nClust; j++){
            //Fill the angles, direction vector, and 4D energy-momentum vector for this particle
            fillIndInfo(0,j);

            //Fills the histograms.
            h_X_HC_XY[0]->Fill(cl_x[j], cl_y[j]);
            h_X_E_theta[0]->Fill(theta[0], cl_E[j]);
            if(vZ[0]<90000){h_X_vZ[0][0]->Fill(vZ[0]);} //Only fill histograms that have at least two points to rebuild the vertex to.
            
            for(Int_t k = j + 1; k < nClust; k++){

                //Fill the angles, direction vector, and 4D energy-momentum vector for this particle
                fillIndInfo(1,k);

                for(Int_t m = k + 1; m < nClust; m++){

                    //Fill the angles, direction vector, and 4D energy-momentum vector for this particle
                    fillIndInfo(2,m);

                    //if(nClust != 3){continue;}
                    searchXEvent_electron(j, k, m);
                    
                }

                if(nClust == 2) searchMollerEvent(j, k);
            }
        }
    }
}

/**
 * @brief Search the current set of clusters for potential X candidates.
 *
 * @param j - the index of the first cluster
 * @param k - the index of the second cluster
 * @param m - the index of the third cluster
 * @note - a helper of the search_events_electrons function
 */
void MassFinder::searchXEvent_electron(Int_t j, Int_t k, Int_t m){
    findMaxAndMinTime(j,k,m);

    //Fill the summed energy and resolution for these 3 particles.
    sumRes = combined_EnergyRes(partEnergies);
    sum = particle_e[0] + particle_e[1] + particle_e[2];

    //Fill the information for the potential X candidates.
    X[0] = particle_e[0] + particle_e[1]; //Particle 0 and 1
    X[1] = particle_e[0] + particle_e[2]; //Particle 0 and 2
    X[2] = particle_e[1] + particle_e[2]; //Particle 1 and 2
    X_eli = {0, 1, 2};
    for(Int_t n = 0; n < 3; n++){
        X_th[n] = findTheta(X[n].Px(), X[n].Py(), X[n].Pz());
        X_phi[n] = findPhi(X[n].Px(), X[n].Py());
    }

    phi_dif[0] = TMath::Abs(X_phi[0]-phi[2]);
    phi_dif[1] = TMath::Abs(X_phi[1]-phi[1]);
    phi_dif[2] = TMath::Abs(X_phi[2]-phi[0]);

    //Cut 1 - Require at least 3 clusters to exist.
    cut = 1;
    fillCutHistos_electron(cut,j,k,m);
                    
    //{"None", "3 or More", "Timing", "Fiducial", "Cluster E", "E Sum", "Coplanarity", "X_E", "1 GEM Match", "2 GEM Match", "Vertex Z"}
    //Cut 2 - Require that the maximum time difference between clusters is 16ns.
    if(dT < 16.0){
        cut++;
        fillCutHistos_electron(cut, j,k,m);

        //Cut 3 - Require that all clusters within the 2nd half of the first crystal layer and the inner edge of the last. 
        if(layer[0] > 0 && layer[0] < 16 && layer[1] > 0 && layer[1] < 16 && layer[2] > 0 && layer[2] < 16){
            cut++;
            fillCutHistos_electron(cut, j, k, m);

            //Cut 4 - Require that all clusters have energy above 70 MeV, and below 0.8 times beam energy.
            if(cl_E[j] > 70.0 && cl_E[j] < Ebeam*0.8 && cl_E[k] > 70.0 && cl_E[k] < Ebeam*0.8 && cl_E[m] > 70.0 && cl_E[m] < Ebeam*0.8){
                cut++;
                fillCutHistos_electron(cut, j, k, m);

                //Cut 5 - Require that the sum of the three clusters be within 5 sigma of beam energy given their respective resolutions added in quadrature.
                if(TMath::Abs(Ebeam-sum.E()) < sigma_E*sumRes){
                    cut++;
                    fillCutHistos_electron(cut, j, k, m);

                    //Cut 6 - Require that the each invariant mass candidate and the e' are coplanar within 10 degrees.
                    for(Int_t d = 0; d < 3; d++){
                        if(TMath::Abs(180-phi_dif[d]) > 15){X_eli.at(d) = -1;}
                    }
                    //Only go deeper into the if statements if there are still elligible candidates.
                    if(!(X_eli.at(0) < 0 && X_eli.at(1) < 0 && X_eli.at(0) < 0)){
                        cut++;
                        fillCutHistos_electron(cut, j, k, m);

                        //Cut 7 - Require that the each invariant mass candidate more than half the beam energy.
                        for(Int_t b = 0; b < 3; b++){
                            if(X_eli.at(b) >= 0 && X[b].E() < Ebeam*0.5){X_eli.at(b) = -1;}
                        }
                        //Only go deeper into the if statements if there are still elligible candidates.
                        if(!(X_eli.at(0) < 0 && X_eli.at(1) < 0 && X_eli.at(0) < 0)){
                            cut++;
                            fillCutHistos_electron(cut, j, k, m);

                            //Cut 8 - Require that every particle had at least one GEM match.
                            if(oneGEM[0] && oneGEM[1] && oneGEM[2]){
                                cut++;
                                fillCutHistos_electron(cut, j, k, m);

                                //Cut 9 - Require that every particle had at a match on both GEM layers.
                                if(twoGEM[0] && twoGEM[1] && twoGEM[2]){
                                        cut++;
                                        fillCutHistos_electron(cut, j, k, m);

                                        //Cut 10 - Require that every particle have a track within 2 meters of the target.
                                        if(TMath::Abs(vZ[0]) < 2000 && TMath::Abs(vZ[1]) < 2000 && TMath::Abs(vZ[2]) < 2000){
                                            cut++;
                                            fillCutHistos_electron(cut, j, k, m);
                                        }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Fills the histograms for the given cut and the indices of each particle.
 *
 * @param j - the index of the first particle to be included in the plots
 * @param k - the index of the second particle to be included in the plots
 * @param m - the index of the third particle to be included in the plots
 */
void MassFinder::fillCutHistos_electron(Int_t c, Int_t j, Int_t k, Int_t m){

    //Fill in the histograms that plot the individual cluster information.
    if(j == 0){ //Make sure to not accidentally double count particles.
        h_X_HC_XY[c]->Fill(cl_x[j], cl_y[j]);
        h_X_E_theta[c]->Fill(theta[0], cl_E[j]);

        if(vZ[0]<90000){
            h_X_vZ[c][0]->Fill(vZ[0]);
        }
    }
    if(k == 1){ //Make sure to not accidentally double count particles.
        h_X_HC_XY[c]->Fill(cl_x[k], cl_y[k]);
        h_X_E_theta[c]->Fill(theta[1], cl_E[k]);

        if(vZ[1]<90000){
            h_X_vZ[c][0]->Fill(vZ[1]);
        }
    }

    h_X_HC_XY[c]->Fill(cl_x[m], cl_y[m]);
    h_X_E_theta[c]->Fill(theta[2], cl_E[m]);
    h_X_sumE[c]->Fill(sum.E());

    //Fill the vertex position for each particle and then separate by minimum, median, and maximum energy particle.
    //TString m[4] = {"All", "Min_{E}", "Med_{E}", "Max_{E}"};
    if(vZ[2]<90000){
        h_X_vZ[c][0]->Fill(vZ[2]);
    }
    Int_t maxE_ind = distance(partEnergies.begin(), max_element(partEnergies.begin(), partEnergies.end()));
    if(vZ[maxE_ind]<90000){h_X_vZ[c][3]->Fill(vZ[maxE_ind]);}

    Int_t minE_ind = distance(partEnergies.begin(), min_element(partEnergies.begin(), partEnergies.end()));
    if(vZ[minE_ind]<90000){h_X_vZ[c][1]->Fill(vZ[minE_ind]);}
    
    Int_t medE_ind = 3 - (maxE_ind + minE_ind);
    if(vZ[medE_ind]<90000){h_X_vZ[c][2]->Fill(vZ[medE_ind]);}

    h_X_maxE[c]->Fill(partEnergies.at(maxE_ind));
    h_X_minE[c]->Fill(partEnergies.at(minE_ind));
    h_X_medE[c]->Fill(partEnergies.at(medE_ind));


    for(Int_t i = 0; i < 3; i++){
        if(X_eli.at(i) < 0) continue;

        h_X_invM[c]->Fill(X[X_eli.at(i)].M());
        h_X_diffPhi[c]->Fill(phi_dif[X_eli.at(i)]);//TMath::Abs(X_phi[X_ind.at(i)]-phi[ep_ind]));
    }

    h_X_Sum_pt[c]->Fill(sum.Pt());
    h_X_Sum_pxVpy[c]->Fill(sum.Px(), sum.Py());
    h_X_timing[c]->Fill(dT);

}

/**
 * Finds the maximum and mimum cluster time of each of the clusters handed in and assigns it to the proper varibales for this object
 * for these clusters, and then gets the difference.
 *
 * @param j - the index of the first cluster
 * @param k - the index of the second cluster
 * @param m - the index of the third cluster
 */
void MassFinder::findMaxAndMinTime(Int_t j, Int_t k, Int_t m){
    vector<Double_t> times = {cl_time[j], cl_time[k], cl_time[m]};

    maxT = *(max_element(times.begin(), times.end()));
    minT = *(min_element(times.begin(), times.end()));

    dT = maxT - minT;
}

/**
 * Finds the Z vertex as defined by the point of closest approach to the z-axis.
 *
 * @param j - the index of the hit to fit.
 */
Utils::Point MassFinder::findVertZ(Int_t j){
    vector<Float_t> x_vec;
    vector<Float_t> y_vec;
    vector<Float_t> z_vec;
    /*x_vec.push_back(cl_x[j]);
    y_vec.push_back(cl_y[j]);
    z_vec.push_back(cl_z[j]);*/ //Currently the HyCal position still needs correction so they are not being used in the line of best fit yet.
    for(Int_t q = 0; q < MAX_GEMS; q++){
        if(mgz[j][q]>0){
            x_vec.push_back(mgx[j][q]);
            y_vec.push_back(mgy[j][q]);
            z_vec.push_back(mgz[j][q]);
        }
    }

    Utils::Point p;
    p.z = 100000;
    if(x_vec.size() > 1){
        Utils::LineOfBestFit line = Utils::FitLine(x_vec, y_vec, z_vec);
        p = Utils::ClosestApproachToZAxis(line);
    }

    return p;
}

/**
 * Fills the information for an individual particle.
 *
 * @param i - the particle index to fill
 * @param j - the cluster index to use for filling
 * @note - A helper for search_events_electrons method
 */
void MassFinder::fillIndInfo(Int_t i, Int_t j){
    theta[i] = findTheta(cl_x[j], cl_y[j], cl_z[j]);
    phi[i] = findPhi(cl_x[j], cl_y[j]);
    dir[i] = makeDirVector(cl_x[j], cl_y[j], cl_z[j]);
    particle_e[i] = make4D_EMomVector_Electron(cl_E[j], dir[i]);
    partEnergies.at(i) = cl_E[j];
    layer[i] = hycal_layer(cl_x[j], cl_y[j]);
    close[i] = findVertZ(j);
    vZ[i] = close[i].z;
    oneGEM[i] = (match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1)) || (match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3));
    twoGEM[i] = ((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) && ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)));
}

/**
 * @brief Performs all cuts for the clusters found at the provided indices.
 *
 * @param j - the index of the first particle of the Moller pair
 * @param k - the index of the second particle of the Moller pair
 */
void MassFinder::searchMollerEvent(Int_t j, Int_t k){
    MsumRes = combined_EnergyRes({cl_E[j], cl_E[k]});
    Mol = particle_e[0] + particle_e[1];
    dT_Mol = TMath::Abs(cl_time[j]-cl_time[k]);
    phi_dif_Mol = TMath::Abs(phi[0] - phi[1]);
    exp_E[0] = ee_ExpectedE(theta[0], Ebeam);
    exp_E[1] = ee_ExpectedE(theta[1], Ebeam);

    //Cut 0 - Require at least 2 clusters to exist.
    //{"2 or More", "Timing", "Fiducial", "Kinematic (E_{exp})", "Coplanarity", "Elasticity", "1 GEM Match", "2 GEM Match", "Vertex Z"}
    Int_t Mcut = 0;
    fillMollerCutHistos(Mcut, j, k);

    //Cut 1 - Require that the maximum time difference between clusters is 16ns.
    if(dT_Mol < 16.0){
        Mcut++;
        fillMollerCutHistos(Mcut, j, k);

        //Cut 2 - Require that all clusters within the 2nd half of the first crystal layer and the inner edge of the last. 
        if(layer[0] > 0 && layer[0] < 16 && layer[1] > 0 && layer[1] < 16){
            Mcut++;
            fillMollerCutHistos(Mcut, j, k);

            //This cut has been removed: Cut 3 - Require that all clusters pass the expected Moller Energy cut for the angle they are at. 
            if(true){//TMath::Abs(cl_E[j]-exp_E[0]) < sigma_E*EnergyRes(exp_E[0]) && TMath::Abs(cl_E[k]-exp_E[1]) < sigma_E*EnergyRes(exp_E[1])){
                Mcut++;
                fillMollerCutHistos(Mcut, j, k);

                //Cut 4 - Require that both clusters are coplanar with on another within 10 degrees.
                if(TMath::Abs(180-phi_dif_Mol) < 10){
                    Mcut++;
                    fillMollerCutHistos(Mcut, j, k);

                    //Cut 5 - Require that both clusters add up to being an elastic event given the Energy resolutions added in quadrature.
                    if(TMath::Abs(cl_E[j] + cl_E[k] - Ebeam - m_e) < sigma_E*MsumRes){
                        Mcut++;
                        fillMollerCutHistos(Mcut, j, k);

                        //Cut 6 - Require that both clusters have a match on at least one of the GEM planes.
                        if(oneGEM[0] && oneGEM[1]){
                            Mcut++;
                            fillMollerCutHistos(Mcut, j, k);

                            //Cut 7 - Require that both clusters have a match on each of the GEM planes.
                            if(twoGEM[0] && twoGEM[1]){
                                Mcut++;
                                fillMollerCutHistos(Mcut, j, k);

                                if(TMath::Abs(vZ[0]) < 2000 && TMath::Abs(vZ[1]) < 2000){
                                    Mcut++;
                                    fillMollerCutHistos(Mcut, j, k);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Fills the Moller Event histograms
 *
 * @param j - the index of the first cluster
 * @param k - the index of the second cluster
 */
void MassFinder::fillMollerCutHistos(Int_t c, Int_t j, Int_t k){
    //Fill in the histograms that plot the individual cluster information.
    if(j == 0){ //Make sure to not accidentally double count particles.
        h_M_HC_XY[c]->Fill(cl_x[j], cl_y[j]);
        h_M_E_theta[c]->Fill(theta[0], cl_E[j]);

        if(vZ[0]<90000){
            h_M_vZ[c]->Fill(vZ[0]);
        }
    }

    h_M_HC_XY[c]->Fill(cl_x[k], cl_y[k]);
    h_M_E_theta[c]->Fill(theta[1], cl_E[k]);
    //Fill the vertex position for each particle if one exists.
    if(vZ[1]<90000){
        h_M_vZ[c]->Fill(vZ[1]);
    }

    h_M_timing[c]->Fill(dT_Mol);
    h_M_invM[c]->Fill(Mol.M());
    h_M_Sum_pt[c]->Fill(Mol.Pt());
    h_M_sumE[c]->Fill(Mol.E());
    h_M_Sum_pxVpy[c]->Fill(Mol.Px(), Mol.Py());
    h_M_diffPhi[c]->Fill(phi_dif_Mol);
}

/**
 * Instantiates each relevant histogram for the Moller event selections.
 */
void MassFinder::setup_Moller_histos(){

    for(Int_t i = 0; i < MCUT_NUM; i++){

        TString xy_Name  = TString::Format("h_M_HC_XY_cut%d", i);
        TString xy_Title = TString::Format("Moller HyCal XY Position Run %d - Cut: %s;x [mm]; y [mm]", runNum, MCUT_NAME[i].Data());
        h_M_HC_XY[i] = new TH2F(xy_Name, xy_Title, HC_XBINS, HC_XMIN, HC_XMAX, HC_YBINS, HC_YMIN, HC_YMAX);

        TString th_Name  = TString::Format("h_M_E_theta_cut%d", i);
        TString th_Title = TString::Format("Moller Energy v. #theta Run %d - Cut: %s;#theta [#circ]; Energy [MeV]", runNum, MCUT_NAME[i].Data());
        h_M_E_theta[i] = new TH2F(th_Name, th_Title, TH_BINS, MIN_TH, MAX_TH, Ebeam+200.0, 0.0 ,Ebeam+200.0);

        TString invM_Name = TString::Format("h_M_invM_cut%d",i);
        TString invM_Title = TString::Format("Moller Invariant Mass of Potential Particle 1/3 MeV per bin Run %d - Cut: %s;Mass [MeV/c^{2}]; Count", runNum, MCUT_NAME[i].Data());
        h_M_invM[i] = new TH1F(invM_Name, invM_Title, INVM_BINS, MIN_INVM, MAX_INVM);

        TString sum_pt_Name = TString::Format("h_M_Sum_pt_cut%d", i);
        TString sum_pt_Title = TString::Format("Moller 2 particle #Sigmap_{t} Run %d - Cut: %s;p_{t} [MeV/c];Count", runNum, MCUT_NAME[i].Data());
        h_M_Sum_pt[i] = new TH1F(sum_pt_Name, sum_pt_Title, PT_BINS, MIN_PT, MAX_PT);

        TString sum_pxVpy_Name = TString::Format("h_M_pxVpy_cut%d", i);
        TString sum_pxVpy_Title = TString::Format("Moller 2 particle #Sigmap_{y} v. #Sigmap_{x} Run %d - Cut: %s;p_{x};p_{y}", runNum, MCUT_NAME[i].Data());
        h_M_Sum_pxVpy[i] = new TH2F(sum_pxVpy_Name, sum_pxVpy_Title, 2*PT_BINS, -1.0*MAX_PT, MAX_PT, 2*PT_BINS, -1.0*MAX_PT, MAX_PT);

        TString diffPhi_Name = TString::Format("h_M_diffPhi_cut%d", i);
        TString diffPhi_Title = TString::Format("Moller #Delta#phi of particle candidate and e' Run %d - Cut: %s;#Delta#phi [#circ];Count/1#circ", runNum, MCUT_NAME[i].Data());
        h_M_diffPhi[i] = new TH1F(diffPhi_Name, diffPhi_Title, PHI_BINS, MIN_PHI, MAX_PHI);

        TString timing_Name = TString::Format("h_M_timing_cut%d", i);
        TString timing_Title = TString::Format("Moller #Deltat (t_{cl_{max}} - t_{cl_{min}}) Run %d - Cut: %s;#Deltat [ns]; Count/ns", runNum, MCUT_NAME[i].Data());
        h_M_timing[i] = new TH1F(timing_Name, timing_Title, TIME_BINS, MIN_TIME, MAX_TIME);
       
        TString vZ_Name = TString::Format("h_M_vZ_cut%d", i);
        TString vZ_Title = TString::Format("Moller V_{z} Run %d - Cut: %s; z [mm]; Count/1cm", runNum, MCUT_NAME[i].Data());
        h_M_vZ[i] = new TH1F(vZ_Name, vZ_Title, VZ_BINS, MIN_VZ, MAX_VZ);

        TString sum_E_Name = TString::Format("h_M_sumE_cut%d", i);
        TString sum_E_Title = TString::Format("Moller 2 particle E_{Sum} Run %d - Cut %s;E_{sum} [MeV]; Count/MeV", runNum, MCUT_NAME[i].Data());
        h_M_sumE[i] = new TH1F(sum_E_Name, sum_E_Title, Ebeam+200.0, 0.0 ,Ebeam+200.0);
    }
    cout<<"Set up Moller histograms." << endl;
}

void MassFinder::delete_histos(){
    for(Int_t i = 0; i < XCUT_NUM; i++){
        delete h_X_HC_XY[i];
        delete h_X_E_theta[i];
        delete h_X_invM[i];
        delete h_X_Sum_pt[i];
        delete h_X_Sum_pxVpy[i];
        delete h_X_diffPhi[i];
        delete h_X_timing[i];
        delete h_X_sumE[i];
        delete h_X_minE[i];
        delete h_X_medE[i];
        delete h_X_maxE[i];

        for(Int_t j = 0; j < 4; j++){
            delete h_X_vZ[i][j];
        }
    }

    for(Int_t j = 0; j < MCUT_NUM; j++){
        delete h_M_HC_XY[j];
        delete h_M_E_theta[j];
        delete h_M_invM[j];
        delete h_M_Sum_pt[j];
        delete h_M_Sum_pxVpy[j];
        delete h_M_diffPhi[j];
        delete h_M_timing[j];
        delete h_M_vZ[j];
        delete h_M_sumE[j];
    }
}