/**
 * Trigger Validation Processor
 * @author Erik Wrightson
 * @version 06.08.2026
 * @creation 04.05.2026
 */

#include <iostream>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>


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

#include "includes/MassFinder.h"
#include "includes/Utils.h"
#include "includes/Physics_Utils.h"

using namespace std;

/**
 * The main function that launches the trigger analysis.
 *
 * @param argc - the number of input arguments
 * @param argv - an array of the different arguments as an array of char* (strings).
 */
int main (int argc, char **argv){

    string fileName;
    string fileListFileName;

    TString outputDirectory = "outfiles/";
    TString rootOutputDirectory = "rootOutfiles/";
    TString fn = "default";

	if (argc<2) {
		cout<<"ERR: Incorrect Arguments: " <<endl;
        Utils::printUsage(argv[0]);
		
		return -1;
	}
    

    // ── Parse command-line ───────────────────────────────────────────────
    int opt;
    while ((opt = getopt(argc, argv, "f:L:o:")) != -1) {
        switch (opt) {
            case 'f': fileName = optarg; break;
            case 'L': fileListFileName = optarg; break;
            case 'o': fn = optarg; break;
            case 'h':
            default: Utils::printUsage(argv[0]); return (opt == 'h') ? 0 : 1;
        }
    }

    struct stat buffer;   
    bool existOne = (stat(fileName.c_str(), &buffer) == 0);
    bool existList = (stat(fileListFileName.c_str(), &buffer) == 0);

    if(!existOne && !existList){
        cerr<<"A single valid input file or a filelist txt file was not provided.\n";
        return -3;
    }

    vector<TString> fileNameVec;
    if(existList){
        fileNameVec = Utils::processFileList(fileListFileName);
    }
    if(existOne){
        fileNameVec.push_back((TString) fileName);
    }

    TChain* fChain = Utils::makeChain(fileNameVec, "recon");
    
    MassFinder mass = new MassFinder(fChain);

    mass.search_events_electrons();
    mass.save_histos(outputDirectory + fn + ".root");
    mass.delete_histos();

    return 0;
}