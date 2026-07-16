/**
 * Contains Utility functions that can be used independently.
 * 
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 06.15.2026
 * @creation 06.15.2026
 */
#include "Utils.h"

/**
 * Prints out the proper usage directions for this program and what each flag means.
 *
 * @param prog - the progam name that is being currently run.
 */
void Utils::printUsage(const char *prog){
    cerr << "Usage: " << prog << " [options]\n"
              << "\t-f <fileName> input root file in the PRad-II/X17 reconstructed format\n"
              << "\t-L <listFileName> input text file containing a list of root files in the PRadII-X17 format\n"
              << "\t-b <runNumber> input the run Number"
              << "\t-h Show this help\n"
              << "\tNOTE: Either option -f or -L are REQUIRED for running properly.\n";
}

/**
 * Processes the file with the input file name and reads out the root file names contained within and adds them to a vector.
 *
 * @param fileListFileName - the name of the file with the list of ROOT file names.
 */
vector<TString> Utils::processFileList(string fileListFileName){
        ifstream file(fileListFileName);
        string line;
        
        vector<TString> list;
        while(getline(file,line)){
            TString l(line);
            list.push_back(l);
        }

        return list;
}

/**
 * Makes a TChain of the entries in a TString vector assuming they are valid paths to ROOT files.
 *
 * @param names - the vector of ROOT file names to be linked in the chain.
 */
TChain* Utils::makeChain(vector<TString> names, TString treeName){
        TChain* chain = new TChain(treeName);//"recon");//"events");
        for(unsigned int i = 0; i < names.size(); i++){
            chain->Add(names.at(i));
        }

        return chain;
}

void Utils::makeNeighborMap(string dbName, map<string, vector<string>> &m, map<Int_t, vector<Int_t>> &m2){
    ifstream file(dbName);
    string line;

    while(getline(file,line)){
        stringstream ss(line);
        string word;
        
        Int_t wordNum = 0;
        string key;
        vector<string> val;

        Int_t key2;
        vector<Int_t> val2;

        while (ss >> word) {
            if(wordNum == 0){
                key = word;
                string temp = word;
                if(temp.front() == 'G'){
                    temp.erase(0,1);
                    key2 = stoi(temp);
                }
                else{
                    temp.erase(0,1);
                    key2 = 1000+stoi(temp);
                }
            }
            if(wordNum != 0){
                if(wordNum == 1){
                    word.erase(0,2);
                    word.erase(word.length()-2, word.length());
                }
                else if(word.back() == ']'){
                    word.erase(0,1);
                    word.erase(word.length()-2, word.length());
                }
                else{
                    word.erase(0,1);
                    word.erase(word.length()-2, word.length());
                }
                val.push_back(word);

                string temp = word;
                if(temp.front() == 'G'){
                    temp.erase(0,1);
                    val2.push_back(stoi(temp));
                }
                else{
                    temp.erase(0,1);
                    val2.push_back(1000+stoi(temp));
                }
            }
            wordNum++;
        }
        m[key] = val;
        m2[key2] = val2;
    }

    cout<<"Made the quick reference neighbor map.\n";
}

void Utils::makeGainMap(string dbName, map<string, Float_t> &m, map<Int_t, Float_t> &m2){
    ifstream file(dbName);
    string line;

    while(getline(file,line)){
        stringstream ss(line);
        string word;
        
        Int_t wordNum = 0;
        string key;
        Float_t val;

        Int_t key2;

        while (ss >> word) {
            if(wordNum == 0){
                key = word;
                string temp = word;
                if(temp.front() == 'G'){
                    temp.erase(0,1);
                    key2 = stoi(temp);
                }
                else{
                    temp.erase(0,1);
                    key2 = 1000+stoi(temp);
                }
            }
            else{
                val = stof(word);
            }
            wordNum++;
        }
        m[key] = val;
        m2[key2] = val;
    }

    cout<<"Made the quick reference gain map.\n";
}

 Utils::LineOfBestFit Utils::FitLine(vector<Float_t> x, vector<Float_t> y, vector<Float_t> z){
    Utils::LineOfBestFit line;

    int n = x.size();

    // TPrincipal does PCA
    TPrincipal pca(3, "");

    Double_t point[3];

    for (int i = 0; i < n; i++) {
        point[0] = x[i];
        point[1] = y[i];
        point[2] = z[i];

        pca.AddRow(point);
    }

    pca.MakePrincipals();

    // Mean point: point on the fitted line
    const TVectorD *mean = pca.GetMeanValues();

    line.x0 = (*mean)[0];
    line.y0 = (*mean)[1];
    line.z0 = (*mean)[2];

    // Eigenvectors: first principal component gives best-fit direction
    const TMatrixD *eig = pca.GetEigenVectors();

    line.vx = (*eig)(0,0);
    line.vy = (*eig)(1,0);
    line.vz = (*eig)(2,0);

    return line;
}

/**
 * Returns the point of the closest approach from the z axis for given the line of best fit handed in.
 *
 * @param line - the line of best fit to find its approach to the z-axis.
 *
 * @return - the point of closest approach to the z axis.
 */
Utils::Point Utils::ClosestApproachToZAxis(Utils::LineOfBestFit line){
    Float_t t = -1.0 * ((line.x0*line.vx)+(line.y0*line.vy))/((line.vx*line.vx) + (line.vy*line.vy));

    Utils::Point p;
    p.x = line.x0 + t*line.vx;
    p.y = line.y0 + t*line.vy;
    p.z = line.z0 + t*line.vz;

    return p;
}

/**
 * Extracts the first continuous set of numbers in the provided TString.
 *
 * @param str - a call by reference to the TString object to search
 *
 * @return - the first continuous set of numbers in the string converted to an Int_t
 */
Int_t Utils::extractFirstInt(const TString& str){
    TString digits = "";

    bool foundDigit = false;

    for (int i = 0; i < str.Length(); i++) {
        if (str[i] >= '0' && str[i] <= '9') {
            digits += str[i];
            foundDigit = true;
        } 
        else if (foundDigit && digits.Length() > 4) {
            break;  // stop after the first continuous block with more than 4 digits
        }
        else if(foundDigit && digits.Length() <= 4){
            foundDigit = false;
            digits = "";
        }
    }

    return foundDigit ? digits.Atoi() : -1;
}