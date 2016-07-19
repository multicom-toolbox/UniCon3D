#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cassert>
#include <functional>
#include <string>
#include <vector>
#include <cfloat>
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mocapy.h"
#include <time.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <unistd.h>
#include "mocapy.h"
#include "vonmises2d/vonmises2dess.h"
#include "vonmises2d/vonmises2ddensities.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "MatVec.h"

#define X 0
#define Y 1
#define Z 2

using namespace mocapy;
using namespace std;
using std::vector;

const double PI = std::atan(1.0)*4;
const double RAD2DEG = 180/PI;
const double DEG2RAD = PI/180;
const double CA2CA = 3.8;
const double START_TEMP = 1000.0;
const double FINAL_TEMP = 298.0;
const double BOLTZMANN_CONSTANT = 0.0019872041;
const double RR_CONTACT_THRESHOLD = 8.0;
const int RR_SEQ_SEP = 1;
const int MIN_FOLDON_LEN = 20;
const int NUM_DAT = 7;
const int NUM_MIS = 6;

// Define a base random number generator and initialize it with a seed.
boost::minstd_rand baseGen(std::time(0));

// Define distribution U[0,1) [double values]
boost::uniform_real<> uniDblUnit(0,1);

// Define a random variate generator using our base generator and distribution
boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uniDblGen(baseGen, uniDblUnit);

// amino acid residue types
const int ALA = 0;
const int CYS = 1;
const int ASP = 2;
const int GLU = 3;
const int PHE = 4;
const int GLY = 5;
const int HIS = 6;
const int ILE = 7;
const int LYS = 8;
const int LEU = 9;
const int MET = 10;
const int ASN = 11;
const int PRO = 12;
const int GLN = 13;
const int ARG = 14;
const int SER = 15;
const int THR = 16;
const int VAL = 17;
const int TRP = 18;
const int TYR = 19;

// amino acide residue code to three letter format
string seq3[20] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

// amino acide residue code to one letter format
string seq[20] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

// secondary structure types
const int ALPHA_HELIX = 0;
const int THREE_HELIX = 1;
const int FIVE_HELIX = 2;
const int ISOLATED_STRAND = 3;
const int EXTENDED_STRAND = 4;
const int HBONDED_TURN = 5;
const int NONHBONDED_BEND = 6;
const int RANDOM_COIL = 7;

// secondary structure code to one letter format
string sec[8] = {"H", "G", "I", "B", "E", "T", "S", "C"};

// point3d object
struct point3d {
    double x;
    double y;
    double z;
};

// pdbInfo object
struct pdbInfo {
    int id;
    int aa;
    point3d ca;
    point3d sc;
};

// poseInfo object
struct poseInfo {
    int id;
    int aa;
    int ss;
    double bca;
    double tao;
    double theta;
    double bsc;
    double phi;
    double delta;
};

// contactInfo object
struct contactInfo {
    int ri;
    int rj;
    double wt;
};

// vector to store the conatcts
vector<contactInfo> rrContact;

//max and min fragment size
int MIN_FRG_LEN = 1;
int MAX_FRG_LEN = 15;

// weights of energy terms
double w_sc_sc = 1.0;
double w_sc_bb = 2.73684;
double w_bb_bb = 0.06833;
double w_ri_rj = 3.0;

// energyInfo object
struct energyInfo {
    double sc_sc;
    double sc_bb;
    double bb_bb;
    double ri_rj;
    double total;
};

// sseInfo object
struct sseInfo {
    int start;
    int end;
    int type;
};

// amino acid string defined by Liwo et. al. 1997 J. Comp. Chem. I
string liwoAminoStr = "CMFILVWYAGTSQNEDHRKP";

// epsilon0 parameter (lower triangle and digonal) for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoEpsilon0[20][20] = {
    {1.050, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.030, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020},
    {1.260, 1.450, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.030, 0.030, 0.030, 0.020, 0.020, 0.020, 0.030, 0.020},
    {1.190, 1.340, 1.270, 0.010, 0.010, 0.010, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.010, 0.010, 0.020, 0.010, 0.010, 0.010},
    {1.300, 1.470, 1.410, 1.580, 0.010, 0.010, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {1.250, 1.510, 1.400, 1.590, 1.550, 0.010, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {1.170, 1.380, 1.310, 1.520, 1.500, 1.400, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {0.990, 1.170, 1.150, 1.210, 1.180, 1.100, 0.970, 0.020, 0.020, 0.020, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020},
    {0.920, 1.150, 1.050, 1.220, 1.180, 1.040, 0.870, 0.810, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {0.980, 1.200, 0.990, 1.240, 1.260, 1.190, 0.770, 0.810, 1.020, 0.010, 0.010, 0.010, 0.030, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {0.980, 1.030, 0.840, 1.060, 1.130, 1.010, 0.710, 0.720, 0.820, 0.560, 0.010, 0.010, 0.020, 0.020, 0.010, 0.020, 0.020, 0.020, 0.000, 0.010},
    {0.800, 0.910, 0.760, 0.980, 0.900, 0.870, 0.560, 0.580, 0.640, 0.550, 0.430, 0.010, 0.020, 0.020, 0.020, 0.010, 0.020, 0.020, 0.000, 0.010},
    {0.770, 0.860, 0.680, 0.900, 0.930, 0.830, 0.510, 0.520, 0.590, 0.470, 0.450, 0.280, 0.020, 0.020, 0.010, 0.020, 0.020, 0.020, 0.000, 0.010},
    {0.810, 1.020, 0.720, 0.950, 0.980, 0.850, 0.590, 0.630, 0.750, 0.330, 0.350, 0.260, -0.280, 0.060, 0.030, 0.020, 0.030, 0.040, 0.010, 0.020},
    {0.730, 0.950, 0.700, 0.870, 1.000, 0.890, 0.600, 0.600, 0.770, 0.490, 0.380, 0.380, 0.530, 0.660, 0.010, 0.040, 0.030, 0.050, 0.000, 0.020},
    {0.640, 0.810, 0.530, 0.880, 0.790, 0.720, 0.520, 0.510, 0.470, -0.060, 0.200, 0.040, -0.230, -0.020, -1.580, 0.100, 0.030, 0.040, 0.040, 0.020},
    {0.680, 0.640, 0.520, 0.790, 0.680, 0.620, 0.530, 0.570, 0.510, 0.230, 0.290, 0.120, -0.120, 0.270, -0.930, -0.660, 0.030, 0.030, 0.050, 0.020},
    {0.910, 1.050, 0.920, 0.940, 0.980, 0.830, 0.820, 0.760, 0.650, 0.560, 0.570, 0.490, 0.380, 0.590, 0.420, 0.620, 0.800, 0.030, 0.000, 0.020},
    {0.580, 0.870, 0.690, 0.960, 0.950, 0.750, 0.660, 0.670, 0.530, 0.380, 0.430, 0.400, 0.360, 0.330, 1.010, 1.000, 0.490, -0.020, 0.090, 0.020},
    {0.590, 0.810, 0.550, 0.960, 0.970, 0.850, 0.500, 0.620, 0.680, -0.010, 0.000, -0.010, -0.020, 0.000, 1.300, 1.090, -0.010, -0.480, -11.960, 0.020},
    {0.820, 0.950, 0.810, 0.980, 1.000, 0.920, 0.770, 0.790, 0.740, 0.690, 0.570, 0.580, 0.620, 0.620, 0.420, 0.420, 0.610, 0.530, 0.560, 0.820}
};

// sigma_0 parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoSigma0[20] = {2.4366, 2.3204, 2.5098, 2.5933, 2.2823, 2.3359, 2.3409, 2.5919, 2.7249, 2.8905, 2.4984, 2.6954, 2.723, 2.6269, 2.3694, 2.4471, 2.6047, 2.7251, 1.6947, 2.1346};

// r_0 parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoR0[20] = {2.1574, 5.7866, 1.4061, 1.2819, 6.3367, 2.5197, 3.357, 4.4859, 0.2712, 3.3121, 3.5449, 0.7634, 3.332, 1.1813, 1.8119, 2.2432, 3.0723, 3.777, 9.2904, 4.8607};

// sigma_double_bar_over_sigma_tee_square parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoSigmaDoubleBarOverSigmaTeeSquare[20] = {1.809, 2.6006, 1.9262, 2.5707, 3.964, 1.0429, 3.6263, 3.2406, 8.0078, 2.3636, 4.4303, 2.0433, 1.7905, 2.6172, 6.6061, 1.6795, 2.2451, 2.0347, 7.5089, 5.9976};

// chi_prime parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoChiPrime[20] = {0.0333, -0.0025, 0.309, 0.4904, 0.0491, 0.2238, 0.1351, 0.0897, 0.579, 0.0749, 0.0968, 0.2732, -0.1105, 0.296, 0.2624, -0.0029, 0.0236, 0.077, 0.0731, 0.1177};

// alpha parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoAlpha[20] = {0.1052, 0.0299, 0.025, -0.0266, 0.0801, -0.1277, 0.0589, 0.0664, 0.0115, 0.1108, 0.0878, -0.0064, -0.019, 0.0505, 0.0062, -0.0348, -0.0264, 0.0679, 0.0549, 0.0438};

// bsc parameter for side chain length by Levitt M. 1976 JMB (Table 1)
double levittBsc[20] = {0.77, 1.38, 1.99, 2.63, 2.97, 0.0, 2.76, 1.83, 2.94, 2.08, 2.34, 1.98, 1.42, 2.58, 3.72, 1.28, 1.43, 1.49, 3.58, 3.36};

//global variables
char jobId[200] = "";
char faFile[200] = "";
char ssFile[200] = "";
char cmFile[200] = "";
char moFile[200] = "";
int numCycles = 100;
int numDecoys = 100;
char natFile[200] = "";
char statFile[200] = "";

//indicators
bool jId = false;
bool fFile = false;
bool sFile = false;
bool cFile = false;
bool mFile = false;
bool nFile = false;

void parseNextItem(int argc, char ** argv, int & i);
void parseCommandLine(int argc, char ** argv);

int getAA(const char * aa);
int getSS(const char * ss);

double getDistance(point3d & p1, point3d &p2);
double getDotProduct(point3d & p1, point3d &p2);
point3d getCrossProduct(point3d & p1, point3d &p2);
double getNorm(point3d & p);
point3d getDifference(point3d & p1, point3d &p2);
point3d getMidpoint(point3d & p1, point3d &p2);
point3d getUnit(point3d & p);
double getAngle(point3d & p1, point3d &p2, point3d &p3);
double getDihedral(point3d & p1, point3d &p2, point3d &p3, point3d &p4);
point3d getCentroid(vector<point3d> &pointCloud);
void setCoordinate(point3d &c0, point3d &c1, point3d &c2, point3d &c3, double alpha, double tau, double normw);

void loadPdb(char *filename, vector<pdbInfo> &pdb);
void loadFasta(char *filename, vector<int> &aa);
void loadSec(char *filename, vector<int> &ss);
void loadCmap(char *filename, MDArray<double> &cm);

void pdb2pose(vector<pdbInfo> &pdb, vector<int> &ss, vector<poseInfo> &pose);
void pose2pdb(vector<poseInfo> &pose, vector<pdbInfo> &pdb);
void sample2pose(MDArray<double> &sample, vector<poseInfo> &pose);
void pose2sample(vector<poseInfo> &pose, MDArray<double> &sample);

bool isHydrophobic(int aa);
void getRrContact(MDArray<double> &cm, int end);
double getUpperBoundHarmonic(double d, double bound);
double getFade(double d, double lb, double ub, double z, double w);
double getLiwoEpsilon0(int ai, int aj);

double getScScEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
double getScBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
double getBbBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
double getRiRjEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
energyInfo getWeightedEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn);
void showEnergy(energyInfo &e);

double getRmsd(vector<pdbInfo> &pdb1, vector<pdbInfo> &pdb2, string mode = "ca");
void getAlignment(vector<vector<double> > &  A, vector<vector<double> > & B, double rot[3][3], double trans[3]);
void quat2mat(double q[4], double mat[3][3]);
void vApply (double m[3][3], const double a[3], double b[3]);

void assembleFoldon(vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow);
void doSimulatedAnneling(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp);
bool acceptMove(energyInfo &prev, energyInfo &curr, double temperature);

void writePdb(vector<pdbInfo> &pdb, char * filename);

/*************************************************************************
 * Name        : main
 * Purpose     : United residue protein folding via stepwise, conditional sampling
 * Arguments   : int argc, char ** argv
 * Return Type : int
 *************************************************************************/
int main(int argc, char ** argv) {
    cout << endl;
    cout << "###########################################################################" << endl;
    cout << "#                                UniCon3D                                 #" << endl;
    cout << "#                              Version: 1.0                               #" << endl;
    cout << "#         De novo protein structure prediction using united-residue       #" << endl;
    cout << "#        conformational search via stepwise, probabilistic sampling       #" << endl;
    cout << "#         Copyright (C) Debswapna Bhattacharya and Jianlin Cheng          #" << endl;
    cout << "#     Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory      #" << endl;
    cout << "#              University of Missouri, Columbia MO 65211                  #" << endl;
    cout << "###########################################################################" << endl;
    
    parseCommandLine(argc, argv);
    
    // load model
    DBN dbn;
    dbn.load(moFile);
    dbn.randomGen->get_rand();
    
    // load amino acid
    vector<int> aa;
    loadFasta(faFile, aa);
    
    // load secondary structure
    vector<int> ss;
    loadSec(ssFile, ss);
    
    if (aa.size() != ss.size()) {
        cout << "Error! size mismatch in fasta, ss" << endl;
        exit(0);
    }
    
    //load contact matrix
    MDArray<double> cm;
    cm.set_shape(aa.size(), aa.size());
    loadCmap(cmFile, cm);
    
    // load native pdb if provided
    vector<pdbInfo> native;
    vector<poseInfo> nativePose;
    if (nFile) {
        loadPdb(natFile, native);
        if (aa.size() != native.size()) {
            cout << "Error! size mismatch in native structure" << endl;
            exit(0);
        }
        pdb2pose(native, ss, nativePose);
    }
		
		MAX_FRG_LEN = min(MAX_FRG_LEN, aa.size());
			
    // extract secondary structural elements
    vector<sseInfo> sse;
    bool sseStarted = false;
    sseInfo sseData;
    for (int i = 0; i < ss.size() - 1; i++) {
        if (!sseStarted) {
            sseData.start = i;
            sseData.type = ss[i];
            sseStarted = true;
        }
        if (sseStarted && ss[i] != ss[i+1]) {
            sseData.end = i;
            sse.push_back(sseData);
            sseStarted = false;
        }
    }
    // for the last sse
    sseData.type = ss[ss.size() - 1];
    sseData.end = ss.size() - 1;
    sse.push_back(sseData);
    
    
    // extract foldon units
    vector<sseInfo> foldon;
    bool foldonStarted = false;
    sseInfo foldonData;
    for (int i = 0; i < sse.size(); i++) {
        if (!foldonStarted) {
            foldonData.start = sse[i].start;
            foldonStarted = true;
        }
        // consider regular secondary structural elements of at least 20 residues long
        if (foldonStarted && (sse[i].end - foldonData.start + 1) > min(MIN_FOLDON_LEN, aa.size() - 1) && (sse[i].type == ALPHA_HELIX || sse[i].type == THREE_HELIX || sse[i].type == FIVE_HELIX || sse[i].type == ISOLATED_STRAND || sse[i].type == EXTENDED_STRAND)) {
	    			foldonData.end = sse[i].end;
            foldonData.type = sse[i].type;
            foldon.push_back(foldonData);
            foldonStarted = false;
        }
    }
    foldon[foldon.size() - 1].end = sse[sse.size()-1].end;
    
    // populate sample corresponding to foldon units
    vector<MDArray<double> > foldonSample;
    for (int i = 0; i < foldon.size(); i++) {
        MDArray<double> foldonSampleData;
        int foldonSize = foldon[i].end + 1;
        foldonSampleData.set_shape(2 * foldonSize, NUM_DAT);
        for (int j = 0; j <= foldon[i].end; j++) {
            // for backbone
            foldonSampleData.set(2 * j, 0, 0);
            foldonSampleData.set(2 * j, 1, 0);
            foldonSampleData.set(2 * j, 2, aa[j]);
            foldonSampleData.set(2 * j, 3, ss[j]);
            foldonSampleData.set(2 * j, 4, CA2CA);
            foldonSampleData.set(2 * j, 5, PI);
            foldonSampleData.set(2 * j, 6, PI);
            // for side chain
            foldonSampleData.set(2 * j + 1, 0, 1);
            foldonSampleData.set(2 * j + 1, 1, 0);
            foldonSampleData.set(2 * j + 1, 2, aa[j]);
            foldonSampleData.set(2 * j + 1, 3, ss[j]);
            foldonSampleData.set(2 * j + 1, 4, levittBsc[aa[j]]);
            foldonSampleData.set(2 * j + 1, 5, 0);
            foldonSampleData.set(2 * j + 1, 6, 0);
        }
        foldonSample.push_back(foldonSampleData);
    }
    
    // assemble foldon units via stepwise addition
    char pdbFoldonFile[200];
    vector<pdbInfo> pdbFoldon;
    vector<poseInfo> poseFoldon;
    

    sprintf(statFile, "%s_stats.txt", jobId);
    
    cout << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    cout << "Job Id                   : " << jobId << endl;
    cout << "Protein Length           : " << aa.size() << endl;
    cout << "Amino Acid Sequence      : ";
    for (int i = 0; i < aa.size(); i++) {
        cout << seq[aa[i]];
    }
    cout << endl;
    cout << "Secondary Structure      : ";
    for (int i = 0; i < ss.size(); i++) {
        cout << sec[ss[i]];
    }
    cout << endl;
    cout << "Number of Decoys         : " << numDecoys << endl;
    cout << "Monte Carlo Cycle Factor : " << numCycles/100 << endl;
    cout << "Native Structure         : ";
    if (nFile) {
        cout << "Available" << endl;
    }
    else {
        cout << "Not Available" << endl;
    }
    cout << "Number of Foldon Units   : " << foldon.size() << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
    
    ofstream stat(statFile);
    if (!stat.is_open()) {
        cout << "Error! folding statistics file can not open " << statFile << endl;
        exit(0);
    }
    // buffer statistics
    char bufStat[1000];

    for (int n = 1; n <= numDecoys; n++) {
        cout << "UniCon3D.io: Generating decoy " << n << " of " << numDecoys << endl;
        for (int i = 0; i < foldon.size(); i++) {
            
            // get residue-residue contacts from native pdb
            getRrContact(cm, foldon[i].end);
            
            if (i == 0) {
                cout << endl <<"UniCon3D.sampling: Sampling foldon unit " << i + 1 << " [residue : 1-" << foldon[i].end + 1 << "]"<< endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            else {
                cout << endl << "UniCon3D.sampling: Conditional sampling foldon unit " << i + 1 << " [residue : " << foldon[i-1].end + 2 << "-" << foldon[i].end + 1 << "]" << endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            cout << "UniCon3D.sampling: Performing simulated annealing energy minimization" << endl;
            assembleFoldon(foldonSample, i, dbn, pdbFoldon, poseFoldon);

            cout << endl;
            //sprintf(pdbFoldonFile, "%s_D%.6d_F%.6d.pdb", jobId, n, i+1);
            //if (saveInd == 0) {
            //    writePdb(pdbFoldon, pdbFoldonFile);
            //}
        }
        sprintf(pdbFoldonFile, "%s_%.6d.pdb", jobId, n);
        writePdb(pdbFoldon, pdbFoldonFile);
        
        cout << endl << "UniCon3D.io: Finished generating decoy " << pdbFoldonFile << endl;
        
        
        double rmsd;
        if (nFile) {
            rmsd = getRmsd(pdbFoldon, native);
            cout << "UniCon3D.io: Ca_rmsd to native structure = " << rmsd << endl;
        }
        // prepare to write to folding statistics file
        energyInfo e = getWeightedEnergy(pdbFoldon, poseFoldon, dbn);
        if (nFile) {
            sprintf(bufStat,"Decoy_name %s E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f Ca_rmsd %8.3f", pdbFoldonFile, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total, rmsd);
            stat << bufStat << endl;
        }
        else {
            sprintf(bufStat,"Decoy_name %s E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f", pdbFoldonFile, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total);
            stat << bufStat << endl;
        }
        cout << "UniCon3D.io: Folding statistics written to " << statFile << endl;
        cout << "---------------------------------------------------------------------------" << endl;
    }
    return EXIT_SUCCESS;
}

/*************************************************************************
 * Name        : parseNextItem
 * Purpose     : parse next item in command line argument
 * Arguments   : int argc, char ** argv, int & i
 * Return Type : void
 *************************************************************************/
void parseNextItem(int argc, char ** argv, int & i) {
    if (strncmp(argv[i], "-i", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No job id provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        jId = true;
        strcpy(jobId, argv[++i]);
    }
    else if (strncmp(argv[i], "-f", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No fasta file provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        fFile = true;
        strcpy(faFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-s", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No secondary structure file provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        sFile = true;
        strcpy(ssFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-c", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No contact map file provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        cFile = true;
        strcpy(cmFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-m", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No model file provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        mFile = true;
        strcpy(moFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-d", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No number of decoys provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        numDecoys = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-x", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No monte carlo cycle factor provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        numCycles = numCycles * atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-n", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No native structure provided" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        nFile = true;
        strcpy(natFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-h", 2) == 0) {
        cout << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    else {
        cout << endl;
        cout << "Error! Invalid option" << endl << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    i++;
}

/*************************************************************************
 * Name        : parseCommandLine
 * Purpose     : parse command line arguments
 * Arguments   : int argc, char ** argv
 * Return Type : void
 *************************************************************************/
void parseCommandLine(int argc, char ** argv) {
    int i = 1;
    while (i < argc)
        parseNextItem(argc, argv, i);
    if (!jId) {
        cout << endl;
        cout << "Error! Job id must be provided" << endl << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    if (!fFile) {
        cout << endl;
        cout << "Error! Fasta file must be provided" << endl << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    else {
        ifstream fin( faFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Fasta file not present" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
    }
    if (!sFile) {
        cout << endl;
        cout << "Error! Secondary structure file must be provided" << endl << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    else {
        ifstream fin( ssFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Secondary structure file not present" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
    }
    if (!cFile) {
        cout << endl;
        cout << "Error! Contact map file must be provided" << endl << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    else {
        ifstream fin( cmFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Contact map file not present" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
    }
    if (!mFile) {
        cout << endl;
        cout << "Error! Model file must be provided" << endl << endl;
        cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
        cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
        cout << "   -i id     : job id                                     : mandatory" << endl;
        cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
        cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
        cout << "   -c cmap   : contact map file                           : mandatory" << endl;
        cout << "   -m model  : model file                                 : mandatory" << endl;
        cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
        cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
        cout << "   -n native : native structure for comparison            : optional" << endl;
        cout << "   -h help   : this message" << endl;
        exit(0);
    }
    else {
        ifstream fin( moFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Model file not present" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
    }
    if (nFile) {
        ifstream fin( natFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Native structure not present" << endl << endl;
            cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
            cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
            cout << "   -i id     : job id                                     : mandatory" << endl;
            cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
            cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
            cout << "   -c cmap   : contact map file                           : mandatory" << endl;
            cout << "   -m model  : model file                                 : mandatory" << endl;
            cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
            cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
            cout << "   -n native : native structure for comparison            : optional" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
    }
}

/*************************************************************************
 * Name        : getAA
 * Purpose     : convert AA name to a numerical code
 * Arguments   : const char * aa
 * Return Type : int
 *************************************************************************/
int getAA(const char * aa) {
    if (strlen(aa) == 3) {
        if (strcmp(aa, "ALA") == 0)
            return (ALA);
        else if (strcmp(aa, "ARG") == 0)
            return (ARG);
        else if (strcmp(aa, "ASN") == 0)
            return (ASN);
        else if (strcmp(aa, "ASP") == 0)
            return (ASP);
        else if (strcmp(aa, "CYS") == 0)
            return (CYS);
        else if (strcmp(aa, "GLN") == 0)
            return (GLN);
        else if (strcmp(aa, "GLU") == 0)
            return (GLU);
        else if (strcmp(aa, "GLY") == 0)
            return (GLY);
        else if (strcmp(aa, "HIS") == 0)
            return (HIS);
        else if (strcmp(aa, "ILE") == 0)
            return (ILE);
        else if (strcmp(aa, "LEU") == 0)
            return (LEU);
        else if (strcmp(aa, "LYS") == 0)
            return (LYS);
        else if (strcmp(aa, "MET") == 0)
            return (MET);
        else if (strcmp(aa, "PHE") == 0)
            return (PHE);
        else if (strcmp(aa, "PRO") == 0)
            return (PRO);
        else if (strcmp(aa, "SER") == 0)
            return (SER);
        else if (strcmp(aa, "THR") == 0)
            return (THR);
        else if (strcmp(aa, "TRP") == 0)
            return (TRP);
        else if (strcmp(aa, "TYR") == 0)
            return (TYR);
        else if (strcmp(aa, "VAL") == 0)
            return (VAL);
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else if (strlen(aa) == 1) {
        if (aa[0] == 'A')
            return ALA;
        else if (aa[0] == 'C')
            return CYS;
        else if (aa[0] == 'D')
            return ASP;
        else if (aa[0] == 'E')
            return GLU;
        else if (aa[0] == 'F')
            return PHE;
        else if (aa[0] == 'G')
            return GLY;
        else if (aa[0] == 'H')
            return HIS;
        else if (aa[0] == 'I')
            return ILE;
        else if (aa[0] == 'K')
            return LYS;
        else if (aa[0] == 'L')
            return LEU;
        else if (aa[0] == 'M')
            return MET;
        else if (aa[0] == 'N')
            return ASN;
        else if (aa[0] == 'P')
            return PRO;
        else if (aa[0] == 'Q')
            return GLN;
        else if (aa[0] == 'R')
            return ARG;
        else if (aa[0] == 'S')
            return SER;
        else if (aa[0] == 'T')
            return THR;
        else if (aa[0] == 'V')
            return VAL;
        else if (aa[0] == 'W')
            return TRP;
        else if (aa[0] == 'Y')
            return TYR;
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else {
        cout << "Error! Invalid amino acid " << aa << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : getSS
 * Purpose     : convert SS name to a numerical code
 * Arguments   : const char * ss
 * Return Type : int
 *************************************************************************/
int getSS(const char * ss) {
    if (strlen(ss) == 1) {
        if (strcmp(ss, "H") == 0)
            return (ALPHA_HELIX);
        else if (strcmp(ss, "G") == 0)
            return (THREE_HELIX);
        else if (strcmp(ss, "I") == 0)
            return (FIVE_HELIX);
        else if (strcmp(ss, "B") == 0)
            return (ISOLATED_STRAND);
        else if (strcmp(ss, "E") == 0)
            return (EXTENDED_STRAND);
        else if (strcmp(ss, "T") == 0)
            return (HBONDED_TURN);
        else if (strcmp(ss, "S") == 0)
            return (NONHBONDED_BEND);
        else if (strcmp(ss, "C") == 0)
            return (RANDOM_COIL);
        else {
            cout << "Error! Invalid secondary structure " << ss << endl;
            exit(0);
        }
    }
    else {
        cout << "Error! Invalid secondary structure " << ss << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadPdb
 * Purpose     : loads a pdb file into a vector of pdbInfo object
 * Arguments   : char *filename, vector<pdbInfo> &pdb
 * Return Type : void
 *************************************************************************/
void loadPdb(char *filename, vector<pdbInfo> &pdb) {
    string line, str;
    string atom ("ATOM ");
    int prevRes = -999999;
    point3d caAtom;
    point3d scAtom;
    vector<point3d> scAtomCloud;
    pdbInfo pdbData;
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            if(line.compare(0, atom.length(), atom)==0) {
                int res = atoi(line.substr(22, 4).c_str());
                int anmae = getAA(line.substr(17, 3).c_str());
                // seek for the next residue
                if (res != prevRes) {
                    if (prevRes != -999999) {
                        // inser the data collected so far
                        pdbData.ca = caAtom;
                        pdbData.sc = getCentroid(scAtomCloud);
                        pdb.push_back(pdbData);
                    }
                    prevRes = res;
                    scAtomCloud.clear();
                    pdbData.id = res;
                    pdbData.aa = anmae;
                }
                // consider the first alternate location id (i.e. A) if present
                if( line.compare(16, 1, " ") == 0 || line.compare(16, 1, "A") == 0 ) {
                    // get the CA atom coordinate
                    if( line.compare(12, 4, "CA  ") == 0 || line.compare(12, 4, " CA ") == 0 || line.compare(12, 4, "  CA") == 0 ) {
                        
                        caAtom.x = atof(line.substr(30, 8).c_str());
                        caAtom.y = atof(line.substr(38, 8).c_str());
                        caAtom.z = atof(line.substr(46, 8).c_str());
                        scAtomCloud.push_back(caAtom);
                    }
                    // get the sidechain heavy atoms
                    if(line.compare(12, 4, "N   ") != 0 && line.compare(12, 4, " N  ") != 0 && line.compare(12, 4, "  N ") != 0 &&
                       line.compare(12, 4, "C   ") != 0 && line.compare(12, 4, " C  ") != 0 && line.compare(12, 4, "  C ") != 0 &&
                       line.compare(12, 4, "O   ") != 0 && line.compare(12, 4, " O  ") != 0 && line.compare(12, 4, "  O ") != 0 &&
                       line.compare(12, 4, "CA  ") != 0 && line.compare(12, 4, " CA ") != 0 && line.compare(12, 4, "  CA") != 0 &&
                       line.compare(77, 1, "H") != 0) {
                        
                        scAtom.x = atof(line.substr(30, 8).c_str());
                        scAtom.y = atof(line.substr(38, 8).c_str());
                        scAtom.z = atof(line.substr(46, 8).c_str());
                        scAtomCloud.push_back(scAtom);
                    }
                }
            }
        }
        fin.close();
        // for the last residue
        pdbData.ca = caAtom;
        pdbData.sc = getCentroid(scAtomCloud);
        pdb.push_back(pdbData);
    }
    else {
        cout << "Error! pdb file can not open " << filename << endl;
        exit(0);
    }
    
}

/*************************************************************************
 * Name        : loadFasta
 * Purpose     : loads a fasta file into a vector of int object
 * Arguments   : char *filename, vector<int> &aa
 * Return Type : void
 *************************************************************************/
void loadFasta(char *filename, vector<int> &aa) {
    string line, str;
    string header (">");
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            if(line.length() != 0 && line.compare(0, header.length(), header) != 0) {
                for (int i = 0; i < line.length(); i++) {
                    int aa_code = getAA(line.substr(i,1).c_str());
                    aa.push_back(aa_code);
                }
            }
        }
        fin.close();
    }
    else {
        cout << "Error! fasta file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadSec
 * Purpose     : loads a secondary structure file into a vector of int object
 * Arguments   : char *filename, vector<int> &ss
 * Return Type : void
 *************************************************************************/
void loadSec(char *filename, vector<int> &ss) {
    string line, str;
    string header (">");
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            if(line.length() != 0 && line.compare(0, header.length(), header) != 0) {
                for (int i = 0; i < line.length(); i++) {
                    int ss_code = getSS(line.substr(i,1).c_str());
                    ss.push_back(ss_code);
                }
            }
        }
        fin.close();
    }
    else {
        cout << "Error! secondary structure file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadCmap
 * Purpose     : loads a contact map into a MDArray object
 * Arguments   : char *filename, vector<int> &sa
 * Return Type : void
 *************************************************************************/
void loadCmap(char *filename, MDArray<double> &cm) {
    ifstream fin (filename);
    double val;
    if (fin.is_open()) {
        for (int i = 0; i < cm.get_shape()[0]; i++) {
            for (int j = 0; j < cm.get_shape()[0]; j++) {
                fin >> val;
                cm.set(i, j, val);
            }
        }
        fin.close();
    }
    else {
        cout << "Error! contact map file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : getDistance
 * Purpose     : gets distance between two points
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : double
 *************************************************************************/
double getDistance(point3d & p1, point3d &p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

/*************************************************************************
 * Name        : getDotProduct
 * Purpose     : gets the dot product (i.e. scalar product) for two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : double
 *************************************************************************/
double getDotProduct(point3d & p1, point3d &p2) {
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

/*************************************************************************
 * Name        : getCrossProduct
 * Purpose     : gets the cross product (i.e. vector product) for two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getCrossProduct(point3d & p1, point3d &p2) {
    point3d p;
    p.x = p1.y * p2.z - p1.z * p2.y;
    p.y = p1.z * p2.x - p1.x * p2.z;
    p.z = p1.x * p2.y - p1.y * p2.x;
    return p;
}

/*************************************************************************
 * Name        : getNorm
 * Purpose     : gets the norm (i.e. length) of a vector
 * Arguments   : point3d & p
 * Return Type : double
 *************************************************************************/
double getNorm(point3d & p) {
    return sqrt( p.x * p.x + p.y * p.y + p.z * p.z);
}

/*************************************************************************
 * Name        : getDifference
 * Purpose     : gets the difference vector between two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getDifference(point3d & p1, point3d &p2) {
    point3d p;
    p.x = p2.x - p1.x;
    p.y = p2.y - p1.y;
    p.z = p2.z - p1.z;
    return p;
}

/*************************************************************************
 * Name        : getMidpoint
 * Purpose     : gets the midpoint vector between two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getMidpoint(point3d & p1, point3d &p2) {
    point3d p;
    p.x = (p2.x + p1.x) / 2.0;
    p.y = (p2.y + p1.y) / 2.0;
    p.z = (p2.z + p1.z) / 2.0;
    return p;
}

/*************************************************************************
 * Name        : getUnit
 * Purpose     : gets the unit vector for a vector
 * Arguments   : point3d & p
 * Return Type : point3d
 *************************************************************************/
point3d getUnit(point3d & p) {
    point3d u;
    double norm = getNorm(p);
    if (norm == 0.0) {
        u.x = 0.0;
        u.y = 0.0;
        u.z = 0.0;
    }
    else {
        u.x = p.x / norm;
        u.y = p.y / norm;
        u.z = p.z / norm;
    }
    return u;
}

/*************************************************************************
 * Name        : getAngle
 * Purpose     : gets the angle between three points
 * Arguments   : point3d & p1, point3d &p2, point3d &p3
 * Return Type : double
 *************************************************************************/
double getAngle(point3d & p1, point3d &p2, point3d &p3) {
    double acc = 0.0;
    acc = (p2.x - p1.x) * (p2.x - p3.x) + (p2.y - p1.y) * (p2.y - p3.y) + (p2.z - p1.z) * (p2.z - p3.z);
    double d1 = getDistance(p1, p2);
    double d2 = getDistance(p3, p2);
    if (d1 == 0 || d2 == 0) {
        return 0;
    }
    acc = acc / (d1 * d2);
    if (acc > 1.0)
        acc = 1.0;
    else if (acc < -1.0)
        acc = -1.0;
    acc = acos(acc);
    return acc;
}

/*************************************************************************
 * Name        : getDihedral
 * Purpose     : gets the dihedral between four points
 * Arguments   : point3d & p1, point3d &p2, point3d &p3, point3d &p4
 * Return Type : double
 *************************************************************************/
double getDihedral(point3d & p1, point3d &p2, point3d &p3, point3d &p4) {
    point3d q, r, s, t, u, v, z;
    double acc, w;
    z.x = z.y = z.z = 0.0;
    q = getDifference(p1, p2);
    r = getDifference(p3, p2);
    s = getDifference(p4, p3);
    t = getCrossProduct(q, r);
    u = getCrossProduct(s, r);
    v = getCrossProduct(u, t);
    w = getDotProduct(v, r);
    acc = getAngle(t, z, u);
    if (w < 0)
        acc = -acc;
    return (acc);
}

/*************************************************************************
 * Name        : getCentroid
 * Purpose     : gets centroid of a cloud of points
 * Arguments   : vector<point3d> &pointCloud
 * Return Type : point3d
 *************************************************************************/
point3d getCentroid(vector<point3d> &pointCloud) {
    point3d centroid;
    centroid.x = 0.0;
    centroid.y = 0.0;
    centroid.z = 0.0;
    for (int i = 0; i < pointCloud.size(); i++) {
        centroid.x += pointCloud[i].x / pointCloud.size();
        centroid.y += pointCloud[i].y / pointCloud.size();
        centroid.z += pointCloud[i].z / pointCloud.size();
    }
    return centroid;
}

/*************************************************************************
 * Name        : setCoordinate
 * Purpose     : set coordinate of a point based on dihedral alpha, bond angle tau, bond length normw
 * Arguments   : point3d &c0, point3d &c1, point3d &c2, point3d &c3, double alpha, double tau, double normw
 * Return Type : void
 *************************************************************************/
void setCoordinate(point3d &c0, point3d &c1, point3d &c2, point3d &c3, double alpha, double tau, double normw) {
    double	u1, u2, u3, v1, v2, v3, norm, pvuv1, pvuv2, pvuv3, pvvuv1, pvvuv2, pvvuv3, nsa, nca, nct;
    u1 = (c1.x - c0.x);
    u2 = (c1.y - c0.y);
    u3 = (c1.z - c0.z);
    v1 = (c2.x - c1.x);
    v2 = (c2.y - c1.y);
    v3 = (c2.z - c1.z);
    norm = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
    v1 /=  norm;
    v2 /=  norm;
    v3 /=  norm;
    pvuv1 = u2 * v3 - u3 * v2;
    pvuv2 = u3 * v1 - u1 * v3;
    pvuv3 = u1 * v2 - u2 * v1;
    norm = sqrt(pvuv1 * pvuv1 + pvuv2 * pvuv2 + pvuv3 * pvuv3);
    pvuv1 /= norm;
    pvuv2 /= norm;
    pvuv3 /= norm;
    pvvuv1 = v2 * pvuv3 - v3 * pvuv2;
    pvvuv2 = v3 * pvuv1 - v1 * pvuv3;
    pvvuv3 = v1 * pvuv2 - v2 * pvuv1;
    norm = sqrt(pvvuv1 * pvvuv1 + pvvuv2 * pvvuv2 + pvvuv3 * pvvuv3);
    pvvuv1 /= norm;
    pvvuv2 /= norm;
    pvvuv3 /= norm;
    nca = cos(alpha);
    nsa = sin(alpha);
    nct = tan(tau - PI/2);
    u1 = nca * (-pvvuv1) + nsa * pvuv1 + v1 * nct;
    u2 = nca * (-pvvuv2) + nsa * pvuv2 + v2 * nct;
    u3 = nca * (-pvvuv3) + nsa * pvuv3 + v3 * nct;
    norm = sqrt(u1 * u1 + u2 * u2 + u3 * u3);
    u1 = u1 * normw/norm;
    u2 = u2 * normw/norm;
    u3 = u3 * normw/norm;
    c3.x = u1 + c2.x;
    c3.y = u2 + c2.y;
    c3.z = u3 + c2.z;
}

/*************************************************************************
 * Name        : pdb2pose
 * Purpose     : converts pdb to pose
 * Arguments   : vector<pdbInfo> &pdb, vector<int> &ss, vector<poseInfo> &pose
 * Return Type : void
 *************************************************************************/
void pdb2pose(vector<pdbInfo> &pdb, vector<int> &ss, vector<poseInfo> &pose) {
    pose.clear();
    poseInfo poseData;
    //define a virtual N terminus
    point3d N;
    N.x = pdb[0].ca.x - cos(PI) * CA2CA;
    N.y = pdb[0].ca.y - sin(PI) * CA2CA;
    N.z = pdb[0].ca.z;
    //define a virtual C terminus
    point3d C;
    C.x = pdb[pdb.size()-1].ca.x + cos(PI) * CA2CA;
    C.y = pdb[pdb.size()-1].ca.y + sin(PI) * CA2CA;
    C.z = pdb[pdb.size()-1].ca.z;
    for (int i = 0; i < pdb.size(); i++) {
        poseData.id = pdb[i].id;
        poseData.aa = pdb[i].aa;
        poseData.ss = ss[i];
        // calculate bca
        if (i == 0) {
            poseData.bca = CA2CA;
        }
        else {
            poseData.bca = getDistance(pdb[i-1].ca, pdb[i].ca);
        }
        // calculate tao
        if (i == 0 || i == pdb.size() - 1 || i == pdb.size() - 2) {
            poseData.tao = 2 * PI;
        }
        else {
            poseData.tao = getDihedral(pdb[i-1].ca, pdb[i].ca, pdb[i+1].ca, pdb[i+2].ca);
        }
        // calculate theta
        if (i == 0 || i == pdb.size() - 1) {
            poseData.theta = 2 * PI;
        }
        else {
            poseData.theta = getAngle(pdb[i-1].ca, pdb[i].ca, pdb[i+1].ca);
        }
        //calculate bsc
        poseData.bsc = getDistance(pdb[i].ca, pdb[i].sc);
        //calculate phi
        if (i == 0) {
            poseData.phi = getDihedral(pdb[i+1].ca, N, pdb[i].ca, pdb[i].sc);
        }
        else if (i == pdb.size() - 1) {
            poseData.phi = getDihedral(C, pdb[i-1].ca, pdb[i].ca, pdb[i].sc);
        }
        else {
            poseData.phi = getDihedral(pdb[i+1].ca, pdb[i-1].ca, pdb[i].ca, pdb[i].sc);
        }
        //calculate delta
        if (i == 0) {
            poseData.delta = getAngle(N, pdb[i].ca, pdb[i].sc);
        }
        else {
            poseData.delta = getAngle(pdb[i-1].ca,  pdb[i].ca, pdb[i].sc);
        }
        // populate the pose
        pose.push_back(poseData);
    }
    
}

/*************************************************************************
 * Name        : pose2pdb
 * Purpose     : converts pose to pdb
 * Arguments   : vector<poseInfo> &pose, vector<pdbInfo> &pdb
 * Return Type : void
 *************************************************************************/
void pose2pdb(vector<poseInfo> &pose, vector<pdbInfo> &pdb) {
    pdb.clear();
    pdbInfo pdbData;
    //define a virtual N terminus
    point3d N;
    N.x = 0.0 - cos(PI) * CA2CA;
    N.y = 0.0 - sin(PI) * CA2CA;
    N.z = 0.0;
    // set the first residue's CA atom at the origin
    pdbData.id = pose[0].id;
    pdbData.aa = pose[0].aa;
    pdbData.ca.x = 0.0;
    pdbData.ca.y = 0.0;
    pdbData.ca.y = 0.0;
    pdbData.sc = pdbData.ca;
    pdb.push_back(pdbData);
    // set the second residue's CA atom by bond angle and bond length
    pdbData.id = pose[1].id;
    pdbData.aa = pose[1].aa;
    pdbData.ca.x = pose[1].bca;
    pdbData.ca.y = 0.0;
    pdbData.ca.y = 0.0;
    pdbData.sc = pdbData.ca;
    pdb.push_back(pdbData);
    // set the third residue's CA atom by dihedral, bond angle and bond length
    pdbData.id = pose[2].id;
    pdbData.aa = pose[2].aa;
    setCoordinate(N, pdb[0].ca, pdb[1].ca, pdbData.ca, pose[1].tao, pose[1].theta, pose[1].bca);
    pdbData.sc = pdbData.ca;
    pdb.push_back(pdbData);
    // set the rest of the residues' CA atoms by dihedral, bond angle and bond length
    for (int i = 3; i < pose.size(); i++) {
        pdbData.id = pose[i].id;
        pdbData.aa = pose[i].aa;
        setCoordinate(pdb[i-3].ca, pdb[i-2].ca, pdb[i-1].ca, pdbData.ca, pose[i-2].tao, pose[i-1].theta, pose[i-1].bca);
        pdbData.sc = pdbData.ca;
        pdb.push_back(pdbData);
    }
    // define a virtual C terminus
    point3d C;
    C.x = pdb[pdb.size()-1].ca.x + cos(PI) * CA2CA;
    C.y = pdb[pdb.size()-1].ca.y + sin(PI) * CA2CA;
    C.z = pdb[pdb.size()-1].ca.z;
    // set all residues' SC atoms based on CA atoms by dihedral, bond angle and bond length
    for (int i = 0; i < pose.size(); i++) {
        if (pose[i].aa != GLY) {
            if (i == 0) {
                setCoordinate(pdb[i+1].ca, N, pdb[i].ca, pdb[i].sc, pose[i].phi, pose[i].delta, pose[i].bsc);
            }
            else if (i == pose.size() - 1) {
                setCoordinate(C, pdb[i-1].ca, pdb[i].ca, pdb[i].sc, pose[i].phi, pose[i].delta, pose[i].bsc);
            }
            else {
                setCoordinate(pdb[i+1].ca, pdb[i-1].ca, pdb[i].ca, pdb[i].sc, pose[i].phi, pose[i].delta, pose[i].bsc);
            }
        }
    }
}

/*************************************************************************
 * Name        : sample2pose
 * Purpose     : converts sample to pose
 * Arguments   : MDArray<double> &sample, vector<poseInfo> &pose
 * Return Type : void
 *************************************************************************/
void sample2pose(MDArray<double> &sample, vector<poseInfo> &pose) {
    pose.clear();
    for (int i = 0; i < sample.get_shape()[0] / 2; i++) {
        poseInfo poseData;
        poseData.id = i + 1;
        poseData.aa = (int)sample.get(2 * i, 2);
        poseData.ss = (int)sample.get(2 * i, 3);
        poseData.bca = sample.get(2 * i, 4);
        poseData.tao = sample.get(2 * i, 5);
        poseData.theta = sample.get(2 * i, 6);
        poseData.bsc = sample.get(2 * i + 1, 4);
        poseData.phi = sample.get(2 * i + 1, 5);
        poseData.delta = sample.get(2 * i + 1, 6);
        pose.push_back(poseData);
    }
}

/*************************************************************************
 * Name        : pose2sample
 * Purpose     : converts pose to sample
 * Arguments   : vector<poseInfo> &pose, MDArray<double> &sample
 * Return Type : void
 *************************************************************************/
void pose2sample(vector<poseInfo> &pose, MDArray<double> &sample) {
    sample.set_shape(pose.size() * 2, NUM_DAT);
    for (int i = 0; i < pose.size(); i++) {
        sample.set(2 * i, 0, 0);
        sample.set(2 * i, 1, 0);
        sample.set(2 * i, 2, pose[i].aa);
        sample.set(2 * i, 3, pose[i].ss);
        sample.set(2 * i, 4, pose[i].bca);
        sample.set(2 * i, 5, pose[i].tao);
        sample.set(2 * i, 5, pose[i].theta);

        sample.set(2 * i + 1, 0, 1);
        sample.set(2 * i + 1, 1, 0);
        sample.set(2 * i + 1, 2, pose[i].aa);
        sample.set(2 * i + 1, 3, pose[i].ss);
        sample.set(2 * i + 1, 4, pose[i].bsc);
        sample.set(2 * i + 1, 5, pose[i].phi);
        sample.set(2 * i + 1, 5, pose[i].delta);
   }
}

/*************************************************************************
 * Name        : isHydrophobic
 * Purpose     : returns true if residue is hydrophobic, false otherwise
 * Arguments   : int aa
 * Return Type : bool
 *************************************************************************/
bool isHydrophobic(int aa) {
    if ( aa == CYS || aa == MET || aa == PHE || aa == ILE || aa == LEU || aa == VAL || aa == TRP || aa == TYR || aa == ALA) {
        return true;
    }
    else {
        return false;
    }
}

/*************************************************************************
 * Name        : getRrContact
 * Purpose     : gets the residue-residue contacts upto a specific residue
 * Arguments   : MDArray<double> &cm, int end
 * Return Type : void
 *************************************************************************/
void getRrContact(MDArray<double> &cm, int end) {
    rrContact.clear();
    for (int i = 0; i < end; i++) {
        for (int j = i + RR_SEQ_SEP; j < end; j++) {
            contactInfo contactData;
            contactData.ri = i;
            contactData.rj = j;
            contactData.wt = cm.get(i, j);
            rrContact.push_back(contactData);
        }
    }
}

/*************************************************************************
 * Name        : getUpperBoundHarmonic
 * Purpose     : gets upper bound harmonic potential
 * Arguments   : double d, double bound
 * Return Type : double
 *************************************************************************/
double getUpperBoundHarmonic(double d, double bound) {
    double potential = 0.0;
    if (d > bound) {
        potential = ((d - bound) * (d - bound));
    }
    return potential;
}

/*************************************************************************
 * Name        : getFade
 * Purpose     : gets FADE potential
 * Arguments   : double d, double lb, double ub, double z, double w
 * Return Type : double
 *************************************************************************/
double getFade(double d, double lb, double ub, double z, double w) {
    double potential = 0.0;
    double lf = lb + z;
    double uf = ub - z;
    if (d < lb || d > ub) {
        potential = 0.0;
    }
    else if (d < lf) {
        potential = -2.0 * pow(((d - lf) / z), 3.0) - 3.0 * pow(((d - lf) / z), 2.0);
    }
    else if (d < uf) {
        potential = -2.0 * pow(((d - uf) / z), 3.0) - 3.0 * pow(((d - uf) / z), 2.0);
    }
    else {
        potential = w;
    }
    return potential;
}

/*************************************************************************
 * Name        : getLiwoEpsilon0
 * Purpose     : gets the epsilon 0 values for amino acid pairs
 * Arguments   : int ai, int aj
 * Return Type : double
 *************************************************************************/
double getLiwoEpsilon0(int ai, int aj) {
    int id1 = liwoAminoStr.find(seq[ai]);
    int id2 = liwoAminoStr.find(seq[aj]);
    if (id1 < 0 || id1 > 19) {
        return 0;
    }
    if (id2 < 0 || id2 > 19) {
        return 0;
    }
    // swap in case id1 is less than id2
    if (id1 < id2)
    {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
    }
    return liwoEpsilon0[id1][id2];
}

/*************************************************************************
 * Name        : getScScEnergy
 * Purpose     : gets the side chain - side chain interaction energy (GBV)
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getScScEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < pdb.size() - 1; i++) {
        for (int j = i + 1; j < pdb.size(); j++) {
            // consider only secondary structural elements
            if (pose[i].ss != RANDOM_COIL && pose[j].ss != RANDOM_COIL) {
            //if (pose[i].ss != RANDOM_COIL && pose[j].ss != RANDOM_COIL) {
                double rij = getDistance(pdb[i].sc, pdb[j].sc);
                double rij0 = liwoR0[pose[i].aa] + liwoR0[pose[j].aa];
                double sigmaij0 = sqrt(liwoSigma0[pdb[i].aa] * liwoSigma0[pdb[i].aa] + liwoSigma0[pdb[j].aa] * liwoSigma0[pdb[j].aa]);
                point3d dij1 = getDifference(pdb[i].ca, pdb[i].sc);
                point3d dij2 = getDifference(pdb[j].ca, pdb[j].sc);
                point3d drij = getDifference(pdb[i].sc, pdb[j].sc);
                point3d uij1 = getUnit(dij1);
                point3d uij2 = getUnit(dij2);
                point3d urij2 = getUnit(drij);
                double omegaij1 = getDotProduct(uij1, urij2);
                double omegaij2 = getDotProduct(uij2, urij2);
                double omegaij12 = getDotProduct(uij1, uij2);
                double chiij1Nr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[i].aa] - 1.0;
                double chiij1Dr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[i].aa] + pow((liwoSigma0[pdb[j].aa]/liwoSigma0[pdb[i].aa]), 2.0);
                double chiij1 = chiij1Nr / chiij1Dr;
                double chiij2Nr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[j].aa] - 1.0;
                double chiij2Dr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[j].aa] + pow((liwoSigma0[pdb[i].aa]/liwoSigma0[pdb[j].aa]), 2.0);
                double chiij2 = chiij2Nr / chiij2Dr;
                double epsilonij0 = getLiwoEpsilon0(pdb[i].aa, pdb[j].aa);
                double epsilonij1Dr = (1.0 - chiij1 * chiij2 * omegaij12);
                double epsilonij1 = sqrt(epsilonij1Dr);
                double epsilonij2Nr = (liwoChiPrime[pdb[i].aa] * omegaij1 * omegaij1 + liwoChiPrime[pdb[j].aa] * omegaij2 * omegaij2 - 2.0 * liwoChiPrime[pdb[i].aa] * liwoChiPrime[pdb[j].aa] * omegaij1 * omegaij2 * omegaij12);
                double epsilonij2Dr = (1.0 - liwoChiPrime[pose[i].aa] * liwoChiPrime[pose[j].aa] * omegaij12 * omegaij12);
                double epsilonij2 = (1.0 - epsilonij2Nr / epsilonij2Dr) * (1.0 - epsilonij2Nr / epsilonij2Dr);
                double epsilonij3Nr = (1.0 - liwoAlpha[pdb[i].aa] * omegaij1 + liwoAlpha[pdb[j].aa] * omegaij2 - 0.5 * (liwoAlpha[pdb[i].aa] + liwoAlpha[pdb[j].aa]) * omegaij12);
                double epsilonij3 = epsilonij3Nr * epsilonij3Nr;
                double sigmaijNr = (chiij1 * omegaij1 * omegaij1 + chiij2 * omegaij2 * omegaij2 - 2.0 * chiij1 * chiij2 * omegaij1 * omegaij2 * omegaij12);
                double sigmaijDr = (1.0 - chiij1 * chiij2 * omegaij12 * omegaij12);
                double sigmaij = sigmaij0 * sqrt(1.0 - sigmaijNr / sigmaijDr);;
                double epsilonij = epsilonij0 * epsilonij1 * epsilonij2 * epsilonij3;
                double xijNr = rij0;
                double xijDr = (rij - sigmaij + rij0);
                double xij = 0.0;
                if (xijDr != 0.0) {
                    xij = xijNr / xijDr;
                }
                energy += 4 * (abs(epsilonij) * pow(xij, 12.0) - epsilonij * pow(xij, 6.0));
            }
        }
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getScBbEnergy
 * Purpose     : gets the side chain - backbone interaction energy (repulsive only)
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getScBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < pdb.size(); i++) {
        for (int j = 0; j < pdb.size() - 1; j++) {
            if (j != i - 1 && j != i) {
                point3d pj = getMidpoint(pdb[j].ca, pdb[j+1].ca);
                double rij = getDistance(pdb[i].sc, pj);
                energy += 0.3 * pow((4.0 / rij), 6.0);
            }
        }
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getBbBbEnergy
 * Purpose     : gets the backnone - backbone interaction energy (electrostatic)
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getBbBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < pdb.size() - 3; i++) {
        for (int j = i + 2; j < pdb.size() - 1; j++) {
            double A = 0.0;
            double B = 0.0;
            double epsilon = 0.0;
            double r0 = 0.0;
            // proline-proline interaction
            if (pdb[i+1].aa == PRO && pdb[j+1].aa == PRO) {
                epsilon = 0.574;
                r0 = 4.48;
                A = 5.13;
                B = 335.0;
            }
            // ordinary-proline interaction
            if ((pdb[i+1].aa == PRO && pdb[j+1].aa != PRO) || (pdb[i+1].aa != PRO && pdb[j+1].aa == PRO)) {
                epsilon = 0.365;
                r0 = 4.54;
                A = 0.0;
                B = 1129.0;
            }
            // ordinary-ordinary interaction
            if (pdb[i+1].aa != PRO && pdb[j+1].aa != PRO) {
                epsilon = 0.305;
                r0 = 4.51;
                A = 3.73;
                B = 1306.0;
            }
            point3d vi = getDifference(pdb[i].ca, pdb[i+1].ca);
            point3d vj = getDifference(pdb[j].ca, pdb[j+1].ca);
            point3d uvi = getUnit(vi);
            point3d uvj = getUnit(vj);
            point3d pi = getMidpoint(pdb[i].ca, pdb[i+1].ca);
            point3d pj = getMidpoint(pdb[j].ca, pdb[j+1].ca);
            point3d pij = getDifference(pi, pj);
            point3d upij = getUnit(pij);
            double rij = getDistance(pi, pj);
            double cosAlphaij = getDotProduct(uvi, uvj);
            double cosBetaij = getDotProduct(uvi, upij);
            double cosGammaij = getDotProduct(uvj, upij);
            double term1 = (A / pow(rij, 3.0)) * (cosAlphaij - 3.0 * cosBetaij * cosGammaij);
            double term2 = (B / pow(rij, 6.0)) * (4.0 + ((cosAlphaij - 3.0 * cosBetaij * cosGammaij) * (cosAlphaij - 3.0 * cosBetaij * cosGammaij)) - 3.0 * (cosBetaij * cosBetaij + cosGammaij * cosGammaij));
            double term3 = epsilon * (pow((r0 / rij), 12.0) - 2.0 * pow((r0 / rij), 6.0));
            double total = term1 - term2 + term3;
            energy += total;
        }
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getRiRjEnergy
 * Purpose     : gets the residue-residue contact energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getRiRjEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < rrContact.size(); i++) {
        contactInfo contactData;
        int ci = rrContact[i].ri;
        int cj = rrContact[i].rj;
        double dij = getDistance(pdb[ci].sc, pdb[cj].sc);
        double pij = rrContact[i].wt;
        if (dij <= RR_CONTACT_THRESHOLD ) {
            energy += -pij;
        }
        else {
            energy += (-pij * exp(-pow((dij - RR_CONTACT_THRESHOLD), 2.0))) + (pij * ((dij - RR_CONTACT_THRESHOLD)/dij));
        }
        //energy += (log(abs(ci - cj)) * pij * (dij - RR_CONTACT_THRESHOLD));
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getWeightedEnergy
 * Purpose     : gets total weighted energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn
 * Return Type : energyInfo
 *************************************************************************/
energyInfo getWeightedEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn) {
    energyInfo e;
    e.sc_sc = getScScEnergy(pdb, pose, 1.0);
    e.sc_bb = getScBbEnergy(pdb, pose, 1.0);
    e.bb_bb = getBbBbEnergy(pdb, pose, 1.0);
    e.ri_rj = getRiRjEnergy(pdb, pose, 1.0);
    e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj;
    return e;
}

/*************************************************************************
 * Name        : showEnergy
 * Purpose     : shows (prints) the energy breakdown and total weighted energy
 * Arguments   : energyInfo
 * Return Type : void
 *************************************************************************/
void showEnergy(energyInfo &e) {
    char buf[300];
    sprintf(buf,"UniCon3D.scoring: E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_total = %8.3f", e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total);
    printf("%s\r",buf);
    fflush(stdout);
    //cout << buf << endl;
}

/*************************************************************************
 * Name        : getRmsd
 * Purpose     : gets the root mean square deviation between two pdb objects
 * Arguments   : vector<pdbInfo> &pdb1, vector<pdbInfo> &pdb2, string mode
 * Return Type : double
 *************************************************************************/
double getRmsd(vector<pdbInfo> &pdb1, vector<pdbInfo> &pdb2, string mode) {
    if (pdb1.size() != pdb2.size()) {
        cout << "Error! The two compared structures are of different length" << endl;
        exit(0);
    }
    double rot[3][3];
    double trans[3];
    vector<vector<double> > A;
    for (int i = 0; i < pdb1.size(); i++) {
        vector<double> a(3);
        //load CA only
        if (mode.compare("ca") == 0) {
            a[0] = pdb1[i].ca.x;
            a[1] = pdb1[i].ca.y;
            a[2] = pdb1[i].ca.z;
            A.push_back(a);
        }
        //load SC only
        else if (mode.compare("sc") == 0) {
            a[0] = pdb1[i].sc.x;
            a[1] = pdb1[i].sc.y;
            a[2] = pdb1[i].sc.z;
            A.push_back(a);
        }
        //load all atoms (CA + SC)
        else if (mode.compare("aa") == 0) {
            a[0] = pdb1[i].ca.x;
            a[1] = pdb1[i].ca.y;
            a[2] = pdb1[i].ca.z;
            A.push_back(a);
            a[0] = pdb1[i].sc.x;
            a[1] = pdb1[i].sc.y;
            a[2] = pdb1[i].sc.z;
            A.push_back(a);
        }
        else {
            cout << "Error! Invalid mode in rmsd calculation" << endl;
            exit(0);
        }
    }
    vector<vector<double> > B;
    for (int i = 0; i < pdb2.size(); i++) {
        vector<double> b(3);
        //load CA only
        if (mode.compare("ca") == 0) {
            b[0] = pdb2[i].ca.x;
            b[1] = pdb2[i].ca.y;
            b[2] = pdb2[i].ca.z;
            B.push_back(b);
        }
        //load SC only
        else if (mode.compare("sc") == 0) {
            b[0] = pdb2[i].sc.x;
            b[1] = pdb2[i].sc.y;
            b[2] = pdb2[i].sc.z;
            B.push_back(b);
        }
        //load all atoms (CA + SC)
        else if (mode.compare("aa") == 0) {
            b[0] = pdb2[i].ca.x;
            b[1] = pdb2[i].ca.y;
            b[2] = pdb2[i].ca.z;
            B.push_back(b);
            b[0] = pdb2[i].sc.x;
            b[1] = pdb2[i].sc.y;
            b[2] = pdb2[i].sc.z;
            B.push_back(b);
        }
        else {
            cout << "Error! Invalid mode in rmsd calculation" << endl;
            exit(0);
        }
    }
    double Bt[A.size()][3];
    getAlignment(A, B, rot, trans);
    // Transform the second chain to optimally align with the first.
    for (int k = 0; k < A.size(); k++) {
        Bt[k][X] = B[k][X] * rot[0][0] + B[k][Y] * rot[1][0] +
        B[k][Z] * rot[2][0] + trans[0];
        Bt[k][Y] = B[k][X] * rot[0][1] + B[k][Y] * rot[1][1] +
        B[k][Z] * rot[2][1] + trans[1];
        Bt[k][Z] = B[k][X] * rot[0][2] + B[k][Y] * rot[1][2] +
        B[k][Z] * rot[2][2] + trans[2];
    }
    double rmsd = 0;
    for (int i = 0; i < A.size(); i++) {
        double a0 = Bt[i][X] - A[i][X];
        double a1 = Bt[i][Y] - A[i][Y];
        double a2 = Bt[i][Z] - A[i][Z];
        
        rmsd += (a0*a0 + a1*a1 + a2*a2);
    }
    return sqrt(rmsd / A.size());
}

/*************************************************************************
 * Name        : getAlignment
 * Purpose     : gets the optimal alignemnt between two sets of point clouds
 * Arguments   : vector<vector<double> > &  A, vector<vector<double> > & B, double rot[3][3], double trans[3]
 * Return Type : void
 *************************************************************************/
void getAlignment(vector<vector<double> > &  A, vector<vector<double> > & B, double rot[3][3], double trans[3]) {
    int i,j;
    double c1[3],c2[3];   /* center of mass for two point collections */
    double v1[3],v2[3];
    double recip;
    double tr;
    double m[4][4], q[4][4];
    double v[4];
    double cov[3][3];
    double aij[3][3];
    double quat[4];
    // find the center of mass for the two collections
    c1[X] = c1[Y] = c1[Z] = 0;
    c2[X] = c2[Y] = c2[Z] = 0;
    
    for (int i = 0; i < A.size(); i++) {
        c1[X] += A[i][X];
        c1[Y] += A[i][Y];
        c1[Z] += A[i][Z];
        
        c2[X] += B[i][X];
        c2[Y] += B[i][Y];
        c2[Z] += B[i][Z];
    }
    recip = 1.0 / A.size();
    c1[X] *= recip;
    c1[Y] *= recip;
    c1[Z] *= recip;
    c2[X] *= recip;
    c2[Y] *= recip;
    c2[Z] *= recip;
    // create the cross-covariance matrix
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            cov[i][j] = 0;
    for (int i = 0; i < A.size(); i++) {
        v1[X] = A[i][X] - c1[X];
        v1[Y] = A[i][Y] - c1[Y];
        v1[Z] = A[i][Z] - c1[Z];
        v2[X] = B[i][X] - c2[X];
        v2[Y] = B[i][Y] - c2[Y];
        v2[Z] = B[i][Z] - c2[Z];
        cov[X][X] += v1[X] * v2[X];
        cov[X][Y] += v1[X] * v2[Y];
        cov[X][Z] += v1[X] * v2[Z];
        cov[Y][X] += v1[Y] * v2[X];
        cov[Y][Y] += v1[Y] * v2[Y];
        cov[Y][Z] += v1[Y] * v2[Z];
        cov[Z][X] += v1[Z] * v2[X];
        cov[Z][Y] += v1[Z] * v2[Y];
        cov[Z][Z] += v1[Z] * v2[Z];
    }
    // aij = cov - transpose(cov)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            aij[i][j] = cov[i][j] - cov[j][i];
    // find the trace of the covariance matrix
    tr = cov[X][X] + cov[Y][Y] + cov[Z][Z];
    m[0][0] = tr;
    m[1][0] = m[0][1] = aij[1][2];
    m[2][0] = m[0][2] = aij[2][0];
    m[3][0] = m[0][3] = aij[0][1];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            m[i+1][j+1] = cov[i][j] + cov[j][i] - (i == j) * tr;
    // find the eigenvector corresponding to the largest eigenvalue of this matrix
    Meigen4(q, v, m);
    if( v[0] > v[1] ) {
        if( v[0] > v[2] ) {
            if( v[0] > v[3] )
                i = 0;
            else
                i = 3;
        }
        else{
            if( v[2] > v[3] )
                i = 2;
            else
                i =3;
        }
    }
    else {
        if( v[1] > v[2] ){
            if( v[1] > v[3] )
                i = 1;
            else
                i = 3;
        }
        else {
            if( v[2] > v[3] )
                i = 2;
            else
                i =3;
        }
    }
    quat[0] = q[0][i];
    quat[1] = q[1][i];
    quat[2] = q[2][i];
    quat[3] = q[3][i];
    quat2mat(quat, rot);
    double c3[3];
    // determine best translation
    vApply (rot, c2, c3);
    trans[0] = c1[0] - c3[0];
    trans[1] = c1[1] - c3[1];
    trans[2] = c1[2] - c3[2];
}

/*************************************************************************
 * Name        : quat2mat
 * Purpose     : convert quat to mat
 * Arguments   : double q[4], double mat[3][3]
 * Return Type : void
 *************************************************************************/
void quat2mat(double q[4], double mat[3][3]) {
    double q00,q01,q02,q03;
    double q11,q12,q13;
    double q22,q23;
    double q33;
    q00 = q[0] * q[0];
    q01 = q[0] * q[1];
    q02 = q[0] * q[2];
    q03 = q[0] * q[3];
    q11 = q[1] * q[1];
    q12 = q[1] * q[2];
    q13 = q[1] * q[3];
    q22 = q[2] * q[2];
    q23 = q[2] * q[3];
    q33 = q[3] * q[3];
    mat[X][X] = q00 + q11 - q22 - q33;
    mat[X][Y] = 2 * (q12 - q03);
    mat[X][Z] = 2 * (q13 + q02);
    mat[Y][X] =  2 * (q12 + q03);
    mat[Y][Y] =  q00 + q22 - q11 - q33;
    mat[Y][Z] = 2 * (q23 - q01);
    mat[Z][X] = 2 * (q13 - q02);
    mat[Z][Y] = 2 * (q23 + q01);
    mat[Z][Z] = q00 + q33 - q11 - q22;
}

/*************************************************************************
 * Name        : vApply
 * Purpose     : apply translation
 * Arguments   : double m[3][3], const double a[3], double b[3]
 * Return Type : void
 *************************************************************************/
void vApply (double m[3][3], const double a[3], double b[3]) {
    int j;
    for (j = 0; j <= 2; j++)
        b[j] = a[0] * m[0][j] + a[1] * m [1][j] +
        a[2] * m[2][j];
}

/*************************************************************************
 * Name        : assembleFoldon
 * Purpose     : assemble foldon unit via simulated anneling
 * Arguments   : vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow
 * Return Type : void
 *************************************************************************/
void assembleFoldon(vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow) {
    int start = 0;
    if (i != 0) {
        start = foldonSample[i-1].get_shape()[0];
    }
    int end = foldonSample[i].get_shape()[0];
    // populate the current foldonSample based on the previos other than the first foldonSample
    if (i != 0) {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0]; j++) {
            foldonSample[i].set(j, 0, foldonSample[i-1].get(j, 0));
            foldonSample[i].set(j, 1, foldonSample[i-1].get(j, 1));
            foldonSample[i].set(j, 2, foldonSample[i-1].get(j, 2));
            foldonSample[i].set(j, 3, foldonSample[i-1].get(j, 3));
            foldonSample[i].set(j, 4, foldonSample[i-1].get(j, 4));
            foldonSample[i].set(j, 5, foldonSample[i-1].get(j, 5));
            foldonSample[i].set(j, 6, foldonSample[i-1].get(j, 6));
        }
    }
    // populate data and mismask from foldonSample
    Sequence data;
    data.set_shape(foldonSample[i].get_shape()[0], NUM_DAT);
    MDArray<eMISMASK> mism;
    mism.set_shape(foldonSample[i].get_shape()[0], NUM_MIS);
    if (i == 0) {
        for (int j = 0; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            mism.set(j * 2, 4, MOCAPY_HIDDEN);
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            mism.set(j * 2, 5, MOCAPY_HIDDEN);
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
            }
        }
    }
    else {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
                mism.set(j * 2, 4, MOCAPY_HIDDEN);
            }
            else {
                mism.set(j * 2, 4, MOCAPY_OBSERVED);
            }
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
                mism.set(j * 2, 5, MOCAPY_HIDDEN);
            }
            else {
                mism.set(j * 2, 5, MOCAPY_OBSERVED);
            }
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
            }
        }
        for (int j = foldonSample[i-1].get_shape()[0] / 2; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            mism.set(j * 2, 4, MOCAPY_HIDDEN);
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            mism.set(j * 2, 5, MOCAPY_HIDDEN);
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
            }
        }
    }
    // setup the sampler
    SampleInfEngineHMM sampler(&dbn, data, mism, 1);
    //MDArray<double> sample = sampler.sample_next();
    //sampler.undo();
    sampler.set_start_end(start, end);
    MDArray<double> sample = sampler.sample_next();

    double Ti = START_TEMP;
    double Tf = FINAL_TEMP;
    
    // perform simulated anneling
    doSimulatedAnneling(sampler, dbn, data, mism, sample, poseLow, MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf);
    // populate the lowest scoring sample into current foldonSample
    for (int j = 0; j < sample.get_shape()[0]; j++) {
        foldonSample[i].set(j, 0, sample.get(j, 0));
        foldonSample[i].set(j, 1, sample.get(j, 1));
        foldonSample[i].set(j, 2, sample.get(j, 2));
        foldonSample[i].set(j, 3, sample.get(j, 3));
        foldonSample[i].set(j, 4, sample.get(j, 4));
        foldonSample[i].set(j, 5, sample.get(j, 5));
        foldonSample[i].set(j, 6, sample.get(j, 6));
    }
    pose2pdb(poseLow, pdbLow);
}

/*************************************************************************
 * Name        : doSimulatedAnneling
 * Purpose     : perform simulated anneling
 * Arguments   : SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp, CImgDisplay &dispLow, CImgDisplay &dispLive
 * Return Type : void
 *************************************************************************/
void doSimulatedAnneling(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp) {
    vector<poseInfo> pose;
    sample2pose(sample, pose);
    vector<pdbInfo> pdb;
    pose2pdb(pose, pdb);
    energyInfo energy = getWeightedEnergy(pdb, pose, dbn);
    
    // assign lowest energy pose to current pose
    poseLow = pose;
    vector<pdbInfo> pdbLow;
    pdbLow = pdb;
    energyInfo energyLow = getWeightedEnergy(pdbLow, poseLow, dbn);
    
    // define integer, discrete uniform distribution for position
    boost::uniform_int<> uniIntPos(0, sample.get_shape()[0]/2 - maxLength);
    // define a random variate generator using our base generator and distribution
    boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenPos(baseGen, uniIntPos);
    
    // define integer, discrete uniform distribution for length
    boost::uniform_int<> uniIntLen(minLength, maxLength);
    // define a random variate generator using our base generator and distribution
    boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenLen(baseGen, uniIntLen);
    
    // cycles and temperature schedule
    int outerCycles = sample.get_shape()[0]/2;
    int innerCycles = numCycles;
    //double initTemp = 1000.0;
    //double finalTemp = 298.0;
    double gamma = pow((finalTemp / initTemp),(double)(1.0 / (outerCycles * innerCycles)));
    double T = initTemp;
    
    // simulated anneling
    for (int i = 0; i < outerCycles; i++) {
        //sampler.set_seq_mismask(sample, mism);
        for (int j = 0; j < innerCycles; j++) {
            
            // randomly sample start and end
            int startPos = uniIntGenPos();
            int endPos;
            int len;
            while(true) {
                len = uniIntGenLen();
                endPos = startPos + len;
                if (endPos < sample.get_shape()[0]/2) {
                    break;
                }
            }
            sampler.set_start_end(startPos * 2, endPos * 2);
            MDArray<double> sampleNext = sampler.sample_next();
            
            for (int k = startPos; k < endPos; k++) {
                if (mism.get(k * 2, 5) == MOCAPY_OBSERVED) {
                	sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
                  sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
                }
            }
            vector<poseInfo> poseNext;
            sample2pose(sampleNext, poseNext);
            vector<pdbInfo> pdbNext;
            pose2pdb(poseNext, pdbNext);
            energyInfo energyNext = getWeightedEnergy(pdbNext, poseNext, dbn);
            // accept the move
            if (acceptMove(energy, energyNext, T)) {
                sampler.set_seq_mismask(sampleNext, mism);
                energy = energyNext;

                //if (((i + 1) * ( j + 1)) % innerCycles == 0) {
                    showEnergy(energyNext);
                //}
                
            }
            // reject the move
            else {
                sampler.undo();
            }
            // in case energy is lower than lowest-energy conformation so far
            if (energyNext.total < energyLow.total) {
                poseLow = poseNext;
                pdbLow = pdbNext;
                energyLow = energyNext;
                sample = sampleNext;
            }
            T = T * gamma;
        }
        
    }
}

/*************************************************************************
 * Name        : acceptMove
 * Purpose     : accept or reject a move
 * Arguments   : energyInfo &prev, energyInfo &curr, double temperature
 * Return Type : bool
 *************************************************************************/
bool acceptMove(energyInfo &prev, energyInfo &curr, double temperature) {
    double diff = prev.total - curr.total;
    if (diff > 0)
        return true;
    else {
        double e = exp(diff / (temperature * BOLTZMANN_CONSTANT));
        double r = uniDblGen();
        if (e > r)
            return true;
        else
            return false;
    }
}

/*************************************************************************
 * Name        : writePdb
 * Purpose     : writes a pdb file in  united residue representation into a file
 * Arguments   : vector<pdbInfo> &pdb, char * filename
 * Return Type : void
 *************************************************************************/
void writePdb(vector<pdbInfo> &pdb, char * filename) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cout << "Error! output pdb file can not open " << filename << endl;
        return;
    }
    
    for (int i=0; i < pdb.size(); i++) {
        // buffer CA atoms
        char bufCa[300];
        sprintf(bufCa,"ATOM %6d  CA  %-3s %5d    %8.3f%8.3f%8.3f  1.00  0.00", 2 * i + 1, seq3[pdb[i].aa].c_str(), pdb[i].id, pdb[i].ca.x, pdb[i].ca.y, pdb[i].ca.z);
        // buffer SC atoms
        char bufSc[300];
        sprintf(bufSc,"ATOM %6d  SC  %-3s %5d    %8.3f%8.3f%8.3f  1.00  0.00", 2 * i + 2, seq3[pdb[i].aa].c_str(), pdb[i].id, pdb[i].sc.x, pdb[i].sc.y, pdb[i].sc.z);
        // write CA and SC atoms to file
        fout << bufCa << endl;
        fout << bufSc << endl;
    }
}
