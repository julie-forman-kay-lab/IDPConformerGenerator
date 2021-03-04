/*
This module contains a copy of the FASPR algorithm for sidechain packing
taken from:

https://github.com/tommyhuangthu/FASPR/tree/f0e6a6d8e8312f34341203f2600cf18df252cab1

Please read the instructions in the original repository for references
to the authors and how to cite the original work. Please, cite FASPR
original articles if you use this module!

The original code is licensed under MIT License and a copy of that
license is provided in the same folder as this file.

For simplicity, the original code was merged into this single file.
Hence, the position of functions and classes in this file does not
represent the authors original layout, instead, this file was structured
by @joaomcteixeira.

The original functions and methods related to data input and output,
which are not needed here, were not included or are included but
commented.

Handling of dun2010bbdeb.bin file has been adapted as well. The bin file
is read once and kept as a global variable, and its state is managed by
faspr_sidechains(). Hence the implementation here and does not reflect
the strategy the FASPR authors adopted to handle this file.

One2Three was rewritten.

Message outputs were disable.

All other parts of the FASPR software are maintained EXACTLY as in the
original code, in both execution and presentation.

The pybind11 adapter and the other adapter functions needed to transfer
data from the Python modules to this module were created by
@joaomcteixeira. These are signed in the docstring accordingly.

pybind11 is available at:
https://github.com/pybind/pybind11

@joaomcteixeira thanks @tdegeus for his tutorials at:
https://github.com/tdegeus/pybind11_examples
*/

//#pragma warning(disable:4018)
//#pragma warning(disable:4305)
//#pragma warning(disable:4101)
//#pragma warning(disable:4244)
//
#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

using namespace std;
namespace py = pybind11;

ifstream infile;

//# contants
#define DEE_THRESHOLD 0.0
#define EPAIR_CUT      2.0
#define TREEWIDTH_CUT  5

#define WGT_HBOND       1.0
#define WGT_SSBOND      6.0
#define CACB_DIST_CUT   2.35
#define RESI_DIST_CUT   4.25
#define SEC_CUT        15.0
//atomic parameters
#define RADIUS_C1       1.78
#define RADIUS_C2       1.40
#define RADIUS_C3       2.30
#define RADIUS_C4       1.91
#define RADIUS_C5       1.96
#define RADIUS_C6       1.73
#define RADIUS_C7       1.43
#define RADIUS_C8       1.99
#define RADIUS_O1       1.48
#define RADIUS_O2       1.44
#define RADIUS_O3       1.40
#define RADIUS_O4       1.43
#define RADIUS_N1       1.42
#define RADIUS_N2       1.69
#define RADIUS_N3       1.56
#define RADIUS_N4       1.70
#define RADIUS_S1       2.15
#define RADIUS_S2       1.74
#define DEPTH_C1        0.25
#define DEPTH_C2        0.14
#define DEPTH_C3        0.30
#define DEPTH_C4        0.37
#define DEPTH_C5        0.48
#define DEPTH_C6        0.38
#define DEPTH_C7        0.07
#define DEPTH_C8        0.38
#define DEPTH_O1        0.22
#define DEPTH_O2        0.27
#define DEPTH_O3        0.07
#define DEPTH_O4        0.12
#define DEPTH_N1        0.08
#define DEPTH_N2        0.24
#define DEPTH_N3        0.46
#define DEPTH_N4        0.48
#define DEPTH_S1        0.44
#define DEPTH_S2        0.40
//VDW
#define DSTAR_MIN_CUT   0.015
#define DSTAR_MAX_CUT   1.90
#define VDW_REP_CUT    10.0
//Hbond energy
#define OPT_HBOND_DIST  2.8
#define MIN_HBOND_DIST  2.6
#define MAX_HBOND_DIST  3.2
#define MIN_HBOND_THETA 90.
#define MIN_HBOND_PHI   90.
//SSbond energy
#define OPT_SSBOND_DIST 2.03
#define MIN_SSBOND_DIST 1.73
#define MAX_SSBOND_DIST 2.53
#define OPT_SSBOND_ANGL 105.
#define MIN_SSBOND_ANGL 75.
#define MAX_SSBOND_ANGL 135.

// weights
#define WGT_CYS  5.5
#define WGT_ASP  2.0
#define WGT_GLU  1.0
#define WGT_PHE  1.5
#define WGT_HIS  3.0
#define WGT_ILE  1.0
#define WGT_LYS  2.0
#define WGT_LEU  2.0
#define WGT_MET  1.5
#define WGT_ASN  2.0
#define WGT_PRO  1.5
#define WGT_GLN  2.5
#define WGT_ARG  1.5
#define WGT_SER  1.5
#define WGT_THR  2.0
#define WGT_VAL  2.0
#define WGT_TRP  3.5
#define WGT_TYR  1.5
#define ROT_PROB_CUT_MIN  0.01
#define ROT_PROB_CUT_ACC  0.97

const float DELTA=1.0e-15;
const float BONDDIST=2.1;
const float PI=3.1415926;
const float DEG2RAD=PI/180.0e0;
const float RAD2DEG=180.0e0/PI;

//# mappings
char Three2One(string aa3)
{
  if(aa3=="ALA")return 'A';
  else if(aa3=="CYS")return 'C';
  else if(aa3=="ASP")return 'D';
  else if(aa3=="GLU")return 'E';
  else if(aa3=="PHE")return 'F';
  else if(aa3=="GLY")return 'G';
  else if(aa3=="HIS")return 'H';
  else if(aa3=="ILE")return 'I';
  else if(aa3=="LYS")return 'K';
  else if(aa3=="LEU")return 'L';
  else if(aa3=="MET")return 'M';
  else if(aa3=="ASN")return 'N';
  else if(aa3=="PRO")return 'P';
  else if(aa3=="GLN")return 'Q';
  else if(aa3=="ARG")return 'R';
  else if(aa3=="SER")return 'S';
  else if(aa3=="THR")return 'T';
  else if(aa3=="VAL")return 'V';
  else if(aa3=="TRP")return 'W';
  else if(aa3=="TYR")return 'Y';

  else if(aa3=="AYA")return 'A';
  else if(aa3=="CAS")return 'C';
  else if(aa3=="CAY")return 'C';
  else if(aa3=="CEA")return 'C';
  else if(aa3=="CME")return 'C';
  else if(aa3=="CMT")return 'C';
  else if(aa3=="CSB")return 'C';
  else if(aa3=="CSD")return 'C';
  else if(aa3=="CSE")return 'C';
  else if(aa3=="CSO")return 'C';
  else if(aa3=="CSP")return 'C';
  else if(aa3=="CSS")return 'C';
  else if(aa3=="CSW")return 'C';
  else if(aa3=="CSX")return 'C';
  else if(aa3=="CYG")return 'C';
  else if(aa3=="CYM")return 'C';
  else if(aa3=="NPH")return 'C';
  else if(aa3=="OCS")return 'C';
  else if(aa3=="OCY")return 'C';
  else if(aa3=="SNC")return 'C';
  else if(aa3=="PYX")return 'C';
  else if(aa3=="SMC")return 'C';
  else if(aa3=="SNN")return 'D';
  else if(aa3=="ASQ")return 'D';
  else if(aa3=="BHD")return 'D';
  else if(aa3=="DOH")return 'D';
  else if(aa3=="PHD")return 'D';
  else if(aa3=="CGU")return 'E';
  else if(aa3=="EHP")return 'F';
  else if(aa3=="ACY")return 'G';
  else if(aa3=="GL3")return 'G';
  else if(aa3=="HSD")return 'H';
  else if(aa3=="HSE")return 'H';
  else if(aa3=="HSP")return 'H';
  else if(aa3=="MHS")return 'H';
  else if(aa3=="NEP")return 'H';
  else if(aa3=="H2P")return 'H';
  else if(aa3=="HIC")return 'H';
  else if(aa3=="HIP")return 'H';
  else if(aa3=="INI")return 'K';
  else if(aa3=="MLY")return 'K';
  else if(aa3=="MLZ")return 'K';
  else if(aa3=="KCX")return 'K';
  else if(aa3=="LLP")return 'K';
  else if(aa3=="LLY")return 'K';
  else if(aa3=="LYZ")return 'K';
  else if(aa3=="M3L")return 'K';
  else if(aa3=="CXM")return 'M';
  else if(aa3=="FME")return 'M';
  else if(aa3=="MHO")return 'M';
  else if(aa3=="MSE")return 'M';
  else if(aa3=="OMT")return 'M';
  else if(aa3=="SME")return 'M';
  else if(aa3=="ASX")return 'N';
  else if(aa3=="MEN")return 'N';
  else if(aa3=="HYP")return 'P';
  else if(aa3=="PRS")return 'P';
  else if(aa3=="GLX")return 'Q';
  else if(aa3=="MGN")return 'Q';
  else if(aa3=="PCA")return 'Q';
  else if(aa3=="AAR")return 'R';
  else if(aa3=="AGM")return 'R';
  else if(aa3=="OPR")return 'R';
  else if(aa3=="MIS")return 'S';
  else if(aa3=="SEP")return 'S';
  else if(aa3=="SVA")return 'S';
  else if(aa3=="AEI")return 'T';
  else if(aa3=="TPO")return 'T';
  else if(aa3=="FTR")return 'W';
  else if(aa3=="HTR")return 'W';
  else if(aa3=="TRF")return 'W';
  else if(aa3=="TRN")return 'W';
  else if(aa3=="TRO")return 'W';
  else if(aa3=="ACE")return 'X';
  else if(aa3=="UNK")return 'X';
  else if(aa3=="FOR")return 'X';
  else if(aa3=="TPQ")return 'Y';
  else if(aa3=="TYI")return 'Y';
  else if(aa3=="TYN")return 'Y';
  else if(aa3=="TYQ")return 'Y';
  else if(aa3=="TYS")return 'Y';
  else if(aa3=="TYY")return 'Y';
  else if(aa3=="YOF")return 'Y';
  else if(aa3=="PAQ")return 'Y';
  else if(aa3=="PTH")return 'Y';
  return 'X';
}

// rewritten
string One2Three(char aa1)
{
  if(aa1 == 'A')return "ALA";
  else if(aa1 == 'C')return "CYS";
  else if(aa1 == 'D')return "ASP";
  else if(aa1 == 'E')return "GLU";
  else if(aa1 == 'F')return "PHE";
  else if(aa1 == 'G')return "GLY";
  else if(aa1 == 'H')return "HIS";
  else if(aa1 == 'I')return "ILE";
  else if(aa1 == 'K')return "LYS";
  else if(aa1 == 'L')return "LEU";
  else if(aa1 == 'M')return "MET";
  else if(aa1 == 'N')return "ASN";
  else if(aa1 == 'P')return "PRO";
  else if(aa1 == 'Q')return "GLN";
  else if(aa1 == 'R')return "ARG";
  else if(aa1 == 'S')return "SER";
  else if(aa1 == 'T')return "THR";
  else if(aa1 == 'V')return "VAL";
  else if(aa1 == 'W')return "TRP";
  else if(aa1 == 'Y')return "TYR";
  return "XXX";
  }

//# types
typedef map<int, set <int> > Graph;
typedef vector<float> FV1;
typedef vector<vector<float> > FV2;
typedef vector<vector<vector<float> > > FV3;
typedef vector<vector<vector<vector<float> > > > FV4;
typedef vector<double> DV1;
typedef vector<vector<double> > DV2;
typedef vector<vector<vector<double> > > DV3;
typedef vector<vector<vector<vector<double> > > > DV4;
typedef vector<int> IV1;
typedef vector<vector<int> > IV2;
typedef vector<vector<vector<int> > > IV3;
typedef vector<vector<vector<vector<int> > > > IV4;
typedef vector<string> SV1;
typedef vector<vector<string> > SV2;
typedef enum{
  Root,
  Inner,
  Leaf,
  None
}BagType;



//# Function prototypes
bool Internal2Cartesian(FV1 &c1,FV1 &c2,FV1 &c3,FV1 &p,FV1 &cc);
bool VectorL2Norm(FV1 &data);
bool VectorNormalization(FV1 &c);
char Three2One(string aa3);
float Angle(FV1 &p1,FV1 &p2,FV1 &p3);
float Dihedral(FV1 &p1,FV1 &p2,FV1 &p3,FV1 &p4);
float Distance(FV1 &p1,FV1 &p2);
float Sign(float c);
float VectorAngle(FV1 &c1,FV1 &c2);
float VectorDotProduct(FV1 &c1,FV1 &c2);
float VectorModulo(FV1 &c);
string One2Three(char aa1);
vector<double> main_faspr(const vector<double>& xyz, const char *residue_labels);
void MatrixByVector(FV2 &mtr,FV1 &vec,FV1 &cc);
//void ShowGraph(Graph &graph);
void Transpose(float **&tab,int i,int j);
void VectorAdd(FV1 &c1,FV1 &c2,FV1 &cc);
void VectorCrossProduct(FV1 &c1,FV1 &c2,FV1 &cc);
void VectorMinus(FV1 &cc);
void VectorMultiply(float par,FV1 &c);
void VectorSubtract(FV1 &c1,FV1 &c2,FV1 &cc);
// There are two function prototypes right before inject_coords()

//# Stand alone functions
//void ShowGraph(Graph &graph)
//{
//  for(Graph::iterator it=graph.begin(); it!=graph.end(); ++it){
//    cout <<it->first<<" =>";
//    for(set<int>::iterator it2=it->second.begin();it2!=it->second.end();++it2){
//      cout<<" "<<*it2;
//    }
//    cout<<endl;
//  }
//}

void Transpose(float **&tab,int i,int j)
{
  float** tmp=new float* [j];
  for(int k=0;k<j;k++){
    tmp[k]=new float [i];
    for(int l=0;l<i;l++){
      tmp[k][l]=tab[l][k];
    }
  }
  tab=tmp;
}


float Distance(FV1 &p1,FV1 &p2){
  return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
}

float VectorAngle(FV1 &c1,FV1 &c2){
  VectorNormalization(c1);
  VectorNormalization(c2);
  float angle=VectorDotProduct(c1,c2);
  if(angle>1.)angle=1.;
  if(angle<-1.)angle=-1.;
  return acos(angle)*RAD2DEG;
}

float Angle(FV1 &p1,FV1 &p2,FV1 &p3){
  FV1 c1,c2;
  VectorSubtract(p1,p2,c1);
  VectorSubtract(p3,p2,c2);
  return VectorAngle(c1,c2);
}

void VectorCrossProduct(FV1 &c1,FV1 &c2,FV1 &cc){
  cc.clear();
  cc.push_back(c1[1]*c2[2]-c1[2]*c2[1]);
  cc.push_back(c1[2]*c2[0]-c1[0]*c2[2]);
  cc.push_back(c1[0]*c2[1]-c2[0]*c1[1]);
}

float VectorDotProduct(FV1 &c1,FV1 &c2){
  return c1[0]*c2[0]+c1[1]*c2[1]+c1[2]*c2[2];
}

void VectorAdd(FV1 &c1,FV1 &c2,FV1 &cc){
  cc.clear();
  cc.push_back(c1[0]+c2[0]);
  cc.push_back(c1[1]+c2[1]);
  cc.push_back(c1[2]+c2[2]);
}

void VectorSubtract(FV1 &c1,FV1 &c2,FV1 &cc){
  cc.clear();
  cc.push_back(c1[0]-c2[0]);
  cc.push_back(c1[1]-c2[1]);
  cc.push_back(c1[2]-c2[2]);
}

void MatrixByVector(FV2 &mtr,FV1 &vec,FV1 &cc){
  cc.clear();
  cc.push_back(VectorDotProduct(mtr[0],vec));
  cc.push_back(VectorDotProduct(mtr[1],vec));
  cc.push_back(VectorDotProduct(mtr[2],vec));
}

bool VectorNormalization(FV1 &c){
  float len=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
  if(len==1.00){
    return true;
  }
  if(len<=DELTA){
    cout<<"warning! the vector is zero-vector"<<endl;
    return false;
  }
  c[0]/=len;
  c[1]/=len;
  c[2]/=len;
  return true;
}

void VectorMultiply(float par,FV1 &c)
{
  c[0]*=par;
  c[1]*=par;
  c[2]*=par;
}

float Dihedral(FV1 &p1,FV1 &p2,FV1 &p3,FV1 &p4){
  FV1 c1,c2,c3;
  VectorSubtract(p1,p2,c1);
  VectorSubtract(p2,p3,c2);
  VectorSubtract(p3,p4,c3);
  FV1 v1,v2,v3;
  VectorCrossProduct(c2,c1,v1);
  VectorCrossProduct(c3,c2,v2);
  VectorCrossProduct(v2,v1,v3);
  return  Sign(VectorDotProduct(v3,c2))*VectorAngle(v1,v2);
}

bool Internal2Cartesian(FV1 &c1,FV1 &c2,FV1 &c3,FV1 &p,FV1 &cc){
  FV1 d2(3,0);
  FV2 mtr(3,d2);
  d2[0]=p[0]*cos(DEG2RAD*p[1]);
  d2[1]=p[0]*cos(DEG2RAD*p[2])*sin(DEG2RAD*p[1]);
  d2[2]=p[0]*sin(DEG2RAD*p[2])*sin(DEG2RAD*p[1]);
  FV1 ab,bc,n;
  VectorSubtract(c2, c1, ab);
  VectorSubtract(c3, c2, bc);
  VectorNormalization(bc);
  VectorCrossProduct(ab,bc,n);
  VectorNormalization(n);
  VectorCrossProduct(n,bc,ab);
  for(int i=0;i<3;i++){
    mtr[i][0]=-bc[i];
    mtr[i][1]=ab[i];
    mtr[i][2]=n[i];
  }
  MatrixByVector(mtr,d2,bc);
  VectorAdd(c3,bc,cc);

  return true;
}

float VectorModulo(FV1 &c)
{
  return sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
}

float Sign(float c)
{
  if(c>0.)
    return 1.;
  else if(c<0.)
    return -1.;
  else
    return 0.;
}

bool VectorL2Norm(FV1 &data)
{
  int i,n=data.size();
  double tot=0;
  for(i=0;i<n;i++){
    if(data[i]!=0){
      tot+=data[i]*data[i];
    }
  }
  if(tot<=0)return false;
  tot=sqrt(tot);
  for(i=0;i<n;i++)
  data[i]/=tot;
  return true;
}

void VectorMinus(FV1 &cc)
{
  cc[0]=-cc[0];
  cc[1]=-cc[1];
  cc[2]=-cc[2];
}


//# Struct objects
struct Residue
{
  string name;
  char chID;
  int pos;
  char ins;
  SV1 atNames;
  string atTypes;
  FV2 xyz;
};

typedef vector<Residue> PV1;
typedef vector<vector<Residue> > PV2;

struct Topology
{
  int nchi;
  SV1 atnames;
  IV2 former3AtomIdx;
  FV2 ic;
};



/* ************************
FASPR Classes
******************* */


//# Class Structure
class Structure
{
public:
  string seq;
  PV1 pdb;
  int nres;

  ~Structure();
  void Pdb2Fas();
// removed from the original code
//  void ReadPDB(string &pdbfile);
//  void OutputPDB(PV1 &pdb);
//  void OutputPDB(PV1 &pdb,string &pdbfile);
//  void WritePDB(string &pdbfile);
};

void Structure::Pdb2Fas()
{
  seq.clear();
  for(int i=0;i<nres;i++){
    seq.push_back(Three2One(pdb[i].name));
  }
}


Structure::~Structure()
{
  nres=0;
  pdb.clear();
  seq.clear();
}

//# Rotamer Builder
class RotamerBuilder:public Structure
{
public:
  ~RotamerBuilder();
  void LoadSeq();
  //void LoadSeq(string &seqfile);    * removed from the original code*
  void LoadParameter();
  void LoadBBdepRotlib2010();
  void BuildSidechain();
  // removed from the original code
  // void RotlibFromBinary2Text(string binlibfile,string &txtlibfile);
  // void RotlibFromText2Binary(string &fulltextlib,string &binlibfile);
  void PhiPsi();
  void AssignSidechainTopology();
  int LoadBackbone(int site);
  void SideChain(int site,int rot,FV2& rxyz);

  IV1 subStat;
  IV1 nrots;
  FV1 phi,psi;
  FV3 chi;
  FV2 probRot;
  FV1 maxProb;
  PV1 stru;
  FV4 sc;
  map<char,float> wRotlib;
  map<char,Topology> sidechainTopo;
};

void RotamerBuilder::LoadSeq()
{
  Pdb2Fas();
  subStat.assign(nres,1);
  PhiPsi();
}

void RotamerBuilder::LoadParameter()
{
  wRotlib['C']=WGT_CYS;
  wRotlib['D']=WGT_ASP;
  wRotlib['E']=WGT_GLU;
  wRotlib['F']=WGT_PHE;
  wRotlib['H']=WGT_HIS;
  wRotlib['I']=WGT_ILE;
  wRotlib['K']=WGT_LYS;
  wRotlib['L']=WGT_LEU;
  wRotlib['M']=WGT_MET;
  wRotlib['N']=WGT_ASN;
  wRotlib['P']=WGT_PRO;
  wRotlib['Q']=WGT_GLN;
  wRotlib['R']=WGT_ARG;
  wRotlib['S']=WGT_SER;
  wRotlib['T']=WGT_THR;
  wRotlib['V']=WGT_VAL;
  wRotlib['W']=WGT_TRP;
  wRotlib['Y']=WGT_TYR;

}

void RotamerBuilder::LoadBBdepRotlib2010()
{
  map<char,int> Chin;//number of chi for the residue
  map<char,int> Rotn;//rotamer number of a residue
  map<char,int> Rotl;//skipped line in reading rotamer library
  Chin['R']=4;Chin['N']=2;Chin['D']=2;Chin['C']=1;Chin['Q']=3;Chin['E']=3;Chin['H']=2;Chin['I']=2;
  Chin['L']=2;Chin['K']=4;Chin['M']=3;Chin['F']=2;Chin['P']=2;Chin['S']=1;Chin['T']=1;Chin['W']=2;
  Chin['Y']=2;Chin['V']=1;

  Rotn['R']=75;Rotn['N']=36;Rotn['D']=18;Rotn['C']=3;Rotn['Q']=108;Rotn['E']=54;Rotn['H']=36;
  Rotn['I']=9;Rotn['L']=9;Rotn['K']=73;Rotn['M']=27;Rotn['F']=18;Rotn['P']=2;Rotn['S']=3;
  Rotn['T']=3;Rotn['W']=36;Rotn['Y']=18;Rotn['V']=3;
  Rotl['R']=0;Rotl['N']=75;Rotl['D']=111;Rotl['C']=129;Rotl['Q']=132;Rotl['E']=240;Rotl['H']=294;
  Rotl['I']=330;Rotl['L']=339;Rotl['K']=348;Rotl['M']=421;Rotl['F']=448;Rotl['P']=466;Rotl['S']=468;
  Rotl['T']=471;Rotl['W']=474;Rotl['Y']=510;Rotl['V']=528;

  int i,j,k;
  int phin,psin,Xn;
  float p,mp,ap;
  string buf;
  FV1 ftmp,stmp;
  FV2 ctmp;
  //fstream infile;
  short sht;
  //string rotfile=(string)"/home/joao/GitHub/FASPR/src/dun2010bbdep.bin";
  //infile.open(rotfile.c_str(),ios::in|ios::binary);
  //if(!infile){
    //cerr<<"error! cannot open rotamer library "<<rotfile<<endl;
    //exit(0);
  //}
  for(i=0;i<nres;i++){
    mp=0.;
    if(subStat[i]==0||seq[i]=='A'||seq[i]=='G'){
      goto CT;
    }

    p=phi[i]+180.;
    phin=(int)floor((p+5.)/10.);
    if(phin>=36){
      phin-=36;
    }
    p=psi[i]+180.;
    psin=(int)floor((p+5.)/10.);
    if(psin>=36){
      psin-=36;
    }
    Xn=Chin[seq[i]];
    infile.seekg((1296*Rotl[seq[i]]+(36*phin+psin)*Rotn[seq[i]])*20,ios::beg);
    ap=0.;
    for(j=0;j<Rotn[seq[i]];j++){
      infile.read((char*)&p,4);
      if(p<ROT_PROB_CUT_MIN)break;
      if(p==0.){
        p=1e-7;
      }
      if(p>mp){
        mp=p;
      }
      stmp.push_back(p);
      for(k=0;k<Xn;k++){
        infile.read((char*)&sht,2);
        ftmp.push_back(((float)sht)/10.);
      }
      infile.seekg((8-Xn)*2,ios::cur);//skip the other fields
      ctmp.push_back(ftmp);
      ftmp.clear();
      ap+=p;
      if(ap>ROT_PROB_CUT_ACC)break;
    }
CT:    chi.push_back(ctmp);
    ctmp.clear();
    nrots.push_back(stmp.size());
    probRot.push_back(stmp);
    stmp.clear();
    maxProb.push_back(mp);
  }
  //infile.close();
  Rotn.clear();
  Rotl.clear();
  Chin.clear();
}

// Removed from the original code
//void RotamerBuilder::RotlibFromBinary2Text(string binlibfile,string &txtlibfile)
//{
//  fstream infile;
//  string pp=PROGRAM_PATH+binlibfile;
//  infile.open(pp.c_str(),ios::in|ios::binary);
//  if(!infile){
//    cerr<<"error! cannot open rotamer library "<<pp<<endl;exit(0);
//  }
//  ofstream outfile;
//  outfile.open(txtlibfile.c_str());
//  outfile<<setiosflags(ios::fixed)<<setprecision(6);
//  while(!infile.eof()){
//    float p;
//    infile.read((char*)&p,4);
//    outfile<<setprecision(6)<<setw(8)<<p;
//    short val;
//    for(int i=0;i<8;++i){
//      infile.read((char*)&val,2);
//      outfile<<setprecision(1)<<setw(7)<<(float)val/10.0;
//    }
//    outfile<<endl;
//  }
//  infile.close();
//  outfile.close();
//}

/**************************************************************
RotlibFromText2Binary() converts an standard Dunbrack library
into a binary rotamer library for fast access
***************************************************************/
//void RotamerBuilder::RotlibFromText2Binary(string &fulltextlib,string &binlibfile)
//{
//  ofstream outfile;
//  outfile.open(binlibfile.c_str(),ios::binary);
//
//  char inpath[2048];
//  sprintf(inpath,"%s%s",PROGRAM_PATH.c_str(),fulltextlib.c_str());
//  FILE* infile=fopen(inpath,"r");
//  if(infile==NULL){
//    cerr<<"error! cannot open rotamer library "<<inpath<<endl;
//    exit(0);
//  }
//  char line[2048];
//  while(fgets(line,2048,infile)){
//    if(line[0]=='#') continue;
//    char aaname[4];
//    int phi,psi;
//    int a,b,c,d,e;
//    float prob,x1,x2,x3,x4,v1,v2,v3,v4;
//    sscanf(line,"%s %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f\n",
//      aaname,&phi,&psi,&a,&b,&c,&d,&e,&prob,&x1,&x2,&x3,&x4,&v1,&v2,&v3,&v4);
//    if(phi==180 || psi==180) continue;
//    outfile.write((char*)&prob,sizeof(float));
//    short s1=(short)(10*x1);
//    short s2=(short)(10*x2);
//    short s3=(short)(10*x3);
//    short s4=(short)(10*x4);
//    short s5=(short)(10*v1);
//    short s6=(short)(10*v2);
//    short s7=(short)(10*v3);
//    short s8=(short)(10*v4);
//    outfile.write((char*)&s1,sizeof(short));
//    outfile.write((char*)&s2,sizeof(short));
//    outfile.write((char*)&s3,sizeof(short));
//    outfile.write((char*)&s4,sizeof(short));
//    outfile.write((char*)&s5,sizeof(short));
//    outfile.write((char*)&s6,sizeof(short));
//    outfile.write((char*)&s7,sizeof(short));
//    outfile.write((char*)&s8,sizeof(short));
//  }
//  fclose(infile);
//  outfile.close();
//}
//

int RotamerBuilder::LoadBackbone(int site)
{
  if(subStat[site]==0)stru[site]=pdb[site];
  else if(seq[site]=='G')return 0;
  else{
    int i;
    SV1 atmp=sidechainTopo[seq[site]].atnames;
    for(i=0;i<atmp.size();i++){
      stru[site].atNames.push_back(atmp[i]);
      stru[site].atTypes.push_back(atmp[i][1]);
    }
    FV1 p,xyz;
    p=sidechainTopo[seq[site]].ic[0];
    Internal2Cartesian(stru[site].xyz[2],stru[site].xyz[0],stru[site].xyz[1],p,xyz);
    stru[site].xyz.push_back(xyz);
  }
  return 0;
}

void RotamerBuilder::SideChain(int site,int rot,FV2& rxyz)
{
  rxyz.clear();
  FV2 ftmp2=stru[site].xyz;
  int i,chipos=0;
  IV1 Inx;
  FV1 p,xyz,dihe=chi[site][rot];
  //starting from i=1 to skip the CB atoms
  for(i=1;i<sidechainTopo[seq[site]].atnames.size();i++){
    Inx=sidechainTopo[seq[site]].former3AtomIdx[i];
    p=sidechainTopo[seq[site]].ic[i];
    if(p[2]==181.0){
      p[2]=dihe[chipos];chipos++;
    }
    Internal2Cartesian(ftmp2[Inx[0]],ftmp2[Inx[1]],ftmp2[Inx[2]],p,xyz);
    ftmp2.push_back(xyz);
    rxyz.push_back(xyz);
    xyz.clear();
  }
}

void RotamerBuilder::PhiPsi()
{
  stru.clear();
  map<char,int> atNum;
  atNum['C']=atNum['S']=6;atNum['P']=atNum['T']=atNum['V']=7;atNum['M']=atNum['N']=atNum['D']=atNum['I']=atNum['L']=8;
  atNum['Q']=atNum['E']=atNum['K']=9;atNum['H']=10;atNum['R']=atNum['F']=11;atNum['Y']=12;atNum['W']=14;atNum['A']=5;
  int i,j,k,l,Inx[4];
  SV1 stmp(4," N  ");stmp[1]=" CA ";stmp[2]=" C  ";stmp[3]=" O  ";
  //default, set phi as -60 and psi as 60 (see SCWRL4 paper)
  phi.assign(nres,-60.);
  psi.assign(nres,60.);
  Residue rtmp;
  for(i=0;i<nres;i++){
    if(pdb[i].atNames.size()<3){
      cout<<"error! backbone is incomplete at residue: "<<pdb[i].name<<pdb[i].pos<<", please check!"<<endl;
      exit(0);
    }
    l=0;
    for(j=0;j<pdb[i].atNames.size();j++){
      for(k=0;k<4;k++){
        if(pdb[i].atNames[j]==stmp[k]){
          Inx[k]=j;
          l++;
        }
      }
      if(l==4)break;
    }
    if(l<4){
      if(l==3){
        l=0;
        for(j=0;j<pdb[i].atNames.size();j++){
          for(k=0;k<3;k++){
            if(pdb[i].atNames[j]==stmp[k]){
              Inx[k]=j;l++;
            }
          }
          if(l==3) goto ADDO;
        }
      }
      cout<<"error: backbone is incomplete at residue: "<<pdb[i].name<<pdb[i].pos<<", please check!"<<endl;
      exit(0);
ADDO:      Inx[3]=pdb[i].atNames.size();
      FV1 p(3,1.234),xyz;p[1]=120;p[2]=180;
      if(i==nres-1){
        Internal2Cartesian(pdb[i].xyz[Inx[0]],pdb[i].xyz[Inx[1]],pdb[i].xyz[Inx[2]],p,xyz);
      }
      else if(Distance(pdb[i].xyz[2],pdb[i+1].xyz[0])>BONDDIST){
        Internal2Cartesian(pdb[i].xyz[Inx[0]],pdb[i].xyz[Inx[1]],pdb[i].xyz[Inx[2]],p,xyz);
      }
      else{
        Internal2Cartesian(pdb[i+1].xyz[0],pdb[i].xyz[Inx[1]],pdb[i].xyz[Inx[2]],p,xyz);
      }
      pdb[i].xyz.push_back(xyz);
    }
    rtmp.chID=pdb[i].chID;rtmp.ins=pdb[i].ins;rtmp.name=pdb[i].name;rtmp.pos=pdb[i].pos;
    rtmp.atNames=stmp;rtmp.atTypes="NCCO";
    for(j=0;j<4;j++){
      rtmp.xyz.push_back(pdb[i].xyz[Inx[j]]);
    }
    stru.push_back(rtmp);
    rtmp.atNames.clear();
    rtmp.atTypes.clear();
    rtmp.xyz.clear();
    if(subStat[i]==0){
      if(atNum[seq[i]]>pdb[i].atNames.size()){
        cout<<"warning! incomplete residue "<<pdb[i].name<<pdb[i].pos<<endl;
        subStat[i]=1;
      }
    }
    for(j=0;j<3;j++){
      if(Distance(stru[i].xyz[j],stru[i].xyz[j+1])>BONDDIST){
        cout<<"warning! chain is discontinuous at residue "<<pdb[i].name<<stru[i].pos<<endl;
      }
    }
    if(i==0)continue;
    if(Distance(stru[i-1].xyz[2],stru[i].xyz[0])>BONDDIST){
      cout<<"warning! chain is discontinuous at ";
      cout<<seq[i-1]<<stru[i-1].pos<<" and "<<seq[i]<<stru[i].pos<<endl;
      continue;
    }
    psi[i-1]=Dihedral(stru[i-1].xyz[0],stru[i-1].xyz[1],stru[i-1].xyz[2],stru[i].xyz[0]);
    phi[i]=Dihedral(stru[i-1].xyz[2],stru[i].xyz[0],stru[i].xyz[1],stru[i].xyz[2]);
  }
  atNum.clear();
}

void RotamerBuilder::AssignSidechainTopology()
{
  Topology ttmp;
  IV1 itmp;
  FV1 ftmp;

  //Ala
  ttmp.nchi=0;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['A']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Arg
  ttmp.nchi=4;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(114.1);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD ");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(111.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NE ");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.461);ftmp.push_back(112);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CZ ");
  itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.33);ftmp.push_back(124.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NH1");
  itmp.push_back(6);itmp.push_back(7);itmp.push_back(8);itmp.push_back(9);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.326);ftmp.push_back(120);ftmp.push_back(0);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NH2");
  itmp.push_back(9);itmp.push_back(7);itmp.push_back(8);itmp.push_back(10);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.326);ftmp.push_back(120);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['R']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Asn
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.516);ftmp.push_back(112.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.231);ftmp.push_back(120.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" ND2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.328);ftmp.push_back(116.5);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['N']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Asp
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.516);ftmp.push_back(112.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.25);ftmp.push_back(118.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OD2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.25);ftmp.push_back(118.5);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['D']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Cys
  ttmp.nchi=1;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" SG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.807);ftmp.push_back(114);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['C']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Gln
  ttmp.nchi=3;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(114.1);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD ");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.516);ftmp.push_back(112.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OE1");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.231);ftmp.push_back(120.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NE2");
  itmp.push_back(7);itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.328);ftmp.push_back(116.5);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['Q']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Glu
  ttmp.nchi=3;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(114.1);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD ");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.516);ftmp.push_back(112.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OE1");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.25);ftmp.push_back(118.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OE2");
  itmp.push_back(7);itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.25);ftmp.push_back(118.5);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['E']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //His
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.5);ftmp.push_back(113.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" ND1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.378);ftmp.push_back(122.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.354);ftmp.push_back(131);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE1");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.32);ftmp.push_back(109.2);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NE2");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);itmp.push_back(9);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.374);ftmp.push_back(107.2);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['H']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Ile
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.546);ftmp.push_back(111.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG1");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.3);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG2");
  itmp.push_back(5);itmp.push_back(1);itmp.push_back(4);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.521);ftmp.push_back(110.5);ftmp.push_back(-122.6);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.516);ftmp.push_back(114);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['I']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Leu
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(116.3);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.521);ftmp.push_back(110.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.521);ftmp.push_back(110.5);ftmp.push_back(122.6);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['L']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Lys
  ttmp.nchi=4;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(114.1);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD ");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(111.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE ");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(111.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NZ ");
  itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.489);ftmp.push_back(112);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['K']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Met
  ttmp.nchi=3;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.52);ftmp.push_back(114.1);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" SD ");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.807);ftmp.push_back(112.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE ");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.789);ftmp.push_back(100.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['M']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Phe
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.5);ftmp.push_back(113.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.391);ftmp.push_back(120.7);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.391);ftmp.push_back(120.7);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE1");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.393);ftmp.push_back(120.7);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE2");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);itmp.push_back(9);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.393);ftmp.push_back(120.7);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CZ ");
  itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);itmp.push_back(10);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.39);ftmp.push_back(120);ftmp.push_back(0);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['F']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Pro
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(103.2);ftmp.push_back(-120);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.495);ftmp.push_back(104.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD ");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.507);ftmp.push_back(105.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['P']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Ser
  ttmp.nchi=1;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.417);ftmp.push_back(110.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['S']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Thr
  ttmp.nchi=1;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.542);ftmp.push_back(111.5);ftmp.push_back(-122);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OG1");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.433);ftmp.push_back(109.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG2");
  itmp.push_back(5);itmp.push_back(1);itmp.push_back(4);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.521);ftmp.push_back(110.5);ftmp.push_back(-120);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['T']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Trp
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.5);ftmp.push_back(113.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.365);ftmp.push_back(126.9);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.433);ftmp.push_back(126.7);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" NE1");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.375);ftmp.push_back(110.2);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE2");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);itmp.push_back(9);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.413);ftmp.push_back(107.2);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE3");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);itmp.push_back(10);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.4);ftmp.push_back(133.9);ftmp.push_back(0);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CZ2");
  itmp.push_back(5);itmp.push_back(7);itmp.push_back(9);itmp.push_back(11);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.399);ftmp.push_back(122.4);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CZ3");
  itmp.push_back(5);itmp.push_back(7);itmp.push_back(10);itmp.push_back(12);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.392);ftmp.push_back(118.7);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CH2");
  itmp.push_back(7);itmp.push_back(9);itmp.push_back(11);itmp.push_back(13);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.372);ftmp.push_back(117.5);ftmp.push_back(0);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['W']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Tyr
  ttmp.nchi=2;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.53);ftmp.push_back(110.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG ");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.511);ftmp.push_back(113.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD1");
  itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.394);ftmp.push_back(120.8);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CD2");
  itmp.push_back(6);itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.394);ftmp.push_back(120.8);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE1");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.392);ftmp.push_back(121.1);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CE2");
  itmp.push_back(4);itmp.push_back(5);itmp.push_back(7);itmp.push_back(9);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.392);ftmp.push_back(121.1);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CZ ");
  itmp.push_back(5);itmp.push_back(6);itmp.push_back(8);itmp.push_back(10);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.385);ftmp.push_back(119.5);ftmp.push_back(0);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" OH ");
  itmp.push_back(6);itmp.push_back(8);itmp.push_back(10);itmp.push_back(11);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.376);ftmp.push_back(119.7);ftmp.push_back(180);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['Y']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
  //Val
  ttmp.nchi=1;
  ttmp.atnames.push_back(" CB ");
  itmp.push_back(2);itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.546);ftmp.push_back(111.5);ftmp.push_back(-122.5);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG1");
  itmp.push_back(0);itmp.push_back(1);itmp.push_back(4);itmp.push_back(5);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.521);ftmp.push_back(110.5);ftmp.push_back(181);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  ttmp.atnames.push_back(" CG2");
  itmp.push_back(5);itmp.push_back(1);itmp.push_back(4);itmp.push_back(6);
  ttmp.former3AtomIdx.push_back(itmp);itmp.clear();
  ftmp.push_back(1.521);ftmp.push_back(110.5);ftmp.push_back(122.6);
  ttmp.ic.push_back(ftmp);ftmp.clear();
  sidechainTopo['V']=ttmp;
  ttmp.atnames.clear();ttmp.former3AtomIdx.clear();ttmp.ic.clear();
}


RotamerBuilder::~RotamerBuilder()
{
  stru.clear();
  sc.clear();
  phi.clear();
  psi.clear();
  chi.clear();
  probRot.clear();
  maxProb.clear();
  nrots.clear();
  subStat.clear();
  sidechainTopo.clear();
}


void RotamerBuilder::BuildSidechain()
{
  LoadParameter();
  LoadBBdepRotlib2010();
  AssignSidechainTopology();
  int i,j;
  for(i=0;i<nres;i++){
    LoadBackbone(i);
  }

  FV2 rxyz;
  FV3 rtmp;
  for(i=0;i<nres;i++){
    if(nrots[i]!=0){
      for(j=0;j<nrots[i];j++){
        SideChain(i,j,rxyz);
        rtmp.push_back(rxyz);
      }
    }
    sc.push_back(rtmp);
    rtmp.clear();
  }  
  sidechainTopo.clear();
  subStat.clear();
}


//# SelfEnergy
class SelfEnergy:public RotamerBuilder
{
public:
  IV1 bestrot;
  FV2 eTableSelf;
  IV2 conMap;
  FV2 radius,depth;
  IV2 atomIdx;

  ~SelfEnergy();
  void AssignConMap();

  void SetVdwPar();
  float VDWType(int a,int b,float rij,float dist);
  float VDWEnergyAtomAndAtom(float ddash,float eij);
  float RotamerPreferenceEnergy(int site,int rot);

  void EnergySidechainAndBackbone(int site,FV1 &ener);
  void EnergyRotamerSidechainAndFixedSidechain(int site,FV1 &ener);
  void CalcSelfEnergy();
  float HbondEnergyAtomAndAtom(FV1 &DBxyz,FV1 &Dxyz,FV1 &Axyz,FV1 &ABxyz,float Dangle,float Aangle);
  float SSbondEnergyAtomAndAtom(FV1 &CA1xyz,FV1 &CB1xyz,FV1 &SG1xyz,FV1 &SG2xyz,FV1 &CB2xyz,FV1 &CA2xyz);
  float EnergyPolarSidechainAndBackbone(int site1,int rot1,int site2);
  float EnergyPolarSidechainAndSidechain(int site1,int rot1,int site2,int rot2);
};



SelfEnergy::~SelfEnergy(){
  bestrot.clear();
  eTableSelf.clear();
  conMap.clear();
  radius.clear();
  depth.clear();
  atomIdx.clear();
}

void SelfEnergy::AssignConMap()
{
  bestrot.assign(nres,-1);
  map<char,float> res_rad;
  res_rad['G']=2.4;
  res_rad['A']=res_rad['C']=res_rad['D']=res_rad['I']=res_rad['L']=res_rad['N']=res_rad['P']=res_rad['S']=res_rad['T']=res_rad['V']=3.2;
  res_rad['Q']=res_rad['E']=res_rad['H']=3.7;
  res_rad['M']=res_rad['F']=4.3;
  res_rad['K']=5.0;
  res_rad['W']=5.3;
  res_rad['Y']=5.7;
  res_rad['R']=6.2;

  IV1 ivtmp;
  conMap.assign(nres,ivtmp);
  for(int i=0;i<nres;i++){
    if(nrots[i]<2){
      if(nrots[i]==1){
        bestrot[i]=0;
      }
      continue;
    }
    for(int j=0;j<nres;j++){
      if(j==i){
        continue;
      }
      int cb=(seq[j]=='G'||seq[j]=='A')?1:4;
      if(Distance(stru[i].xyz[4],stru[j].xyz[cb])<res_rad[seq[i]]+res_rad[seq[j]]+RESI_DIST_CUT &&
        Distance(stru[i].xyz[4],stru[j].xyz[cb])<Distance(stru[i].xyz[1],stru[j].xyz[1])+CACB_DIST_CUT){
        ivtmp.push_back(j);
      }
    }
    conMap[i]=ivtmp;
    ivtmp.clear();
  }
  res_rad.clear();
}


void SelfEnergy::SetVdwPar()
{
  //set vdw parameters: [X][0]->radius; [X][1]->well-depth
  FV1 vtmp(2,2.0);
  FV2 vdwpar(18,vtmp);//18 atom types in CHARMM19
  vdwpar[ 0][0]=RADIUS_C1; vdwpar[ 0][1]=DEPTH_C1;//1  mainchain CA
  vdwpar[ 1][0]=RADIUS_C2; vdwpar[ 1][1]=DEPTH_C2;//2  mainchain C
  vdwpar[ 2][0]=RADIUS_C3; vdwpar[ 2][1]=DEPTH_C3;//3  CH1
  vdwpar[ 3][0]=RADIUS_C4; vdwpar[ 3][1]=DEPTH_C4;//4  CH2
  vdwpar[ 4][0]=RADIUS_C5; vdwpar[ 4][1]=DEPTH_C5;//5  CH3
  vdwpar[ 5][0]=RADIUS_C6; vdwpar[ 5][1]=DEPTH_C6;//6  Aromatic CH1 and C (Phe/Tyr/Trp)
  vdwpar[ 6][0]=RADIUS_C7; vdwpar[ 6][1]=DEPTH_C7;//7  CO,COO,NCNN
  vdwpar[ 7][0]=RADIUS_C8; vdwpar[ 7][1]=DEPTH_C8;//8  Cys CB
  vdwpar[ 8][0]=RADIUS_N1; vdwpar[ 8][1]=DEPTH_N1;//9  mainchain NH
  vdwpar[ 9][0]=RADIUS_N2; vdwpar[ 9][1]=DEPTH_N2;//10 His/Arg/Trp NH, Asn/Gln/Arg NH2, Lys NH3
  vdwpar[10][0]=RADIUS_N3; vdwpar[10][1]=DEPTH_N3;//11 His C=N-C
  vdwpar[11][0]=RADIUS_N4; vdwpar[11][1]=DEPTH_N4;//12 Pro N
  vdwpar[12][0]=RADIUS_O1; vdwpar[12][1]=DEPTH_O1;//13 mainchain O
  vdwpar[13][0]=RADIUS_O2; vdwpar[13][1]=DEPTH_O2;//14 sidechain C=O
  vdwpar[14][0]=RADIUS_O3; vdwpar[14][1]=DEPTH_O3;//15 sidechain COO
  vdwpar[15][0]=RADIUS_O4; vdwpar[15][1]=DEPTH_O4;//16 Ser/Thr/Tyr OH
  vdwpar[16][0]=RADIUS_S1; vdwpar[16][1]=DEPTH_S1;//17 Cys S
  vdwpar[17][0]=RADIUS_S2; vdwpar[17][1]=DEPTH_S2;//18 Met S

  map<char,IV1> vAtomIdx;//vector of atom index
  IV1 inx(4,9);
  inx[1]=1;inx[2]=2;inx[3]=13;vAtomIdx['G']=inx;//4 atoms
  inx.push_back(5);vAtomIdx['A']=inx;//5 atoms
  inx[4]=8;inx.push_back(17); vAtomIdx['C']=inx;//6 atoms
  inx[4]=4;inx[5]=16; vAtomIdx['S']=inx;//6 atoms
  inx[5]=4;inx.push_back(4);inx[0]=12; vAtomIdx['P']=inx;//7 atoms
  inx[0]=9;inx[4]=3;inx[5]=16;inx[6]=5; vAtomIdx['T']=inx;//7 atoms
  inx[5]=5; vAtomIdx['V']=inx;//7 atoms
  inx[4]=4;inx[5]=4;inx[6]=18;inx.push_back(5);vAtomIdx['M']=inx;//8 atoms
  inx[5]=7;inx[6]=14;inx[7]=10;  vAtomIdx['N']=inx;//8 atoms
  inx[6]=15;inx[7]=15; vAtomIdx['D']=inx;//8 atooms
  inx[5]=3;inx[6]=5;inx[7]=5;  vAtomIdx['L']=inx;//8 atoms
  inx[4]=3;inx[5]=4; vAtomIdx['I']=inx;//8 atoms
  inx[4]=4;inx[6]=7;inx[7]=14;inx.push_back(10);vAtomIdx['Q']=inx;//9 atoms
  inx[7]=15;inx[8]=15; vAtomIdx['E']=inx;//9 atoms
  inx[6]=4;inx[7]=4;inx[8]=10; vAtomIdx['K']=inx;//9 atoms
  //for His, index 6 stands for ND1, index 9 stands for NE2
  inx[5]=6;inx[6]=10;inx[7]=6;inx[8]=6;inx.push_back(11);vAtomIdx['H']=inx; //10 atoms
  inx[5]=4;inx[6]=4;inx[7]=10;inx[8]=7;inx[9]=10;inx.push_back(10);vAtomIdx['R']=inx;//11 atoms
  inx[5]=6;inx[6]=6;inx[7]=6;inx[8]=6;inx[9]=6;inx[10]=6;vAtomIdx['F']=inx;//11 atoms
  inx.push_back(16);vAtomIdx['Y']=inx;//12 atoms
  inx[8]=10;inx[11]=6;inx.push_back(6);inx.push_back(6);vAtomIdx['W']=inx;//14 atoms
  for(int i=0;i<nres;i++){
    FV1 radtmp,depthtmp;
    IV1 ati_tmp=vAtomIdx[seq[i]];
    for(int j=0;j<ati_tmp.size();j++){
      //get the atom type index of atom j on residue i
      radtmp.push_back(vdwpar[ati_tmp[j]-1][0]);
      depthtmp.push_back(vdwpar[ati_tmp[j]-1][1]);      
    }
    radius.push_back(radtmp);
    depth.push_back(depthtmp);
    atomIdx.push_back(ati_tmp);
    radtmp.clear();
    depthtmp.clear();
    ati_tmp.clear();
  }
}

float SelfEnergy::VDWType(int a,int b,float rij,float dist)
{
  int i1=a,i2=b;
  if(i1>i2) swap(i1,i2);
  if(i1==9 && (i2==11 || i2==13 || i2==14 || i2==15 || i2==16)){
    if(dist>MIN_HBOND_DIST && dist<MAX_HBOND_DIST){
      return OPT_HBOND_DIST;
    }
  }
  else if(i1==10 && (i2==11 || i2==13 || i2==14 || i2==15 || i2==16)){
    if(dist>MIN_HBOND_DIST && dist<MAX_HBOND_DIST){
      return OPT_HBOND_DIST;
    }
  }
  else if(i1==11 && i2==16){
    if(dist>MIN_HBOND_DIST && dist<MAX_HBOND_DIST){
      return OPT_HBOND_DIST;
    }
  }
  else if((i1==13 || i1==14 || i1==15) && i2==16){
    if(dist>MIN_HBOND_DIST && dist<MAX_HBOND_DIST){
      return OPT_HBOND_DIST;
    }
  }
  else if(i1==16 && i2==16){
    if(dist>MIN_HBOND_DIST && dist<MAX_HBOND_DIST){
      return OPT_HBOND_DIST;
    }
  }
  else if(i1==17 && i2==17){//disulfide bond
    if(dist>MIN_SSBOND_DIST && dist<MAX_SSBOND_DIST){
      return OPT_SSBOND_DIST;
    }
  }
  return rij;
}


//dstar=dij/rij
float SelfEnergy::VDWEnergyAtomAndAtom(float dstar,float eij)
{
  if(dstar>DSTAR_MAX_CUT) return 0.;
  else if(dstar>1){
    double energy = 4*eij*(pow(1/dstar,12)-pow(1/dstar,6));
    return energy;
  }
  else if(dstar>DSTAR_MIN_CUT){
    return VDW_REP_CUT*(dstar-1)/(DSTAR_MIN_CUT-1);
  }
  else{
    return VDW_REP_CUT;
  }
}

float SelfEnergy::HbondEnergyAtomAndAtom(FV1 &DBxyz,FV1 &Dxyz,FV1 &Axyz,FV1 &ABxyz,float Dangle,float Aangle){
  float dist=Distance(Dxyz,Axyz);
  if(dist>MAX_HBOND_DIST || dist<MIN_HBOND_DIST) return 0.;
  float theta=Angle(DBxyz,Dxyz,Axyz);
  if(theta<MIN_HBOND_THETA) return 0.;
  float phi=Angle(Dxyz,Axyz,ABxyz);
  if(phi<MIN_HBOND_PHI) return 0.;
  float fdist=5.*pow(OPT_HBOND_DIST/dist,12)-6.*pow(OPT_HBOND_DIST/dist,10);
  float ftheta=cos((theta-Dangle)*DEG2RAD)*cos((theta-Dangle)*DEG2RAD);
  float fphi=cos((phi-Aangle)*DEG2RAD)*cos((phi-Aangle)*DEG2RAD);
  float energy = fdist*ftheta*fphi;
  energy *= WGT_HBOND;
  if(energy>0.) energy=0.;
  return energy;
}

float SelfEnergy::SSbondEnergyAtomAndAtom(FV1 &CA1xyz,FV1 &CB1xyz,FV1 &SG1xyz,FV1 &SG2xyz,FV1 &CB2xyz,FV1 &CA2xyz){
  float dist=Distance(SG1xyz,SG2xyz);
  if(dist>MAX_SSBOND_DIST || dist<MIN_SSBOND_DIST) return 0.;
  float ang1=Angle(CB1xyz,SG1xyz,SG2xyz);
  if(ang1>MAX_SSBOND_ANGL || ang1<MIN_SSBOND_ANGL) return 0;
  float ang2=Angle(CB2xyz,SG2xyz,SG1xyz);
  if(ang2>MAX_SSBOND_ANGL || ang2<MIN_SSBOND_ANGL) return 0.;
  float torsion1=Dihedral(CB1xyz,SG1xyz,SG2xyz,CB2xyz);
  float energy=100.*(dist-OPT_SSBOND_DIST)*(dist-OPT_SSBOND_DIST)-4.
    +0.01*(ang1-OPT_SSBOND_ANGL)*(ang1-OPT_SSBOND_ANGL)-2.
    +0.01*(ang2-OPT_SSBOND_ANGL)*(ang2-OPT_SSBOND_ANGL)-2.
    +2.*cos(2.*torsion1*DEG2RAD);
  energy *= WGT_SSBOND;
  if(energy>0.) energy=0.;
  return energy;
}


float SelfEnergy::RotamerPreferenceEnergy(int site,int rot)
{
  float elib=-1.*log(probRot[site][rot]/maxProb[site]);
  if(elib>5.){
    elib=5.;
  }
  return wRotlib[seq[site]]*elib;
}

float SelfEnergy::EnergyPolarSidechainAndBackbone(int site1,int rot1,int site2){
  float energy=0.;
  if(site1-site2<=1 && site1-site2>=-1) return energy;
  if(seq[site1]=='D'){
    energy = HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
  }
  else if(seq[site1]=='E'){
    energy = HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
  }
  else if(seq[site1]=='K'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],stru[site2].xyz[3],stru[site2].xyz[2],109.5,120.);
  }
  else if(seq[site1]=='R'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.)+
      HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.)+
      HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.);
  }
  else if(seq[site1]=='W'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.);
  }
  else if(seq[site1]=='H'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.);
  }
  else if(seq[site1]=='N'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.);
  }
  else if(seq[site1]=='Q'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.);
  }
  else if(seq[site1]=='S'){
    energy = HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],stru[site2].xyz[3],stru[site2].xyz[2],109.5,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5);
  }
  else if(seq[site1]=='T'){
    energy = HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],stru[site2].xyz[3],stru[site2].xyz[2],109.5,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5);
  }
  else if(seq[site1]=='Y'){
    energy = HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],stru[site2].xyz[3],stru[site2].xyz[2],120.,120.)+
      HbondEnergyAtomAndAtom(stru[site2].xyz[1],stru[site2].xyz[0],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.);
  }
  return energy;
}


float SelfEnergy::EnergyPolarSidechainAndSidechain(int site1,int rot1,int site2,int rot2){
  float energy=0.;
  if(seq[site1]=='C'){
    if(seq[site2]=='C'){
      energy = SSbondEnergyAtomAndAtom(stru[site1].xyz[1],stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][0],stru[site2].xyz[4],stru[site2].xyz[1]);
    }
  }
  else if(seq[site1]=='D'){
    if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][1],sc[site1][rot1][0],109.5,120.)+
        HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][2],sc[site1][rot1][0],109.5,120.);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][1],sc[site1][rot1][0],109.5,120.)+
        HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][2],sc[site1][rot1][0],109.5,120.);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
  }
  else if(seq[site1]=='E'){
    if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][3],sc[site1][rot1][1],109.5,120.);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.)+
        HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][3],sc[site1][rot1][1],109.5,120.);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.)+
        HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][3],sc[site1][rot1][1],109.5,120.);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][3],sc[site1][rot1][1],120.,120.);
    }
  }
  else if(seq[site1]=='K'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][1],sc[site2][rot2][0],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][0],109.5,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][1],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][3],sc[site2][rot2][1],109.5,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][4],sc[site2][rot2][2],109.5,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][1],sc[site2][rot2][0],109.5,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][1],109.5,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][0],stru[site2].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][0],stru[site2].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][2],sc[site1][rot1][3],sc[site2][rot2][6],sc[site2][rot2][5],109.5,120.);
    }
  }
  else if(seq[site1]=='R'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][4],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][3],sc[site1][rot1][5],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][2],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.);
    }
  }
  else if(seq[site1]=='W'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.);
    }
  }
  else if(seq[site1]=='H'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][4],sc[site1][rot1][2],109.5,120.);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][4],sc[site1][rot1][2],120.,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][4],sc[site1][rot1][2],120.,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][1],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][4],sc[site1][rot1][2],120.,120.);
    }
  }
  else if(seq[site1]=='N'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],109.5,120.);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][1],sc[site1][rot1][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][1],sc[site1][rot1][0],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][1],sc[site1][rot1][0],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][1],sc[site1][rot1][0],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][0],sc[site1][rot1][2],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.);
    }
  }
  else if(seq[site1]=='Q'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][2],sc[site1][rot1][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][0],stru[site2].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][2],sc[site1][rot1][1],109.5,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][1],sc[site1][rot1][3],sc[site2][rot2][6],sc[site2][rot2][5],120.,120.);
    }
  }
  else if(seq[site1]=='S'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][0],stru[site1].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][4],sc[site2][rot2][2],109.5,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][1],sc[site2][rot2][0],109.5,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][2],sc[site2][rot2][1],109.5,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][0],stru[site1].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][0],stru[site1].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][0],stru[site1].xyz[4],120.0,109.5);
    }
  }
  else if(seq[site1]=='T'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][0],stru[site1].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][4],sc[site2][rot2][2],109.5,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][1],sc[site2][rot2][0],109.5,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][0],stru[site1].xyz[4],120.,109.5)+
        HbondEnergyAtomAndAtom(stru[site1].xyz[4],sc[site1][rot1][0],sc[site2][rot2][2],sc[site2][rot2][1],109.5,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][0],stru[site1].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][0],stru[site1].xyz[4],109.5,109.5);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][0],stru[site1].xyz[4],120.0,109.5);
    }
  }
  else if(seq[site1]=='Y'){
    if(seq[site2]=='D'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][2],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='E'){
      energy = HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][3],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='K'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][2],sc[site2][rot2][3],sc[site1][rot1][6],sc[site1][rot1][5],109.5,120.);
    }
    else if(seq[site2]=='R'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][4],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][3],sc[site2][rot2][5],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][2],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.);
    }
    else if(seq[site2]=='W'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.);
    }
    else if(seq[site2]=='H'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][1],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][4],sc[site2][rot2][2],120.,120.);
    }
    else if(seq[site2]=='N'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][0],sc[site2][rot2][2],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][1],sc[site2][rot2][0],120.,120.);
    }
    else if(seq[site2]=='Q'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][1],sc[site2][rot2][3],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.)+
        HbondEnergyAtomAndAtom(sc[site1][rot1][5],sc[site1][rot1][6],sc[site2][rot2][2],sc[site2][rot2][1],120.,120.);
    }
    else if(seq[site2]=='S'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][6],sc[site1][rot1][5],109.5,120.);
    }
    else if(seq[site2]=='T'){
      energy = HbondEnergyAtomAndAtom(stru[site2].xyz[4],sc[site2][rot2][0],sc[site1][rot1][6],sc[site1][rot1][5],109.5,120.);
    }
    else if(seq[site2]=='Y'){
      energy = HbondEnergyAtomAndAtom(sc[site2][rot2][5],sc[site2][rot2][6],sc[site1][rot1][6],sc[site1][rot1][5],120.,120.);
    }
  }
  return energy;
}

void SelfEnergy::EnergySidechainAndBackbone(int site,FV1 &ener)
{
  ener.assign(nrots[site],0.);
  //residue internal energy
  for(int j=0;j<4;j++){
    for(int k=6;k<stru[site].atTypes.size();k++){
      float eij=sqrt(depth[site][k]*depth[site][j]);
      for(int r=0;r<nrots[site];r++){
        float dist=Distance(stru[site].xyz[j],sc[site][r][k-5]);
        //if(dist>VDW_DIST_CUT) continue;
        float rij=VDWType(atomIdx[site][k],atomIdx[site][j],radius[site][k]+radius[site][j],dist);
        ener[r]+=VDWEnergyAtomAndAtom(dist/rij,eij);
      }
    }
  }

  for(int i=0;i<conMap[site].size();i++){
    int ipair=conMap[site][i];
    //side-chain vs. main-chain vdW interactions
    int n=seq[ipair]=='G'?4:5;
    for(int j=0;j<n;j++){
      for(int k=5;k<stru[site].atTypes.size();k++){
        float eij=sqrt(depth[site][k]*depth[ipair][j]);
        for(int r=0;r<nrots[site];r++){
          float dist=Distance(stru[ipair].xyz[j],sc[site][r][k-5]);
          //if(dist>VDW_DIST_CUT) continue;
          float rij=VDWType(atomIdx[site][k],atomIdx[ipair][j],radius[site][k]+radius[ipair][j],dist);
          ener[r]+=VDWEnergyAtomAndAtom(dist/rij,eij);
        }
      }
    }
    //side-chain vs. main-chain Hbond
    for(int r=0;r<nrots[site];r++){
      ener[r]+=EnergyPolarSidechainAndBackbone(site,r,ipair);
    }
  }
}

void SelfEnergy::EnergyRotamerSidechainAndFixedSidechain(int site,FV1 &ener)
{
  for(int i=0;i<conMap[site].size();i++){
    int ipair=conMap[site][i];
    if(bestrot[ipair]==-1){
      continue;
    }
    for(int j=5;j<stru[ipair].atTypes.size();j++){
      for(int k=5;k<stru[site].atTypes.size();k++){
        float eij=sqrt(depth[site][k]*depth[ipair][j]);
        for(int r=0;r<nrots[site];r++){
          float dist=Distance(sc[ipair][0][j-5],sc[site][r][k-5]);
          //if(dist>VDW_DIST_CUT) continue;
          float rij=VDWType(
        atomIdx[site][k],
        atomIdx[ipair][j],
        radius[site][k]+radius[ipair][j],
        dist);
          ener[r]+=VDWEnergyAtomAndAtom(dist/rij,eij);
        }
      }
    }
    //side-chain vs. side-chain polar energy
    for(int r=0;r<nrots[site];r++){
      ener[r]+=EnergyPolarSidechainAndSidechain(site,r,ipair,0);
    }
  }
}

void SelfEnergy::CalcSelfEnergy()
{
  SetVdwPar();
  AssignConMap();
  
  FV1 ftmp;
  eTableSelf.assign(nres,ftmp);
  //calculate self-energy
  for(int i=0;i<nres;i++){
    if(nrots[i]<2){
      continue;
    }
    EnergySidechainAndBackbone(i,eTableSelf[i]);
  }

  int nPosFixedNonAlaGly=0;
  for(int i=0;i<nres;i++){
    if(nrots[i]==1){
      nPosFixedNonAlaGly++;
    }
    if(nrots[i]<2){
      continue;
    }
    float emin=1e8;
    for(int j=0;j<nrots[i];j++){
      if(eTableSelf[i][j]<emin){
        emin=eTableSelf[i][j];
      }
    }
    for(int j=0;j<nrots[i];j++){
      if(eTableSelf[i][j]>emin+SEC_CUT){
        nrots[i]--;
        eTableSelf[i].erase(eTableSelf[i].begin()+j);
        probRot[i].erase(probRot[i].begin()+j);
        chi[i].erase(chi[i].begin()+j);
        sc[i].erase(sc[i].begin()+j);
        j--;
      }
    }
    if(nrots[i]==1){
      //cout<<"fixed residue index: "<<i<<endl;
      bestrot[i]=0;
      nPosFixedNonAlaGly++;
    }
  }
  //cout<<"#residues fixed during self-energy-check: "<<nPosFixedNonAlaGly<<endl;
  
  //update the self-energy table
  for(int i=0;i<nres;i++){
    if(nrots[i]<2){
      continue;
    }
    for(int j=0;j<nrots[i];j++){
      eTableSelf[i][j]+=RotamerPreferenceEnergy(i,j);
    }
    EnergyRotamerSidechainAndFixedSidechain(i,eTableSelf[i]);
  }
}


// Pair Energy

class PairEnergy:public SelfEnergy{
public:
  float**** eTablePair;

  ~PairEnergy();
  void CalcPairEnergy();
  bool EnergyRotamerSidechainAndRotamerSidechain(int site1,int site2,float **&tab);
  //void ShowPairEnergy();
  //void ShowPairEnergy(int site1, int site2);
};



PairEnergy::~PairEnergy()
{
  int i,j,k;
  for(i=0;i<nres;i++){
    for(j=0;j<nres;j++){
      if(eTablePair[i][j]!=NULL){
        for(k=0;k<nrots[i];k++){
          delete [] eTablePair[i][j][k];
        }
      }
      delete [] eTablePair[i][j];
    }
    delete [] eTablePair[i];
  }
  delete [] eTablePair;
}



bool PairEnergy::EnergyRotamerSidechainAndRotamerSidechain(int site1,int site2,float **&tab)
{
  tab= new float* [nrots[site1]];
  for(int i=0;i<nrots[site1];i++){
    tab[i]=new float [nrots[site2]];
    for(int j=0;j<nrots[site2];j++){
      tab[i][j]=0;
    }
  }
  
  float aij,eij,dist;
  for(int i=5;i<stru[site1].atTypes.size();i++){
    for(int j=5;j<stru[site2].atTypes.size();j++){
      float eij=sqrt(depth[site1][i]*depth[site2][j]);
      for(int k=0;k<nrots[site1];k++){
        for(int l=0;l<nrots[site2];l++){
          float dist=Distance(sc[site1][k][i-5],sc[site2][l][j-5]);
          float rij=VDWType(atomIdx[site1][i],atomIdx[site2][j],radius[site1][i]+radius[site2][j],dist);
          tab[k][l]+=VDWEnergyAtomAndAtom(dist/rij,eij);
        }
      }
    }
  }

  //side-chain vs. main-chain polar energy
  for(int k=0;k<nrots[site1];k++){
    for(int l=0;l<nrots[site2];l++){
      tab[k][l]+=EnergyPolarSidechainAndSidechain(site1,k,site2,l);
    }
  }
  
  for(int i=0;i<nrots[site1];i++){
    for(int j=0;j<nrots[site2];j++){
      if(tab[i][j]!=0){
        return true;
      }
    }
  }
  
  return false;
}

void PairEnergy::CalcPairEnergy()
{
  int i,j,k,l;
  eTablePair= new float***[nres];
  for(i=0;i<nres;i++){
    eTablePair[i]= new float** [nres];
    for(j=0;j<nres;j++)
    eTablePair[i][j]=NULL;
  }
  float **tab;
  int ip;
  for(i=0;i<nres-1;i++){
    if(nrots[i]<2) continue;
    for(j=0;j<conMap[i].size();j++){
      ip=conMap[i][j];
      if(nrots[ip]<2) continue;
      if(ip<i) continue;
      if(EnergyRotamerSidechainAndRotamerSidechain(i,ip,tab)){
        eTablePair[i][ip]=tab;
        Transpose(tab,nrots[i],nrots[ip]);
        eTablePair[ip][i]=tab;
      }
    }
  }
}

// disabled from original,
//void PairEnergy::ShowPairEnergy()
//{
//  for(int i=0;i<nres-1;i++){
//    for(int j=i+1;j<nres;j++){
//      if(eTablePair[i][j] != NULL){
//        cout<<"pairwise energies between site "<<i<<" and "<<j<<":";
//        for(int k=0;k<nrots[i];k++){
//          for(int s=0;s<nrots[j];s++){
//            cout<<" "<<eTablePair[i][j][k][s];
//          }
//        }
//        cout<<endl;
//      }
//    }
//  }
//}



//# template classes
template <class Elem> class RankElement{
public:
  RankElement(int a,Elem b):idx(a),element(b){}
  int idx;
  Elem element;
  bool operator < (const RankElement &m)const {
    return element<m.element;
  }
  bool operator > (const RankElement &m)const {
    return element>m.element;
  }
};


//sort elements according to the size
template <class Elem> class RankSize{
public:
  RankSize(int a,Elem b):idx(a),element(b){}
  int idx;
  Elem element;
  bool operator < (const RankSize &m)const {
    return element.size()<m.element.size();
  }
  bool operator > (const RankSize &m)const {
    return element.size()>m.element.size();
  }
};


//# Bag Class

class Bag
{
public:
  set<int> left;
  set<int> right;
  set<int> total;
  int parentBagIdx;
  set<int> childBagIdx;
  int type;
  int childCounter;

  /* data structure to record solution*/
  IV1 lsites;
  IV1 rsites;
  IV1 tsites;
  IV2 lrots;
  IV2 rrots;
  FV2 Etcom;
  IV2 Rtcom;
  IV1 indices;
  bool deployFlag;

  /* for backtrack*/
  float EGMEC;
  IV1 RLGMEC;
  IV1 RRGMEC;

  Bag(){
    type=Leaf;
    deployFlag=false;
  }
  //void ShowBag();
};

//void Bag::ShowBag()
//{
//  cout<<"(";
//  for(set<int>::iterator it=left.begin();it!=left.end();++it){
//    cout<<*it<<" ";
//  }
//  cout<<"|";
//  for(set<int>::iterator it=right.begin();it!=right.end();++it){
//    cout<<" "<<*it;
//  }
//  cout<<")";
//}


//# TreeDecomposition Class

class TreeDecomposition
{
public:
  vector <Bag> bags;
  vector <Bag> connBags;

  void Subgraph2TreeDecomposition(int index,Graph &graph);
  void MergeBags(int depth);
  int CheckTreewidth();

};


void TreeDecomposition::Subgraph2TreeDecomposition(int index,Graph &graph)
{
  //sort the vertices in the cluster by the number of edges
  int ver=1;
  set <int> neibs;
  RankSize <set <int> > nbsize(ver, neibs);
  vector<RankSize< set<int> > > sortgraph;
  for(Graph::iterator it=graph.begin(); it!=graph.end(); ++it){
    ver=it->first;
    neibs=it->second;
    nbsize.idx=ver;
    nbsize.element=neibs;
    sortgraph.push_back(nbsize);
  }
  stable_sort(sortgraph.begin(), sortgraph.end(), less <RankSize<set<int> > >());
  //cout<<endl<<"*subgraph "<<index<<" sorted (low-to-high degree):"<<endl;
  for(int i=0;i<sortgraph.size();++i){
    RankSize<set<int> > &st=sortgraph[i];
    //cout <<st.idx<<" => ";
    //for(set<int>::iterator it2=st.element.begin();it2!=st.element.end();++it2){
    //  cout<<" "<<*it2;
    //}
    //cout<<endl;
  }

  //construct the bags using a loop
  while (sortgraph.size()>0){
    //add to bag
    Bag newbag;
    RankSize< set<int> > &st=sortgraph[0];
    newbag.right.insert(st.idx);
    newbag.total.insert(st.idx);
    for(set<int>::iterator it2=st.element.begin();it2!=st.element.end();++it2){
      newbag.left.insert(*it2);
      newbag.total.insert(*it2);
    }
    bags.push_back(newbag);

    //1.remove edges from the sorted graph
    sortgraph.erase(sortgraph.begin());
    for(set<int>::iterator it3=newbag.left.begin();it3!=newbag.left.end();++it3){
      for(vector<RankSize<set<int> > >::iterator it4=sortgraph.begin();it4!=sortgraph.end();++it4){
        if(*it3==it4->idx){
          for(set<int>::iterator it5=it4->element.begin();it5!=it4->element.end();){
            if(*(newbag.right.begin())==*it5){
              it4->element.erase(it5++);
              if(it4->element.empty()){
                break;
              }
            }
            else{
              it5++;
            }
          }
          if(it4->element.empty()){
            sortgraph.erase(it4++);
            break;
          }
        }
      }
    }

    //2.add edges for any two residues in the left bag
    if(newbag.left.size()>1){
      for(set<int>::iterator it3=newbag.left.begin();it3!=newbag.left.end();++it3){
        bool exist=false;
        for(int it5=0;it5<sortgraph.size();++it5){
          if(sortgraph[it5].idx==*it3){
            exist=true;
            break;
          }
        }
        if(exist==false){
          set <int> neibs;
          RankSize< set<int> > newedge(*it3,neibs);
          sortgraph.push_back(newedge);
        }
      }
      for(set<int>::iterator it3=newbag.left.begin();it3!=newbag.left.end();++it3){
        for(set<int>::iterator it4=it3;it4!=newbag.left.end();++it4){
          if(it4==it3) continue;
          for(int it5=0;it5<sortgraph.size();++it5){
            if(sortgraph[it5].idx==*it3){
              sortgraph[it5].element.insert(*it4);
              break;
            }
          }

          for(int it5=0;it5<sortgraph.size();++it5){
            if(sortgraph[it5].idx==*it4){
              sortgraph[it5].element.insert(*it3);
              break;
            }
          }
        }
      }
    }
    
    //3.resort the subgraph
    stable_sort(sortgraph.begin(), sortgraph.end(), less <RankSize<set<int> > >());
  }

  //show all of the bags on the tree
  //if(true){
  //  cout<<"@list bags (nodes) of the tree decomposition:"<<endl;
  //  for(int i=0; i<bags.size(); ++i){
  //    bags[i].ShowBag();
  //    cout<<endl;
  //  }
  //}

  //connect the bags into a tree
  set<int> vtsOnTree;
  int counter=0;
  while(bags.size()>0){
    Bag bg=bags[bags.size()-1];
    bags.pop_back();
    if(counter==0){
      //for the root bag, set the parent Bag Index as -1
      bg.parentBagIdx=-1;
      bg.type=Root;
      connBags.push_back(bg);
      //add vertices from bag to vtsOnTree;
      set_union(vtsOnTree.begin(),vtsOnTree.end(),
        bg.total.begin(),bg.total.end(),
        inserter(vtsOnTree,vtsOnTree.begin())
        );
    }
    else{
      set<int> intersection;
      set_intersection(vtsOnTree.begin(),vtsOnTree.end(),
        bg.total.begin(),bg.total.end(),
        inserter(intersection,intersection.begin())
        );
      for(int i=0; i<connBags.size();++i){
        Bag nbg=connBags[i];
        set<int> intersection2;
        set_intersection(nbg.total.begin(),nbg.total.end(),
          intersection.begin(),intersection.end(),
          inserter(intersection2,intersection2.begin())
          );
        if(intersection2==intersection){
          //bag nbg is the parent of the current bag
          //modify the attributes of parent bag
          bg.parentBagIdx=i;
          if(connBags[i].type!=Root){
            connBags[i].type=Inner;
          }
          connBags[i].childBagIdx.insert(connBags.size());
          connBags.push_back(bg);
          //add vertices from bag to vtsOnTree;
          set_union(vtsOnTree.begin(),vtsOnTree.end(),
            bg.total.begin(),bg.total.end(),
            inserter(vtsOnTree,vtsOnTree.begin())
            );
          break;
        }
      }
    }
    counter++;
  }
  //set the number of children for each bag
  for(int i=0; i<connBags.size(); ++i){
    connBags[i].childCounter=connBags[i].childBagIdx.size();
  }

  //show the connected bags, parent and children
  //if(false){
  //  //cout<<"@list connected bags of the tree decomposition:"<<endl;
  //  for(int i=0; i<connBags.size(); ++i){
  //    connBags[i].ShowBag();
  //    if(connBags[i].type==Root){
  //      //cout<<":"<<endl<<"index = "<<setw(4)<<i<<", type =  root, parent = none";
  //      //cout<<", children =";
  //      for(set<int>::iterator it1=connBags[i].childBagIdx.begin(); it1!=connBags[i].childBagIdx.end();++it1){
  //        cout<<" "<<*it1;
  //      }
  //      cout<<endl;
  //    }
  //    else if(connBags[i].type==Inner){
  //      cout<<":"<<endl<<"index = "<<setw(4)<<i<<", type = inner, parent = "<<setw(4)<<connBags[i].parentBagIdx;
  //      cout<<", children =";
  //      for(set<int>::iterator it1=connBags[i].childBagIdx.begin(); it1!=connBags[i].childBagIdx.end();++it1){
  //        cout<<" "<<*it1;
  //      }
  //      cout<<endl;
  //    }
  //    else{
  //      cout<<":"<<endl<<"index = "<<setw(4)<<i<<", type =  leaf, parent = "<<setw(4)<<connBags[i].parentBagIdx;
  //      cout<<", children = none"<<endl;
  //    }
  //  }
  //}

}

void TreeDecomposition::MergeBags(int depth)
{
  if(depth<connBags.size()){
    for(int i=0;i<connBags.size();i++){
      if(connBags[i].type==Leaf){
        //merge the leaf to its parent
        //cout<<"Merge leaf "<<i<<" with its parent "<<connBags[i].parentBagIdx<<": ";
        set_union(connBags[i].total.begin(),connBags[i].total.end(),
          connBags[connBags[i].parentBagIdx].total.begin(),connBags[connBags[i].parentBagIdx].total.end(),
          inserter(connBags[connBags[i].parentBagIdx].total,connBags[connBags[i].parentBagIdx].total.begin())
          );
        //for(set<int>::iterator it=connBags[connBags[i].parentBagIdx].total.begin();it!=connBags[connBags[i].parentBagIdx].total.end();++it){
        //  cout<<*it<<" ";
        //}
        //cout<<endl;
        connBags[i].type=None;
        connBags[connBags[i].parentBagIdx].childBagIdx.erase(i);
        if(connBags[connBags[i].parentBagIdx].childBagIdx.size()==0 && connBags[connBags[i].parentBagIdx].type != Root){
          connBags[connBags[i].parentBagIdx].type=Leaf;
        }
        depth++;
      }
    }
    MergeBags(depth);
  }
}

int TreeDecomposition::CheckTreewidth()
{
  int width=0;
  for(int i=0;i<connBags.size();i++){
    if(width<connBags[i].total.size()){
      width=connBags[i].total.size();
    }
  }
  width++;
  return width;
}


//# Solution Class

class Solution:public PairEnergy
{
public:
  ~Solution();
  IV1 unfixres;
  bool DEESearch(IV1 &pos);
  int DEEGoldstein(IV1& pos);
  int DEEsplit(IV1& pos);
  void Pick(int site,int rot);

  vector< Graph > graphs;
  void ConstructAdjMatrix(int nunfix,IV2 &adjMatrix);
  void ConstructSubgraphs(int nunfix,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix);
  void FindSubgraphByDFS(Graph &graph,int u,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix,stack<int> &vertices);
  //void ShowGraphs();
  void GraphEdgeDecomposition(IV2 &adjMatrix,float threshold);

  TreeDecomposition tree;
  void TreeDecompositionBottomToTopCalcEnergy();
  void CalcRightBagRotamerCombinationEnergy(Bag &leafbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Ercom,IV2 &Rrcom);
  void CalcLeftBagRotamerCombinationEnergy(Bag &rootbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Elcom,IV2 &Rlcom);
  void GetLeftBagRotamerCombination(Bag &leafbag,int depth,IV1 &Rtmp,IV2 &Rlcom);
  void BagDeploySites(Bag &leafbag);
  void LeafBagCalcEnergy(Bag &leafbag,IV2 &Rlcom);
  void CombineChildIntoParentBag(Bag &leafbag,Bag &parbag,IV2 &Rclcom);
  void SubsetCheck(IV1 &subset,IV1 &fullset, IV1 &indices);
  void RootBagFindGMEC(Bag &rootbag);
  void TreeDecompositionTopToBottomAssignRotamer(Bag &parbag,Bag &childbag);
  void TreeDecompositionRelease();
  void Search();
};


Solution::~Solution()
{
  unfixres.clear();
}


bool Solution::DEESearch(IV1 &pos)
{
  //DEEGoldstein
  IV1 fixres;
  int ndeadDEE=1,iterDEE=1;
  unfixres.clear();
  while(ndeadDEE!=0){
    ndeadDEE=DEEGoldstein(pos);
    //cout<<"iter "<<iterDEE<<" DEEgoldstein eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  //cout<<"#residues fixed after DEE-Goldstein: "<<fixres.size()<<endl;
  //cout<<"#residues unfixed after DEE-Goldstein: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  //DEEsplit
  pos.clear();
  fixres.clear();
  unfixres.clear();
  ndeadDEE=1,iterDEE=1;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  while(ndeadDEE!=0){
    ndeadDEE=DEEsplit(pos);
    //cout<<"iter "<<iterDEE<<" DEEsplit eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  //cout<<"#residues fixed after DEE-split: "<<fixres.size()<<endl;
  //cout<<"#residues unfixed after DEE-split: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  //DEEGoldstein
  pos.clear();
  fixres.clear();
  unfixres.clear();
  ndeadDEE=1,iterDEE=1;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  while(ndeadDEE!=0){
    ndeadDEE=DEEGoldstein(pos);
    //cout<<"iter "<<iterDEE<<" DEEgoldstein eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  //cout<<"#residues fixed after DEE-Goldstein: "<<fixres.size()<<endl;
  //cout<<"#residues unfixed after DEE-Goldstein: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  //DEEsplit
  pos.clear();
  fixres.clear();
  unfixres.clear();
  ndeadDEE=1,iterDEE=1;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  while(ndeadDEE!=0){
    ndeadDEE=DEEsplit(pos);
    //cout<<"iter "<<iterDEE<<" DEEsplit eliminates "<<ndeadDEE<<" rotamers"<<endl;
    iterDEE++;
  };
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    int n=0;
    int rot;
    for(int j=0;j<nrots[ip];j++){
      if(eTableSelf[ip][j]<999.){
        rot=j;
        n++;
      }
    }
    if(n==1){
      nrots[ip]=1;
      bestrot[ip]=rot;
      fixres.push_back(ip);
    }
    else if(n>1){
      unfixres.push_back(ip);
    }
  }
  //cout<<"#residues fixed after DEE-split: "<<fixres.size()<<endl;
  //cout<<"#residues unfixed after DEE-split: "<<unfixres.size()<<endl;
  if(unfixres.size()==0) return false;
  for(int i=0;i<unfixres.size();i++){
    int ipos=unfixres[i];
    for(int j=0;j<fixres.size();j++){
      int jpos=fixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      int rot=bestrot[jpos];
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]<999.){
          eTableSelf[ipos][k] += eTablePair[ipos][jpos][k][rot];
        }
      }
    }
  }

  return true;
}

int Solution::DEEGoldstein(IV1& pos)
{
  int elimination=0;
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    for(int s=0;s<nrots[ip];s++){
      if(eTableSelf[ip][s]>999.) continue;
      for(int r=0;r<nrots[ip];r++){
        if(r==s) continue;
        else if(eTableSelf[ip][r]>999.) continue;
        float ex=eTableSelf[ip][s]-eTableSelf[ip][r];
        for(int j=0;j<pos.size();j++){
          if(j==i) continue;
          int jp=pos[j];
          if(eTablePair[ip][jp]==NULL) continue;
          float ey=1e8;
          for(int t=0;t<nrots[jp];t++){
            if(eTableSelf[jp][t]>999.) continue;
            float em=eTablePair[ip][jp][s][t]-eTablePair[ip][jp][r][t];
            if(em<ey) ey=em;
          }
          ex+=ey;
        }
        if(ex>DEE_THRESHOLD){
          eTableSelf[ip][s]=1000.;
          elimination++;
          break;
        }
      }
    }
  }
  return elimination;
}

int Solution::DEEsplit(IV1& pos)
{
  int elimination=0;
  for(int i=0;i<pos.size();i++){
    int ip=pos[i];
    for(int s=0;s<nrots[ip];s++){
      if(eTableSelf[ip][s]>999.) continue;
      FV1 storeYj(pos.size(),0.);
      FV2 storeYjr(nrots[ip],storeYj);
      for(int r=0;r<nrots[ip];r++){
        if(r==s) continue;
        else if(eTableSelf[ip][r]>999.) continue;
        for(int j=0;j<pos.size();j++){
          if(j==i) continue;
          int jp=pos[j];
          if(eTablePair[ip][jp]==NULL) continue;
          float ey=1e8;
          for(int t=0;t<nrots[jp];t++){
            if(eTableSelf[jp][t]>999.) continue;
            float em=eTablePair[ip][jp][s][t]-eTablePair[ip][jp][r][t];
            if(em<ey) ey=em;
          }
          storeYjr[r][j]=ey;
        }
      }

      IV1 elim;
      for(int k=0;k<pos.size();k++){
        if(k==i) continue;
        int kp=pos[k];
        if(eTablePair[ip][kp]==NULL) continue;
        for(int v=0;v<nrots[kp];v++){
          elim.push_back(0);
        }
        for(int r=0;r<nrots[ip];r++){
          if(r==s) continue;
          else if(eTableSelf[ip][r]>999.) continue;

          float ex=eTableSelf[ip][s]-eTableSelf[ip][r];
          for(int j=0;j<pos.size();j++){
            if(j==i || j==k) continue;
            int jp=pos[j];
            if(eTablePair[ip][jp]==NULL) continue;
            ex+=storeYjr[r][j];
          }
          for(int v=0;v<nrots[kp];v++){
            if(eTableSelf[kp][v]>999.) continue;
            if(ex+eTablePair[ip][kp][s][v]-eTablePair[ip][kp][r][v]>DEE_THRESHOLD){
              elim[v]=1;
            }
          }
        }
        bool allelim=true;
        for(int v=0;v<nrots[kp];v++){
          if(elim[v]==0){
            allelim=false;
            break;
          }
        }
        if(allelim==true){
          eTableSelf[ip][s]=1000.;
          elimination++;
          goto FLAG_SPLIT;
        }
      }
FLAG_SPLIT:;
    }
  }
  return elimination;
}

void Solution::Pick(int site,int rot)
{
  for(int i=5;i<stru[site].atTypes.size();i++){
    stru[site].xyz.push_back(sc[site][rot][i-5]);
  }
}

void Solution::ConstructAdjMatrix(int nunfix,IV2 &adjMatrix)
{
  //cout<<endl<<"construct adjacent matrix ..."<<endl;
  //cout<<"remove edges (residue-residue interactions) with small energy values ... ";
  int n_removed_edges=0;
  for(int i=0;i<nunfix-1;i++){
    int ipos=unfixres[i];
    for(int j=i+1;j<nunfix;j++){
      int jpos=unfixres[j];
      if(eTablePair[ipos][jpos]==NULL)continue;
      for(int k=0;k<nrots[ipos];k++){
        if(eTableSelf[ipos][k]>999.)continue;
        for(int l=0;l<nrots[jpos];l++){
          if(eTableSelf[jpos][l]>999.)continue;
          if(eTablePair[ipos][jpos][k][l]>EPAIR_CUT || eTablePair[ipos][jpos][k][l]<-1.*EPAIR_CUT){
            adjMatrix[i][j]=1;
            adjMatrix[j][i]=1;
            goto FLAG1;
          }
        }
      }
      eTablePair[ipos][jpos]=NULL;
      //cout<<"residue "<<ipos<<" and "<<jpos<<" does not interact"<<endl;
      n_removed_edges++;
FLAG1: continue;
    }
  }
  //cout<<"#edges removed: "<<n_removed_edges<<endl;

  //cout<<"remove the residues that have no interaction edge ... ";
  int n_res_noedge=0;
  for(int i=0;i<nunfix;i++){
    bool allzeros=true;
    for(int j=0;j<nunfix;j++){
      if(adjMatrix[i][j]==1){
        allzeros=false;
        break;
      }
    }
    if(allzeros==true){
      float emin=1e8;
      int rot;
      for(int r=0;r<nrots[unfixres[i]];r++){
        if(eTableSelf[unfixres[i]][r]<emin){
          emin=eTableSelf[unfixres[i]][r];
          rot=r;
        }
      }
      nrots[unfixres[i]]=1;
      bestrot[unfixres[i]]=rot;
      n_res_noedge++;
    }
  }
  //cout<<"#residues without edges removed: "<<n_res_noedge<<endl;
}

void Solution::ConstructSubgraphs(int nunfix,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix)
{
  visited.assign(nunfix,0);
  for(int i=0;i<nunfix;i++){
    if(!visited[i]){
      Graph newsg;
      stack<int> vertices;
      FindSubgraphByDFS(newsg,i,visited,adjMatrix,flagMatrix,vertices);
      if(!newsg.empty()){
        graphs.push_back(newsg);
      }

      if(newsg.size()==1){
        int site=unfixres[i];
        int rot=0;
        float eMin=1000.;
        for(int j=0;j<nrots[site];j++){
          if(eTableSelf[site][j]>999.)continue;
          if(eTableSelf[site][j]<eMin){
            eMin=eTableSelf[site][j];
            rot=j;
          }
        }
        bestrot[site]=rot;
        nrots[site]=1;
      }
    }
  }

  for(vector<int>::iterator it=unfixres.begin();it!=unfixres.end();){
    if(nrots[*it]==1){
      it=unfixres.erase(it);
    }
    else{
      ++it;
    }
  }


}

void Solution::FindSubgraphByDFS(Graph &graph,int u,IV1 &visited,IV2 &adjMatrix,IV2 &flagMatrix,stack<int> &vertices)
{
  visited[u]=1;
  vertices.push(u);
  while(!vertices.empty()){
    bool hasEdge=false;
    int w=-1;
    for(w=0;w<adjMatrix[u].size();++w){
      if(adjMatrix[u][w]==1 && flagMatrix[u][w]==0){
        hasEdge=true;
        break;
      }
    }

    if(hasEdge==false){
      u=vertices.top();
      vertices.pop();
    }
    else{
      visited[w]=1;
      vertices.push(w);
      flagMatrix[u][w]=1;
      flagMatrix[w][u]=1;
      graph[unfixres[u]].insert(unfixres[w]);
      graph[unfixres[w]].insert(unfixres[u]);
    }
  }
}

//void Solution::ShowGraphs()
//{
//  for(int i1=0; i1<graphs.size(); i1++){
//    //cout<<"*subgraph "<<i1<<":"<<endl;
//    Graph &graph=graphs[i1];
//    ShowGraph(graph);
//  }
//}

void Solution::GetLeftBagRotamerCombination(Bag &leafbag,int depth,IV1 &Rtmp,IV2 &Rlcom)
{
  if(depth<leafbag.lsites.size()){
    for(int i=0;i<leafbag.lrots[depth].size();i++){
      int roti=leafbag.lrots[depth][i];
      Rtmp.push_back(roti);
      GetLeftBagRotamerCombination(leafbag,depth+1,Rtmp,Rlcom);
      Rtmp.pop_back();
    }
  }
  else{
    Rlcom.push_back(Rtmp);
  }
}

void Solution::CalcLeftBagRotamerCombinationEnergy(Bag &rootbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Elcom,IV2 &Rlcom)
{
  if(depth<rootbag.lsites.size()){
    int site=rootbag.lsites[depth];
    for(int i=0;i<rootbag.lrots[depth].size();i++){
      int roti=rootbag.lrots[depth][i];
      float Eold=Etmp;
      float Enew=eTableSelf[site][roti];
      for(int k=0;k<depth;k++){
        if(eTablePair[rootbag.lsites[k]][site]==NULL) continue;
        Enew+=eTablePair[rootbag.lsites[k]][site][Rlcom[Rlcom.size()-1][k]][roti];
      }
      Etmp+=Enew;
      Rtmp.push_back(roti);
      CalcLeftBagRotamerCombinationEnergy(rootbag,depth+1,Etmp,Rtmp,Elcom,Rlcom);
      Rtmp.pop_back();
      Etmp=Eold;
    }
  }
  else{
    Rlcom.push_back(Rtmp);
    Elcom.push_back(Etmp);
  }
}

void Solution::CalcRightBagRotamerCombinationEnergy(Bag &leafbag,int depth,float &Etmp,IV1 &Rtmp,FV1 &Ercom,IV2 &Rrcom)
{
  if(depth<leafbag.rsites.size()){
    int site=leafbag.rsites[depth];
    for(int i=0;i<leafbag.rrots[depth].size();i++){
      int roti=leafbag.rrots[depth][i];
      float Eold=Etmp;
      float Enew=eTableSelf[site][roti];
      for(int k=0;k<depth;k++){
        if(eTablePair[leafbag.rsites[k]][site]==NULL) continue;
        Enew+=eTablePair[leafbag.rsites[k]][site][Rrcom[Rrcom.size()-1][k]][roti];
      }
      Etmp+=Enew;
      Rtmp.push_back(roti);
      CalcRightBagRotamerCombinationEnergy(leafbag,depth+1,Etmp,Rtmp,Ercom,Rrcom);
      Rtmp.pop_back();
      Etmp=Eold;
    }
  }
  else{
    Rrcom.push_back(Rtmp);
    Ercom.push_back(Etmp);
  }
}

void Solution::BagDeploySites(Bag &leafbag)
{
  if(leafbag.deployFlag==false){
    leafbag.deployFlag=true;
    for(set<int>::iterator it=leafbag.left.begin();it!=leafbag.left.end();it++){
      leafbag.lsites.push_back(*it);
      leafbag.tsites.push_back(*it);
    }
    for(set<int>::iterator it=leafbag.right.begin();it!=leafbag.right.end();it++){
      leafbag.rsites.push_back(*it);
      leafbag.tsites.push_back(*it);
    }

    IV1 rots;
    for(int idx2=0;idx2<leafbag.lsites.size();idx2++){
      int site=leafbag.lsites[idx2];
      for(int j=0;j<nrots[site];j++){
        if(eTableSelf[site][j]>999.) continue;
        rots.push_back(j);
      }
      leafbag.lrots.push_back(rots);
      rots.clear();
    }
    for(int idx2=0;idx2<leafbag.rsites.size();idx2++){
      int site=leafbag.rsites[idx2];
      for(int j=0;j<nrots[site];j++){
        if(eTableSelf[site][j]>999.) continue;
        rots.push_back(j);
      }
      leafbag.rrots.push_back(rots);
      rots.clear();
    }

    int lcount=1;
    int counter=0;
    while(counter<leafbag.lrots.size()){
      lcount *=leafbag.lrots[counter].size();
      counter++;
    }

    int rcount=1;
    counter=0;
    while(counter<leafbag.rrots.size()){
      rcount *= leafbag.rrots[counter].size();
      counter++;
    }
    for(int i=0;i<lcount;i++){
      IV1 Rtmp(rcount,-1);
      FV1 Etmp(rcount,0.);
      leafbag.Etcom.push_back(Etmp);
      leafbag.Rtcom.push_back(Rtmp);
    }
  }
}

void Solution::LeafBagCalcEnergy(Bag &leafbag,IV2 &Rlcom)
{
  FV1 Ercom;
  float Ertmp=0.;
  IV2 Rrcom;
  IV1 Rrtmp, Rltmp;
  CalcRightBagRotamerCombinationEnergy(leafbag,0,Ertmp,Rrtmp,Ercom,Rrcom);
  GetLeftBagRotamerCombination(leafbag,0,Rltmp,Rlcom);
  for(int idx2=0;idx2<Rlcom.size();idx2++){
    float emin=1e8;
    IV1 Rrmin;
    for(int idx3=0;idx3<Rrcom.size();idx3++){
      float eval=0.;
      for(int idx4=0;idx4<Rlcom[idx2].size();idx4++){
        for(int idx5=0;idx5<Rrcom[idx3].size();idx5++){
          if(eTablePair[leafbag.lsites[idx4]][leafbag.rsites[idx5]] != NULL){
            eval+=eTablePair[leafbag.lsites[idx4]][leafbag.rsites[idx5]][Rlcom[idx2][idx4]][Rrcom[idx3][idx5]];
          }
        }
      }
      eval+=Ercom[idx3];
      if(leafbag.Etcom.size() != 0){
        //add the pre-stored energy
        eval+=leafbag.Etcom[idx2][idx3];
      }
      if(eval<emin){
        Rrmin=Rrcom[idx3];
        emin=eval;
      }
    }
    leafbag.Etcom[idx2][0]=emin;

    //record the best rotamer combinations
    IV1 Rttmp;
    for(int idx4=0;idx4<Rlcom[idx2].size();idx4++){
      Rttmp.push_back(Rlcom[idx2][idx4]);
    }
    for(int idx5=0;idx5<Rrmin.size();idx5++){
      Rttmp.push_back(Rrmin[idx5]);
    }
    leafbag.Rtcom[idx2]=Rttmp;
    Rttmp.clear();
  }

  Ercom.clear();
  Rrcom.clear();
  Rrtmp.clear();
  Rltmp.clear();
}

void Solution::CombineChildIntoParentBag(Bag &leafbag,Bag &parbag,IV2 &Rclcom)
{
  FV1 Ercom;
  float Ertmp=0.;
  IV2 Rrcom, Rlcom;
  IV1 Rrtmp, Rltmp;
  CalcRightBagRotamerCombinationEnergy(parbag,0,Ertmp,Rrtmp,Ercom,Rrcom);
  GetLeftBagRotamerCombination(parbag,0,Rltmp,Rlcom);
  SubsetCheck(leafbag.lsites,parbag.tsites,leafbag.indices);
  for(int idx2=0;idx2<Rlcom.size();idx2++){
    for(int idx3=0;idx3<Rrcom.size();idx3++){
      IV1 ppartrot,ppartsite;
      for(int j=0;j<leafbag.indices.size();j++){
        if(leafbag.indices[j]>=parbag.lsites.size()){
          ppartrot.push_back(Rrcom[idx3][leafbag.indices[j]-parbag.lsites.size()]);
          ppartsite.push_back(parbag.rsites[leafbag.indices[j]-parbag.lsites.size()]);
        }
        else{
          ppartrot.push_back(Rlcom[idx2][leafbag.indices[j]]);
          ppartsite.push_back(parbag.lsites[leafbag.indices[j]]);
        }
      }

      for(int j=0;j<Rclcom.size();j++){
        if(ppartsite==leafbag.lsites && ppartrot==Rclcom[j]){
          parbag.Etcom[idx2][idx3]+=leafbag.Etcom[j][0];
        }
      }
    }
  }
}

void Solution::SubsetCheck(IV1 &leaflsites,IV1 &partsites, IV1 &indices)
{
  indices.clear();
  for(int i=0;i<leaflsites.size();i++){
    indices.push_back(-1);
    for(int j=0;j<partsites.size();j++){
      if(leaflsites[i]==partsites[j]){
        indices[i]=j;
      }
    }
  }
}

void Solution::RootBagFindGMEC(Bag &rootbag)
{
  FV1 Ercom,Elcom;
  float Ertmp=0., Eltmp=0.;
  IV2 Rrcom, Rlcom;
  IV1 Rrtmp, Rltmp;
  CalcRightBagRotamerCombinationEnergy(rootbag,0,Ertmp,Rrtmp,Ercom,Rrcom);
  CalcLeftBagRotamerCombinationEnergy(rootbag,0,Eltmp,Rltmp,Elcom,Rlcom);

  rootbag.EGMEC=1e8;
  for(int idx2=0;idx2<Rlcom.size();idx2++){
    for(int idx3=0;idx3<Rrcom.size();idx3++){

      float energy=Elcom[idx2];
      for(int idx4=0;idx4<Rlcom[idx2].size();idx4++){
        for(int idx5=0;idx5<Rrcom[idx3].size();idx5++){
          if(eTablePair[rootbag.lsites[idx4]][rootbag.rsites[idx5]] != NULL){
            energy+=eTablePair[rootbag.lsites[idx4]][rootbag.rsites[idx5]][Rlcom[idx2][idx4]][Rrcom[idx3][idx5]];
          }
        }
      }
      energy+=Ercom[idx3];
      if(rootbag.Etcom.size() != 0){
        energy+=rootbag.Etcom[idx2][idx3];
      }

      if(energy<rootbag.EGMEC){
        rootbag.EGMEC=energy;
        rootbag.RLGMEC.assign(Rlcom[idx2].begin(),Rlcom[idx2].end());
        rootbag.RRGMEC.assign(Rrcom[idx3].begin(),Rrcom[idx3].end());
      }

    }
  }

  //set optimal rotamer index
  for(int i=0;i<rootbag.lsites.size();i++){
    nrots[rootbag.lsites[i]]=1;
    bestrot[rootbag.lsites[i]]=rootbag.RLGMEC[i];
  }
  for(int i=0;i<rootbag.rsites.size();i++){
    nrots[rootbag.rsites[i]]=1;
    bestrot[rootbag.rsites[i]]=rootbag.RRGMEC[i];
  }

}

void Solution::TreeDecompositionBottomToTopCalcEnergy()
{
  if(tree.connBags.size()==1){
    BagDeploySites(tree.connBags[0]);
  }
  else{
    while(true){
      int leafCount=0;
      for(int idx1=0;idx1<tree.connBags.size();idx1++){
        if(tree.connBags[idx1].type==Leaf){
          IV2 Rclcom;
          BagDeploySites(tree.connBags[idx1]);
          LeafBagCalcEnergy(tree.connBags[idx1],Rclcom);

          BagDeploySites(tree.connBags[tree.connBags[idx1].parentBagIdx]);
          CombineChildIntoParentBag(tree.connBags[idx1],tree.connBags[tree.connBags[idx1].parentBagIdx],Rclcom);

          tree.connBags[idx1].type=None;
          tree.connBags[tree.connBags[idx1].parentBagIdx].childCounter--;
          if(tree.connBags[tree.connBags[idx1].parentBagIdx].childCounter==0 && 
            tree.connBags[tree.connBags[idx1].parentBagIdx].type != Root){
              tree.connBags[tree.connBags[idx1].parentBagIdx].type=Leaf;
          }
          leafCount++;
        }
      }
      if(leafCount==0){
        break;
      }
    }
  }

}


void Solution::TreeDecompositionTopToBottomAssignRotamer(Bag &parbag,Bag &childbag)
{
  bool find=false;
  IV1 ppartrot,ppartsite;
  for(int j=0;j<childbag.indices.size();j++){
    if(childbag.indices[j]>=parbag.lsites.size()){
      ppartrot.push_back(parbag.RRGMEC[childbag.indices[j]-parbag.lsites.size()]);
      ppartsite.push_back(parbag.rsites[childbag.indices[j]-parbag.lsites.size()]);
    }
    else{
      ppartrot.push_back(parbag.RLGMEC[childbag.indices[j]]);
      ppartsite.push_back(parbag.lsites[childbag.indices[j]]);
    }
  }

  for(int j=0;j<childbag.Rtcom.size();j++){
    IV1 Rcltmp,Rcrtmp;
    Rcltmp.assign(childbag.Rtcom[j].begin(),childbag.Rtcom[j].begin()+childbag.lsites.size());
    Rcrtmp.assign(childbag.Rtcom[j].begin()+childbag.lsites.size(),childbag.Rtcom[j].end());
    if(ppartsite==childbag.lsites && ppartrot==Rcltmp){
      childbag.RLGMEC=Rcltmp;
      childbag.RRGMEC=Rcrtmp;
      for(int p=0;p<childbag.rsites.size();p++){
        nrots[childbag.rsites[p]]=1;
        bestrot[childbag.rsites[p]]=childbag.RRGMEC[p];
      }
      find=true;
      break;
    }
  }

  for(set<int>::iterator it=childbag.childBagIdx.begin();it!=childbag.childBagIdx.end();it++){
    TreeDecompositionTopToBottomAssignRotamer(childbag,tree.connBags[*it]);
  }
}


void Solution::TreeDecompositionRelease()
{
  tree.bags.clear();
  tree.connBags.clear();
}


void Solution::GraphEdgeDecomposition(IV2 &adjMatrix,float threshold)
{
  int nedgeremoved=0;
  for(int i=0;i<adjMatrix.size()-1;i++){
    int k=unfixres[i];
    for(int j=i+1;j<adjMatrix.size();j++){
      int l=unfixres[j];
      if(adjMatrix[i][j]==1 && eTablePair[k][l]!=NULL){
        int counterm=0, countern=0;
        for(int m=0;m<nrots[k];m++){
          if(eTableSelf[k][m]>999.0) continue;
          counterm++;
        }
        for(int n=0;n<nrots[l];n++){
          if(eTableSelf[l][n]>999.0) continue;
          countern++;
        }

        float abar=0.;
        for(int m=0;m<nrots[k];m++){
          if(eTableSelf[k][m]>999.0) continue;
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.0) continue;
            abar+=eTablePair[k][l][m][n];
          }
        }
        abar/=(2.*counterm*countern);

        FV1 ak,bl;
        ak.assign(nrots[k],1000.);
        bl.assign(nrots[l],1000.);
        for(int m=0;m<nrots[k];m++){
          float temp=0.;
          if(eTableSelf[k][m]>999.0) continue;
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.0) continue;
            temp+=eTablePair[k][l][m][n];
          }
          temp/=countern;
          ak[m]=temp-abar;
        }

        for(int n=0;n<nrots[l];n++){
          float temp=0.;
          if(eTableSelf[l][n]>999.0) continue;
          for(int m=0;m<nrots[k];m++){
            if(eTableSelf[k][m]>999.0) continue;
            temp+=eTablePair[k][l][m][n];
          }
          temp/=counterm;
          bl[n]=temp-abar;
        }

        //estimate max deviation
        float maxdev=-1e8;
        for(int m=0;m<nrots[k];m++){
          if(eTableSelf[k][m]>999.0) continue;
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.0) continue;
            float dev=abs(eTablePair[k][l][m][n]-ak[m]-bl[n]);
            if(dev>maxdev){
              maxdev=dev;
            }
          }
        }

        if(maxdev<=threshold){
          adjMatrix[i][j]=adjMatrix[j][i]=0;
          for(int m=0;m<nrots[k];m++){
            if(eTableSelf[k][m]>999.) continue;
            eTableSelf[k][m]+=ak[m];
          }
          for(int n=0;n<nrots[l];n++){
            if(eTableSelf[l][n]>999.) continue;
            eTableSelf[l][n]+=bl[n];
          }
          //cout<<"edge between site "<<k<<" and "<<l<<" is decomposed with threshold "<<threshold<<endl;
          nedgeremoved++;
        }
      }
    }
  }
  //cout<<"edge decomposition with threshold "<<threshold<<"... #edges removed: "<<nedgeremoved<<endl;

  int resremoved=0;
  for(int i=0;i<adjMatrix.size();i++){
    int k=unfixres[i];
    int nonzero=0;
    for(int j=0;j<adjMatrix.size();j++){
      if(adjMatrix[i][j]!=0){
        nonzero++;
        break;
      }
    }
    if(nonzero==0){
      int best;
      float emin=1e8;
      for(int j=0;j<nrots[k];j++){
        if(eTableSelf[k][j]>999.0) continue;
        if(eTableSelf[k][j]<emin){
          emin=eTableSelf[k][j];
          best=j;
        }
      }
      bestrot[k]=best;
      nrots[k]=1;
      //cout<<"residue "<<k<<" eliminated by edge decomposition"<<endl;
      resremoved++;
    }
  }
  //cout<<"remove residues with no edge ... "<<"#residues removed: "<<resremoved<<endl;
}


void Solution::Search()
{
  IV1 pos;
  for(int i=0;i<nres;i++){
    if(nrots[i]>1) pos.push_back(i);
  }
  if(DEESearch(pos)==true){
    bool hardmode=false;
    float threshold=0.5;
LOOP:
    int nunfix=unfixres.size();
    IV1 visited(nunfix,0);
    IV2 adjMatrix(nunfix,visited);
    IV2 flagMatrix(nunfix,visited);
    //construct subgraphs
    ConstructAdjMatrix(nunfix,adjMatrix);
    if(hardmode==true){
      hardmode=false;
      GraphEdgeDecomposition(adjMatrix,threshold);
      threshold*=2;
    }
    ConstructSubgraphs(nunfix,visited,adjMatrix,flagMatrix);
    //ShowGraphs();

    for(int i1=0; i1<graphs.size(); i1++){
      Graph &graph=graphs[i1];
      tree.Subgraph2TreeDecomposition(i1,graph);
      int treewidth=tree.CheckTreewidth();
      //cout<<"#treewidth = "<<treewidth<<endl;
      if(treewidth<=TREEWIDTH_CUT){
        TreeDecompositionBottomToTopCalcEnergy();
        RootBagFindGMEC(tree.connBags[0]);
        for(set<int>::iterator it=tree.connBags[0].childBagIdx.begin();it!=tree.connBags[0].childBagIdx.end();it++){
          TreeDecompositionTopToBottomAssignRotamer(tree.connBags[0],tree.connBags[*it]);
        }
        for(Graph::iterator it2=graph.begin();it2!=graph.end();it2++){
          int site=it2->first;
          for(vector<int>::iterator it=unfixres.begin();it!=unfixres.end();it++){
            if(*it==site){
              it=unfixres.erase(it);
              break;
            }
          }
        }
        //cout<<"current tree has been solved"<<endl;
      }
      else{
        hardmode=true;
        //cout<<"current tree is hard, to be solved later"<<endl;
      }
      TreeDecompositionRelease();
      graph.clear();
    }
    graphs.clear();
    if(hardmode) goto LOOP;
  }

  //cout<<"optimal rotamer index:"<<endl;
  for(int i=0;i<nres;i++){
    if(nrots[i]==1){
      Pick(i,bestrot[i]);
      //cout<<i<<"\t"<<optRotIdx[i]<<endl;
    }
  }

  pdb=stru;
  stru.clear();
  sc.clear();
}


//# Adpator functions
// adaptor functions were created by @joaomcteixeira

void extract_coordinates(vector<double> &coords, PV1 &pdb, int nres);
void inject_coords(const vector<double> &coords, const char *labels, PV1 &pdb);

void inject_coords(const vector<double> &coords, const char *labels, PV1 &pdb)
{
    /*
    Injects input coordinates and labels to Structure.pdb in the form of
    Residues.

    This function replaces, and is based upon, the original Structure::ReadPDB.

    Author: @joaomcteixeira, Feb 2021
    */

    Residue residue;
    FV1 xyz(3);

    for(int i = 0; i < strlen(labels); i++)
    {
        residue.name = One2Three(labels[i]);
        residue.chID = 'A';
        residue.pos = i + 1;
        residue.ins = ' ';
        residue.atNames.push_back(" N  ");
        residue.atNames.push_back(" CA ");
        residue.atNames.push_back(" C  ");
        residue.atNames.push_back(" O  ");
        residue.atTypes = "NCCO";

        for(int atom = 0; atom < 4; atom++)
        {
            // NOTE: residue.xyz is FV2
            xyz.push_back(coords[i * 12 + 0 + 3 * atom]);
            xyz.push_back(coords[i * 12 + 1 + 3 * atom]);
            xyz.push_back(coords[i * 12 + 2 + 3 * atom]);
            residue.xyz.push_back(xyz);
            xyz.clear();
        }

        pdb.push_back(residue);  // sends residue to protein
        residue.xyz.clear();
        residue.atNames.clear();
        residue.atTypes.clear();
    }
}


void extract_coordinates(vector<double> &coords, PV1 &pdb, int nres)
{
    /*
    Extracts coordinates from Structure.pdb.residue(s) to vector of doubles.

    Author: @joaomcteixeira, Feb 2021
    */
    for(int residue = 0; residue < nres; residue++)
    {
        for(int atom = 0; atom < pdb[residue].atNames.size(); atom++)
        {
            for (int i = 0; i < 3; i++)
            {
                coords.push_back(pdb[residue].xyz[atom][i]);
            }
        }
    }
}


vector<double> main_faspr(const vector<double>& xyz, const char *residue_labels)
{
    /*

    Reproduces the original workflow of FASPR.
    Code related to PDB input and output was removed from the original.

    Author: @joaomcteixeira, Feb 2021
    */
    Solution faspr;
    inject_coords(xyz, residue_labels, faspr.pdb);
    faspr.nres = faspr.pdb.size(); // avoids sending faspr to inject_coords
    faspr.LoadSeq();
    faspr.BuildSidechain();
    faspr.CalcSelfEnergy();
    faspr.CalcPairEnergy();
    faspr.Search();

    // transfers coordinates in Structure to a vector
    // this vector will be returned to Python
    vector<double> coordinates;
    extract_coordinates(coordinates, faspr.pdb, faspr.nres);
    return coordinates;
}


// beautiful pybind11 :-)
py::array faspr_sidechains( \
        py::array_t<double, py::array::c_style | py::array::forcecast> conf_xyz, \
        const char *reslabel, \
        const string &dunbbdepfile)
{
    /*
    Links Python with this c++ module.

    Parameters
    ----------
    conf_xyz : (N, 3) float array
        The coordinates of ONLY N, CA, C, O backbone atoms. Atoms must
        be sorted in the indicated fashion and residues must be sorted.
        Failure to comply may result in unexpected errors or
        Segmentation Fault. For efficiency, several input checks
        performed originally by FASPR were desactivate here.

    reslabel : UTF-8 Python string
        The 1-letter residue code sequence of the protein.

    dunbbdepfile : UTF-8 Python string
        The absolute path to the dun2010bbdeb.bin file.

    written by @joaomcteixeira
    */

    // The original FASPR reads the dun2010bbdeb.bin file on every
    // execution. In this module, the bin file is kept open as a global
    // variable and every call to the module restarts the file position.
    if(infile.is_open())
    {
        infile.clear();
        infile.seekg(0);
    }
    else  // in the first call, the file is openned
    {
        infile.open(dunbbdepfile.c_str(), ifstream::in|ifstream::binary);
    }

    // allocate std::vector (to pass to the C++ function)
    std::vector<double> xyz(conf_xyz.size());

    // copy py::array -> std::vector
    std::memcpy(xyz.data(), conf_xyz.data(), conf_xyz.size() * sizeof(double));

    // call pure C++ function
    std::vector<double> coords_with_sidechains = main_faspr(xyz, reslabel);

    ssize_t              ndim    = 2;
    const int            scsize  = coords_with_sidechains.size();
    std::vector<ssize_t> shape   = {scsize / 3, 3};
    std::vector<ssize_t> strides = {sizeof(double) * 3 ,sizeof(double)};

    // return 2-D NumPy array
    return py::array(py::buffer_info(
    coords_with_sidechains.data(),           /* data as contiguous array  */
    sizeof(double),                          /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    ndim,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
    ));
}

// wrap as Python module
// this must have the name of the .cpp file also
// @joaomcteixeira
PYBIND11_MODULE(idpcpp, m)
{
  m.doc() = "CPP intensive executions for IDPConfGen.";

  m.def( \
    "faspr_sidechains", \
    &faspr_sidechains, \
    "Compute sidechains using FASPR methodology.");
}
