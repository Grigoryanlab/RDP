// a function that creates an ensemble of randomly docked structures, optionally biased towards certain binding residues for either binding partner, and outputs a distribution of their RMSDs when compared to the correct docked conformation

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include <cmath>
/*#include <time.h>
#include <chrono>
#include <sys/time.h>
#include <ctime>*/
#include <iostream>
#include "msttypes.h"
#include "mstoptions.h"
#include "mstrotlib.h"
#include "msttransforms.h"
#include "mstcondeg.h"
#include "mstfasst.h"
#include "mstsequence.h"

using std::cout; using std::endl;
/*using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::system_clock;*/

using namespace std;
using namespace MST;

void lilWriteTest(string fileName, string toWrite) {
  fstream myFile;
  myFile.open(fileName, ios::out);
  if( !myFile ) {
     cerr << "Error: file could not be opened" << endl;
     exit(1);
  }
  myFile << toWrite << endl;
  myFile.close();
}

void printMapIntVecInt(map<int,vector<int>> myMap) {
  for(auto it = myMap.begin(); it != myMap.end(); ++it) {
    std::cout << "key: " << it->first << endl;
    vector<int> valuesVect = myMap[it->first];
    for (int i = 0; i < valuesVect.size(); i++) {
      cout << "value: " << valuesVect[i] << endl;
    }
  }
}

AtomPointerVector getBBatoms(AtomPointerVector startingAtoms) {
  AtomPointerVector bbAtoms;
  for (int i = 0; i < startingAtoms.size(); i++) {
    Atom* currAtom = startingAtoms[i];
    string cAtomName = currAtom->getName();
    if ((cAtomName == "CA") or (cAtomName == "N") or (cAtomName == "O") or (cAtomName == "C")) {
      bbAtoms.push_back(currAtom);
    }
  }
  return bbAtoms;
}

AtomPointerVector getHeavyAtoms(AtomPointerVector startingAtoms) {
  AtomPointerVector heavyAtoms;
  for (int i = 0; i < startingAtoms.size(); i++) {
    Atom* currAtom = startingAtoms[i];
    string cAtomName = currAtom->getName();
    char firstChar = cAtomName[0];
    if ((firstChar != 'H')) {
      heavyAtoms.push_back(currAtom);
    }
  }
  return heavyAtoms;
}

vector<vector<int>> getABindexes(AtomPointerVector aAtoms,AtomPointerVector bAtoms,mstreal conDist) {
  ProximitySearch psA(aAtoms, 15);
  vector<int> aIdx = {};
  vector<int> bIdx = {};

  for (int i = 0; i < bAtoms.size(); i++) {
    vector <int> abContactPts = psA.getPointsWithin(bAtoms[i], 0.0, conDist, false); // checks if the given atom in B is w/in [conDist] Angstroms from any in A
    if (abContactPts.size() > 0) {
      bIdx.push_back(i);
      for (int ii = 0; ii < abContactPts.size(); ii++) {
        int aPoint = abContactPts[ii];
        if (std::find(aIdx.begin(), aIdx.end(), aPoint) == aIdx.end()) {
          aIdx.push_back(aPoint);
        }
      }
    }
  }
  vector<vector<int>> abIdxVector = {aIdx,bIdx};
  return abIdxVector;
}

// gets all the residue contacts between a & b as a dict, with the B indexes as keys and A indexes as a list of values
map<int,vector<int>> getABcontactDict(AtomPointerVector aAtoms,AtomPointerVector bAtoms,mstreal conDist, ProximitySearch psA) {

  map<int,vector<int>> abDict;

  for (int i = 0; i < bAtoms.size(); i++) {
    vector <int> abContactPts = psA.getPointsWithin(bAtoms[i], 0.0, conDist, false); // checks if the given atom in B is w/in 5 Angstroms from any in A
    if (abContactPts.size() > 0) {
      for (int ii = 0; ii < abContactPts.size(); ii++) {

        int currResB = bAtoms[i]->getResidue()->getResidueIndex();
        int currResA = aAtoms[abContactPts[ii]]->getResidue()->getResidueIndex();

        auto it = abDict.find(currResB);
        if (it != abDict.end()) {
          vector <int> currBcons = abDict[currResB];
          bool found = (std::find(currBcons.begin(), currBcons.end(), currResA) != currBcons.end());
          if (!found) {
            abDict[currResB].push_back(currResA);
          }
        }
        else {
          abDict[currResB] = {currResA};
        }
      }
    }
  }
  return abDict;
}

// get DOCKQ
mstreal getDOCKQ(Structure* s2evalP, Structure* s2evalAP, AtomPointerVector s2evalAbb, Structure* s2evalBP, AtomPointerVector s2evalBbb, AtomPointerVector realA, AtomPointerVector realAbb, AtomPointerVector realB, AtomPointerVector realBbb, map<int,vector<int>> interfaceIdxMap5, int totResidueSize, vector<int> interfaceIdx10A, vector<int> interfaceIdx10B, AtomPointerVector oriInterface, RMSDCalculator rc, ProximitySearch ps) {

  Structure s2eval = *s2evalP;
  Structure s2evalA = *s2evalAP;
  Structure s2evalB = *s2evalBP;
  AtomPointerVector s2evalAatoms = s2evalA.getAtoms();
  AtomPointerVector s2evalBatoms = s2evalB.getAtoms();

  map<int,vector<int>> s2evalABdict = getABcontactDict(s2evalAatoms, s2evalBatoms, 5.0, ps);

  mstreal totRecoveredSize = 0;
  for (auto it = interfaceIdxMap5.cbegin(); it != interfaceIdxMap5.cend(); ++it) { //go through the original contacts map
    auto it2 = s2evalABdict.find(it->first); //see if the contact key is in the structure to evaluate's map
    if (it2 != s2evalABdict.end()) { //if so
      for (auto it3 = it->second.cbegin(); it3 != it->second.cend(); ++it3) { //go through the original contact's list...
        if (std::find(it2->second.begin(), it2->second.end(), *it3) != it2->second.end()) {
          totRecoveredSize++;
        }
      }
    }
  }

  // ps is for contacts in s2evalA / s2evalB

  mstreal fnat = totRecoveredSize / totResidueSize;

  // get interface RMSD
  AtomPointerVector evalInterface;
  for (int i = 0; i < interfaceIdx10A.size(); i++) { // add supposed-to-be-interface atoms to the evalInterface
    evalInterface.push_back(s2evalAatoms[interfaceIdx10A[i]]);
  }
  for (int i = 0; i < interfaceIdx10B.size(); i++) { // add supposed-to-be-interface atomsto the evalInterface
    evalInterface.push_back(s2evalBatoms[interfaceIdx10B[i]]);
  }

  mstreal iRMSD = rc.bestRMSD(evalInterface,oriInterface);
  mstreal iRMSDscaled = 1 / (1 + pow((iRMSD / 1.5),2));

  // get scaled ligand RMSD
  mstreal ligRMSD = 0;
  if (realAbb.size() > realBbb.size()) {

    rc.align(s2evalAbb,realAbb,*s2evalP);
    ligRMSD = rc.rmsd(s2evalBbb,realBbb);
  }
  else {
    rc.align(s2evalBbb,realBbb,*s2evalP);
    ligRMSD = rc.rmsd(s2evalAbb,realAbb);

  }

  mstreal ligRMSDscaled = 1 / (1 + pow((ligRMSD / 8.5),2));

  //get DOCKQ~!
  mstreal dockQscore = (fnat + iRMSDscaled + ligRMSDscaled) / 3.0;
  return dockQscore;
}

// quickly gets a random rotation and orientation matrix - see 1992 Arvo paper for details https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.53.1357&rep=rep1&type=pdf
Transform getRandRotOrient() {  

    // make Z pole rotation matrix...
    mstreal x1 = MstUtils::randUnit(0,1);
    mstreal zRot[3][3] = {
    {cos(2*M_PI*x1), sin(2*M_PI*x1), 0},
    {-sin(2*M_PI*x1), cos(2*M_PI*x1), 0},
    {0, 0, 1},
    } ;

    // make v (unit vector paralel to zp, the line between the top of the Z pole and the point to rotate to)
    mstreal x2 = MstUtils::randUnit(0,1);
    mstreal x3 = MstUtils::randUnit(0,1);
    mstreal v[3] = {
    cos(2*M_PI*x2) * sqrt(x3),
    sin(2*M_PI*x2) * sqrt(x3),
    sqrt(1 - x3)
    } ;

    mstreal houseHolder [3][3] = { 
    { 1 - (2*v[0]*v[0]) , 0 - (2*v[0]*v[1]) , 0 - (2*v[0]*v[2]) },
    { 0 - (2*v[1]*v[0]) , 1 - (2*v[1]*v[1]) , 0 - (2*v[1]*v[2]) },
    { 0 - (2*v[2]*v[0]) , 0 - (2*v[2]*v[1]) , 1 - (2*v[2]*v[2]) }   
    } ;

    // final matrix M = -(H @ R) (Householder on rotation, negative to flip and preserve left-rightness). Maxtrix multiplication goes along rows (counted by i) of the left matrix (H) and columns (counted by j) of the right matrix (R)
    
    vector<vector<mstreal> > finalMvalues;
    for (int i = 0; i < 3; i++) {

        vector<mstreal> rowValues;
        for (int j = 0; j < 3; j++) {

            mstreal nextValue = 0;
            for (int k = 0; k < 3; k++) {
                nextValue += -1 * houseHolder[i][k] * zRot[k][j];
            }
            rowValues.push_back(nextValue);
        }

        finalMvalues.push_back(rowValues);
    }

    Transform finalM(finalMvalues);
    return(finalM);
}

int main(int argc, char** argv) {
  /*time_t start_time = time(NULL);
  cout << "start time: " << ctime(&start_time) << endl;*/
  // Setup and get the input arguments
  MstOptions op;
  op.setTitle("Generates randomly docked pairs of structures, with optional binding residue constraints, and computes their RMSDs to create a distribution of randomly docked structures. The canonical structure and simulated structure you want to evaluate must have the same atom indexing, so before generating your random structures with this script, use atomMatcher on the canonical & simulated structures");
  op.addOption("m", "use model binding partners for random docking instead of crystal binding partners", false);
  op.addOption("r", "instead of using complex RMSD, use the sum of ligand-RMSD & receptor-RMSD", false);
  op.addOption("dq", "instead of using the sum of complex RMSD, do DOCKQ as the underlying measurement", false);
  op.addOption("mbl", "there are multiple binding locations for a; this flag will take paths to pdb files containing backbone-matched alternative binding locations, separated out by commas", false);
  op.addOption("ta", "one half of the correctly docked structure, in pdb form, after having its backbone paired with that of the simulated docked structure via match_backbones.py; if using antibody mode, this should be the antibody ", true);
  op.addOption("ma", "one half of the simulated docked structure, in pdb form, after having its backbone paired with that of the correctly docked structure via match_backbones.py (see -ta); if using antibody mode, this should be the antibody ", true);
  op.addOption("tb", "the other half of the correctly docked structure, in pdb form, after having its backbone paired with that of the simulated docked structure via match_backbones.py", true);
  op.addOption("mb", "the other half of the simulated docked structure, in pdb form, after having its backbone paired with that of the correctly docked structure via match_backbones.py (see -tb); if using antibody mode, this should be the antibody ", true);
  op.addOption("n", "the number of valid docking positions to generate; defaults to 1 million, which should take around 1 hour per every 500 residues in the complex");
  op.addOption("i", "the number of valid interaction residues required to accept the structure; defaults to 0");
  op.addOption("cla", "the number of clashes allowed in an accepted structure; defaults to 0");
  op.addOption("sd", "the standard deviation of the normal distribution used to pull the docking partners apart from each other; defaults to 2.0; in angstroms");
  op.addOption("al", "an optional list of binding residues for the first docking partner (see -ta); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like 'A,100, ;A,100,A;A,100,B'. If you're only giving binding residues for one of the two partners, it must be this one (i.e. you cannot give bl without giving al)");
  op.addOption("bl", "an optional list of binding residues for the second docking partner (see -tb); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like A,100,;A,100,A;A,100,B.");
  op.addOption("abm", "an optional antibody / TCR mode, which will treat all the loops of the antibody as the binding residues; the antibody structure must use IMGT numbering and its structure must be entered using the -ta argument not the -tb argument.");
  op.addOption("cdr3", "limit the antibody / TCR mode to using CDR3 loop contacts only");
  op.addOption("limA", "Must be used with al or abm, and bl. Optionally limit the range of angles of accepted dockings, wherein each binding partner has a line drawn between its geometric center and the geometric center of its binding residues; one line is used as an axis, and the other is used to calculate the angle of tilt relative to that axis. Primarily useful for TCRpMHC random dockings.");
  op.addOption("cache", "Use a cached distribution of underlying distance metrics instead; supply the path to that csv file here.");
  op.addOption("j", "Just print the underlying distance metric between the two structures instead of doing comparisons.");
  op.addOption("o", "the output file to save the distribution of RMSDs, as a comma separated list", true);
  op.addOption("t", "create testing files; provide the output directory to save the testing files. This also sets the number of dockings to just 10 so you aren't inundated with files on accident :)");
  op.addOption("pbs", "save partner B's structures, for visualization purposes");
  op.setOptions(argc, argv);
  MstUtils::seedRandEngine();

  if (!op.isGiven("al") && !op.isGiven("abm") && op.isGiven("bl")) {
    cout << "you cannot give bl without giving al or abm";
    exit(1);
  }

  if (op.isGiven("limA") && !op.isGiven("bl")) {
    cout << "you cannot give limA without giving (al or abm) and bl, as there's no 'binding angle' between two patches of binding residues, without binding residues defined for each docking partner";
    exit(1);
  }

  int numberDockingsRequired;
  if (op.isGiven("t")) {
    numberDockingsRequired = 10;
  }
  else {
    if (op.isGiven("n")) {
      numberDockingsRequired = stoi(op.getString("n"));
    }
    else {
      numberDockingsRequired = 1000000;
    }
  }

  int contactsRequired;
  if (op.isGiven("i")) {
    contactsRequired = stoi(op.getString("i"));
  }
  else {
    contactsRequired = 1;
  }

  int clashesAllowed;
  if (op.isGiven("cla")) {
    clashesAllowed = stoi(op.getString("cla"));
  }
  else {
    clashesAllowed = 0;
  }

  mstreal normalDistBase;
  if (op.isGiven("sd")) {
    normalDistBase = stof(op.getString("sd"));
  }
  else {
    normalDistBase = 2.0;
  }

  // load backbone or full-atom structures; make sure just working w/ heavy atoms

  Structure TAfull(op.getString("ta"));
  AtomPointerVector TAheavyatoms = getHeavyAtoms(TAfull.getAtoms());
  Structure TAheavy(TAheavyatoms);

  Structure TBfull(op.getString("tb"));
  AtomPointerVector TBheavyatoms = getHeavyAtoms(TBfull.getAtoms());
  Structure TBheavy(TBheavyatoms);

  Structure MAfull(op.getString("ma"));
  AtomPointerVector MAheavyatoms = getHeavyAtoms(MAfull.getAtoms());
  Structure MAheavy(MAheavyatoms);

  Structure MBfull(op.getString("mb"));
  AtomPointerVector MBheavyatoms = getHeavyAtoms(MBfull.getAtoms());
  Structure MBheavy(MBheavyatoms);

  // make a full structure for T, which is connected to TAheavy / TBheavy rather than a copy of them
  Structure T;
  for (int i = 0; i < TAheavy.chainSize(); i++) {
    T.appendChain(&TAheavy.getChain(i));
  }
  for (int i = 0; i < TBheavy.chainSize(); i++) {
    T.appendChain(&TBheavy.getChain(i));
  }

  //make copies of the true struccture, just A from the true structure, and just B from the true structure, to do RMSD calculations with
  Structure Tcopy(T); 
  AtomPointerVector realAtoms = getHeavyAtoms(Tcopy.getAtoms());
  Structure TACopy(TAheavy);
  AtomPointerVector realAtomsTA = TACopy.getAtoms();
  Structure TBCopy(TBheavy);
  AtomPointerVector realAtomsTB = TBCopy.getAtoms();

  ProximitySearch psA(realAtomsTA, 15);
  map<int,vector<int>> abResCons = getABcontactDict(realAtomsTA,realAtomsTB,5.0,psA);

  int conSize = 0;
  for (auto it = abResCons.cbegin(); it != abResCons.cend(); ++it) {
    conSize += it->second.size();
  }

  // get the indexes and then interface structure for A & B, by a 10 Angstrom cutoff, in case DOCKQ is being used as an underlying distance metric

  vector<vector<int>> ab10realIndexes = getABindexes(realAtomsTA,realAtomsTB,10.0);
  vector<int> a10interfaceIndexesAll = ab10realIndexes[0];
  vector<int> b10interfaceIndexesAll = ab10realIndexes[1];
  vector<int> a10interfaceIndexesBB;
  vector<int> b10interfaceIndexesBB;
  AtomPointerVector ab10interface;
  for (int i = 0; i < a10interfaceIndexesAll.size(); i++) {
    Atom* checkingAtom = realAtomsTA[a10interfaceIndexesAll[i]];
    if ((checkingAtom->getName() == "CA") or (checkingAtom->getName() == "C") or (checkingAtom->getName() == "N") or (checkingAtom->getName() == "O")) {
      ab10interface.push_back(checkingAtom);
      a10interfaceIndexesBB.push_back(a10interfaceIndexesAll[i]);
    }
  }

  for (int i = 0; i < b10interfaceIndexesAll.size(); i++) {
    Atom* checkingAtom = realAtomsTB[b10interfaceIndexesAll[i]];
    if ((checkingAtom->getName() == "CA") or (checkingAtom->getName() == "C") or (checkingAtom->getName() == "N") or (checkingAtom->getName() == "O")) {
      ab10interface.push_back(checkingAtom);
      b10interfaceIndexesBB.push_back(b10interfaceIndexesAll[i]);
    }
  }

  // make a full structure for M (the model), which is connected to MAheavy / MBheavy rather than a copy of them
  Structure M;
  for (int i = 0; i < MAheavy.chainSize(); i++) {
    M.appendChain(&MAheavy.getChain(i));
  }
  for (int i = 0; i < MBheavy.chainSize(); i++) {
    M.appendChain(&MBheavy.getChain(i));
  }

  // if you're using model-structure partners for docking, make A & B out of them (MAheavy & MB heavy); otherwise stick with the crystal partners (TAheavy & TBheavy) which is default

  Structure A;
  Structure B;

  vector<AtomPointerVector> altBindingLocs;
  vector<AtomPointerVector> altBindingLocsBB;
  vector<AtomPointerVector> altBindingLocsAonly;
  vector<AtomPointerVector> altBindingLocsAbbOnly;
  Structure TAalt;
  Structure TBCOPY;

  if (op.isGiven("mbl")) {
    vector<string> altBindingAstructs = MstUtils::split(op.getString("mbl"), ",");
    for (int i = 0; i < altBindingAstructs.size(); i++) {
      TBCOPY = Structure(TBheavy);
      string taAltPath = altBindingAstructs[i];
      TAalt = Structure(taAltPath);
      AtomPointerVector taAltAtoms = getHeavyAtoms(TAalt.getAtoms());

      altBindingLocsAonly.push_back(taAltAtoms);
      altBindingLocsAbbOnly.push_back(getBBatoms(taAltAtoms));

      for (int ii = 0; ii < TBCOPY.chainSize(); ii++) {
        TAalt.appendChain(&TBCOPY.getChain(ii));
      }
      AtomPointerVector cAltAtoms = getHeavyAtoms(TAalt.getAtoms());
      
      altBindingLocs.push_back(cAltAtoms);
      altBindingLocsBB.push_back(getBBatoms(cAltAtoms));
    }
  }

  // if alt binding locations: make the indexes for the interface, by a 5 angstrom cutoff, in case DOCKQ is being used as an underlying distance metric 

   vector<map<int,vector<int>>> abResConsAlts;
   vector<int> conSizeAlts;
   for (int i = 0; i < altBindingLocsAonly.size(); i++) {
    ProximitySearch psAlt(altBindingLocsAonly[i], 15);

    map<int,vector<int>> abResConsAlt = getABcontactDict(altBindingLocsAonly[i],realAtomsTB,5.0,psAlt);

    int conSizeAlt = 0;
    for (auto it = abResConsAlt.cbegin(); it != abResConsAlt.cend(); ++it) {
      conSizeAlt += it->second.size();
    }

    abResConsAlts.push_back(abResConsAlt);
    conSizeAlts.push_back(conSizeAlt);
   }

   //plus if alt binding locations, get the indexes / interface Structure for a 10 angstrom cutoff - also in case DOCKQ is used 
  vector<vector<int>> a10interfaceIndexesBBAlts;
  vector<vector<int>> b10interfaceIndexesBBAlts;
  vector<AtomPointerVector> ab10interfaceAlts;
  for (int i = 0; i < altBindingLocsAonly.size(); i++) {
    vector<vector<int>> ab10realIndexesAlt = getABindexes(altBindingLocsAonly[i],realAtomsTB,10.0);
    vector<int> a10interfaceIndexesAllAlt = ab10realIndexes[0];
    vector<int> b10interfaceIndexesAllAlt = ab10realIndexes[1];
    vector<int> a10interfaceIndexesBBAlt;
    vector<int> b10interfaceIndexesBBAlt;
    AtomPointerVector ab10interfaceAlt;  

    for (int ii = 0; ii < a10interfaceIndexesAllAlt.size(); ii++) {

      int altInterfaceIndx = a10interfaceIndexesAllAlt[ii];

      Atom* checkingAtom = altBindingLocsAonly[i][altInterfaceIndx];

      if ((checkingAtom->getName() == "CA") or (checkingAtom->getName() == "C") or (checkingAtom->getName() == "N") or (checkingAtom->getName() == "O")) {
        ab10interfaceAlt.push_back(checkingAtom);
        a10interfaceIndexesBBAlt.push_back(altInterfaceIndx);
      }
    } 

    for (int ii = 0; ii < b10interfaceIndexesAllAlt.size(); ii++) {
      Atom* checkingAtom = realAtomsTB[b10interfaceIndexesAllAlt[ii]];
      if ((checkingAtom->getName() == "CA") or (checkingAtom->getName() == "C") or (checkingAtom->getName() == "N") or (checkingAtom->getName() == "O")) {
        ab10interfaceAlt.push_back(checkingAtom);
        b10interfaceIndexesBBAlt.push_back(b10interfaceIndexesAllAlt[ii]);
      }
    }

    a10interfaceIndexesBBAlts.push_back(a10interfaceIndexesBBAlt);
    b10interfaceIndexesBBAlts.push_back(b10interfaceIndexesBBAlt);
    ab10interfaceAlts.push_back(ab10interfaceAlt);
  }

  if (op.isGiven("m")) {
    A = Structure(MAheavy); // this duplicates / does not connect A & MAheavy
    B = Structure(MBheavy);
  }
  else {
    A = Structure(TAheavy);
    B = Structure(TBheavy);
  }

  // make a full structure for W (our "working" structure, made up of partners A and B that will dock a million times), which is connected to A / B rather than a copy of them
  Structure W;
  for (int i = 0; i < A.chainSize(); i++) {
    W.appendChain(&A.getChain(i));
  }
  for (int i = 0; i < B.chainSize(); i++) {
    W.appendChain(&B.getChain(i));
  }

  vector <Residue*> Areses = A.getResidues();
  AtomPointerVector Aatoms = A.getAtoms();

  vector <Residue*> Breses = B.getResidues();
  AtomPointerVector Batoms = B.getAtoms();

  vector <Residue*> MAreses = MAheavy.getResidues();
  AtomPointerVector MAatoms = MAheavy.getAtoms();

  // for DOCKQ to get fnat later on

  ProximitySearch psDQ2(MAatoms, 15);

  vector <Residue*> MBreses = MBheavy.getResidues();
  //AtomPointerVector MBatoms = getHeavyAtoms(MBheavy.getAtoms());
  AtomPointerVector MBatoms = MBheavy.getAtoms();

  AtomPointerVector AbbAtoms = getBBatoms(Aatoms);
  AtomPointerVector BbbAtoms = getBBatoms(Batoms);
  AtomPointerVector MAbbAtoms = getBBatoms(MAatoms);
  AtomPointerVector MBbbAtoms = getBBatoms(MBatoms);
  AtomPointerVector realAtomsCAbb = getBBatoms(realAtomsTA);
  AtomPointerVector realAtomsCBbb = getBBatoms(realAtomsTB);

  RMSDCalculator rc;
  vector <mstreal> rmsdList {};

  if ((!op.isGiven("cache")) && (!op.isGiven("j"))) {

    // position 1
    if (op.isGiven("t")) {
      W.writePDB(op.getString("t") + "position1" + ".pdb");
    }

    // translate each to the origin

    Transform TCA0 = TransformFactory::translate(-Aatoms.getGeometricCenter());
    TCA0.apply(A);

    Transform TCB0 = TransformFactory::translate(-Batoms.getGeometricCenter());
    TCB0.apply(B);

    // position 2
    if (op.isGiven("t")) {
      W.writePDB(op.getString("t") + "position2" + ".pdb");
    }

    // if binding residues are given for A (or A has binding residues by virtue of this being run in antibody-antigen binding mode) get them; same for B

    AtomPointerVector AbinderAtoms;
    AtomPointerVector BbinderAtoms;

    vector<tuple <string, int, char>> aResesDetails; // stores details for binding residues for A
    vector<tuple <string, int, char>> bResesDetails; // stores details for binding residues for B

    if (op.isGiven("abm")) {
      for (int i = 0; i < Areses.size(); i++) {
        Residue* currR = Areses[i];
        int resNum = currR->getNum();
        bool isLoop = false;

        //check if loop by normal numbering, unless mm mode, in which case check if loop by crystal numbering

        if (op.isGiven("cdr3") && (105 <= resNum && resNum <= 117)) {
          isLoop = true;
        }

        else if ((27 <= resNum && resNum <= 38) || (56 <= resNum && resNum <= 65) || (105 <= resNum && resNum <= 117))
        {
          isLoop = true;
        }

        if (isLoop == true) {
          AtomPointerVector taResAtoms = currR->getAtoms();
          for (int ca = 0; ca < taResAtoms.size(); ca++) {
            Atom* currA = taResAtoms[ca];
            string cResAtomName = currA->getName();
            if (cResAtomName == "CA") { // only get indexes for CAs, as contacts will be defined as inter-A-distance of 8 angstroms or less
              AbinderAtoms.push_back(currA);
            } 
          }
        }
      }
      //rotate the binding residues in A up
      CartesianPoint geoCenterA = AbinderAtoms.getGeometricCenter();
      Transform TZ = TransformFactory::alignVectorWithZAxis(geoCenterA);
      TZ.apply(A);

      if (op.isGiven("t")) {
        W.writePDB(op.getString("t") + "position3" + ".pdb");
      }
    }

    else if (op.isGiven("al")) {

      // get the details for A's binding residues
      string aBinders = op.getString("al");
      vector<string> aResSplit1 = MstUtils::split(aBinders, ";");
      for (int i = 0; i < aResSplit1.size(); i++) {
          vector<string> aResSplit2 = MstUtils::split(aResSplit1[i], ",");
          string resChain = aResSplit2[0];
          int resNum = std::stoi(aResSplit2[1]);
          string resIcodeBase = aResSplit2[2];
          char resIcode = resIcodeBase[0];
          auto aResDetails = std::make_tuple (resChain, resNum, resIcode);
          aResesDetails.push_back(aResDetails);
      }

      // go over each residue, and if its chain / number / ID code match, add it to the binding residues list! Also make a list of the binder indexes

      for (int i = 0; i < Areses.size(); i++) {
        Residue* currR = Areses[i];
        string resChain = currR->getChainID();
        int resNum = currR->getNum();
        char resIcode = currR->getIcode();
        for (int ii = 0; ii < aResesDetails.size(); ii++) {
          auto aResTuple = aResesDetails[ii];
          bool isbindingRes = false;

          //check if binding residue by normal numbering, unless mm mode, in which case check if binding residue by crystal numbering

          if ((get<0>(aResTuple) == resChain) && (get<1>(aResTuple) == resNum) && (get<2>(aResTuple) == resIcode)) {
            isbindingRes = true;
          }
          //}

          if (isbindingRes) {
            AtomPointerVector taResAtoms = currR->getAtoms();
            for (int ca = 0; ca < taResAtoms.size(); ca++) {
              Atom* currA = taResAtoms[ca];
              string cResAtomName = currA->getName();
              if (cResAtomName == "CA") {
                AbinderAtoms.push_back(currA);
              }
            }
          }
        }
      }
      //rotate the binding residues in A up

      CartesianPoint geoCenterA = AbinderAtoms.getGeometricCenter();
      Transform TZ = TransformFactory::alignVectorWithZAxis(geoCenterA);
      TZ.apply(A);

      if (op.isGiven("t")) {
        W.writePDB(op.getString("t") + "position3" + ".pdb");
      }
    }

    //***
    else {
      for (int a = 0; a < Aatoms.size(); a++) {
        Atom* currAtom = Aatoms[a];
        string currName = currAtom->getName();
        if (currName == "CA") { 
          AbinderAtoms.push_back(currAtom);
        }
      }
    }

    // if binding residues are given for B, pre-store the details to iterate over easily

    // *** CHECK OVER THIS WHOOOOLE SECTION AND THE ABOVE - MAKE SURE IT'S ALL RIGHT, THEN MAKE SURE THE "i" CHECKING IS BETWEEN DISCRETE RESIDUES ONLY

    if (op.isGiven("bl")) {
      string bBinders = op.getString("bl");
      vector<string> bResSplit1 = MstUtils::split(bBinders, ";");
      for (int i = 0; i < bResSplit1.size(); i++) {
          vector<string> bResSplit2 = MstUtils::split(bResSplit1[i], ",");
          string resChain = bResSplit2[0];
          int resNum = std::stoi(bResSplit2[1]);
          char resIcode = bResSplit2[2][0];
          auto bResDetails = std::make_tuple (resChain, resNum, resIcode);
          bResesDetails.push_back(bResDetails);
      }

      for (int i = 0; i < Breses.size(); i++) {
        Residue* currR = Breses[i];
        string resChain = currR->getChainID();
        int resNum = currR->getNum();
        char resIcode = currR->getIcode();
        for (int ii = 0; ii < bResesDetails.size(); ii++) {
          auto bResTuple = bResesDetails[ii];

          bool isbindingRes = false;

          //check if binding residue by normal numbering, unless mm mode, in which case check if binding residue by crystal numbering

          if ((get<0>(bResTuple) == resChain) && (get<1>(bResTuple) == resNum) && (get<2>(bResTuple) == resIcode)) {
            isbindingRes = true;
          }
          //}

          if (isbindingRes) {
            AtomPointerVector taResAtoms = currR->getAtoms();
            for (int ca = 0; ca < taResAtoms.size(); ca++) {
              Atom* currAtom = taResAtoms[ca];
              string currName = currAtom->getName();
              if (currName == "CA") {
                BbinderAtoms.push_back(currAtom);
              }
            }
          }
        }
      }
    }

    else {
      vector<int> bBinderIndexes; // store the atom indexes of B binding atoms
      for (int i = 0; i < Batoms.size(); i++) { // go over all atoms in B...
        Atom* currAtom = Batoms[i];
        string currAtomName = currAtom->getName();
        if (!(op.isGiven("dq")) && (currAtomName != "CA")) { // if not DOCKQ mode, and not a A atom...
          continue;
        }
        Residue* currRes = currAtom->getParent();
        string currChain = currRes->getChainID();
        int currResNum = currRes->getNum();
        char currResIcode = currRes->getIcode();
        for (int ii = 0; ii < bResesDetails.size(); ii++) {
          auto bResTuple = bResesDetails[ii];
          if ((get<0>(bResTuple) == currChain) && (get<1>(bResTuple) == currResNum) && (get<2>(bResTuple) == currResIcode)) { // if the residue is a binding residue...
            bBinderIndexes.push_back(i);
          }
        }
      }
    }

    // let's get started on the random docking now! Set up variables (d = number of random dockings required, t = a variable used to count the initial 10 structures to be saved for testing mode)

    int d = 0;
    int t = 0;
    int fails = 0;

    // set up important vectors...

    // keep a list of indexes of previous clashes
    vector <Atom*> recentClashes;
    vector <int> recentClashIndexes;
    // and have a list of current clashes
    vector <Atom*> currClashes;
    vector <int> currClashIndexes;

    // set the curr coords of A to be alt cords, so that they can be reset every new random docking

    for (int a = 0; a < Aatoms.size(); a++) {
      Atom* currAtom = Aatoms[a];
      currAtom->clearAlternatives();
      currAtom->addAlternative(currAtom);
    }

    // set the current coords of B to be alternate coords, so that they can be reset every new random docking - also fill out which indexes are possible to make contacts in B, and get their atoms

    for (int a = 0; a < Batoms.size(); a++) {
      Atom* currAtom = Batoms[a];
      currAtom->clearAlternatives();
      currAtom->addAlternative(currAtom);
      // as well as get indexes for the potential contact residues, if not already gotten from abm / al / bl flags being given
      if (op.isGiven("bl")) {
        continue;
      }
      else {
        string currName = currAtom->getName();
        if (currName == "CA") {
          BbinderAtoms.push_back(currAtom);
        }
      }
    }

    // set up xyz high & low values to be used later in the loop, for checking xLow & xHigh when pulling out proteins

    mstreal axLow;
    mstreal axHigh;
    mstreal ayLow;
    mstreal ayHigh;
    mstreal azLow;
    mstreal azHigh;

    mstreal bxLow;
    mstreal bxHigh;
    mstreal byLow;
    mstreal byHigh;
    mstreal bzLow;
    mstreal bzHigh;

    //& docking requirements

    mstreal clashDistance = 3.0;
    mstreal contactDistance = 8.0;



    // set up proximity search with A, since B will be the partner moving

    ProximitySearch psBB(AbbAtoms, 15);

    // for DOCKQ
    ProximitySearch psDQ(Aatoms, 15);

    // and a proximity search with just A's A binding atoms - this is a ProximitySearch*

    ProximitySearch psAbinders(AbinderAtoms, 15);


    mstreal angleLim = 180;

    if (op.isGiven("limA")) {
      mstreal angleLim1 = stoi(op.getString("limA"));
      if (angleLim1 < 180) {
        angleLim = angleLim1;
      }
      else {
        cout << "limA argument must be a number less than 180 (180 degrees would mean no limit)";
        exit(0);
      }
    }

    // find the widest angle between two points in the binding atoms of A, and the geometric center of A; limit rotations to be at that angle or less, to skew away from useless dockings with no contacts between binding residues in A and B

    if (op.isGiven("al") || op.isGiven("abm")) {

      mstreal largestAngle = 0;
      CartesianPoint geoCenterA = Aatoms.getGeometricCenter();

      for (int ai = 0; ai < AbinderAtoms.size(); ai++) {

        Atom* firstAtom = AbinderAtoms[ai];

        for (int aii = ai + 1; aii < AbinderAtoms.size(); aii++) {

          Atom* secondAtom = AbinderAtoms[aii];
          mstreal angleOff = CartesianGeometry::angle(firstAtom,geoCenterA,secondAtom);

          if (angleOff > largestAngle) {
            largestAngle = angleOff;
          }
        }
      }

      if (largestAngle < angleLim) {
        angleLim = largestAngle;
      }

    }

    mstreal polCos;
    mstreal polSin;
    mstreal azCos;
    mstreal azSin;
    
    mstreal intAnegLim = cos(angleLim*M_PI/180);
    mstreal intAposLim = 1;

    mstreal randPolar;
    mstreal randAz;
    mstreal randX;
    mstreal randY;
    mstreal randZ;
    int backNforth;

    while (d < numberDockingsRequired) { // while the number of docked structures is still insufficient...

      if (op.isGiven("t")) {
        t++;
      }

      if (op.isGiven("t")) {
        W.writePDB(op.getString("t") + "position4" + to_string(d) + ".pdb");
      }

      // reset to alt coors

      if (op.isGiven("r")) {
        for (int i = 0; i < Areses.size(); i++) {
          Areses[i]->makeAlternativeMain(0); // reset to origin position
        }
      }

      for (int i = 0; i < Breses.size(); i++) {
        Breses[i]->makeAlternativeMain(0); // reset to origin position
      }

      //position 4

      if (op.isGiven("t")) {
        W.writePDB(op.getString("t") + "position5" + to_string(d) + ".pdb");
      }

      // rotate B randomly using Fast Random Rotation Matrixes

      Transform randRot = getRandRotOrient();
      randRot.apply(B);

        //position 5
      if ((op.isGiven("t")) && (t == 0)) {
          W.writePDB(op.getString("t") + "position6" + ".pdb");
      }

      // generate a random unit vector to pull B away by - if limiting binding angles between binding partners, implement that when getting polCos (the cosine of the polar angle is pre-computed above!)

      polCos = MstUtils::randUnit(intAnegLim,intAposLim);
      polSin = sqrt(1 - (polCos*polCos));

      randAz = MstUtils::randUnit(0,2*M_PI);
      randX = polSin * cos(randAz);
      randY = polSin * sin(randAz);
      randZ = polCos;
      
      /*
      randX = MstUtils::randNormal(0,1.0);
      while ((randX < -1.0) || (randX > 1.0)) {
        randX = MstUtils::randNormal(0,1.0);
      }
      mstreal randY = MstUtils::randNormal(0,1.0);
      while ((randY < -1.0) || (randY > 1.0)) {
        randY = MstUtils::randNormal(0,1.0);
      }
      mstreal randZ = MstUtils::randNormal(0,1.0);
      while ((randZ < -1.0) || (randZ > 1.0)) {
        randZ = MstUtils::randNormal(0,1.0);
      }
      mstreal denom = sqrt((randX*randX) + (randY*randY) + (randZ*randZ));
      mstreal randXstepScale = randX/denom;
      mstreal randYstepScale = randY/denom;
      mstreal randZstepScale = randZ/denom;*/

      // sort the atom pointers of B according to the largest value out of randX / randY /randZ, so you can check them intelligently along the vector it's being most pulled out along

      /*
      char maxAbsAxis;
      int negatorInt;
      if (abs(randX) > abs(randY)) {
        if (abs(randX) > abs(randZ)) {
          maxAbsAxis = 'X';
          if (randX <= 0) {
            negatorInt = -1;
          }
          else {
            negatorInt = 1;
          }
        }
        else {
          maxAbsAxis = 'Z';
          if (randZ < 0) {
            negatorInt = -1;
          }
          else {
            negatorInt = 1;
          }
        }
      } 
      else {
        if (abs(randY) > abs(randZ)) {
          maxAbsAxis = 'Y';
          if (randY < 0) {
            negatorInt = -1;
          }
          else {
            negatorInt = 1;
          }
        }
        else {
          maxAbsAxis = 'Z';
          if (randZ < 0) {
            negatorInt = -1;
          }
          else {
            negatorInt = 1;
          }
        }
      }

      
      vector<tuple<mstreal,Atom*>> axSortBatomsTuple;
      AtomPointerVector axSortBatoms;
      for (int axi = 0; axi < Batoms.size(); axi++) {

        Atom* currAtomCheck = Batoms[axi];
        mstreal currAxval;

        if (maxAbsAxis == 'X') {
          currAxval = currAtomCheck->getX() * negatorInt;
        }
        if (maxAbsAxis == 'Y') {
          currAxval = currAtomCheck->getY() * negatorInt;
        }
        if (maxAbsAxis == 'Z') {
          currAxval = currAtomCheck->getZ() * negatorInt;
        }
        
        axSortBatomsTuple.push_back(std::make_tuple(currAxval,currAtomCheck));
      }
      sort(axSortBatomsTuple.begin(), axSortBatomsTuple.end());
      for (int axi = 0; axi < Batoms.size(); axi++) {
        axSortBatoms.push_back(get<1>(axSortBatomsTuple[axi]));
      }*/

      // while a proper docked stage has not been reached, but there are still opporunitites for success

      bool absoluteSuccess = false;
      bool absoluteFailure = false;

      /// reset recent clashes, if left-over from the last round
      recentClashes.resize(0);

      mstreal contactsCount;
      mstreal zRand;

      while ((absoluteSuccess == false) && (absoluteFailure == false)) {

        // reset current clashes to be empty
        //currClashes.resize(0);

        // pull B atoms out along the X axis, in small steps with a normal distrubtion of ~1.0 angstrom and mew of 0, except always positive so it's pulling away from A - can try variations on the 1.0 and adjust

        mstreal randStep = abs(MstUtils::randNormal(0,normalDistBase));

        for (int aa = 0; aa < Batoms.size(); aa++) {
            Atom* currAtom = Batoms[aa];
            mstreal xStart = currAtom->getX();
            currAtom->setX(xStart + (randX*randStep));
            mstreal yStart = currAtom->getY();
            currAtom->setY(yStart + (randY*randStep));
            mstreal zStart = currAtom->getZ();
            currAtom->setZ(zStart + (randZ*randStep));
        }

        //position 8
        if (op.isGiven("t") && (t == 0)) {
          W.writePDB(op.getString("t") + "position7" + ".pdb");
        }

        // reject / continue depending on if there are too many clashes or not

        int clashCount = 0;
        currClashes.resize(0);
        currClashIndexes.resize(0);

        // first check positions in B that previously had clashes... recentClashes holds Atom*s and their index numbers in Batoms, for atoms in B that clashed previously

        for (int aa = 0; aa < recentClashes.size(); aa++) {
          Atom* currBatom = recentClashes[aa];
          CartesianPoint currBcoor = currBatom->getCoor();
          vector <int> currClashPts = psBB.getPointsWithin(currBcoor, 0, clashDistance, false);
          int currBsize = currClashPts.size();

          if (currBsize > 0) {
            currClashes.push_back(currBatom);
            currClashIndexes.push_back(recentClashIndexes[aa]);
            clashCount+=currBsize;
          }

          if (clashCount > clashesAllowed) {
            break;
          }
        }

        if (clashCount > clashesAllowed) {
          continue;
        }

        // then if you're not over the clash limit, check the rest of the valid residues, and add to recentClashes
        
        int a;
        backNforth = 0;
        if (recentClashIndexes.size() > 0) {
          a = recentClashIndexes.back() + 1;
        }
        else {
          a = floor(BbbAtoms.size() / 2);
        }

        while ((a > 0) and (a < BbbAtoms.size() - 1)) { //a goes back n forth around the most recent clashes to check nearby atoms, or if no recent clashes, starts in the middle then goes back n forth

          a += backNforth;

          if (backNforth > 0) {
            backNforth = (backNforth + 1)*-1;
          }
          else if (backNforth < 0) {
            backNforth = (backNforth - 1)*-1;
          }
          else {
            backNforth += 1;
          }

          if (std::find(recentClashIndexes.begin(), recentClashIndexes.end(), a) != recentClashIndexes.end()) {
            continue; // don't double-check, if already checked above!
          }

          Atom* currBatom = BbbAtoms[a];
          CartesianPoint currBcoor = currBatom->getCoor();
          vector <int> currClashPts = psBB.getPointsWithin(currBcoor, 0, clashDistance, false);
          int currBsize = currClashPts.size();
          if (currBsize > 0) {
            currClashes.push_back(currBatom);
            currClashIndexes.push_back(a);
            clashCount+=currBsize;
          }

          if (clashCount > clashesAllowed) {
            break;
          }
        }

        if ((clashCount <= clashesAllowed) && (a > 0)) { // if still few clashes, and the "end" of the atompointervector was hit, but the start of it still needs to be reached...
          
          if (a >= BbbAtoms.size()) {
            a = BbbAtoms.size() - 1;
          }

          //a is currently = Batoms.size() - 1, aka the last index; backNforth is negative

          a += backNforth;

          while (a >= 0) {

            if (std::find(recentClashIndexes.begin(), recentClashIndexes.end(), a) != recentClashIndexes.end()) {
              a--;
              continue; // don't double-check, if already checked above!
            }

            Atom* currBatom = BbbAtoms[a];
            CartesianPoint currBcoor = currBatom->getCoor();
            vector <int> currClashPts = psBB.getPointsWithin(currBcoor, 0, clashDistance, false);
            int currBsize = currClashPts.size();
            if (currBsize > 0) {
              currClashes.push_back(currBatom);
              currClashIndexes.push_back(a);
              clashCount+=currBsize;
            }
            if (clashCount > clashesAllowed) {
              break;
            }
            a--;
          }
        }

        if ((clashCount <= clashesAllowed) && (a < BbbAtoms.size())) { // if still few clashes, and the "start" of the atompointervector was hit, but the end of it still needs to be reached...

          if (a < 0) {
            a = 0;
          }

          //a is currently = 0, aka the first index; backNforth is positive
          a += backNforth;

          while (a < BbbAtoms.size()) {

            if (std::find(recentClashIndexes.begin(), recentClashIndexes.end(), a) != recentClashIndexes.end()) {
              a++;
              continue; // don't double-check, if already checked above!
            }

            Atom* currBatom = BbbAtoms[a];
            CartesianPoint currBcoor = currBatom->getCoor();
            vector <int> currClashPts = psBB.getPointsWithin(currBcoor, 0, clashDistance, false);
            int currBsize = currClashPts.size();
            if (currBsize > 0) {
              currClashes.push_back(currBatom);
              currClashIndexes.push_back(a);
              clashCount+=currBsize;
            }
            if (clashCount > clashesAllowed) {
              break;
            }
            a++;
          }
        }

        //change recentClashes to hold the value of currClashes, and the same for their respective indexes

        recentClashes.resize(0);
        recentClashIndexes.resize(0);
        for (int r = 0; r < currClashes.size(); r++) {
          recentClashes.push_back(currClashes[r]);
          recentClashIndexes.push_back(currClashIndexes[r]);
        }

        if (clashCount > clashesAllowed) {
          continue;
        }

        //position 9
        if (op.isGiven("t") && (t == 0)) {
          W.writePDB(op.getString("t") + "position8" + ".pdb");
        }

        // we've only reached this point if there are an aceeptable number of clashes, so accept if there are at least the required number of contacts - possibly limited by antibody binding mode, a binding residues, and/or b binding residues

        if (op.isGiven("i")) {

          contactsCount = 0;
          for (int b = 0; b < BbinderAtoms.size(); b++) { // go over all CA atoms in B possible to be part of contacts

            Atom* currBatom = BbinderAtoms[b];
            CartesianPoint currBcoor = currBatom->getCoor();
            vector <int> currContactPts = psAbinders.getPointsWithin(currBcoor, 0.0, contactDistance, false); // find all CA atoms in A it's in conract with
            contactsCount += currContactPts.size();

            if (contactsCount >= contactsRequired) {
              absoluteSuccess = true;
              break;
            }
          }
          if (contactsCount < contactsRequired) {
            absoluteFailure = true;
          }
        }
        else {
          absoluteSuccess = true;
        }
      }

      if (absoluteFailure) {
        fails++;
        if (op.isGiven("t")) {
          W.writePDB(op.getString("t") + "position9_failedWithContacts" + to_string(contactsCount) + ".pdb");
        }
        continue;
      }
      else {
        d++; // accepted! :M

        if (op.isGiven("t")) {
          W.writePDB(op.getString("t") + "dock" + to_string(d) + "withCon" + to_string(contactsCount) + ".pdb");
        }

        if (op.isGiven("pbs")) {
          B.writePDB(op.getString("pbs") + to_string(d) + "dock.pdb");
        }
      }

      /*if (d%10000 == 0) {
        cout << d << " total docks accepted..." << endl;
      }*/

      // for those accepted, compare to the correct structure to calculate the RMSD

      mstreal simulatedRealRMSDorDOCKQ;

      if (op.isGiven("r")) {
        rc.align(Aatoms,realAtomsTA,W);
        if (op.isGiven("t") && (t == 0)) {
          W.writePDB(op.getString("t") + "Tposition10.pdb");
          TACopy.writePDB(op.getString("t") + "TAposition10.pdb");
        }
        mstreal bRMSD = rc.rmsd(Batoms,realAtomsTB); // gets RMSD without alignment
        rc.align(Batoms,realAtomsTB,W); // aligns B with the original location of B, but transforms the whole structure W
        if (op.isGiven("t") && (t == 0)) {
          W.writePDB(op.getString("t") + "Tposition11.pdb");
          TBCopy.writePDB(op.getString("t") + "CBTposition11.pdb");
        }
        mstreal aRMSD = rc.rmsd(Aatoms,realAtomsTA); // gets RMSD without alignment
        simulatedRealRMSDorDOCKQ = (aRMSD + bRMSD);

        if (op.isGiven("t")) {
          if (simulatedRealRMSDorDOCKQ < 5) {
            W.writePDB(op.getString("t") + "lowDock" + to_string(simulatedRealRMSDorDOCKQ) + ".pdb");
          }
        }

        if (op.isGiven("mbl")) {
          for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations

            AtomPointerVector altBindingAtoms = altBindingLocsAonly[abl];
            mstreal aaRMSD = rc.rmsd(Aatoms,altBindingAtoms); // gets RMSD without alignment - because W's already aligned by origianl location of B

            rc.align(Aatoms,altBindingAtoms,W); // aligns A with the alt original location of A, but transforms the whole structure W
            mstreal abRMSD = rc.rmsd(Batoms,realAtomsTB); // gets RMSD without alignment

            mstreal bestArmsd = min(aaRMSD, aRMSD);
            mstreal bestBrmsd = min(abRMSD, bRMSD);
            mstreal altSimRealRMSD = (bestArmsd + bestBrmsd);
            if (op.isGiven("t")) {
              cout << "simulatedRealRMSDorDOCKQ was: " << simulatedRealRMSDorDOCKQ << endl;
              cout << "best simulatedRealRMSDorDOCKQ taking into acct alt A locations was: " << altSimRealRMSD << endl;
            }
            if (altSimRealRMSD < simulatedRealRMSDorDOCKQ) {
              simulatedRealRMSDorDOCKQ = altSimRealRMSD;
            }
            if (op.isGiven("t")) {
              cout << "lower RMSD picked for comparison was: " << simulatedRealRMSDorDOCKQ << endl;
            }
          }
        }
      }

      else if (op.isGiven("dq")) {

        if (op.isGiven("t")) {
          cout << "testing AbbAtoms[0]" << AbbAtoms[0] << endl;
          cout << "testing getBBatoms(Aatoms)[0]" << getBBatoms(Aatoms)[0] << endl;
        }

        mstreal currDOCKQ = getDOCKQ(&W, &A, AbbAtoms, &B, BbbAtoms, realAtomsTA, realAtomsCAbb, realAtomsTB, realAtomsCBbb, abResCons, conSize, a10interfaceIndexesBB, b10interfaceIndexesBB, ab10interface, rc, psDQ);
        simulatedRealRMSDorDOCKQ = currDOCKQ;

        if (op.isGiven("t")) {
          cout << "complex DOCKQ was: " << currDOCKQ << endl;
        }
        if (op.isGiven("mbl")) {
          for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations

            mstreal altDOCKQ = getDOCKQ(&W, &A, AbbAtoms, &B, BbbAtoms, altBindingLocsAonly[abl], altBindingLocsAbbOnly[abl], realAtomsTB, realAtomsCBbb, abResConsAlts[abl], conSizeAlts[abl], a10interfaceIndexesBBAlts[abl], b10interfaceIndexesBBAlts[abl], ab10interfaceAlts[abl], rc, psDQ);
            if (op.isGiven("t")) {
              cout << "altDOCKQ was: " << altDOCKQ << endl;
              cout << "altBindingLocsAonly[abl][0] was: " << altBindingLocsAonly[abl][0] << endl;
              cout << "realAtomsTA[0] was: " << realAtomsTA[0] << endl;
              cout << "*realAtomsTA[0] was: " << *realAtomsTA[0] << endl;
            }
            if (altDOCKQ > simulatedRealRMSDorDOCKQ) {
              simulatedRealRMSDorDOCKQ = altDOCKQ;
            }
            if (op.isGiven("t")) {
              cout << "best DOCKQ picked for comparison was: " << simulatedRealRMSDorDOCKQ << endl;
            }
          }
        }
      }

      else {
        AtomPointerVector Catoms = W.getAtoms();
        simulatedRealRMSDorDOCKQ = rc.bestRMSD(Catoms,realAtoms);

        if (op.isGiven("t")) {
          cout << "complex RMSD was: " << simulatedRealRMSDorDOCKQ << endl;
        }

        if (op.isGiven("mbl")) {
          for (int abl = 0; abl < altBindingLocs.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations
            AtomPointerVector altBindingAtoms = altBindingLocs[abl];
            mstreal altSimRealRMSD = rc.bestRMSD(Catoms,altBindingAtoms);
            if (op.isGiven("t")) {
              cout << "simulatedRealRMSDorDOCKQ was: " << simulatedRealRMSDorDOCKQ << endl;
              cout << "altSimRealRMSD was: " << altSimRealRMSD << endl;
            }
            if (altSimRealRMSD < simulatedRealRMSDorDOCKQ) {
              simulatedRealRMSDorDOCKQ = altSimRealRMSD;
            }
            if (op.isGiven("t")) {
              cout << "lower RMSD picked for comparison was: " << simulatedRealRMSDorDOCKQ << endl;
            }
          }
        }
      }

      rmsdList.push_back(simulatedRealRMSDorDOCKQ);

      // if in testing mode, save the pdb files of the first 10 random things accepted, with their scores in their titles

      if (op.isGiven("t")) {
        if (t < 10) {
          W.writePDB(op.getString("t") + "position9_succeededWithContacts" + to_string(contactsCount) + to_string(t) + ".pdb");
        }
      }
    }

    // save the distribution of RMSDs as just a list of all the RMSDs, in a text file, comma separated - also compare to comparison RMSD (best correct vs model RMSD)

    cout << "# failed: " << fails << " # worked: " << d << endl;
  
  }

  else if (op.isGiven("cache")) { //instead of doing random docking, load a csv file of cached distance distributions
    string line;
    ifstream distFile;
    int lineCount = 0;
    distFile.open(op.getString("cache"));

   if(!distFile.is_open()) {
      perror("Error opening distribution file");
      exit(EXIT_FAILURE);
    }
    
    while(getline(distFile, line)) {
      if (lineCount < 3) { //skip over the saved previous distance evaluation & P-value
        lineCount++;
        continue;
      }
      rmsdList.push_back(MstUtils::toReal(line));
    }
  }

  mstreal comparisonRMSDorDOCKQ;

  if (op.isGiven("r")) {
    if (op.isGiven("t")) {
      mstreal OLDbRMSD = rc.rmsd(MBheavy.getAtoms(),realAtomsTB); // gets RMSD without alignment
      cout << "testing OLDbRMSD: " << OLDbRMSD << endl;
    }
    rc.align(MAheavy.getAtoms(),realAtomsTA,M); // aligns MAheavy with the original location of A, but transforms the whole structure M - if I'm doing this right lol
    if (op.isGiven("t")) {
      MAheavy.writePDB(op.getString("t") + "MAposition14.pdb");
      TACopy.writePDB(op.getString("t") + "TAposition14.pdb");
      M.writePDB(op.getString("t") + "Mposition14.pdb");
      Tcopy.writePDB(op.getString("t") + "CTposition14.pdb");
    }
    mstreal bRMSD = rc.rmsd(MBheavy.getAtoms(),realAtomsTB); // gets RMSD without alignment
    rc.align(MBheavy.getAtoms(),realAtomsTB,M); // aligns B with the original location of B, but transforms the whole structure W - if I'm doing this right lol
    if (op.isGiven("t")) {
      M.writePDB(op.getString("t") + "Mposition15.pdb");
      Tcopy.writePDB(op.getString("t") + "CTposition15.pdb");
    }
    mstreal aRMSD = rc.rmsd(MAheavy.getAtoms(),realAtomsTA); // gets RMSD without alignment
    comparisonRMSDorDOCKQ = (aRMSD + bRMSD);

    if (op.isGiven("t")) {
      cout << "testing aRMSD... " << aRMSD << endl;
      cout << "testing bRMSD... " << bRMSD << endl;
      cout << "testing comparisonRMSDorDOCKQ... " << comparisonRMSDorDOCKQ << endl;
    }

    if (op.isGiven("mbl")) {
      for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations

        AtomPointerVector altBindingAtoms = altBindingLocsAonly[abl];
        mstreal altcomparisonRMSDorDOCKQa = rc.rmsd(MAheavy.getAtoms(),altBindingAtoms); // gets RMSD without alignment

        rc.align(MAheavy.getAtoms(),altBindingAtoms,M); // aligns A with the alt original location of A, but transforms the whole structure W
        mstreal altcomparisonRMSDorDOCKQb = rc.rmsd(MBheavy.getAtoms(),realAtomsTB); // gets RMSD without alignment

        mstreal bestArmsdComp = min(altcomparisonRMSDorDOCKQa, aRMSD);
        mstreal bestBrmsdComp = min(altcomparisonRMSDorDOCKQb, bRMSD);

        mstreal altcomparisonRMSDorDOCKQ = (bestArmsdComp + bestBrmsdComp);
        if (op.isGiven("t")) {
          cout << "real comparisonRMSDorDOCKQ was: " << comparisonRMSDorDOCKQ << endl;
          cout << "real altcomparisonRMSDorDOCKQ was: " << altcomparisonRMSDorDOCKQ << endl;
        }
        if (altcomparisonRMSDorDOCKQ < comparisonRMSDorDOCKQ) {
          comparisonRMSDorDOCKQ = altcomparisonRMSDorDOCKQ;
        }
      }
    }

    
    if (op.isGiven("j")) {
      cout << "receptor + ligand RMSD between model and true structure is: " << to_string(comparisonRMSDorDOCKQ) << endl;
      exit(0);
    }
  }

  else if (op.isGiven("dq")) {

    cout << "getting DOCKQ between crystal structure and model..." << endl;

    comparisonRMSDorDOCKQ = getDOCKQ(&M, &MAheavy, MAbbAtoms, &MBheavy, MBbbAtoms, realAtomsTA, realAtomsCAbb, realAtomsTB, realAtomsCBbb, abResCons, conSize, a10interfaceIndexesBB, b10interfaceIndexesBB, ab10interface, rc, psDQ2);

      if (op.isGiven("mbl")) {
      for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) {

        mstreal altcomparisonRMSDorDOCKQ = getDOCKQ(&M, &MAheavy, MAbbAtoms, &MBheavy, MBbbAtoms, altBindingLocsAonly[abl], altBindingLocsAbbOnly[abl], realAtomsTB, realAtomsCBbb, abResConsAlts[abl], conSizeAlts[abl], a10interfaceIndexesBBAlts[abl], b10interfaceIndexesBBAlts[abl], ab10interfaceAlts[abl], rc, psDQ2); 
        if (op.isGiven("t")) {
          cout << "real comparisonRMSDorDOCKQ was: " << comparisonRMSDorDOCKQ << endl;
          cout << "real altcomparisonRMSDorDOCKQ was: " << altcomparisonRMSDorDOCKQ << endl;
        }
        if (altcomparisonRMSDorDOCKQ > comparisonRMSDorDOCKQ) {
          comparisonRMSDorDOCKQ = altcomparisonRMSDorDOCKQ;
        }
      }
    }

    if (op.isGiven("j")) {
      cout << "DOCKQ between model and true structure is: " << to_string(comparisonRMSDorDOCKQ) << endl;
      exit(0);
    }
  }

  else {
    AtomPointerVector Datoms = M.getAtoms();
    comparisonRMSDorDOCKQ = rc.bestRMSD(realAtoms,Datoms);

    if (op.isGiven("t")) {
      cout << "testing complex comparisonRMSDorDOCKQ... " << comparisonRMSDorDOCKQ << endl;
    }

    if (op.isGiven("mbl")) {
      for (int abl = 0; abl < altBindingLocs.size(); abl++) {
        AtomPointerVector altBindingAtoms = altBindingLocs[abl];
        mstreal altcomparisonRMSDorDOCKQ = rc.bestRMSD(altBindingAtoms,Datoms);

        if (op.isGiven("t")) {
          cout << "real comparisonRMSDorDOCKQ was: " << comparisonRMSDorDOCKQ << endl;
          cout << "real altcomparisonRMSDorDOCKQ was: " << altcomparisonRMSDorDOCKQ << endl;
        }
        if (altcomparisonRMSDorDOCKQ < comparisonRMSDorDOCKQ) {
          comparisonRMSDorDOCKQ = altcomparisonRMSDorDOCKQ;
        }
      }
    }

    if (op.isGiven("j")) {
      cout << "complex RMSD between model and true structure is: " << to_string(comparisonRMSDorDOCKQ) << endl;
      exit(0);
    }
  }

  mstreal pValueCount = 1;
  mstreal totCount = rmsdList.size();
  if (op.isGiven("dq")) {
    cout << "DOCKQ between correct structure & structure to evaluate was: " << comparisonRMSDorDOCKQ << endl;

    sort(rmsdList.begin(), rmsdList.end());
    std::reverse(rmsdList.begin(),rmsdList.end());
    for (int i = 0; i < rmsdList.size(); i++) {
      if (rmsdList[i] >= comparisonRMSDorDOCKQ) {
        pValueCount++;
      }
      else {
        break;
      }
    }
  }
  else {
    if (op.isGiven("r")) {
      cout << "R+L RMSD between correct structure & structure to evaluate was: " << comparisonRMSDorDOCKQ << endl;
    }
    else {
      cout << "complex RMSD between correct structure & structure to evaluate was: " << comparisonRMSDorDOCKQ << endl;
    }

    sort(rmsdList.begin(), rmsdList.end());
    for (int i = 0; i < rmsdList.size(); i++) {
      if (rmsdList[i] <= comparisonRMSDorDOCKQ) {
        pValueCount++;
      }
      else {
        break;
      }
    }
  }

  mstreal pValue = pValueCount/totCount;
  cout << "pValueCount is: " << pValueCount << endl;
  cout << "totCount is: " << totCount << endl;
  cout << "p value is: " << pValue << endl;

  string outFileName = op.getString("o");
  std::ofstream outFile(outFileName);
  outFile << comparisonRMSDorDOCKQ << "\n";
  outFile << pValue << "\n\n";

  for (int i = 0; i < rmsdList.size(); i++) {
    if (op.isGiven("t")) {
      cout << i << endl;
      cout << rmsdList[i] << endl;
    }
    outFile << rmsdList[i] << "\n";
  }
  outFile.close();
  return(0);
}