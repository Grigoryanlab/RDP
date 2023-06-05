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

// ***
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
  //AtomPointerVector s2evalAatoms = getHeavyAtoms(s2evalA.getAtoms());
  //AtomPointerVector s2evalBatoms = getHeavyAtoms(s2evalB.getAtoms());
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
  op.addOption("ca", "one half of the correctly docked structure, in pdb form, after having its backbone paired with that of the simulated docked structure via match_backbones.py; if using antibody mode, this should be the antibody ", true);
  op.addOption("da", "one half of the simulated docked structure, in pdb form, after having its backbone paired with that of the correctly docked structure via match_backbones.py (see -ca); if using antibody mode, this should be the antibody ", true);
  op.addOption("cb", "the other half of the correctly docked structure, in pdb form, after having its backbone paired with that of the simulated docked structure via match_backbones.py", true);
  op.addOption("db", "the other half of the simulated docked structure, in pdb form, after having its backbone paired with that of the correctly docked structure via match_backbones.py (see -cb); if using antibody mode, this should be the antibody ", true);
  op.addOption("n", "the number of valid docking positions to generate; defaults to 1 million, which should take around 1 hour per every 500 residues in the complex");
  op.addOption("i", "the number of valid interaction residues required to accept the structure; defaults to 1");
  op.addOption("cla", "the number of clashes allowed in an accepted structure; defaults to 3");
  op.addOption("sd", "the standard deviation of the normal distribution used to pull the docking partners apart from each other; defaults to 1.0; in angstroms");
  op.addOption("al", "an optional list of binding residues for the first docking partner (see -ca); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like 'A,100, ;A,100,A;A,100,B'. If you're only giving binding residues for one of the two partners, it must be this one (i.e. you cannot give bl without giving al)");
  op.addOption("bl", "an optional list of binding residues for the second docking partner (see -cb); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like A,100,;A,100,A;A,100,B.");
  op.addOption("abm", "an optional antibody / TCR mode, which will treat all the loops of the antibody as the binding residues; the antibody structure must use IMGT numbering and its structure must be entered using the -ca argument not the -cb argument.");
  op.addOption("cdr3", "limit the antibody / TCR mode to using CDR3 loop contacts only");
  op.addOption("q", "an optional quick mode for the --al or --abm flags, wherein the docking distribution is skewed towards conformations involving the binding residues given on the A side, to make the calculation faster");
  op.addOption("limA", "Must be used with al or abm, and bl. Optionally limit the range of angles of accepted dockings, wherein each binding partner has a line drawn between its geometric center and the geometric center of its binding residues; one line is used as an axis, and the other is used to calculate the angle of tilt relative to that axis. Primarily useful for TCRpMHC random dockings.");
  op.addOption("cache", "Use a cached distribution of underlying distance metrics instead; supply the path to that csv file here.");
  op.addOption("j", "Just print the underlying distance metric between the two structures instead of doing comparisons.");
  op.addOption("o", "the output file to save the distribution of RMSDs, as a comma separated list", true);
  op.addOption("t", "create testing files; provide the output directory to save the testing files. This also sets the number of dockings to just 10 so you aren't inundated with files on accident :)");
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

  if (op.isGiven("q") && (!(op.isGiven("al")) && !(op.isGiven("abm")))) {
    cout << "you cannot give q without giving al or abm";
    exit(1);
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
    clashesAllowed = 3;
  }

  mstreal normalDistBase = 1.0;

  // load backbone or full-atom structures

  Structure CA00(op.getString("ca"));
  AtomPointerVector CA0heavy = getHeavyAtoms(CA00.getAtoms());
  Structure CA0(CA0heavy);

  Structure CB00(op.getString("cb"));
  AtomPointerVector CB0heavy = getHeavyAtoms(CB00.getAtoms());
  Structure CB0(CB0heavy);

  Structure DA0(op.getString("da"));
  AtomPointerVector DAheavy = getHeavyAtoms(DA0.getAtoms());
  Structure DA(DAheavy);

  Structure DB0(op.getString("db"));
  AtomPointerVector DBheavy = getHeavyAtoms(DB0.getAtoms());
  Structure DB(DBheavy);

  // make a full structure for C0, which is connected to CA0 / CB0 rather than a copy of them
  Structure C0;
  for (int i = 0; i < CA0.chainSize(); i++) {
    C0.appendChain(&CA0.getChain(i));
  }
  for (int i = 0; i < CB0.chainSize(); i++) {
    C0.appendChain(&CB0.getChain(i));
  }

  Structure CC(C0); // a copy for doing RMSD calculations with
  AtomPointerVector realAtoms = getHeavyAtoms(CC.getAtoms());
  Structure CAC(CA0);
  //AtomPointerVector realAtomsCA = getHeavyAtoms(CAC.getAtoms());
  AtomPointerVector realAtomsCA = CAC.getAtoms();
  Structure CBC(CB0);
  //AtomPointerVector realAtomsCB = getHeavyAtoms(CBC.getAtoms());
  AtomPointerVector realAtomsCB = CBC.getAtoms();

  ProximitySearch psA(realAtomsCA, 15);
  map<int,vector<int>> abResCons = getABcontactDict(realAtomsCA,realAtomsCB,5.0,psA);

  int conSize = 0;
  for (auto it = abResCons.cbegin(); it != abResCons.cend(); ++it) {
    conSize += it->second.size();
  }

  // get the indexes and then interface structure for A & B, by a 10 Angstrom cutoff, in case DOCKQ is being used as an underlying distance metric

  vector<vector<int>> ab10realIndexes = getABindexes(realAtomsCA,realAtomsCB,10.0);
  vector<int> a10interfaceIndexesAll = ab10realIndexes[0];
  vector<int> b10interfaceIndexesAll = ab10realIndexes[1];
  vector<int> a10interfaceIndexesBB;
  vector<int> b10interfaceIndexesBB;
  AtomPointerVector ab10interface;
  for (int i = 0; i < a10interfaceIndexesAll.size(); i++) {
    Atom* checkingAtom = realAtomsCA[a10interfaceIndexesAll[i]];
    if ((checkingAtom->getName() == "CA") or (checkingAtom->getName() == "C") or (checkingAtom->getName() == "N") or (checkingAtom->getName() == "O")) {
      ab10interface.push_back(checkingAtom);
      a10interfaceIndexesBB.push_back(a10interfaceIndexesAll[i]);
    }
  }

  for (int i = 0; i < b10interfaceIndexesAll.size(); i++) {
    Atom* checkingAtom = realAtomsCB[b10interfaceIndexesAll[i]];
    if ((checkingAtom->getName() == "CA") or (checkingAtom->getName() == "C") or (checkingAtom->getName() == "N") or (checkingAtom->getName() == "O")) {
      ab10interface.push_back(checkingAtom);
      b10interfaceIndexesBB.push_back(b10interfaceIndexesAll[i]);
    }
  }

  // make a full structure for D, which is connected to DA / DB rather than a copy of them
  Structure D;
  for (int i = 0; i < DA.chainSize(); i++) {
    D.appendChain(&DA.getChain(i));
  }
  for (int i = 0; i < DB.chainSize(); i++) {
    D.appendChain(&DB.getChain(i));
  }

  // if you're using model-structure partners for docking, make CA & CB out of them; otherwise stick with the crystal partners (CA0 & CB0) which is default

  Structure CA;
  Structure CB;

  vector<AtomPointerVector> altBindingLocs;
  vector<AtomPointerVector> altBindingLocsBB;
  vector<AtomPointerVector> altBindingLocsAonly;
  vector<AtomPointerVector> altBindingLocsAbbOnly;
  Structure CALT;
  Structure CBCOPY;

  if (op.isGiven("mbl")) {
    vector<string> altBindingAstructs = MstUtils::split(op.getString("mbl"), ",");
    for (int i = 0; i < altBindingAstructs.size(); i++) {
      CBCOPY = Structure(CB0);
      string caltPath = altBindingAstructs[i];
      CALT = Structure(caltPath);
      AtomPointerVector caAltAtoms = getHeavyAtoms(CALT.getAtoms());

      altBindingLocsAonly.push_back(caAltAtoms);
      altBindingLocsAbbOnly.push_back(getBBatoms(caAltAtoms));

      for (int ii = 0; ii < CBCOPY.chainSize(); ii++) {
        CALT.appendChain(&CBCOPY.getChain(ii));
      }
      AtomPointerVector cAltAtoms = getHeavyAtoms(CALT.getAtoms());
      
      altBindingLocs.push_back(cAltAtoms);
      altBindingLocsBB.push_back(getBBatoms(cAltAtoms));
    }
  }

  // if alt binding locations: make the indexes for the interface, by a 5 angstrom cutoff, in case DOCKQ is being used as an underlying distance metric 

   vector<map<int,vector<int>>> abResConsAlts;
   vector<int> conSizeAlts;
   for (int i = 0; i < altBindingLocsAonly.size(); i++) {
    ProximitySearch psAlt(altBindingLocsAonly[i], 15);

    map<int,vector<int>> abResConsAlt = getABcontactDict(altBindingLocsAonly[i],realAtomsCB,5.0,psAlt);

    int conSizeAlt = 0;
    for (auto it = abResConsAlt.cbegin(); it != abResConsAlt.cend(); ++it) {
      conSizeAlt += it->second.size();
    }

    abResConsAlts.push_back(abResConsAlt);
    conSizeAlts.push_back(conSizeAlt);
   }

   //plus if alt binding locations, get the indexes / interface Structure for a 10 angstrom cutoff - also in case DOCKQ is used ***
  vector<vector<int>> a10interfaceIndexesBBAlts;
  vector<vector<int>> b10interfaceIndexesBBAlts;
  vector<AtomPointerVector> ab10interfaceAlts;
  for (int i = 0; i < altBindingLocsAonly.size(); i++) {
    vector<vector<int>> ab10realIndexesAlt = getABindexes(altBindingLocsAonly[i],realAtomsCB,10.0);
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
      Atom* checkingAtom = realAtomsCB[b10interfaceIndexesAllAlt[ii]];
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
    CA = Structure(DA); // this duplicates / does not connect CA & DA
    CB = Structure(DB);
  }
  else {
    CA = Structure(CA0);
    CB = Structure(CB0);
  }

  // make a full structure for C, which is connected to CA / CB rather than a copy of them
  Structure C;
  for (int i = 0; i < CA.chainSize(); i++) {
    C.appendChain(&CA.getChain(i));
  }
  for (int i = 0; i < CB.chainSize(); i++) {
    C.appendChain(&CB.getChain(i));
  }

  //lilWriteTest("/dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/dockingDistDebugging.txt","testing 3?");

  vector <Residue*> CAreses = CA.getResidues();
  //AtomPointerVector CAatoms = getHeavyAtoms(CA.getAtoms());
  AtomPointerVector CAatoms = CA.getAtoms();

  vector <Residue*> CBreses = CB.getResidues();
  //AtomPointerVector CBatoms = getHeavyAtoms(CB.getAtoms());
  AtomPointerVector CBatoms = CB.getAtoms();

  vector <Residue*> DAreses = DA.getResidues();
  //AtomPointerVector DAatoms = getHeavyAtoms(DA.getAtoms());
  AtomPointerVector DAatoms = DA.getAtoms();

  // for DOCKQ to get fnat later on

  ProximitySearch psDQ2(DAatoms, 15);

  vector <Residue*> DBreses = DB.getResidues();
  //AtomPointerVector DBatoms = getHeavyAtoms(DB.getAtoms());
  AtomPointerVector DBatoms = DB.getAtoms();

  AtomPointerVector CAbbAtoms = getBBatoms(CAatoms);
  AtomPointerVector CBbbAtoms = getBBatoms(CBatoms);
  AtomPointerVector DAbbAtoms = getBBatoms(DAatoms);
  AtomPointerVector DBbbAtoms = getBBatoms(DBatoms);
  AtomPointerVector realAtomsCAbb = getBBatoms(realAtomsCA);
  AtomPointerVector realAtomsCBbb = getBBatoms(realAtomsCB);

  RMSDCalculator rc;
  vector <mstreal> rmsdList {};

  if ((!op.isGiven("cache")) && (!op.isGiven("j"))) {

    // position 1
    if (op.isGiven("t")) {
      C.writePDB(op.getString("t") + "position1" + ".pdb");
    }

    // translate each to the origin

    Transform TCA0 = TransformFactory::translate(-CAatoms.getGeometricCenter());
    TCA0.apply(CA);

    Transform TCB0 = TransformFactory::translate(-CBatoms.getGeometricCenter());
    TCB0.apply(CB);

    // position 2
    if (op.isGiven("t")) {
      C.writePDB(op.getString("t") + "position2" + ".pdb");
    }

    // if binding residues are given for A (or A has binding residues by virtue of this being run in antibody-antigen binding mode) get them; same for B

    AtomPointerVector CAbinderAtoms;
    AtomPointerVector CBbinderAtoms;

    vector<tuple <string, int, char>> aResesDetails; // stores details for binding residues for A
    vector<tuple <string, int, char>> bResesDetails; // stores details for binding residues for B

    Structure aBinderS;
    if (op.isGiven("abm")) {
      for (int i = 0; i < CAreses.size(); i++) {
        Residue* currR = CAreses[i];
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
          aBinderS.addResidue(currR);
          AtomPointerVector cResAtoms = currR->getAtoms();
          for (int ca = 0; ca < cResAtoms.size(); ca++) {
            Atom* currA = cResAtoms[ca];
            string cResAtomName = currA->getName();
            if (!(op.isGiven("dq")) && (cResAtomName != "CA")) { // if not DOCKQ mode, and not a CA atom...
              continue;
            } // only get indexes for CAs, as contacts will be defined as inter-CA-distance of 10 angstroms or less, unless using DOCKQ mode in which case it's all atoms 5 angstroms or less
            CAbinderAtoms.push_back(currA);
          }
        }
      }
      if (op.isGiven("q")) {
        //AtomPointerVector aBinderAtoms = getHeavyAtoms(aBinderS.getAtoms());
        AtomPointerVector aBinderAtoms = aBinderS.getAtoms();
        CartesianPoint geoCenterA = aBinderAtoms.getGeometricCenter();
        Transform TZ = TransformFactory::alignVectorWithXAxis(geoCenterA);
        TZ.apply(CA);

        if (op.isGiven("t")) {
          C.writePDB(op.getString("t") + "position3" + ".pdb");
        }
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

      for (int i = 0; i < CAreses.size(); i++) {
        Residue* currR = CAreses[i];
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
            aBinderS.addResidue(currR);
            AtomPointerVector cResAtoms = currR->getAtoms();
            for (int ca = 0; ca < cResAtoms.size(); ca++) {
              Atom* currA = cResAtoms[ca];
              string cResAtomName = currA->getName();
              if (!(op.isGiven("dq")) && (cResAtomName != "CA")) { // if not DOCKQ mode, and not a CA atom...
                continue;
              }
              CAbinderAtoms.push_back(currA);
            }
          }
        }
      }
      if (op.isGiven("q")) {
        AtomPointerVector aBinderAtoms = aBinderS.getAtoms();
        CartesianPoint geoCenterA = aBinderAtoms.getGeometricCenter();
        Transform TZ = TransformFactory::alignVectorWithXAxis(geoCenterA);
        TZ.apply(CA);

        if (op.isGiven("t")) {
          C.writePDB(op.getString("t") + "position3" + ".pdb");
        }
      }
    }

    // if binding residues are given for B, pre-store the details to iterate over easily

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

      for (int i = 0; i < CBreses.size(); i++) {
        Residue* currR = CBreses[i];
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
            AtomPointerVector cResAtoms = currR->getAtoms();
            for (int ca = 0; ca < cResAtoms.size(); ca++) {
              string cResAtomName = cResAtoms[ca]->getName();
              if (cResAtomName == "CA") {
                int cResIndx = cResAtoms[ca]->getIndex();
                CBbinderAtoms.push_back(cResAtoms[ca]);
              }
            }
          }
        }
      }
    }

    vector<int> bBinderIndexes; // store the atom indexes of CA binding atoms
    for (int i = 0; i < CBatoms.size(); i++) { // go over all atoms in CB...
      Atom* currAtom = CBatoms[i];
      string currAtomName = currAtom->getName();
      if (!(op.isGiven("dq")) && (currAtomName != "CA")) { // if not DOCKQ mode, and not a CA atom...
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

    // let's get started on the random docking now! Set up variables (d = number of random dockings required, t = a variable used to count the initial 10 structures to be saved for testing mode)

    int d = 0;
    int t = 0;
    int t2 = 0;
    //int t3 = 0;
    int fails = 0;

    // set up important vectors...

    // keep a list of indexes of previous clashes
    vector <int> recentClashes;
    // and have a list of current clashes
    vector <int> currClashes;
    // as well as a list of indexes possible to be involved in contacts (CA atoms)
    vector <int> aPossibleContacts;
    vector <int> bPossibleContacts;

    // set the curr coords of CA to be alt cords, so that they can be reset every new random docking

    for (int a = 0; a < CAatoms.size(); a++) {
      Atom* currAtom = CAatoms[a];
      currAtom->clearAlternatives();
      currAtom->addAlternative(currAtom);
    }

    // set the current coords of CB to be alternate coords, so that they can be reset every new random docking - also fill out which indexes are possible to make contacts in B, and get their atoms

    for (int a = 0; a < CBatoms.size(); a++) {
      Atom* currAtom = CBatoms[a];
      currAtom->clearAlternatives();
      currAtom->addAlternative(currAtom);
      // as well as get indexes for the potential contact residues, if not already gotten from abm / al / bl flags being given
      if (op.isGiven("bl")) {
        continue;
      }
      else {
        string currName = currAtom->getName();
        if (op.isGiven("dq")) {
          bPossibleContacts.push_back(a);
          CBbinderAtoms.push_back(currAtom);
        }
        if (!(op.isGiven("dq")) && (currName == "CA")) {
          bPossibleContacts.push_back(a);
          CBbinderAtoms.push_back(currAtom);
        }
      }
    }

    // and fill out which ones are possible to make contacts in A / get their atoms

    for (int a = 0; a < CAatoms.size(); a++) {
      Atom* currAtom = CAatoms[a];
      if (op.isGiven("al") || op.isGiven("abm")) {
        continue;
      }
      else {
        string currName = currAtom->getName();
        if (!(op.isGiven("dq")) && (currName != "CA")) { // if not DOCKQ mode, and not a CA atom...
          continue;
        }
        aPossibleContacts.push_back(a);
        CAbinderAtoms.push_back(currAtom);
      }
    }

    if (op.isGiven("bl")) {
      bPossibleContacts = bBinderIndexes;
    }

    // make structures out of the A binder atoms & B binder atoms, for using with calculateExtent below to get Y & Z dimensions to scale the Y / Z translations by

    //Structure CAbinderS(CAbinderAtoms);
    //Structure CBbinderS(CBbinderAtoms);

    // set up xyz high & low values to be used later in the loop, for checking xLow & xHigh when pulling out proteins

    mstreal xLow;
    mstreal xHigh;
    mstreal yLow;
    mstreal yHigh;
    mstreal zLow;
    mstreal zHigh;

    // and variables for calculating the Y & Z extents for the Y & Z translations

    mstreal ayExtent;
    mstreal azExtent;
    mstreal axLow;
    mstreal axHigh;
    mstreal ayLow;
    mstreal ayHigh;
    mstreal azLow;
    mstreal azHigh;
    mstreal ayLen;
    mstreal azLen;

    mstreal byExtent;
    mstreal bzExtent;
    mstreal bxLow;
    mstreal bxHigh;
    mstreal byLow;
    mstreal byHigh;
    mstreal bzLow;
    mstreal bzHigh;
    mstreal byLen;
    mstreal bzLen;

    //& docking requirements

    mstreal clashDistance = 3.0;
    
    mstreal contactDistance;
    if (op.isGiven("dq")) { // if using all atoms
      contactDistance = 5.0;
    }
    else { // if using CA atoms
      contactDistance = 8.0;
    }

    mstreal angleOff;

    while (d < numberDockingsRequired) { // while the number of docked structures is still insufficient...

      // reset to alt coors

      for (int i = 0; i < CAreses.size(); i++) {
        CAreses[i]->makeAlternativeMain(0); // reset to origin position
      }

      for (int i = 0; i < CBreses.size(); i++) {
        CBreses[i]->makeAlternativeMain(0); // reset to origin position
      }

      //position 4

      if ((op.isGiven("t")) && (t == 0)) {
        C.writePDB(op.getString("t") + "position4" + ".pdb");
      }

      // rotate CA randomly along the x, y, and z axes, if no binding residues given for A, or limited to prevent A's binding residues from facing too far away from B if binding residues specified

      if (!op.isGiven("q")) {
        mstreal xaAngle = MstUtils::randUnit(0,360);
        Transform TXA = TransformFactory::rotateAroundX(xaAngle);
        TXA.apply(CA);
        mstreal yaAngle = MstUtils::randUnit(0,360);
        Transform TYA = TransformFactory::rotateAroundY(yaAngle);
        TYA.apply(CA);
        mstreal zaAngle = MstUtils::randUnit(0,360);
        Transform TZA = TransformFactory::rotateAroundZ(zaAngle);
        TZA.apply(CA);

        //position 5
        if ((op.isGiven("t")) && (t == 0)) {
          C.writePDB(op.getString("t") + "position5" + ".pdb");
        }
      }

      // rotate CB randomly along the x, y, and z axes

      mstreal xbAngle = MstUtils::randUnit(0,360);
      Transform TXB = TransformFactory::rotateAroundX(xbAngle);
      TXB.apply(CB);
      mstreal ybAngle = MstUtils::randUnit(0,360);
      Transform TYB = TransformFactory::rotateAroundY(ybAngle);
      TYB.apply(CB);
      mstreal zbAngle = MstUtils::randUnit(0,360);
      Transform TZB = TransformFactory::rotateAroundZ(zbAngle);
      TZB.apply(CB);

      //position 6
      if ((op.isGiven("t")) && (t == 0)) {
        C.writePDB(op.getString("t") + "position6" + ".pdb");
      }

      // make standard deviations for random translations in the Y & Z directions, based on the Y & Z lengths of CB - if binding residues are given for a

      if (((op.isGiven("al")) || (op.isGiven("abm"))) && op.isGiven("q")) {

        ProximitySearch::calculateExtent(CAbinderAtoms,axLow,ayLow,azLow,axHigh,ayHigh,azHigh);
        ayExtent = ayHigh - ayLow;
        azExtent = azHigh - azLow;

        ProximitySearch::calculateExtent(CBatoms,bxLow,byLow,bzLow,bxHigh,byHigh,bzHigh);
        byExtent = byHigh - byLow;
        bzExtent = bzHigh - bzLow;

        //ProximitySearch::calculateExtent(CBbinderAtoms,bxLow,byLow,bzLow,bxHigh,byHigh,bzHigh);
        //byExtent = byHigh - byLow;
        //bzExtent = bzHigh - bzLow;

        //mstreal yLowExtent = std::min(byExtent, ayExtent);
        //mstreal zLowExtent = std::min(bzExtent, azExtent);
        //mstreal ySig = (1.0 / 6.0) * yLowExtent;
        //mstreal zSig = (1.0 / 6.0) * zLowExtent;


        //mstreal yRand = MstUtils::randNormal(0,yLowExtent);
        //mstreal zRand = MstUtils::randNormal(0,zLowExtent);
        
        //mstreal yRand = MstUtils::randNormal(0,(1.0 / 6.0) * (ayExtent + byExtent));
        //mstreal zRand = MstUtils::randNormal(0,(1.0 / 6.0) * (azExtent + bzExtent));
        
        mstreal yRand = MstUtils::randUnit(0,(1.0 / 2.0) * (ayExtent + byExtent));
        mstreal zRand = MstUtils::randUnit(0,(1.0 / 2.0) * (azExtent + bzExtent));


        // use those to randomly translate CB in the Y & Z axes so that the angle of binding can vary, with a normal distribution with standard deviation of 1/6 the distance between the furthest Y points / Z points (so 1/3 the distance from the middle and furthest Y / Z point, which means 0.13% will be beyond the range of being able to connect, which is fine)

        for (int a = 0; a < CBatoms.size(); a++) {
          mstreal yStart = CBatoms[a]->getY();
          CBatoms[a]->setY(yStart + yRand);
          mstreal zStart = CBatoms[a]->getZ();
          CBatoms[a]->setZ(zStart + zRand);
        }

        //position 7
        if ((op.isGiven("t")) && (t == 0)) {
          C.writePDB(op.getString("t") + "position7" + ".pdb");
        }

        // if there are B binder atoms, check if they're outside the Y or Z ranges of the A binder atoms, and so could never meet; go back to start without attempting to dock if so

        if (op.isGiven("bl")) {

          ProximitySearch::calculateExtent(CBbinderAtoms,bxLow,byLow,bzLow,bxHigh,byHigh,bzHigh);
          
          if ((byLow > ayHigh) || (byHigh < ayLow) || (bzLow > azHigh) || (bzHigh < azLow)) {
            continue;
          }
        }

        // if limiting binding angles between binding partners, check what the binding angle is, and go back to start if outside the limits 

        if (op.isGiven("limA")) {
          CartesianPoint geoCenterB = CBatoms.getGeometricCenter();
          CartesianPoint geoCenterBbinders = CBbinderAtoms.getGeometricCenter();
          CartesianPoint geoCenterBprime(geoCenterB.getX() - 10, geoCenterB.getY(), geoCenterB.getZ());
          angleOff = CartesianGeometry::angle(geoCenterBprime,geoCenterB,geoCenterBbinders);
          if (angleOff > stoi(op.getString("limA"))) {
            continue;
          }
        }
      }

      // set up proximity search with A, since B will be the partner moving (if you move an object you need to re-make a prox-search or it won't work / will segfault)

      ProximitySearch ps;

      if (op.isGiven("dq")) {
        ps = ProximitySearch(CAbbAtoms, 15);
      }
      else {
        ps = ProximitySearch(CAatoms, 15);
      }

      // for DOCKQ
      ProximitySearch psDQ1(CAatoms, 15);

      // and a proximity search with just A's CA binding residues - this is a ProximitySearch*

      ProximitySearch ps2(CAbinderAtoms, 15);

      // while a proper docked stage has not been reached, but there are still opporunitites for success

      bool absoluteSuccess = false;
      bool absoluteFailure = false;

      /// reset recent clashes, if left-over from the last round
      recentClashes.resize(0);

      mstreal contactsCount;
      mstreal xRand;

      while ((absoluteSuccess == false) && (absoluteFailure == false)) {

        // reset current clashes to be empty
        currClashes.resize(0);

        // pull B atoms out along the X axis, in small steps with a normal distrubtion of ~1.0 angstrom and mew of 0, except always positive so it's pulling away from A - can try variations on the 1.0 and adjust

        xRand = fabs(MstUtils::randNormal(0,normalDistBase));

        for (int a = 0; a < CBatoms.size(); a++) {
          mstreal xStart = CBatoms[a]->getX();
          CBatoms[a]->setX(xStart + xRand);
        }

        //position 8
        if (op.isGiven("t") && (t == 0)) {
          C.writePDB(op.getString("t") + "position8" + ".pdb");
        }

        // reject / continue depending on if there are too many clashes or not

        int clashCount = 0;
        currClashes.resize(0);

        // first check positions in B that previously had clashes...

        for (int a = 0; a < recentClashes.size(); a++) {
          int currBcheck = recentClashes[a];
          Atom* currBatom = CBatoms[currBcheck];
          CartesianPoint currBcoor = currBatom->getCoor();
          vector <int> currClashPts = ps.getPointsWithin(currBcoor, 0, clashDistance, false);
          int currBsize = currClashPts.size();
          if (currBsize > 0) {
            currClashes.push_back(currBcheck);
            clashCount+=currBsize;
          }
          if (clashCount > clashesAllowed) {
            break;
          }
        }

        if (clashCount > clashesAllowed) {
          recentClashes = currClashes;
          continue;
        }

        // then if you're not over the clash limit, check the rest of the valid residues

        for (int a = 0; a < CBatoms.size(); a++) {
          int currBcheck = a;
          if (std::find(recentClashes.begin(), recentClashes.end(), currBcheck) != recentClashes.end()) {
            continue; // don't double-count if already checked as part of recent clashes
          }

          Atom* currBatom = CBatoms[currBcheck];
          CartesianPoint currBcoor = currBatom->getCoor();
          vector <int> currClashPts = ps.getPointsWithin(currBcoor, 0, clashDistance, false);
          int currBsize = currClashPts.size();
          if (currBsize > 0) {
            currClashes.push_back(currBcheck);
            clashCount+=currBsize;
          }
          if (clashCount > clashesAllowed) {
            break;
          }
        }

        if (clashCount > clashesAllowed) {
          recentClashes = currClashes;
          continue;
        }

        //position 9
        if (op.isGiven("t") && (t == 0)) {
          C.writePDB(op.getString("t") + "position9" + ".pdb");
        }

        // we've only reached this point if there are an aceeptable number of clashes, so accept if there are at least the required number of contacts - possibly limited by antibody binding mode, a binding residues, and/or b binding residues

        contactsCount = 0;
        vector<pair<Residue*,Residue*>> pastContacts = {};

        for (int a = 0; a < bPossibleContacts.size(); a++) { // go over all CA atoms in B possible to be part of contacts
          int currBcheck = bPossibleContacts[a];
          Residue* currBres = CBatoms[currBcheck]->getParent();
          CartesianPoint currBcoor = CBatoms[currBcheck]->getCoor();
          vector <int> currContactPts = ps2.getPointsWithin(currBcoor, 0.0, contactDistance, false); // ps2 is only checking for CA atoms in A possible to be part of contacts
          if (op.isGiven("dq")) {
            for (int aa = 0; aa < currContactPts.size(); aa++) {
            Residue* currAres = CAbinderAtoms[currContactPts[aa]]->getParent();
            if (std::find(pastContacts.begin(), pastContacts.end(), make_pair(currBres,currAres)) != pastContacts.end()) {
              continue;
            }
            else {
              pastContacts.push_back(make_pair(currBres,currAres));
              contactsCount += 1;
            }
          }
          }
          else {
            contactsCount += currContactPts.size();
          }
          if (contactsCount >= contactsRequired) {
            absoluteSuccess = true;
            break;
          }
        }

        //check if all atoms in CB are further positive than all atoms in CA by over 8 angstroms... b/c if so the docking is irrecoverable

        ProximitySearch::calculateExtent(CB,xLow,yLow,zLow,xHigh,yHigh,zHigh);
        mstreal bX = xLow;
        ProximitySearch::calculateExtent(CA,xLow,yLow,zLow,xHigh,yHigh,zHigh);
        mstreal aX = xHigh;

        if (bX - aX > 8.0) {
          absoluteFailure = true;
        }
      }

      if (absoluteFailure) {
        fails++;
        if (op.isGiven("t") && (t2 == 0)) {
          C.writePDB(op.getString("t") + "position11_failedWithContacts" + to_string(contactsCount) + ".pdb");
          t2++;
        }
        continue;
      }
      else {
        d++; // accepted! :D
        if (op.isGiven("t") && (t2 == 0)) {
          C.writePDB(op.getString("t") + "position11_succeededWithContacts" + to_string(contactsCount) + ".pdb");
          t2++;
        }
      }

      if (d%10000 == 0) {
        cout << d << " total docks accepted..." << endl;
      }

      // for those accepted, compare to the correct structure to calculate the RMSD

      mstreal simulatedRealRMSDorDOCKQ;

      if (op.isGiven("r")) {
        rc.align(CA.getAtoms(),realAtomsCA,C);
        if (op.isGiven("t") && (t == 0)) {
          C.writePDB(op.getString("t") + "Cposition12.pdb");
          CAC.writePDB(op.getString("t") + "CACposition12.pdb");
        }
        mstreal bRMSD = rc.rmsd(CB.getAtoms(),realAtomsCB); // gets RMSD without alignment
        rc.align(CB.getAtoms(),realAtomsCB,C); // aligns CB with the original location of CB, but transforms the whole structure C
        if (op.isGiven("t") && (t == 0)) {
          C.writePDB(op.getString("t") + "Cposition13.pdb");
          CBC.writePDB(op.getString("t") + "CBCposition13.pdb");
        }
        mstreal aRMSD = rc.rmsd(CA.getAtoms(),realAtomsCA); // gets RMSD without alignment
        simulatedRealRMSDorDOCKQ = (aRMSD + bRMSD);

        if (op.isGiven("mbl")) {
          for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations

            AtomPointerVector altBindingAtoms = altBindingLocsAonly[abl];
            mstreal aaRMSD = rc.rmsd(CA.getAtoms(),altBindingAtoms); // gets RMSD without alignment - because C's already aligned by origianl location of CB

            rc.align(CA.getAtoms(),altBindingAtoms,C); // aligns CA with the alt original location of CA, but transforms the whole structure C
            mstreal abRMSD = rc.rmsd(CB.getAtoms(),realAtomsCB); // gets RMSD without alignment

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
          cout << "testing CAbbAtoms[0]" << CAbbAtoms[0] << endl;
          cout << "testing getBBatoms(CA.getAtoms())[0]" << getBBatoms(CA.getAtoms())[0] << endl;
        }

        mstreal currDOCKQ = getDOCKQ(&C, &CA, CAbbAtoms, &CB, CBbbAtoms, realAtomsCA, realAtomsCAbb, realAtomsCB, realAtomsCBbb, abResCons, conSize, a10interfaceIndexesBB, b10interfaceIndexesBB, ab10interface, rc, psDQ1);
        simulatedRealRMSDorDOCKQ = currDOCKQ;

        if (op.isGiven("t")) {
          cout << "complex DOCKQ was: " << currDOCKQ << endl;
        }
        if (op.isGiven("mbl")) {
          for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations

            mstreal altDOCKQ = getDOCKQ(&C, &CA, CAbbAtoms, &CB, CBbbAtoms, altBindingLocsAonly[abl], altBindingLocsAbbOnly[abl], realAtomsCB, realAtomsCBbb, abResConsAlts[abl], conSizeAlts[abl], a10interfaceIndexesBBAlts[abl], b10interfaceIndexesBBAlts[abl], ab10interfaceAlts[abl], rc, psDQ1);
            if (op.isGiven("t")) {
              cout << "altDOCKQ was: " << altDOCKQ << endl;
              cout << "altBindingLocsAonly[abl][0] was: " << altBindingLocsAonly[abl][0] << endl;
              cout << "realAtomsCA[0] was: " << realAtomsCA[0] << endl;
              cout << "*realAtomsCA[0] was: " << *realAtomsCA[0] << endl;
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
        AtomPointerVector Catoms = C.getAtoms();
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
          C.writePDB(op.getString("t") + "position11_succeededWithContacts" + to_string(contactsCount) + to_string(t) + ".pdb");
          //cout << "angleOff was: " << endl;
          //cout << angleOff << endl;
          t++;
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
      mstreal OLDbRMSD = rc.rmsd(DB.getAtoms(),realAtomsCB); // gets RMSD without alignment
      cout << "testing OLDbRMSD: " << OLDbRMSD << endl;
    }
    rc.align(DA.getAtoms(),realAtomsCA,D); // aligns DA with the original location of CA, but transforms the whole structure D - if I'm doing this right lol
    if (op.isGiven("t")) {
      DA.writePDB(op.getString("t") + "DAposition14.pdb");
      CAC.writePDB(op.getString("t") + "CACposition14.pdb");
      D.writePDB(op.getString("t") + "Dposition14.pdb");
      CC.writePDB(op.getString("t") + "CCposition14.pdb");
    }
    mstreal bRMSD = rc.rmsd(DB.getAtoms(),realAtomsCB); // gets RMSD without alignment
    rc.align(DB.getAtoms(),realAtomsCB,D); // aligns CB with the original location of CB, but transforms the whole structure C - if I'm doing this right lol
    if (op.isGiven("t")) {
      D.writePDB(op.getString("t") + "Dposition15.pdb");
      CC.writePDB(op.getString("t") + "CCposition15.pdb");
    }
    mstreal aRMSD = rc.rmsd(DA.getAtoms(),realAtomsCA); // gets RMSD without alignment
    comparisonRMSDorDOCKQ = (aRMSD + bRMSD);

    if (op.isGiven("t")) {
      cout << "testing aRMSD... " << aRMSD << endl;
      cout << "testing bRMSD... " << bRMSD << endl;
      cout << "testing comparisonRMSDorDOCKQ... " << comparisonRMSDorDOCKQ << endl;
    }

    if (op.isGiven("mbl")) {
      for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) { // altBindingLocs is a list of atom lists representing alternative complex conformations

        AtomPointerVector altBindingAtoms = altBindingLocsAonly[abl];
        mstreal altcomparisonRMSDorDOCKQa = rc.rmsd(DA.getAtoms(),altBindingAtoms); // gets RMSD without alignment

        rc.align(DA.getAtoms(),altBindingAtoms,D); // aligns CA with the alt original location of CA, but transforms the whole structure C
        mstreal altcomparisonRMSDorDOCKQb = rc.rmsd(DB.getAtoms(),realAtomsCB); // gets RMSD without alignment

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

    comparisonRMSDorDOCKQ = getDOCKQ(&D, &DA, DAbbAtoms, &DB, DBbbAtoms, realAtomsCA, realAtomsCAbb, realAtomsCB, realAtomsCBbb, abResCons, conSize, a10interfaceIndexesBB, b10interfaceIndexesBB, ab10interface, rc, psDQ2);

      if (op.isGiven("mbl")) {
      for (int abl = 0; abl < altBindingLocsAonly.size(); abl++) {

        mstreal altcomparisonRMSDorDOCKQ = getDOCKQ(&D, &DA, DAbbAtoms, &DB, DBbbAtoms, altBindingLocsAonly[abl], altBindingLocsAbbOnly[abl], realAtomsCB, realAtomsCBbb, abResConsAlts[abl], conSizeAlts[abl], a10interfaceIndexesBBAlts[abl], b10interfaceIndexesBBAlts[abl], ab10interfaceAlts[abl], rc, psDQ2); 
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
    AtomPointerVector Datoms = D.getAtoms();
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
    cout << "RMSD between correct structure & structure to evaluate was: " << comparisonRMSDorDOCKQ << endl;

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

  cout << "testing" << endl;
  return(0);
}

//example command: ./dockingDistribution --ca "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAb_matchedbb/crystal/6yu8_bef6732d-1f2b-41e2-960b-d5c68c4250b2_0.9186929628704293.pdb" --cb "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAg_matchedbb/crystal/6yu8_4acaabcc-e98e-4ddc-a8c1-cece355b44c3_2.555964177325445.pdb" --da "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAb_matchedbb/af/6yu8_bef6732d-1f2b-41e2-960b-d5c68c4250b2_0.9186929628704293.pdb" --db "/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgTERMs/alphafold8080ColabResultsAg_matchedbb/af/6yu8_4acaabcc-e98e-4ddc-a8c1-cece355b44c3_2.555964177325445.pdb"  --o "/dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/ColabAf6yu8DD.txt" --n 10000 --abm