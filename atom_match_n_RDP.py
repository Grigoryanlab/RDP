# Matches the atoms of backbones between two structures, so they can have RMSDs done to them - assumes that if something is on one complete chain in one structure, it's on one complete chain in another. Does not allow mismatches; each chain must be at least 80% a match to the other (so up to 20% gaps)

# For underlying distance metrics (RMSDs and DOCKQ) to be computed, atoms need to be paired between all structures, wherein each pair is the correct position of the atom vs its predicted position. For this our script can match either all backbone atoms or all heavy atoms between a true structure \& its predicted model. Only the chains for each binding partner in the crystal structure need to be specified, so that many different models with different chain naming standards can all be easily backbone-matched in batch and then scored in batch, which is useful for evaluating models produced by different methods.

import argparse
import sys
from Bio import pairwise2
from Bio.Seq import Seq
import itertools
from operator import itemgetter
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-t", dest="t", required=True, help="the path to the true complex structure")
parser.add_argument("-m", dest="m", required=True, help="the path to the model complex structure")
parser.add_argument("-ta", dest="ta", required=False, help="the chains making up the first docking partner out of two; if there are multiple binding locations for one of the partners, this must be that partner. If antibody-binding-mode RDP is being run, this must be the antibody. Enter chains separated by commas, with alternate binding locations specified by colons (e.g a two-chain partner with multiple binding locations, either chains A&B or C&D, would be entered like 'A:C,B:D'. If there are only two chains in the crystal structure, -ta and -tb are not necessary; if there are multiple chains you can only provide one and all other chains will be presumed to go to the other flag.")
parser.add_argument("-tb", dest="tb", required=False, help="the chains making up the second docking partner out of two; if there are multiple binding locations for one of the partners, this cannot be that partner. If antibody-binding-mode RDP is being run, this must be the antigen. Enter chains split by ',' (e.g a two-chain partner would be entered like 'E,F'. If there are only two chains in the crystal structure, -ta and -tb are not necessary; if there are multiple chains you can only provide one and all other chains will be presumed to go to the other flag.")
parser.add_argument("-d", dest="d", required=True, help="the underlying distance metric used for the distribution the P-value is taken from; accepts 'c' for complex RMSD, 'r' for receptor + ligand RMSD, and 'dq' for DOCKQ. R+L RMSD is the default.")
parser.add_argument("-o", dest="o", required=True, help="the out path for the text file containing the RDP results'")
parser.add_argument("-ml", dest="ml", required=False, help="the path to MST's lib subdir, like '/path/to/MST/lib'; not required if the library path is where it's supposed to be, as the makefile will save the path and this script will automatically load it.")
parser.add_argument("-e", dest="explicit", required=False, help="explicitly designate which chains to match between the true structure and the model, instead of having this be automatically done. Enter without spaces, chains to pair separated by colons (true first, mdodel second) and sets of chains separated by commas, like 'A:H,B:L,E:E,F:F'. If multiple binding locations are allowed for the true structure, use the first one entered (e.g. for -ta 'A:C,B:D' you would use A and B for the true structures' chains)")
parser.add_argument("-mod", dest="mod", required=False, help="model binding partners for random docking instead of crystal binding partners")
parser.add_argument("-test", dest="test", required=False, help="test mode - enter the pathway to start off all test output files for this flag")
parser.add_argument("-mbl", dest="mbl", required=False, help="there are multiple binding locations for a; this flag will take the path / paths to pdb files containing backbone-matched alternative binding locations, separated out by commas if multiple")
parser.add_argument("-mult", dest="mult", required=False, help="allow matching of multiple chains in the true structure, to the same chain in the model (useful for cases like SnugDock or AbAdapt models, wherein all antigen chains are labeled as just chain A)")
parser.add_argument("-n", dest="n", required=False, help="the number of valid docking positions to generate; defaults to 1 million, which should take around 1 hour per every 500 residues in the complex")
parser.add_argument("-i", dest="i", required=False, help="the number of valid interaction residues required to accept the structure; defaults to 3")
parser.add_argument("-cla", dest="cla", required=False, help="the number of clashes to allow; defaults to 0")
parser.add_argument("-sd", dest="sd", required=False, help="the standard deviation of the random step size; defaults to 2.0 Angstroms")
parser.add_argument("-abm", dest="abm", required=False, help="an optional antibody / TCR mode, which will treat all the loops of the antibody as the binding residues; the antibody structure must use IMGT numbering and its structure must be entered using the -ca argument not the -cb argument.")
parser.add_argument("-al", dest="al", required=False, help="an optional list of binding residues for the first docking partner (A); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like 'A,100, ;A,100,A;A,100,B'. If you're only giving binding residues for one of the two partners, it must be this one (i.e. you cannot give bl without giving al)")
parser.add_argument("-bl", dest="bl", required=False, help="an optional list of binding residues for the second docking partner (B); this will bias the random docking towards conformations including those residues in the binding site. Should be a list of tuples, where each tuple has the chain followed by the residue number followed by the residue insertion code (or ' ' if no insertion code). Separate each member of the tuple with a comma, and each tuple with a semi-colon, like 'A,100, ;A,100,A;A,200, '")
parser.add_argument("-q", dest="q", required=False, help="an optional quick mode for the --al or --abm flags, wherein the docking distribution is skewed towards conformations involving the binding residues given on the A side, to make the calculation much faster")
parser.add_argument("-limA", dest="limA", required=False, help="Must be used with al or abm, and bl. Optionally limit the range of angles of accepted dockings, wherein each binding partner has a line drawn between its geometric center and the geometric center of its binding residues; one line is used as an axis, and the other is used to calculate the angle of tilt relative to that axis. Primarily useful for TCRpMHC random dockings.")
parser.add_argument("-cache", dest="cache", required=False, help="Use a cached distribution of underlying distance metrics instead; supply the path to that csv file here.")
parser.add_argument("-j", dest="j", required=False, help="Just print the underlying distance metric between the two structures instead of doing comparisons.")
arguments = parser.parse_args()

if arguments.q and ((not arguments.al) and (not arguments.abm)):
    print("-q can only be used with -al or -abm")
    quit()

if arguments.bl and ((not arguments.al) and (not arguments.abm)):
    print("-bl can only be used with -al or -abm")
    quit()

if arguments.limA and ((not arguments.al) and (not arguments.abm)):
    print("-limA can only be used with both (-al or -abm) and -bl")
    quit()

if arguments.limA and (not arguments.bl):
    print("-limA can only be used with both -al and -bl")
    quit()

if (arguments.d != 'dq') and (arguments.d != 'r') and (arguments.d != 'c'):
    print("-m must be 'dq' for DOCKQ, 'c' for complex RMSD, or 'r' for receptor + ligand RMSD")
    quit()

if not arguments.ml:
    f = open("mstObjsDir.txt","r")
    mlBase = f.readline()
    arguments.ml = mlBase.replace("objs","lib").strip()

sys.path.append(arguments.ml)
import mstpython as mst

fullAtoms = False 
if arguments.d == 'dq':
    fullAtoms = True

outPathP1 = "matchedBackbones/P1"
outPathP2 = "matchedBackbones/P2"

# make a dictionary for paired chains, if they're explicity specified

if arguments.explicit:
    chainMappingDict = {}
    chainMappingSplit = arguments.explicit.split(",")
    for i in chainMappingSplit:
        chainMappingSet = i.split(":")
        if (len(chainMappingSet) != 2) and (not arguments.mult):
            print("chains must be mapped one to one")
            quit()
        chainMappingDict[chainMappingSet[0]] = chainMappingSet[1]

# function for getting a dictionary of chains & their sequences, from a structure and the list of chain IDs

def resDictMaker(struct, chainList):
    resesDict = {}
    for chID in chainList:
        ch = struct.getChainByID(chID)
        reses = ch.getResidues()
        chSeq = ''
        chNums = []
        for i in range(len(reses)):
            res = reses[i]
            rName = res.name
            rNum = res.num
            rCode = res.iCode
            rIndex = res.getResidueIndex()
            try:
                rSname = mst.SeqTools.tripleToSingle(rName,"+")
            except:
                rSname = 'X'
            if rSname == '?':
                rSname = 'X'
            chSeq += rSname
            chNums.append(rNum)
        resesDict[chID] = [chSeq,chNums] #,seq]
    return(resesDict)

tPath = arguments.t

try:
    structA = mst.Structure(tPath, "QUIET")
except:
    print("loading structure for " + tPath + " failed for some reason; quitting")
    exit()

aChains = []
aChains1 = []
aChains2 = []
i = 0
while i < structA.chainSize():
    aChain = structA.getChain(i)
    aChainID = aChain.id
    aChains.append(aChainID)
    i += 1

if arguments.ta and not arguments.tb:
    aChains1 = arguments.ta.split(",")
    for aC in aChains:
        if aC in aChains1:
            continue 
        aChains2.append(aC)
elif arguments.tb and not arguments.ta:
    aChains2 = arguments.tb.split(",")
    for aC in aChains:
        if aC in aChains2:
            continue 
        aChains1.append(aC)
elif arguments.ta and arguments.tb:
    aChains1 = arguments.ta.split(",")
    aChains2 = arguments.tb.split(",")
else:
    if len(aChains) == 2:
        aChains1.append(aChains[0])
        aChains2.append(aChains[1])
    else:
        print("either -ta or -tb is required if there are more than two chains in the crystal structure")
        exit()

aaCdict = {} #a alt chain Dict
aChainsBaseChains = []
for aC in aChains1:
    if ":" in aC:
        aCplit = aC.split(":")
        aaC1 = aCplit[0]
        aaC2list = []
        for i in range(1,len(aCplit)):
            aaC2 = aC.split(":")[1]
            aaC2list.append(aaC2)
        aaCdict[aaC1] = aaC2list
        aChainsBaseChains.append(aaC1)
    else:
        aChainsBaseChains.append(aC)

aChains1base = aChainsBaseChains.copy()

for aC in aChains2:
    aChainsBaseChains.append(aC)

mPath = arguments.m
mStart = mPath.split("/")[-1].split(".pdb")[0]

try:
    structB = mst.Structure(mPath, "QUIET")
except:
    print("loading structure for " + mPath + " failed for some reason; quitting")
    exit()

aResesDict = resDictMaker(structA,aChainsBaseChains)

bSize = structB.chainSize()
bChains = []
bSeqs = []
i = 0
while i < bSize:
    bChain = structB.getChain(i)
    bChainID = bChain.id
    bChains.append(bChainID)
    i += 1

bResesDict = resDictMaker(structB,bChains)

a2bAlignments = {}

for keyA in aResesDict:

    a2bAlignments[keyA] = []

    for keyB in bResesDict:

        if arguments.explicit:
            if keyA in chainMappingDict:
                if chainMappingDict[keyA] != keyB:
                    continue
            else:
                print("chains in mapping not found in structures; quitting")
                quit()

        seqA = aResesDict[keyA][0].strip("X")
        numsA = aResesDict[keyA][1]
        seqB = bResesDict[keyB][0].strip("X")
        numsB = bResesDict[keyB][1]
        aLen = float(len(seqA))
        bLen = float(len(seqB))
        if (aLen*10 < bLen) or (bLen*10 < aLen): # takes care of issues where very small seqs are paired against very large ones, which takes forever with pairwise2
            continue

        #aXcount = 0
        #for char in seqA:
        #    if char == 'X':
        #        aXcount += 1
        #aXper = aXcount/aLen

        #bXcount = 0
        #for char in seqB:
        #    if char == 'X':
        #        bXcount += 1
        #bXper = bXcount/bLen

        highestLen = max(len(seqA),len(seqB))
        #mismatchPenalty = -10 * highestLen

        alignmentList = pairwise2.align.localms(seqA, seqB, 1, -100000, -0.1, -0.01) # identical chars get 1 point, mismatches get -100000, gaps get -0.1, and gap extensions -0.01

        for a in alignmentList:
            abScore = a.score
            scorePercent1 = (float(abScore)/len(seqA))
            scorePercent2 = (float(abScore)/len(seqB))

            if (scorePercent1 >= 0.8) or (scorePercent2 >= 0.8):

                seqAalign = a.seqA
                seqBalign = a.seqB

                aChain = structA.getChainByID(keyA)
                aReses = aChain.getResidues()
                bChain = structB.getChainByID(keyB)
                bReses = bChain.getResidues()
                aAtoms = []
                bAtoms = []
                altaAtoms = []

                altResesDict = {}
                if keyA in aaCdict:
                    altChainsList = aaCdict[keyA]
                    for altChainID in altChainsList:
                        altChain = structA.getChainByID(altChainID)
                        altReses = altChain.getResidues()
                        altResesDict[altChainID] = altReses
                        altaAtoms.append([])

                # parse down atoms to only those shared between structures, and add to a list of

                i = 0
                aIndex = 0
                bIndex = 0

                while (i < len(seqAalign)) and (i < len(seqBalign)):
                    resA = seqAalign[i]
                    resB = seqBalign[i]
                    if ((resA != '-') and (resB == '-')):
                        aIndex += 1
                    elif ((resA == '-') and (resB != '-')):
                        bIndex += 1
                    elif ((resA == '-') and (resB == '-')):
                        pass
                    elif (resA == resB):
                        # get C residue
                        currAres = aReses[aIndex]
                        currAnum = numsA[aIndex]
                        # get D residue
                        currBres = bReses[bIndex]
                        currBnum = numsB[bIndex]

                        if fullAtoms:
                            # append /all/ those atoms shared by A & B to atom lists for A & B

                            aAtomlist = currAres.getAtoms()
                            bAtomlist = currBres.getAtoms()


                            #get any alt C residue bbs *** works, prints a list containing AtomLists (just one with this test)
                            altAtomsList = []
                            for key in altResesDict:
                                altRes = altResesDict[key][aIndex]
                                altAtomlist = altRes.getAtoms()
                                altAtomsList.append(altAtomlist)
                                #print("testing altAtomlist:")
                                #print(altAtomlist)

                            #print('testing altAtomsList:')
                            #print(altAtomsList)

                            ii = 0
                            while (ii < len(aAtomlist)):
                                aAtom = aAtomlist[ii]
                                aAtomName = aAtom.name

                                foundBatom = None
                                bb = 0
                                while (bb < len(bAtomlist)):
                                    possBatom = bAtomlist[bb]
                                    possBatomName = possBatom.name
                                    if aAtomName == possBatomName:
                                        foundBatom = possBatom
                                        break
                                    bb += 1
                                if foundBatom == None:
                                    ii += 1
                                    continue

                                possibleAltAatoms = []
                                if len(altAtomsList) != 0:
                                    aa = 0
                                    while aa < len(altAtomsList): # a list of atom lists
                                        checkingList = altAtomsList[aa] # go through each list to find the atom
                                        cl = 0
                                        foundAltAtom = None
                                        while cl < len(checkingList):
                                            checkingAltAtom = checkingList[cl]
                                            checkingAltAtomName = checkingAltAtom.name
                                            if checkingAltAtomName == aAtomName:
                                                foundAltAtom = checkingAltAtom
                                                break
                                            cl += 1
                                        if foundAltAtom == None: # if not found in even one list, break, so you can iterate then continue as below in a few lines
                                            break
                                        else:
                                            possibleAltAatoms.append(foundAltAtom)
                                        aa += 1

                                    if foundAltAtom == None: # if not found in even one list, iterate then continue
                                        ii += 1
                                        continue

                                aAtoms.append(aAtom)
                                bAtoms.append(foundBatom)
                                aal = 0
                                #*** works - right number of atoms in altaAtoms[0]!!
                                while aal < len(possibleAltAatoms):
                                    altaAtoms[aal].append(possibleAltAatoms[aal])
                                    aal += 1
                                ii += 1

                        else:
                            # append those atoms shared by A & B to atom lists for A & B

                            aBBlist = mst.RotamerLibrary.getBackbone(currAres,False)
                            bBBlist = mst.RotamerLibrary.getBackbone(currBres,False)

                            #get any alt C residue bbs *** works, prints a list containing AtomLists (just one with this test)
                            altBBsList = []
                            for key in altResesDict:
                                altRes = altResesDict[key][aIndex]
                                altBBlist = mst.RotamerLibrary.getBackbone(altRes,False)
                                altBBsList.append(altBBlist)

                            ii = 0
                            while ii < 4:
                                aBBatom = aBBlist[ii]
                                bBBatom = bBBlist[ii]

                                # get alt atoms *** works, prints a list containing a single atom when there's only one alt position given
                                possibleAltAatoms = []
                                for altBB in altBBsList:
                                    altAtom = altBB[ii]
                                    possibleAltAatoms.append(altAtom)

                                #print('testing possibleAltAatoms')
                                #print(possibleAltAatoms)

                                if (aBBatom != None) and (bBBatom != None) and (None not in possibleAltAatoms):
                                    aAtoms.append(aBBatom)
                                    bAtoms.append(bBBatom)
                                    aal = 0
                                    #*** works - right number of atoms in altaAtoms[0]!!
                                    while aal < len(possibleAltAatoms):
                                        altaAtoms[aal].append(possibleAltAatoms[aal])
                                        aal += 1
                                ii += 1

                        aIndex += 1
                        bIndex += 1
                    else: # a rare mismatch case that shouldn't show up in practice
                        aIndex += 1
                        bIndex += 1
                    i += 1

                # add atom list to the set of alignments

                a2bAlignments[keyA].append([keyB,aAtoms,bAtoms,altaAtoms])

            else: #don't go over remaining scores below 0.8
                break

# for all keyAs that only have one match, remove their keyBs from consideration in other keyAs
if not arguments.mult:
    for key in a2bAlignments:
        aSet = a2bAlignments[key]
        if len(aSet) == 1:
            keyB = aSet[0][0]
            for key2 in a2bAlignments:
                if key2 == key: # if the same position, skip - otherwise, remove keyB from consideration from other sets
                    continue
                aSet2 = a2bAlignments[key2]
                i = 0
                while i < len(aSet2):
                    aSetList = aSet2[i]
                    if aSetList[0] == keyB:
                        aSet2.pop(i)
                    else:
                        i += 1

# if a chain requested has no match, error and print which chain had no match:
for key in a2bAlignments:
    aSet = a2bAlignments[key]
    if len(aSet) == 0:
        print("no match found for chain: " + key)
        quit()

aChainList = []
for key in a2bAlignments:
    aChainList.append([]) # make a list to hold possible alignments for the current A chain
    aSet = a2bAlignments[key]
    for aAlign in aSet:
        # make a list for a specific alignment formatted like: A chain, B chain, A atoms, B atoms, A alt atoms list; so if there are multiple possible alignments, you get mulitple lists starting w/ A, for example
        aAlignList = [key]
        aAlignList.extend(aAlign)
        aChainList[-1].append(aAlignList)

allCombos = list(itertools.product(*aChainList))

if not arguments.mult:
    i = 0
    while i < len(allCombos):
        currCombo = allCombos[i] #side note: currCombo is a tuple, wherein each entry is a list, containing A chain, B chain, A atoms list, B atoms list, and alt A atoms list of lists
        ii = 0
        bList2 = []
        broke = False
        while ii < len(currCombo):
            currB = currCombo[ii][1] # get B chain for this A chain's match-up
            if currB in bList2:
                allCombos.pop(i)
                broke = True
                break
            bList2.append(currB)
            ii += 1
        if not broke:
            i += 1

rc = mst.RMSDCalculator()

comboRMSDlist = []
for combo in allCombos:
    aAtomsTogether = []
    bAtomsTogether = []
    aAltAtomsTogether = []
    ci = 0
    while ci < len(combo[0][4]):
        aAltAtomsTogether.append([])
        ci += 1
    for chainAlignment in combo:
        aAtomsTogether.extend(chainAlignment[2])
        bAtomsTogether.extend(chainAlignment[3])
        cii = 0
        while cii < len(chainAlignment[4]):
            #print("testing chainAlignment[4]") #looks good - is a list of atom lists, just one cuz just one alt location w/ this test
            aAltAtomsTogether[cii].extend(chainAlignment[4][cii])
            cii += 1
    aAPV = mst.AtomPointerVector(aAtomsTogether)
    bAPV = mst.AtomPointerVector(bAtomsTogether)
    currRMSD = rc.apvBestRMSD(aAPV,bAPV)
    comboRMSDlist.append([currRMSD, aAtomsTogether, bAtomsTogether,aAltAtomsTogether]) # and now aAltAtomsTogether is a list of atom lists!

# chose the alignment w/ the best RMSD

sortedComboRMSD = sorted(comboRMSDlist, key=itemgetter(0))
bestCombo = sortedComboRMSD[0]
bestRMSD = bestCombo[0]
bestAatoms = bestCombo[1]
bestBatoms = bestCombo[2]
bestAltAatoms = bestCombo[3]


fileEnd = "_" + str(bestRMSD) + ".pdb"
#uniqueID = uuid.uuid4()
partner1Apath = outPathP1 + "/crystal" #A = crystal, B = model
partner2Apath = outPathP2 + "/crystal"
partner1Bpath = outPathP1 + "/model"
partner2Bpath = outPathP2 + "/model"
if fullAtoms:
    partner1Apath = outPathP1 + "/crystalFull" #A = crystal, B = model
    partner2Apath = outPathP2 + "/crystalFull"
    partner1Bpath = outPathP1 + "/modelFull"
    partner2Bpath = outPathP2 + "/modelFull"
if not os.path.exists(partner1Apath):
    os.mkdir(partner1Apath)
if not os.path.exists(partner2Apath):
    os.mkdir(partner2Apath)
if not os.path.exists(partner1Bpath):
    os.mkdir(partner1Bpath)
if not os.path.exists(partner2Bpath):
    os.mkdir(partner2Bpath)

fileCrP1 = partner1Apath + "/" + mStart + fileEnd 
fileModP1 = partner1Bpath + "/" + mStart + fileEnd 
fileCrP2 = partner2Apath + "/" + mStart + fileEnd
fileModP2 = partner2Bpath + "/" + mStart + fileEnd
structP1_crystal = mst.emptyStructure()
structP2_crystal = mst.emptyStructure()
structP1_model = mst.emptyStructure()
structP2_model = mst.emptyStructure()

aa = 0
while aa < len(bestAatoms):
    aaRes = bestAatoms[aa].getAtomParent()
    aaID = aaRes.getChainID(True)
    if aaID in aChains1base: # partner A
        #print("testing P1 found")
        structP1_crystal.addAtom(bestAatoms[aa]) #A = crystal, B = model
        structP1_model.addAtom(bestBatoms[aa])
    elif aaID in aChains2:
        #print('testing p2 found')
        structP2_crystal.addAtom(bestAatoms[aa]) #A = crystal, B = model
        structP2_model.addAtom(bestBatoms[aa])
    aa += 1

structP1_crystal.writePDB(fileCrP1,"QUIET")
structP1_model.writePDB(fileModP1,"QUIET")
structP2_crystal.writePDB(fileCrP2,"QUIET")
structP2_model.writePDB(fileModP2,"QUIET")
print(mStart + " complex RMSD was: " + str(bestRMSD))

aaIndex = 0
for aa in bestAltAatoms:
    structAalt_out = mst.emptyStructure()
    structAalt_out.addAtoms(aa)
    altPath = outPathP1 + "/crystalAlt"
    if fullAtoms:
        altPath = outPathP1 + "/crystalFullAlt"
    if not os.path.exists(altPath):
        os.mkdir(altPath)
    fileAalt = altPath + "/" + mStart + "_" + str(aaIndex) + fileEnd
    structAalt_out.writePDB(fileAalt,"QUIET")
    aaIndex += 1

jobCommand = 'bin/dockingDistribution --ta "' + fileCrP1 + '" --tb "' + fileCrP2 + '" --ma "' +  fileModP1 + '" --mb "' + fileModP2 + '"  --o "' + arguments.o + '"'

# add on options to jobCommand 
if arguments.mod:
    jobCommand += ' --m '

if arguments.d == 'r':
    jobCommand += ' --r 1'

if arguments.d == 'dq':
    jobCommand += ' --dq 1'

if arguments.test:
    jobCommand += ' --t ' + arguments.test

if arguments.mbl:
    jobCommand += ' --mbl '

if arguments.n:
    jobCommand += ' --n ' + arguments.n

if arguments.i:
    jobCommand += ' --i ' + arguments.i

if arguments.cla:
    jobCommand += ' --cla ' + arguments.cla

if arguments.sd:
    jobCommand += ' --sd ' + arguments.sd

if arguments.abm:
    jobCommand += ' --abm 1'

if arguments.al:
    jobCommand += " --al '" + arguments.al + "'" 

if arguments.bl:
    jobCommand += " --bl '" + arguments.bl + "'" 

if arguments.q:
    jobCommand += ' --q '

if arguments.limA:
    jobCommand += ' --limA ' + arguments.limA

if arguments.cache:
    jobCommand += ' --cache ' + arguments.cache + ' '

if arguments.j:
    jobCommand += ' --j '

print(jobCommand)

proc = subprocess.Popen(jobCommand, stdout=subprocess.PIPE, shell=True)
(bOut, bErr) = proc.communicate()
jobOut = bOut.decode("utf-8")
if bErr:
    jobErr = bErr.decode("utf-8")
    print("Error running dockingDistribution binary: " + jobErr)
else:
    print(jobOut)
