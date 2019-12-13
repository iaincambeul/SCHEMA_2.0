import sys
import re
import os
import csv
import pandas as pd
from Bio.PDB import *
from Bio import SeqIO
from Bio import SeqUtils

global letterToNumber_dict
letterToNumber_dict = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5, "G": 6, "H": 7, "I": 8, "J": 9, "K": 10, "L": 11, "M": 12, "N" : 13, "O" : 14, "P" : 15, "Q" : 16, "R": 17, "S": 18, "T": 19, "U": 20, "V": 21, "W": 22, "X": 23, "Y": 24, "Z": 25}

"""Note: Alignment must have the sequence that stands in for parent 1 at the top."""

"""Sometimes, when a chain lacks structural information, the indices of fdResPartnerObjects will not match the index of the actual Fd residues, leading to errors when
I try to pass those indices in FdAtomToResDict.  Need lengthen fdResPartnerObjects"""

"""Note: There are SEQ IO atom objects and then there are the atom object class I have created.  
My object class is really only better in the sense that it has atom indexes."""


#Object Class Atom
class atom:
  def __init__(self, resIndex, chainID, atomOB, atomInd):
    self.resIndex=resIndex #The residue index that the atom is in
    self.chainID=chainID #The chain that the atom is in
    self.atomOB=atomOB #The SEQIO Atom object
    self.atomInd=atomInd #The index of the atom relative to all atoms in the structure

#Object class contactedRes.  Represents atoms that have indices and calculated contacts.  Requires a list of contacting atoms to instantiate.
class contactedRes:
  def __init__(self, resType, index, contacts):
    self.resType=SeqUtils.seq1(resType)
    self.index=index
    self.contacts=contacts #List of indexes for residues that this residue contacts

class parent:
  def __init__(self, name, residues):
    self.name=name
    self.residues=residues  
    rawSequence=[str(x.resType) for x in residues] #Generate sequence based on residue objects
    condSequence="".join(rawSequence) #Condense list of one letter AA to one string
    self.sequence=condSequence #Initialize condensed sequence as parent sequence

  def spitSequence(self):
    residueList=self.residues
    indexList=[x.index for x in residueList]
    residueMax=max(indexList) #Find maximum index of amino acids
  
    sortingList=[[] for x in range(residueMax)]
    for res in residueList:
      sortingList[res.index-1]=res #Dump residue objects into their appropriate place in the list.  Subtract one to get pythonic indexing from human indexing.

    return residueMax, sortingList

class chimera:
  def __init__(self, gapSequence, coords, coordType, p1, p2, schemaP1, schemaP2):
    self.gapSequence=gapSequence #sequence of chimera with gaps from parent's alignment
    self.coords=coords #coords of chimera crossover
    self.coordType=coordType # 1 or 2, which correspond to which parent the two coordinates of the double crossover map onto    
    self.p1=p1 #First parent object
    self.p2=p2 #Second parent object
    self.schemaP1=schemaP1 #SCHEMA score of chimera relative to first parent
    self.schemaP2=schemaP2 #SCHEMA score of chimera relative to second parent


    self.schema_min=min(schemaP1,schemaP2)
    noGapSequence=[x for x in gapSequence if x != "-"] #Remove all gaps from sequence
    self.sequence=noGapSequence

    self.test1=(1 for a,b in zip(gapSequence,p1.alignSeq) if a!=b)
    self.hamP1=sum(1 for a,b in zip(gapSequence,p1.alignSeq) if a!=b) #Add up and sum 1 for each difference between the aligned residues of the chimera sequence and p1 sequence
    self.hamP2=sum(1 for a,b in zip(gapSequence,p2.alignSeq) if a!=b) #Add up and sum 1 for each difference between the aligned residues of the chimera sequence and p2 sequence

def plate_Grabber(file_name):
  plate=pd.read_csv(file_name, index_col=0)
  return plate

def chain_dict_maker(pdb_file, pName):
  parser= PDBParser()
  structure = parser.get_structure("test", pdb_file)
  chains = Selection.unfold_entities(structure, 'C')
  chain_list=[chain.get_full_id()[2] for chain in chains]
  chain_dict = {}

  #Set baseline values for returned variables in case no chain monikers are detected
  Fd_choice=0
  parentName =""
  #Skip over this if no chains are found
  if len(chain_list) != 0:
    print(chain_list)
    Fd_choice = input("Which uppercase letter corresponds to the Fd in the list (A, B, C, etc.) for " + str(pName) +"?")
    try:
      Fd_choice = letterToNumber_dict[Fd_choice]
    except KeyError:
      print("\nTry again and please input an uppercase single letter that matches the ferredoxin.\n")
      sys.exit()
    print(Fd_choice)
  
  return chain_dict, Fd_choice

def pdb_parser(pdb_file):
  parser= PDBParser()
  structure = parser.get_structure("test", pdb_file)
  chains = Selection.unfold_entities(structure, 'C')
  #We will populate the below list with lists of atoms from the respective chains
  chain_atoms = []
  chain_res = []
  for ind, chain in enumerate(chains):
    working_chain = chains[ind]
    working_res = []
    for res in working_chain.get_list():
      if res.get_id()[0][0]!=" ":#This should skip over hetero-atoms and water.  Organic residues have a blank string in this field
        continue
      else:
        working_res.append(res)
    chain_res.append(working_res)
    working_atoms = []
    for res in working_res:
      working_atoms.append(res.get_list())
    working_atoms = [item for sublist in working_atoms for item in sublist]#This is a list expresion that turns working_atoms into a single list from a list of lists
    #Now to grab PDB sequence
    chain_atoms.append(working_atoms)

  return chain_atoms, chain_res

def reIndexDictMaker(parentResidues, alignmentSequence):
  """Need to re-index both residues and contacts."""
  gapInds=[]
  for i, aa in enumerate(alignmentSequence):    
    if aa == "-":
      gapInds.append(i+1)

  maxInd=(parentResidues[-1].get_full_id()[3][1])  
  minInd=(parentResidues[0].get_full_id()[3][1])
  pdbLength=maxInd-minInd
  
  pdb_Inds= [x+minInd for x in range(pdbLength+1)] #Add min pdb ind so we in agreement with pdb values.  Add 1 to range to make inclusive.
  align_Inds=[x+1 for x in range(pdbLength+1)] #Add one so that we are in line with the human numbering scheme of the alignment.  Add 1 to range to make inclusive.
  #These two lists should have the same length.
  #Add in scalars to match them to their original respective numbering schemes.
  #Add one to align inds if there is a gap detected in the sequence. (Boosts everything after it up)
  #They can then be used to create a dictionary that will change sequence indexes into alignment indices

  for gap in gapInds:
    align_Inds=[x+1 if x >= gap else x for x in align_Inds]#Add +1 to each ind value if it comes on or after a gap

  indDict={}
  for i,key in enumerate(pdb_Inds):
    indDict[key]=align_Inds[i]
  return indDict

def atomSorter(chain_atoms):
  #Need to combine all lists of atoms
  allAtoms = [item for sublist in chain_atoms for item in sublist] #Collapses all sublists into a single list.  This puts all atoms from non-Fd chains into a single list
  #Make atom objects
  atomIndexDict={} #Dictionary relating atom index (relative to all atoms in structure) to 
  for i, atomSeqIO in enumerate(allAtoms):
    resIndex=atomSeqIO.get_full_id()[3][1]
    chainID=letterToNumber_dict[atomSeqIO.get_full_id()[2]]  #Get the chain Letter code and translate to a number with the letterToNumber_dict.
    currentAtom=atom(resIndex, chainID, atomSeqIO, i)
    atomIndexDict[i]=currentAtom
  
  atomToResDict={} #Dictionary that relates atom index to residue Seq IO Object
  i=0
  for key in atomIndexDict.keys():
    atomSeqIO=atomIndexDict[key].atomOB
    residueSeqIO=Selection.unfold_entities([atomSeqIO],"R") #Have to make the atom seq IO object a list so it can be subscripted.  Very annoying, I know.
    atomToResDict[key]=residueSeqIO
  return atomIndexDict, atomToResDict

def distCheck(atom1, atom2, cutoff):
  if atom2-atom1 <= cutoff:
    return True
  else:
    return False

def contact_maker_test(Fd_chainID, parentName, res_list, indDict, atomIndexDict, atomToResDict):
  cutoff=4.5

  fdAtomIndexDict={} #Create dictionary for Fd atoms only

  Fd_atoms = []
  for key in atomIndexDict.keys():
    atom=atomIndexDict[key]
    if atom.chainID == Fd_chainID:
      atomOB=atom.atomOB
      Fd_atoms.append(atomOB) 
      
      fdAtomIndexDict[key]=atom #Start dumping into Fd dictionary


  fdAtomPartnerObjects = {} # This is the atom-based dictionary we will put the Fd partner contacts into
  #Making a contact map with a 4.5 angstrom cutoff

  #for i, atom in enumerate(Fd_atoms): #Going through all Fd atoms
  for key in fdAtomIndexDict.keys(): #Going through all Fd atoms
    atom=fdAtomIndexDict[key]
    res_contact_ids= []
    
    atomSeqIO=atom.atomOB
    center = atomSeqIO.get_coord()

    contactingAtoms=[]
    for atomOther in Fd_atoms:
      if distCheck(atomSeqIO, atomOther, cutoff) == True:
        #print("Yeehaw")
        contactingAtoms.append(atomOther)

    #print(len(contactingAtoms))
    residue_list = Selection.unfold_entities(contactingAtoms, 'R')
    #print(residue_list)
    
    #neighbors = neighbor_searcher.search(center, 4.5) # 4.5 angstroms from center of selected atom
    #residue_list = Selection.unfold_entities(neighbors, 'R') #Finds residues of all neighbors (non-redundantly). R for residues
    
    for item in residue_list:
      if item.get_id()[0][0]!=" ":#This should skip over hetero-atoms and water
        continue
      else:
        res_contact_ids.append(item)
    
    fdAtomPartnerObjects[key]=res_contact_ids


  #Pulling out the res index of the last atom in the selected Fd chain
  #firstFdResIndex=chain_atoms[Fd_ind][0].get_full_id()[3][1]
  #lastFdResIndex = chain_atoms[Fd_ind][-1].get_full_id()[3][1]

  fdResPartnerObjects={} # The keys will be the alignment residue indices and the values will be the SEQIO objects that contact them.
  fdResPartnerIndexs={} #The keys will be the alignment residue indices and the values will be the alignment index of the residues that contact them.  
  
  #This nested for loop will add partner residue objects to the appropriate
  #Fd residue indices using the FdAtomToResDict dictionary as a key
  for atomIndex in fdAtomPartnerObjects.keys():
    contactingResidues=fdAtomPartnerObjects[atomIndex]

    myAtomObject=atomIndexDict[atomIndex] #Get my version of the relevant atom's object
    homeResidue=indDict[myAtomObject.resIndex] #Transform the atoms pdb residuce index to alignment index
    resObjects=[]
    contactingResidueIndices=[]
    
    for res in contactingResidues: #Cycle through residues, dumping either their object or index into a box, omitting duplicates
      pdbResIndex=res.get_full_id()[3][1]
      alignResIndex=indDict[pdbResIndex]
      
      if res not in resObjects:
        resObjects.append(res)
      
      if alignResIndex not in contactingResidueIndices:
        contactingResidueIndices.append(alignResIndex)
    
    fdResPartnerObjects[homeResidue]=resObjects
    fdResPartnerIndexs[homeResidue]=contactingResidueIndices

  parentResidues=res_list[Fd_chainID]
  crs=[]

  if len(parentResidues) != len(fdResPartnerIndexs):
    print("Something wrong with align sequence versus pdb sequence.  Fix.")
    sys.exit()
  
  indKeys=list(fdResPartnerIndexs.keys())

  for i, res in enumerate(parentResidues):
    alignmentIndex=indKeys[i]

    contactIndices=fdResPartnerIndexs[alignmentIndex]
    if alignmentIndex in contactIndices: 
      contactIndices=[x for x in contactIndices if x is not alignmentIndex] #Remove the residue's self index from the list of residues it contacts.  Don't care if a residue contacts itself.

    resName=res.get_resname()

    workingResidue=contactedRes(resName, alignmentIndex, contactIndices)# Create contacted Res objects
    crs.append(workingResidue)
  
  p=parent(parentName, crs)
  return p


def contact_maker(Fd_chainID, parentName, res_list, indDict, atomIndexDict, atomToResDict):

  fdAtomIndexDict={} #Create dictionary for Fd atoms only

  Fd_atoms = []
  for key in atomIndexDict.keys():
    atom=atomIndexDict[key]
    if atom.chainID == Fd_chainID:
      atomOB=atom.atomOB
      Fd_atoms.append(atomOB) 
      
      fdAtomIndexDict[key]=atom #Start dumping into Fd dictionary
  

  #Creates NeighborSearch object with Fd_atoms list. Generates the object only.  Using only the Fd atoms here.  
  #Can use nonFdAtoms if wanting to tests contacts with other chains.
  neighbor_searcher = NeighborSearch(Fd_atoms)
  

  fdAtomPartnerObjects = {} # This is the atom-based dictionary we will put the Fd partner contacts into
  #Making a contact map with a 4.5 angstrom cutoff

  #for i, atom in enumerate(Fd_atoms): #Going through all Fd atoms
  for key in fdAtomIndexDict.keys(): #Going through all Fd atoms
    atom=fdAtomIndexDict[key]
    res_contact_ids= []
    
    atomSeqIO=atom.atomOB
    center = atomSeqIO.get_coord()

    neighbors = neighbor_searcher.search(center, 4.5) # 4.5 angstroms from center of selected atom
    residue_list = Selection.unfold_entities(neighbors, 'R') #Finds residues of all neighbors (non-redundantly). R for residues
    for item in residue_list:
      if item.get_id()[0][0]!=" ":#This should skip over hetero-atoms and water
        continue
      else:
        res_contact_ids.append(item)
    
    fdAtomPartnerObjects[key]=res_contact_ids


  #Pulling out the res index of the last atom in the selected Fd chain
  #firstFdResIndex=chain_atoms[Fd_ind][0].get_full_id()[3][1]
  #lastFdResIndex = chain_atoms[Fd_ind][-1].get_full_id()[3][1]

  fdResPartnerObjects={} # The keys will be the alignment residue indices and the values will be the SEQIO objects that contact them.
  fdResPartnerIndexs={} #The keys will be the alignment residue indices and the values will be the alignment index of the residues that contact them.  
  
  #This nested for loop will add partner residue objects to the appropriate
  #Fd residue indices using the FdAtomToResDict dictionary as a key
  for atomIndex in fdAtomPartnerObjects.keys():
    contactingResidues=fdAtomPartnerObjects[atomIndex]

    myAtomObject=atomIndexDict[atomIndex] #Get my version of the relevant atom's object
    homeResidue=indDict[myAtomObject.resIndex] #Transform the atoms pdb residuce index to alignment index
    resObjects=[]
    contactingResidueIndices=[]
    
    for res in contactingResidues: #Cycle through residues, dumping either their object or index into a box, omitting duplicates
      pdbResIndex=res.get_full_id()[3][1]
      alignResIndex=indDict[pdbResIndex]
      
      if res not in resObjects:
        resObjects.append(res)
      
      if alignResIndex not in contactingResidueIndices:
        contactingResidueIndices.append(alignResIndex)
    
    fdResPartnerObjects[homeResidue]=resObjects
    fdResPartnerIndexs[homeResidue]=contactingResidueIndices

  parentResidues=res_list[Fd_chainID]
  crs=[]

  if len(parentResidues) != len(fdResPartnerIndexs):
    print("Something wrong with align sequence versus pdb sequence.  Fix.")
    sys.exit()
  
  indKeys=list(fdResPartnerIndexs.keys())

  for i, res in enumerate(parentResidues):
    alignmentIndex=indKeys[i]

    contactIndices=fdResPartnerIndexs[alignmentIndex]
    if alignmentIndex in contactIndices: 
      contactIndices=[x for x in contactIndices if x is not alignmentIndex] #Remove the residue's self index from the list of residues it contacts.  Don't care if a residue contacts itself.

    resName=res.get_resname()

    workingResidue=contactedRes(resName, alignmentIndex, contactIndices)# Create contacted Res objects
    crs.append(workingResidue)
  
  p=parent(parentName, crs)
  return p


#Takes correctly indexed parents and removes the contacts from all amino acids which will not be changed by recombination.
def redundancyRemover(indexedParents, alignmentSequence1, alignmentSequence2):
  redundantIndices=[] #List that will hold the index of amino acids that are the same in each parent.
  redDict={}
  for i, aa in enumerate(alignmentSequence1):    
    if alignmentSequence1[i] == alignmentSequence2[i]:
      redundantIndices.append(i+1)#Need to add one to go from pythonic indices to pdb indices.
      redDict[i]=aa


  for res in indexedParents[0].residues:#Get rid of contacts in all redundant residues in parent 1
    if res.index in redundantIndices:
      res.contacts=[] 

  for res in indexedParents[1].residues:#Get rid of contacts in all redundant residues in parent 2
    if res.index in redundantIndices:
      res.contacts=[] 

  return indexedParents

#This script will cycle through all possible double crossover recombinants of the two parents
#calculate schema scores and generate chimera objects that store schema, sequence, and hamming distance information. 
def schemaCalc(parent1, parent2):
  residueMax1, sortingListP1=parent1.spitSequence()
  residueMax2, sortingListP2=parent2.spitSequence()
  
  chimerap1List=[] #List to hold all chimeras made by inserting a chunk of P2 into P1
  coordType=1
  parent1SchemaCount=0 #Calculate SCHEMA score relative to parent1, which was inserted with a piece of parent2.
  parent1Contacts=[]
  for res in sortingListP1[:0]+sortingListP1[1:]:
    a=1
    b=2
    tCount=0
    try:
      tCount=sum(1 for con in res.contacts if not a > con and not con >= b) #Add one if the contact is not less than beginning of insert or greater/equal to the index right after insert.
      parent1Contacts+=[con for con in res.contacts if not a > con and not con >= b]      
      for con in res.contacts: 
        if not a > con and not con >= b:
          print(res.resType)
          print(res.index)
          print(con)
    except AttributeError:
      continue        
    parent1SchemaCount+=tCount 
  
  """for x in range(residueMax1):
    for y in range(x+1,residueMax1): #After this, all processes go towards making a single chimera object.
      
      a=x+1 #First index for pdb indexes (always +1 pythonic indexes)
      b=y+1 #Second index for pdb indexes (always +1 pythonic indexes)
      
      parent2SchemaCount=0 #Calculate SCHEMA score relative to parent2, which was inserted into parent 1.
      parent2Contacts=[]
      for res in sortingListP2[x:y]:
        tCount=0
        try: 
          tCount=sum(1 for con in res.contacts if not a <= con < b)
          parent2Contacts+=[con for con in res.contacts if not a <= con < b]
        except AttributeError:
          continue
        parent2SchemaCount+=tCount

      parent1SchemaCount=0 #Calculate SCHEMA score relative to parent1, which was inserted with a piece of parent2.
      parent1Contacts=[]
      for res in sortingListP1[:0]+sortingListP1[1:]:
        tCount=0
        try:
          tCount=sum(1 for con in res.contacts if not a > con and not con >= b) #Add one if the contact is not less than beginning of insert or greater/equal to the index right after insert.
          parent1Contacts+=[con for con in res.contacts if not a > con and not con >= b]      
          for con in res.contacts: 
            if not a > con and not con >= b:
              print(res.resType)
              print(res.index)
              print(con)
        except AttributeError:
          continue        
        parent1SchemaCount+=tCount 
      
      cGapSequence= parent1.alignSeq[:x] + parent2.alignSeq[x:y] + parent1.alignSeq[y:] #Make double crossover chimera sequence by taking bits of both parents
      workingChim=chimera(cGapSequence, [x,y], coordType, parent1, parent2, parent1SchemaCount,parent2SchemaCount)#Make chimera and add to first list
      workingChim.p2Cons=parent2Contacts
      workingChim.p1Cons=parent1Contacts
      chimerap1List.append(workingChim)
  chimerap2List=[] #List to hold all chimeras made by inserting a chunk of P1 into P2
  coordType=2
  for x in range(residueMax1):
    for y in range(x+1,residueMax1): #After this, all processes go towards making a single chimera object.
      
      a=x+1 #First index for pdb indexes (always +1 pythonic indexes)
      b=y+1 #Second index for pdb indexes (always +1 pythonic indexes)
      
      parent1SchemaCount=0 #Calculate SCHEMA score relative to parent1, which was inserted into parent 2.
      parent1Contacts=[]
      for res in sortingListP1[x:y]:
        tCount=0
        try:
          tCount=sum(1 for con in res.contacts if not a <= con < b)
          parent1Contacts+=[con for con in res.contacts if not a <= con < b]
        except AttributeError:
          continue         
        parent1SchemaCount+=tCount

      parent2SchemaCount=0 #Calculate SCHEMA score relative to parent2, which was inserted with a piece of parent1.
      parent2Contacts=[]
      for res in sortingListP2[:x]+sortingListP2[y:]:
        tCount=0
        try:
          tCount=sum(1 for con in res.contacts if not a > con and not con >= b)
          parent2Contacts+=[con for con in res.contacts if not a > con and not con >= b]

        except AttributeError:
          continue        
        parent2SchemaCount+=tCount 

      cGapSequence= parent1.alignSeq[:x] + parent2.alignSeq[x:y] + parent1.alignSeq[y:] #Make double crossover chimera sequence by taking bits of both parents
      
      workingChim=chimera(cGapSequence, [x,y], coordType, parent1, parent2,parent1SchemaCount,parent2SchemaCount)#Make chimera and add to first list
      workingChim.p2Cons=parent2Contacts
      workingChim.p1Cons=parent1Contacts
      chimerap2List.append(workingChim)
  return chimerap1List,chimerap2List"""

def file_writer(chim1List, chim2List):
  final_name = "allChimeras.csv"
  with open(final_name,'w', newline="") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["Sequence", "Coords", "Coord Type", "Parent 1 E", "Parent 2 E", "Min. E", "Parent 1 Hamming", "Parent 2 Hamming", "Parent 1 Contacts", "Parent 2 Contacts"])
    for chim in chim1List:
      writer.writerow([chim.sequence, chim.coords, chim.coordType, chim.schemaP1,chim.schemaP2, chim.schema_min, chim.hamP1, chim.hamP2, chim.p1Cons, chim.p2Cons])
    for chim in chim2List:
      writer.writerow([chim.sequence, chim.coords, chim.coordType, chim.schemaP1, chim.schemaP2, chim.schema_min, chim.hamP1, chim.hamP2, chim.p1Cons, chim.p2Cons])
 
def main():
  pdbFiles=["1rfk.pdb","newPssm2.pdb"]
  #pdbFiles=["pssm.pdb","1rfk.pdb"]
  pNames=["1rfk","pssm"]
  alignmentFile="parentAlignEd.txt"
  chimFile="chimIndices.csv"
  chims=plate_Grabber(chimFile)
  
  choices=[1,0]

  parents=[]
  
  #Getting the sequences of the aligned parents
  align=SeqIO.parse(alignmentFile, "fasta")
  alignSequences=[]
  for record in align:
    alignSequences.append(record.seq)
  print(alignSequences)

  indDicts=[]
  #Loop through several functions to add a parent object to the list.
  for i, file in enumerate(pdbFiles):
    pName=pdbFiles[i].replace(".csv","")#Get structure name
    
    if "choices" not in locals(): #Use biopython to query user about chain of interest if not already assigned.
      c_List, Fd = chain_dict_maker(file, pName)
    else:
      Fd=choices[i]
    
    c_atoms, c_res = pdb_parser(file) #Pulls out all atoms and residues from pdb file, each into single lists
    Fd_res = c_res[Fd] #Selects Fd chain from all residue chains
    indDict=reIndexDictMaker(Fd_res,alignSequences[i]) #Create dictionary relating non-aligned residue to aligned residues.
    indDicts.append(indDict)

    atomIndexDict, FdAtomToResDict=atomSorter(c_atoms)

    parentOb = contact_maker_test(Fd, pName, c_res, indDict, atomIndexDict, FdAtomToResDict) #Finds all residues that contact selected chain residues
    parents.append(parentOb)
    with open(str(pNames[i])+".csv", "w", newline="") as csvFile:
      csvWriter= csv.writer(csvFile)
      csvWriter.writerow(["Res Number", "Contacts", "Number of Contacts"])
      resList=parentOb.residues
      for res in resList:
        csvWriter.writerow([res.index, res.contacts, len(res.contacts)])


    #completeFdRes = chain_filler(Fd_res) #Fills in blank lists where non-structured Fd should be in the chain


  p1AlignSeq=list(str(alignSequences[0]))
  p2AlignSeq=list(str(alignSequences[1]))
  #print(p1AlignSeq)
  #print(p2AlignSeq)
  nonRedParents= redundancyRemover(parents,p1AlignSeq,p2AlignSeq)

  """p1.alignSeq=p1AlignSeq
  p2.alignSeq=p2AlignSeq
  schemaCalc(nonRedParents[0], nonRedParents[1])
  #p1Chims, p2Chims=schemaCalc(nonRedParents[0], nonRedParents[1])"""

  """file_writer(p1Chims,p2Chims)"""
  

if __name__ == '__main__':
    main() 