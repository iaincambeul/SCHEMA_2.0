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

def chain_dict_maker(pdb_file):
  pdb_open = open(pdb_file, "r")
  chain_list = []
  chain_list = re.findall("DBREF\s+\w+\s+(\w)\s+\w+\s+\w+\s+\w+\s+\w+\s+(\w+)", pdb_open.read())
  chain_dict = {}
  for chain in chain_list:
    chain_dict[chain[0]] = chain[1]
  pdb_open.close()
  
  #Set baseline values for returned variables in case no chain monikers are detected
  Fd_choice=0
  parentName =""
  
  #Skip over this if no chains are found
  if len(chain_list) != 0:
    print(chain_dict)
    Fd_choice = input("Which uppercase letter corresponds to the Fd in the list (A, B, C, etc.)? ")
    parentName=chain_dict[Fd_choice]
    try:
      Fd_choice = letterToNumber_dict[Fd_choice]
    except KeyError:
      print("\nTry again and please input an uppercase single letter that matches the ferredoxin.\n")
      sys.exit()
    print(Fd_choice)
  
  return chain_dict, Fd_choice, parentName

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

  j=1#Gap counter. This will used to offset the incremental 
  # discrepancies in the sortingList index caused by the gaps.
  gapInds=[]
  for i, aa in enumerate(alignmentSequence):    
    if aa == "-":         
      gapInds.append(i)

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
    chainID=atomSeqIO.get_full_id()[2]
    currentAtom=atom(resIndex, chainID, atomSeqIO, i)
    atomIndexDict[i]=currentAtom
  
  atomToResDict={} #Dictionary that relates atom index to residue Seq IO Object
  i=0
  for key in atomIndexDict.keys():
    atomSeqIO=atomIndexDict[key].atomOB
    residueSeqIO=Selection.unfold_entities([atomSeqIO],"R") #Have to make the atom seq IO object a list so it can be subscripted.  Very annoying, I know.
    atomToResDict[key]=residueSeqIO
  return atomIndexDict, atomToResDict

def contact_maker(Fd_ind, parentName, res_list, indDict, atomIndexDict, atomToResDict):

  #Need to combine all lists of not selected Fd atoms
  nonFdAtoms = []
  for ind, chain in enumerate(chain_atoms):
    if ind == Fd_ind:
      Fd_atoms = chain_atoms[ind] #Selecting Fd atoms from chain atoms list
    else:
      nonFdAtoms.append(chain_atoms[ind]) #Add non-Fd atoms to non-Fd atoms list

  nonFdAtoms = [item for sublist in nonFdAtoms for item in sublist] #Collapses all sublists into a single list.  This puts all atoms from non-Fd chains into a single list


  #Creates NeighborSearch object with Fd_atoms list. Generates the object only.  Using only the Fd atoms here.  
  #Can use nonFdAtoms if wanting to tests contacts with other chains.
  neighbor_searcher = NeighborSearch(Fd_atoms)
  
  #Creating dictionary so that I can sort contacting partner residues 
  #once I have them in them in the atom indexed list
  atomToResDict = {}

  allAtoms= [item for sublist in chain_atoms for item in sublist]
  #Make a dictionary relating the index of all atoms to their residue number
  for ind, atom in enumerate(allAtoms):
    atomToResDict[ind] = atom.get_full_id()[3][1]-1 #Have to subtract one to get into Pythonic indices
  print(atomToResDict)
  sys.exit()
  
  fdAtomPartnerObjects = {} # This is the atom-based dictionary we will put the Fd partner contacts into
  #Making a contact map with a 4.5 angstrom cutoff
  for i, atom in enumerate(Fd_atoms): #Going through all Fd atoms
    res_contacts = []
    res_contact_ids= []
    
    center = Fd_atoms[i].get_coord()
    
    neighbors = neighbor_searcher.search(center, 4.5) # 4.5 angstroms from center of selected atom
    residue_list = Selection.unfold_entities(neighbors, 'R') #Finds residues of all neighbors (non-redundantly). R for residues
    for item in residue_list:
      if item.get_id()[0][0]!=" ":#This should skip over hetero-atoms and water
        continue
      else:
        res_contacts.append(item.id[1])
        res_contact_ids.append(item)
    
    fdAtomPartnerObjects.append(res_contact_ids)


  #Pulling out the res index of the last atom in the selected Fd chain
  firstFdResIndex=chain_atoms[Fd_ind][0].get_full_id()[3][1]
  lastFdResIndex = chain_atoms[Fd_ind][-1].get_full_id()[3][1]
  fdResPartnerObjects={} # The keys will be the pdb residue indices and the values with the SEQIO objects that contact them.

  #fdResPartnerObjects = [[] for a in range(0, lastFdResIndex+1)]# This is the residue-based list we will put the Fd partner contacts into
  
  #This nested for loop will add partner residue objects to the appropriate
  #Fd residue indices using the FdAtomToResDict dictionary as a key
  #The result will be a list of lists, in which each sublist contains the SEQIO residue objects that contact 
  #the atom with the index in the protein equal to the index in the list of lists.
  for i, atom in enumerate(fdAtomPartnerObjects):
    print(fdAtomPartnerObjects)
    sys.exit()
    for contactingRes in atom:
      print(atom)
      sys.exit()
      FdAtomToResDict[index]
      if contactingRes not in fdResPartnerObjects[FdAtomToResDict[index]]:
        
        fdResPartnerObjects[0] = contactingRes

  print(fdResPartnerObjects)

  parentResidues=res_list[Fd_ind]
  crs=[]
  residueIndices=range(firstFdResIndex, lastFdResIndex+1)

  for res_ind in residueIndices:
    res= parentResidues[res_ind]

    contactIndices=[x.get_full_id()[3][1]-1 for x in fdResPartnerObjects[res_ind]] #Pull out only the indices of the contacting residues for the residue in question (from the box of contacting residues).

    if res_ind in contactIndices: 
      contactIndices=[x for x in contactIndices if x is not res_ind] #Remove the residue's self index from the list of residues it contacts.  Don't care if a residue contacts itself.

    resName=res.get_resname()
    resIndex=res.get_full_id()[3][1]
    resIndex=indDict[resIndex] #Re-index residues to correspond to alignment
    contactIndices=[indDict[x] for x in contactIndices] #Re-index contacting residue to correspond to alignment
 
    workingResidue=contactedRes(resName, resIndex, contactIndices)# Create contacted Res objects
    crs.append(workingResidue)
  
  p=parent(parentName, crs)
  print(p.sequence)
  sys.exit()
  return p

#May want to write script that will fill in holes in Fd residue chain
def chain_filler(resObjects):
  lastFdResIndex = resObjects[-1].get_full_id()[3][1]
  filledChainList = [[] for a in range(0, lastFdResIndex)]
  for res in resObjects:
    Fd_index = res.get_full_id()[3][1] -1 #Must subtract one to get in with pythonic indices
    filledChainList[Fd_index] = res
  return filledChainList
  
#This function scans through the alignment sequence to look for gaps 
# and add +1 to the index of the amino acid at that position
# and each amino acid after within the parent object.
def reIndexer(parent, alignmentSequence):
  residueMax, sortingList= parent.spitSequence()

  """Need to re-index both residues and contacts.
  """
  #print([x.index for x in sortingList])
  print([x.index for x in sortingList])

  j=1#Gap counter. This will used to offset the incremental 
  # discrepancies in the sortingList index caused by the gaps.
  gapInds=[]
  for i, aa in enumerate(alignmentSequence):    
    if aa == "-":         
      gapInds.append(i)

  for i in gapInds:
    for res in sortingList[i+1-j:]: #Increase the index of every residue object after a gap is seen in the alignment.
      try:
        res.index+=1
      except TypeError:
        continue
      j+=1
  
  print([x.index for x in sortingList])
  print(sortingList[1].contacts)
  
  j=1
  for i in gapInds:
    for res in sortingList:
      for con in res.contacts:
        if con > i:
          con +=1
    j+=1    
  print(sortingList[1].contacts)



  """print([x.index for x in sortingList])
  print(parent.residues[0].index)
  print(parent.residues[0].resType)"""

  reIndexedRes=[] 
  for item in sortingList:
    if item != []:
      reIndexedRes.append(item)
  
  parent.residues=reIndexedRes

  return parent #Spit out parent with reindexed residues.

#Takes correctly indexed parents and removes the contacts from all amino acids which will not be changed by recombination.
def redundancyRemover(indexedParents, alignmentSequence1, alignmentSequence2):
  residueMax1, sortingList1= indexedParents[0].spitSequence()
  residueMax2, sortingList2= indexedParents[1].spitSequence()
  redundantIndices=[] #List that will hold the index of amino acids that are the same in each parent.
  redDict={}
  for i, aa in enumerate(alignmentSequence1):    
    if alignmentSequence1[i] == alignmentSequence2[i]:
      redundantIndices.append(i+1)#Need to add one to go from pythonic indices to pdb indices.
      redDict[i]=aa
  print(redundantIndices)
  print(redDict)
  
  for res in indexedParents[0].residues:#Get rid of contacts in all redundant residues in parent 1
    if res.index in redundantIndices:
      #print(res.resType)
      #print(res.index)
      res.contacts=[] 

  for res in indexedParents[1].residues:#Get rid of contacts in all redundant residues in parent 2
    if res.index in redundantIndices:
      #print(res.resType)
      #print(res.index)
      res.contacts=[] 

  #print(indexedParents[0].residues)
  """print(indexedParents[0].residues[0].resType)
  print(indexedParents[0].residues[0].contacts)
  print(indexedParents[0].residues[0].index)
  print(indexedParents[0].residues[-2].resType)
  print(indexedParents[0].residues[-2].contacts)
  print(indexedParents[0].residues[-2].index)"""  

  #print([x.index for x in sortingList])
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
  pdbFiles=["1rfk.pdb","pssm.pdb"]
  alignmentFile="parentAlignEd.txt"
  chimFile="chimIndices.csv"
  chims=plate_Grabber(chimFile)
  print(chims)
  
  parents=[]
  
  #Getting the sequences of the aligned parents
  align=SeqIO.parse(alignmentFile, "fasta")
  alignSequences=[]
  for record in align:
    alignSequences.append(record.seq)
  print(alignSequences)

  #Loop through several functions to add a parent object to the list.
  for i, file in enumerate(pdbFiles):
    c_dict, Fd, pName = chain_dict_maker(file) #Uses Regex to make dictionary relating chain ID to datebase name, aka {C: CYFD001} and queries user to designate Chain of interest
    c_atoms, c_res = pdb_parser(file) #Pulls out all atoms and residues from pdb file, each into single lists
    Fd_res = c_res[Fd] #Selects Fd chain from all residue chains
    #print(Fd_res)
    indDict=reIndexDictMaker(Fd_res,alignSequences[i]) #Create dictionary relating non-aligned residue to aligned residues.

    atomIndexDict, FdAtomToResDict=atomSorter(c_atoms)

    parentOb = contact_maker(Fd, pName, c_res, indDict, atomIndexDict, FdAtomToResDict) #Finds all residues that contact selected chain residues
    sys.exit()
    parents.append(parentOb)
    
    #completeFdRes = chain_filler(Fd_res) #Fills in blank lists where non-structured Fd should be in the chain
  print(len(parents))
  p1=parents[0]
  print(p1.name)
  print(p1.sequence)
  rm, sL = p1.spitSequence()
  print(len(sL))
  #for con in sL:
  #  print(con.contacts)
  sys.exit()


  
  """p1AlignSeq=list(str(sequences[0]))
  p2AlignSeq=list(str(sequences[1]))

  p1=reIndexer(parents[0],p1AlignSeq)"""

  """p2=reIndexer(parents[1],p2AlignSeq)
  nonRedParents= redundancyRemover([p1,p2],p1AlignSeq,p2AlignSeq)

  p1.alignSeq=p1AlignSeq
  p2.alignSeq=p2AlignSeq
  schemaCalc(nonRedParents[0], nonRedParents[1])
  #p1Chims, p2Chims=schemaCalc(nonRedParents[0], nonRedParents[1])"""

  """file_writer(p1Chims,p2Chims)"""
  

if __name__ == '__main__':
    main() 