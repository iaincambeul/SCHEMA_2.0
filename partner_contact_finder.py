import sys
import re
import os
import csv
from Bio.PDB import *
from Bio import SeqIO

global letterToNumber_dict
letterToNumber_dict = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5, "G": 6, "H": 7, "I": 8, "J": 9, "K": 10, "L": 11, "M": 12, "N" : 13, "O" : 14, "P" : 15, "Q" : 16, "R": 17, "S": 18, "T": 19, "U": 20, "V": 21, "W": 22, "X": 23, "Y": 24, "Z": 25}

"""Sometimes, when a chain lacks structural information, the indices of fdResPartnerObjects will not match the index of the actual Fd residues, leading to errors when
I try to pass those indices in FdAtomToResDict.  Need lengthen fdResPartnerObjects"""

def chain_dict_maker(pdb_file):
  pdb_open = open(pdb_file, "rU")
  chain_list = []
  chain_list = re.findall("DBREF\s+\w+\s+(\w)\s+\w+\s+\w+\s+\w+\s+\w+\s+(\w+)", pdb_open.read())
  print(chain_list)
  chain_dict = {}
  for chain in chain_list:
    chain_dict[chain[0]] = chain[1]
  print(chain_dict)
  pdb_open.close()
  Fd_choice = input("Which uppercase letter corresponds to the Fd in the list (A, B, C, etc.)? ")
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

def contact_maker(Fd_ind, chain_atoms, chain_dict, res_list):
  #Need to combine all lists of not selected Fd atoms
  nonFdAtoms = []
  for ind, chain in enumerate(chain_atoms):
    if ind == Fd_ind:
      print("Fd_ind:" + str(ind))
    else:
      nonFdAtoms.append(chain_atoms[ind])
      print(ind)
  nonFdAtoms = [item for sublist in nonFdAtoms for item in sublist] #Collapses all sublists into a single list.  This puts all atoms from non-Fd chains into a single list

  #Creates NeighborSearch object with Fd_atoms list. Generates the object only.  Want to put in the atoms that Fd atoms will be contacting
  neighbor_searcher = NeighborSearch(nonFdAtoms)
  
  ind=0 #This is just to help me keep track of which atom the machine is on and can be deleted later
  Fd_atoms = chain_atoms[Fd_ind] #Selecting Fd atoms from chain atoms list
  
  #Creating dictionary so that I can sort contacting partner residues 
  #once I have them in them in the atom indexed list
  FdAtomToResDict = {}
  for ind, atom in enumerate(Fd_atoms):
    FdAtomToResDict[ind] = atom.get_full_id()[3][1]-1 #Have to subtract one to get into Pythonic indices

  fdAtomPartnerObjects = [] # This is the atom-based list we will put the Fd partner contacts into
  #Making a contact map with a 4.5 angstrom cutoff
  for i in range(0, len(Fd_atoms)): #Going through all Fd atoms
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
  lastFdResIndex = chain_atoms[Fd_ind][-1].get_full_id()[3][1]
  fdResPartnerObjects = [[] for a in range(0, lastFdResIndex)]# This is the residue-based list we will put the Fd partner contacts int
  print("Fd Res count:" + str(len(fdResPartnerObjects)))
    
  index = -1
  #This nested for loop will add partner residue objects to the appropriate
  #Fd residue indices using the FdAtomToResDict dictionary as a key
  for atom in fdAtomPartnerObjects:
    index += 1
    for item in atom:
      #print "FdAtomToResDict[index]:" + str(FdAtomToResDict[index])
      #print "fdResPartnerObjects[FdAtomToResDict[index]]:" + str(fdResPartnerObjects[FdAtomToResDict[index]])
      if item not in fdResPartnerObjects[FdAtomToResDict[index]]:
        fdResPartnerObjects[FdAtomToResDict[index]].append(item)

  return fdResPartnerObjects

#May want to write script that will fill in holes in Fd residue chain
def chain_filler(resObjects):
  lastFdResIndex = resObjects[-1].get_full_id()[3][1]
  filledChainList = [[] for a in range(0, lastFdResIndex)]
  for res in resObjects:
    Fd_index = res.get_full_id()[3][1] -1 #Must subtract one to get in with pythonic indices
    filledChainList[Fd_index] = res
  return filledChainList
  
#Parses object list data from partner_contact_finder into a csv printable format
def object_parser(partnerContacts, ferredoxinRes, chain_dict):
  res_length = len(partnerContacts) #Should get the number of residues in the Fd chain, including those that don't have pdb data
  bothRes = [[] for a in range(0, res_length)]  
  for res in ferredoxinRes:
    if res != []:
      fdResIndex = res.get_full_id()[3][1]-1 #pulls out the residue index for each Fd residue
      bothRes[fdResIndex] = ([ferredoxinRes[fdResIndex], partnerContacts[fdResIndex]]) #Places Fd residues and contacting parnter residues in the right indices

  #Using list comprehension to remove all residue pairs that are: 1. an empty list and 2. Have an empty list in the second element (No contacting residues)
  cleanedBothRes = [x for x in bothRes if x!= [] and x[1] != []]
  
  #Splitting up pairs in the cases where one Fd residue contacts multiple partners
  splitCleanedBothRes= []
  for pair in cleanedBothRes:
    if len(pair[1]) > 1:
      for ele in pair[1]:
        splitCleanedBothRes.append([pair[0], ele])
    else:
      splitCleanedBothRes.append([pair[0]] + pair[1])#Add lists and make Fd object a list to make list depth congruent with multiple partner case

  textLines = []
  for pair in splitCleanedBothRes:
    fd_number = str(pair[0].get_full_id()[3][1])
    partner_number = str(pair[1].get_full_id()[3][1])
    fd_chain_id = str(pair[0].get_full_id()[2])
    partner_chain_id = str(pair[1].get_full_id()[2])
    fd_res_type = str(pair[0].get_resname())
    partner_res_type = str(pair[1].get_resname())
    fd_chain_name = str(chain_dict[fd_chain_id])
    partner_chain_name = str(chain_dict[partner_chain_id])
    textLines.append([fd_chain_name+"-"+fd_chain_id+"-"+fd_number, partner_chain_name+"-"+partner_chain_id+"-"+partner_number, fd_res_type, partner_res_type]) 
  return textLines

def file_writer(name, data):
  final_name = name.replace(".pdb","") + "_fd_partner_contacts.csv"
  with open(final_name,'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["Structure:" + name])
    writer.writerow(["Interactions:"])
    writer.writerow(["Fd Res#", "Partner Res#", "Fd Res Type", "Partner Res Type"])
    for interaction in data:
      writer.writerow(interaction)

def main():
  c_dict,Fd = chain_dict_maker(sys.argv[1]) #Uses Regex to make dictionary relating chain ID to datebase name, aka {C: CYFD001} and queries user to designate Chain of interest
  c_atoms, c_res = pdb_parser(sys.argv[1]) #Pulls out all atoms and residues from pdb file, each into single lists
  fdPartnerContacts = contact_maker(Fd, c_atoms, c_dict, c_res) #Finds all residues that contact selected chain residues
  Fd_res = c_res[Fd] #Selects Fd chain from all residue chains
  completeFdRes = chain_filler(Fd_res) #Fills in blank lists where non-structured Fd should be in the chain
  fileName = str(sys.argv[1])
  
  printable_data = object_parser(fdPartnerContacts, completeFdRes, c_dict)
  file_writer(fileName, printable_data)
  
if __name__ == '__main__':
    main() 