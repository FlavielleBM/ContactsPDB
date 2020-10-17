import numpy as np
import itertools
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB
from Bio.PDB import *

from evaluation import precision,coverage

def getPairs(short_range,medium_range,long_range):
	pairs = []

	pairs =  short_range +  medium_range + long_range

	print(pairs)

	return pairs


def getNumberContactPairs(short_range,medium_range,long_range):

	numberContact =  len(short_range) + len(medium_range) + len(long_range) 

	print(numberContact)

	return numberContact

def getNumberContact(lengthAA):
	L = lengthAA

	print(L, round(L/2 - 0.5),round(L/5 - 0.5), round(L/10- 0.5))

	return L, round(L/2 - 0.5),round(L/5 - 0.5), round(L/10- 0.5)

def contactRangeNativ(structure_pdb, diff):

	short_range =[]
	medium_range =[]
	long_range =[]

	parser = PDBParser()
	structure = parser.get_structure("structure", structure_pdb)  

	model = structure[0]

	residues = [r for r in model.get_residues() if r.get_id()[0] == " "]
	for each in itertools.combinations(residues, 2):
		if each[0].get_resname() == 'GLY':
			atom1 = each[0]['CA'].get_coord()	
		#	print("atom1",atom1)
		else:
			atom1 =each[0]['CB'].get_coord()
		#	print("atom1",atom1)
		
		if each[1].get_resname() == 'GLY':
			atom2 = each[1]['CA'].get_coord()	
		#	print("atom2",atom2)
		
		else:
			atom2 =each[1]['CB'].get_coord()
		#	print("atom2",atom2)

		res1 = each[0].get_id()[1]
		res2 = each[1].get_id()[1]

		if diff != 1:
			res1 = res1 - diff
			res2 = res2 - diff

		
		residues_pairs = res1,res2
		
		distance = np.linalg.norm(atom1-atom2)
		if distance < 8:
			if res2-res1 < 6:
				short_range.append(residues_pairs)
			if res2-res1 > 11 and res2-res1 < 24:
				medium_range.append(residues_pairs)
			if res2-res1 > 24:
				long_range.append(residues_pairs)

	#print("short", short_range)
	#print("medium", medium_range)
	#print("long", long_range)

	return short_range,medium_range,long_range

              

def getPositionAA(protein):
	positionAA = [] 

	parser = PDBParser()
	structure = parser.get_structure("structure", protein)   

	model = structure[0]

	residues = [r for r in model.get_residues() if r.get_id()[0] == " "]
	for each in itertools.combinations(residues, 1):
		
		positionAA.append(each[0].get_id()[1])

	lengthAA = len(positionAA)

	#menos um para o resíduo começar em 1 e não em 0
	diff = positionAA[0] - 1
	print(positionAA,lengthAA,diff)

	return positionAA, lengthAA, diff


def contactRangePredicted(pairs):
	short_range =[]
	medium_range =[]
	long_range =[]

	for i in pairs:
		distance =  i[1]- i[0]
	#	print(distance)
		if distance < 6:
			short_range.append(i)
		if distance > 11 and distance < 24:
			medium_range.append(i)
		if distance > 24:
			long_range.append(i)

#	print("short_predicted", short_range)
#	print("medium_predicted", medium_range)
#	print("long_predicted", long_range)

	return short_range,medium_range,long_range

def getNumberContactPredicted(short_range,medium_range,long_range):
	numberContact = len(short_range) +  len(medium_range) + len(long_range) 

	print("number contact predicted",numberContact)

	return numberContact


def get_contact(file_contact):
	pairs = []
	
	file = open(file_contact,"r")

	for i in file:
		information = str(i)
		values = information.split()
		res1=int(values[0])
		res2=int(values[1])
		probability= float(values[4])
		if probability > 0.3:
			residues_pairs = res1,res2
			pairs.append(residues_pairs)

	#print("Pares", pairs)	
	return pairs


def contactNative(structure_pdb):
	#get pairs contact native
	positionAA, lengthAA, diff = getPositionAA(structure_pdb)

	L,L_2,L_5,L_10 = getNumberContact(lengthAA)

	short_range_nativ, medium_range_nativ, long_range_nativ = contactRangeNativ(structure_pdb, diff)

	numberContactNative = getNumberContactPairs(short_range_nativ, medium_range_nativ, long_range_nativ)

	pairsNative = getPairs(short_range_nativ,medium_range_nativ,long_range_nativ)


	return pairsNative, numberContactNative

def contactPredicted(file_contact):

	#get pairs contact predicted

	pairsPredicted = get_contact(file_contact)

	short_range_pred, medium_range_pred, long_range_pred = contactRangePredicted(pairsPredicted)

	numberContactPredicted = getNumberContactPairs(short_range_pred, medium_range_pred, long_range_pred)

	return pairsPredicted, numberContactPredicted

def run():

	structure_pdb = "/home/acer-linux/Documentos/CReF_Set/1zdd.pdb"

	file_contact = '/home/acer-linux/Documentos/ContactsPDB/1zdd.rr'
	
	sequence = 'FNMQCQRRFYEALHDPNLNEEQRNAKIKSIRDDC'
	
	pairsNative, numberContactNative = contactNative(structure_pdb)

	pairsPredicted, numberContactPredicted = contactPredicted(file_contact)


	p = precision(pairsPredicted, numberContactPredicted, pairsNative, numberContactNative)

	cov = coverage(pairsPredicted,pairsNative,numberContactNative)




run()