import numpy as np



#Precisão: P=TP/TP+FP | TP: número de contatos preditos corretamente | TP+FP: número total de contatos preditos
def precision(pairsPredicted, numberContactPredicted,pairsNative, numberContactNative):

	pairCommon = set(pairsPredicted).intersection(pairsNative)

	TP = len(set(pairsPredicted).intersection(pairsNative))

	precision = TP/numberContactPredicted

	print("predicted",pairsPredicted)
	print("native",pairsNative)

	print("pairCommon",pairCommon)
	
	print("precision", round(precision,2))

	return round(precision,2)

#Cobertura: | Cov = 100 x TP / número de contatos nativos
def coverage(pairsPredicted,pairsNative,numberContactNative):
	
	TP = len(set(pairsPredicted).intersection(pairsNative))

	cov = (100 * TP)/numberContactNative

	print("cov",round(cov,2))

	return round(cov,2)