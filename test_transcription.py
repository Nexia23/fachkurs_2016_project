import molecules as mol
import transcription as trsc
import random



dic={}
nuc=['A','C','G','T']

def rand_genedic(dic):
	for gene in range(10):
		id=random.randint(1,100)
		name="Gene_{}".format(gene)
		seq=''
		for i in range(id):
			n=random.randint(0,3)
			seq+=nuc[n]
		dic.update({id: mol.Gene(id, name, seq, 1)})
	return dic

dic= rand_genedic(dic)


#mygene=mol.Gene(9875569, 'supergene', 'ACGGGGGTGTGACC',1)
mypol=mol.RNAPolymeraseII(2830, 'bigpol', 42)

trnsc=trsc.Transcription(894, 'perfect_transc', mypol)
rna_pool=[]

for s in range(1000):
	#print(mygene.sequence[s])
	rna_pool= trnsc.update(dic, rna_pool)

print(len(rna_pool))

for r in rna_pool:
	print(r.name+": "+r.sequence)

#print('\n')
#for g in dic.keys():
#	print(dic[g].name+" :"+dic[g].sequence)
#print(mygene.sequence_binding)


