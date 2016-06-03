import molecules as mol
import transcription as trsc
import random



dic={}
nuc=['A','C','G','T']

def rand_genedic(dic):
	for gene in range(20):
		id=random.randint(20,100)
		#id=10
		name="Gene_{}".format(gene)
		seq=''
		for i in range(id):
			n=random.randint(0,3)
			seq+=nuc[n]
		dic.update({id: mol.Gene(id, name, 'chr1', seq, (234,85443), 45, random.randint(1,2))})
	return dic

dic= rand_genedic(dic)


#mygene=mol.Gene(9875569, 'supergene', 'ACGGGGGTGTGACC',1)
mypol=mol.RNAPolymeraseII(2830, 'bigpol', 42)
nucleotides=mol.NucleotidPool(235,'nuc', 100000)

trnsc=trsc.Transcription(894, 'perfect_transc', mypol, nucleotides)
rna_pool=[]

for s in range(1):
	#print(mygene.sequence[s])
	rna_pool= trnsc.update(dic, rna_pool)

print(len(rna_pool))

#for r in rna_pool:
#	print(r.name+": "+r.sequence)

#print('\n')
for g in dic.keys():
	print(dic[g].name+" :"+str(dic[g].count))
#print(mygene.sequence_binding)

#print(nucleotides.count_nuc)

#chr_list = createchromosomes()
#chr=chr_list[0]
#chr.binding_molecules=[[(234,347),(348,956),(1039,1299)],['pol1','pol2','pol3','replicase']]
#print(chr.chromosome_bound([957,1022]))