import processes 
import molecules
import random
import numpy as np
from matplotlib import pyplot as plt


################### TO DO !!!!!! ########################
#respect gene count in binding probability  ! (Lea)

#look to data: elongation-processes per 1s and transcripion-rates (Huyen)

#get information how many mRNAs are transcribed at the moment

#include correct gene-class: with gene start, end and strand instead of sequence ! (Paula)

#at the moment: only one polymerase transcribes an RNA ! (Paula)

#include rRNA and tRNA - transcription

#bigger polymerase, occupying more than only the transcribed position ! (Lea)

#######################################################

#mRNA half lifes

########################################################

class Transcription(processes.Process):
	"Transcription process needs to be described."


	def __init__(self, id, name, polymerase):
		""" All new Transcription variables need to be implemented."""
		
		super().__init__(id, name)
		#initialized through set_states in model: 
		#id_enzymes= one single RNA Polymerase -> not necessary as polymerase is an argument of Transcription 
		# -> multiple transcription processes with different polymerases possible
		#id_substrates= a dictionary with all genes of the genome: {geneid1: Geneobject, geneid2: GeneObject, ...} 
		self.mypolymerase=polymerase

				#polymerase=RNAPolymeraseII(needs to be specified) -> initialization can be done in transcription or model
				#expect polymerase to be unbound -> check in model.py

		#####for visualization of selected genes
		#self.allgenes=[[],[]]
		#########################################


	def __repr__(self):
		# todo: each process class should have something like this
		pass

	def __str__(self):
		# todo: each process class should have something like this
		pass



	#def update(self, model):
	def update(self, genedic, rna_pool):
		#genedic=model.genes
		#rna_pool=model.states

		update_per_s=5000

		for steps in range(update_per_s):
			rna = self.onestep(genedic)
			if isinstance(rna, molecules.MRNA):
				rna_pool.append(rna)
				#if rna.name in model.states:
    			#	model.states[rna.name].append(rna)
    			#else:
    			#	model.states[rna.name] = [rna]

    	####visualization of selected genes ######
		#plt.plot(range(len(self.allgenes[0])),self.allgenes[1])
		#plt.xlabel(self.allgenes[0])
		#plt.savefig('tests/count_histogram_genes.pdf')
		###########################################

		return rna_pool


			
	def onestep(self,genedic):			#expect a dictionary of genes 

		transc_gene=self.select_gene(genedic)

		#testprints
		#print(transc_gene.name)
		#print(transc_gene.sequence)	
		#print(transc_gene.pol_on_gene)

		######for visualization of selected genes ################
		#if transc_gene.name in self.allgenes[0]:
		#	self.allgenes[1][self.allgenes[0].index(transc_gene.name)]+=1
		#else:
		#	self.allgenes[0].append(transc_gene.name)
		#	self.allgenes[1].append(1)
		#########################################################

		if not transc_gene.pol_on_gene:
			self.initiate(transc_gene)
		else:
			#if gene very long - likely to have a second pol binding there
				#self.initiate(transc_gene)
			pol_position=random.choice(transc_gene.pol_on_gene)
			mrna = self.transcribe(transc_gene, pol_position)
			if isinstance(mrna, molecules.MRNA):
				#print('this is an mrna')
				#if mrna.name in model.states:
				#	model.states[mrna.name].append(mrna)
				#else:
				#	model.states[mrna.name] = [mrna]
				return mrna


	def select_gene(self, genedic):

		#trandscribed gene identified
		#include parameter gene copies
		copies=[]
		for g in genedic.keys():
			copies.append(genedic[g].count)
		copies=np.array(copies)
		copie_probs=copies/sum(copies)

		rand_index=np.random.choice(range(len(copies)),p=copie_probs)

		gene_ids=list(genedic.keys())	#no random selection, but weighted by copies of each gene
		transc_gene=genedic[gene_ids[rand_index]]

		return transc_gene


	def initiate(self,gene):

		""" for one Polymerase: look for a gene to transcribe. when ORF is found: initiate """

		#### from the data group: we expect an self.genes-dictionary in model.py containing all genes

		#print(gene.sequence_binding)

		if gene.sequence_binding[0]==0 and self.mypolymerase.count>0:
		
		#optional: add stochastical condition for binding, combination with gene transcription rates 
			gene.pol_on_gene.append(0)
			gene.sequence_binding[0]=molecules.MRNA("mRNA_{}".format(gene.mid), "mRNA_{0}".format(gene.name.split("_")[-1]), '',)		
			#later: add bigger polymerase -> sequence_binding[0+i]=1
			self.mypolymerase.count+=-1


	def transcribe(self, gene, position):

		""" elongate mRNA for given ORF for only one step. if ORF ends: call terminate-function """

		pos=position
		index=gene.pol_on_gene.index(position)
		mrna=gene.sequence_binding[pos]

		#appending the correct codon from the coding strand
		nuc=gene.sequence[pos]
		if nuc=='T':
			mrna.sequence+='U'
		else:
			mrna.sequence+=nuc

		if pos+1<len(gene.sequence_binding):

			if gene.sequence_binding[pos+1]==0:
				gene.pol_on_gene[index]+=1
				gene.sequence_binding[pos+1]=mrna
				gene.sequence_binding[pos]=0
				return 0
		else:
			return self.terminate(gene, position)
		


	def terminate(self, gene, position):
		""" separate mRNA-DNA-Polymerase-complex. release and store new mRNA """
		
		self.mypolymerase.count+=1

		mrna = gene.sequence_binding[position]
		gene.sequence_binding[position]=0
		del gene.pol_on_gene[gene.pol_on_gene.index(position)]

		return mrna


		

