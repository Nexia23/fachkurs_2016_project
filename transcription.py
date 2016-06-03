import processes 
import molecules
import random
import numpy as np
from matplotlib import pyplot as plt


################### TO DO !!!!!! ########################
#look to data: elongation-processes per 1s and transcripion-rates (Huyen)

#include real data and global data-construction

#get information how many mRNAs are transcribed at the moment

#at the moment: only one polymerase transcribes an RNA ! (Paula)

#include rRNA and tRNA - transcription

#######################################################

#mRNA half lifes

########################################################

class Transcription(processes.Process):
	"Transcription process needs to be described."


	def __init__(self, id, name):
		""" All new Transcription variables need to be implemented."""
		
		super().__init__(id, name)
		#initialized through set_states in model: 
		#id_enzymes= one single RNA Polymerase -> not necessary as polymerase is an argument of Transcription 
		# -> multiple transcription processes with different polymerases possible
		#id_substrates= a dictionary with all genes of the genome: {geneid1: Geneobject, geneid2: GeneObject, ...} 

		#self.mypolymerase=polymerase
		#self.my_nucleotides=nucleotide_pool

				#polymerase=RNAPolymeraseII(needs to be specified) -> initialization can be done in transcription or model
				#expect polymerase to be unbound -> check in model.py

		#length of occupied sequence by RNA-Polymerase II on DNA: already validated with data!
		self.polymerase_size=17

		#####for visualization of selected genes
		#self.allgenes=[[],[]]
		#########################################


	def __repr__(self):
		# todo: each process class should have something like this
		pass

	def __str__(self):
		# todo: each process class should have something like this
		pass



	def update(self, model):
	#def update(self, genedic, rna_pool):

		genedic=model.genes
		#rna_pool=model.mrnas
		self.mypolymerase=self.enzyme_ids
		self.my_nucleotides=model.nucleotides

		#number still needed: empirically! 
		update_per_s=500

		gene_id, weights=self.make_weights(genedic)


		for steps in range(update_per_s):
			#rna = self.onestep(genedic)

			#print('New step: \n')	

			bound_gene=False

			#elongation
			for g in gene_id:
				g=genedic[g]
				for i in g.pol_on_gene[0]:
					bound_gene=True
					if random.randint(1,10)<9:
						#print('elongate!')
						rna = self.elongation(g, model.chromosomes, i)
						
						if isinstance(rna, molecules.MRNA):
					#rna_pool.append(rna)
							if rna.name in model.states:
								model.states[rna.name].append(rna)
							else:
								model.states[rna.name] = [rna]

			transc_gene=self.select_gene(genedic)
			while len(transc_gene.location) >1:
				transc_gene=self.select_gene(genedic)

			if not transc_gene.pol_on_gene[0]: 

				if bound_gene==False:
					self.intialization(transc_gene, model.chromosomes, weights, gene_id)
				elif bound_gene==True and random.randint(1,40)==20:
					self.intialization(transc_gene, model.chromosomes, weights, gene_id)

			##2. initialization?
			else:
				self.intialization(transc_gene, model.chromosomes, weights, gene_id)

			#print('Step ended \n')
				

    	####visualization of selected genes ######
		#plt.plot(range(len(self.allgenes[0])),self.allgenes[1])
		#plt.xlabel(self.allgenes[0])
		#plt.savefig('tests/count_histogram_genes.pdf')
		###########################################

		#return rna_pool

	def elongation(self, transc_gene, all_chromosomes, pos):

		if int(transc_gene.location[0][0])< int(transc_gene.location[0][1]):
			strand='+'
		else:
			strand='-'

			#get chromosome-position: on which chromosome are we? 
		for chr in list(all_chromosomes.keys()):
			if transc_gene.chr==all_chromosomes[chr].id:
				mychr=all_chromosomes[chr]
				break

		pol_position=pos
		mrna = self.transcribe(transc_gene, pol_position, mychr, strand)
		if isinstance(mrna, molecules.MRNA):
			return mrna

			
	#def onestep(self,genedic):			#expect a dictionary of genes 
	def intialization(self, transc_gene, all_chromosomes, weights, gene_ids):
		chr_found=False

		for chr in list(all_chromosomes.keys()):
			if transc_gene.chr==all_chromosomes[chr].id:
				mychr=all_chromosomes[chr]
				break
		
			

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
		try:
			if mychr.chromosome_bound((int(transc_gene.location[0][0])-self.polymerase_size,int(transc_gene.location[0][0])+self.polymerase_size)) is None and self.mypolymerase.count>0:
				#print('Initiate!')
				mychr.bind_to_chrom( (int(transc_gene.location[0][0])-self.polymerase_size,int(transc_gene.location[0][0])+self.polymerase_size), self.mypolymerase)
				""" for one Polymerase: look for a gene to transcribe. when ORF is found: initiate """

				#### from the data group: we expect an self.genes-dictionary in model.py containing all genes

		
				if isinstance(self.mypolymerase, molecules.RNAPolymeraseII): 
					transc_gene.pol_on_gene[0].append(int(transc_gene.location[0][0]))
					transc_gene.pol_on_gene[1].append(molecules.MRNA("mRNA_{}".format(transc_gene.mid), "mRNA_{0}".format(transc_gene.name), '',))
			
			
				self.mypolymerase.count+=-1
	
		except UnboundLocalError:
			pass




				

	def make_weights(self, genedic):
		transc_rate=[]
		gene_ids=list(genedic.keys())
		for g in gene_ids:
			transc_rate.append(genedic[g].transrate)
		transc_rate=np.array(transc_rate)

		weights=transc_rate
		weights=weights/sum(weights)
		return (gene_ids, weights)



	def select_gene(self, genedic):

		#transcribed gene identified: no random selection, but weighted by copies of each gene

		copies=[]
		transc_rate=[]
		weight_for_binding=[]
		gene_ids=list(genedic.keys())
		for g in gene_ids:
			copies.append(genedic[g].count)
			transc_rate.append(genedic[g].transrate)
			if genedic[g].pol_on_gene:
				weight_for_binding.append(1000000000)
			else:
				weight_for_binding.append(1)
		copies=np.array(copies)
		transc_rate=np.array(transc_rate)

		weights=copies*transc_rate*weight_for_binding
		#weights=copies*transc_rate
		weights=weights/sum(weights)

		rand_index=np.random.choice(range(len(weights)),p=weights)

			
		transc_gene=genedic[gene_ids[rand_index]]

		return transc_gene


	#def initiate(self,gene):
#moved to initialization!!!!!!

	def transcribe(self, gene, pos, mychr, strand):

		""" elongate mRNA for given ORF for only one step. if ORF ends: call terminate-function """

		for i in range(20):
	
			index=gene.pol_on_gene[0].index(pos)

			mrna=gene.pol_on_gene[1][index]

			#appending the correct codon from the coding strand
			
			#polymerase moves forward! in which direction to move? (strand-dependend)
			if strand=='+':

			#if we are not on the end of the ORF-string
				if pos+1<int(gene.location[0][1]):

					if mychr.chromosome_bound(pos+self.polymerase_size+1) is None:
						nuc=gene.sequence[pos-int(gene.location[0][0])]

						if nuc=='T':
							mrna.sequence+='U'
							self.my_nucleotides.count_nuc['U']+=-1
						else:
							mrna.sequence+=nuc
							self.my_nucleotides.count_nuc[nuc]+=-1

						mychr.del_on_chrom((pos-self.polymerase_size, pos+self.polymerase_size))
						mychr.bind_to_chrom((pos-self.polymerase_size+1, pos+self.polymerase_size+1), self.mypolymerase)
						gene.pol_on_gene[0][index]+=1
						pos+=1

				else:
					nuc=gene.sequence[int(gene.location[0][0])-pos]


					if nuc=='T':
						mrna.sequence+='U'
						self.my_nucleotides.count_nuc['U']+=-1
					else:
						mrna.sequence+=nuc
						self.my_nucleotides.count_nuc[nuc]+=-1

					mychr.del_on_chrom((pos-self.polymerase_size, pos+self.polymerase_size))
					return self.terminate(gene, pos)

			else:
				if pos-1>int(gene.location[0][1]):

					if mychr.chromosome_bound(pos-self.polymerase_size-1) is None:
						nuc=gene.sequence[int(gene.location[0][0])-pos]


						if nuc=='T':
							mrna.sequence+='U'
							self.my_nucleotides.count_nuc['U']+=-1
						else:
							mrna.sequence+=nuc
							self.my_nucleotides.count_nuc[nuc]+=-1

						mychr.del_on_chrom((pos-self.polymerase_size, pos+self.polymerase_size))
						mychr.bind_to_chrom((pos-self.polymerase_size-1, pos+self.polymerase_size-1), self.mypolymerase)
						gene.pol_on_gene[0][index]+=-1
						pos +=-1
				else:

					nuc=gene.sequence[int(gene.location[0][0])-pos]


					if nuc=='T':
						mrna.sequence+='U'
						self.my_nucleotides.count_nuc['U']+=-1
					else:
						mrna.sequence+=nuc
						self.my_nucleotides.count_nuc[nuc]+=-1

					mychr.del_on_chrom((pos-self.polymerase_size, pos+self.polymerase_size))
					return self.terminate(gene, pos)
			
		return 0



	def terminate(self, gene, position):
		""" separate mRNA-DNA-Polymerase-complex. release and store new mRNA """
		#print(gene.pol_on_gene)
		self.mypolymerase.count+=1
		rna = gene.pol_on_gene[1][gene.pol_on_gene[0].index(position)]
		
		del gene.pol_on_gene[1][gene.pol_on_gene[0].index(position)]
		del gene.pol_on_gene[0][gene.pol_on_gene[0].index(position)]
		

		return rna


	def rand_distr(self, gene, gene_ids, weights): #if gene_ids and weight are global we don't need to send them in the function 
		ran=random.gauss(0.1,0.05)
		index_id = gene_ids.index(gene.mid)
		if weights[index_id]>= ran:
			return 1
		else:
			return 0

