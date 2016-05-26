import processes 
import molecules
import random


################### TO DO !!!!!! ########################
#respect gene count in binding probability 

#implement an update function which refers to 1s of cell cycle
# -> e.g. bigger update function iterating transcribe-steps for the genes
# -> e.g. look to data: elongation-processes per 1s and transcripion-rates

#get information how many mRNAs are transcribed at the moment



########################################################

class Transcription(processes.Process):
	"Transcription process needs to be described."


	def __init__(self, id, name, polymerase):
		""" All new Transcription variables need to be implemented."""
		
		super().__init__(id, name)
		#initialized through set_states in model: 
		#id_enzymes= one single RNA Polymerase
		#id_substrates= a dictionary with all genes of the genome: {geneid1: Geneobject, geneid2: GeneObject, ...} 
		self.mypolymerase=polymerase
		self.position=None
		self.gene=0
				#polymerase=RNAPolymeraseII(needs to be specified) -> initialization can be done in transcription or model
				#expect polymerase to be unbound -> check in model.py


	#def update(self, model):		#expect an dictionary genes 
	def update(self,genedic):
		
		#check for polymerase-dna-binding
		if self.position==None:
			#allgenes=model.genes

			#trandscribed gene identified
			#gene_ids=model.genes.keys()
			gene_ids=list(genedic.keys())
			rand_index=random.randint(0,len(gene_ids))

			#transc_gene=model.genes[rand_index]
			transc_gene=genedic[gene_ids[rand_index]]
			print(transc_gene.sequence)


			self.initiate(transc_gene)
		else:
			mrna = self.transcribe()
			if isinstance(mrna, molecules.MRNA):
				#if mrna.name in model.states:
				#	model.states[mrna.name].append(mrna)
				#else:
				#	model.states[mrna.name] = [mrna]
				return mrna



	def initiate(self,gene):

		""" for one Polymerase: look for a gene to transcribe. when ORF is found: initiate """

		#### from the data group: we expect an self.genes-dictionary in model.py containing all genes

		if gene.sequence_binding[0]==0 and self.mypolymerase.count>0:
		
		#optional: add stochastical condition for binding
			self.position=0
			self.gene=gene
			gene.sequence_binding[0]=molecules.MRNA("mRNA_{}".format(gene.mid), gene.name, '',)		#later: add bigger polymerase -> sequence_binding[0+i]=1
			self.mypolymerase.count+=-1


	def transcribe(self):

		""" elongate mRNA for given ORF for only one step. if ORF ends: call terminate-function """

		pos=self.position
		gene=self.gene
		mrna=gene.sequence_binding[pos]

		#appending the correct codon from the coding strand
		nuc=gene.sequence[pos]
		if nuc=='T':
			mrna.sequence+='U'
		else:
			mrna.sequence+=nuc

		if pos+1<len(gene.sequence_binding):

			if gene.sequence_binding[pos+1]==0:
				self.position+=1
				gene.sequence_binding[pos+1]=mrna
				gene.sequence_binding[pos]=0
				return 0
		else:
			return self.terminate()
		




	def terminate(self):
		""" separate mRNA-DNA-Polymerase-complex. release and store new mRNA """
		
		
		self.mypolymerase.count+=1

		mrna = self.gene.sequence_binding[self.position]
		self.gene.sequence_binding[self.position]=0
		self.position=None
		self.gene=0	

		return mrna


		

