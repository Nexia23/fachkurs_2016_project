from molecules import BioMoleculeCount
from molecules import BioMolecule
from processes import Process


class Chromosome:
    def __init__(self, name, sequence):
        self.sequence = sequence
        self.name = name
        self.replication_ori_bound = False
#BioM class PM, He
# wenn He gebunden, PM bindet (if bound...)


#mehrere Origins
#definition i?

class Helicase(BioMoleculeCount):
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)

    def __str__(self):
        return 'helicases: {0}'.format(self.count)
        
    #Bindung?
    # bindet an Stelle X
    
class Polymerase(BioMoleculeCount):
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)

    def __str__(self):
        return 'polymerase: {0}'.format(self.count, self.name)
  

class Replication(Process):
    def __init__(self, id, name):
        super().__init__(id, name)
        self.chromosomes = {}
        self.duplication = {}
        #ruft hel+PM auf
        #self._init_helicase()
        #self._init_polymerase()
    
    # hel+PM werden definiert

    
    def update(self, model):
        #call state (> transcription?)
        self.chromosome_names = self.substrate_ids
        print(self.chromosome_names)
        self.polymerase =  model.states['Polymerase3']
        self.helicase = model.states['DnaB']
        
        self.old_chromosomes = [model.states[chromsome_name] for chromsome_name in self.chromosome_names]

        for old_chromosome in self.old_chromosomes:
            if not old_chromosome.replication_ori_bound and not self.duplication[chrom.name] and self.polymerase.count > 0 and self.helicase.count > 0:
                self.initiate(old_chromosome)
            elif old_chromosome.replication_ori_bound:
                self.elongate(old_chromosome)
            elif (self.polymerase.count == 0 or self.helicase.count == 0) and not old_chromosome.replication_ori_bound:
                continue
            else:
                raise NotImplementedError
        print(self.helicase.count)
        print(self.polymerase.count)
        print(self.chromosomes['chrom1'].sequence)
        print(len(self.chromosomes['chrom1'].sequence))
    
    def initiate(self, chrom: Chromosome):
        print('initiation')
            #Zeitpunkt der Replikation?
        
        

            #wenn nichts gebunden ist, dann:
        
                #wird ein neues Chromosom erstellt, ist noch leer
             #Listenform für Sequenz?
        self.chromosomes[chrom.name]=Chromosome(chrom.name,[])
        self.polymerase.count -= 1 #jew count -1 für jedes gebundene molekül
        self.helicase.count -= 1
        chrom.replication_ori_bound = True
        self.duplication[chrom.name] = False
  
    def elongate(self, old_chrom: Chromosome):
        print('elongation')
        #mehrere Schritte (Zeitschritte) >listeneinträge new_chrom  >>vgl alte sequ
            # dictionary 'chromosomes' enthält 'Chromosomname':'Eigenschaften(Name, Sequenz')
        
    
        #definition i?
        new_chrom = self.chromosomes[old_chrom.name]
        

        if len(new_chrom.sequence) == 0:
            sequence_to_replicate = old_chrom.sequence[0:100]
        elif len(old_chrom.sequence) - len(new_chrom.sequence) > 100:
            sequence_to_replicate = old_chrom.sequence[len(new_chrom.sequence)+1:len(new_chrom.sequence)+100]
        elif len(old_chrom.sequence) - len(new_chrom.sequence) == 0:
            return self.terminate(new_chrom.sequence)
        else:
            sequence_to_replicate = old_chrom.sequence[len(new_chrom.sequence) + 1:]

              

        for x in sequence_to_replicate:
            if x == 'A':
                    #soll die Sequence von new_chrom befüllen
                new_chrom.sequence.append('T')      #>>>???
            if x == 'T':
                new_chrom.sequence.append('A')
            if x == 'G':
                new_chrom.sequence.append('C')
            if x == 'C':
                new_chrom.sequence.append('G')

        
    def terminate(self, chrom, i):

            #wenn die Länge der neuen Sequence der länge der alten entspricht:    
        
        self.polymerase += 1
        self.helicase += 1
        self.duplication[chrom.name] = True
        chrom.replication_ori_bound = False

        return new_chrom


        pass

if __name__ == '__main__':
    rep = Replication('replication', 'replication')
    print(rep.helicase)
    print (rep.polymerase)
    chrom1 = Chromosome('chrom1', ['A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A''A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A'])
    print (len(chrom1.sequence))
    rep.initiate(chrom1)
    rep.elongate(chrom1)
