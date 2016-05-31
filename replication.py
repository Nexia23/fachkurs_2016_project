from molecules import BioMoleculeCount
from molecules import BioMolecule
from processes import Process

# generell:
# start replication?
# replication fork
# origins of Replication? ARS , Origin recognition complex bindet an A-rich DNA
# Wo bindet Helikase, wo Polymerase?
# ATP Verbrauch?
# mehrere Chromosome gleichzeitig? 
# Gegenstrang 
# vllt noch Reparaturmodul für Mutation
# koordination mit Transkription

class Chromosome:
    def __init__(self, name, sequence):
        self.sequence = sequence
        self.name = name
        self.replication_ori_bound = False # beginnt ungebunden, 
        # da in der Replikation erstmal die Helikase und die Polymerase binden müssen 
        # > Initiationsbedingung


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
        self.chromosomes = {} # dict für das neue Chromosom
        self.duplication = {}
		# def. old_chromosomes > states
        #> zeigt, welche Chromosome schon doppelt vorliegen
        #ruft hel+PM auf
        #self._init_helicase()
        #self._init_polymerase()
    
    # hel+PM werden definiert

    
    def update(self, model):
        #call state (> transcription?)
        # transcription blockiert die Stelle für replication
        #for #key in #
        if self.duplication == {}:
	        self.duplication = {chrom: None for chrom in self.substrate_ids} # dict für schon replizierte Chromosome 

        self.chromosome_names = self.substrate_ids #Chromosom wird aufgerufen
        
        self.polymerase =  model.states['Polymerase3']
        self.helicase = model.states['DnaB']
        self.old_chromosomes = [model.states[chromsome_name] for chromsome_name in self.chromosome_names]


        for i, old_chromosome in enumerate(self.old_chromosomes): # wenn das aufgerufene Chromosom folgende Bedingungen erfüllt wird es zur entprechenden Phase weitergeleitet
            if not old_chromosome.replication_ori_bound and not self.duplication[self.chromosome_names[i]] and self.polymerase.count > 0 and self.helicase.count > 0:
            #and not transcription an aktueller Stelle:
                self.initiate(old_chromosome) 
            elif old_chromosome.replication_ori_bound: #and not transcription an der Stelle:
                self.elongate(old_chromosome)
            elif (self.polymerase.count == 0 or self.helicase.count == 0) and not old_chromosome.replication_ori_bound and not self.duplication[chrom.name]:
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
            # Chromosom wird aufgerufen, ein neues mit gleichem Namen aber noch leerer Sequenz wird erstellt:
        self.chromosomes[chrom.name]=Chromosome(chrom.name,[])# legt das dict 'Chromosomes' an, Name bleibt erhalten, Sequenz ändert sich
                # wenn die sequenz repliziert wird haben wir einen Gegenstrang, müsste der Formhalber nicht wieder davon der Gegenstrang genommen werden, 
                # damit wir eine konsistente Speicherung der Stänge haben?
        self.polymerase.count -= 1 #jew count -1 für jedes gebundene molekül
        self.helicase.count -= 1
        chrom.replication_ori_bound = True #setzt das Chromosom auf "gebunden"
        #> damit geht es im nächsten timestep einfach direkt bei der Elongation weiter
        self.duplication[chrom.name] = False # Chromosom ist noch nicht dupliziert
  
    def elongate(self, old_chrom: Chromosome):
        print('elongation')
        #mehrere Schritte (Zeitschritte) >listeneinträge new_chrom  >>vgl alte sequ
            # dictionary 'chromosomes' enthält 'Chromosomname':'Eigenschaften(Name, Sequenz')

        
    
        #definition i?
        new_chrom = self.chromosomes[old_chrom.name]
        #import von initiation als old_chrom, hier wird das New_chrom erstellt >neue sequenz
        
        # replikation erfolgt mit 100 bp/s also in jedem timestep werden nur 100 nukleotide der Sequnenz gelesen und repliziert, 
        # wenn die restlichen nukleotide weniger als 100 sind läuft die PM einfach bis zum Schluss der Sequenz, 
        #leitet die Termination ein und returned die neue Sequence

        if len(new_chrom.sequence) == 0:
            sequence_to_replicate = old_chrom.sequence[0:100]
            # Bindung checken für 100-200 usw.
        elif len(old_chrom.sequence) - len(new_chrom.sequence) > 100:
            sequence_to_replicate = old_chrom.sequence[len(new_chrom.sequence)+1:len(new_chrom.sequence)+100]
        elif len(old_chrom.sequence) - len(new_chrom.sequence) == 0:
            return self.terminate(new_chrom.sequence)
        else:
            sequence_to_replicate = old_chrom.sequence[len(new_chrom.sequence) + 1:]
         # transcription berücksichtigen

         # hier kommt noch die Mutation hin: (bisher nur reiner Komplementärstrang)     
      

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
        self.duplication[chrom.name] = True# Eintrag als dupliziertes Chromosom
        chrom.replication_ori_bound = False# wieder ungebunden 

        return new_chrom
         #! jedes Chromosom sollte nur einmal dupliziert werden!


        pass

if __name__ == '__main__':
    rep = Replication('replication', 'replication')
    print(rep.helicase)
    print (rep.polymerase)
    chrom1 = Chromosome('chrom1', ['A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A''A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A'])
    #chrom1 = chr_list[0]
    #chrom2 = chr_list[1]
    print (len(chrom1.sequence))

    rep.initiate(chrom1)
    rep.elongate(chrom1)
   

