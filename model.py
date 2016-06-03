import modeldata
import molecules as mol
import translation
import replication as rep
import transcription

class Output:
    """
    class for handling the simulation results of the different species types
    """

    def __init__(self, model):
        self.meta = {}
        self.model = model
        self.timecourses = {state: SimulationResult(model.states[state]) for state in model.states}

    def add_timepoint(self, species):
        """
        add a simulation time point for one species
        @param species: mol.BioMolecule
        @return: None
        """
        if isinstance(self.model.states[species], mol.Polymer):
            pass  # TODO: implement a useful method for Polymers
        elif isinstance(self.model.states[species], mol.BioMoleculeCount):
            self.timecourses[species].add_timepoint(self.model.states[species].count, self.model.timestep)


class SimulationResult:
    """
    handles and stores a simulation result for one species
    """

    def __init__(self, species):
        self.name = species.name
        self.value = []
        self.time = []

    def add_timepoint(self, time, value):
        self.value.append(value)
        self.time.append(time)


class Model:
    """
    Initializes the states and processes for the model and lets the processes update their corresponding states.
    """

    def __init__(self):
        self.states = {}        #dictionary with all molecules {Rib_name: Rib_object, mrna_ids: mrna_object, mrna2_id: ...}
        self.processes = {}     #dictionary filled with all active processes
        self.timestep = 0
        self.mrnas = {}  # all selfs should be initialized in the constructor
        self.ribosomes = {} #dictionary will be filled with 10 Ribosomes
        self.helicases = {}
        self.polymerases = {}
        self.chromosomes = {}
        self.volume = 1
        self.db = modeldata.ModelData()
        self.chromosomes=modeldata.ModelData.createchromosomes()    #list with chromosomes
        self.genes=modeldata.ModelData.creategenes()                #dictionary with genes

        self.__initialize_macromolecules()
        self.__initialize_states()
        self.__initialize_processes()
        #self.results = Output(self)  


    def __initialize_macromolecules(self):
        self.ribosomes = {'Ribosomes': mol.Ribosome('Ribos', 'Ribosomes', 187000)}
        self.polymerase2= mol.RNAPolymeraseII('Pol2', 'Polymerase2', 100000000)
        self.nucleotides= mol.NucleotidPool('Nucs','Nucleotides', 100000)
        self.helicases = {'DnaB': rep.Helicase("Helicase", "DnaB", 100)}
        self.polymerases = {'Polymerase3' :rep.Polymerase("Polymerase", "Polymerase3", 100)}
        self.chromosomes = {x.id:x for x in modeldata.ModelData.createchromosomes()}
        
        #for i, mrna in enumerate(self.db.get_states(mol.MRNA)):
         #   mid, name, sequence = mrna
          #  sequence=list(sequence)
           # sequence[0:3]='AUG'
            #sequence[12:15]='UGA'
            #sequence=''.join(sequence)
            #self.mrnas[mid] = [mol.MRNA(mid, name, sequence)]
            
        
    def __initialize_states(self):
        """
        initialize the different states
        """
        self.states.update(self.ribosomes)  #adding dictionaries to self.states
        self.states.update(self.helicases)
        self.states.update(self.polymerases)
        self.states.update(self.chromosomes)
        self.states.update(self.mrnas)
        self.states["Nucleotides"] = self.nucleotides


    def __initialize_processes(self):
        """
        initialize processes
        """
        
        # transcription
        trsc = transcription.Transcription(0, 'Transcription')
        trsc.set_states(self.genes.keys(), self.polymerase2)
        self.processes["Transkription"] = trsc

        # translation
        trsl = translation.Translation(1, "Translation")
        trsl.set_states(self.mrnas.keys(), self.ribosomes.keys())           #states in Process are keys: Rib_name, mrna_name?!
        self.processes["Translation"] = trsl

        # replication
        #repl =rep.Replication(2, "Replication")
        #replication_enzyme_ids= list(self.helicases.keys()).extend(list(self.polymerases.keys()))
        #repl.set_states(list(self.chromosomes.keys()), replication_enzyme_ids)
        #self.processes.update({"Replication":repl})


    def step(self):
        """
        Do one update step for each process.

        """
        for p in self.processes:
            self.processes[p].update(self)

        #for state in self.states:
         #   self.results.add_timepoint(state)

        self.timestep += 1

    def simulate(self, steps, log=True):
        """
        Simulate the model for some time.
        """
        for s in range(steps):
            self.step()
            if log:  # This could be an entry point for further logging
                #print all states
                print('\r{}'.format([self.states[x] for x in self.states.keys() ], end=''))


if __name__ == "__main__":
    c = Model()
    c.simulate(300, log=True)
