class Process(object): #master-process definition
    """
    Parent for all cellular processes.
    """

    def __init__(self, id, name):
        self.id = id
        self.name = name

        self.enzyme_ids = []	#data input for process
        self.substrate_ids = []

    def set_states(self, substrate_ids, enzyme_ids):		#set-method for process-input
        self.enzyme_ids = enzyme_ids
        self.substrate_ids = substrate_ids

    def update(self, model):
        """
        Has to be implemented by child class.
        """
        pass
