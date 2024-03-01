from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


class Part:
    """Class represents the most fundamental componenet of a plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Sequence: String (only ACTG and no spaces)
    Role: String, gives function for further info and serialiastion
    Description: String, add some information about the part"""

    def __init__(self, name, unique_id, sequence, role, description=""):
        self.name = name
        self.unique_id = unique_id
        self.sequence = sequence
        self.role = role
        self.description = description

    @property
    def name(self):
        """gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """gets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """gets the name attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def sequence(self):
        """gets the sequence attribute"""
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """sets the sequence attribute, type: String, Only contains A,C,T or G with no spaces"""
        if not isinstance(sequence, str):
            raise TypeError("Name must be a string")
        if not all(c in "ACTG" for c in sequence) or " " in sequence:
            raise ValueError("Sequence must only contain A,C,T,G with no spaces")
        self._sequence = sequence

    @property
    def role(self):
        """gets the role attribute"""
        return self._role

    @role.setter
    def role(self, role):
        """sets the role attribute, type: String, must only be certain roles also:
            "gene",
            "abr",
            "ori",
            "terminator",
            "promoter",
            "scar",
            "re",
            "oriT",
            "cargo","""
        

        if not isinstance(role, str):
            raise TypeError("Role must be a string")
        if role not in [
            "gene",
            "abr",
            "ori",
            "terminator",
            "promoter",
            "scar",
            "re",
            "oriT",
            "cargo",
        ]:
            raise ValueError(
                "Role must be either gene, abr, ori, terminator, promoter, scar, re, oriT, cargo"
            )

        self._role = role

    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description


class Module:
    """Class represents the most module componenet of a plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Module Type: String (either Cargo, Marker, Origin or Invariant)
    Structure: List of Part objects
    Description: String, add some information about the part
    Fetched: True or False, allows for different structure construction if sequence is imported
    """


    def __init__(
        self, name, unique_id, module_type, description="", structure=[], fetched=False
    ):
        self.name = name
        self.unique_id = unique_id
        self.module_type = module_type
        self.description = description
        self.fetched = fetched
        self.structure = structure

    @property
    def name(self):
        """gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """gets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """sets the unique_id attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id 
        
      
        
    @property
    def module_type(self):
        """gets the module_type attribute"""
        return self._module_type

    @module_type.setter
    def module_type(self, module):
        """sets the unique_id attribute"""
        if module not in ["cargo", "marker", "origin", "invariant"]:
            raise ValueError("Module must be either a cargo, marker, origin, invariant")
        self._module_type = module


    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description

    @property
    def structure(self):
        """gets the structure attribute"""
        return self._structure

    @structure.setter
    def structure(self, structure_list):
        """Function which sets the structure
        Input: list of part(s) objects
        Depending on module_type a predefined structure with corresponding flanking RE sites are assigned
        """

        if not isinstance(structure_list, list):
            raise TypeError("Input structure must be of Type: list")
        for parts in structure_list:
            if not isinstance(parts, Part):
                raise TypeError("List should contain Part types only")

        if (
            self.module_type == "cargo"
        ):  # Creates Part objects for flanking RE sites dynamically
            self.cargo = None
            self.PacI = Part("PacI", "paci", "TTAATTAA", "re")
            self.SpeI = Part("SpeI", "spei", "ACTAGT", "re")
            self._structure = [self.PacI, self.cargo, self.SpeI]
            # The input part list is inserted by a slice to remove the placeholder
            self._structure[1:2] = (
                structure_list  
            )

        if self.module_type == "origin":
            self.replication = None
            self.FseI = Part("FseI", "fsei", "GGCCGGCC", "re")
            self.AscI = Part("AscI", "asci", "GGCGCGCC", "re")
            self._structure = [self.FseI, self.replication, self.AscI]
            # The input part list is inserted by a slice to remove the placeholder
            self._structure[1:2] = (
                structure_list  
            )

        if self.module_type == "marker" and self.fetched == False:
            self.marker = None
            self.PshAI = Part(
                "pshai", "1", "GACGTC", "re"
            )  #Need to sort this logic out a bit more, however SEVA synbiohub is down for some reason so I can't access the sequence logic requried
            self.SwaI = Part("swai", "1", "ATTTAAAT", "re")
            self._structure = [self.SwaI, self.marker, self.PshAI]
            # The input part list is inserted by a slice to remove the placeholder
            self._structure[1:2] = (
                structure_list  
            )

        if self.module_type == "marker" and self.fetched == True:
            self._structure = structure_list  # If sequence is imported I want to assign PshAI dynamically as certain SEVA plasmids contain different sequence

        if self.module_type == "invariant":
            self._structure = structure_list




    def get_sequence(self):
        """Obtains a sequence by concatenating sequence of parts together"""
        sequence = ""
        for i in self.structure:
            sequence += i.sequence
        return sequence


class Plasmid:
    """Class which builds final plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Module Type: String (either Cargo, Marker, Origin or Invariant)
    Structure: List of Module objects
    Description: String, add some information about the part
    """

    def __init__(self, name, unique_id, description=""):
        self.name = name
        self.unique_id = unique_id
        self.description = description
        self.components = []

        # Three placeholder variables to just provide template structure
        self.cargo_module = None
        self.marker_module = None  
        self.origin_module = None

        # Creating the non-variable modules on construction of an object
        self.t1_part = Part(
            "t1_part",
            "t1_part",
            "CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT",
            "terminator",
        )
        self.t1 = Module("t1", "t1_module", "invariant", structure=[self.t1_part])

        self.t0_part = Part(
            "t0_part",
            "t0_part",
            "CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG",
            "terminator",
        )
        self.SanDI = Part("SanDI_part", "sandi", "GGGTCCC", "re")
        self.t0 = Module(
            "t0", "t0_module", "invariant", structure=[self.t0_part, self.SanDI]
        )

        self.orit_part = Part(
            "orit_part",
            "orit_part",
            "CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA",
            "oriT",
        )
        self.orit = Module(
            "orit", "orit_module", "invariant", structure=[self.orit_part]
        )

        # Skip first validation on initilisation so ._
        self._structure = [
            self.cargo_module,
            self.t0,
            self.marker_module,
            self.orit,
            self.origin_module,
            self.t1,
        ]

    @property
    def name(self):
        """Gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """Sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """Sets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """Sets the unique_id attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description


    @property
    def structure(self):
        """Gets the structure attribute"""
        return self._structure

    @structure.setter
    def structure(self, structure):
        """Sets the name attribute
        Input: A list of three individual seperate modules, type Cargo,Marker,Origin"""
        if not isinstance(structure, list):
            raise TypeError("Input structure must be of Type: list")
        for modules in structure:
            if not isinstance(modules, Module):
                raise TypeError("List should contain Module types only")

        # Creates a list containing the module types and compares it with predefined structure to see if its correct
        correct_val = ["cargo", "marker", "origin"]
        order = [module.module_type for module in structure]

        if correct_val != order:
            raise ValueError(
                "Please ensure the modules are in the correct order: [Cargo,Marker,Origin]"
            )
        
        #For use in the future construction methods
        self.components = structure 

        iterVal = iter(structure)
        for parts, element in enumerate(self._structure):
            if (
                element is None
            ):  # Finds which variables are "None" in structure AKA placeholders and then swaps them for each corresponding module
                self._structure[parts] = next(iterVal)

    def get_sequence(
        self,
    ):  # it iterates through each module in structure and then calls their get_sequence function and concatenates all together
        plasmid_sequence = ""
        for module in self.structure:
            plasmid_sequence += module.get_sequence()
        return plasmid_sequence
    

    def insert_scar(self, scar, module=None, position=None):
        """ "Method to insert scars:
        Modules: T1,T0,OriT
        Position: before or after"""
        module_dict = {"T1": self.t1, "T0": self.t0, "OriT": self.orit}

        module = module_dict.get(
            module
        )  # Get the module object based on the module name

        if (
            module
        ):  # Inserts the scar depending on if it comes before or after the module
            if position == "before":
                module.structure.insert(0, scar)
            elif position == "after":
                module.structure.append(scar)

    def gather_data(self):
        """Method which gathers the data of each part in the plasmid and puts it all into a dictionary with key being name and values being other data attributes and locations"""
        sequence = self.get_sequence()
        positions = []
        end_point = []
        my_dict = {}

        for module in self.structure:
            for part in module.structure:
                position = sequence.find(
                    part.sequence
                )  # Finds where the part occurs in the sequence and adds it to dictionary along with other data attributes
                positions.append(position)
                my_dict[part.name] = [position, part.unique_id, part.role, part.description]


        for i in range(
            len(positions) - 1
        ):  # Looks at position of where next part is and then -1 as that will be the endpoint
            end_point.append(positions[i + 1] - 1)
        end_point.append(len(sequence))

        for part, (key, value) in enumerate(my_dict.items()):
            # Update each key in my_dict with its corresponding end point
            if part < len(end_point):
                my_dict[key].append(end_point[part])

        return my_dict

    def serialise_genbank(self, path):
        """Serialises in genbank format: Need input file path"""
        sequence = self.get_sequence()
        sequence_r = SeqRecord(
            Seq(sequence),
            id=self.unique_id,
            name=self.name,
            description=self.description,
        )
        print(sequence_r)
        data_dict = (
            self.gather_data()
        )  # Obtains data dictionary and uses this to add features of each part

        conversion = {
            "re": "misc_feature",  # Dictionary to convert the role into the GenBank accepted role
            "cargo": "misc_feature",
            "abr": "misc_feature",
            "ori": "rep_origin",
        }

        for part, attribute in data_dict.items():  # Role Conversion
            if attribute[2] in conversion:
                attribute[2] = conversion[attribute[2]]

        sequence_r.annotations["molecule_type"] = (
            "DNA"  # Annotate the sequence to show that its DNA
        )

        # Going through the data dictionary an using it to annotate each part of the sequence
        for part, attribute in data_dict.items():
            feature = SeqFeature(
                FeatureLocation(attribute[0], attribute[4]),
                type=attribute[2],
                qualifiers={
                    "gene": [part],  # Gene name
                    "locus_tag": [attribute[1]],  # A unique identifier for the gene
                    "note": [attribute[3]],  # Any additional notes
                },
            )
            sequence_r.features.append(feature)

        # Write the SeqRecord to a GenBank file
        with open(path, "w") as output_handle:
            SeqIO.write(sequence_r, output_handle, "genbank")

        print(f"GenBank file saved at {path}")
        print(sequence_r)  # Help show that the serialisation has worked
