from plasmid import Part, Module, Plasmid
from utils import get_module, get_scars


class Construction():

    def __init__(self):
        self.construction_dictionary = {}

    

    def linear_construction(self,plasmid1, plasmid2, desired_plasmid_name, desired_plasmid):
        """A function which models the construction of a desired plasmid from two input plasmids
        The two input plasmids are of Type(Plasmid) and desired plasmid should be a list of three modules
        which are from the two input plasmids"""

        #Input validation
        if not isinstance(desired_plasmid, list):
            raise TypeError("Desired Plasmid must be of Type: list")
        for modules in desired_plasmid:
            if not isinstance(modules, Module):
                raise TypeError("List should contain Module types only")
            
        if not isinstance(desired_plasmid_name, str):
            raise TypeError("Desired Plasmid name must be of Type List")

            # Creates a list containing the module types and compares it with predefined structure to see if its correct
        correct_val = ["cargo", "marker", "origin"]
        order = []
        for module in desired_plasmid:
            order.append(module.module_type)
        if correct_val != order:
            raise ValueError(
                "Please ensure the modules are in the correct order: [Cargo,Marker,Origin]"
            )

            # Check if the desired plasmid name already exists in the dictionary
        if desired_plasmid_name not in self.construction_dictionary:
            self.construction_dictionary[desired_plasmid_name] = []
        else:
            raise ValueError("Desired Plasmid has already been constructed")


        #Creates list of module unique_ids from each plasmid for identification
        desired_list = [module.unique_id for module in desired_plasmid]

        plasmid1_list = [module.unique_id for module in plasmid1.components]
        
        plasmid2_list = [module.unique_id for module in plasmid2.components]
    
        #Set creation to model which plasmid is going to be parent/donor
        desired_set = set(desired_list)
        plasmid1_set = set(plasmid1_list)
        plasmid2_set = set(plasmid2_list)

        #As only inserting one module into the plasmid, the desired plasmid must contain either one or two modules
        contribution1 = len(desired_set.intersection(plasmid1_set))
        contribution2 = len(desired_set.intersection(plasmid2_set))
        if contribution1 in (0,3) or contribution2 in (0,3):
            raise ValueError("Linear Construction is not possible")

        if contribution1 > contribution2:
            parent_plasmid = plasmid1_list
            parent_name = plasmid1.name
            donor_name = plasmid2.name
        else:
            parent_plasmid = plasmid2_list
            parent_name = plasmid2.name
            donor_name = plasmid1.name

        modules = ["Cargo", "Marker", "Origin"]

        # Check which modules need to be updated and add messages to the dictionary
        for i, module in enumerate(modules):
            if parent_plasmid[i] != desired_list[i]:
                print(f"{module} module from {donor_name} has been inserted into the parent plasmid: {parent_name}")
                self.construction_dictionary[desired_plasmid_name].append((module, donor_name, parent_name))
                break 



seq = "TGAACGTGATACCACCATGCCGGCAGCAATGGCGACCACCCTGCGTAAACTGCTGACGGGTGAGCTGCTGACCCTGGCAAGCCGCCAGCAACTGATTGATTGGATGGAAGCGGATAAAGTGGCGGGTCCGCTGCTGCGTAGCGCGCTGCCGGCTGGCTGGTTTATTGCGGATAAAAGCGGTGCGGGCGAACGTGGCAGCCGTGGCATTATTGCGGCGCTGGGCCCGGATGGTAAACCGAGCCGTATTGTGGTGATTTATACCACCGGCAGCCAGGCGACGATGGATGAACGTAACCGTCAGATTGCGGAAATTGGCGCGAGCCTGATTAAACATTGGTAAACCGATACAATTAAAGGCTCCTTTTGGAGCCTTTTTTTTTGGACGACCCTTGTCCTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTAGGCCGGCCGATCTGAAGATCAGCAGTTCAACCTGTTGATAGTACGTACTAAGCTCTCATGTTTCACGTACTAAGCTCTCATGTTTAACGTACTAAGCTCTCATGTTTAACGAACTAAACCCTCATGGCTAACGTACTAAGCTCTCATGGCTAACGTACTAAGCTCTCATGTTTCACGTACTAAGCTCTCATGTTTGAACAATAAAATTAATATAAATCAGCAACTTAAATAGCCTCTAAGGTTTTAAGTTTTATAAGAAAAAAAAGAATATATAAGGCTTTTAAAGCCTTTAAGGTTTAACGGTTGTGGACAACAAGCCAGGGATGTAACGCACTGAGAAGCCCTTAGAGCCTCTCAAAGCAATTTTGAGTGACACAGGAACACTTAACGGCTGACATGGGGCGCGCCCAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCTTTAATTAAAGCGGATAACAATTTCACACAGGAGGCCGCCTAGGCCGCGGCCGCGCGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGCGGCCGCGTCGTGACTGGGAAAACCCTGGCGACTAGTCTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAGGGGTCCCCAATAATTACGATTTAAATTAGTAGCCCGCCTAATGAGCGGGCTTTTTTTTAATTCCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGCATTCAGCATTTTCGTGTGGCGCTGATTCCGTTTTTTGCGGCGTTTTGCCTGCCGGTGTTTGCGCATCCGGAAACCCTGGTGAAAGTGAAAGATGCGGAAGATCAACTGGGTGCGCGCGTGGGCTATATTGAACTGGATCTGAACAGCGGCAAAATTCTGGAATCTTTTCGTCCGGAAGAACGTTTTCCGATGATGAGCACCTTTAAAGTGCTGCTGTGCGGTGCGGTTCTGAGCCGTGTGGATGCGGGCCAGGAACAACTGGGCCGTCGTATTCATTATAGCCAGAACGATCTGGTGGAATATAGCCCGGTGACCGAAAAACATCTGACCGATGGCATGACCGTGCGTGAACTGTGCAGCGCGGCGATTACCATGAGCGATAACACCGCGGCGAACCTGCTGCTGACGACCATTGGCGGTCCGAAAGAACTGACCGCGTTTCTGCATAACATGGGCGATCATGTGACCCGTCTGGATCGTTGGGAACCGGAACTGAACGAAGCGATTCCGAACGA"
fetched_cargo = get_module(seq, "fetched_cargo", "cargo","cargo")
fetched_marker = get_module(seq, "fetched_marker", "marker","marker")
fetched_origin = get_module(seq, "fetched_origin", "origin","origin")
check_scars = get_scars(seq)
scar = Part("scar1","scar111","CAATAATTACG","scar") #Check the scars previously and this can be inserted manually.
plasmid2 = Plasmid("plasmid2","plasmid222")
plasmid2.structure = [fetched_cargo,fetched_marker,fetched_origin]
plasmid2.insert_scar(scar,"T0","after") #insert scar part object previously made



promoter2 = Part("promoter", "1", "GGG", "promoter")
cds2 = Part("cds", "2", "TTT", "gene")
terminator2 = Part("terminator", "3", "AAA", "terminator")
marker_part = Part("marker_part","4","GGG","abr")
origin_part = Part("origin_part","5","CCC","ori")

cargo = Module("cargo_name", "6", "cargo")
marker = Module("marker_name", "7", "marker", structure = [marker_part])
origin = Module("origin_name", "8", "origin", structure = [origin_part])

cargo_list = [promoter2, cds2, terminator2]
cargo.structure = cargo_list

plasmid1 = Plasmid("plasmid1", "9")
plasmid_list = [cargo, marker, origin]
plasmid1.structure = plasmid_list


hello = Construction()
hello.linear_construction(plasmid1,plasmid2,"New_Plasmid", [cargo, fetched_marker, fetched_origin])
hello.linear_construction(plasmid1,plasmid2,"New_Plasmid1", [fetched_cargo, marker, fetched_origin])
hello.linear_construction(plasmid1,plasmid2,"New_Plasmid2", [cargo, marker, fetched_origin])
print(hello.construction_dictionary)