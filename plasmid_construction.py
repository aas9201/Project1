from plasmid import Part, Module, Plasmid
from utils import get_module, get_scars
from sbol_serialisation import serialise_sbol
from sbol2 import Document, setHomespace
from utils import get_module, get_scars


seq = "TGAACGTGATACCACCATGCCGGCAGCAATGGCGACCACCCTGCGTAAACTGCTGACGGGTGAGCTGCTGACCCTGGCAAGCCGCCAGCAACTGATTGATTGGATGGAAGCGGATAAAGTGGCGGGTCCGCTGCTGCGTAGCGCGCTGCCGGCTGGCTGGTTTATTGCGGATAAAAGCGGTGCGGGCGAACGTGGCAGCCGTGGCATTATTGCGGCGCTGGGCCCGGATGGTAAACCGAGCCGTATTGTGGTGATTTATACCACCGGCAGCCAGGCGACGATGGATGAACGTAACCGTCAGATTGCGGAAATTGGCGCGAGCCTGATTAAACATTGGTAAACCGATACAATTAAAGGCTCCTTTTGGAGCCTTTTTTTTTGGACGACCCTTGTCCTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTAGGCCGGCCGATCTGAAGATCAGCAGTTCAACCTGTTGATAGTACGTACTAAGCTCTCATGTTTCACGTACTAAGCTCTCATGTTTAACGTACTAAGCTCTCATGTTTAACGAACTAAACCCTCATGGCTAACGTACTAAGCTCTCATGGCTAACGTACTAAGCTCTCATGTTTCACGTACTAAGCTCTCATGTTTGAACAATAAAATTAATATAAATCAGCAACTTAAATAGCCTCTAAGGTTTTAAGTTTTATAAGAAAAAAAAGAATATATAAGGCTTTTAAAGCCTTTAAGGTTTAACGGTTGTGGACAACAAGCCAGGGATGTAACGCACTGAGAAGCCCTTAGAGCCTCTCAAAGCAATTTTGAGTGACACAGGAACACTTAACGGCTGACATGGGGCGCGCCCAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCTTTAATTAAAGCGGATAACAATTTCACACAGGAGGCCGCCTAGGCCGCGGCCGCGCGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGCGGCCGCGTCGTGACTGGGAAAACCCTGGCGACTAGTCTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAGGGGTCCCCAATAATTACGATTTAAATTAGTAGCCCGCCTAATGAGCGGGCTTTTTTTTAATTCCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGCATTCAGCATTTTCGTGTGGCGCTGATTCCGTTTTTTGCGGCGTTTTGCCTGCCGGTGTTTGCGCATCCGGAAACCCTGGTGAAAGTGAAAGATGCGGAAGATCAACTGGGTGCGCGCGTGGGCTATATTGAACTGGATCTGAACAGCGGCAAAATTCTGGAATCTTTTCGTCCGGAAGAACGTTTTCCGATGATGAGCACCTTTAAAGTGCTGCTGTGCGGTGCGGTTCTGAGCCGTGTGGATGCGGGCCAGGAACAACTGGGCCGTCGTATTCATTATAGCCAGAACGATCTGGTGGAATATAGCCCGGTGACCGAAAAACATCTGACCGATGGCATGACCGTGCGTGAACTGTGCAGCGCGGCGATTACCATGAGCGATAACACCGCGGCGAACCTGCTGCTGACGACCATTGGCGGTCCGAAAGAACTGACCGCGTTTCTGCATAACATGGGCGATCATGTGACCCGTCTGGATCGTTGGGAACCGGAACTGAACGAAGCGATTCCGAACGA"
fetched_cargo = get_module(seq, "fetched_cargo", "cargo","cargo")
fetched_marker = get_module(seq, "fetched_marker", "marker","marker")
fetched_origin = get_module(seq, "fetched_origin", "origin","origin")
check_scars = get_scars(seq)
scar = Part("scar1","scar111","CAATAATTACG","scar")
plasmid2 = Plasmid("plasmid2","plasmid222")
plasmid2.structure = [fetched_cargo,fetched_marker,fetched_origin]
plasmid2.insert_scar(scar,"T0","after")





#How you would make it from scratch
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



def linear_construction(plasmid1,plasmid2,desired_plasmid):
    if not isinstance(desired_plasmid, list):
            raise TypeError("Desired Plasmid must be of Type: list")
    for modules in desired_plasmid:
        if not isinstance(modules, Module):
            raise TypeError("List should contain Module types only")

        # Creates a list containing the module types and compares it with predefined structure to see if its correct
    correct_val = ["cargo", "marker", "origin"]
    order = []
    for module in desired_plasmid:
        order.append(module.module_type)
    if correct_val != order:
        raise ValueError(
            "Please ensure the modules are in the correct order: [Cargo,Marker,Origin]"
        )      
    
    desired_list = []
    for module in desired_plasmid:
         desired_list.append(module.unique_id)
    
    plasmid1_list = []
    for module in plasmid1.components:
        plasmid1_list.append(module.unique_id)

    plasmid2_list = []
    for module in plasmid2.components:
        plasmid2_list.append(module.unique_id)

    desired_set = set(desired_list)
    plasmid1_set = set(plasmid1_list)
    plasmid2_set = set(plasmid2_list)
    
    contribution1 = len(desired_set.intersection(plasmid1_set))
    if contribution1 == (0 or 3):
         raise ValueError("Linear Construction is not possible")
    
    contribution2 = len(desired_set.intersection(plasmid2_set))
    if contribution2 == (0 or 3):
         raise ValueError("Linear Construction is not possible")
    
    if contribution1 > contribution2:
         parent_plasmid = plasmid1_list
         if parent_plasmid[0] != desired_list[0]:
              print(f"Cargo module from {plasmid2.name} has been inserted into the parent plasmid: {plasmid1.name}")
         if parent_plasmid[1] != desired_list[1]:
              print(f"Marker module from {plasmid2.name} has been inserted into the parent plasmid: {plasmid1.name}")
         if parent_plasmid[2] != desired_list[2]:
              print(f"Origin module from {plasmid2.name} has been inserted into the parent plasmid: {plasmid1.name}")
    
    if contribution1 < contribution2:
         parent_plasmid = plasmid2_list
         if parent_plasmid[0] != desired_list[0]:
              print(f"Cargo module from {plasmid1.name} has been inserted into the parent plasmid: {plasmid2.name}")
         if parent_plasmid[1] != desired_list[1]:
              print(f"Marker module from {plasmid1.name} has been inserted into the parent plasmid: {plasmid2.name}")
         if parent_plasmid[2] != desired_list[2]:
              print(f"Origin module from {plasmid1.name} has been inserted into the parent plasmid: {plasmid2.name}")



# top_list = ['6', 'marker', 'origin']
# list1 = ['6', '7', '8']
# list2 = ['cargo', 'marker', 'origin']

# # Convert lists to sets for intersection operation
# top_set = set(top_list)
# set1 = set(list1)
# set2 = set(list2)

# # Calculate intersection sizes
# intersection_size1 = len(top_set.intersection(set1))
# intersection_size2 = len(top_set.intersection(set2))

#     print(desired_list)
#     print(plasmid1_list)
#     print(plasmid2_list)



linear_construction(plasmid1,plasmid2,[cargo,fetched_marker,fetched_origin])










