from plasmid import Part, Module, Plasmid
from utils import get_module, get_scars
from sbol_serialisation import serialise_sbol
from sbol2 import Document, setHomespace
from utils import get_module, get_scars
setHomespace('http://project7.com')
doc = Document()
doc.clear()


#This example has an input sequence and uses this to create modules/plasmid
seq = "TGAACGTGATACCACCATGCCGGCAGCAATGGCGACCACCCTGCGTAAACTGCTGACGGGTGAGCTGCTGACCCTGGCAAGCCGCCAGCAACTGATTGATTGGATGGAAGCGGATAAAGTGGCGGGTCCGCTGCTGCGTAGCGCGCTGCCGGCTGGCTGGTTTATTGCGGATAAAAGCGGTGCGGGCGAACGTGGCAGCCGTGGCATTATTGCGGCGCTGGGCCCGGATGGTAAACCGAGCCGTATTGTGGTGATTTATACCACCGGCAGCCAGGCGACGATGGATGAACGTAACCGTCAGATTGCGGAAATTGGCGCGAGCCTGATTAAACATTGGTAAACCGATACAATTAAAGGCTCCTTTTGGAGCCTTTTTTTTTGGACGACCCTTGTCCTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTAGGCCGGCCGATCTGAAGATCAGCAGTTCAACCTGTTGATAGTACGTACTAAGCTCTCATGTTTCACGTACTAAGCTCTCATGTTTAACGTACTAAGCTCTCATGTTTAACGAACTAAACCCTCATGGCTAACGTACTAAGCTCTCATGGCTAACGTACTAAGCTCTCATGTTTCACGTACTAAGCTCTCATGTTTGAACAATAAAATTAATATAAATCAGCAACTTAAATAGCCTCTAAGGTTTTAAGTTTTATAAGAAAAAAAAGAATATATAAGGCTTTTAAAGCCTTTAAGGTTTAACGGTTGTGGACAACAAGCCAGGGATGTAACGCACTGAGAAGCCCTTAGAGCCTCTCAAAGCAATTTTGAGTGACACAGGAACACTTAACGGCTGACATGGGGCGCGCCCAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCTTTAATTAAAGCGGATAACAATTTCACACAGGAGGCCGCCTAGGCCGCGGCCGCGCGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGCGGCCGCGTCGTGACTGGGAAAACCCTGGCGACTAGTCTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAGGGGTCCCCAATAATTACGATTTAAATTAGTAGCCCGCCTAATGAGCGGGCTTTTTTTTAATTCCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGCATTCAGCATTTTCGTGTGGCGCTGATTCCGTTTTTTGCGGCGTTTTGCCTGCCGGTGTTTGCGCATCCGGAAACCCTGGTGAAAGTGAAAGATGCGGAAGATCAACTGGGTGCGCGCGTGGGCTATATTGAACTGGATCTGAACAGCGGCAAAATTCTGGAATCTTTTCGTCCGGAAGAACGTTTTCCGATGATGAGCACCTTTAAAGTGCTGCTGTGCGGTGCGGTTCTGAGCCGTGTGGATGCGGGCCAGGAACAACTGGGCCGTCGTATTCATTATAGCCAGAACGATCTGGTGGAATATAGCCCGGTGACCGAAAAACATCTGACCGATGGCATGACCGTGCGTGAACTGTGCAGCGCGGCGATTACCATGAGCGATAACACCGCGGCGAACCTGCTGCTGACGACCATTGGCGGTCCGAAAGAACTGACCGCGTTTCTGCATAACATGGGCGATCATGTGACCCGTCTGGATCGTTGGGAACCGGAACTGAACGAAGCGATTCCGAACGA"
fetched_cargo = get_module(seq, "fetched_cargo", "cargo","cargo")
fetched_marker = get_module(seq, "fetched_marker", "marker","marker")
fetched_origin = get_module(seq, "fetched_origin", "origin","origin")
check_scars = get_scars(seq)
scar = Part("scar1","scar111","CAATAATTACG","scar")
plasmid12 = Plasmid("plasmid2","plasmid222")
plasmid12.structure = [fetched_cargo,fetched_marker,fetched_origin]
plasmid12.insert_scar(scar,"T0","after")
print(plasmid12.get_sequence())


#Use this to serialise in genbank and create plasmid you can visualise in Benchling etc
#plasmid12.serialise_genbank("path")


# #How you would make it from scratch
# promoter2 = Part("promoter", "1", "GGG", "promoter")
# cds2 = Part("cds", "2", "TTT", "gene")
# terminator2 = Part("terminator", "3", "AAA", "terminator")
# marker_part = Part("marker_part","4","GGG","abr")
# origin_part = Part("origin_part","5","CCC","ori")

# cargo = Module("cargo_name", "6", "cargo")
# marker = Module("marker_name", "7", "marker", structure = [marker_part])
# origin = Module("origin_name", "8", "origin", structure = [origin_part])

# cargo_list = [promoter2, cds2, terminator2]
# cargo.structure = cargo_list

# Plasmid1 = Plasmid("plasmid", "9")
# plasmid_list = [cargo, marker, origin]
# Plasmid1.structure = plasmid_list
# plasmid_sequence = Plasmid1.get_sequence()
# print(plasmid_sequence)

# #Serialise in SBOL
# #serialise_sbol(plasmid12,"file_path")