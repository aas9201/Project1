from plasmid import Part, Module, Plasmid
import unittest

class TestPart(unittest.TestCase):
    test_part = Part("name", "id", "AAA", "gene", "description")

    def test_part_name_initialisation(self):
        self.assertEqual(TestPart.test_part.name, "name")

    def test_part_unique_id_initialisation(self):
        self.assertEqual(TestPart.test_part.unique_id, "id")

    def test_part_sequence_initialisation(self):
        self.assertEqual(TestPart.test_part.sequence, "AAA")

    def test_part_role_initialisation(self):
        self.assertEqual(TestPart.test_part.role, "gene")

    def test_part_description_initialisation(self):
        self.assertEqual(TestPart.test_part.description, "description")

    def test_part_instance(self):
        self.assertIsInstance(TestPart.test_part, Part)

    def test_name_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Part(1, "unique_id", "ACTG", "gene")

    def test_unique_id_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Part("name", 1, "ACTG", "gene")

    def test_sequence_incorrect_initialisation_type_error(self):
        with self.assertRaises(TypeError):
            Part("name", "unique_id", 1, "gene")

    def test_sequence_incorrect_initialisation_character_error(self):
        with self.assertRaises(ValueError):
            Part("name", "unique_id", "AXTG", "gene")

    def test_sequence_incorrect_initialisation_spaces_error(self):
        with self.assertRaises(ValueError):
            Part("name", "unique_id", "A CTG", "gene")

    def test_role_incorrect_initialisation_type_error(self):
        with self.assertRaises(TypeError):
            Part("name", "unique_id", "ACTG", None)

    def test_role_incorrect_initialisation_invalid_role(self):
        with self.assertRaises(ValueError):
            Part("name", "unique_id", "ACTG", "invalid_role")

    def test_description_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Part("name", "unique_id", "ACTG", "gene", 1)



class TestModule(unittest.TestCase):

    test_module = Module("name", "id", "cargo", "description", [], False)

    def test_module_name_initialisation(self):
        self.assertEqual(TestModule.test_module.name, "name")

    def test_module_unique_id_initialisation(self):
        self.assertEqual(TestModule.test_module.unique_id, "id")

    def test_module_type_initialisation(self):
        self.assertEqual(TestModule.test_module.module_type, "cargo")

    def test_module_description_initialisation(self):
        self.assertEqual(TestModule.test_module.description, "description")

    def test_module_fetched_initialisation(self):
        test_module = Module("name", "id", "cargo", "description", [], True)
        self.assertTrue(test_module.fetched)


    def test_module_cargo_structure_initialisation(self):
        part = Part("name1", "id1", "AAA", "gene")
        test_module = Module("name", "id", "cargo", "description", [part], False)
        self.assertEqual(test_module.structure, [test_module.PacI,part,test_module.SpeI])
        self.assertIsInstance(test_module.PacI, Part)
        self.assertIsInstance(test_module.SpeI, Part)
        self.assertEqual(test_module.PacI.sequence, "TTAATTAA")
        self.assertEqual(test_module.SpeI.sequence, "ACTAGT")


    def test_module_origin_structure_initialisation(self):
        part = Part("name1", "id1", "AAA", "gene")
        test_module = Module("name", "id", "origin", "description", [part], False)
        self.assertEqual(test_module.structure, [test_module.FseI,part,test_module.AscI])
        self.assertIsInstance(test_module.FseI, Part)
        self.assertIsInstance(test_module.AscI, Part)
        self.assertEqual(test_module.FseI.sequence, "GGCCGGCC")
        self.assertEqual(test_module.AscI.sequence, "GGCGCGCC")


    def test_module_marker_structure_initialisation(self):
        part = Part("name1", "id1", "AAA", "gene")
        test_module = Module("name", "id", "marker", "description", [part], False)
        self.assertEqual(test_module.structure, [test_module.SwaI,part,test_module.PshAI])
        self.assertIsInstance(test_module.PshAI, Part)
        self.assertIsInstance(test_module.SwaI, Part)
        self.assertEqual(test_module.PshAI.sequence, "GACGTC")
        self.assertEqual(test_module.SwaI.sequence, "ATTTAAAT")


    def test_fetched_module_marker_structure_initialisation(self):
        part = Part("name1", "id1", "AAA", "gene")
        test_module = Module("name", "id", "marker", "description", [part], True)
        self.assertEqual(test_module.structure, [part])


    def test_module_invariant_structure_initialisation(self):
        part = Part("name1", "id1", "AAA", "gene")
        test_module = Module("name", "id", "invariant", "description", [part], False)
        self.assertEqual(test_module.structure, [part])
  
  
    def test_name_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Module(1, "unique_id", "cargo", "description")

    def test_unique_id_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Module("name", 1, "cargo", "description")

    def test_module_type_incorrect_initialisation(self):
        with self.assertRaises(ValueError):
            Module("name", "unique_id", "invalid_type", "description")

    def test_description_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Module("name", "unique_id", "cargo", 1)


    def test_get_sequence(self):
        part = Part("name1", "id1", "AAA", "gene")
        test_module = Module("name", "id", "marker", "description", [part], False)
        self.assertEqual(test_module.get_sequence(),"ATTTAAATAAAGACGTC")



class TestPlasmid(unittest.TestCase):
    test_plasmid1 = Plasmid("plasmid_name", "plasmid_id1","plasmid_description")

    test_plasmid2 = Plasmid("plasmid_name", "plasmid_id2","plasmid_description")
    test_part1 = Part("name", "id1", "AAA", "gene", "description1")
    test_part2 = Part("name", "id2", "CCC", "gene", "description2")
    test_part3 = Part("name", "id3", "TTT", "gene", "description3")
    test_cargo_module = Module("name", "id4", "cargo", "description4", [test_part1], False)
    test_marker_module = Module("name", "id5", "marker", "description5", [test_part2], False)
    test_origin_module = Module("name", "id6", "origin", "description6", [test_part3], False)
    test_plasmid2.structure = [test_cargo_module,test_marker_module,test_origin_module]


    def test_plasmid_name_initialisation(self):
        self.assertEqual(TestPlasmid.test_plasmid1.name, "plasmid_name")

    def test_plasmid_unique_id_initialisation(self):
        self.assertEqual(TestPlasmid.test_plasmid1.unique_id, "plasmid_id1")

    def test_plasmid_description_initialisation(self):
        self.assertEqual(TestPlasmid.test_plasmid1.description, "plasmid_description")

    def test_plasmid_name_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Plasmid(1, "plasmid_id1","plasmid_description")
    
    def test_plasmid_unique_id_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Plasmid("plasmid_name", 1, "plasmid_description")

    def test_plasmid_description_incorrect_initialisation(self):
        with self.assertRaises(TypeError):
            Plasmid("plasmid_name", "plasmid_id1", 1)

    def test_plasmid_t1(self):
        self.assertIsInstance(TestPlasmid.test_plasmid1.t1_part, Part)
        self.assertIsInstance(TestPlasmid.test_plasmid1.t1, Module)

    def test_plasmid_t0(self):
        self.assertIsInstance(TestPlasmid.test_plasmid1.t0_part, Part)
        self.assertIsInstance(TestPlasmid.test_plasmid1.t0, Module)  
        self.assertIsInstance(TestPlasmid.test_plasmid1.SanDI, Part) 
        self.assertEqual(TestPlasmid.test_plasmid1.t0.structure, [TestPlasmid.test_plasmid1.t0_part,TestPlasmid.test_plasmid1.SanDI])
 
    def test_plasmid_orit(self):
        self.assertIsInstance(TestPlasmid.test_plasmid1.orit_part, Part)
        self.assertIsInstance(TestPlasmid.test_plasmid1.orit, Module)

    def test_plasmid_intial_structure(self):
        self.assertEqual(TestPlasmid.test_plasmid1.structure, [None, TestPlasmid.test_plasmid1.t0, None, TestPlasmid.test_plasmid1.orit, None, TestPlasmid.test_plasmid1.t1])


    def test_plasmid_set_list_structure(self):
        with self.assertRaises(TypeError):
            TestPlasmid.test_plasmid2.structure = TestPlasmid.test_cargo_module,TestPlasmid.test_marker_module,TestPlasmid.test_origin_module

    def test_plasmid_set_module_type_structure(self):
        with self.assertRaises(TypeError):
            TestPlasmid.test_plasmid2.structure = [1,TestPlasmid.test_marker_module,TestPlasmid.test_origin_module]

    def test_plasmid_set_module_list_structure(self):
        with self.assertRaises(ValueError):
            TestPlasmid.test_plasmid2.structure = [TestPlasmid.test_marker_module,TestPlasmid.test_cargo_module,TestPlasmid.test_origin_module]

    def test_plasmid_set_structure(self):
        self.assertEqual(TestPlasmid.test_plasmid2.structure, [TestPlasmid.test_cargo_module, TestPlasmid.test_plasmid2.t0, TestPlasmid.test_marker_module, TestPlasmid.test_plasmid2.orit, TestPlasmid.test_origin_module, TestPlasmid.test_plasmid2.t1])
    
    def test_plasmid_components(self):
        self.assertEqual(TestPlasmid.test_plasmid2.components,[TestPlasmid.test_cargo_module,TestPlasmid.test_marker_module,TestPlasmid.test_origin_module])


    def test_plasmid_sequence(self):
        self.assertEqual(TestPlasmid.test_plasmid2.get_sequence(),"TTAATTAAAAAACTAGTCTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAGGGGTCCCATTTAAATCCCGACGTCCTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTAGGCCGGCCTTTGGCGCGCCCAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT")

    def test_plasmid_insert_scar_before_T1(self):
        scar1 = Part("Scar","Scar1","AAA","scar")
        test_plasmid3 = Plasmid("plasmid3", "plamid3id", "description")
        test_plasmid3.insert_scar(scar1,"T1","before")
        self.assertEqual(test_plasmid3.t1.structure, [scar1, test_plasmid3.t1_part])

    def test_plasmid_insert_scar_after_T1(self):
        scar1 = Part("Scar","Scar1","AAA","scar")
        test_plasmid3 = Plasmid("plasmid3", "plamid3id", "description")
        test_plasmid3.insert_scar(scar1,"T1","after")
        self.assertEqual(test_plasmid3.t1.structure, [test_plasmid3.t1_part, scar1])

    def test_plasmid_insert_scar_before_T0(self):
        scar1 = Part("Scar","Scar1","AAA","scar")
        test_plasmid3 = Plasmid("plasmid3", "plamid3id", "description")
        test_plasmid3.insert_scar(scar1,"T0","before")
        self.assertEqual(test_plasmid3.t0.structure, [scar1, test_plasmid3.t0_part, test_plasmid3.SanDI])


    def test_plasmid_insert_scar_after_T0(self):
        scar1 = Part("Scar","Scar1","AAA","scar")
        test_plasmid3 = Plasmid("plasmid3", "plamid3id", "description")
        test_plasmid3.insert_scar(scar1,"T0","after")
        self.assertEqual(test_plasmid3.t0.structure, [test_plasmid3.t0_part,test_plasmid3.SanDI,scar1,])


    def test_plasmid_insert_scar_before_orit(self):
        scar1 = Part("Scar","Scar1","AAA","scar")
        test_plasmid3 = Plasmid("plasmid3", "plamid3id", "description")
        test_plasmid3.insert_scar(scar1,"OriT","before")
        self.assertEqual(test_plasmid3.orit.structure, [scar1, test_plasmid3.orit_part])


    def test_plasmid_insert_scar_after_orit(self):
        scar1 = Part("Scar","Scar1","AAA","scar")
        test_plasmid3 = Plasmid("plasmid3", "plamid3id", "description")
        test_plasmid3.insert_scar(scar1,"OriT","after")
        self.assertEqual(test_plasmid3.orit.structure, [test_plasmid3.orit_part,scar1])


    def test_plasmid_gather_data(self):
        test_dict = TestPlasmid.test_plasmid2.gather_data()
        reference_dict = {'PacI': [0, 'paci', 're', '', 5], 'name': [67, 'id3', 'gene', 'description3', 10], 'SpeI': [11, 'spei', 're', '', 16], 't0_part': [17, 't0_part', 'terminator', '', 119], 'SanDI_part': [120, 'sandi', 're', '', 126], 'swai': [127, '1', 're', '', 123], 'pshai': [138, '1', 're', '', 137], 'orit_part': [144, 'orit_part', 'oriT', '', 143], 'FseI': [390, 'fsei', 're', '', 389], 'AscI': [401, 'asci', 're', '', 66], 't1_part': [409, 't1_part', 'terminator', '', 400]}
        self.assertEqual(test_dict,reference_dict)


if __name__ == '__main__':
    unittest.main()