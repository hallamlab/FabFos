import os
import unittest


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.test_data_dir = "test_data"
        self.backbone = os.path.join(self.test_data_dir, "trim_sequences.fasta")
        self.fabfos_db_path = os.path.join(self.test_data_dir, "FabFos_DB")
        return

    # TODO: flesh out basic integrative tests
    def test_fabfos_se(self):
        from src.FabFos import fabfos_main
        fabfos_main(["--miffed", "test_data/test_SE_miffed.csv",
                     "--background", self.backbone,
                     "--reads", self.test_data_dir,
                     "--assembler", "spades",
                     "--fabfos_path", self.fabfos_db_path,
                     "--threads", str(4),
                     "-p", "se",
                     "--force"])
        self.assertEqual(True, True)
        return

    def test_fabfos_pe(self):
        from src.FabFos import fabfos_main
        fabfos_main("-m tests/test_data/test_PE_miffed.csv "
                    "-b tests/test_data/trim_sequences.fasta "
                    "-r tests/test_data/fwd/ -2 tests/test_data/rev/ "
                    "-a megahit --fabfos_path tests/test_data/FabFos_DB/ -T 4 --force".split())
        self.assertTrue(True)
        return

    # TODO: test end mapping workflow
    def test_get_assembly_nx(self):
        from src.FabFos import get_assembly_nx
        test_fasta = {"s1": "A"*100, "s2": "C"*200, "s3": "G"*400, "s4": "T"*800}
        nx_stats = get_assembly_nx(test_fasta)
        self.assertEqual(11, len(nx_stats))
        self.assertEqual(800, nx_stats[0])
        self.assertEqual(800, nx_stats[0.5])
        self.assertEqual(400, nx_stats[0.6])
        self.assertEqual(100, nx_stats[1.0])
        return

    def test_validate_dependency_versions(self):
        from src.FabFos import validate_dependency_versions
        valid = validate_dependency_versions({"samtools": "1.11", "trimmomatic.jar": "0.40"})
        self.assertTrue(valid)
        return


if __name__ == '__main__':
    unittest.main()
