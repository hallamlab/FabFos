import os
import unittest


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.test_data_dir = "test_data"
        self.backbone = os.path.join(self.test_data_dir, "trim_sequences.fasta")
        self.fabfos_db_path = os.path.join(self.test_data_dir, "FabFos_DB")
        self.fq_interleaved = os.path.join(self.test_data_dir, "fwd",
                                           "52414.2.321701.AACCGTTC-GAACGGTT.filter-METAGENOME.fosmid.fastq")
        return

    def test_fabfos_se(self):
        from src.FabFos import fabfos_main
        fabfos_main(["--miffed", self.test_data_dir + os.sep + "test_SE_miffed.csv",
                     "--background", self.backbone,
                     "--reads", self.test_data_dir,
                     "--assembler", "spades_sc",
                     "--fabfos_path", self.fabfos_db_path,
                     "--threads", str(4),
                     "-p", "se",
                     "--force", "--overwrite"])
        self.assertEqual(True, True)
        return

    def test_fabfos_pe(self):
        from src.FabFos import fabfos_main
        fabfos_main(["-m", self.test_data_dir + os.sep + "test_PE_miffed.csv",
                     "-b", self.backbone,
                     "-r", self.test_data_dir + os.sep + "fwd/", "-2", self.test_data_dir + os.sep + "rev/",
                     "-a", "spades_meta",
                     "--fabfos_path", self.fabfos_db_path] +
                    "-T 4 --force --overwrite".split())
        self.assertTrue(True)
        return

    def test_fabfos_main_interleaved(self):
        from src.FabFos import fabfos_main
        fabfos_main(["-m", self.test_data_dir + os.sep + "miffed_monoaromatics.csv",
                     "-b", self.backbone,
                     "-r", self.test_data_dir + os.sep + "fwd",
                     "-a", "megahit",
                     "--fabfos_path", self.fabfos_db_path,
                     "--interleaved", "--overwrite", "-T", str(4), "--verbose"])
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

    def test_summarize_dependency_versions(self):
        from src.FabFos import summarize_dependency_versions
        summarize_dependency_versions({})
        return

    def test_deinterleave_fastq(self):
        from src.FabFos import deinterleave_fastq, read_fastq_to_dict
        num_interleaved = len(read_fastq_to_dict(self.fq_interleaved))
        fwd_fq, rev_fq = deinterleave_fastq(self.fq_interleaved, output_dir=self.test_data_dir)
        num_split = sum([len(read_fastq_to_dict(fwd_fq)), len(read_fastq_to_dict(rev_fq))])
        self.assertEqual(num_interleaved, num_split)
        return


if __name__ == '__main__':
    unittest.main()
