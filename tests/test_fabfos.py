import unittest


class MyTestCase(unittest.TestCase):
    # TODO: flesh out basic integrative tests
    def test_something(self):
        self.assertEqual(True, False)

    # ./FabFos.py -m tests/test_data/test_SE_miffed.csv -b tests/test_data/trim_sequences.fasta -r tests/test_data/ -a spades --fabfos_path tests/test_data/FabFos_DB/ -T 4 -p se --force
    # ./FabFos.py -m tests/test_data/test_PE_miffed.csv -b tests/test_data/trim_sequences.fasta -r tests/test_data/fwd/ -2 tests/test_data/rev/ -a megahit --fabfos_path tests/test_data/FabFos_DB/ -T 4 --force

    # TODO: test end mapping workflow


if __name__ == '__main__':
    unittest.main()
