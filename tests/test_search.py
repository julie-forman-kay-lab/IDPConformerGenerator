"""Test ProteinSearch."""

from idpconfgen.conformer_generator import ProteinSearch


class TestSearchNoMismatch:

    def test_class(self):
        proteins = ProteinSearch()
        proteins.start_search("ACD", "ACD")
        assert proteins.proteins_sequences[-1] == "ACD"

        proteins.start_search("KED", "KED")
        assert proteins.proteins_sequences[-1] == "KED"

    def test_edgecase1(self):
        proteins = ProteinSearch()
        proteins.start_search("", "ACD")
        assert len(proteins.results[-1]) == 0

    def test_edgecase2(self):
        proteins = ProteinSearch()
        proteins.start_search("", "")
        assert len(proteins.results) == 0
        assert len(proteins.proteins_sequences) == 0

        proteins.start_search("ACD", None)
        assert len(proteins.results) == 0
        assert len(proteins.proteins_sequences) == 0

    def test_edgecase3(self):
        proteins = ProteinSearch()
        proteins.start_search("ACD", "ACD")
        assert len(proteins.results[-1]) == 1

    def test_edgecase4(self):
        proteins = ProteinSearch()
        proteins.start_search("ACDED", "ACDED")
        assert len(proteins.results[-1]) == 6

        proteins.start_search("ACDEDHPCTRWNK", "ACDED")
        assert len(proteins.results[-1]) == 6


class TestSearchWithMismatch:

    def test_mismatch1(self):
        proteins = ProteinSearch()
        proteins.start_search("ACDDKC", "ACDMKC", 3, 10)
        # should never find any mismatches since the minimum percentage
        # for this minimum chunk would be 25%
        assert len(proteins.results[-1]) == 1

    def test_mismatch2(self):
        proteins = ProteinSearch()
        proteins.start_search("ACDDKC", "ACDMKC", 3, 25)
        print(proteins.results[-1])
        assert len(proteins.results[-1]) == 2

    def test_mismatch3(self):
        proteins = ProteinSearch()
        proteins.start_search("ACDDK", "ACDMK", 3, 50)
        assert len(proteins.results[-1]) == 3

    def test_mismatch4(self):
        proteins = ProteinSearch()
        proteins.start_search("ACDDKMACD", "ACDDKOKOKOKMACD", 5, 0)
        assert len(proteins.results[-1]) == 2

        proteins.start_search("ACDDKMACD", "ACDDKOKOKOKMACD", 5, 45)
        assert len(proteins.results[-1]) == 6
