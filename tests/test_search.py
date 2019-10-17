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