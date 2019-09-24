from collections import defaultdict
import logging


class Protein_Search:

    def __init__(self):
        # both lists should be the same size
        self.proteins_sequences = []
        self.results = []

    def start_search(self, input_pattern, protein_bank):
        """ search every sequential combination of
         sequences in pattern, inside of word """

        def recursive_search(pattern_index, pattern_size):
            """ recursively increase bracket_size
             and search through the protein_bank"""

            # base case
            if pattern_index + pattern_size > input_pattern_length:
                return

            current_pattern = input_pattern[pattern_index:pattern_index + pattern_size]
            previous_pattern = current_pattern[:-1]

            # optimization using dynamic programming
            if current_pattern in result:
                recursive_search(pattern_index, pattern_size + 1)

            for index in result[previous_pattern]:

                try:
                    character = protein_bank[index + pattern_size - 1]
                except IndexError:  # index is at the end of the sequence
                    continue

                # found a match
                if current_pattern[-1] == character:
                    if index not in result[current_pattern]:
                        result[current_pattern].append(index)
                    recursive_search(pattern_index, pattern_size + 1)

        try:
            protein_bank_length = len(protein_bank)
            assert protein_bank_length != 0
            input_pattern_length = len(input_pattern)
        except TypeError:
            logging.exception("One of the parameters is None")
            return
        except AssertionError:
            logging.exception("The parameter protein_bank is an empty String")
            return
        
        result = defaultdict(lambda: [])
        protein_bank_dict = defaultdict(lambda: [])
        minimum_sequence_length = 3

        # splitting the protein_bank sequence and
        # initializing protein_bank_dict to contain minimum brackets
        for bracket_index in range(protein_bank_length - minimum_sequence_length + 1):
            minimum_sequence = protein_bank[bracket_index: bracket_index + minimum_sequence_length]
            protein_bank_dict[minimum_sequence].append(bracket_index)

        for bracket_index in range(input_pattern_length - 2):
            minimum_sequence = input_pattern[bracket_index: bracket_index + minimum_sequence_length]

            # has to has been seen before in the protein_bank
            if minimum_sequence in protein_bank_dict:
                indices = protein_bank_dict[minimum_sequence]

                # has not been seen before
                if minimum_sequence not in self.results:
                    result[minimum_sequence] = indices
                
                recursive_search(bracket_index, minimum_sequence_length + 1)

        self.proteins_sequences.append(protein_bank)
        self.results.append(result)
        

class TestSearch:

    def test_class(self):
        proteins = Protein_Search()
        proteins.start_search("ACD", "ACD")
        assert proteins.proteins_sequences[-1] == "ACD"

        proteins.start_search("KED", "KED")
        assert proteins.proteins_sequences[-1] == "KED"

    def test_edgecase1(self):
        proteins = Protein_Search()
        proteins.start_search("", "ACD")
        assert len(proteins.results[-1]) == 0

    def test_edgecase2(self):
        proteins = Protein_Search()
        proteins.start_search("", "")
        assert len(proteins.results) == 0
        assert len(proteins.proteins_sequences) == 0

        proteins.start_search("ACD", None)
        assert len(proteins.results) == 0
        assert len(proteins.proteins_sequences) == 0

    def test_edgecase3(self):
        proteins = Protein_Search()
        proteins.start_search("ACD", "ACD")
        assert len(proteins.results[-1]) == 1

    def test_edgecase4(self):
        proteins = Protein_Search()
        proteins.start_search("ACDED", "ACDED")
        assert len(proteins.results[-1]) == 6

        proteins.start_search("ACDEDHPCTRWNK", "ACDED")
        assert len(proteins.results[-1]) == 6

# -------------- FOR MANUAL TESTING ONLY USE BELOW --------------------

# def string_generator():
#     string = "GAYHLRHEVKEDMALAAAAAAAPNSLWWITVICIHPCTRWNKTAPEMAQVWWWGVQERMFAWSPMLLLKFFAPLFKGLCFTMRAQPGRAMKQNLVDGEWTRKYVCIRDRKSSNHRSWVSMFNYQLMCSEGATSIRIRPFHWQWPYEWMYHDNYHKSWWLSVEEWCHWFNRNDNGSLNHTMPCRDYCINDASHHKIWHGVILATTTTTTKIHEWHMSVNAWMGLAKLAIMLHHFSNVMGLTQMETLCQHATHEIMFVMWVIDFMNGPAQQHDKDFEDVSEYPKLFQNDPVSPEQYWIGEPDDQLKYTYDITVVFAMNDNMCVTGNPDHIVHHHFKNWANHVNHNPFEQNDDDASYGIKHEQWCYNVDMDVHKGWRWPSFKKYYTTNAEAPNFSDNFISEDQMLEELVDIQD"
#     randomize = "GAYHLRHEVKEDMNSLWWITVICIHPCTRWNKTAPEMAQVWWWGVQERMFAWSPMLLLKFFAPLFKGLCFTMRAQPGRAMKQNLVDGEWTRKYVCIRDRKSSNHRSWVSMFNYQLMCSEGATSIRIRPFHWQWPYEWMYHDNYHKSWWLSVEEWCHWFNRNDNGSLNHTMPCRDYCINDASHHKIWHGVILATTTTTTKIHEWHMSVNAWMGLAKLAIMLHHFSNVMGLTQMETLCQHATHEIMFVMWVIDFMNGPAQQHDKDFEDVSEYPKLFQNDPVSPEQYWIGEPDDQLKYTYDITVVFAMNDNMCVTGNPDHIVHHHFKNWANHVNHNPFEQNDDDASYGIKHEQWCYNVDMDVHKGWRWPSFKKYYTTNAEAPNFSDNFISEDQMLEELVDIQD"
#     for _ in range(28000):
#         string+=(''.join(random.sample(randomize,len(randomize))))
#     return string

# def stopWatch(value):
#     '''From seconds to Days;Hours:Minutes;Seconds'''

#     valueD = (((value/365)/24)/60)
#     Days = int (valueD)

#     valueH = (valueD-Days)*365
#     Hours = int(valueH)

#     valueM = (valueH - Hours)*24
#     Minutes = int(valueM)

#     valueS = (valueM - Minutes)*60
#     Seconds = int(valueS)
#     print(Hours, " Hours", Minutes, " Minutes" ,Seconds, " seconds")

# if __name__ == '__main__':
#     # 400 residues
#     input_seq = "YCSMLILDQNWALAAAAAAPSEGSSLRCCGCSNMMLFASFMMCKEYHQQYNWREWKSTNTACCHMWHISKLHGYRVQLDKCDLYHQDDDRTDRYWDIINLADEEIRSRGSPDQKSWCQSMGRYQVWSMQSAWGCLEPSPRKLMEPPYYRGMSKWEKDSKYSNRTMFAPCADRESRYWEMYQRQKIEPKRLDFAQTCRWLTVASRIFTQWICWWPKKQDMPVKISLQMGKGGISEWAVLSALFGQGMNQVKGKMPVKQYCKYNCGWEAYNKSPTWVVKGMAMMQGLTRTGSNDYGHAMWDTLKWEAVFCFRQTRQWRWWIECIGKYWHYWVFYWDVQRLAVTLTRYGRRDYEPGAGWSNQDAVINRFDAYHATFYPAIGEYPSLGGGWAMNQDRQWWQMCFSSCPQGHEM"

#     print("generating random sequence..")
#     # start = time.time()
#     # protein_bank = string_generator()
#     # middle = time.time()
#     # stopWatch(middle-start)

#     search = Protein_Search()
#     print("searching..")
#     search.start_search(input_seq, "")
#     # end = time.time()
#     # stopWatch(end-middle)
#     print("done")
