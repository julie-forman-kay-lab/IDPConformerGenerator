from collections import defaultdict
import logging
import math


class ProteinSearch:

    def __init__(self):
        # both lists should be the same size
        self.proteins_sequences = []
        self.results = []

    def start_search(
            self,
            input_pattern,
            primary_seq_db,
            min_seq_chunk_size=3,
            max_mismatch=0
            ):
        """ search every sequential combination of
         sequences in pattern, inside of word """

        def recursive_search(pattern_index, pattern_size):
            """ recursively increase bracket_size
             and search through the primary_seq_db"""

            # base case
            if pattern_index + pattern_size > input_pattern_length:
                return
            
            current_pattern = \
                input_pattern[pattern_index:pattern_index + pattern_size]

            previous_pattern = current_pattern[:-1]

            # optimization using dynamic programming
            if current_pattern in result:
                recursive_search(pattern_index, pattern_size + 1)

            for index, prev_mismatch_num in result[previous_pattern]:

                try:
                    character = primary_seq_db[index + pattern_size - 1]
                except IndexError:  # index is at the end of the sequence
                    continue
                
                mismatch = math.ceil(((prev_mismatch_num + 1)
                                        / (pattern_size)) * 100)

                # check if the characters are the same,
                # if they're not make sure we're still
                # under the max_mismatch specified
                if current_pattern[-1] == character or mismatch <= max_mismatch:

                    same_char_bool = 0 if current_pattern[- 1] == character else 1
                    new_mismatch = prev_mismatch_num + same_char_bool
                    new_result = result[current_pattern]
                    if not new_result or index not in list(zip(*new_result))[0]:
                        result[current_pattern].append((index, new_mismatch))

                    # have we JUST reached the max because we found a mismatch?
                    if mismatch == max_mismatch and new_mismatch != prev_mismatch_num:
                        continue

                    recursive_search(pattern_index, pattern_size + 1)

        try:
            primary_seq_db_length = len(primary_seq_db)
            assert primary_seq_db_length != 0
            input_pattern_length = len(input_pattern)
        except TypeError:
            logging.exception("One of the parameters is None")
            return
        except AssertionError:
            logging.exception("The parameter primary_seq_db is an empty String")
            return
        
        result = defaultdict(lambda: [])
        primary_seq_db_dict = defaultdict(lambda: [])

        # splitting the primary_seq_db sequence and
        # initializing primary_seq_db_dict to contain minimum brackets
        maximum_allowed_index = primary_seq_db_length - min_seq_chunk_size
        for bracket_index in range(maximum_allowed_index + 1):
            minimum_sequence = \
                primary_seq_db[bracket_index:bracket_index + min_seq_chunk_size]
            primary_seq_db_dict[minimum_sequence].append(bracket_index)

        maximum_allowed_index = input_pattern_length - min_seq_chunk_size
        for bracket_index in range(maximum_allowed_index + 1):
            minimum_sequence = \
                input_pattern[bracket_index: bracket_index + min_seq_chunk_size]

            # has to has been seen before in the primary_seq_db
            if minimum_sequence in primary_seq_db_dict:
                indices = primary_seq_db_dict[minimum_sequence]

                # has not been seen before
                if minimum_sequence not in self.results:
                    # initiallize all to 0 mismatches
                    result[minimum_sequence] = [(index, 0) for index in indices]
                
                recursive_search(bracket_index, min_seq_chunk_size + 1)

        self.proteins_sequences.append(primary_seq_db)
        self.results.append(result)

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

#     input_seq = "GAYHLRHEVKEDMAMGCIKSKGKDSLSDDGVDLKTQPVRNTERTIYVRDPTSNKQQRPVPESQLLPGQRFQTKDPE"

#     print("generating random sequence..")
#     start = time.time()
#     primary_seq_db = string_generator()
#     middle = time.time()
#     stopWatch(middle-start)

#     search = ProteinSearch()
#     print("searching..")
#     search.start_search(input_seq, primary_seq_db, 5, 25)
#     end = time.time()
#     stopWatch(end-middle)
#     print("done")
#     print(len(search.results[-1]))
#     search.start_search(input_seq, primary_seq_db, 5, 0)
#     print(len(search.results[-1]))

    # print("\ndifference:\n")
    # value = { k : search.results[-2][k] for k in set(search.results[-2]) - set(search.results[-1]) }
    # print(value)

# ----------------NOTES------------------
# We don't find every combination of mismatches. we only go far enough until the mismatch has gone
# over the max_mismatch and stop, we don't consider situations where there might be a match after,
# and therefore the percentage of mismatches will lower. Reason for this is we shouldn't stop
# at an alpha-helix, we should start and end inside of a loop
