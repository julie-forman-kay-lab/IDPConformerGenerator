import os
from idpconfgen import conformer_generator
from collections import defaultdict

def read_data_into_dict(folder_path):
    """
    Return the data in the following format
    {
        "L": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...},
        "H": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...},
        "S": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...},
    }
    """


    result_dict = {
                    "L": {},
                    "H": {},
                    "S": {},
                    "NEW": {},
                    }

    with os.scandir(folder_path) as it:
        for entry in it:
            if entry.name.endswith(".data") and entry.is_file():
                data = open(entry.path, "r").readlines()
                data =  [line.strip().split(",") for line in data]
                prev_SS = "NEW"
                prev_fragment = ""
                prev_fragment_data = []
                for line in data:
                    cur_SS = line[1]
                    if cur_SS == prev_SS: #we're still in the same fragment
                        prev_fragment += line[0]
                        prev_fragment_data.extend(line[2:])
                    else: #new fragment found

                        # save the old fragment
                        result_dict[prev_SS][prev_fragment] = prev_fragment_data

                        # make new fragment
                        prev_SS = cur_SS
                        prev_fragment = line[0]
                        prev_fragment_data = line[2:]

                # save the old fragment
                result_dict[prev_SS][prev_fragment] = prev_fragment_data
    del result_dict["NEW"]

    return result_dict

def find_aa_patterns(input_pattern, database, min_seq_chunk_size=None, max_mismatch=None):
        """
        Finds patterns of input_pattern in the data using the ProteinSearch algorithm.

        The input data has to be of the following form:
        {
            "L": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...},
            "H": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...},
            "S": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...}
        }

        The output data will be of the following form:
        {
            "L" : {match_seq1 : [ [x,y,z,phi,psi,omega,chi1], [x,y,z,phi,psi,omega,chi1] ], match_seq2: ...},
            "H" : {match_seq1 : [ [x,y,z,phi,psi,omega,chi1], [x,y,z,phi,psi,omega,chi1] ], match_seq2: ...},
            "S" : {match_seq1 : [ [x,y,z,phi,psi,omega,chi1], [x,y,z,phi,psi,omega,chi1] ], match_seq2: ...},
        }
        """
        result = {}
        protein_search = conformer_generator.ProteinSearch()

        for ss, sequences in database.items():
            result[ss] = defaultdict(lambda: [])
            for fragment, fragment_data in sequences.items():
                protein_search.start_search(
                                            input_pattern,
                                            fragment,
                                            min_seq_chunk_size=min_seq_chunk_size,
                                            max_mismatch=max_mismatch,
                                            )
                # most recent result
                search_result = protein_search.results[-1]
                for matched_seq, positions in search_result.items():
                    for position, _ in positions:
                        # grab the fragment data that we found
                        result[ss][matched_seq].append(fragment_data[position*7:(position+len(matched_seq))*7])
                
        return result
                

class Filter:
    """
    This class has been designed such that 
    the output of any of its' functions can 
    be the input to any of its' functions.
    This will allow the user to use a 
    combination of different filters togther.
    """
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.data = read_data_into_dict(folder_path)
        self.filtered_data = self.data

    def return_ss(self, ss, lower_bound=None, upper_bound=None):
        """
        This function returns all the ss fragments
        that have lower_bound<=size<=upper_bound.
        Where ss can be loops, helixes or sheets
        """
        if ss not in self.filtered_data:
            return {}

        if not lower_bound and not upper_bound:
            #return all sizes
            self.filtered_data = {ss :self.filtered_data[ss]}
        self.filtered_data = {ss :self.__return_ss_ranges(ss, lower_bound, upper_bound)}

    def __return_ss_ranges(self, ss, lower_bound, upper_bound):
        ss_data = self.filtered_data[ss]
        return dict(filter(lambda fragment: len(fragment[0]) >= lower_bound and len(fragment[0]) <= upper_bound, ss_data.items()))
    
    def clear_filter(self):
        self.filtered_data = self.data

    def return_filter(self):
        return self.filtered_data


if __name__ == "__main__":
    filter_step = Filter("/Users/alaashamandy/IDPCalcPDBDownloader/alphas/data/")
    res = find_aa_patterns("YCV",  filter_step.return_filter(),3, 0)
    print(res)


