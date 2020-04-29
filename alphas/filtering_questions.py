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

def read_data_into_sanwiched_dict(folder_path, ss_sandwich, min_sanwich_length, ss_burger, min_burger_length):
    """
    This function adds a "sanwiched" secondary structure
    to the data. For example if ss_sandwich is "L" with
    min length 10 and ss_burger is "H" with min length 4 
    then an acceptable aa sequence will be 
    "LLLLLLLLLLHHHHHLLLLLLLLLLL".

    Returns the data in the following format
    {
        ss_sandwich + ss_burger: {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1], ...}
    }
    """
    data_dict = defaultdict(lambda: [])
    with os.scandir(folder_path) as it:
        for entry in it:
            if entry.name.endswith(".data") and entry.is_file():
                data = open(entry.path, "r").readlines()
                data =  [line.strip().split(",") for line in data]
                fragment_ss = ""
                fragment_aa = ""
                fragment_data = []

                for line in data:
                    ss = line[1]
                    if ss == ss_sandwich or ss == ss_burger:
                        fragment_ss += ss
                        fragment_aa += line[0]
                        fragment_data.append(line[2:])
                    else:
                        data_dict[(fragment_ss, fragment_aa)] = fragment_data
                        fragment_ss = ""
                        fragment_aa = ""
                        fragment_data = []
    result = remove_illegal_fragments(data_dict, ss_sandwich, min_sanwich_length, ss_burger, min_burger_length)
    return {ss_sandwich + " + " + ss_burger: result}
    

def remove_illegal_fragments(data, ss_sandwich, min_sanwich_length, ss_burger, min_burger_length):
    
    results = defaultdict(lambda: [])
    for fragment, fragment_data in data.items():

        fragment_ss = fragment[0]
        fragment_aa = fragment[1]

        chopped = fragment_ss.split(ss_sandwich*min_sanwich_length)
        cur_fragment_aa = ""
        cur_fragment_index = 0
        start_index = 0

        # Taking care of the first sub_fragment
        if chopped[0].find(ss_burger) >= 0:
            # illegal fragment
            fragment_data = fragment_data[len(chopped[0])*7:]
            fragment_aa = fragment_aa[len(chopped[0]):]
            chopped[0] = ""
            
        for sub_fragment in chopped[:-1]:
            if sub_fragment:
                new_sub_fragment = sub_fragment
                index_of_burger = sub_fragment.find(ss_burger) 
                if ss_burger != 0: #check if there is any ss_sandwich residue in the beginning
                    #add and start from ss_burger
                    cur_fragment_aa += fragment_aa[cur_fragment_index:cur_fragment_index+index_of_burger]
                    cur_fragment_index += index_of_burger
                    new_sub_fragment = sub_fragment[index_of_burger:]
                
                # only contains ss_burger, check if within legal size
                len_ss_burger = len(new_sub_fragment)
                if len_ss_burger >= min_burger_length:
                    cur_fragment_aa += fragment_aa[cur_fragment_index:cur_fragment_index+len_ss_burger]
                    cur_fragment_index += len_ss_burger
                else:
                    results[cur_fragment_aa].append(fragment_data[start_index*7: (start_index+len(cur_fragment_aa)*7)])
                    cur_fragment_index += len(new_sub_fragment)
                    start_index = cur_fragment_index
                    cur_fragment_aa = ""

            # add the next min_sanwich_length amino acids
            # because we know for sure they are ss_sandwich
            cur_fragment_aa += fragment_aa[cur_fragment_index: cur_fragment_index+min_sanwich_length]
            cur_fragment_index+=min_sanwich_length

        # Taking care of the last fragment
        if fragment_aa:
            last_sub_fragment = chopped[-1]
            for ss in last_sub_fragment:
                if ss != ss_sandwich:
                    break
                cur_fragment_aa += fragment_aa[cur_fragment_index]
                cur_fragment_index+=1
            results[cur_fragment_aa] = fragment_data[start_index*7: (start_index+len(cur_fragment_aa)*7)]
    return results

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

    def apply_ss(self, ss, lower_bound=None, upper_bound=None):
        """
        This function filters the data to have ss
        fragments that have lower_bound<=size<=
        upper_bound. Where ss can be loops, helices
        or sheets.
        """
        if ss not in self.filtered_data:
            return {}

        if not lower_bound or not upper_bound:
            #return all sizes
            self.filtered_data = {ss :self.filtered_data[ss]}
        self.filtered_data = {ss :self.__return_ss_ranges(ss, lower_bound, upper_bound)}

    def __return_ss_ranges(self, ss, lower_bound, upper_bound):
        ss_data = self.filtered_data[ss]
        return dict(filter(lambda fragment: len(fragment[0]) >= lower_bound and len(fragment[0]) <= upper_bound, ss_data.items()))

    def filter_ss_sandwich(self, ss_sandwich, min_sanwich_length, ss_burger, min_burger_length, keep_others=True):
        """
        NOTE: this function will add to a new secondary structure
        (ss_sandwich + ss_burger) to self.filtered_data, therefore
        there may be some overlap between the secondary structures.
        keep_others is a boolean that either keeps this overlap
        or removes it by making this sandwich the only ss in self.filtered.
        """
        data = read_data_into_sanwiched_dict(self.folder_path, ss_sandwich, min_sanwich_length, ss_burger, min_burger_length)
        new_ss = ss_sandwich + " + " + ss_burger

        if not keep_others:
            self.clear_filter()
        self.filtered_data[new_ss] = data[new_ss]

    def clear_filter(self):
        self.filtered_data = dict()

    def reset_filter(self):
        self.filtered_data = self.data

    def return_filter(self):
        return self.filtered_data




