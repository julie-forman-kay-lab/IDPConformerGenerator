import os
from idpconfgen import conformer_generator

def read_data_into_dict(folder_path):
    """
    Return the data in the following format
    {
        "L": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1]}, ...,
        "H": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1]}, ...,
        "S": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1]}, ...,
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

class Filter:
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.data = read_data_into_dict(folder_path)

    def return_ss(self, ss, lower_bound=None, upper_bound=None):
        """
        This function returns all the loop fragments
        that have lower_bound<=size<=upper_bound.
        """
        # if ss not in self.data:
            # TODO: raise exception

        if not lower_bound and not upper_bound:
            #return all sizes
            return {ss :self.data[ss]}
        return {ss :self.__return_ss_ranges(ss, lower_bound, upper_bound)}

    def __return_ss_ranges(self, ss, lower_bound, upper_bound):
        ss_data = self.data[ss]
        return dict(filter(lambda fragment: len(fragment[0]) >= lower_bound and len(fragment[0]) <= upper_bound, ss_data.items()))
    
    # def find_aa_patterns(self, input_pattern, min_seq_chunk_size=None, max_mismatch=None, data=None):
    #     """
    #     Finds patterns of input_pattern in the data using the ProteinSearch algorithm.

    #     The input data has to be of the following form:
    #     {
    #         "L": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1]}, ...,
    #         "H": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1]}, ...,
    #         "S": {seq1 : [x,y,z,phi,psi,omega,chi1], seq2: [x,y,z,phi,psi,omega,chi1]}, ...,
    #     }
    #     """
    #     database = self.data
    #     if data:
    #         database = data

    #     protein_search = conformer_generator.ProteinSearch()
    #     for _, sequences in database.items():
    #         for fragment, fragment_data in sequences.items():
    #             protein_search.start_search(
    #                                             input_pattern,
    #                                             fragment,
    #                                             min_seq_chunk_size=min_seq_chunk_size,
    #                                             max_mismatch=max_mismatch
    #                                         )
                
    #     return protein_search
                

    


if __name__ == "__main__":
    filter_step = Filter("/Users/alaashamandy/IDPCalcPDBDownloader/alphas/data/")
    print(filter_step.find_aa_patterns("YCV", 3, 0).results)
    # print(filter_step.return_ss("H", 3, 9))

    # protein_search = conformer_generator.ProteinSearch()
    # protein_search.start_search("AAA", "AAA")
    # print(protein_search.results[0]["AAA"])

