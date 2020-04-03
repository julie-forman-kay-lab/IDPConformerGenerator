import os

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


if __name__ == "__main__":
    data = read_data_into_dict("/Users/alaashamandy/IDPCalcPDBDownloader/alphas/data/")
    print(data.items())