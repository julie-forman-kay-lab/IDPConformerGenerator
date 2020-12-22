
# from cli_build.py
# vie 18 dic 2020 16:34:40 EST

ss = [
    [
        np.where(residue_numbers == _resnum)[0][4:],
        np.full((_natoms - 4, 3), NAN, dtype=np.float64),
        ]
    for _resnum, _natoms in sorted(Counter(residue_numbers).items())
    ]

# special cases to the above implementation
# the first residues contains H1 H2 H3 atoms instead of H
ss[0][0] = np.where(residue_numbers == residue_numbers[0])[0][6:]
ss[0][1] = np.full((len(ss[0][0]), 3), NAN, dtype=np.float64)

# the last residue will contain an extra atom OXT, which needs to be
# removed
ss[-1][0] = ss[-1][0][:-1]
ss[-1][1] = ss[-1][1][:-1]
assert ss[-1][0].size == ss[-1][1].shape[0]



# NOT USED
def calc_outer_sum_upper_diagonal(data):
    """
    Calculate the sig_ij for all atom pairs.

    See Amber20 manual equation 14.13 on page 256.
    """
    # require
    assert data.ndim == 1, 'Array should have only one dimension.'
    assert data.dtype in (np.int, np.float), data.dtype

    indices = np.triu_indices(data.size, k=+1)

    result = (data[:, None] + data)[indices]

    # ensure
    _expected_size = (data.size * data.size - data.size) // 2
    assert result.size == _expected_size, f'Sizes differ from expectations: {(result.size, _expected_size)}'
    assert abs(result[0] - (data[0] + data[1])) < 0.0000001, 'Values are not as expected.'
    assert abs(result[-1] - (data[-2] + data[-1])) < 0.0000001, 'Values are not as expected.'

    return result


# NOT USED
def generate_ij_pairs_without_repetition(*data):
    """."""
    assert all(len(data[i]) > 0 for i in range(len(data))), \
        f'Data index {i} is empty, we expect some data in.'
    assert sum(len(d) for d in data) // len(data) == len(data[0]), \
        'All datasets must have the same length'

    ij_pairs = [[] for i in range(len(data))]
    ij_appends = [l.append for l in ij_pairs]

    # the first two loops represent the upper matrix diagonal
    len_ = len(data[0])

    for i in range(len_ - 1):
        for j in range(i + 1, len_):

            # this loops adds the pairs to the multiple data lists
            for append_func, series in zip(ij_appends, data):
                append_func((series[i], series[j]))

    assert ij_pairs[-1][-1] == (data[-1][-2], data[-1][-1])

    assert isinstance(ij_pairs, list)
    return ij_pairs


# NOT USED
def create_atom_pair_connectivity_masks(
        res_names_ij,
        res_num_ij,
        atom_names_ij,
        connectivity_intra,
        connectivity_inter,
        ):
    """
    To generate bonds masks we need the residue numbers and the atom names.

    Depends
    -------
    `are_connected`
    """
    assert len(res_names_ij) == len(res_num_ij) == len(atom_names_ij), \
        'Sizes of `ij` vectors differ.'

    zipit = zip(res_names_ij, res_num_ij, atom_names_ij)

    # boolean matrix arrays
    # all are false until otherwise stated
    connectivity_mask = np.full(len(res_names_ij), False)
    counter = 0
    for (rn1, rn2), (n1, n2), (a1, a2) in zipit:

        found_connectivity = are_connected(
            int(n1),
            int(n2),
            rn1,
            a1,
            a2,
            connectivity_intra,
            connectivity_inter,
            )

        if found_connectivity:
            connectivity_mask[counter] = True

        counter += 1

    assert isinstance(connectivity_mask, np.ndarray)
    return connectivity_mask


# NOT USED
def calc_outer_multiplication_upper_diagonal(data):
    """Calculate electrostatic interactions for ij pair."""
    assert data.ndim == 1, \
        'Data has {data.ndim} dimensions, should have only 1.'

    indices = np.triu_indices(data.size, k=+1)

    result = np.outer(data, data)[indices]

    # ensure
    assert result.size == (data.size * data.size - data.size) // 2
    assert abs(result[0] - data[0] * data[1]) < 0.0000001

    return result


