def write_tuples_in_bed(tuples, filename, append=False):
    file_mode = 'a' if append else 'w'
    with open(filename, file_mode) as bed:
        for tupl in tuples:
            print(*tupl, file=bed, sep="\t")


def write_positions_in_bed(positions, filename, label=None, append=False):
    # chr is pos[0], loc is pos[1]
    if label is None:
        write_tuples_in_bed(map(lambda pos: (pos[0], pos[1]-1, pos[1]), positions),
                            filename, append)
    else:
        write_tuples_in_bed(map(lambda pos: (pos[0], pos[1]-1, pos[1], label), positions),
                            filename, append)
