def bind_to_chrom(tuple, binding_molecules, object):
    inserted=False
    if binding_molecules:
        for i in range(len(binding_molecules[0])):
            if int(tuple[0]) < int(binding_molecules[0][i][0]):
                binding_molecules[0].insert(i, tuple)
                binding_molecules[1].insert(i, object)
                inserted=True
    if inserted==False:
        binding_molecules[0].append(tuple)
        binding_molecules[1].append(object)
    return binding_molecules


list=[[(234,238), (456,576), (598, 627)],['bla', 'bli', 'blub']]

print(bind_to_chrom((45435,300065), [[],[]], 'blaaaab'))