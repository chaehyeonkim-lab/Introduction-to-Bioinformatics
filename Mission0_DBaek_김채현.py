# Mission0. ch22의 A,G,C,T 개수와 비율 측정
# 2022-11001 김채현

# counting nucleotide number function
def seq_count(gene_Seq):
    
    nt_Dict = {}
    
    for seq in gene_Seq:
        if seq == 'N' or seq == '\n':
            continue
        elif seq in nt_Dict:
            nt_Dict[seq] += 1
        elif seq not in nt_Dict:
            nt_Dict[seq] = 1
                    
    return nt_Dict


# Toy example
example = 'TTNNAAGGAAGG'
print('Toy example for seq:', example, seq_count(example))


# chromosome 22
handle = open("chr22.fa", "r")
handle.readline()
ch22 = handle.read()
handle.close()

ch22 = ch22.upper()
result = seq_count(ch22)

    
# nucleotide frequency calculation
bp_length = 0
ratio_sum = 0

for nt in result:
    bp_length += result[nt]

for nt in result:
    nt_ratio = result[nt]/bp_length
    ratio_sum += nt_ratio

    print(nt, result[nt], nt_ratio)

print('Sum', bp_length, ratio_sum)
