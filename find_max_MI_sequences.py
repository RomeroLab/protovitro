from sys import argv
import numpy

filename = argv[1]

# get the query sequence (the first seq)
with open(filename) as infile:
    for line in infile:
        if line[0]!='#' and len(line)>10:
            name,queryseq = line.split()
            break

# the positions where the query doesn't have a gap
querypos = [i for i in range(len(queryseq)) if queryseq[i]!='-']
queryseq = ''.join([queryseq[i] for i in querypos])




coverage_cut = 0.75 # want sequences to be at least 75%
identity_cut = 0.20

sequences = []
names = []
with open(filename) as infile:
    for line in infile:
        if len(line)>10 and line[0]!='#':
            name,seq = line.split()
            seq = ''.join([seq[i] for i in querypos])
            cov = 1 - seq.count('-')/len(seq)
            iden = len([p for p in range(len(queryseq)) if queryseq[p]==seq[p]])/(cov*len(queryseq))
            if cov>coverage_cut and iden>identity_cut:
                sequences.append(seq)
                names.append(name)



# This is an MSA.  It is represented as a list of lists
# the first dimension is MSA postions, the second dimension is over homolgous sequences
# you should do any required filtering (gaps, iden, etc) before
# I've been using jackhmmer to generate
align = tuple(zip(*sequences))



# generate solution space, omit gaps
SS = [tuple([a for a in set(p) if a!='-']) for p in align] 


# calculate the terms (positon,AA)
terms = []
for i,pos in enumerate(SS):
    for aa in pos:
        terms.append((i,aa))


# encode sequences in a binary X matrix 
# here rows correspond to sequences and columns to terms
Xall = []
for seq in sequences:
    Xall.append([1 if seq[term[0]]==term[1] else 0 for term in terms])


# Set your intital sequences.  Possibly sequences you have already characterized 
# this could be empty
indices = [0]
X = [Xall[i] for i in indices] 


# greedy algorithm to select 10 sequences
for i in range(10):
    print(i)
    H = []
    for x in Xall:
        Xnew = numpy.array(X+[x])
        # covariance matrix
        S = numpy.dot(Xnew,Xnew.T) +  1e-10*numpy.eye(Xnew.shape[0])
        # the MI between a set of sequences and the entire set is equal to the entropy (H) plus a constant
        # H is proportional to the log determinant of the covariance matrix
        LD = 2*sum(numpy.log(numpy.diag(numpy.linalg.cholesky(S)))) # this is how you calculate log(det(A))
        H.append(LD) # H is proportinal to LD

    ind = H.index(max(H))
    indices.append(ind)
    X = X+[Xall[ind]]


chosen_seqs = [sequences[s] for s in indices]
chosen_seqs = [''.join(s) for s in zip(*[p for p in zip(*chosen_seqs) if set(p)!=set('-')])] # remove all positions that are only gaps


open('output_seqs.txt','w').write('\n'.join(chosen_seqs))
