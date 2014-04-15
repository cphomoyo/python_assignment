
import read_fasta_file
import sys
import reverse_complement
import find_orfs
import get_genes_by_ORF
import translate
import get_fasta

path_file=raw_input('enter fasta file path:  ')

dna=read_fasta_file(open(sys.argv[1], 'r').readlines())

for i in dna:
	rc=reverse.complement(i)
frame=find_ORFs(dna)
#if dna does not have any ORFs then the find_ORFs function will be used to obtain ORFs on the reverse complement.
frame=find_ORFs(rc)

gene=get_gene_by_ORF(frame,0)

proteinsequence=translate(gene)

fasta_format=get_fasta(id,proteinsequence)

protein_file=raw_input('enter protein file path: ')

protein=open(protein_file,'w')

protein.write=(fasta_format)

protein.close()
		
