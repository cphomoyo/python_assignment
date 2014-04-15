def read_fasta_file(filename):

	fastafile=open(filename, 'r').readlines()
	fasta=[]
	index=0
	sequence = ''
	
	for line in fastafile:
		if line.startswith ('>'):
			if index>=1:
				fasta=fasta+line[1:]
			index =index+1
			name=line[:-1]
			s=fasta(name,sequence)
		else:
			sequence=sequence+line[:-1]
			s=fasta(name, sequence)
		
	return s
print read_fasta_file('Home/sipho/genomics/Elf1-500.fa')


def reverse_complement(s):
	s=s.upper()
	s = s[::-1]
	rc=[]
	for i in s:
		if i=='A':
			rc.append('T')
		elif i=='T':
			rc.append('A')
		elif i=='C':
			rc.append('G')
		elif i=='G':
			rc.append('C')
	return ''.join(rc)	
print reverse_complement('ATTATGCCCTTGCG')


def find_orfs(sequence):
   orf=''
   for i in range(0,len(sequence),3):
      codon=sequence[i:i+3]
      if codon == 'ATG':
            orf=codon 
      if codon=='TAA' or codon=='TGA' or codon=='TAG':
            orf=orf+codon
            break
      return 0
   else:
      return -1
   
         
   for i in range(1,len(sequence),3):
         codon=sequence[i:i+3]
         if codon == 'ATG':
               orf=codon 
         if codon=='TAA' or codon=='TGA' or codon=='TAG':
               orf=orf+codon
               break
         return 1
   else:
         -1
   for i in range(2,len(sequence),3):
            codon=sequence[i:i+3]
            if codon == 'ATG':
                  orf=codon 
            if codon=='TAA' or codon=='TGA' or codon=='TAG':
                  orf=orf+codon
                  break
            return 2
   else:
         return -1

      
print find_orfs('ATTATGCCCTAGCG')

def get_genes_by_ORF(sequence, ORF):
    gene=''
    result=''
    #To find the first start codon within the frame
    for i in range (0,len(sequence),3):
        codon=sequence[i:i+3]
        if codon == 'ATG':
            result=codon
        if codon=='CCC' or codon=='ATA' or codon=='ATC' or codon=='ATT' or codon=='ACA' or codon=='ACC' or codon=='ACG' or codon=='ACT' or codon=='AAC' or codon=='AAT' or codon=='AAA' or codon=='AAG' or codon=='AGC' or codon=='AGT' or codon=='AGA' or codon=='AGG' or codon=='CTA' or codon=='CTC' or codon=='CTG' or codon=='CTT' or codon=='CCA' or codon=='CCC'or codon=='CCG' or codon=='CCT' or codon=='CAC' or codon=='CAT' or codon=='CAA' or codon=='CAG' or codon=='CGA' or codon=='CGC' or codon=='CGG'or codon=='CGT' or codon=='GTA' or codon=='GTC' or codon=='GTG' or codon=='GTT'or codon=='GCA' or codon=='GCC' or codon=='GCG' or codon=='GCT' or codon=='GAC' or codon=='GAT' or codon=='GAA' or codon=='GAG' or codon=='GGA' or codon=='GGC'or codon=='GGG'or codon=='GGT' or codon=='TCA' or codon=='TCC' or codon=='TCG' or codon=='TCT'or codon=='TTC' or codon=='TTT' or codon=='TTA'or codon=='TTG' or codon=='TAC' or codon=='TAT' or codon=='TGC' or codon=='TGT' or codon=='TGG':
            result=result+codon
        if codon=='TAA' or codon=='TGA' or codon=='TAG':
            result=result+codon
            break
    
    return''.join(result)
print get_genes_by_ORF('TGCATGCCCTAGTGC',0)

def translate(sequence):

	sequence=sequence.upper()

	codon_table={
	    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	    'TAC':'Y', 'TAT':'Y', 'TAA':'stop', 'TAG':'stop',
	    'TGC':'C', 'TGT':'C', 'TGA':'stop', 'TGG':'W',
	     }

	proteinsequence = ''
	for i in range(0,len(sequence),3):
		codon=sequence[i:i+3]
		if codon=='TAA' or codon=='TGA' or codon=='TAG':
			proteinsequence=proteinsequence + codon_table[codon]
			break
		proteinsequence = proteinsequence + codon_table[codon]
	

	return proteinsequence
print translate('ATGCCCTAGTGC')

def get_fasta(id, sequence):
    s=sequence.split()
    line=''
    for i in s:
        if len(line) + len(i)<=60:
            line=line + '%s' %i
        else:
            line = '%s' %i
    
    fasta_format='>'+id+'\n'+line+'\n'
    return fasta_format
print get_fasta('protein_seq001', 'MPstop')



    

