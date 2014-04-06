#this program outputs the number of non-nucleotide bases in an input sequence
s=raw_input("enter sequence:")
seq=list(s)
non_nucleotide=['b','d','e','f','h','i','j','k','l','m','n','o','p','q','r','s','u','v','w','x','y','z']
nucleotide=['a','c','t','g']
counter=0
for base in seq: 
	if base in nucleotide:
		counter=counter
	if base in non_nucleotide:
		counter=counter+1

print counter



