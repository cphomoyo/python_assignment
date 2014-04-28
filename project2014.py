#Sipho Moyo
#G10m3876
#Multiple Sequence Analyzer

multiple_alignment = str(None)
file_upload = False
option=0000

#Print menu that dispalys the options to be selected
print '********************************************************************************'
print '* MULTIPLE ALIGNMENT ANALYZER                                                  *'
print '********************************************************************************'
print '* Select an option from below:                                                 *'
print '*                                                                              *'

menu = {}
menu['*      (1)']=" Open a Multiple Alignment File                   (0)               *" 
menu['*      (2)']=" Alignment Information                            (A)               *"
menu['*      (3)']=" Alignment Explorer                               (E)               *"
menu['*      (4)']=" Information per Sequence                         (I)               *"
menu['*      (5)']=" Analysis of Glycosylation Signatures             (S)               *"
menu['*      (6)']=" Export to Fasta                                  (X)               *"
menu['*      (7)']=" Sequence stability analysis                      (U)               *"
menu['*      (8)']=" Exit                                             (Q)               *"
while option!='': 
  options=menu.keys()
  options.sort()
  print '********************************************************************************'
  print '* MULTIPLE ALIGNMENT ANALYZER                                                  *'
  print '********************************************************************************'
  print '* Select an option from below:                                                 *'
  print '*                                                                              *' 
  for entry in options:                                                                 
      print entry, menu[entry]
  print '*                                                 File: ' + multiple_alignment +'                  *' 
  print '********************************************************************************'  

  option=raw_input("Please select option:  ")
  #extracts the file data and stores it in memory
  if option =='1' or option=='O': 
    while file_upload!=True: 
      multiple_alignment=raw_input('enter file path: ')
      import os.path
      file_check=os.path.exists(multiple_alignment) #allows the program to check if the file exists anywhere in memory so that the program doesnt crash if the file cannot be found
      if file_check == True:#if the file exists the program will begin to extract the data
	file_upload=True
	aligned_file = open(multiple_alignment,'r')
	f=aligned_file.readline()
	if f[0]!='CLUSTAL': 
	  print 'Error! This is not a ClustalW file'
	else:
	
	  multiple_align={}
	  alignment_type =[]
	  sequences=[]
	  IDs=[] 
	  matches=[]
	  mat=[]
	  a=0
	  for line in aligned_file:
		  if line.startswith('CLUSTAL'):
			  x=5000
			  alignment_type.append(line[0:-1])
			  
		  elif line=='\n': #calculates the number of empty lines to be used further in the program
			  if a < x:
				  a+=1
		  elif line[0]!=' ': 
			  x=a
			  i=line.find(' ')
			  IDs.append(line[0:i])
			  n=line[i:-1]
			  j=len(n)
			  seq=n.strip()
			  k=len(seq)
			  l=j-k
			  sequences.append(seq)
		  else:
			  b=line.find(' ')
			  c=line[b:-1]
			  d=len(c)
			  e=c.strip()
			  f=len(e)
			  g=d-f
			  matches.append(line.strip())
			  z=''.join(matches)
	  index=0
	  aligned_file.seek(0)	
	  for name in IDs:
		  if name in multiple_align.keys():
			  multiple_align[name]=multiple_align[name]+sequences[index]
		  else:
			  multiple_align[name]=sequences[index]
		  index+=1    	
	  print 'The file', multiple_alignment, 'has been loaded successfully'
      else:
	  print 'File does not exist!'
	  multiple_alignment=raw_input('enter file path: ')
	
  elif option=='2' or option=='A':
    if file_upload!=True:
      print 'No file uploaded!'
    else:
      aligned_file = open(multiple_alignment,'r')
      multiple_align={}
      alignment_type =[]
      sequences=[]
      IDs=[] 
      matches=[]
      mat=[]
      a=0
      for line in aligned_file:
	      if line.startswith('CLUSTAL'):
		      x=5000
		      alignment_type.append(line[0:-1])
	      elif line=='\n':
		      if a < x:
			      a+=1
	      elif line[0]!=' ':
		      x=a
		      i=line.find(' ')
		      IDs.append(line[0:i])
		      n=line[i:-1]
		      j=len(n)
		      seq=n.strip()
		      k=len(seq)
		      l=j-k
		      sequences.append(seq)
	      #calculates the number of spaces between the match and line[0]
	      else:
		      b=line.find(' ')
		      c=line[b:-1]
		      d=len(c)
		      e=c.strip()
		      f=len(e)
		      g=d-f
		      matches.append(line.strip())
		      z=''.join(matches)
      index=0
      aligned_file.seek(0)	
      for name in IDs:
	      if name in multiple_align.keys():
		      multiple_align[name]=multiple_align[name]+sequences[index]
	      else:
		      multiple_align[name]=sequences[index]
	      index+=1        
      print 'Filename:     ' + multiple_alignment
      print 'Sequences:    ' + ', '.join(multiple_align.keys())
      for i in multiple_align.values():
	y = 'Length:       ' + str(len(i)) 
      print y
      for n in multiple_align.values():
	for i in n:
	  value=[]
	  if n=='A'or n=='T'or n=='C'or n=='G':
			value.append(n)
    
      align=[]
      for x in z:
	if x=='*':
	  align.append(x)
      match_length = float(len(align))
      
      q=round((match_length/len(z))*100,3)
      
      print '% of matches:         ' + str(q) + '%'
      print 'Press [enter] to display the menu again: '
      
  elif option=='3' or option=='E':
    if file_upload!=True:
      print 'No file uploaded!'
    else:    
      
      h=g-l
      multiple_align.update({h*' ':z})
      start = raw_input('Enter the start of segment:  ')
      stop = raw_input('Enter the stop of segment:   ')
      y='['+ start +'-'+ stop +']'
      print 'CLUSTAL segment' + y + 'of the', multiple_alignment, 'alignment'
      print (a/2)*'\n'
      line=''
      for n in multiple_align.keys():
	      if len(multiple_align[n][int(start):int(stop)+1])<=60:
		      print n, l*' ', multiple_align[n][int(start):int(stop)+1]
	      
      print 'Press [enter] to display the menu again: '	
		
  elif option=='4' or option=='I':
    name=raw_input('enter sequence ID:  ')
    s=[]
    for i in multiple_align:
	if name==multiple_align.keys():
	    print 'ID:   ' + name
	elif name in multiple_align:
	    multiple_align[name].replace('-','')
	    sname=len(multiple_align[name].replace('-',''))
    print 'ID:    ' + name
    print 'Length:   ' + str(sname)
    print 'Base Frequency:    '
    print '             [A]:    ' + str(multiple_align[name].count('A'))
    print '             [T]:    ' + str(multiple_align[name].count('T'))
    print '             [C]:    ' + str(multiple_align[name].count('C'))
    print '             [G]:    ' + str(multiple_align[name].count('G'))
    print 'Sequence:    ' + multiple_align[name].replace('-','')
    print 'Press [enter] to display the menu again: '
  elif option=='5' or option=='S':
    import re
    if file_upload!=True:
	  print 'No file uploaded!'
    else:
      aligned_file = open(multiple_alignment,'r')
      multiple_align={}
      alignment_type =[]
      sequences=[]
      IDs=[] 
      matches=[]
      mat=[]
      a=0
      for line in aligned_file:
	      if line.startswith('CLUSTAL'):
		      x=5000
		      alignment_type.append(line[0:-1])
	      elif line=='\n':
		      if a < x:
			      a+=1
	      elif line[0]!=' ':
		      x=a
		      i=line.find(' ')
		      IDs.append(line[0:i])
		      n=line[i:-1]
		      j=len(n)
		      seq=n.strip()
		      k=len(seq)
		      l=j-k
		      sequences.append(seq)
	      #calculates the number of spaces between the match and line[0]
	      else:
		      b=line.find(' ')
		      c=line[b:-1]
		      d=len(c)
		      e=c.strip()
		      f=len(e)
		      g=d-f
		      matches.append(line.strip())
		      z=''.join(matches)
      index=0
      aligned_file.seek(0)	
      for name in IDs:
	      if name in multiple_align.keys():
		      multiple_align[name]=multiple_align[name]+sequences[index]
	      else:
		      multiple_align[name]=sequences[index]
	      index+=1              
      prot=[]
      my_dict={}
      seq=multiple_align.values()
      for i in seq:
		i=i.replace('-','')
		x=len(i)%3 
		if x==0:
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
			                                    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
			                                    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
			                                     }
		      proteinsequence=''
		      for n in range(0,len(i),3):
			    codon = i[n:n+3]
			    if codon=='TAA' or codon=='TGA' or codon=='TAG':
				  proteinsequence=proteinsequence + codon_table[codon]
				  break	 
			    proteinsequence = proteinsequence + codon_table[codon]
			    prot.append(proteinsequence)
			    regexobj=re.compile('N[^P](S|T)+')
			    		    
						  
		elif x==1:
		      i=i[:-1]
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
			                                        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
			                                        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
			                                         }
				
		      proteinsequence = ''
		      for n in range(0,len(i),3):
			    codon = i[n:n+3]
			    if codon=='TAA' or codon=='TGA' or codon=='TAG':
				  proteinsequence=proteinsequence + codon_table[codon]
				  break	 
			    proteinsequence = proteinsequence + codon_table[codon]
			    prot.append(proteinsequence)
			    regexobj=re.compile('N[^P](S|T)+')
			    			    
		else:
		      i=i[:-2]
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
			                                              'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
			                                              'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
			                                               }
				      
		      proteinsequence = ''
		      for n in range(0,len(i),3):	
			    codon = i[n:n+3]
			    if codon=='TAA' or codon=='TGA' or codon=='TAG':
				  proteinsequence=proteinsequence + codon_table[codon]
				  break	 
			    proteinsequence = proteinsequence + codon_table[codon]
			    prot.append(proteinsequence)
			    regexobj=re.compile('N[^P](S|T)+')
			    
		print 'Press [enter] to display the menu again: '	      
  elif option=='6' or option=='X':
    if file_upload!=True:
	  print 'No file uploaded!'
    else:
      aligned_file = open(multiple_alignment,'r')
      multiple_align={}
      alignment_type =[]
      sequences=[]
      IDs=[] 
      matches=[]
      mat=[]
      a=0
      for line in aligned_file:
	      if line.startswith('CLUSTAL'):
		      x=5000
		      alignment_type.append(line[0:-1])
	      elif line=='\n':
		      if a < x:
			      a+=1
	      elif line[0]!=' ':
		      x=a
		      i=line.find(' ')
		      IDs.append(line[0:i])
		      n=line[i:-1]
		      j=len(n)
		      seq=n.strip()
		      k=len(seq)
		      l=j-k
		      sequences.append(seq)
	      else:
		      b=line.find(' ')
		      c=line[b:-1]
		      d=len(c)
		      e=c.strip()
		      f=len(e)
		      g=d-f
		      matches.append(line.strip())
		      z=''.join(matches)
      index=0
      aligned_file.seek(0)	
      for name in IDs:
	      if name in multiple_align.keys():
		      multiple_align[name]=multiple_align[name]+sequences[index]
	      else:
		      multiple_align[name]=sequences[index]
	      index+=1        
      
      file_option=raw_input('If you want to create Multiple files(one per sequence) press M or press S for a single file that includes all the sequences:   ')
      if file_option=='S':
	seq_input=raw_input('Do you want to save the DNA sequences(input D) or the proteins (input P): ')
	if seq_input=='D':
	  save_file = raw_input('Please input the name of the file to save:    ')
	  new_file = open(save_file, 'w')
	  seq_names = multiple_align.keys()
	  for i in seq_names:
	    line = '>'+ i + '\n'+ multiple_align[i]
	    new_file.write(line)
	elif seq_input=='P':
	  save_file = raw_input('Please input the name of the file to save:    ')
	  f=open(save_file, 'w')
	  for i in prot:
	    f.write('>'+i)
      elif file_option== 'M':
	seq_input=raw_input('Do you want to save the DNA sequences(input D) or the proteins (input P): ')
	if seq_input=='D':
	  for i in multiple_align.keys():
	    save_file = raw_input('Please input the name of the file to save:    ')
	    open(save_file, 'w').write('>'+ i + '\n' + multiple_align[i])
	if seq_input=='P':
	  for n in prot:
	    save_file = raw_input('Please input the name of the file to save:    ')
	    open(save_file, 'w').write('>'+ '\n' + prot[n])	    
	print 'Press [enter] to display the menu again: '
  elif option=='7' or option=='U':
    if file_upload!=True:
	      print 'No file uploaded!' 
    else:
      index=0
      count_result=[]
      for i in multiple_align.keys():
	      
	      y = multiple_align[i].count('C')
	      u = multiple_align[i].count('G')
	      w=y+u
	      count_result+=[(i,w)]
      counter_result=sorted(count_result, key=lambda a:a[1])#sorts list of list from the smallest value to the largest
      print count_result[-1][0] + ' is the most stable sequence with ' + str(count_result[-1][1]) + ' G+C residues'    
      print 'Press [enter] to display the menu again: '
  elif option=='8' or option=='Q':
    option=''
    
    
      
    
    
	  
      
    
    
    
    
    
     
