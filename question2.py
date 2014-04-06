#the program outputs the first 3 and last 3 bases of the input sequence

s=raw_input("enter sequence: ")
if len(s)<6:
	print 'error'
elif len(s)>=6:
	s= s[:3] + s[-3:]
	print s 


