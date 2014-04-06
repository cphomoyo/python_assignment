#the program does point extension on input sequence

s=raw_input("enter sequence: ")
s=s.replace('t','U') 
s=s.replace('T','U')	
s = s[::-1]
print s
