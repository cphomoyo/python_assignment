#the program computes the average length of a string list

my_list=['TGAC', 'aattggcc', 'ccCCcc']
length = [len(i) for i in my_list]
print (float(sum(length)) / len(length))
