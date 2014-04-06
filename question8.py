#the program computes prime numbers

n=input("enter number: ")
# 0 and 1 are not primes
if n<2:
	print n,"is not a prime number"
# 2 is the only even prime number
if n == 2:
	print "2 is a prime number"
elif n%2==0 & n%3==0 & n%5==0:
	print n,"is a prime number"
else:
	print n, "is not a prime number"



