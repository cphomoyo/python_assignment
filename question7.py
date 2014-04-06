#the program computes multiples of 5 from a fibonacci sequence

five_fib = 0 #the multiples of five 
s = 0
a = 1
b = 0
fib = []
while s <=4000:
    s = a + b
    b = a
    a = s
    if s % 5 == 0:
        five_fib += s
        fib.append(s)
print(fib)
print(five_fib)
