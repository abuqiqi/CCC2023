import cmath

n=1024
lim=10

ind=[]
for i in range(n):
    t = 0
    for j in range(lim):
        if ((i >> j) & 1):
            t |= (1 << (lim - j - 1))
    if i < t:
        ind.append(i)
        ind.append(t)
print(ind)
print(len(ind))
print(len(ind)*2*16)

# omg = [cmath.exp(-2j * cmath.pi * i / n) for i in range(n//2)]
# for x in omg:
#     print("{",float(x.real),",",float(x.imag),"},",end='')
