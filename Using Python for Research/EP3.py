def IsPrime(n):
    rr=True
    for j in range(2,n):
        if (n% j) == 0:
            rr=False
            break
    return rr




no=13195
Stop=IsPrime(no)

t=no//3;
while Stop==False and t>1:
    if no%t == 0 and IsPrime(t):
        Stop=True
        result=t
    t=t-1
print(result)