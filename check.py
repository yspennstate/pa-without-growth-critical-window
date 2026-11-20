from fractions import Fraction as F
from itertools import combinations
import sys
def paths(n,m,a):
    # enumerate all simple sequential paths; compute P (simple process) and Q (urn, ordered edges, unordered endpoints)
    out=[]
    def rec(E,deg,P,Qp,k):
        if k==m: out.append((P,Qp)); return
        w=[deg[v]+a for v in range(n)]; W=sum(w)
        non=[(u,v) for u,v in combinations(range(n),2) if (u,v) not in E]
        Zs=sum(w[u]*w[v] for u,v in non)
        S=sum(x*x for x in w); Qs=sum(w[u]*w[v] for u,v in E)
        b=F(S+W+2*Qs, W*(W+1))
        # identity check: sum over nonedges of urn probs = 1-b
        assert sum(F(2*w[u]*w[v],W*(W+1)) for u,v in non)==1-b
        for u,v in non:
            p=F(w[u]*w[v],Zs); q=F(2*w[u]*w[v],W*(W+1))
            assert q==p*(1-b)
            d=list(deg); d[u]+=1; d[v]+=1
            rec(E|{(u,v)},d,P*p,Qp*q,k+1)
    rec(frozenset(),[0]*n,F(1),F(1),0)
    return out
for n,m,a in [(4,3,F(1)),(5,4,F(1)),(5,3,F(1,4)),(4,4,F(4))]:
    o=paths(n,m,a); sP=sum(p for p,q in o); sQ=sum(q for p,q in o)
    tv=sum(abs(p-q/sQ) for p,q in o)/2
    print(n,m,a,"paths",len(o),"sumP",sP,"Q(A)",sQ,"TV",tv)
# NB factorial moments: E(D)_k=(a)_k/(a+1)^k
for a in [F(1,4),F(1),F(4)]:
    mu=a/(a+1); kap=(a+2)/(a+1)
    f1=mu; f2=a*(a+1)/(a+1)**2; f3=a*(a+1)*(a+2)/(a+1)**3
    ED2=f2+f1; ED3=f3+3*f2+f1; eta=mu*ED3-ED2**2
    print(a,"f2==mu",f2==mu,"f3==mu*kap",f3==mu*kap,"eta==mu^2kap",eta==mu*mu*kap)
# ODE solutions symbolic check at random points
for a in [F(1,4),F(1),F(4)]:
    c=(a+1)/a
    for t in [F(1,3),F(2)]:
        x=a+2*t; s=c*x*x-x; q=t*c*c*x*x
        ds=2*c*x*2-2; dq=c*c*x*x+t*c*c*2*x*2
        assert ds==4*s/x+2 and dq==4*q/x+(s/x+1)**2
        assert (s+x+2*q)/x**2==c+2*t*c*c
    assert c*a*a-a==a*a
print("ODE ok")
# Dirichlet-multinomial exact mean of sum d(d-1) over 2m labels, n vertices, weight a each
for a in [F(1,4),F(1),F(4)]:
    n=1000; m=n*a/(2*(a+1))
    L=2*m; ES=n*L*(L-1)*a*(a+1)/(n*a*(n*a+1))
    nu=ES/L; print(a,"nu-1 =",float(nu-1),"(order 1/n)")
