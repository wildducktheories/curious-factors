import math
import cmath

from math import gcd
from ordered_set import OrderedSet
from src import *

class PrecisionError(ValueError):
    def __init__(self, p):
        super().__init__(f"precision error: p={p} is not an integer")

class Q:
    PRECISION=12

    def lift(m, max_depth=10):
        if isinstance(m, Q):
            return m
        else:
            if abs(m) == 0:
                return Q(0,1)
            elif abs(m) < 1:
                r= Q.lift(1/m,max_depth-1)**-1
                return r
            else:
                try:
                    m1=round(m,Q.PRECISION)
                    return Q(m1, 1)
                except PrecisionError as e:
                    p=Q(int(round(m)),1)
                    if max_depth > 0:
                        frac=Q.lift(m-p.p, max_depth-1) #sometimes unbounded
                        return p+frac
                    else:
                        return p

    def __init__(self, p,q):

        if not int(p) == p:
            raise PrecisionError(p)

        if not int(q) == q:
            raise PrecisionError(q)

        p=int(p)
        q=int(q)

        if q < 0:
            p=p*-1
            q=q*-1

        if q == 0:
            p,q=int(0),int(0)
        else:
            p,q=int(p/gcd(p,q)),int(q/gcd(p,q))

        self.p = p
        self.q = q


    def is_undefined(self):
        return self.q == 0

    def __str__(self):
        if self.q==0:
            return "\u221e"
        elif self.q==1:
            return f"{self.p}"
        else:
            return f"{self.p}/{self.q}"

    def __repr__(self):
        return self.__str__()

    def __add__(self, pq):
        pq=Q.lift(pq)
        return Q(self.p*pq.q+pq.p*self.q, self.q*pq.q)

    def __sub__(self, pq):
        pq=Q.lift(pq)
        return self+(pq*-1)

    def __mul__(self, pq):
        pq=Q.lift(pq)
        return Q(self.p*pq.p, self.q*pq.q)

    def __pow__(self, pow):
        if pow < 0:
            return Q(self.q**-pow, self.p**-pow)
        elif pow == 0:
            return Q(0,1)
        else:
            return Q(self.p**pow, self.q**pow)

    def __truediv__(self, pq):
        pq=Q.lift(pq)
        return self * pq**-1

    def __eq__(self, pq):
        pq=Q.lift(pq)
        return self.p == pq.p and self.q == pq.q

    def __lt__(self, pq):
        pq=Q.lift(pq)
        pq = self - pq
        return pq.p < 0

    def __gt__(self, pq):
        pq=Q.lift(pq)
        pq = self - pq
        return pq.p > 0

    def __ge__(self, pq):
        return not self < pq

    def __hash__(self):
        return self.p^self.q

    def float(self):
        return self.p/self.q

Q0=Q(0,1)
Q1=Q(1,1)


