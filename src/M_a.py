from .Q import *

class InvalidBasisValueError(ValueError):
    def __init__(self):
        super("you must specify an 'a' parameter")

class M_a():
    """
        A system of numbers that correspond to vectors in the complex plane that are Q
        multiples of Q unit vectors, where a Q unit vector is a unit vector whose
        argument is 2pi/(p+a) for some Q values of r, k, p and a.

            r.e^{2\pi{i}\frac{k}{p+a}} = r.M_a(p)^k

        It can be shown that M_a(p).M_a(q) = M_a(0) for all iff p and q are complementary factors of a^2.

        In particular, for a=2, then M_2(0) = -1.

        That this means is that for all p in {1,2,4}

            M_2(p).M_2(4/p) = M_2(0) = -1

        which happens to mirror a analogous relationship in the integer domain:

            1/(p+2) + 1/(4/p+2) = 1/2

        In otherwords the same values of p such that a convex polygon of p+2 sides can successfully
        tessellate the plane.

        In sense, where p={1,2,4} are factors of 4 in the integer domain. M_2(p), M_2(q) are factors of M_2(0)
        in the complex plane.
    """

    def check_a(a):
        if a is None:
            raise InvalidBasisValueError()

    def __init__(self, p, a, k=1, r=1):
        try:
            a = Q.lift(a)
            r = Q.lift(r)
            p = Q.lift(p)
            k = Q.lift(k)

            kn=k/(p+a)
            a=Q.lift(a.p)
            p=Q.lift(kn.q)-a
            k=Q.lift(kn.p)

            if p.q != 1:
                k=k*p.q
                p=Q.lift(p.p)+a*(p.q-1)

            if k.q != 1:
                p=p*k.q+a*(k.q-1)
                k=Q.lift(k.p)

            if p+a < 0:
                k=k*-1
                p=(p+a)*-1-a

            if p.p+a.p != 0:
                k = Q.lift(k.p % (p.p+a.p))

            if k.p == 0:
                p.p = 0

            self.a=a
            self.r=r
            self.p=p
            self.k=k
        except ZeroDivisionError as e:
            print("div by zerp", k,p,a,r)
            raise e

    def from_ratio(pq, a=None):
        M_a.check_a(a)
        a = Q.lift(a)
        pq = Q.lift(pq)
        p=(Q_1-pq*a)/pq

        return M_a(p=p, a=a)

    def lift(p, a=None, k=1, r=1):
        """lifts a number of another kind into an M_a"""
        if isinstance(p, M_a):
            if a and a != p.a:
                return M_a(p=p.p+p.a-a,a=a,k=p.k*k,r=p.r*r)
            return p
        else:
            M_a.check_a(a)
            if isinstance(p, Q):
                return M_a(p, a, k=k, r=r)
            elif isinstance(p, int) or isinstance(p, float):
                return M_a(Q.lift(p), a, k=k, r=1)
            elif isinstance(p, complex):
                a = Q.lift(a)
                r, phi = cmath.polar(p)
                phi = math.fmod(phi+math.pi*2, math.pi*2)
                n=Q.lift((2*math.pi)/phi, max_depth=2)
                p=Q.lift(n.p)-a
                k=Q.lift(n.q)
                return M_a(p, a, k=k, r=r)
            else:
                raise ValueError(f"{p} cannot be lifted to an M_a")

    def conjugate(self):
        """creates an object that represents the conjugate of the receiver"""
        return M_a(p=self.p, a=self.a, k=self.k*-1, r=self.r)

    def complex(self, precision=Q.PRECISION):
        """creates a complex number that approximates the receiver"""
        kn=self.k_n().float()
        c=cmath.rect(self.r.float(), 2*math.pi*kn)
        if precision:
            c = complex(round(c.real, precision), round(c.imag, precision))
        return c

    def ratio(self):
        return self.k/(self.p+self.a)

    def __str__(self):
        k=f"{self.k}"
        if self.k.q > 1:
            k=f"({k})"
        frac=f"{self.k}/({self.p}+{self.a})"
        if not self.r == 1:
            frac=f"{frac}x{self.r}"
        if not self.a == 2:
            frac=f"{frac}@{self.a}"
        return f"{frac}"

    def __repr__(self):
        return f"(a,r,p,k)={(self.a, self.r, self.p, self.k)}={str(self)}={self.complex()}"

    def latex(self, use_negative=True):
        k=self.k
        if use_negative:
            n=self.a+self.p
            if k*2 > n:
                k=(k-n)
        out=f"e^{{2\\pi\\frac{{{k}}}{{{self.p}+{self.a}}}{{i}}}}"
        if not self.r == 1:
            out =f"{self.r}{out}"
        return f"${out}$"

    def latex_symbolic(self, use_negative=True):
        k=self.k
        if use_negative:
            n=self.a+self.p
            if k*2 > n:
                k=(k-n)
        out=f"M_{{{self.a}}}({self.p})^{{{k}}}"
        if not self.r == 1:
            out =f"{self.r}{out}"
        return f"${out}$"

    def k_n(self):
        return self.k/(self.p+self.a)

    def __eq__(self, pq):
        pq=M_a.lift(pq)
        return self.k_n() == pq.k_n()

    def __hash__(self):
        kn=self.k_n()
        return kn.p^kn.q

    def __lt__(self, pq):
        pq=M_a.lift(pq)
        return self.k_n() < pq.k_n()

    def __mul__(self, q):
        q=M_a.lift(q)
        r=self.r*q.r
        kn=self.k_n()+q.k_n()
        k=kn.p
        p=Q.lift(kn.q)-self.a
        return M_a(p=p, k=k, a=self.a, r=r)

    def __pow__(self, pow):
        pow = Q.lift(pow)
        k = self.k*pow
        return M_a(a=self.a, p=self.p, k=k, r=self.r)
