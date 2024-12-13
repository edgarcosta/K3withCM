AttachSpec("spec");
R3<x,y,z> := PolynomialRing(Integers(), 3);
f := x* y *z *(x^3 - 3 *x^2* z - 3 *x *y^2 - 3 *x *y *z + y^3 + 9* y^2 *z + 6 *y *z^2 + z^3);
_<x> := PolynomialRing(Rationals());
// http://www.lmfdb.org/NumberField/6.0.419904.1
K<a> := NumberField(x^6 + 6*x^4 + 9*x^2 + 1);
F := Subfields(K,3)[1,1];
degree2cmsubfields := [L[1] : L in Subfields(K, 2) | Signature(L[1]) eq 0];
assert #degree2cmsubfields eq 1;
Qi<I> := Polredabs(degree2cmsubfields[1]);
// we are assuming L = Q(i)
assert I^2 eq -1;

ZK := Integers(K);
p2 := Factorisation(2*ZK)[1][1];
assert 1 eq #Factorisation(2*ZK);
p3 := Factorisation(3*ZK)[1][1];


HX := HeckeCharacterGroup(p2*p2);
GCX := Grossencharacter(HX ! 1,[[2,0],[1,1],[1,1]]);

LX := LSeries(GCX);

primes := [p :  p in PrimesUpTo(1000) | not p in {2, 3}];
// check the first few primes
time assert &and[EulerFactor(LX, p : Integral:=true) eq Reverse(WeilPolynomialOfDegree2K3Surface(f, p) div (x - p)) : p in primes[2..9]];

// other weight 2 factors

// Known factors
GC2 := Grossencharacter(HeckeCharacterGroup(1*ZK)!1, [[2,0], [0,2], [1,1]]);
LF := TateTwist(LSeries(F), -1);
LK := LSeries(GC2);

// the 3 fold
C := HyperellipticCurve(Polynomial(Reverse([1,0,6,0,9,0,1,0])));
HA := HeckeCharacterGroup(p2^4);
exps := [1,0,2];
GCA := Grossencharacter( &*[HA | HA.i^e : i->e in exps], [[0,1], [0,1], [0,1]]);
LA := LSeries(GCA);
// check the first primes
time assert &and[EulerFactor(LA, p : Integral:=true) eq LPolynomial(ChangeRing(C, GF(p))) : p in primes[2..30]];


// Wedge square
LWedgeA := [<p, AlternatingSquareCharacteristicPolynomial(Lp)> where Lp := EulerFactor(LA, p : Integral:=true) : p in primes];


// check they match
LWedgeA eq [<p, - &*[EulerFactor(L, p : Integral:=true) : L in [LX, LF, LK]]> : p in primes];
