AttachSpec("spec");
R3<x,y,z> := PolynomialRing(Integers(), 3);
f := x*y*z*(7* x^3 - 7 *x^2* y + 49 *x^2* z - 21* x* y* z + 98 *x* z^2 + y^3 - 7 *y^2* z + 49 *z^3);


// http://www.lmfdb.org/NumberField/6.0.153664.1
_<x> := PolynomialRing(Rationals());
K<a> := NumberField(x^6 + 5*x^4 + 6*x^2 + 1);
F := Subfields(K,3)[1,1];
degree2cmsubfields := [L[1] : L in Subfields(K, 2) | Signature(L[1]) eq 0];
assert #degree2cmsubfields eq 1;
Qi<I> := Polredabs(degree2cmsubfields[1]);
// we are assuming L = Q(i)
assert I^2 eq -1;


ZK := Integers(K);

p2 := Factorisation(2*ZK)[1][1];
assert 1 eq #Factorisation(2*ZK);
p7 := Factorisation(7*ZK)[1][1];
assert 1 eq #Factorisation(7*ZK);

primes := [p :  p in PrimesUpTo(1000) | not p in {2, 7}];

// Euler Factors ok K3
HX := HeckeCharacterGroup(p2 * p2 * p7);
GCX := Grossencharacter(HX.3^4,[[0,2],[1,1],[1,1]]);
LX := LSeries(GCX);

// check the first few primes
time assert &and[EulerFactor(LX, p : Integral:=true) eq Reverse(WeilPolynomialOfDegree2K3Surface(f, p) div (x - p)) : p in primes[2..9]];


// other weight 2 factors

// Known factors
GC2 := Grossencharacter(HeckeCharacterGroup(1*ZK)!1, [[0,2], [0,2], [1,1]]);
LF := TateTwist(LSeries(F), -1);
LK := LSeries(GC2);

// the 3 fold
C := HyperellipticCurve(Polynomial(Reverse([1,0,7,0,14,0,7,0])));
HA := HeckeCharacterGroup(p2^3*p7);
exps := [1,1,4];
GCA := Grossencharacter( &*[HA | HA.i^e : i->e in exps], [[0,1], [0,1], [0,1]]);
LA := LSeries(GCA);
// check the first primes
time assert &and[EulerFactor(LA, p : Integral:=true) eq LPolynomial(ChangeRing(C, GF(p))) : p in primes[2..20]];
// Wedge square
LWedgeA := [<p, AlternatingSquareCharacteristicPolynomial(Lp)> where Lp := EulerFactor(LA, p : Integral:=true) : p in primes];


// check they match
LWedgeA eq [<p, - &*[EulerFactor(L, p : Integral:=true) : L in [LX, LF, LK]]> : p in primes];
