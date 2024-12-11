AttachSpec("spec");
R3<x,y,z> := PolynomialRing(Integers(), 3);
f := x *y* z* (49*x^3 - 304*x^2*y + 570*x^2*z + 361*x*y^2 - 2793*x*y*z + 2033*x*z^2 + 361*y^3 + 2888*y^2*z - 5415*y*z^2 + 2299*z^3);

//https://www.lmfdb.org/NumberField/6.0.8340544.1
_<x,y,z> := PolynomialRing(Rationals(), 3);
//f3 := x *y* z* (49*x^3 - 304*x^2*y + 570*x^2*z + 361*x*y^2 - 2793*x*y*z+ 2033*x*z^2 + 361*y^3 + 2888*y^2*z - 5415*y*z^2 + 2299*z^3);
// the lines factor over the cubic field
_<x> := PolynomialRing(Rationals());
K := NumberField(x^6 + 13*x^4 + 50*x^2 + 49);
F := Subfields(K,3)[1,1];
degree2cmsubfields := [L[1] : L in Subfields(K, 2) | Signature(L[1]) eq 0];
assert #degree2cmsubfields eq 1;
Qi<I> := Polredabs(degree2cmsubfields[1]);
// we are assuming L = Q(i)
assert I^2 eq -1;


ZK := Integers(K);

p2 := Factorisation(2*ZK)[1][1];
p19 := Factorisation(19*ZK)[1][1];


// Euler Factors of K3
HX := HeckeCharacterGroup(p2^2 * p19);
GCX := Grossencharacter(HX.3^30,[[0,2],[1,1],[1,1]]);
primes := [p :  p in PrimesUpTo(1000) | not p in {2, 7, 11, 19}];
LX := LSeries(GCX);
// check the first few primes
time assert &and[EulerFactor(LX, p : Integral:=true) eq Reverse(WeilPolynomialOfDegree2K3Surface(f, p) div (x - p)) : p in primes[2..9]];

// other weight 2 factors

// Known factors
GC2 := Grossencharacter(HeckeCharacterGroup(1*ZK)!1, [[2,0], [0,2], [1,1]]);
LF := TateTwist(LSeries(F), -1);
LK := LSeries(GC2);

// the 3 fold
C := HyperellipticCurve(Polynomial(Reverse([1,0,1786,0,44441,0,278179,0])));
HA := HeckeCharacterGroup(p2^3*p19);
GCA := Grossencharacter(HA.1*HA.3^90, [[0,1], [0,1], [0,1]]);
LA := LSeries(GCA);
// check the first primes
assert &and[EulerFactor(LA, p : Integral:=true) eq LPolynomial(ChangeRing(C, GF(p))) : p in primes[2..20]];

// Wedge square
LWedgeA := [<p, AlternatingSquareCharacteristicPolynomial(Lp)> where Lp := EulerFactor(LA, p : Integral:=true) : p in primes];


// check they match
LWedgeA eq [<p, - &*[EulerFactor(L, p : Integral:=true) : L in [LX, LF, LK]]> : p in primes];

