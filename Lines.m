intrinsic FactorIntoLinear(f::RngMPolElt) -> SeqEnum[RngMPolElt]
  {return the lines defined by f}
  res := Factorisation(f);
  if #res eq Degree(f) then // factors completely over ground field
    return res;
  end if;
  Pf := Parent(f);
  require Rank(Pf) eq 3 : "Argument 1 must be over a polynomial ring with rank = 3";
  P2 := ProjectiveSpace(Pf);
  OL := BaseRing(Pf);
  L := OL;
  if IsIsomorphic(OL, RationalsAsNumberField()) then
    L := RationalsAsNumberField();
  end if;
  for elt in res do
      g := elt[1];
      require elt[2] eq 1 :  "f must be square free";
      K, _ := FieldOfGeometricIrreducibility(Curve(P2, [g] ));
      L := CompositeFields(L, K)[1];
  end for;
  res := [elt[1] : elt in Factorisation(ChangeRing(f, L))];
  require #res eq Degree(f): "f must be a product of linear factors over the algebraic closure";
  Sort(~res, func<x, y | Degree(sub<L | Coefficients(x)>) - Degree(sub<L | Coefficients(y)>)>);
  return res;
end intrinsic;

intrinsic BadPrimesOfDegree2K3Surface(f::RngMPolElt) -> SeqEnum[RngInt]
 {return the primes for which the surface w^2 = f(x,y,z) has more than ADE singularities}
  R := Generic(Parent(f));
  // maybe this would all work over a Dedekind domain
  require
    (BaseRing(R) cmpeq Integers()) and
    (Rank(R) eq 3) and
    (VariableWeights(R) eq [1,1,1]) and
    IsHomogeneous(f) and
    (TotalDegree(f) eq 6)
    :
    "Argument 1 must be a homogeneous degree 6 polynomial in a uniweighted 3-variable polynomial ring over the integers";
  K := FieldOfFractions(BaseRing(R));
  P2 := ProjectiveSpace(R);
  require Degree(ChangeRing(SingularSubscheme(Scheme(P2, f)), K)) eq 15 : "Argument 1 must be a product of linear factors over the algebraic closure";
  R2 := PolynomialRing(Integers(), 2);
  res := {};
  for v in [1..3] do
    fv := hom<R -> R2 | [i eq v select 1 else R2.(i mod 2 + 1) : i in [1..3]]>(f);
    ss := SingularSubscheme(Scheme(AffineSpace(R2), fv));
    SS := SingularSubscheme(ss);
    GB := GroebnerBasis(Ideal(SS));
    res join:=SequenceToSet(PrimeDivisors(Integers()!GB[#GB]));
  end for;
  return res;
end intrinsic;



intrinsic HyperplaneToPoint(h::Sch) -> Pt
  {}
  P := AmbientSpace(h);
  require Degree(h) eq 1 : "h must be of degree 1.";
  f1 := DefiningEquation(h);
  R := Parent(f1);
  return P![Coefficient(f1, R.i, 1) : i in [1..Rank(R)]];
end intrinsic;

intrinsic PointToHyperplane(pt::Pt) -> Sch
  {}
  P := Scheme(pt);
  if not IsAmbient(P) then
    P := Ambient(P);
  end if;
  R := CoordinateRing(P);
  c := Coordinates(pt);
  return Scheme(P, &+[c[i]*R.i : i in [1..Rank(R)]]);
end intrinsic;

intrinsic MoveLines(original::SeqEnum[Sch], target::SeqEnum[Sch]) -> MatAlgElt
  {return a linear map such that original are mapped to target}
  require #original eq #target: "target and original should have the same length";
  originalpts := [HyperplaneToPoint(p) : p in original];
  targetpts := [HyperplaneToPoint(p) : p in target];
  A := Matrix([Coordinates(p) : p in originalpts]);
  B := Matrix([Coordinates(p) : p in targetpts]);
  k := BaseRing(A);
  V := VectorSpace(k, NumberOfColumns(A));
  require Rank(A) eq #original: "original is not generic";
  require Rank(B) eq #original: "target is not generic";
  extraA := [Eltseq(v) : v in ExtendBasis(RowSpace(A), V)[#originalpts + 1..NumberOfColumns(A)]];
  A := Matrix([Coordinates(p) : p in originalpts] cat extraA);
  extraB := [Eltseq(v) : v in ExtendBasis(RowSpace(B), V)[#targetpts + 1..NumberOfColumns(B)]];
  B := Matrix([Coordinates(p) : p in targetpts] cat extraB);

  P := AmbientSpace(original[1]);
  R := CoordinateRing(P);
  gens := Matrix(Rank(R), 1, [R.i : i in [1..Rank(R)]]);
  m := map<P->P| Eltseq(ChangeRing(B^-1*A, R)*gens)>;
  assert [m(elt) : elt in original] eq target;
  return m;
end intrinsic;


/*
//TODO we could perhaps handle l1*l2, l3*l4*l5*l6
intrinsic LinesToEllipticFibration(l1::Sch, l2::Sch, f4::Sch) -> Any //FIXME
  {}
  P2 := AmbientSpace(l1);
  k := BaseRing(P2);
  e1, p := IsFinite(k);
  inputs := [l1,l2,f4];
  require #{AmbientSpace(elt) : elt in inputs} eq 1 : "the 3 schemes must be in the same projective plane";
  require Dimension(P2) eq 2 : "Branch locus must be in P^2.";
  require IsField(k): "Base ring must be a field";
  require [Degree(elt) : elt in inputs] eq [1, 1, 4] : "Branch locus must be of degree 6 = 1 + 1 + 4.";
  require not e1 or p gt 2 : "Characteristic 2 case is not covered";
  require &and[IsReduced(elt) : elt in inputs]: "Branch locus must be reduced";

  R3<x,y,z> := CoordinateRing(P2);
  // move l1 -> x and l2 ->z
  m := MoveLines(inputs[1..2], [Scheme(P2, x), Scheme(P2, z)]);
  inputs := [ m(elt) : elt in inputs];
  l1, l2, f4 := Explode(inputs);
  assert DefiningEquation(l1) eq x;
  assert DefiningEquation(l2) eq z;

  // We identify k(P1) with K<t>
  K<t> := FunctionField(k);
  P1 := ProjectiveSpace(k, 1);
  R<X> := PolynomialRing(K);

  // Now construct the morphism defining the elliptic fibration.
  phi := map<P2->P1 | [x,z]>;
  // Don't forget the twist to account for the lines x, z.
  F := t * DefiningEquation(f4);

  C := HyperellipticCurve(F);
    O := C ! [-Coefficient(branch[1],0)/Coefficient(branch[1],1), 0];
    E := EllipticCurve(C, O);
  //HERE
  return inputs;
end intrinsic;
*/

