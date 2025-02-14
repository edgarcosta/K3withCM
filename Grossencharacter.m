import !"Ring/FldNum/grossenchar.m" : check_1modI_is_OK, make_grossencharacter, get_chi, IsCMField, ensure_field;
declare attributes FldNum: MatchGaloisWithRootsandootypes;


function MatchGaloisWithRootsandootypes(K)
  if not assigned K`MatchGaloisWithRootsandootypes then
    G, A, toA := AutomorphismGroup(K);
    assert #G eq Degree(K);
    z := K.1;
    b, cc := HasComplexConjugate(K);
    assert b;
    ccs := [g : g in G | toA(g)(z) eq cc(z)];
    assert #ccs eq 1;
    ccelt := ccs[1];

    // match InfinitePlaces with G
    Gs := [g : g in G];
    rootsK := [toA(g)(z) : g in Gs];
    // we need to make a choice
    rootsCC := [Evaluate(r, InfinitePlaces(K)[1]) : r in rootsK];
    roots_paired := [];
    G_paired := [];
    for pl in InfinitePlaces(K) do
      zCC := Evaluate(z, pl);
      _, i := Min([Abs(zCC - rCC) : rCC in rootsCC]);
      r := rootsK[i];
      g := Gs[i];
      roots_paired cat:= [r, cc(r)];
      G_paired  cat:= [g, ccelt*g];
    end for;
    assert roots_paired eq [toA(g)(K.1) : g in G_paired];
    // could have used Sym(n)![i1, i2..., in];
    Sn := SymmetricGroup(#G);
    Gsigma := AssociativeArray();
    for g in Gs do
      Gsigma[g] := Sn![Index(roots_paired, m(r)) : r in roots_paired] where m:=toA(g);
    end for;
    flat := func<x | &cat x>;
    unflat := func<x | [[x[i], x[i+1]] : i in [1..#x by 2]]>;
    permoo := func<oo, g | unflat(PermuteSequence(flat(oo), Gsigma[g]))>;
    K`MatchGaloisWithRootsandootypes := <G, A, toA, permoo>;
  end if;
  return K`MatchGaloisWithRootsandootypes;
end function;

intrinsic GaloisSignatureComplex(psi::GrossenChar : Precision:=false) -> Any
  { returns a signature for the Grossencharacter psi such that two Grossencharacters give the same L-function iff their Galois signatures match }
  // there should be a cheap way to do this algebraically
  if Precision cmpeq false then
    Precision := GetPrecision(0.5);
  end if;
  N, oo := Conductor(psi);
  K := NumberField(Order(N));
  //if H`issubgroup then H:=H`ambient; end if;
  R, RtoI := RayClassGroup(N, oo);
  // generator representatives of the ray class group
  gens := [RtoI(g) : g in Generators(R)];
  //G, A, toA := AutomorphismGroup(K);
  G, A, toA, permoo := Explode(MatchGaloisWithRootsandootypes(K));
  return <N, oo,
    {[ComplexField(Precision) | '@'(toA(g)(gen), psi : Precision:=2*Precision) : gen in gens] : g in G},
    {permoo(psi`type, g) : g in G}
    >;
end intrinsic;


intrinsic EulerFactorsIntegralityCheck(chi::GrossenChar, ps::SeqEnum[RngIntElt] : Precision:=false) -> BoolElt, RngIntElt
{ check if the Euler factors at given primes are integral }
  if Precision cmpeq false then
    Precision := GetPrecision(0.5);
  end if;
  eps := 10^(-Precision*0.5);
  for p in ps do
    ef := EulerFactor(chi, p : Integral:=false, Precision:=Precision);
    if exists(notused){ 1 : c in Coefficients(ef) | Abs(Round(Real(c)) - c) gt eps}then
      return false, p;
    end if;
  end for;
  return true, _;
end intrinsic;


intrinsic GaloisSignatureRational(psi::GrossenChar : Precision:=false, B:=100) -> Any
  { returns a signature for the Grossencharacter psi such that two Grossencharacters give the same L-function iff their Galois signatures match }
  if Precision cmpeq false then
    Precision := GetPrecision(0.5);
  end if;
  require EulerFactorsIntegralityCheck(psi, PrimesUpTo(B) : Precision:=Precision) : "The character did not pass the integrality test";

  N, oo := Conductor(psi);
  K := NumberField(Order(N));
  //if H`issubgroup then H:=H`ambient; end if;
  R, RtoI := RayClassGroup(N, oo);
  // generator representatives of the ray class group
  gens := [RtoI(g) : g in Generators(R)];
  //G, A, toA := AutomorphismGroup(K);
  G, A, toA, permoo := Explode(MatchGaloisWithRootsandootypes(K));
  if #gens gt 0 then
    b := false;
    for tries in [1..5] do
      // we don't want to redefine K
      Kiso := NumberFieldExtra(DefiningPolynomial(K) : prec:=Precision);
      sig_cc := Matrix([[ComplexField(Precision) | '@'(toA(g)(gen), psi : Precision:=2*Precision) : gen in gens] : g in G]);
      b, sig_alg := AlgebraizeMatrixExtra(sig_cc, Kiso);
      if b then break; end if;
      Precision *:=Precision;
    end for;
    require b : Sprintf("Failed to algebraize elements while working with maximal precision %o", Precision);
    sig_alg := {[K | Eltseq(x) : x in Eltseq(r)] : r in Rows(sig_alg)};
  else
    sig_alg := {[]};
  end if;
  return <N, oo,
    sig_alg,
    {permoo(psi`type, g) : g in G}
    >;
end intrinsic;



intrinsic ConductorBoundRationalHeckeCharacters(K::FldNum, Support::SetEnum) -> RngOrdIdl
  { returns the maximum level for an Hecke character chi such that its values lie in a subfield of K }
  // some unecessary assumption, but simplifies the code for the 4 fields of interest
  ntor := #TorsionSubgroup(UnitGroup(K));
  assert PrimeDivisors(ntor) eq [2];
  ZK := Integers(K);
  if Universe(Support) cmpeq Integers() then
    Support := [ pK : pK in PrimeIdealsOverPrime(K, p), p in Support ];
  end if;
  // deal with the odd part
  max_level := &*[Universe(Support) | p : p in Support | Norm(p) mod 2 eq 1];
  // deal with the even part
  even_primes := [p : p in Support | Norm(p) mod 2 eq 0 ];
  for p in even_primes do
    powerp := p;
    inv_factors := 1+Degree(p)*RamificationIndex(p); // comes from the torsion units
    // things get messy if one tries to do this with HeckeCharacterGroup/RayClassGroup when the ClassGroup is not trivial
    // making sure that the image of Units(K) is as large as necessary
    while #[1 : elt in InvariantFactors(AbelianGroup(UnitGroup(quo< ZK | powerp>))) | elt mod ntor eq 0] ne inv_factors do
      powerp *:=p;
    end while;
    // "Q", InvariantFactors(AbelianGroup(UnitGroup(quo< ZK | powerp>)));
    // "H", InvariantFactors(AbelianGroup(HeckeCharacterGroup(powerp)));
    max_level *:= powerp;
  end for;
  return max_level;
end intrinsic;



intrinsic RationalHeckeCharacters(H::GrpHecke :  UpToGalois:=true) -> SeqEnum[GrpHeckeElt]
  { returns the all the Hecke characters chi such that its values lie in a subfield of K }
  K := NumberField(Order(Modulus(H)));
  ntor := #TorsionSubgroup(UnitGroup(K));
  G, GtoH := AbelianGroup(H);

  // the subgroup such that chi image is contained in K
  S, StoG := sub<G | [(Order(g) div GCD(ntor, Order(g)))*g : g in Generators(G)]>;

  if not UpToGalois then
    res := [GtoH(StoG(g)) : g in S];
  else
    n := Exponent(AbelianGroup(S));
    coprime_to_n := [m : m in [1 .. n] | IsCoprime(m, n)];
    // this seems to be faster than trying to be smart
    Selts_up_to_galois := {{k*g : k in coprime_to_n} : g in S};
    res := [GtoH(StoG(Rep(g))) : g in Selts_up_to_galois];
  end if;
  return res;
end intrinsic;




function possible_infinity_types(ootype)
  // there must be a smarter way to do this
  ootype_to_mset := func<ootype | {* {* x : x in t *} : t in ootype *}>;
  mset := ootype_to_mset(ootype);
  possibilities := ootype;
  possibilities cat:=[[elt[2], elt[1]] : elt in ootype | elt[2] ne elt[1]];
  return { [x : x in elt] : elt in CartesianPower(possibilities, #ootype) | ootype_to_mset(elt) eq mset };
end function;



intrinsic GrossencharacterSearch(
  psis::SeqEnum[GrpHeckeElt],
  ootype::SeqEnum[SeqEnum[RngIntElt]],
  euler_factors::SeqEnum[Tup]
  : UpToGalois:=true, Primitive:=true, Jobs:=1, B:=100) -> SeqEnum[GrpHeckeElt]
{ Given a list of }
  N := Universe(psis)`Modulus;
  K := NumberField(Order(N));
  assert IsCMField(K);
  assert Degree(K) eq 2*#ootype;
  assert #Set([&+t : t in ootype]) eq 1; // "All oo-type parts must have same weight";
  ensure_field(K);

  //psis := RationalHeckeCharacters(K, Support : UpToGalois:=UpToGalois);
  ooall := possible_infinity_types(ootype);
  dirich := AssociativeArray();
  DG:=DirichletGroup(N,[]);
  for oo in ooall do
    if not check_1modI_is_OK(N, oo) then continue; end if;
    b, chi := get_chi(N, oo);
    if b then dirich[oo] := DG!chi; end if;
  end for;
  function match_psi(psi)
    gcs := [];
    for oo->chi in dirich do
      // this skips over check_1modI_is_OK and get_chi, which were done above
      gc := make_grossencharacter(chi, psi, N, &+oo[1], oo);
      if not exists(notused){ 1 : elt in euler_factors |
        EulerFactor(gc, p : Integral:=true) ne ef where p, ef := Explode(elt)} then
        Append(~gcs, gc);
      end if;
    end for;
    if Primitive then // the primitive check is more expensive than the EulerFactor check
      gcs := [g : g in gcs | IsPrimitive(g)];
    end if;
    res := [elt`type : elt in gcs];
    return res;
  end function;

  b, p := IsIntrinsic("ParallelCall");
  if b and Jobs gt 1 then
    // The function block is necessary to prevent Magma from getting confused when forking
    out := function ()
      return p(Jobs, match_psi, [<psi> : psi in psis], 1);
    end function();
    assert &and[elt[1] : elt in out];
    valid_oo := [elt[2,1] : elt in out];
  else
    valid_oo := [match_psi(psi) : psi in psis];
  end if;

  gcs := [make_grossencharacter(dirich[oo], psi, N, &+oo[1], oo) :  oo in valid_oo[i], i->psi in psis];
  if UpToGalois then
    signatures := [GaloisSignatureRational(g) : g in gcs | EulerFactorsIntegralityCheck(g, PrimesUpTo(B)) ];
    gcs := [gcs[i] : i in {Index(signatures, sig) : sig in signatures}];
  end if;
  return gcs;
end intrinsic;





/*
function possible_GC(K, Support, ootype)
  assert IsCMField(K);
  assert Degree(K) eq 2*#ootype;
  assert #Set([&+t : t in ootype]) eq 1; // "All oo-type parts must have same weight";
  ensure_field(K);
  // we want all the characters with the same modulus
  time psis := possible_hecke_characters(K, Support);
  N := Universe(psis)`Modulus;
  oo4, oo2, ooall := possible_infinity_types(K, ootype);
  dirich := AssociativeArray();
  for oo in ooall do
    if not check_1modI_is_OK(N, oo) then continue; end if;
    b, chi := get_chi(N, oo);
    if b then dirich[oo] := chi; end if;
  end for;
  oo4 meet:= Keys(dirich);
  oo2 meet:= Keys(dirich);
  ooall meet:= Keys(dirich);
  gcs := [];
  for i->psi in psis do
    oos := Order(chi) le 2 select oo2 else (Order(chi) eq 4 select oo4 else ooall);
    // this skips over check_1modI_is_OK and get_chi, which were done above
    gcs cat:= [make_grossencharacter(dirich[oo], psi, N, &+oo[1], oo) : oo in oos];
  end for;
  ps := [p : p in PrimesUpTo(1000) | not p in Support and #PrimeIdealsOverPrime(K, p) eq Degree(K)][1..10];
  return filter_integral(gcs, ps);
end function;

function match_GC(K, Support, ootype, euler_factors)
  time GCS := SequenceToSet(possible_GC(K, Support, ootype));
  for pef in euler_factors do
    p, ef := Explode(pef);
    GCS := { g : g in GCS | EulerFactor(g, p : Integral:=true) eq ef};
    p, #GCS;
  end for;
  return GCS;
end function;
*/
