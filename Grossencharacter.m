import !"Ring/FldNum/grossenchar.m" : check_1modI_is_OK, make_grossencharacter, get_chi, IsCMField, ensure_field, initialise_grossenideal;
declare attributes FldNum: MatchGaloisWithRootsandootypes;
declare verbose GCSearch, 3;

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


function is_primitive_abstract_generic(CG)
  assert Type(CG) in [GrpDrchNF, GrpHecke];
  constructor := Type(CG) eq GrpDrchNF select DirichletGroup else HeckeCharacterGroup;
  N, oo := Modulus(CG);
  max_divs := [ <N/elt[1], oo> : elt in Factorization(N)] cat [ <N, [t : i->t in oo | i ne j]> : j in [1..#oo]];
  max_subgps := [Extend(constructor(elt[1], elt[2]), CG) : elt in max_divs];
  return function(g)
      return forall(notused){1 : G in max_subgps | not g in G`AbGrp};
  end function;
end function;

intrinsic IsPrimitiveAbstract(H::GrpHecke) -> UserProgram
{ return a function that replicates IsPrimitive on the abstract abelian group}
  return is_primitive_abstract_generic(H);
end intrinsic;

intrinsic IsPrimitiveAbstract(D::GrpDrchNF) -> UserProgram
{ " } //"
  return is_primitive_abstract_generic(D);
end intrinsic;


// elements of the abelian group that have order dividing O
function bounded_order_subgroup(G, O)
  return sub<G | [(Order(g) div GCD(O, Order(g)))*g : g in Generators(G)]>;
end function;

function rational_characters_abstract(C, O, UpToGalois)
  // the subgroup, as an abstract abelian group, such that chi image is contained in K
  S, StoG := bounded_order_subgroup(C, O);
  if not UpToGalois then
    res := [StoG(g) : g in S];
  else
    n := Exponent(AbelianGroup(S));
    coprime_to_n := [m : m in [1 .. n] | IsCoprime(m, n)];
    // this seems to be faster than trying to be smart
    Selts_up_to_galois := {{k*g : k in coprime_to_n} : g in S};
    res := [StoG(Rep(g)) : g in Selts_up_to_galois];
  end if;
  return res;
end function;

function rational_characters(CG, UpToGalois, Abstract)
  ntor := #TorsionUnitGroup(K) where K := NumberField(Order(Modulus(CG)));
  G, GtoH := AbelianGroup(CG);
  resG := rational_characters_abstract(G, ntor, UpToGalois);
  if Abstract then
    return resG;
  else
    return [GtoH(g) : g in resG];
  end if;
end function;




intrinsic RationalHeckeCharacters(H::GrpHecke :  UpToGalois:=false, Abstract:=false) -> SeqEnum[GrpHeckeElt]
  { returns the all the characters chi such that its values lie in a subfield of K }
  return rational_characters(H, UpToGalois, Abstract);
end intrinsic;

intrinsic RationalDirichletCharacters(D::GrpDrchNF :  UpToGalois:=false, Abstract:=false) -> SeqEnum[GrpDrchNFElt]
  { " } //"
  return rational_characters(D, UpToGalois, Abstract);
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

  ooall := possible_infinity_types(ootype);
  vprintf GCSearch, 2: "#psis = %o #ooall = %o\n", #psis, #ooall;
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
  // filter by integrality
  gcs := [g : g in gcs | EulerFactorsIntegralityCheck(g, PrimesUpTo(B))];
  vprintf GCSearch, 2: "#matches = %o\n", #gcs;
  if UpToGalois then
    signatures := [GaloisSignatureRational(g) : g in gcs ];
    gcs := [gcs[i] : i in {Index(signatures, sig) : sig in signatures}];
    vprintf GCSearch, 2: "#matches up to Galois = %o\n", #gcs;
  end if;
  return gcs;
end intrinsic;


intrinsic CandidateRationalGrossenCharacters(
  N::RngOrdIdl,
  ootype::SeqEnum[SeqEnum[RngIntElt]]: Primitive:=true
  ) -> SeqEnum[GrossenChar]
{ Return all the Grossencharacters such that Q(Restriction(gr`hecke)/gr`dirich) < K }
  K := NumberField(Order(N));
  ensure_field(K);
  // call internal functions to avoid try/catch
  if not check_1modI_is_OK(N, ootype) then return []; end if;
  b, chi0 := get_chi(N, ootype);
  if not b then return []; end if;
  HG := HeckeCharacterGroup(N, []);
  HGAb, toHG := AbelianGroup(HG);
  DG := DirichletGroup(N, []);
  DGAb, toDG := AbelianGroup(DG);

  toHGAb := hom<HG -> HGAb | x:->HGAb!Eltseq(x)>;
  toDGAb := hom<DG -> DGAb | x:->DGAb!Eltseq(x)>;
  // chi0 by default is in the subgroup of DG generated by itself
  chi0_ab := DGAb ! Eltseq(DG!chi0);
  UGAb := AbelianGroup(UnitTrivialSubgroup(DG));
  // the line below is equivalent to Grossencharacter(HG.0, ootype)
  // but this skips over check_1modI_is_OK and get_chi, which were done above
  // and creates the principal quasi-character
  psi_princ := make_grossencharacter(chi0, HG.0, N, &+ootype[1], ootype);

  // compute the Dirichlet characters chi/chi0 such that Q(chi/chi0) < K
  ntor := #TorsionUnitGroup(K);
  // bounded_order_subgroup gives us the subgroup of characters {chi' = chi/chi0} such that Order(chi') | ntor in DGAb
  // note that we avoid RationalDirichletCharacters or rational_characters, as we want to work with the abstract group
  S, StoDGAb := bounded_order_subgroup(DGAb, ntor);
  // we further want:
  // - chi to be liftable to a Hecke character, i.e., chi is trivial on the units
  // - (If Primitive) chi' = chi/chi0 to be primitive, as this is equivalent that the quasi-character is also primitive

  if Primitive then
    N, oo := Modulus(DG);
    max_divs := [ <N/elt[1], oo> : elt in Factorization(N)] cat [ <N, [t : i->t in oo | i ne j]> : j in [1..#oo]];
    max_subgps := [Extend(DirichletGroup(elt[1], elt[2]), DG)`AbGrp : elt in max_divs];
  else
    max_subgps := [];
  end if;

  // we want to take the interesection
  // S meet &meet [Complement(G) : G in max_subgps] meet chi0_ab^{-1} U
  // we first try to reduce S by interesecting it with the sugroup generated by chi0_ab^{-1} U
  Sprime := sub<DGAb | [-chi0_ab] cat [u-chi0_ab : u in Generators(UGAb)]> meet S;

  // we now could loop over Sprime, but our answer will come in cosets
  // UGAb meet &meet [Gi : Gi in max_subgps] meet Sprime
  // thus we quotient everything by it
  B := UGAb meet S;
  if not IsEmpty(max_subgps) then
    B meet:= &meet [G : G in max_subgps];
  end if;
  assert B subset Sprime;
  max_subgps_Sprime := [G meet Sprime : G in max_subgps];
  is_primitive := func<g | forall(notused){1 : G in max_subgps_Sprime | not g in G}>;

  chis_cosets := [chi0_ab + StoDGAb(g) : g in Transversal(Sprime, B) | is_primitive(g) and chi0_ab + StoDGAb(g) in UGAb];
  chis := [g + b : b in B, g in chis_cosets];


  // we now lift the chis, by looping over the preimages of the restriction
  restriction := hom<HGAb -> DGAb | [<toHGAb(g), toDGAb(DG!DirichletRestriction(g))> : g in Generators(HG)]>;
  HCelts := [HGAb | g : g in AbelianGroup(HilbertCharacterSubgroup(HG)) ];
  // equivalent to psi0s := &cat[ [HG | psi*toHG(elt) :  elt in HCelts] where psi := HeckeLift(toDG(chi)) : chi in chis];
  psi0s := [HG | toHG(psi + elt) :  elt in HCelts, psi in chis @@ restriction];
  // we now twist by the Hecke characters
  return [psi_princ*psi0 : psi0 in psi0s ];
end intrinsic;


intrinsic GrossencharacterSearch(
  N::RngOrdIdl,
  ootype::SeqEnum[SeqEnum[RngIntElt]],
  euler_factors::SeqEnum[Tup]
  :
  UpToGalois:=true,
  Primitive:=true,
  Jobs:=1,
  B:=100,
  Permuteoo:=true
  ) -> SeqEnum[GrpHeckeElt]
{ Given a list of }
  vprintf GCSearch, 3: "%o, %o, %o\n", LMFDBLabel(N), StripWhiteSpace(Sprint(ootype)), euler_factors;
  vprintf GCSearch, 3: "%o, %o, %o %o\n", UpToGalois, Primitive, Jobs, Permuteoo;
  ensure_field(NumberField(Order(N)));
  has_ParallelCall, ParallelCall := IsIntrinsic("ParallelCall");
  ooall := Permuteoo select possible_infinity_types(ootype) else [ootype];
  Ns := Primitive select [N] else Divisors(N);


  function match_psi(psi)
    return not exists(notused){ 1 : elt in euler_factors |
        EulerFactor(psi, p : Integral:=true) ne ef where p, ef := Explode(elt)};
  end function;


  if has_ParallelCall and Jobs gt 1 then
    // fill up cache, so we can make consistent choices, and thus serialize the elements
    for I in Ns do
      _ := DirichletGroup(I);
      _ := HeckeCharacterGroup(I);
      _ := initialise_grossenideal(I);
    end for;

    norms := [-Norm(elt) : elt in Ns];
    ParallelSort(~norms, ~Ns);
    nopermute := func<N, oo | [<Eltseq(Parent(psi)`ambient!psi) : psi in [* g`hecke, g`dirich *]> cat <g`type> : g in GrossencharacterSearch(N, oo, euler_factors : UpToGalois:=UpToGalois, Jobs:=1, Primitive:=true, B:=B, Permuteoo:=false)]>;
    inputs := [<N, oo> : oo in ooall, N in Ns];
    // The function block is necessary to prevent Magma from getting confused when forking
    out := function()
      return ParallelCall(Jobs, nopermute, inputs, 1);
    end function();
    assert &and[elt[1] : elt in out];
    gcs := &cat [[Grossencharacter(HeckeCharacterGroup(inputs[i,1])!x[1], DirichletGroup(inputs[i,1])!x[2], x[3]) : x in elt[2,1]] : i->elt in out];
    assert &and [match_psi(psi) : psi in gcs];

    vprintf GCSearch, 2: "global #matches = %o\n", #gcs;
    if UpToGalois then
      signatures := [GaloisSignatureRational(g) : g in gcs ];
      gcs := [gcs[i] : i in {Index(signatures, sig) : sig in signatures}];
      vprintf GCSearch, 2: "global #matches up to Galois = %o\n", #gcs;
    end if;
    return gcs;
  end if;


  vprintf GCSearch, 2: "Assembling candidates list...";
  vtime GCSearch, 2:
  psis := &cat[ CandidateRationalGrossenCharacters(N, oo : Primitive:=true) : oo in ooall, N in Ns];
  vprintf GCSearch, 2: "#psis = %o\n", #psis;


  // filter by match and integrality
  gcs := [psi : i->psi in psis | match_psi(psi) and EulerFactorsIntegralityCheck(psi, PrimesUpTo(B))];
  vprintf GCSearch, 2: "#matches = %o\n", #gcs;
  if UpToGalois then
    signatures := [GaloisSignatureRational(g) : g in gcs ];
    gcs := [gcs[i] : i in {Index(signatures, sig) : sig in signatures}];
    vprintf GCSearch, 2: "#matches up to Galois = %o\n", #gcs;
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
