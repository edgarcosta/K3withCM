import !"Ring/FldNum/grossenchar.m" : check_1modI_is_OK, make_grossencharacter, get_chi, IsCMField, ensure_field;


//TODO: many of these should be converted to intrinsics

function StringToReal(s) // (s::MonStgElt) -> RngIntElt
//{ Converts a decimal string (like 123.456 or 1.23456e40 or 1.23456e-10) to a real number at default precision. }
  if #s eq 0 then return 0.0; end if;
  if "e" in s then
    t := Split(s,"e");
    //require #t eq 2: "Input should have the form 123.456e20 or 1.23456e-10";
      return $$(t[1])*10.0^StringToInteger(t[2]);
  end if;
  t := Split(s,".");
  //require #t le 2: "Input should have the form 123 or 123.456 or 1.23456e-10";
  n := StringToInteger(t[1]);  s := t[1][1] eq "-" select -1 else 1;
  return #t eq 1 select RealField()!n else RealField()!n + s*RealField()!StringToInteger(t[2])/10^#t[2];
end function;

function bound_conductor(K, conductor_support)
    // some unecessary assumption, but simplifies the code for the 4 fields of interest
    ntor := #TorsionSubgroup(UnitGroup(K));
    assert PrimeDivisors(ntor) eq [2];
    ZK := Integers(K);
    if Universe(conductor_support) cmpeq Integers() then
        conductor_support := [ pK : pK in PrimeIdealsOverPrime(K, p), p in conductor_support ];
    end if;
    // deal with the odd part
    max_level := &*[Universe(conductor_support) | p : p in conductor_support | Norm(p) mod 2 eq 1];
    inv_factors := Degree(K) - (Degree(K) div 2) - 1; // 4 = 6 - 2 for the 4 sextic CM fields
    // deal with the even part
    even_primes := [p : p in conductor_support | Norm(p) mod 2 eq 0 ];
    for p in even_primes do
        powerp := 1*ZK;
        // increase the power of two until we get enough generators
        while #InvariantFactors(AbelianGroup(HeckeCharacterGroup(powerp))) ne inv_factors do
            powerp *:=p;
        end while;
        // increase the power of two until we get the torsion part
        while #[1 : elt in InvariantFactors(AbelianGroup(HeckeCharacterGroup(powerp))) | elt mod ntor eq 0] ne inv_factors do
            powerp *:=p;
        end while;
        max_level *:= powerp;
    end for;
    return max_level;
end function;

function possible_hecke_characters(K, conductor_support)
    N := bound_conductor(K, conductor_support);
    H := HeckeCharacterGroup(N);
    ntor := #TorsionSubgroup(UnitGroup(K));
    chis := { chi : chi in Elements(H) | ntor mod Order(chi) eq 0 };

    // exclude Galois conjugates
    n := Exponent(AbelianGroup(H));
    coprime_to_n := [m : m in [1 .. n] | IsCoprime(m, n)];
    chis_up_to_galois := {};
    while #chis gt 0 do
      chi := Rep(chis);
      Include(~chis_up_to_galois, chi);
      for m in coprime_to_n do
        Exclude(~chis, chi^m);
      end for;
    end while;
    return chis_up_to_galois;
end function;

function possible_hecke_characters(K, conductor_support)
    N := bound_conductor(K, conductor_support);
    H := HeckeCharacterGroup(N);
    ntor := #TorsionSubgroup(UnitGroup(K));
    G, GtoH := AbelianGroup(H);

    // the subgroup such that chi image is contained in K
    S, StoG := sub<G | [(Order(g) div GCD(ntor, Order(g)))*g : g in Generators(G)]>;

    //time Selts := {g : g in S};;
    n := Exponent(AbelianGroup(S));
    coprime_to_n := [m : m in [1 .. n] | IsCoprime(m, n)];
    // this is faster than trying to be smart
    Selts_up_to_galois := {{k*g : k in coprime_to_n} : g in S};
    return [GtoH(StoG(Rep(g))) : g in Selts_up_to_galois];
end function;

function filter_integral(chis, ps)
    res := [];
    for chi in chis do
        for p in ps do
            efz := EulerFactor(chi, p : Integral:=true);
            ef := EulerFactor(chi, p : Integral:=false, Precision:=500);
            if Max([Abs(c) : c in Coefficients(efz - ef)]) gt 10^-100 then
                continue chi;
            end if;
        end for;
        Append(~res, chi);
    end for;
    return res;
end function;



function galois_conjugate_ootypes(K)
    // Getting the galois group acting on CM places
    // there is way to do this the Magma way, I don't know it
    G, aut, fromG := AutomorphismGroup(K);
    z := K.1;
    b, cc := HasComplexConjugate(K);
    assert b;
    ccs := [g : g in G | fromG(g)(z) eq cc(z)];
    assert #ccs eq 1;
    ccelt := ccs[1];

    // match InfinitePlaces with G
    rootsK := [fromG(g)(z) : g in G];
    // we need to make a choice
    rootsCC := [Evaluate(r, InfinitePlaces(K)[1]) : r in rootsK];
    roots_paired := [];

    for pl in InfinitePlaces(K) do
        zCC := Evaluate(z, pl);
        _, i := Min([Abs(zCC - rCC) : rCC in rootsCC]);
        r := rootsK[i];
        roots_paired cat:= [r, cc(r)];
    end for;
    assert SequenceToSet(roots_paired) eq SequenceToSet(rootsK);

    Is := [elt[1] : elt in Roots(Polynomial([1, 0, 1]), K)];
    I := #Is eq 2 select Is[1] else K!1; // we don't need to worry about characters of order 4 otherwise
    GRR := sub<G | [g : g in G | fromG(g)(I) eq I]>;



    // could have used Sym(n)![i1, i2..., in];
    // GRR only looks at the galois action in the real places, i.e., Order(chi) eq 4
    // GCC = G, but acting on the 00-types, and what is needed if Order(chi) le 2
    GRR, GCC := Explode([
        [[Index(roots_paired, m(r)) : r in roots_paired] where m:=fromG(g): g in gp]:
        gp in [* GRR, G*]
        ]);
    flat := func<x | &cat x>;
    unflat := func<x | [[x[i], x[i+1]] : i in [1..#x by 2]]>;
    permoo := func<oo, p | unflat([oof[i] : i in p ]) where oof := flat(oo)>;
    // the first return value is the conjugate ootypes for characters of order 4
    // the second return value is the conjugate ootypes for characters of order 2 and 1
    return func<oo | [permoo(oo, p) : p in GRR]>, func<oo | [permoo(oo, p) : p in GCC]>;
end function;

function possible_infinity_types(K, ootype)
    // there must be a smarter way to do this
    ootype_to_mset := func<ootype | {* {* x : x in t *} : t in ootype *}>;
    mset := ootype_to_mset(ootype);
    possibilities := ootype;
    possibilities cat:=[[elt[2], elt[1]] : elt in ootype | elt[2] ne elt[1]];
    ootypes := { [x : x in elt] : elt in CartesianPower(possibilities, #ootype) | ootype_to_mset(elt) eq mset };
    conj4, conj2 := galois_conjugate_ootypes(K);
    // pick orbit representatives
    ootypes4 := {};
    ootypes2 := {};
    for oo in ootypes do
        oo4 := SequenceToSet(conj4(oo));
        oo2 := SequenceToSet(conj2(oo));
        if IsDisjoint(ootypes4, oo4) then
            Include(~ootypes4, oo);
        end if;
        if IsDisjoint(ootypes2, oo2) then
            Include(~ootypes2, oo);
        end if;
    end for;
    return ootypes4, ootypes2, ootypes;
end function;
function possible_GC(K, conductor_support, ootype)
    assert IsCMField(K);
    assert Degree(K) eq 2*#ootype;
    assert #Set([&+t : t in ootype]) eq 1; // "All oo-type parts must have same weight";
    ensure_field(K);
    // we want all the characters with the same modulus
    time psis := possible_hecke_characters(K, conductor_support);
    N := Universe(psis)`Modulus;
    "#chis = ", #psis;
    oo4, oo2, ooall := possible_infinity_types(K, ootype);
    #oo4, #oo2, #ooall;
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
    ps := [p : p in PrimesUpTo(1000) | not p in conductor_support and #PrimeIdealsOverPrime(K, p) eq Degree(K)][1..10];
    return filter_integral(gcs, ps);
end function;

function match_GC(K, conductor_support, ootype, euler_factors)
    time GCS := SequenceToSet(possible_GC(K, conductor_support, ootype));
    for pef in euler_factors do
        p, ef := Explode(pef);
        GCS := { g : g in GCS | EulerFactor(g, p : Integral:=true) eq ef};
        p, #GCS;
    end for;
    return GCS;
end function;


// TODO a parallel version that checks for matching in parallel
// i.e., parallelizes on the hecke characters
