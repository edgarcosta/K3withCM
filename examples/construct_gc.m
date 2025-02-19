load "examples/data.m";
if assigned target then
    target := [ StringToInteger(target) ];
else
    target := [1,2,3,4];
end if;
if assigned psi then
  psi_target := [ psi ];
else
  psi_target := [ "A", "X" ];
end if;
if assigned jobs then
  jobs := StringToInteger(jobs);
else
  jobs := 12;
end if;
GCs, Ls := function ()
  GCs_X := [];
  oo := [[0,2], [1,1],[1,1]];
  for i->c in conductorsX do
    if not i in target then continue; end if;
    match := GrossencharacterSearch(c, oo, efX[i] : Jobs:=jobs);
    // match := GrossencharacterSearch(RationalHeckeCharacters(HeckeCharacterGroup(c)), oo, efX[i] : Jobs:=jobs);
    assert #match eq 1;
    GCs_X[i] := match[1];
  end for;
  LXs := [LSeries(g : Integral:=true) : g in GCs_X];

  GCs_A := [];
  oo := [[0,1], [0,1], [0,1]];
  for i->c in conductorsC do
      if not i in target then continue; end if;
      match := GrossencharacterSearch(c, oo, efC[i] : Jobs:=jobs);
      // match := GrossencharacterSearch(RationalHeckeCharacters(HeckeCharacterGroup(c)), oo, efC[i] : Jobs:=jobs);
      assert #match eq 1;
      GCs_A[i] := match[1];
  end for;
  LAs := [LSeries(g : Integral:=true) : g in GCs_A];

  GCs_prime := [];
  oo := [[0,2], [0, 2],[1,1]];
  for i->K in Ks do
    if not i in target then continue; end if;
    c := 1*Integers(K);
    match := GrossencharacterSearch(c, oo, efPrime[i] : Jobs:=jobs);
    assert #match eq 1;
    // match := GrossencharacterSearch(RationalHeckeCharacters(HeckeCharacterGroup(c)), oo, efPrime[i] : Jobs:=jobs);
    assert #match eq 1;
    GCs_prime[i] := match[1];
  end for;
  Lprimes := [LSeries(g : Integral:=true) : g in GCs_prime];
return [GCs_X, GCs_A, GCs_prime], [LXs, LAs, Lprimes];
end function();
