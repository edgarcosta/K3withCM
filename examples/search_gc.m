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


verbose:=assigned verbose;

GCs := function ()

GCs_X := [];
if "X" in psi_target then
  oo := [[0,2], [1,1],[1,1]];
  for i->X in Xs do
    if not i in target then continue; end if;
    K := Ks[i];
    c := ConductorBoundRationalHeckeCharacters(K, badps[i]);
    if verbose then "i=", i, ", Conductor bound: for psi_X", fac(c); end if;
    ef := [<p, euler_factor_X(X, p)> : p in PrimesUpTo(200) | not p in badps[i] and #Factorization(p*Integers(K)) eq Degree(K)];
    matches := [];
    for d in Divisors(c) do
      psis := RationalHeckeCharacters(HeckeCharacterGroup(d));
      if verbose then "i=", i, "trying d=", fac(d), "#psis=", #psis; end if;
      match := GrossencharacterSearch(psis, oo, ef : Jobs:=jobs);
      matches cat:= match;
      if #match gt 0 and verbose then "i=", i, "found ", #match, "character(s)"; end if;
    end for;
    if verbose then "i=", i, ",", #matches, " match(es) for psi_X of conductor", Join([fac(Conductor(elt)) : elt in matches], " "); end if;
    GCs_X[i] := matches;
  end for;
end if;

GCs_A := [];
if "A" in psi_target then
  oo := [[0,1], [0,1], [0,1]];
  for i->C in Cs do
    if not i in target then continue; end if;
    K := Ks[i];
    c := ConductorBoundRationalHeckeCharacters(K, badpsC[i]);
    if verbose then "i=", i, ", Conductor bound: for psi_A", fac(c); end if;
    badpmodel := PrimeDivisors(Integers()!Discriminant(C));
    ef := [<p, EulerFactor(C, p)> : p in PrimesUpTo(200) | not p in badpmodel and #Factorization(p*Integers(K)) eq Degree(K)];
    matches := [];
    for d in Divisors(c) do
      psis := RationalHeckeCharacters(HeckeCharacterGroup(d));
      if verbose then "i=", i, "trying d=", fac(d), "#psis=", #psis; end if;
      match := GrossencharacterSearch(psis, oo, ef : Jobs:=jobs);
      matches cat:= match;
      if #match gt 0 and verbose then "i=", i, "found ", #match, "character(s)"; end if;
    end for;
    if verbose then "i=", i, ",", #matches, " match(es) for psi_A of conductor", Join([fac(Conductor(elt)) : elt in matches], " "); end if;
    GCs_A[i] := matches;
  end for;
end if;


return [GCs_X, GCs_A];
end function();
exit;
