// for better printing
SetColumns(0);
SetAutoColumns(false);
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

pipe:= jobs gt 12;

if assigned verbose then
  verbose := StringToInteger(verbose);
else
  verbose := 0;
end if;
SetVerbose("GCSearch", verbose);

GCs := function ()

GCs_X := [];
if "X" in psi_target then
  oo := [[0,2], [1,1],[1,1]];
  for i->X in Xs do
    if not i in target then continue; end if;
    K := Ks[i];
    c := ConductorBoundRationalHeckeCharacters(K, badps[i]);
    Sprintf("i=%o Conductor bound: for psi_X with label \"%o\" = %o", i, LMFDBLabel(c), fac(c));
    ef := [<p, euler_factor_X(X, p)> : p in PrimesUpTo(250) | not p in badps[i] and #Factorization(p*Integers(K)) eq Degree(K)];
    matches := GrossencharacterSearch(c, oo, ef : Primitive:=false, Jobs:=jobs, Pipe:=pipe);
    Sprintf("i=%o %o  match(es) for psi_X of conductor %o", i,  #matches, Join([LMFDBLabel(Conductor(elt)) : elt in matches], ", "));
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
    Sprintf("i=%o Conductor bound: for psi_A with label \"%o\" = %o", i, LMFDBLabel(c), fac(c));
    badpmodel := PrimeDivisors(Integers()!Discriminant(C));
    ef := [<p, EulerFactor(C, p)> : p in PrimesUpTo(250) | not p in badpmodel and #Factorization(p*Integers(K)) eq Degree(K)];
    matches := GrossencharacterSearch(c, oo, ef : Primitive:=false, Jobs:=jobs, Pipe:=pipe);
    Sprintf("i=%o %o  match(es) for psi_A of conductor %o", i,  #matches, Join([LMFDBLabel(Conductor(elt)) : elt in matches], ", "));
    GCs_A[i] := matches;
  end for;
end if;


return [GCs_X, GCs_A];
end function();
exit;
