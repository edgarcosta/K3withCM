SetColumns(0); SetAutoColumns(false);
example, p := Explode(GetScriptArguments()); // to be run with magma -I
try
  known := [Split(elt, ":")[1..2] : elt in Split(Read("examples.txt"),"\n")];
catch e
  known := [];
end try;
if [example, p] in known then
  exit; //nothing to do
end if;
example := StringToInteger(example);
p := StringToInteger(p);
if not IsPrime(p) then
  exit; // nothing to do
end if;

_<x,y,z> := PolynomialRing(Integers(), 3);
fs := [
  x* y *z *(x^3 - 3 *x^2* z - 3 *x *y^2 - 3 *x *y *z + y^3 + 9* y^2 *z + 6 *y *z^2 + z^3),
  x*y*z*(7* x^3 - 7 *x^2* y + 49 *x^2* z - 21* x* y* z + 98 *x* z^2 + y^3 - 7 *y^2* z + 49 *z^3),
  x *y* z* (49*x^3 - 304*x^2*y + 570*x^2*z + 361*x*y^2 - 2793*x*y*z + 2033*x*z^2 + 361*y^3 + 2888*y^2*z - 5415*y*z^2 + 2299*z^3)
];
f := fs[example];
try
  w1, w2 := WeilPolynomialOfDegree2K3Surface(f, p);
catch e
  w1 := Polynomial([0]);
  // nothing to do
end try;
w1<t> := w1;
if Degree(w1) eq 7 then
  StripWhiteSpace(Join(
    [
    Sprint(example),
    Sprint(p),
    Sprint(Eltseq(w1 div (t-p)))
    ], ":"));
end if;
