SetColumns(0); SetAutoColumns(false);
line := GetScriptArguments()[1]; // to be run with magma -I
p, l := Explode(Split(line, ":"));
p := StringToInteger(p);
l := eval l;
F := GF(p);
P2<x,y,z> := ProjectiveSpace(F, 2);
l := [[1,0,0], [0,1,0], [0,0,1], [1,1,1]] cat l;
w1<t>, w2 := WeilPolynomialOfDegree2K3Surface(&*Eltseq(Vector([x,y,z])*Transpose(Matrix(CoordinateRing(P2), l))));
if Degree(w1) eq 7 then
  StripWhiteSpace(line cat ":" cat Sprint(Eltseq(w1 div (t-p))));
end if;
