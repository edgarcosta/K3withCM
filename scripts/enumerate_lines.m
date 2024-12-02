AttachSpec("~/projects/CHIMP/CHIMP.spec");
function normalize_lfl(lfl)
  col := [[MonomialCoefficient(f,Parent(f).i) : i in [1..3]] : f in lfl];

  VV := VectorSpace(Parent(col[1][1]),3);

  M1 := Matrix(3,&cat [col[i] : i in [1..3] ]);
  bool, w := IsConsistent(M1,VV!col[4]);
  M2 := DiagonalMatrix([w[i] : i in [1..3]])*M1;

  co_f4 := [[Parent(col[1][1])!1,0,0],[0,1,0],[0,0,1],[1,1,1]];
  M1a := Matrix(3,&cat [co_f4[i] : i in [1..3] ]);
  bool, w := IsConsistent(M1a,VV!co_f4[4]);
  M2a := DiagonalMatrix([w[i] : i in [1..3]])*M1a;

  trm :=  (M2^-1) * M2a;

  lfl2 := [a^trm : a in lfl];
  lfl2 := [a / LeadingCoefficient(a) : a in lfl2];
  return lfl2, trm;
end function;
function GenerateLines(p)
    F := GF(p);
    P2<x,y,z> := ProjectiveSpace(F, 2);
    // where the lines x, y, z, and x + y + z intersect
    intersection_points := [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, -1], [1, 0, -1], [1, -1, 0]];
    M4 := Transpose(Matrix(F, intersection_points));
    valid_line4 := func<l | &and [elt ne 0 : elt in Eltseq(l*M4)]>;
    // now we generate the lines that do not pass through these points
    P2points := [VectorSpace(F, 3) | Eltseq(elt) : elt in Points(P2)];
    lines := Sort([Eltseq(p) : p in P2points | valid_line4(p)]);
    G := SymmetricGroup(3);
    Gelements := [x : x in G];
    res := {};
    for i->l5 in lines do
        a, b, c := Explode(l5);
        new_intersection_points := [[b-c, c-a, a-b], [0, c, -b], [c, 0, -a], [b, -a, 0]];
        M5 := Transpose(Matrix(F, new_intersection_points));
        valid_line_l6 := func<l | &and [elt ne 0 : elt in Eltseq(l*M5)]>;
        l6_lines := [Eltseq(p) : j->p in lines | i lt j and valid_line_l6(Vector(p))];
        for g in G do
            gl5 := PermuteSequence(l5, g);
            // already_there = lambda elt: ((perm_l5, perm(elt)) if perm_l5 < perm(elt) else (perm(elt), perm_l5)) in res
            l6_lines := [ l : l in l6_lines | not {gl5, PermuteSequence(l5, g)} in res];
        end for;
        res join:= {{l5, l6} : l6 in l6_lines};
    end for;
    // Now we consider the choice of 4 of the 6 lines to normalize
    res_extended := [
        Eltseq(Vector([x,y,z])*Transpose(Matrix(CoordinateRing(P2), [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]] cat SetToSequence(r)))) : r in res
    ];
    triple := func<lin | [MonomialCoefficient(lin, m) : m in [x, y, z]]>;
    min_rep := func<r | [triple(elt) : elt in Sort([Sort(normalize_lfl(PermuteSequence(r, g))[5..6]) : g in Sym(6)])[1] ] >;
    minimal_reps := function()
        return {elt[2,1] : elt in ParallelCall(128, min_rep, [<x> : x in res_extended], 1)};
    end function();
    //minimal_reps := {min_rep(r) : r in res_extended};
    minimal_reps := Sort(SetToSequence(minimal_reps));
    return minimal_reps;
end function;
