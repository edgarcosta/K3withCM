AttachSpec("spec");

R3<x,y,z> := PolynomialRing(Integers(), 3, "lex");
Xs := [
  x* y *z *(x^3 - 3 *x^2* z - 3 *x *y^2 - 3 *x *y *z + y^3 + 9* y^2 *z + 6 *y *z^2 + z^3),
  x*y*z*(7* x^3 - 7 *x^2* y + 49 *x^2* z - 21* x* y* z + 98 *x* z^2 + y^3 - 7 *y^2* z + 49 *z^3),
  x *y* z* (49*x^3 - 304*x^2*y + 570*x^2*z + 361*x*y^2 - 2793*x*y*z + 2033*x*z^2 + 361*y^3 + 2888*y^2*z - 5415*y*z^2 + 2299*z^3)
];



Ks := [NumberField(elt) : elt in [ PolynomialRing(Integers()) |
[ 1, 0, 9, 0, 6, 0, 1 ], // http://www.lmfdb.org/NumberField/6.0.419904.1
[ 1, 0, 6, 0, 5, 0, 1 ], // http://www.lmfdb.org/NumberField/6.0.153664.1
[ 49, 0, 50, 0, 13, 0, 1 ], //https://www.lmfdb.org/NumberField/6.0.8340544.1
[ 64, 0, 116, 0, 21, 0, 1 ] // http://www.lmfdb.org/NumberField/6.0.59105344.1
]];

Klabels := ["6.0.419904.1", "6.0.153664.1", "6.0.8340544.1", "6.0.59105344.1"];
Fs := [Subfields(K,3)[1,1] : K in Ks];
LFs := [TateTwist(LSeries(F), -1) : F in Fs];

Cs := [HyperellipticCurve(elt): elt in
[PolynomialRing(Integers()) | [0,1,0,9,0,6,0,1],
[0,7,0,14,0,7,0,1],
[0,278179,0,44441,0,1786,0,1],
[0,1832265664,0,-3694084,0,961,0,1]
]
];


badps := [BadPrimesOfDegree2K3Surface(X) : X in Xs];
//badpsC := [SequenceToSet(PrimeDivisors(Integers()!Conductor(C))) : C in Cs];
badpsC := [{2, 3}, {2, 7}, {2, 19}, {2, 31}]; // just to avoid the annoying warning
ps := [17,13,37,29];


efX := [ PowerSequence(car<Integers(), PolynomialRing(Integers())>) |
  [<17, [1, -6, 255, 3468, 73695, -501126, 24137569]>],
  [<13, [1, -2, 247, 676, 41743, -57122, 4826809]>],
  [<37, [1, 14, -185, -38332, -253265, 26238254, 2565726409]>],
  [<29, [1, -38, -261, 43732, -219501, -26876678, 594823321]>]
];


efC := [ PowerSequence(car<Integers(), PolynomialRing(Integers())>) |
[<17, [1, -6, 15, -52, 255, -1734, 4913]>],
[<13, [1, 4, 7, 40, 91, 676, 2197]>],
[<37, [1, 4, 15, -152, 555, 5476, 50653]>],
[<29, [1, 4, 51, 216, 1479, 3364, 24389]>]
];

efPrime := [ PowerSequence(car<IntegerRing(), PolynomialRing(Integers())>) |
[<17, [1, 42, 1023, 19244, 295647, 3507882, 24137569]> ],
[ <13, [1, 34, 631, 8476, 106639, 971074, 4826809]> ],
[ <37, [1, 82, 4423, 201724, 6055087, 153681202, 2565726409]> ],
[ <29, [1, 74, 3067, 94772, 2579347, 52338794, 594823321]> ]
];


conductorsX := [* LMFDBIdeal(Ks[i], lbl) : i->lbl in ["64.1", "3136.1", "23104.1", "61504.13" ] *];
conductorsC := [* LMFDBIdeal(Ks[i], lbl) : i->lbl in ["4096.1", "25088.1", "184832.1", "3936256.41"] *];

euler_factor_X := func<X, p | Reverse(WeilPolynomialOfDegree2K3Surface(X, p) div (PolynomialRing(Integers()).1 - p))>;
fac := func<N | StripWhiteSpace(Sprint([<LMFDBLabel(elt[1]), elt[2]> : elt in Factorization(N)]))>;

procedure check_dataCX()
    assert &and [&and [<p, euler_factor_X(X, p) > eq elt where p, efp := Explode(elt) : elt in efX[i]] : i->X in Xs];
    assert &and [&and [<p, EulerFactor(C, p)> eq elt where p, efp := Explode(elt) : elt in efC[i]] : i->C in Cs];
    assert efC eq [[<p, EulerFactor(Cs[i], p)>] : i->p in ps];
    assert &and [
            &and [
              <p, -AlternatingSquareCharacteristicPolynomial(efC[i,j,2]) div (EulerFactor(LFs[i], p) * efX[i,j,2])> eq elt
              where p, efp := Explode(elt)
            : j->elt in efPrime[i]]
            : i->C in Cs];
end procedure;

//check_data();


