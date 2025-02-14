load "examples/construct_gc.m";
GCs_X, GCs_A, GCs_prime := Explode(GCs);
LXs, LAs, Lprimes := Explode(Ls);
function check_construct(B)
    return &and [ &and [
    -AlternatingSquareCharacteristicPolynomial(EulerFactor(LAs[i], p)) eq &*[EulerFactor(L, p : Integral:=true)
    : L in [LXs[i], LFs[i], Lprimes[i]]] : p in PrimesUpTo(B) | not p in (badps[i] join badpsC[i])]  : i in [1..4]];
end function;
