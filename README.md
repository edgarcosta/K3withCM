# K3withCM

This depends on several utils available in [CHIMP](https://github.com/edgarcosta/CHIMP), which in examples below I am assuming that is loaded at startup via the environment variable `MAGMA_USER_SPEC`.

# Constructing the quasi-characters and associated L-functions

```
echo "GCs;" | time magma -b  examples/construct_gc.m
Loading "examples/data.m"
[
[
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair (1,$.1) with modulus of norm 64 over Number Field with defining polynomial x^6 + 6*x^4 + 9*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair ($.1,$.1) with modulus of norm 3136 over Number Field with defining polynomial x^6 + 5*x^4 + 6*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair ($.1,$.1) with modulus of norm 23104 over Number Field with defining polynomial x^6 + 13*x^4 + 50*x^2 + 49 over the Rational Field,
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair ($.1^3,$.1) with modulus of norm 61504 over Number Field with defining polynomial x^6 + 21*x^4 + 116*x^2 + 64 over the Rational Field
],
[
Grossencharacter of type [[ 0, 1 ],[ 1, 0 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1,$.1^3) with modulus of norm 4096 over Number Field with defining polynomial x^6 + 6*x^4 + 9*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 1, 0 ],[ 1, 0 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1^3,$.1) with modulus of norm 25088 over Number Field with defining polynomial x^6 + 5*x^4 + 6*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 0, 1 ],[ 1, 0 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1,$.1) with modulus of norm 184832 over Number Field with defining polynomial x^6 + 13*x^4 + 50*x^2 + 49 over the Rational Field,
Grossencharacter of type [[ 1, 0 ],[ 1, 0 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1^3,$.1) with modulus of norm 3936256 over Number Field with defining polynomial x^6 + 21*x^4 + 116*x^2 + 64 over the Rational Field
],
[
Grossencharacter of type [[ 0, 2 ],[ 2, 0 ],[ 1, 1 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 6*x^4 + 9*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 1, 1 ],[ 0, 2 ],[ 0, 2 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 5*x^4 + 6*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 0, 2 ],[ 2, 0 ],[ 1, 1 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 13*x^4 + 50*x^2 + 49 over the Rational Field,
Grossencharacter of type [[ 1, 1 ],[ 0, 2 ],[ 0, 2 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 21*x^4 + 116*x^2 + 64 over the Rational Field
]
]
magma -b examples/construct_gc.m  36.04s user 1.40s system 250% cpu 14.938 total
```

## We can also verify that the L-functions match up to some bound
```
echo "time check_construct(1000);" | time magma -b  examples/check_gc.m
Loading "examples/construct_gc.m"
Loading "examples/data.m"
true
Time: 6.510
magma -b examples/check_gc.m  45.74s user 1.41s system 212% cpu 22.212 total
```


# Search for the quasi-characters and associated L-function


Some of these searches take a long time, so we suggest to run them with many jobs and one search at the time.

Here is an example of searching for psi_X for the first example:

```
time magma -b verbose:=2 target:=1 psi:=X jobs:=1 < examples/search_gc.m
Loading "examples/data.m"
i= 1 , Conductor bound: for psi_X [<"8.1",7>,<"9.1",1>]
Assembling candidates list...Time: 6.670
#psis = 6144
#matches = 6
#matches up to Galois = 1
i= 1 , 1  match(es) for psi_X of conductor [<"8.1",2>]
magma -b target:=1 psi:=X verbose:=1 jobs:=1 < examples/search_gc.m  50.06s user 0.65s system 99% cpu 50.742 total
```



## Search timings

### psi_A

```
for i in {1..4}; do /usr/bin/time magma -b verbose:=1 target:=$i psi:=A jobs:=512 < examples/search_gc.m ; done
Loading "examples/data.m"
i= 1 , Conductor bound: for psi_A [<"8.1",7>,<"9.1",1>]
i= 1 , 1  match(es) for psi_A of conductor [<"8.1",4>]
133.05user 12.83system 0:13.64elapsed 1069%CPU (0avgtext+0avgdata 61072maxresident)k
0inputs+10448outputs (0major+2070336minor)pagefaults 0swaps
Loading "examples/data.m"
i= 2 , Conductor bound: for psi_A [<"8.1",7>,<"49.1",1>]
i= 2 , 1  match(es) for psi_A of conductor [<"8.1",3>,<"49.1",1>]
139.95user 13.06system 0:12.40elapsed 1233%CPU (0avgtext+0avgdata 56692maxresident)k
0inputs+10464outputs (0major+2183451minor)pagefaults 0swaps
Loading "examples/data.m"
i= 3 , Conductor bound: for psi_A [<"8.1",7>,<"361.1",1>]
i= 3 , 1  match(es) for psi_A of conductor [<"8.1",3>,<"361.1",1>]
398.78user 11.95system 0:17.09elapsed 2402%CPU (0avgtext+0avgdata 43380maxresident)k
0inputs+10448outputs (0major+2135937minor)pagefaults 0swaps
Loading "examples/data.m"
i= 4 , Conductor bound: for psi_A [<"2.1",7>,<"2.2",7>,<"2.3",7>,<"961.1",1>]
i= 4 , 1  match(es) for psi_A of conductor [<"2.1",4>,<"2.2",4>,<"2.3",4>,<"961.1",1>]
44415.75user 5718.03system 1:03:39elapsed 1312%CPU (0avgtext+0avgdata 820212maxresident)k
0inputs+980608outputs (0major+412751930minor)pagefaults 0swaps
```



### psi_X
