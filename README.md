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
for i in {1..4}; do time magma -b target:=$i psi:=A jobs:=12 < examples/search_gc.m ; done
Loading "examples/data.m"
i=1 Conductor bound: for psi_A with label "18874368.1" = [<"8.1",7>,<"9.1",1>]
i=1 1  match(es) for psi_A of conductor 4096.1
magma -b target:=$i psi:=A jobs:=12 < examples/search_gc.m  68.36s user 0.50s system 696% cpu 9.886 total
Loading "examples/data.m"
i=2 Conductor bound: for psi_A with label "102760448.1" = [<"8.1",7>,<"49.1",1>]
i=2 1  match(es) for psi_A of conductor 25088.1
magma -b target:=$i psi:=A jobs:=12 < examples/search_gc.m  56.72s user 0.59s system 636% cpu 9.002 total
Loading "examples/data.m"
i=3 Conductor bound: for psi_A with label "757071872.1" = [<"8.1",7>,<"361.1",1>]
i=3 1  match(es) for psi_A of conductor 184832.1
magma -b target:=$i psi:=A jobs:=12 < examples/search_gc.m  76.81s user 0.63s system 724% cpu 10.684 total
Loading "examples/data.m"
i=4 Conductor bound: for psi_A with label "2015363072.113" = [<"2.1",7>,<"2.2",7>,<"2.3",7>,<"961.1",1>]
i=4 1  match(es) for psi_A of conductor 3936256.41
magma -b target:=$i psi:=A jobs:=12 < examples/search_gc.m  9692.18s user 44.60s system 1003% cpu 16:10.06 total
```



### psi_X


```
for i in {1..2}; do time magma -b target:=$i psi:=X jobs:=12 < examples/search_gc.m ; done
Loading "examples/data.m"
i=1 Conductor bound: for psi_X with label "18874368.1" = [<"8.1",7>,<"9.1",1>]
i=1 1  match(es) for psi_X of conductor 64.1
magma -b target:=$i psi:=X jobs:=12 < examples/search_gc.m  61.89s user 0.81s system 328% cpu 19.087 total
Loading "examples/data.m"
i=2 Conductor bound: for psi_X with label "102760448.1" = [<"8.1",7>,<"49.1",1>]
i=2 1  match(es) for psi_X of conductor 3136.1
magma -b target:=$i psi:=X jobs:=12 < examples/search_gc.m  54.66s user 0.78s system 307% cpu 18.035 total


# on a different machine
time magma -b target:=3 psi:=X jobs:=512 < examples/search_gc.m
i=3 Conductor bound: for psi_X with label "157790721460674756608.45" = [<"8.1",7>,<"49.1",1>,<"49.2",1>,<"49.3",1>,<"121.1",1>,<"121.2",1>,<"121.3",1>,<"361.1",1>]
i=3 1  match(es) for psi_X of conductor 23104.1

real    106m41.505s
user    11142m56.536s
sys     28m10.660s
```
