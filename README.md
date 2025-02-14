# K3withCM

This depends on several utils available in [CHIMP](https://github.com/edgarcosta/CHIMP), which in examples below I am assuming that is loaded at startup via the environment variable `MAGMA_USER_SPEC`.

# Constructing the quasi-characters and associated L-functions

```
time magma -b examples/construct_gc.m
Loading "examples/data.m"
> GCs;
[
[
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair (1,$.1) with modulus of norm 64 over Number Field with defining polynomial x^6 + 6*x^4 + 9*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair ($.1,$.1) with modulus of norm 3136 over Number Field with defining polynomial x^6 + 5*x^4 + 6*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair ($.1,$.1) with modulus of norm 23104 over Number Field with defining polynomial x^6 + 13*x^4 + 50*x^2 + 49 over the Rational Field,
Grossencharacter of type [[ 2, 0 ],[ 1, 1 ],[ 1, 1 ]] for Hecke-Dirichlet pair ($.1,$.1) with modulus of norm 61504 over Number Field with defining polynomial x^6 + 21*x^4 + 116*x^2 + 64 over the Rational Field
],
[
Grossencharacter of type [[ 0, 1 ],[ 1, 0 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1,$.1^3) with modulus of norm 4096 over Number Field with defining polynomial x^6 + 6*x^4 + 9*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 1, 0 ],[ 1, 0 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1^3,$.1) with modulus of norm 25088 over Number Field with defining polynomial x^6 + 5*x^4 + 6*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 0, 1 ],[ 0, 1 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1^3,$.1^3) with modulus of norm 184832 over Number Field with defining polynomial x^6 + 13*x^4 + 50*x^2 + 49 over the Rational Field,
Grossencharacter of type [[ 1, 0 ],[ 0, 1 ],[ 0, 1 ]] for Hecke-Dirichlet pair ($.1,$.1^3) with modulus of norm 3936256 over Number Field with defining polynomial x^6 + 21*x^4 + 116*x^2 + 64 over the Rational Field
],
[
Grossencharacter of type [[ 0, 2 ],[ 2, 0 ],[ 1, 1 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 6*x^4 + 9*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 1, 1 ],[ 0, 2 ],[ 0, 2 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 5*x^4 + 6*x^2 + 1 over the Rational Field,
Grossencharacter of type [[ 0, 2 ],[ 2, 0 ],[ 1, 1 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 13*x^4 + 50*x^2 + 49 over the Rational Field,
Grossencharacter of type [[ 1, 1 ],[ 0, 2 ],[ 0, 2 ]] for Hecke-Dirichlet pair (1,1) with modulus of norm 1 over Number Field with defining polynomial x^6 + 21*x^4 + 116*x^2 + 64 over the Rational Field
]
]
>
magma -b examples/construct_gc.m  54.84s user 1.31s system 128% cpu 43.648 total
```

## We can also verify that the L-functions match up to some bound
```
time magma -b examples/check_gc.m
Loading "examples/construct_gc.m"
Loading "examples/data.m"
> time check_construct(1000);
true
Time: 5.990
>
magma -b examples/check_gc.m  57.16s user 1.37s system 135% cpu 43.320 total
```


# Search for the quasi-characters and associated L-function


Some of these searches take a long time, so we suggest to run them with many jobs and one search at the time.

Here is an example of searching for psi_X for the first example:

```
time magma -b verbose:= target:=1 psi:=X jobs:=1 < examples/search_gc.m
Loading "examples/data.m"
i= 1 , Conductor bound: for psi_X [<"8.1",7>,<"9.1",1>]
i= 1 trying d= [] #psis= 1
i= 1 trying d= [<"8.1",1>] #psis= 1
i= 1 trying d= [<"8.1",2>] #psis= 4
i= 1 found  1 character(s)
i= 1 trying d= [<"8.1",3>] #psis= 4
i= 1 trying d= [<"8.1",4>] #psis= 20
i= 1 trying d= [<"8.1",5>] #psis= 40
i= 1 trying d= [<"8.1",6>] #psis= 72
i= 1 trying d= [<"8.1",7>] #psis= 136
i= 1 trying d= [<"9.1",1>] #psis= 2
i= 1 trying d= [<"8.1",1>,<"9.1",1>] #psis= 2
i= 1 trying d= [<"8.1",2>,<"9.1",1>] #psis= 12
i= 1 trying d= [<"8.1",3>,<"9.1",1>] #psis= 12
i= 1 trying d= [<"8.1",4>,<"9.1",1>] #psis= 72
i= 1 trying d= [<"8.1",5>,<"9.1",1>] #psis= 144
i= 1 trying d= [<"8.1",6>,<"9.1",1>] #psis= 272
i= 1 trying d= [<"8.1",7>,<"9.1",1>] #psis= 528
i= 1 , 1  match(es) for psi_X of conductor [<"8.1",2>]
magma -b verbose:= target:=1 psi:=X jobs:=1 < examples/search_gc.m  80.70s user 0.98s system 99% cpu 1:21.69 total
```
