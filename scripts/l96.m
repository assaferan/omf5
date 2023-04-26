ap_dict := [AssociativeArray() : i in [1,2]];

B := 300;
primes := PrimesUpTo(B);
disc := 96;

good_aps1 := [-4,0,-4,-44,4,36,8,60,184,-300,36,-196,-512,-20,476,36,516,-1032,228,-312,-92,1316,212,140,624,-20,1428,836,1352,-1148,-892,-84,-3236,-3056,2180,788,104,2380,308,100,-4608,308,4428,-1168,1236,-1208,916,-2524,-2044,9360,-5388,3116,-2172,1800,2780,8632,-6460,-700,3836,-916];

good_aps2 := [-20, -80, 88, -116, -28];

good_primes := [p :  p in primes | disc mod p ne 0];
bad_primes := [p : p in primes | p notin good_primes];

for i in [1..#good_primes] do
    ap_dict[1][good_primes[i]] := good_aps1[i];
end for;

for i in [1..#good_aps2] do
    ap_dict[2][good_primes[i]] := good_aps2[i];
end for;

for p in bad_primes do
    ap_dict[1][p] := 0;
    ap_dict[2][p] := 0;
end for;

function local_factor(p,d)
    _<x> := PowerSeriesRing(Integers());
    ap := ap_dict[1][p];
    if (d le 1) then
	return 1-ap*x+O(x^2);
    end if;
    ap2 := ap_dict[2][p];
    return 1 - ap * x + p*(ap2+1+p^2) * x^2 - ap*p^3*x^3 + p^6*x^4;
end function;

lser := LSeries(4, [-1,0,0,1], disc, local_factor : Precision := 22);

_<x> := EulerFactor(lser,5);

lp2 := 8*x^2 + 4*x+1;
lp3 := (-3*x+1)*(27*x^2+8*x+1);

lser := ChangeEulerFactor(lser, 2, lp2);
lser := ChangeEulerFactor(lser, 3, lp3);

CFENew(lser);
// 0.0000000000000000000000
