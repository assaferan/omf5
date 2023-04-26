ap_dict := [AssociativeArray() : i in [1,2]];

B := 300;
primes := PrimesUpTo(B);
disc := 117;

good_aps1 := [-2,12,2,-42,36,18,112,-124,-130,32,-416,128,-206,304,430,428,210,-594,-916,68,258,400,-1340,408,544,-548,-160,1356,1800,-364,-1288,660,1900,2250,928,902,-2670,584,1496,-5008,-5828,612,-2876,-2520,-44,-1714,2522,6968,1344,-834,3532,352,7600,404,-1308,42,-6556,5960,-7664,2172];

good_aps2 := [-1,-4,-64,32,-172];

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

// lp3 := 27*x^2 + ?*x + 1;
// lp13 := (1+13*x)*(1+?*x+13^3*x^2);

lser := ChangeEulerFactor(lser, 3, lp3);
lser := ChangeEulerFactor(lser, 13, lp13);

CFENew(lser);
