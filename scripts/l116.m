ap_dict := [AssociativeArray() : i in [1,2]];

B := 300;
primes := PrimesUpTo(B);
disc := 116;

good_aps1 := [-9,9,-6,29,-49,26,-12,-48,-37,-90,-308,-87,191,-335,572,-222,652,130,-226,85,-336,-2028,476,1118,-334,22,1569,480,-24,996,1702,140,-2673,-1014,180,-4731,-826,-488,1702,427,-2652,-468,1376,-6330,4645,574,-3222,4132,2717,-136,-3089,1127,5511,-713,1538,-3771,-2204,-1795,740,896];

good_aps2 := [8,-18,-28,-24,-126,-48];

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

// lp2 := 8*x^2 + ?*x + 1;
// lp29 := (1+29*x)*(1+?*x+29^3*x^2);

lser := ChangeEulerFactor(lser, 2, lp2);
lser := ChangeEulerFactor(lser, 29, lp29);

CFENew(lser);
