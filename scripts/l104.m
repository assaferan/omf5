ap_dict := [AssociativeArray() : i in [1,2]];

B := 300;
primes := PrimesUpTo(B);
disc := 104;

good_aps1 := [-7,9,-11,34,-33,34,100,-32,-36,-219,402,-209,-7,-158,-638,510,-878,-73,256,484,396,-24,1176,-118,1576,1232,-1915,-1116,-1584,1475,-138,-1453,-372,-469,-884,-864,-624,362,3139,-4274,202,-326,-797,-618,-901,-305,3188,2995,3535,-1987,1452,-3084,-653,1068,-4802,-1991,5322,1832,4072,57];

good_aps2 := [4,-30,-56,-48,-342];

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

lser := LSeries(4, [-1,0,0,1], disc, local_factor : Precision := 23);

_<x> := EulerFactor(lser,5);

lp2 := 8*x^2 + 3*x + 1;
lp13 := (1+13*x)*(1+34*x+13^3*x^2);

lser := ChangeEulerFactor(lser, 2, lp2);
lser := ChangeEulerFactor(lser, 13, lp13);

CFENew(lser);
// 0.00000000000000000000000
