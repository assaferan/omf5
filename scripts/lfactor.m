function getBadPrimesList(p, disc)
    ZZx<x> := PolynomialRing(Integers());
    possible_euler := [1, 1-x, 1+x];
    a := p^2;
    possible_euler cat:= [1+a*x, 1-a*x];
    max_abs_a := Floor(2*SquareRoot(p)^3);
    euler2 := [1 + a*x + p^3*b*x^2 : a in [-max_abs_a..max_abs_a], b in [-1..1]];
    possible_euler cat:= euler2;
    BadPrimesList := [[]];
    for e in possible_euler do
	Append(~BadPrimesList, [<p, e>]);
    end for;
    return BadPrimesList;
end function;

function test_coeff(lser, p, deg, a_range)
    euler_factor<x> := EulerFactor(lser, p);
    euler_factor -:= Coefficient(euler_factor, deg) * x^deg;
    cfes := [CFENew(ChangeEulerFactor(lser, p, a*x^deg + euler_factor)) : a in a_range];
    local_min := [i : i in [2..#a_range-1] | cfes[i-1] ge cfes[i] and cfes[i+1] ge cfes[i] ];
    cands := [a_range[i+delta]*x^deg + euler_factor : i in local_min, delta in [-2..2]];
    return cands;
end function;
