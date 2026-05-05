function dadt = rhs_coeff_RK4(a, c)
    N = length(a)-1;
    b = b_Chebyshev_coeff_RK4(a, N);
    dadt = -c * b;
end