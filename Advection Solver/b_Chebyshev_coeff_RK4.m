function b = b_Chebyshev_coeff_RK4(a, N)
    b = zeros(1,N+1);

    b(N+1) = 0;
    b(N)   = 2*N*a(N+1);  %in theory b(N-1) = 2NaN

    for k = N-1:-1:1
        b(k) = b(k+2) + 2*k*a(k+1);
    end

    b(1) = b(1)/2;
end