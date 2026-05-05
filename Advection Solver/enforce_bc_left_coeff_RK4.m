function a = enforce_bc_left_coeff_RK4(a)
    N = length(a)-1;
    s = 0;
    for n = 0:N-1
        s = s + ((-1)^n) * a(n+1);
    end
    a(N+1) = ((-1)^(N+1)) * s; %last a term
end