function u = cheb_eval_series_RK4(a, x)
    % a: 1x(N+1) or (N+1)x1, coefficients a0..aN
    % x: vector of points in [-1,1]
    a = a(:).';                 % row
    N = length(a)-1;

    T0 = ones(size(x));
    % if N==0
    %     u = a(1)*T0;
    %     return
    % end
    T1 = x;

    u = a(1)*T0 + a(2)*T1;

    for n = 1:N-1
        T2 = 2*x.*T1 - T0;       % T_{n+1}
        u = u + a(n+2)*T2;
        T0 = T1;
        T1 = T2;
    end
end