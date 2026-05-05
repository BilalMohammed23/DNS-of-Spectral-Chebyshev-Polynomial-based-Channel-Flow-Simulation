function a = cheb_coeff_RK4(u, theta)
    N = length(u)-1;
    a = zeros(1,N+1);

    for n = 0:N
        summ = 0;
        for j = 0:N
            c_j = 1;
            if (j==0 || j==N) 
                c_j = 2; 
            end
            summ = summ + ( (1/c_j) * u(j+1) * cos(n*theta(j+1)) );
        end
        c_n = 1;
        if (n==0 || n==N) 
            c_n = 2; 
        end
        a(n+1) = (2/(c_n*N)) * summ;
    end
end