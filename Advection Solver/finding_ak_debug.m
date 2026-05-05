xxxr4eclear 
close all
clc

%%
    
% f = x.^3;
% 
% c = ones(1,N+1);
% c(1) = 0.5; c(end) = 0.5;
% 
% a = zeros(1,N+1);
% for k = 0:N
%     a(k+1) = (2/N) * sum( c .* f .* cos(k*theta) );
% end
% 
% a(1)   = a(1)/2;
% a(end) = a(end)/2;
% disp(a)

N = 15;
j = 0:N;
theta = j*pi/N;
x = cos(theta);
f = x.^2;
a = zeros(1,N+1);

for n = 0:N
    summ = 0;
    for i = 0:N
        c_j = 1;
        if(i==0 || i==N)
            c_j = 2;
        end
        summ = summ + ((1/c_j) .* f(i+1) .* cos(n*theta(i+1)));
    end 
    Cn = 1;
    if (n==0 || n==N)
        Cn = 2;
    end
    a(n+1) = (2/(Cn*N)) * summ;
end
disp(a)

