function [Z, z] = multivariate_gauss(rho, m, n, N, z)

  generateNewZ = isempty(z);
  if generateNewZ
    disp('generating new z');
    z = randn(m,n,N);
  else
    disp('using old z');
  end

    
    sigma = eye(m,m) + (ones(m,m)-eye(m,m))*rho;

    if rho ~= 1
        a = chol(sigma,'lower');
    else
        a = [ones(m,1) zeros(m,m-1)];
    end
    
    Z = ones(m,n,N);
    
    for i=1:N
        Z(:,:,i) = a*z(:,:,i);
    end
end
