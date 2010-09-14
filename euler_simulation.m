function S = euler_simulation(s0, T, sigma, Z, r)
    [n N] = size(Z);
    
    dt = T/N;

    S = zeros(size(Z));
    S(:,1) = s0;
    for i= 1:N
        vols = sigma(i*dt, S(:, i));
        S(:,i+1) = S(:,i) .* (1 + r*dt + vols .* sqrt(dt) .* Z(:,i));
    end
end