function [u] = prob_mod_gen_CFD1D(mu, sigma, f, a, b, ga, gb, N)
    h = (b-a)/N;
    mesh = linspace(a,b,N+1)';          % N+1 elementi
    mesh_in = mesh(2:end-1);            % N-1 elementi
    mesh_aux = mesh(1:end-1) + h/2;     % N elementi
    a = mu(mesh_aux);
    b = -1 * a(2:end-1) / (h*h);
    d = (a(1:end-1)+a(2:end)) / (h*h) + sigma(mesh_in);
    e = f(mesh_in);
    e(1) = e(1) + ga*a(1)/(h*h);
    e(end) = e(end) + gb*a(end)/(h*h);
    u = thomas_simmetrico(d,b,e);
end

