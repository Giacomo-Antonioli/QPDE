function u = solver_Diffusion_generic(f, grids, A, N, dx, dt)
    dim = numel(grids);
    N_vecs = N;
    f_values = reshape(f, N_vecs);
    f_h = fftn(f_values);
    
    k_vecs = cell(1, dim);
    for d = 1:dim
        k_idx = [0:N_vecs(d)/2-1, -N_vecs(d)/2:-1];
        k_vecs{d} = 2i * pi * k_idx / (N_vecs(d) * dx);
    end
    
    K_grids = cell(1, dim);
    [K_grids{:}] = ndgrid(k_vecs{:});
    
    denom = ones(N_vecs); 
    for i = 1:dim
        for j = 1:dim
            if A(i,j) ~= 0
                 denom = denom - dt * A(i,j) * K_grids{i} .* K_grids{j};
            end
        end
    end
    % denom(1)=1;
    u_h = f_h ./ denom;
    u = real(ifftn(u_h));
end