function u = solver_Elliptic_generic(f, grids, A, N, dx)
    % 1. Detect dimension directly from the cell array length
    dim = numel(grids);
    
    % 2. Assign inputs
    N_vecs = N; % Already a vector [Nx, Ny, ...]

    % 3. Reshape RHS (f is a flat vector, map back to N-D grid)
    % FIX: Do not call f(grids), just reshape the data
    f_values = reshape(f, N_vecs);
  
    % 4. Forward N-D FFT
    f_h = fftn(f_values);
    
    % 5. Generate Spectral Wave Numbers (2*pi*i*k)
    k_vecs = cell(1, dim);
    for d = 1:dim
        % Standard MATLAB FFT frequencies: [0, 1... N/2-1, -N/2 ... -1]
        k_idx = [0:N_vecs(d)/2-1, -N_vecs(d)/2:-1];
        
        % Scale by domain size (L = N*dx)
        k_d = 2i * pi * k_idx / (N_vecs(d) * dx);
        
        k_d(1) = 1; % Avoid division by zero (handle mean mode later)
        k_vecs{d} = k_d;
    end
  
    % Create N-D spectral grid
    K_grids = cell(1, dim);
    [K_grids{:}] = ndgrid(k_vecs{:});
    % 6. Build Elliptic Operator Tensor: sum( A_ij * ki * kj )
    denom = zeros(N_vecs); 
    for i = 1:dim
        for j = 1:dim
            if A(i,j) ~= 0
                 denom = denom + A(i,j) * K_grids{i} .* K_grids{j};
            end
        end
    end
    
    % 7. Solve in spectral space
    u_h = f_h ./ denom;
    u_h(1) = 0; % Enforce zero mean (remove artifact from k_d(1)=1)
  
    % 8. Inverse N-D FFT
    u = ifftn(u_h);
end