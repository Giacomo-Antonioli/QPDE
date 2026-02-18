function GenericDiffusion_QPDE(f_handle, u_init_handle, A, N, d, dt, steps)
% GENERICDIFFUSION_QPDE Generic Quantum Diffusion Solver (2D/3D).
% Uses implicit time-stepping for classical verification.

    %% 1. Setup & Grid Generation
    L = 1; 
    dx = L / N;
    x_vec = (0:N-1) * dx;
 
    % Generate Grids
    ndgrid_inputs = repmat({x_vec'}, 1, d);
    grids = cell(1, d);
    [grids{:}] = ndgrid(ndgrid_inputs{:});
    
    % Evaluate Source (f) and Initial Condition (u)
    f_vals = f_handle(grids{:});
    u_init = u_init_handle(grids{:});
    
    % N vector for solver (e.g. [32, 32])
    N_vecs = repmat(N, 1, d);

    %% 2. Operator Construction
    n = log2(N);
    fprintf('Building Quantum Operator for d=%d, N=%d...\n', d, N);
    
    % A) Quantum Operator (Propagator)
    % Ensure QPDE_Generator.m is in your path
    Q_Op = QPDE_Generator_Diffusion(A, n, dt);
    
    % B) Compute Normalization Alpha
    % We compute alpha by running the solver once on a dummy impulse or 
    % by estimating the max eigenvalue of the inverse operator.
    % For diffusion (I - dt*L)^-1, the max value is 1 (at k=0).
    alpha = 1.0; 

    %% 3. Time Evolution Loop
    u_class = u_init;
    u_quant = u_init;
    % Initialize Energy History
    E_class_history = zeros(1, steps);
    E_quant_history = zeros(1, steps);
    fprintf('Running %d steps...\n', steps);
    
    for t = 1:steps
        % --- Classical Step ---
        % Solving: (I - dt * Div(A * Grad)) u_new = u_old - dt * f
        rhs = u_class - dt * f_vals;
        u_class = solver_Diffusion_generic(rhs, grids, A, N_vecs, dx, dt);
        
        % --- Quantum Step ---
        % 1. Apply Source: v = u - dt * f
        v = u_quant - dt * f_vals; 
        
        % 2. Apply Operator: u_new = Q_Op * v
        u_vec_next = (Q_Op * v(:)) * alpha;
        
        % 3. Reshape
        u_quant = reshape(real(u_vec_next), N_vecs);

        % --- Compute Energy for Both ---
        E_class_history(t) = energy(u_class, A, N, dx,d);
        E_quant_history(t) = energy(u_quant, A, N, dx,d);
    end
    
    %% 4. Visualization & Metrics
  

    visualize_simulation_results(u_class, u_quant, d, N)
    visualize_energy(E_class_history, E_quant_history,dt)
end