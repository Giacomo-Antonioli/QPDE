function [grids, N_dir, dx] = GenericElliptic_QPDE(f, A, N, d)
    x_lb = 0;
    x_rb = 1;
    L  = x_rb - x_lb;
    dx = L / N;
    n=log2(N);
    x_vec = x_lb + (0:N-1) * dx;
    
    % Prepare generic inputs for ndgrid
    ndgrid_inputs = repmat({x_vec}, 1, d);
    
    % Initialize 'grids' as a cell array (container for generic data)
    grids = cell(1, d);
    
    % Generate N-dimensional grids
    [grids{:}] = ndgrid(ndgrid_inputs{:});
    
    % N_dir is now a simple vector [N, N, N...]
    N_dir = repmat(N, 1, d);   

    f_vals = f(grids{:});
    f_flat = f_vals(:);
    
   
    
    u_generic = solver_Elliptic_generic(f_flat, grids, A, N_dir, dx);
    op=QPDE_Generator(A,n);
    size(f_flat)
    size(op)
    u_quantum=op*f_flat;
    
    figure;
  
    if d == 2
        % --- 2D Visualization ---
        u_quantum = reshape(real(u_quantum), N, N);
        
        subplot(1,3,1)
        imagesc(u_generic'); 
        axis square; colorbar;
        title('Classical Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
        
        subplot(1,3,2)
        imagesc(u_quantum'); 
        axis square; colorbar;
        title('Quantum Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
        
        subplot(1,3,3)
        imagesc(abs(u_generic - u_quantum)'); 
        axis square; colorbar;
        title('Absolute Error');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
        
    elseif d == 3
      % --- 3D Visualization ---
        u_quantum = reshape(real(u_quantum), N, N, N);
        err = abs(u_generic - u_quantum);
        
        nx = size(u_generic, 1);
        ny = size(u_generic, 2);
        nz = size(u_generic, 3);

        % 1. Classical Solution
        subplot(1,3,1)
        slice(u_generic, [], [], 1:nz)
        hold on
        slice(u_generic, 1:nx, [], [])
        slice(u_generic, [], 1:ny, [])
        shading interp
        axis equal tight
        colorbar
        title('Classical Solution')
        xlabel('x'); ylabel('y'); zlabel('z')
        view(3)
        hold off

        % 2. Quantum Solution
        subplot(1,3,2)
        slice(u_quantum, [], [], 1:nz)
        hold on
        slice(u_quantum, 1:nx, [], [])
        slice(u_quantum, [], 1:ny, [])
        shading interp
        axis equal tight
        colorbar
        title('Quantum Solution')
        xlabel('x'); ylabel('y'); zlabel('z')
        view(3)
        hold off

        % 3. Error
        subplot(1,3,3)
        slice(err, [], [], 1:nz)
        hold on
        slice(err, 1:nx, [], [])
        slice(err, [], 1:ny, [])
        shading interp
        axis equal tight
        colorbar
        title('Absolute Error')
        xlabel('x'); ylabel('y'); zlabel('z')
        view(3)
        hold off
    end
end


