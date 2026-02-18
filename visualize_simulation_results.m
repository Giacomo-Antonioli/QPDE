function visualize_simulation_results(u_class, u_quant, d, N)
    figure('Position', [100, 100, 1500, 500]);

    if d == 2
        % --- 2D Visualization ---
        u_class = reshape(real(u_class), N, N);
        u_quant = reshape(real(u_quant), N, N);
        
        subplot(1,3,1)
        imagesc(u_class'); 
        axis square; colorbar;
        title('Classical Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
        
        subplot(1,3,2)
        imagesc(u_quant'); 
        axis square; colorbar;
        title('Quantum Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
        
        subplot(1,3,3)
        imagesc(abs(u_class - u_quant)'); 
        axis square; colorbar;
        title('Absolute Error');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
        
    elseif d == 3
        % --- 3D Visualization ---
        u_class = reshape(real(u_class), N, N, N);
        u_quant = reshape(real(u_quant), N, N, N);
        err = abs(u_class - u_quant);
        
        nx = N; ny = N; nz = N;
        
        % 1. Classical Solution
        subplot(1,3,1)
        slice(u_class, [], [], 1:nz)
        hold on
        slice(u_class, 1:nx, [], [])
        slice(u_class, [], 1:ny, [])
        shading interp
        axis equal tight
        colorbar
        title('Classical Solution')
        xlabel('x'); ylabel('y'); zlabel('z')
        view(3)
        hold off
        
        % 2. Quantum Solution
        subplot(1,3,2)
        slice(u_quant, [], [], 1:nz)
        hold on
        slice(u_quant, 1:nx, [], [])
        slice(u_quant, [], 1:ny, [])
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