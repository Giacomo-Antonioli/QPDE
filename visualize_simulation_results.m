function visualize_simulation_results(u_generic, u_quantum, d, N,ground_truth)
    figure('Position', [100, 100, 1500, 500]);
    flag=true;

    if nargin < 5

           flag=false;
    end

   if d == 3
    % --- 3D Visualization ---
    u_class = reshape(real(u_generic), N, N, N);
    u_quant = reshape(real(u_quantum), N, N, N);
    nx = N; ny = N; nz = N;

    if flag
        ground_truth_3d = reshape(real(ground_truth), N, N, N);

        % 1. Ground Truth
        subplot(2, 3, 1)
        slice(ground_truth_3d, [], [], 1:nz)
        hold on
        slice(ground_truth_3d, 1:nx, [], [])
        slice(ground_truth_3d, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Ground Truth')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

        % 2. Classical Solution
        subplot(2, 3, 2)
        slice(u_class, [], [], 1:nz)
        hold on
        slice(u_class, 1:nx, [], [])
        slice(u_class, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Classical Solution')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

        % 3. Error: Ground Truth vs Classical
        err_class = abs(ground_truth_3d - u_class);
        subplot(2, 3, 3)
        slice(err_class, [], [], 1:nz)
        hold on
        slice(err_class, 1:nx, [], [])
        slice(err_class, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Classical Error')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

        % 5. Quantum Solution
        subplot(2, 3, 5)
        slice(u_quant, [], [], 1:nz)
        hold on
        slice(u_quant, 1:nx, [], [])
        slice(u_quant, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Quantum Solution')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

        % 6. Error: Ground Truth vs Quantum
        err_quant = abs(ground_truth_3d - u_quant);
        subplot(2, 3, 6)
        slice(err_quant, [], [], 1:nz)
        hold on
        slice(err_quant, 1:nx, [], [])
        slice(err_quant, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Quantum Error')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

    else
        % 1. Classical Solution
        subplot(1, 3, 1)
        slice(u_class, [], [], 1:nz)
        hold on
        slice(u_class, 1:nx, [], [])
        slice(u_class, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Classical Solution')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

        % 2. Quantum Solution
        subplot(1, 3, 2)
        slice(u_quant, [], [], 1:nz)
        hold on
        slice(u_quant, 1:nx, [], [])
        slice(u_quant, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Quantum Solution')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off

        % 3. Absolute Error
        err = abs(u_class - u_quant);
        subplot(1, 3, 3)
        slice(err, [], [], 1:nz)
        hold on
        slice(err, 1:nx, [], [])
        slice(err, [], 1:ny, [])
        shading interp; axis equal tight; colorbar
        title('Absolute Error')
        xlabel('x'); ylabel('y'); zlabel('z'); view(3); hold off
    end

else
    % --- 2D Visualization ---
    u_quantum = reshape(real(u_quantum), N, N);

    if flag
        % 1. Ground Truth
        subplot(2, 3, 1)
        imagesc(ground_truth');
        axis square; colorbar;
        title('Ground Truth');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

        % 2. Classical Solution
        subplot(2, 3, 2)
        imagesc(u_generic');
        axis square; colorbar;
        title('Classical Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

        % 3. Error: Ground Truth vs Classical
        subplot(2, 3, 3)
        imagesc(abs(ground_truth - u_generic)');
        axis square; colorbar;
        title('Classical Error');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

        % 5. Quantum Solution
        subplot(2, 3, 5)
        imagesc(u_quantum');
        axis square; colorbar;
        title('Quantum Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

        % 6. Error: Ground Truth vs Quantum
        subplot(2, 3, 6)
        imagesc((ground_truth - u_quantum)');
        axis square; colorbar;
        title('Quantum Error');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

    else
        % 1. Classical Solution
        subplot(1, 3, 1)
        imagesc(u_generic');
        axis square; colorbar;
        title('Classical Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

        % 2. Quantum Solution
        subplot(1, 3, 2)
        imagesc(u_quantum');
        axis square; colorbar;
        title('Quantum Solution');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

        % 3. Absolute Error
        subplot(1, 3, 3)
        imagesc(abs(u_generic - u_quantum)');
        axis square; colorbar;
        title('Absolute Error');
        xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');
    end
end
end