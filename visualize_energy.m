function visualize_energy(energy_true, energy_num, dt)
% PLOT_LOG_ENERGY Plots the log-energy evolution over time.
% Replicates the provided Python/Matplotlib logic.

    % 1. Create Time Vector
    % time = dt * np.arange(0, len(energy_spec))
    time = dt * (0:length(energy_true)-1);
    
    % 2. Compute Logs
    % log_energy_spec = np.log(energy_spec)
    % log_energy_FG = np.log(energy_FG)
    log_energy_true = log(energy_true);
    log_energy_num  = log(energy_num);
    
    % 3. Plotting
    figure;
    
    % plt.plot(time, log_energy_FG, label="$Energy_{num}$")
    plot(time, log_energy_num, 'LineWidth', 1.5, ...
         'DisplayName', '$Energy_{quantum}$');
    hold on;
    
    % plt.plot(time, log_energy_spec, label="$Energy_{true}$")
    plot(time, log_energy_true, '--', 'LineWidth', 1.5, ...
         'DisplayName', '$Energy_{classical}$');
         
    hold off;
    
    % 4. Labels and Legend
    % plt.ylabel("log(Energy)")
    ylabel('log(Energy)');
    
    % plt.xlabel("Time")
    xlabel('Time');
    
    % plt.legend()
    legend('Interpreter', 'latex', 'Location', 'best');
    
    % 5. Save Figure
    % plt.savefig("Heat2D_energy.png", bbox_inches="tight")
    % saveas(gcf, 'Heat2D_energy.png');
end