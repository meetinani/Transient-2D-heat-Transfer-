% Parameters
Lx = 0.1; % Length of the device in x-direction (m)
Ly = 0.05; % Length of the device in y-direction (m)
Nx = 30; % Number of grid points in x-direction
Ny = 20; % Number of grid points in y-direction
T_initial = 25; % Initial temperature (°C)
T_source = 100; % Temperature of heat source (°C)
alpha = 1e-4; % Thermal diffusivity (m^2/s)
dt = 0.01; % Time step (s)
t_final = 10; % Final time (s)

% Grid spacing
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);

% Initialize temperature matrix
T = T_initial * ones(Ny, Nx);

% Define heat source location
source_x = (Nx / 2);
source_y = (Ny / 2);

% Main time-stepping loop
num_steps = (t_final / dt);
for k = 1:num_steps
    % Boundary conditions (assuming adiabatic walls)
    T(:, 1) = 150; % Left boundary
    T(:, end) = 150; % Right boundary
    T(1, :) = 75; % Top boundary
    T(end, :) = 75; % Bottom boundary
    
    % Apply heat source
    T(source_y, source_x) = T_source;
    
    % Compute temperature at next time step using finite difference
    T_new = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            T_new(j, i) = T(j, i) + alpha * dt * ((T(j, i+1) - 2*T(j, i) + T(j, i-1)) / dx^2 + ...
                                                    (T(j+1, i) - 2*T(j, i) + T(j-1, i)) / dy^2);
        end
    end
    
    % Update temperature matrix for next time step
    T = T_new;
    
    % Plot temperature distribution
    imagesc(T);
    colorbar;
    title(['Transient Heat Transfer (Time = ', num2str(k*dt), 's)']);
    xlabel('x');
    ylabel('y');
    axis equal;
    drawnow;
end