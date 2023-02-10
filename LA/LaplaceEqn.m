clear;
close all;

%%%%%%%%%%%%%%%%%%%%
%%%   Geometry   %%%
%%%%%%%%%%%%%%%%%%%%

Lx = 1;         Ly = 1;         % rectangle dimensions
Nx = 100;       Ny = 100;       % # of intervals
nx = Nx + 1;    ny = Ny + 1;    % # of gridpoints in x and y, w/ boundaries
dx = Lx/Nx;     dy = Ly/Ny;     % grid spacing
x = (0:Nx)*dx;  y = (0:Ny)*dy;  % x and y positions on grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Interation and Starting Conditions   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cc = 1.e-6;     % convergence criteria
in_x = 2:nx - 1;   in_y = 2:ny - 1;   % internal grid points
V = zeros(nx, ny);  % solution with BC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Boundary Conditions   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V(:,1) = 0;      % Bottom BC
V(:,ny) = 0;     % Top BC

V(1,:) = 1;      % Left BC
V(nx,:) = 0;     % Right BC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Jacobi Iteration   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_old = V;
error = 2*cc;
it_count = 0;

while (error > cc)
    
    % Iteration Counter
    it_count = it_count + 1;
    
    % Interation Solution
%     V(in_x, in_y) = 0.25*(V(in_x-1, in_y) + V(in_x, in_y-1) + V(in_x+1, in_y) +V(in_x, in_y+1));
    V = imboxfilt(V, 3);
    % Reseting BC
    V(:,1) = V(:,2);      % Bottom BC
    V(:,ny) = V(:,ny-1);     % Top BC

    V(1,:) = 1;      % Left BC
    V(nx,:) = 0;     % Right BC

    %Error Calculation
    error = max(abs(V(:) - V_old(:)));

    %Plots in 3D
    if mod(it_count, 250) == 0
        subplot(1,2,1)
        surf(V')
        
        subplot(1,2,2)
        [Ex, Ey] = gradient(V);
        quiver(-Ey',-Ex', 15)

%         subplot(1,3,3)
%         grad = imboxfilt(V, 3);
%         imshowpair(grad);
        pause(0.05)
    end
    
    %Checks if interation diverges by looking for non-numbers or infinities
    if any(isnan(V(:))) || any(isinf(V(:)))
        fprintf('Iterations Diverge\n');
    end

    V_old = V;
end

fprintf('%g\n', it_count);

% [Ex, Ey] = gradient(V);
% figure
% quiver(-Ey',-Ex', 0.8)


%%%%%%%%%%%%%%%%%%%%
%%%   Plotting   %%%
%%%%%%%%%%%%%%%%%%%%

% [xx, yy] = meshgrid(x, y);
% v = [0.8 0.6 0.4 0.2 0.1 0.05 0.01];
% contour(xx, yy, V', v, 'ShowText', 'on');
% axis equal;
% 
% Tick = [0 0.2 0.4 0.6 0.8 1];
% set(gca, 'YTick', Tick);
% set(gca, 'XTick', Tick);
% 
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
% title('Interative Solution of the Laplace Equation', 'Interpreter', 'latex', 'FontSize', 16);