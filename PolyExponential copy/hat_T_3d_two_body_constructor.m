function [T] = hat_T_3d_two_body_constructor(Nx,m)
%HAT_T_3D_TWO_BODY_CONSTRUCTOR [hat_T]ij = d_ij (p^2 / 2m)
%   INPUT: Nx: an integer specifying the amount of points to split up each
%               direction (nx, ny, nz, mx, my, mz)
%          m: a scalar representing the mass of the particles (?)

%   OUTPUT: T: an Nx^3*Nx^3 matrix

% First, create a matrix for i and j: [Nx^3, 3].
% The first row represents the value of nx/mx for the i/j-th index, the
% second represents the value of ny/my, the last represents the value of
% nz/mz.
long_nx = repmat((1:Nx), 1, Nx^2);    % create a [1, Nx^3] vector
long_ny = repelem((1:Nx), Nx);        % create a [1, Nx^2] vector
long_ny = repmat(long_ny, 1, Nx);     % repeat it to make it [1, Nx^3]
long_nz = repelem((1:Nx), Nx^2);      % create a [1, Nx^3] vector

n_indices = [long_nx; long_ny; long_nz];    % this can be used to find n values at i or m values at j


% Next, create a vector for p and q: [Nx, 1]
% The index i will return the prx/pry/prz-th element for nx/ny/nz, and same
% for q and m. 
% pr = (2 pi / Nx) (nx - Nx/2)
pr = (2 * pi() / Nx) * ((1:Nx) - (Nx / 2));

% Now begin constructing T. It is diagonal, so start with all zeros
T = zeros(Nx^3);

% Populate the diagonal
i = 1;
while i <= Nx^3
    x = n_indices(1, i);
    y = n_indices(2, i);
    z = n_indices(3, i);
    %disp(append("For i = ", string(i), ":"))
    %disp("x, y, z:")
    %disp([x, y, z])
    %disp(append("Momentum is ", string(pr(x)^2 + pr(y)^2 + pr(z)^2), " / 2m"))
    momentum = (pr(x)^2 + pr(y)^2 + pr(z)^2) / (2 * m);
    %disp(append("Momentum is ", string(momentum)))
    T(i, i) = momentum;

    i = i + 1;
end

end

