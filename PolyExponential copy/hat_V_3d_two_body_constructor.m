function [V] = hat_V_3d_two_body_constructor(Nx, g)
%HAT_V_3D_TWO_BODY_CONSTRUCTOR [hat_V]ij = -g
%   INPUT: Nx: an integer specifying the amount of points to split up each
%               direction (nx, ny, nz, mx, my, mz)
%          g: a scalar representing an interaction constant

%   OUTPUT: V: an Nx^3*Nx^3 matrix

V = - ones(Nx^3, Nx^3) * g;

end

