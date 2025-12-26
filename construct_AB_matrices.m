function [A_bar, B] = construct_AB_matrices(k, x, we_hat, id_hat)
% CONSTRUCT_AB_MATRICES Constructs the A_bar and B matrices for the SDRE controller.
%   [A_bar, B] = construct_AB_matrices(k, x, we_hat, id_hat)
%   Constructs the state-dependent A_bar matrix and the constant B matrix
%   based on motor coefficients (k), current state (x), desired speed (we_hat),
%   and desired d-axis current (id_hat).

A_bar = [ -k(2),                 (k(1) + k(11)*id_hat),   k(11)*x(2);
          -(k(5) + k(10)*x(3)),  -k(4),                   -k(10)*we_hat;
          k(9)*x(2),             k(9)*we_hat,             -k(7) ];
          
B = [ 0,    0;
      k(6), 0;
      0,    k(8) ];

end
