function u_pu = calculate_control_output(R,B,P,x,x_hat,u_hat)
inv_R = diag([1/R(1,1), 1/R(2,2)]);
K = (inv_R * B') * P;

% --- E. Calculate Feedback & Total Input ---
x_tilde = x - x_hat;
u_tilde = -K * x_tilde;

% Total output voltage
u_pu_raw = u_hat + u_tilde;

% --- F. Output Saturation ---
limit = 1.1;
vq_pu = max(min(u_pu_raw(1), limit), -limit);
vd_pu = max(min(u_pu_raw(2), limit), -limit);

u_pu = [vq_pu; vd_pu];

end

