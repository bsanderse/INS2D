O_intrusive = [Diff_intrusive/options.fluid.Re Conv_intrusive];
hat_A = [As; vectorwise_kron(As)];
A_dot_intrusive = O_intrusive*hat_A;
norm(A_dot_intrusive - A_dots)