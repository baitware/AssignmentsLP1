A = [2 1 -2; 1 -1 -3];
b = [0; 0];

x1 = A\b;

A = [4.2 0 1.3; 1 3.5 pi/2; 0 2.2 3];
b = [pi; 1; 2];

x2 = A\b;

A = [4 2 5; 6 0 6; 1 2 3];
b = [6; 3; 2];

x3 = A\b;

%%
% Sym Toolbox practice
% 1:

syms x1 x2 x3 u1 u2 real

x = [x1; x2; x3];
u = [u1; u2];

f = [-x1^2-x2*x3*tanh(u1); -x1*x2*u1*u2; u1*cos(x1)*sin(u2)];


A = jacobian(f,x);
B = jacobian(f,u);

x0 = [1; 1; 1];
u0 = [1; 1];

dx = 1e-3*randn(size(x0));
du = 1e-3*randn(size(u0));

A0 = double( subs(A, [x;u], [x0;u0]));
B0 = double( subs(B, [x;u], [x0;u0]));

delta_x_dot = @(x,u) A0*(x - x0) + B0*(u - u0);

f0 = double( subs(f, [x;u], [x0;u0]));

f_lin = f0 + delta_x_dot(x0 + dx,u0 + du)
f_true = double( subs(f, [x;u], [x0+dx; u0+du]))

%%
% 2:

syms x

f = sin(2*x) + cos(x)*exp(2*x);

f_func = matlabFunction(f, 'Vars', x);

T1_sym = taylor(f, x, 0, 'Order', 1);
T2_sym = taylor(f, x, 0, 'Order', 2);
T3_sym = taylor(f, x, 0, 'Order', 3);
T4_sym = taylor(f, x, 0, 'Order', 4);

T1 = matlabFunction(T1_sym, 'Vars', x)
T2 = matlabFunction(T2_sym, 'Vars', x)
T3 = matlabFunction(T3_sym, 'Vars', x)
T4 = matlabFunction(T2_sym, 'Vars', x)

xn = -2:0.01:2;

plot(xn, f_func(xn), 'k', ...
 xn, T1(xn), '--', ...
 xn, T2(xn), '--', ...
 xn, T3(xn), '--', ...
 xn, T4(xn), '--')

