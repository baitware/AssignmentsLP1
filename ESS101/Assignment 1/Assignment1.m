%% Assignment 1 - Part 1a)
clear; clc;

syms m1 m2 L g real
syms p1_x p1_y p1_z theta phi v real
syms p1_x_dot p1_y_dot p1_z_dot theta_dot phi_dot v_dot real
syms p1_x_dotdot p1_y_dotdot p1_z_dotdot theta_dotdot phi_dotdot real
syms ux uy uz real

u = [ux; uy; uz];

p1 = [p1_x; p1_y; p1_z];
p1_dot = [p1_x_dot; p1_y_dot; p1_z_dot];
p1_dotdot = [p1_x_dotdot; p1_y_dotdot; p1_z_dotdot]; 
p2 = p1 + L * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)];

q = [p1; theta; phi];
q_dot = [p1_dot; theta_dot; phi_dot];
q_dotdot = [p1_dotdot; theta_dotdot; phi_dotdot];
    
p2_dot = jacobian(p2, q) * q_dot;

T = 0.5*m1*(p1_dot.'*p1_dot) + 0.5*m2*(p2_dot.'*p2_dot);
V1 = g*m1*[0 0 1]*p1;
V2 = g*m2*[0 0 1]*p2;
V = V1 + V2;

Lag = T-V;

Q = [u; 0; 0];

dLdq = jacobian(Lag, q).';
dLdq_dot = jacobian(Lag, q_dot).';
ddt_dLdq_dot = jacobian(dLdq_dot, [q; q_dot]) * [q_dot; q_dotdot];

E = simplify(ddt_dLdq_dot - dLdq - Q); % = 0
    
M_1a = simplify(jacobian(E, q_dotdot));
b_1a = simplify(-(E-M_1a*q_dotdot));

LHS_1a = M_1a * q_dotdot;

%% Part 1b)

syms Z zero real
syms p2_x p2_y p2_z real
syms p2_x_dot p2_y_dot p2_z_dot real
syms p2_x_dotdot p2_y_dotdot p2_z_dotdot real

p2 = [p2_x; p2_y; p2_z];
p2_dot = [p2_x_dot; p2_y_dot; p2_z_dot];
p2_dotdot = [p2_x_dotdot; p2_y_dotdot; p2_z_dotdot];

q = [p1; p2];
q_dot = [p1_dot; p2_dot];
q_dotdot = [p1_dotdot; p2_dotdot];

e = p1 - p2;

C = 0.5 * (e.'*e - L^2);

V = g * (m1*p1_z + m2*p2_z);
T = 0.5*m1*(p1_dot.'*p1_dot) + 0.5*m2*(p2_dot.'*p2_dot);

Lag = T - V - Z.' * C;

gradLdq = jacobian(Lag, q).';
gradLdq_dot = jacobian(Lag, q_dot).';
ddt_gradLdq_dot = jacobian(gradLdq_dot, [q; q_dot]) * [q_dot; q_dotdot];

Q = [u; 0; 0; 0];

dCdq = jacobian(C, q);
ddq_dCdq = jacobian(dCdq * q_dot, q);

E = simplify(ddt_gradLdq_dot - gradLdq - Q);

M_1b = simplify(jacobian(E, q_dotdot));
b_1b = simplify(-(E-M_1b*q_dotdot));

C_dotdot = dCdq * q_dotdot + ddq_dCdq * q_dot;

%% Part 2a)

a = jacobian(C, q).';

cTop = -simplify(E - M_1b * q_dotdot - jacobian(E, Z)*Z);

cBot = -simplify(C_dotdot - jacobian(C_dotdot, q_dotdot)*q_dotdot);

LHS_mat = [M_1b a; a.' 0];

LHS = LHS_mat*[q_dotdot; Z];

c = [cTop; cBot];

%% Part 2b)


x = LHS_mat \ c;

LHS_mat_inv = inv(LHS_mat);

RHS = simplify(LHS_mat_inv * c);

%% === 2b Simulation & Plot =================================================
% Export numeric functions
Mfun    = matlabFunction(M_1b,    'Vars',{[p1_x; p1_y; p1_z; p2_x; p2_y; p2_z], m1, m2});
afun    = matlabFunction(a,       'Vars',{[p1_x; p1_y; p1_z; p2_x; p2_y; p2_z]});
cTopfun = matlabFunction(cTop,    'Vars',{[p1_x; p1_y; p1_z; p2_x; p2_y; p2_z], ...
                                         [p1_x_dot; p1_y_dot; p1_z_dot; p2_x_dot; p2_y_dot; p2_z_dot], ...
                                         [ux; uy; uz], m1, m2, L, g});
cBotfun = matlabFunction(cBot,    'Vars',{[p1_x; p1_y; p1_z; p2_x; p2_y; p2_z], ...
                                         [p1_x_dot; p1_y_dot; p1_z_dot; p2_x_dot; p2_y_dot; p2_z_dot]});

% Parameters
m1v = 1.5; m2v = 1.0; Lv = 1.2; gv = 9.81;

% Control: gentle lateral motion, thrust balances weight
ufun = @(t,q,v) [ 3*sin(0.6*t);        % ux(t)
                  2*cos(0.4*t);        % uy(t)
                  (m1v+m2v)*gv + cos(0.8*t) ];      % uz ~ hover

% Consistent initial conditions (C(q0)=0, a(q0)^T v0 = 0):
p1_0 = [0; 0; 1.5];
p2_0 = p1_0 - [0; 0; Lv];   % vertical cable of length L
q0   = [p1_0; p2_0];
v0   = zeros(6,1);
y0   = [q0; v0];

% ODE: y=[q; v], qdot=v, vdot=qdd from KKT solve
rhs = @(t,y) kkt_rhs(t,y,Mfun,afun,cTopfun,cBotfun,ufun,m1v,m2v,Lv,gv);
tspan = [0 15];
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,Y] = ode45(rhs, tspan, y0, opts);

% Extract trajectories
Q  = Y(:,1:6);
P1 = Q(:,1:3);
P2 = Q(:,4:6);

%% Pretty 3D plot for the report
figure('Color','w'); hold on; grid on; box on;
plot3(P1(:,1),P1(:,2),P1(:,3),'LineWidth',2);                 % helicopter
plot3(P2(:,1),P2(:,2),P2(:,3),'--','LineWidth',2);            % payload
% cable snapshots
idx = round(linspace(1,size(P1,1),12));
for k = idx
    plot3([P1(k,1) P2(k,1)], [P1(k,2) P2(k,2)], [P1(k,3) P2(k,3)], 'LineWidth',1);
end
% start/end markers
plot3(P1(1,1),P1(1,2),P1(1,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
plot3(P2(1,1),P2(1,2),P2(1,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
plot3(P1(end,1),P1(end,2),P1(end,3),'s','MarkerFaceColor','k','MarkerEdgeColor','k');
plot3(P2(end,1),P2(end,2),P2(end,3),'s','MarkerFaceColor','k','MarkerEdgeColor','k');

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Helicopter (m_1) and Payload (m_2) Trajectories');
legend({'helicopter m_1','payload m_2','cable snapshots'}, 'Location','best');
axis equal; view(35,20);

% Optional: top view and save
figure('Color','w'); plot(P1(:,1),P1(:,2),'LineWidth',2); hold on;
plot(P2(:,1),P2(:,2),'--','LineWidth',2); axis equal; grid on;
xlabel('x'); ylabel('y'); title('Top view (x–y)'); legend('m_1','m_2');
% saveas(gcf,'kkt_helicopter_payload.png');

%% ---------- local function ----------
function ydot = kkt_rhs(t,y,Mfun,afun,cTopfun,cBotfun,ufun,m1v,m2v,Lv,gv)
    q = y(1:6); v = y(7:12);
    u = ufun(t,q,v);
    M = Mfun(q, m1v, m2v);       % 6x6 (blkdiag(m1*I3, m2*I3))
    A = afun(q);                 % 6x1  ([-d; d] with d=p1-p2)
    cTop = cTopfun(q, v, u, m1v, m2v, Lv, gv);   % 6x1
    cBot = cBotfun(q, v);                          % 1x1
    K = [M, A; A.', 0];
    rhs = [cTop; cBot];
    sol = K \ rhs;                % [qdd; z]
    qdd = sol(1:6);
    ydot = [v; qdd];
end
