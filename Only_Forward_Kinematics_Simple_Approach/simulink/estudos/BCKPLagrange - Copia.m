% function [...
%     M1, M2, M3, M4, M5, M6, M7,...
%     T1, T2, T3, T4, T5, T6, T7,...
%     G1, G2, G3, G4, G5, G6, G7,...
%     ...
%     ] = Lagrange()

g = 9.80665; % m/s^2


%LR - nao usa
I1=[48.1312 3.6966 2.6109; 3.6966 52.4583 1.5051; 2.6109 1.5051 84.0654];

%

I2=[7.8993 0.0146 -0.1057; 0.0146 7.1818 -0.4610; -0.1057 -0.4610 2.3295];  % Kg*m^2 
I3=[6.9751 0.0015 -0.0012; 0.0015 6.9735 0.4844; -0.0012 0.4844 0.7234];
I4=[2.0212 -0.1779 -0.2049; -0.1779 3.7773 0.2234; -0.2049 0.2234 2.7220];
I5=[0.0581 -0.0073 0.0056; -0.0073 0.1652 -0.0003; 0.0056 -0.0003 0.1908];
I6=[0.0116 -0.0018 0; -0.0018 0.0155 0; 0 0 0.0163];
I7=[0.0009 0 0; 0 0.0008 0; 0 0 0.0010];

m1=309.0488;
m2=97.9263;
m3=62.47998;
m4=62.1287;
m5=9.673;
m6=3.6457;
m7=0.4753;

m=[m1 m2 m3 m4 m5 m6 m7];

syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7

q=[q1;q2;q3;q4;q5;q6;q7];
dq=[dq1;dq2;dq3;dq4;dq5;dq6;dq7];
ddq=[ddq1;ddq2;ddq3;ddq4;ddq5;ddq6;ddq7];

 l1=0.644;
 l2=0.250;
 l3=0.550;
 l4=0.850;
 l5=0.145;
 excentricidade=l5;
 l6=0.820;
 l7=0.272;



% Transformações de cinemática direta
T_01=translacao('Z',l1)*translacao('Y',q1);
T_12=rotacao('Z',q2)*translacao('Z',l2)*translacao('X',l3);
T_23=rotacao('Y',q3)*translacao('X',l4);
T_34=rotacao('Y',q4)*translacao('Z',excentricidade)*translacao('X',l6);
T_45=rotacao('X',q5);
T_56=rotacao('Y',q6);
T_67=rotacao('X',q7);
T_78=translacao('X',l7);
T_02=T_01*T_12;
T_03=T_02*T_23;
T_04=T_03*T_34;
T_05=T_04*T_45;
T_06=T_05*T_56;
T_07=T_06*T_67;



% Centro de massa de cada elo em relação ao referencial de cada elo. Esses
% dados foram retirados do programa SolidWorks
CG1=[-0.2509; -0.1666; -0.2277; 1];
CG2=[0.0251;-0.0172;-0.3222;1];
CG3=[-0.0004;-0.2293;-0.5496;1];
CG4=[-0.4497;0.0585;0.0154;1];
CG5=[-0.2220; -0.0100; -0.0093;1];
CG6=[-0.1002 ;0.0173; 0;1];
CG7=[-0.1080;0;0;1];



% Jacobiano para o Braço
P0_G1 = T_02*CG2;
JL0_2 = jacobian(P0_G1(1:3),q);

P0_G2 = T_03*CG3;
JL0_3 = jacobian(P0_G2(1:3),q);

P0_G3 = T_04*CG4;
JL0_4 = jacobian(P0_G3(1:3),q);

% Jacobiano para o Pulso
P0_G4 = T_05*CG5;
JL0_5 = jacobian(P0_G4(1:3),q);

P0_G5 = T_06*CG6;
JL0_6 = jacobian(P0_G5(1:3),q);

P0_G6 = T_07*CG7;
JL0_7 = jacobian(P0_G6(1:3),q);

% Derivada do centro de massa no referencial inercial
dP0_G1 = JL0_2*dq;
dP0_G2 = JL0_3*dq;
dP0_G3 = JL0_4*dq;
dP0_G4 = JL0_5*dq;
dP0_G5 = JL0_6*dq;
dP0_G6 = JL0_7*dq;

% Determinando o Jacobiano Angular
% JA0_1 = [[0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_2 = [[0 0 0]', T_02(1:3,3), [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_3 = [[0 0 0]', T_02(1:3,3), T_03(1:3,2), [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_4 = [[0 0 0]', T_02(1:3,3), T_03(1:3,2), T_04(1:3,2), [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_5 = [[0 0 0]', T_02(1:3,3), T_03(1:3,2), T_04(1:3,2), T_05(1:3,1), [0 0 0]', [0 0 0]'];
JA0_6 = [[0 0 0]', T_02(1:3,3), T_03(1:3,2), T_04(1:3,2), T_05(1:3,1), T_06(1:3,2), [0 0 0]'];
JA0_7 = [[0 0 0]', T_02(1:3,3), T_03(1:3,2), T_04(1:3,2), T_05(1:3,1), T_06(1:3,2), T_07(1:3,1)];

% Velocidade Angular do Trilho
% W0_1 = JA0_1*dq;

% Velocidade Angular do Braço 
W0_2 = JA0_2*dq;
W0_3 = JA0_3*dq;
W0_4 = JA0_4*dq;

% Velocidade Angular do Pulso
W0_5 = JA0_5*dq;
W0_6 = JA0_6*dq;
W0_7 = JA0_7*dq;

% Aceleração da gravidade
gravidade = [0;0;-g];

% Energias

% EC1 = ((m1*dP0_G1')*dP0_G1 + W0_2'*T0_1(1:3,1:3)*T1*T0_1(1:3,1:3)'*W0_2)/2; % Trilho
EC2 = ((m(2)*dP0_G1')*dP0_G1 + W0_2'*((T_02(1:3,1:3)*I2*T_02(1:3,1:3)')*W0_2))/2; % Braço
EC3 = ((m(3)*dP0_G2')*dP0_G2 + W0_3'*((T_03(1:3,1:3)*I3*T_03(1:3,1:3)')*W0_3))/2;
EC4 = ((m(4)*dP0_G3')*dP0_G3 + W0_4'*((T_04(1:3,1:3)*I4*T_04(1:3,1:3)')*W0_4))/2;
EC5 = ((m(5)*dP0_G4')*dP0_G4 + W0_5'*((T_05(1:3,1:3)*I5*T_05(1:3,1:3)')*W0_5))/2; % Pulso
EC6 = ((m(6)*dP0_G5')*dP0_G5 + W0_6'*((T_06(1:3,1:3)*I6*T_06(1:3,1:3)')*W0_6))/2;
EC7 = ((m(7)*dP0_G6')*dP0_G6 + W0_7'*((T_07(1:3,1:3)*I7*T_07(1:3,1:3)')*W0_7))/2;
% EP1 = (-m1*gravidade'*(P0_G1(1:3,1))); % Trilho
EP2 = (-m(2)*gravidade'*(P0_G1(1:3,1))); % Braço
EP3 = (-m(3)*gravidade'*(P0_G2(1:3,1))); 
EP4 = (-m(4)*gravidade'*(P0_G3(1:3,1)));
EP5 = (-m(5)*gravidade'*(P0_G4(1:3,1))); % Pulso
EP6 = (-m(6)*gravidade'*(P0_G5(1:3,1)));
EP7 = (-m(7)*gravidade'*(P0_G6(1:3,1)));


L = (EC2 + EC3 + EC4 + EC5 + EC6 + EC7) - (EP2 + EP3 + EP4 + EP5 + EP6 + EP7);

dL_dq1 = diff(L,q(1)); % Trilho
dL_dq2 = diff(L,q(2));
dL_dq3 = diff(L,q(3));
dL_dq4 = diff(L,q(4));
dL_dq5 = diff(L,q(5));
dL_dq6 = diff(L,q(6));
dL_dq7 = diff(L,q(7));


dL_ddq1 = diff(L,dq(1)); % Trilho
dL_ddq2 = diff(L,dq(2));
dL_ddq3 = diff(L,dq(3));
dL_ddq4 = diff(L,dq(4));
dL_ddq5 = diff(L,dq(5));
dL_ddq6 = diff(L,dq(6));
dL_ddq7 = diff(L,dq(7));    

syms q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q7(t)
q1 = q1(t);
q2 = q2(t);
q3 = q3(t);
q4 = q4(t);
q5 = q5(t);
q6 = q6(t);
q7 = q7(t);

dq1 = diff(q1,t);
dq2 = diff(q2,t);
dq3 = diff(q3,t);
dq4 = diff(q4,t);
dq5 = diff(q5,t);
dq6 = diff(q6,t);
dq7 = diff(q7,t);

ddq1 = diff(dq1,t);
ddq2 = diff(dq2,t);
ddq3 = diff(dq3,t);
ddq4 = diff(dq4,t);
ddq5 = diff(dq5,t);
ddq6 = diff(dq6,t);
ddq7 = diff(dq7,t);

q_L = [q1; q2; q3; q4; q5; q6; q7];
dq_L = [dq1; dq2; dq3; dq4; dq5; dq6; dq7];
ddq_L = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6; ddq7];

dL_ddq1 = subs(dL_ddq1);
dL_ddq2 = subs(dL_ddq2);
dL_ddq3 = subs(dL_ddq3);
dL_ddq4 = subs(dL_ddq4);
dL_ddq5 = subs(dL_ddq5);
dL_ddq6 = subs(dL_ddq6);
dL_ddq7 = subs(dL_ddq7);

dL_ddq1_dt = diff(dL_ddq1,t);
dL_ddq2_dt = diff(dL_ddq2,t);
dL_ddq3_dt = diff(dL_ddq3,t);
dL_ddq4_dt = diff(dL_ddq4,t);
dL_ddq5_dt = diff(dL_ddq5,t);
dL_ddq6_dt = diff(dL_ddq6,t);
dL_ddq7_dt = diff(dL_ddq7,t);

T1 = dL_ddq1_dt ;
T2 = dL_ddq2_dt - dL_dq2;
T3 = dL_ddq3_dt - dL_dq3;
T4 = dL_ddq4_dt - dL_dq4;
T5 = dL_ddq5_dt - dL_dq5;
T6 = dL_ddq6_dt - dL_dq6;
T7 = dL_ddq7_dt - dL_dq7;

T1 = subs(T1 ,[q_L dq_L ddq_L],[q dq ddq]);
T2 = subs(T2 ,[q_L dq_L ddq_L],[q dq ddq]);
T3 = subs(T3 ,[q_L dq_L ddq_L],[q dq ddq]);
T4 = subs(T4 ,[q_L dq_L ddq_L],[q dq ddq]);
T5 = subs(T5 ,[q_L dq_L ddq_L],[q dq ddq]);
T6 = subs(T6 ,[q_L dq_L ddq_L],[q dq ddq]);
T7 = subs(T7 ,[q_L dq_L ddq_L],[q dq ddq]);

% M1
M11 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(1));
M12 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(2));
M13 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(3));
M14 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(4));
M15 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(5));
M16 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(6));
M17 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(7));
M1 = [M11 M12 M13 M14 M15 M16 M17];
% M2
M21 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(1));
M22 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(2));
M23 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(3));
M24 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(4));
M25 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(5));
M26 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(6));
M27 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(7));
M2 = [M21 M22 M23 M24 M25 M26 M27];
% M3
M31 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(1));
M32 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(2));
M33 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(3));
M34 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(4));
M35 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(5));
M36 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(6));
M37 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(7));
M3 = [M31 M32 M33 M34 M35 M36 M37];
% M4
M41 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(1));
M42 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(2));
M43 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(3));
M44 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(4));
M45 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(5));
M46 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(6));
M47 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(7));
M4 = [M41 M42 M43 M44 M45 M46 M47];
% M5
M51 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(1));
M52 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(2));
M53 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(3));
M54 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(4));
M55 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(5));
M56 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(6));
M57 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(7));
M5 = [M51 M52 M53 M54 M55 M56 M57];
% M6
M61 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(1));
M62 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(2));
M63 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(3));
M64 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(4));
M65 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(5));
M66 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(6));
M67 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(7));
M6 = [M61 M62 M63 M64 M65 M66 M67];
% M7
M71 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(1));
M72 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(2));
M73 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(3));
M74 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(4));
M75 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(5));
M76 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(6));
M77 = equationsToMatrix(subs(T7,[q dq ddq],[q dq ddq]),ddq(7));
M7 = [M71 M72 M73 M74 M75 M76 M77];

G1 = subs(T1 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G2 = subs(T2 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G3 = subs(T3 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G4 = subs(T4 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G5 = subs(T5 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G6 = subs(T6 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G7 = subs(T7 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));

T1 = matlabFunction(T1);
T2 = matlabFunction(T2);
T3 = matlabFunction(T3);
T4 = matlabFunction(T4);
T5 = matlabFunction(T5);
T6 = matlabFunction(T6);
T7 = matlabFunction(T7);
M1 = matlabFunction(M1);
M2 = matlabFunction(M2);
M3 = matlabFunction(M3);
M4 = matlabFunction(M4);
M5 = matlabFunction(M5);
M6 = matlabFunction(M6);
M7 = matlabFunction(M7);
G1 = matlabFunction(G1);
G2 = matlabFunction(G2);
G3 = matlabFunction(G3);
G4 = matlabFunction(G4);
G5 = matlabFunction(G5);
G6 = matlabFunction(G6);
G7 = matlabFunction(G7);
% 
% end