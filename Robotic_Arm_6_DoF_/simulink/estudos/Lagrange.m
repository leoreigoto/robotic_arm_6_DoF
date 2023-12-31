% function [...
%     M1, M2, M3, M4, M5, M6, M7,...
%     T1, T2, T3, T4, T5, T6, T7,...
%     G1, G2, G3, G4, G5, G6, G7,...
%     ...
%     ] = Lagrange()

g = 9.80665; % m/s^2




%LR - nao usa
I1=[6.45927 0.00073 -0.42267; 0.00073 8.55120 0.00155; -0.42267 0.00155 13.15421];

%

I2=[11.08697 3.39167 5.38787; 3.39167 37.87049 0.88684; 5.38787 0.88684 35.38487];  % Kg*m^2 
I3=[6.70839 -7.30351 -0.12339; -7.30351 78.03742 0.03995; -0.12339 0.03995 78.85343];
I4=[2.40106 1.64045 -1.15010; 1.64045 48.35669 -0.11971; -1.15010 -0.11971 47.88720];
I5=[0.02668 0.00251 0.00000; 0.00251 0.11156 0.00000; 0.00000 0.00000 0.10998];
I6=[0.26923 -0.04332 -0.00012; -0.04332 0.22283 -0.00011; -0.00012 -0.00011 0.37601];
I7=[0.05118 -0.0005 -0.00002; -0.00005 0.01547 0.00003; -0.00002 0.00003 0.06098]; %flange+ferramenta
%Itarefa=[0.22649 -0.00006 -0.00005;-0.00006 0.17109 -0.00013;-0.00005 -0.00013 0.23614];
syms Itarefaxx Itarefaxy Itarefaxz Itarefayx Itarefayy Itarefayz Itarefazx Itarefazy Itarefazz;
Itarefa=[Itarefaxx Itarefaxy Itarefaxz ;Itarefayx Itarefayy Itarefayz ;Itarefazx Itarefazy Itarefazz];


% m1=207.61613;
% m2=346.81729;
% m3=366.65763;
% m4=243.34094;
% m5=10.76884;
% m6=27.61521;
% m7=3.40236; %flange+ferramenta
%mtarefa=8.90300;
syms m1 m2 m3 m4 m5 m6 m7 mtarefa;
massa_robo=[m1 m2 m3 m4 m5 m6 m7 mtarefa];


m=[m1 m2 m3 m4 m5 m6 m7 mtarefa];

syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7

q=[q1;q2;q3;q4;q5;q6];
dq=[dq1;dq2;dq3;dq4;dq5;dq6];
ddq=[ddq1;ddq2;ddq3;ddq4;ddq5;ddq6];

%  L1=0.590;
%  L2=1.350;
%  L3=1.600;
%  L4=0.215;
%  xf=0.22485;
%  E1=0.75;
%  E2=-0.041;
%  yf=0;
%  zf=0;
syms L1 L2 L3 L4 xf E1 E2 yf zf
size_robo=[L1 L2 L3 L4 E1 E2 xf yf zf];

% Transforma��es de cinem�tica direta
R1_0 = [cos(q1) -sin(q1) 0 0; sin(q1) cos(q1) 0 0; 0 0 1 0; 0 0 0 1];
D1_0 = [1 0 0 E1; 0 1 0 0; 0 0 1 L1; 0 0 0 1];
T1_0 = R1_0*D1_0;
T_01=T1_0;
R2_1 = [cos(q2) 0 sin(q2) 0; 0 1 0 0; -sin(q2) 0 cos(q2) 0; 0 0 0 1];
D2_1 = [1 0 0 L2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
T2_1 = R2_1*D2_1;
T2_0 = T1_0*T2_1;
T_02=T2_0;
R3_2 = [cos(q3) 0 sin(q3) 0; 0 1 0 0; -sin(q3) 0 cos(q3) 0; 0 0 0 1];
D3_2 = [1 0 0 L3; 0 1 0 0; 0 0 1 E2; 0 0 0 1];
T3_2 = R3_2*D3_2;
T3_0 = simplify(T2_0*T3_2);
T_03=T3_0;
R4_3 = [1 0 0 0; 0 cos(q4) -sin(q4) 0; 0 sin(q4) cos(q4) 0; 0 0 0 1];
T4_3 = R4_3;
T4_0 = simplify(T3_0*T4_3);
T_04=T4_0;
R5_4 = [cos(q5) 0 sin(q5) 0; 0 1 0 0; -sin(q5) 0 cos(q5) 0; 0 0 0 1];
T5_4 = R5_4;
T5_0 = simplify(T4_0*T5_4);
T_05=T5_0;
R6_5 = [1 0 0 0; 0 cos(q6) -sin(q6) 0; 0 sin(q6) cos(q6) 0; 0 0 0 1];
T6_5 = R6_5;
T6_0 = simplify(T5_0*T6_5);
T_06=T6_0;
D7_6 = [1 0 0 L4+xf; 0 1 0 yf; 0 0 1 zf; 0 0 0 1];
T7_6 = D7_6;
T7_0 = simplify(T6_0*T7_6);
T_07=T7_0;


% T_01=translacao('Z',l1)*translacao('Y',q1);
% T_12=rotacao('Z',q2)*translacao('Z',l2)*translacao('X',l3);
% T_23=rotacao('Y',q3)*translacao('X',l4);
% T_34=rotacao('Y',q4)*translacao('Z',excentricidade)*translacao('X',l6);
% T_45=rotacao('X',q5);
% T_56=rotacao('Y',q6);
% T_67=rotacao('X',q7);
% T_78=translacao('X',l7);
% T_02=T_01*T_12;
% T_03=T_02*T_23;
% T_04=T_03*T_34;
% T_05=T_04*T_45;
% T_06=T_05*T_56;
% T_07=T_06*T_67;



% Centro de massa de cada elo em rela��o ao referencial de cada elo. Esses
% dados foram retirados do programa SolidWorks
CG1=[-0.04937; 0; -0.47820; 1];
CG2=[-0.57627;0.03772;-0.11203;1];
CG3=[-0.82480;-0.24674;-0.00603;1];
CG4=[-0.90603;-0.01830;0.01229;1];
CG5=[-0.11647; 0.00143; 0.00000;1];
CG6=[-0.14129 ;0.02399; 0.00006;1];
CG7=[-0.18620;0.00038;0.00002;1]; %flange + ferramenta
CGtarefa=[-0.0047;-0.00130;0.01830;1];

% Jacobiano para o Bra�o
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

P0_Gtarefa=T_07*CGtarefa;
JL0_tarefa=jacobian(P0_Gtarefa(1:3),q);

% Derivada do centro de massa no referencial inercial
dP0_G1 = JL0_2*dq;
dP0_G2 = JL0_3*dq;
dP0_G3 = JL0_4*dq;
dP0_G4 = JL0_5*dq;
dP0_G5 = JL0_6*dq;
dP0_G6 = JL0_7*dq;
dP0_Gtarefa=JL0_tarefa*dq;

% Determinando o Jacobiano Angular
JA0_1 = [T_01(1:3,3), [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_2 = [T_01(1:3,3), T_02(1:3,2), [0 0 0]', [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_3 = [T_01(1:3,3), T_02(1:3,2), T_03(1:3,2), [0 0 0]', [0 0 0]', [0 0 0]'];
JA0_4 = [T_01(1:3,3), T_02(1:3,2), T_03(1:3,2), T_04(1:3,1), [0 0 0]', [0 0 0]'];
JA0_5 = [T_01(1:3,3), T_02(1:3,2), T_03(1:3,2), T_04(1:3,1), T_05(1:3,2), [0 0 0]'];
JA0_6 = [T_01(1:3,3), T_02(1:3,2), T_03(1:3,2), T_04(1:3,1), T_05(1:3,2), T_06(1:3,1)];

% Velocidade Angular do Trilho
% W0_1 = JA0_1*dq;

% Velocidade Angular do Bra�o 
W0_1 = JA0_1*dq;
W0_2 = JA0_2*dq;
W0_3 = JA0_3*dq;

% Velocidade Angular do Pulso
W0_4 = JA0_4*dq;
W0_5 = JA0_5*dq;
W0_6 = JA0_6*dq;

% Acelera��o da gravidade
gravidade = [0;0;-g];

% Energias
EC1 = ((m(2)*dP0_G1')*dP0_G1 + W0_1'*((T_02(1:3,1:3)*I2*T_02(1:3,1:3)')*W0_1))/2; % Bra�o
EC2 = ((m(3)*dP0_G2')*dP0_G2 + W0_2'*((T_03(1:3,1:3)*I3*T_03(1:3,1:3)')*W0_2))/2;
EC3 = ((m(4)*dP0_G3')*dP0_G3 + W0_3'*((T_04(1:3,1:3)*I4*T_04(1:3,1:3)')*W0_3))/2;
EC4 = ((m(5)*dP0_G4')*dP0_G4 + W0_4'*((T_05(1:3,1:3)*I5*T_05(1:3,1:3)')*W0_4))/2; % Pulso
EC5 = ((m(6)*dP0_G5')*dP0_G5 + W0_5'*((T_06(1:3,1:3)*I6*T_06(1:3,1:3)')*W0_5))/2;
EC6 = ((m(7)*dP0_G6')*dP0_G6 + W0_6'*((T_07(1:3,1:3)*I7*T_07(1:3,1:3)')*W0_6))/2;
ECtarefa=((m(8)*dP0_Gtarefa')*dP0_Gtarefa+W0_6'*((T_07(1:3,1:3)*Itarefa*T_07(1:3,1:3)')*W0_6))/2;
EP1 = (-m(2)*gravidade'*(P0_G1(1:3,1))); % Bra�o
EP2 = (-m(3)*gravidade'*(P0_G2(1:3,1))); 
EP3 = (-m(4)*gravidade'*(P0_G3(1:3,1)));
EP4 = (-m(5)*gravidade'*(P0_G4(1:3,1))); % Pulso
EP5 = (-m(6)*gravidade'*(P0_G5(1:3,1)));
EP6 = (-m(7)*gravidade'*(P0_G6(1:3,1)));
EPtarefa=(-m(8)*gravidade'*(P0_Gtarefa(1:3,1)));


L = (EC1 + EC2 + EC3 + EC4 + EC5 + EC6+ ECtarefa) - (EP1 + EP2 + EP3 + EP4 + EP5 + EP6+EPtarefa);

dL_dq1 = diff(L,q(1));
dL_dq2 = diff(L,q(2));
dL_dq3 = diff(L,q(3));
dL_dq4 = diff(L,q(4));
dL_dq5 = diff(L,q(5));
dL_dq6 = diff(L,q(6));



dL_ddq1 = diff(L,dq(1));
dL_ddq2 = diff(L,dq(2));
dL_ddq3 = diff(L,dq(3));
dL_ddq4 = diff(L,dq(4));
dL_ddq5 = diff(L,dq(5));
dL_ddq6 = diff(L,dq(6));


syms q1(t) q2(t) q3(t) q4(t) q5(t) q6(t)
q1 = q1(t);
q2 = q2(t);
q3 = q3(t);
q4 = q4(t);
q5 = q5(t);
q6 = q6(t);


dq1 = diff(q1,t);
dq2 = diff(q2,t);
dq3 = diff(q3,t);
dq4 = diff(q4,t);
dq5 = diff(q5,t);
dq6 = diff(q6,t);


ddq1 = diff(dq1,t);
ddq2 = diff(dq2,t);
ddq3 = diff(dq3,t);
ddq4 = diff(dq4,t);
ddq5 = diff(dq5,t);
ddq6 = diff(dq6,t);


q_L = [q1; q2; q3; q4; q5; q6];
dq_L = [dq1; dq2; dq3; dq4; dq5; dq6];
ddq_L = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6];

dL_ddq1 = subs(dL_ddq1);
dL_ddq2 = subs(dL_ddq2);
dL_ddq3 = subs(dL_ddq3);
dL_ddq4 = subs(dL_ddq4);
dL_ddq5 = subs(dL_ddq5);
dL_ddq6 = subs(dL_ddq6);


dL_ddq1_dt = diff(dL_ddq1,t);
dL_ddq2_dt = diff(dL_ddq2,t);
dL_ddq3_dt = diff(dL_ddq3,t);
dL_ddq4_dt = diff(dL_ddq4,t);
dL_ddq5_dt = diff(dL_ddq5,t);
dL_ddq6_dt = diff(dL_ddq6,t);


T1 = dL_ddq1_dt-dL_dq1 ;
T2 = dL_ddq2_dt - dL_dq2;
T3 = dL_ddq3_dt - dL_dq3;
T4 = dL_ddq4_dt - dL_dq4;
T5 = dL_ddq5_dt - dL_dq5;
T6 = dL_ddq6_dt - dL_dq6;


T1 = subs(T1 ,[q_L dq_L ddq_L],[q dq ddq]);
T2 = subs(T2 ,[q_L dq_L ddq_L],[q dq ddq]);
T3 = subs(T3 ,[q_L dq_L ddq_L],[q dq ddq]);
T4 = subs(T4 ,[q_L dq_L ddq_L],[q dq ddq]);
T5 = subs(T5 ,[q_L dq_L ddq_L],[q dq ddq]);
T6 = subs(T6 ,[q_L dq_L ddq_L],[q dq ddq]);


% M1
M11 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(1));
M12 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(2));
M13 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(3));
M14 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(4));
M15 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(5));
M16 = equationsToMatrix(subs(T1,[q dq ddq],[q dq ddq]),ddq(6));
M1 = [M11 M12 M13 M14 M15 M16];
% M2
M21 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(1));
M22 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(2));
M23 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(3));
M24 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(4));
M25 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(5));
M26 = equationsToMatrix(subs(T2,[q dq ddq],[q dq ddq]),ddq(6));
M2 = [M21 M22 M23 M24 M25 M26];
% M3
M31 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(1));
M32 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(2));
M33 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(3));
M34 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(4));
M35 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(5));
M36 = equationsToMatrix(subs(T3,[q dq ddq],[q dq ddq]),ddq(6));
M3 = [M31 M32 M33 M34 M35 M36];
% M4
M41 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(1));
M42 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(2));
M43 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(3));
M44 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(4));
M45 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(5));
M46 = equationsToMatrix(subs(T4,[q dq ddq],[q dq ddq]),ddq(6));
M4 = [M41 M42 M43 M44 M45 M46];
% M5
M51 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(1));
M52 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(2));
M53 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(3));
M54 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(4));
M55 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(5));
M56 = equationsToMatrix(subs(T5,[q dq ddq],[q dq ddq]),ddq(6));
M5 = [M51 M52 M53 M54 M55 M56];
% M6
M61 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(1));
M62 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(2));
M63 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(3));
M64 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(4));
M65 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(5));
M66 = equationsToMatrix(subs(T6,[q dq ddq],[q dq ddq]),ddq(6));
M6 = [M61 M62 M63 M64 M65 M66];

G1 = subs(T1 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G2 = subs(T2 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G3 = subs(T3 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G4 = subs(T4 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G5 = subs(T5 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));
G6 = subs(T6 ,[dq ddq],([0 0; 0 0; 0 0; 0 0; 0 0; 0 0]));


T1 = matlabFunction(T1);
T2 = matlabFunction(T2);
T3 = matlabFunction(T3);
T4 = matlabFunction(T4);
T5 = matlabFunction(T5);
T6 = matlabFunction(T6);
M1 = matlabFunction(M1);
M2 = matlabFunction(M2);
M3 = matlabFunction(M3);
M4 = matlabFunction(M4);
M5 = matlabFunction(M5);
M6 = matlabFunction(M6);
G1 = matlabFunction(G1);
G2 = matlabFunction(G2);
G3 = matlabFunction(G3);
G4 = matlabFunction(G4);
G5 = matlabFunction(G5);
G6 = matlabFunction(G6);
% 
% end
disp('Entrada na din�mica');
%[M, H, G,T]=Dinamica(q,dq,ddq,massa_robo,size_robo,Itarefa);

% disp('bloco M');
% matlabfunctionblock('novamontagemrobo/M',M);
% disp('bloco G');
% matlabfunctionblock('novamontagemrobo/G',G);
% disp('bloco T');
% matlabfunctionblock('novamontagemrobo/T',T);
