clear all, clc, close all
%% 
%   Rotina para cálculo da cinemática direta e inversa de um manipulador
%   antropomórfico alterada para robô da tarefa do estoque
%
%   Autor:  Luciano
%   Data:   09-11-2020
%
%   Revisão: Grupo
%   Data:    06-01-2021

%% Definição de variáveis

syms q1 q2 q3 q4 q5 q6

q = [q1; q2; q3; q4; q5; q6];

syms dq1 dq2 dq3 dq4 dq5 dq6

dq = [dq1; dq2; dq3; dq4; dq5; dq6];

syms d2q1 d2q2 d2q3 d2q4 d2q5 d2q6

d2q = [d2q1; d2q2; d2q3; d2q4; d2q5; d2q6];

syms L1 L2 L3 L4 

syms E1 E2                          %excentricidades e1=750; e2=-41

syms xf yf zf

syms xG1 yG1 zG1 xG2 yG2 zG2 xG3 yG3 zG3 xG4 yG4 zG4 xG5 yG5 zG5 xG6 yG6 zG6 xGf yGf zGf

syms xC yC zC qxC qyC qzC

syms x7_C y7_C z7_C alfa beta phi

%% Análise do Elo 1
R1_0 = [cos(q1) -sin(q1) 0 0; sin(q1) cos(q1) 0 0; 0 0 1 0; 0 0 0 1];
D1_0 = [1 0 0 E1; 0 1 0 0; 0 0 1 L1; 0 0 0 1];
T1_0 = R1_0*D1_0;
p1_0 = T1_0(1:3,4);
JL1_0 = jacobian(p1_0,q);
v1_0 = JL1_0*dq;
dJL1_0 = [jacobian(JL1_0(:,1),q)*dq jacobian(JL1_0(:,2),q)*dq jacobian(JL1_0(:,3),q)*dq jacobian(JL1_0(:,4),q)*dq jacobian(JL1_0(:,5),q)*dq jacobian(JL1_0(:,6),q)*dq];
a1_0 = JL1_0*d2q + dJL1_0*dq;

JA1_0 = [T1_0(1:3,3) zeros(3,5)];
w1_0 = JA1_0*dq;
dJA1_0 = [jacobian(JA1_0(:,1),q)*dq jacobian(JA1_0(:,2),q)*dq jacobian(JA1_0(:,3),q)*dq jacobian(JA1_0(:,4),q)*dq jacobian(JA1_0(:,5),q)*dq jacobian(JA1_0(:,6),q)*dq];
dw1_0 = JA1_0*d2q + dJA1_0*dq;

pG1_1 = [xG1; yG1; zG1; 1];
aux = T1_0*pG1_1;
pG1_0 = aux(1:3,1);
JLG1_0 = jacobian(pG1_0,q);
vG1_0 = JLG1_0*dq;
dJLG1_0 = [jacobian(JLG1_0(:,1),q)*dq jacobian(JLG1_0(:,2),q)*dq jacobian(JLG1_0(:,3),q)*dq jacobian(JLG1_0(:,4),q)*dq jacobian(JLG1_0(:,5),q)*dq jacobian(JLG1_0(:,6),q)*dq];
aG1_0 = JLG1_0*d2q + dJLG1_0*dq;

%% Análise do Elo 2
R2_1 = [cos(q2) 0 sin(q2) 0; 0 1 0 0; -sin(q2) 0 cos(q2) 0; 0 0 0 1];
D2_1 = [1 0 0 L2; 0 1 0 0; 0 0 1 E2; 0 0 0 1];
T2_1 = R2_1*D2_1;
T2_0 = T1_0*T2_1;
p2_0 = T2_0(1:3,4);
JL2_0 = jacobian(p2_0,q);
v2_0 = JL2_0*dq;
dJL2_0 = [jacobian(JL2_0(:,1),q)*dq jacobian(JL2_0(:,2),q)*dq jacobian(JL2_0(:,3),q)*dq jacobian(JL2_0(:,4),q)*dq jacobian(JL2_0(:,5),q)*dq jacobian(JL2_0(:,6),q)*dq];
a2_0 = simplify(JL2_0*d2q + dJL2_0*dq);

JA2_0 = [T1_0(1:3,3) T2_0(1:3,2) zeros(3,4)];
w2_0 = JA2_0*dq;
dJA2_0 = [jacobian(JA2_0(:,1),q)*dq jacobian(JA2_0(:,2),q)*dq jacobian(JA2_0(:,3),q)*dq jacobian(JA2_0(:,4),q)*dq jacobian(JA2_0(:,5),q)*dq jacobian(JA2_0(:,6),q)*dq];
dw2_0 = JA2_0*d2q + dJA2_0*dq;

pG2_2 = [xG2; yG2; zG2; 1];
aux = T2_0*pG2_2;
pG2_0 = aux(1:3,1);
JLG2_0 = jacobian(pG2_0,q);
vG2_0 = JLG2_0*dq;
dJLG2_0 = [jacobian(JLG2_0(:,1),q)*dq jacobian(JLG2_0(:,2),q)*dq jacobian(JLG2_0(:,3),q)*dq jacobian(JLG2_0(:,4),q)*dq jacobian(JLG2_0(:,5),q)*dq jacobian(JLG2_0(:,6),q)*dq];
aG2_0 = JLG2_0*d2q + dJLG2_0*dq;

%% Análise do Elo 3
R3_2 = [cos(q3) 0 sin(q3) 0; 0 1 0 0; -sin(q3) 0 cos(q3) 0; 0 0 0 1];
D3_2 = [1 0 0 L3; 0 1 0 0; 0 0 1 0; 0 0 0 1];
T3_2 = R3_2*D3_2;
T3_0 = simplify(T2_0*T3_2);
p3_0 = T3_0(1:3,4);
JL3_0 = jacobian(p3_0,q);
v3_0 = JL3_0*dq;
dJL3_0 = [jacobian(JL3_0(:,1),q)*dq jacobian(JL3_0(:,2),q)*dq jacobian(JL3_0(:,3),q)*dq jacobian(JL3_0(:,4),q)*dq jacobian(JL3_0(:,5),q)*dq jacobian(JL3_0(:,6),q)*dq];
a3_0 = simplify(JL3_0*d2q + dJL3_0*dq);

JA3_0 = [T1_0(1:3,3) T2_0(1:3,2) T3_0(1:3,2) zeros(3,3)];
w3_0 = JA3_0*dq;
dJA3_0 = [jacobian(JA3_0(:,1),q)*dq jacobian(JA3_0(:,2),q)*dq jacobian(JA3_0(:,3),q)*dq jacobian(JA3_0(:,4),q)*dq jacobian(JA3_0(:,5),q)*dq jacobian(JA3_0(:,6),q)*dq];
dw3_0 = JA3_0*d2q + dJA3_0*dq;

pG3_3 = [xG3; yG3; zG3; 1];
aux = T3_0*pG3_3;
pG3_0 = aux(1:3,1);
JLG3_0 = jacobian(pG3_0,q);
vG3_0 = JLG3_0*dq;
dJLG3_0 = [jacobian(JLG3_0(:,1),q)*dq jacobian(JLG3_0(:,2),q)*dq jacobian(JLG3_0(:,3),q)*dq jacobian(JLG3_0(:,4),q)*dq jacobian(JLG3_0(:,5),q)*dq jacobian(JLG3_0(:,6),q)*dq];
aG3_0 = JLG3_0*d2q + dJLG3_0*dq;

%% Análiose do Elo 4
R4_3 = [1 0 0 0; 0 cos(q4) -sin(q4) 0; 0 sin(q4) cos(q4) 0; 0 0 0 1];
T4_3 = R4_3;
T4_0 = simplify(T3_0*T4_3);
p4_0 = p3_0;
JL4_0 = JL3_0;
v4_0 = v3_0;
dJL4_0 = dJL3_0;
a4_0 = a3_0;

JA4_0 = [T1_0(1:3,3) T2_0(1:3,2) T3_0(1:3,2) T4_0(1:3,1) zeros(3,2)];
w4_0 = JA4_0*dq;
dJA4_0 = [jacobian(JA4_0(:,1),q)*dq jacobian(JA4_0(:,2),q)*dq jacobian(JA4_0(:,3),q)*dq jacobian(JA4_0(:,4),q)*dq jacobian(JA4_0(:,5),q)*dq jacobian(JA4_0(:,6),q)*dq];
dw4_0 = JA4_0*d2q + dJA4_0*dq;

pG4_4 = [xG4; yG4; zG4; 1];
aux = T4_0*pG4_4;
pG4_0 = aux(1:3,1);
JLG4_0 = jacobian(pG4_0,q);
vG4_0 = JLG4_0*dq;
dJLG4_0 = [jacobian(JLG4_0(:,1),q)*dq jacobian(JLG4_0(:,2),q)*dq jacobian(JLG4_0(:,3),q)*dq jacobian(JLG4_0(:,4),q)*dq jacobian(JLG4_0(:,5),q)*dq jacobian(JLG4_0(:,6),q)*dq];
aG4_0 = JLG4_0*d2q + dJLG4_0*dq;

%% Análise do Elo 5
R5_4 = [cos(q5) 0 sin(q5) 0; 0 1 0 0; -sin(q5) 0 cos(q5) 0; 0 0 0 1];
T5_4 = R5_4;
T5_0 = simplify(T4_0*T5_4);
p5_0 = p4_0;
JL5_0 = JL4_0;
v5_0 = v4_0;
dJL5_0 = dJL4_0;
a5_0 = a4_0;

JA5_0 = [T1_0(1:3,3) T2_0(1:3,2) T3_0(1:3,2) T4_0(1:3,1) T5_0(1:3,2) zeros(3,1)];
w5_0 = JA5_0*dq;
dJA5_0 = [jacobian(JA5_0(:,1),q)*dq jacobian(JA5_0(:,2),q)*dq jacobian(JA5_0(:,3),q)*dq jacobian(JA5_0(:,4),q)*dq jacobian(JA5_0(:,5),q)*dq jacobian(JA5_0(:,6),q)*dq];
dw5_0 = JA5_0*d2q + dJA5_0*dq;

pG5_5 = [xG5; yG5; zG5; 1];
aux = T5_0*pG5_5;
pG5_0 = aux(1:3,1);
JLG5_0 = jacobian(pG5_0,q);
vG5_0 = JLG5_0*dq;
dJLG5_0 = [jacobian(JLG5_0(:,1),q)*dq jacobian(JLG5_0(:,2),q)*dq jacobian(JLG5_0(:,3),q)*dq jacobian(JLG5_0(:,4),q)*dq jacobian(JLG5_0(:,5),q)*dq jacobian(JLG5_0(:,6),q)*dq];
aG5_0 = JLG5_0*d2q + dJLG5_0*dq;

%% Análise do Elo 6
R6_5 = [1 0 0 0; 0 cos(q6) -sin(q6) 0; 0 sin(q6) cos(q6) 0; 0 0 0 1];
T6_5 = R6_5;
T6_0 = simplify(T5_0*T6_5);
T6_3 = simplify(T4_3*T5_4*T6_5);
p6_0 = p5_0;
JL6_0 = JL5_0;
v6_0 = v5_0;
dJL6_0 = dJL5_0;
a6_0 = a5_0;

JA6_0 = [T1_0(1:3,3) T2_0(1:3,2) T3_0(1:3,2) T4_0(1:3,1) T5_0(1:3,2) T6_0(1:3,1)];
w6_0 = JA6_0*dq;
dJA6_0 = [jacobian(JA6_0(:,1),q)*dq jacobian(JA6_0(:,2),q)*dq jacobian(JA6_0(:,3),q)*dq jacobian(JA6_0(:,4),q)*dq jacobian(JA6_0(:,5),q)*dq jacobian(JA6_0(:,6),q)*dq];
dw6_0 = JA6_0*d2q + dJA6_0*dq;

pG6_6 = [xG6; yG6; zG6; 1];
aux = T6_0*pG6_6;
pG6_0 = aux(1:3,1);
JLG6_0 = jacobian(pG6_0,q);
vG6_0 = JLG6_0*dq;
dJLG6_0 = [jacobian(JLG6_0(:,1),q)*dq jacobian(JLG6_0(:,2),q)*dq jacobian(JLG6_0(:,3),q)*dq jacobian(JLG6_0(:,4),q)*dq jacobian(JLG6_0(:,5),q)*dq jacobian(JLG6_0(:,6),q)*dq];
aG6_0 = JLG6_0*d2q + dJLG6_0*dq;

%% Análise do Ferramenta
D7_6 = [1 0 0 L4+xf; 0 1 0 yf; 0 0 1 zf; 0 0 0 1];
T7_6 = D7_6;
T7_0 = simplify(T6_0*T7_6);
p7_0 = T7_0(1:3,4);
JL7_0 = jacobian(p7_0,q);
v7_0 = simplify(JL7_0*dq);
dJL7_0 = [jacobian(JL7_0(:,1),q)*dq jacobian(JL7_0(:,2),q)*dq jacobian(JL7_0(:,3),q)*dq jacobian(JL7_0(:,4),q)*dq jacobian(JL7_0(:,5),q)*dq jacobian(JL7_0(:,6),q)*dq];
a7_0 = simplify(JL7_0*d2q + dJL7_0*dq);

w7_0 = w6_0;
dw7_0 = dw6_0;

pG7_7 = [xGf; yGf; zGf; 1];
aux = T7_0*pG7_7;
pG7_0 = aux(1:3,1);
JLG7_0 = jacobian(pG7_0,q);
vG7_0 = JLG7_0*dq;
dJLG7_0 = [jacobian(JLG7_0(:,1),q)*dq jacobian(JLG7_0(:,2),q)*dq jacobian(JLG7_0(:,3),q)*dq jacobian(JLG7_0(:,4),q)*dq jacobian(JLG7_0(:,5),q)*dq jacobian(JLG7_0(:,6),q)*dq];
aG7_0 = JLG7_0*d2q + dJLG7_0*dq;

%% Parametrização da tarefa x ferramenta

% Definição do ambiente de trabalho
DC_0 = [1 0 0 xC; 0 1 0 yC; 0 0 1 zC; 0 0 0 1];
Rz_qzC = [cos(qzC) -sin(qzC) 0 0; sin(qzC) cos(qzC) 0 0; 0 0 1 0; 0 0 0 1];
Ry_qyC = [cos(qyC) 0 sin(qyC) 0; 0 1 0 0; -sin(qyC) 0 cos(qyC) 0; 0 0 0 1];
Rx_qxC = [1 0 0 0; 0 cos(qxC) -sin(qxC) 0; 0 sin(qxC) cos(qxC) 0; 0 0 0 1];
RC_0 = Rz_qzC*Ry_qyC*Rx_qxC;
TC_0 = DC_0*RC_0;
% 
% Tarefa
D7_C = [1 0 0 x7_C; 0 1 0 y7_C; 0 0 1 z7_C; 0 0 0 1];
Rz_alfa = [cos(alfa) -sin(alfa) 0 0; sin(alfa) cos(alfa) 0 0; 0 0 1 0; 0 0 0 1];
Ry_beta = [cos(beta) 0 sin(beta) 0; 0 1 0 0; -sin(beta) 0 cos(beta) 0; 0 0 0 1];
Rx_phi = [1 0 0 0; 0 cos(phi) -sin(phi) 0; 0 sin(phi) cos(phi) 0; 0 0 0 1];
R7_C = Rz_alfa*Ry_beta*Rx_phi;
T7_C = D7_C*R7_C;

% Juntar tarefa e ambiente
T7_0f = TC_0*T7_C;
R7_0f = RC_0*R7_C;

% % xC = 1; yC = 1; zC = 0.7; qxC = 0; qyC = 0; qzC = 0; 
% % x7_C = 0; y7_C = 0; z7_C = 0; alfa = 0; beta = pi/2; phi = 0;
% % aux = eval(T7_0f)
% 
%% Definição da posição do pulso no referencial da base do robô - {0}
p7_0f = T7_0f(1:3,4);
p3_7_7f = [-L4-xf; -yf; -zf; 1];
aux = R7_0f*p3_7_7f;
p3_0f = p7_0f + aux(1:3,1);
% 
% return
% %% Simulações


%% INVERSA || PARTE 1 q possiveis

%%%% Variaveis a receber %%%%%%%%%
%P3_0f <- numerico
%T3_0 <- ou criar dentro 
%T6_3 <- ou criar dentro
%T7_0f <- ou criar dentro ||  Problema se manter L4+E2 no T7_6  
%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%q_possiveis=zeros(10,6); Testar se necessita ou pode apagar
X=p3_0f(1);
Y=p3_0f(2);
Z=p3_0f(3);
ind_k=0;
q1=atan2(Y,X);
if isreal(q1)==0
    error("Região fora da area de trabalho");
end
%if q1>deg2rad(185) || q1<deg2rad(-185)
%    error("Região fora da area de trabalho");
%end

% if sin(q1)<0.1 && sin(q1)>-0.1
%     aux=X/cos(q1);
%     
% else
%     aux=Y/sin(q1);
%     
% end

aux=piecewise(sin(q1)<0.1,X/cos(q1),sin(q1)>-0.1,X/cos(q1),Y/sin(q1));


%- 2*L1*zC  + E2^2+ L2^2 + aux^2 + E1^2 + L1^2 + zC^2 + 2*L2*zC*sin(q2) - (2*E1)*aux + 2*E1*L2*cos(q2) + 2*E2*L1*cos(q2) + 2*E1*E2*sin(q2) - 2*L1*L2*sin(q2) - 2*E2*zC*cos(q2) - (2*L2*cos(q2)*aux) - (2*E2*sin(q2))*aux-L3^2

           



A_q2=2*E1*L2+ 2*E2*L1- 2*E2*zC -(2*L2*aux);
B_q2=2*L2*zC+ 2*E1*E2 - 2*L1*L2- (2*E2)*aux;
C_q2=- 2*L1*zC+ E2^2+ L2^2 + aux^2 + E1^2 + L1^2 + zC^2- (2*E1)*aux-L3^2;

delta_q2= 4*A_q2^2*C_q2^2-4*(A_q2^2+B_q2^2)*(C_q2^2-B_q2^2);
% if delta_q2<0
%     error("Fora da região de trabalho");
% end

 
 q2a=zeros(1,4);
 argumanq2a1=((-2*A_q2*C_q2+sqrt(delta_q2))/(2*(A_q2^2+B_q2^2)));
 argumanq2a1 = max(min(argumanq2a1,1),-1);
 
 argumanq2a2=((-2*A_q2*C_q2-sqrt(delta_q2))/(2*(A_q2^2+B_q2^2)));
 argumanq2a2 = max(min(argumanq2a2,1),-1);
 
 q2a(1)=acos(argumanq2a1);
 q2a(2)=acos(argumanq2a2);
 q2a(3)=-q2a(1);
 q2a(4)=-q2a(2);

%%%%%% PARTE NOVA %%%%%%%%%%%%%%%%%%%%
q3a=zeros(1,2);
for i=1:1:4 %valores q2
    q2=q2a(i);
    if q2>deg2rad(-120) && q2<deg2rad(70) && isreal(q2)==1
        %colocar bloco de restrição de altura se a tarefa requirir | if altura > blablabla continue%
        soma_dos_angulos=acos( (aux-E1-L2*cos(q2)-E2*sin(q2)) /L3 ); 
        q3a(1)=soma_dos_angulos-q2;
        q3a(2)=-(soma_dos_angulos)-q2;
        for i2=1:1:2 %valores q3
            q3=q3a(i2);
            if q3>deg2rad(-120) && q3<deg2rad(155) && isreal(q3)==1
                T3_0num=subs(T3_0);
                if norm(T3_0num(1,4)-X)<0.0001 && norm(T3_0num(2,4)-Y)<0.0001 && norm(T3_0num(3,4)-Z) <0.0001 %opcional
                    T6_0fnum=subs(T7_0f)*inv(T7_6); %%%% SE VIER DE OUTRO BLOCO VIRA NUMERICA
                    T6_3num=inv(T3_0num)*T6_0fnum; %%% COMPARAR COM T6_3
                    q4=atan2(T6_3num(2,1),-T6_3num(3,1)); %% q4 estara sempre no intervalo [-350;350] graus devido a funcao atan2
                    q6=atan2(T6_3num(1,2),T6_3num(1,3)); %% iden a q4 sobre o intervalo
                    q5a=zeros(1,2);
                    q5a(1)=acos(T6_3num(1,1));
                    q5a(2)=-q5a(1);
                    if isreal(q4) == 1 && isreal(q6) == 1 %q4 e q6 reais
                        for i3=1:1:2 %valores q5
                            if q5a(i3)>deg2rad(-125) && q5a(i3)<deg2rad(125) && isreal(q5a) ==1 %q5 valido
                                q5=q5a(i3);
                                T7_0num=subs(T7_0);
                                if norm(T7_0num(1,4)-p7_0f(1))<0.0001 && norm(T7_0num(2,4)-p7_0f(2))<0.0001 && norm(T7_0num(3,4)-p7_0f(3))<0.0001
                                    ind_k=ind_k+1;
                                    q_possiveis(ind_k,1)=q1;
                                    q_possiveis(ind_k,2)=q2;
                                    q_possiveis(ind_k,3)=q3;
                                    q_possiveis(ind_k,4)=q4;
                                    q_possiveis(ind_k,5)=q5;
                                    q_possiveis(ind_k,6)=q6;
                                end %pos braço valida
                            end %if q5 valido
                        end %valores q5                
                    end %q4 e q6 validos
                end %checar posicao braco
            end %q3 validos
        end %valores q3
    end % q2 validos
end %valores q2

% PARTE 2 , pegar combinacao de q mais proxima da anterior
coder.extrinsic('evalin')
coder.extrinsic('assigin')
menor_desloc=0;
if t==0
    Qs=[0 0 0 0 0 0]; %%%POSIÇÃO INICIAL ESTÁ CERTA ?
else
    Qs=evalin('caller','Qs');
end
valor_q_inv=zeros(1,6);
%fprintf("Realizando inversa do ponto %d\n",ind_ponto)
%[q_possiveis_calc,ind_k]=inversa_q_possiveis(mat71);
if ind_k<1
   error("ERRO: Naoo foi possivel realizar uma movimentacao para o ponto")
end
for z=1:ind_k
    vetor_deslocamento_ang=q_possiveis(z,:)-Qs;
    if z==1
        menor_desloc=norm(vetor_deslocamento_ang);
        valor_q_inv=q_possiveis_calc(z,:);
    else
        if norm(vetor_deslocamento_ang) < menor_desloc
            valor_q_inv=q_possiveis_calc(z,:);
            menor_desloc=norm(vetor_deslocamento_ang);
        end
    end
end
assigin('base','Qs',valor_q_inv);
 
 
 
 
 
 
 
 matlabFunctionBlock('montagemrobo2/inversa',valor_q_inversa)
 
 
 
 
 
 
 



    
    
