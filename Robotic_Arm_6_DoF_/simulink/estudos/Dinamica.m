 function [M, H, G,T] = Dinamica(q, dq, ddq,massa_robo,size_robo,Itarefa)
%     M1, M2, M3, M4, M5, M6, M7,...
%     T1, T2, T3, T4, T5, T6, T7,...
%     G1, G2, G3, G4, G5, G6, G7,...
%     ...
%     )

coder.extrinsic('evalin');

M1 = evalin('caller','M1');
M2 = evalin('caller','M2');
M3 = evalin('caller','M3');
M4 = evalin('caller','M4');
M5 = evalin('caller','M5');
M6 = evalin('caller','M6');


T1 = evalin('caller','T1');
T2 = evalin('caller','T2');
T3 = evalin('caller','T3');
T4 = evalin('caller','T4');
T5 = evalin('caller','T5');
T6 = evalin('caller','T6');

 
G1 = evalin('caller','G1'); 
G2 = evalin('caller','G2');
G3 = evalin('caller','G3');
G4 = evalin('caller','G4');
G5 = evalin('caller','G5');
G6 = evalin('caller','G6');





ddq1 = ddq(1);
ddq2 = ddq(2);
ddq3 = ddq(3);
ddq4 = ddq(4);
ddq5 = ddq(5);
ddq6 = ddq(6);

dq1 = dq(1);
dq2 = dq(2);
dq3 = dq(3);
dq4 = dq(4);
dq5 = dq(5);
dq6 = dq(6);

q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);
q5 = q(5);
q6 = q(6);


m1=massa_robo(1);
m2=massa_robo(2);
m3=massa_robo(3);
m4=massa_robo(4);
m5=massa_robo(5);
m6=massa_robo(6);
m7=massa_robo(7);
mtarefa=massa_robo(8);


L1=size_robo(1);
L2=size_robo(2);
L3=size_robo(3);
L4=size_robo(4);
E1=size_robo(5);
E2=size_robo(6);
xf=size_robo(7);
yf=size_robo(8);
zf=size_robo(9);

Itarefaxx=Itarefa(1,1);
Itarefaxy=Itarefa(1,2);
Itarefaxz=Itarefa(1,3);
Itarefayx=Itarefa(2,1);
Itarefayy=Itarefa(2,2);
Itarefayz=Itarefa(2,3);
Itarefazx=Itarefa(3,1);
Itarefazy=Itarefa(3,2);
Itarefazz=Itarefa(3,3);



T1 = T1(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,m2,m3,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
%T2 = T2(Itarefa,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,mtarefa,q1,q2,q3,q4,q5,q6);
T2 = T2(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,m2,m3,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
T3 = T3(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,m3,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
T4 = T4(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
T5 = T5(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
T6 = T6(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq1,dq2,dq3,dq4,dq5,dq6,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);

M1 = M1(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,m2,m3,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
M2 = M2(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,m2,m3,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
M3 = M3(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,m3,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
M4 = M4(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,m4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
M5 = M5(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,m5,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
M6 = M6(E1,E2,Itarefaxx,Itarefaxy,Itarefayx,Itarefaxz,Itarefayy,Itarefazx,Itarefayz,Itarefazy,Itarefazz,L2,L3,L4,m6,m7,mtarefa,q1,q2,q3,q4,q5,q6,xf,yf,zf);
G1 = G1();
G2 = G2(E2,L2,L3,L4,m2,m3,m4,m5,m6,m7,mtarefa,q2,q3,q4,q5,q6,xf,yf,zf);
G3 = G3(E2,L3,L4,m3,m4,m5,m6,m7,mtarefa,q2,q3,q4,q5,q6,xf,yf,zf);
G4 = G4(L4,m4,m5,m6,m7,mtarefa,q2,q3,q4,q5,q6,xf,yf,zf);
G5 = G5(L4,m5,m6,m7,mtarefa,q2,q3,q4,q5,q6,xf,yf,zf);
G6 = G6(m6,m7,mtarefa,q2,q3,q4,q5,q6,yf,zf);



M = [M1; M2; M3; M4; M5; M6];
G = [G1; G2; G3; G4; G5; G6];
T = [T1; T2; T3; T4; T5; T6];
H = T-M*ddq-G;

end