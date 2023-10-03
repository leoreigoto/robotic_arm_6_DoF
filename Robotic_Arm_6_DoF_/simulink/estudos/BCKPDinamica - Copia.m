 function [M, H, G,T] = Dinamica(q, dq, ddq)
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
M7 = evalin('caller','M7');

T1 = evalin('caller','T1');
T2 = evalin('caller','T2');
T3 = evalin('caller','T3');
T4 = evalin('caller','T4');
T5 = evalin('caller','T5');
T6 = evalin('caller','T6');
T7 = evalin('caller','T7');
 
G1 = evalin('caller','G1'); 
G2 = evalin('caller','G2');
G3 = evalin('caller','G3');
G4 = evalin('caller','G4');
G5 = evalin('caller','G5');
G6 = evalin('caller','G6');
G7 = evalin('caller','G7');




ddq1 = ddq(1);
ddq2 = ddq(2);
ddq3 = ddq(3);
ddq4 = ddq(4);
ddq5 = ddq(5);
ddq6 = ddq(6);
ddq7 = ddq(7);
dq1 = dq(1);
dq2 = dq(2);
dq3 = dq(3);
dq4 = dq(4);
dq5 = dq(5);
dq6 = dq(6);
dq7 = dq(7);
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);
q5 = q(5);
q6 = q(6);
q7 = q(7);

T1 = T1(ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,dq2,dq3,dq4,dq5,dq6,q2,q3,q4,q5,q6);
T2 = T2(ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,q2,q3,q4,q5,q6,q7);
T3 = T3(ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,q2,q3,q4,q5,q6,q7);
T4 = T4(ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,q2,q3,q4,q5,q6,q7);
T5 = T5(ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,q2,q3,q4,q5,q6,q7);
T6 = T6(ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,q2,q3,q4,q5,q6,q7);
T7 = T7(ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq2,dq3,dq4,dq5,dq6,dq7,q2,q3,q4,q5,q6,q7);
M1 = M1(q2,q3,q4,q5,q6);
M2 = M2(q2,q3,q4,q5,q6,q7);
M3 = M3(q2,q3,q4,q5,q6,q7);
M4 = M4(q2,q3,q4,q5,q6,q7);
M5 = M5(q2,q3,q4,q5,q6,q7);
M6 = M6(q2,q3,q4,q5,q6,q7);
M7 = M7(q2,q3,q4,q5,q6,q7);
G1 = G1();
G2 = G2();
G3 = G3(q3,q4,q5,q6);
G4 = G4(q3,q4,q5,q6);
G5 = G5(q3,q4,q5,q6);
G6 = G6(q3,q4,q5,q6);
G7 = G7();


M = [M1; M2; M3; M4; M5; M6; M7];
G = [G1; G2; G3; G4; G5; G6; G7];
T = [T1; T2; T3; T4; T5; T6; T7];
H = T-M*ddq-G;

end