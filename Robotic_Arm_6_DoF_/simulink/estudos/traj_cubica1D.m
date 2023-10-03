function [traj_pos_funcao,traj_vel_funcao,traj_acel_funcao]= traj_cubica1D(Trajetoria_funcao ...
    ,temp_desloc,aux_X,aux_Y,aux_Z,r_X,r_Y,r_Z,centro_X_func,centro_Y_func,centro_Z_func)

syms t cub_ang
tamanho_traj=size(Trajetoria_funcao,1);
temp_desloc=temp_desloc/(tamanho_traj-1);
switch aux_X
    case 1
       aux_X=0;
    case 2
       aux_X=cos(cub_ang);
    case 3
       aux_X=sin(cub_ang);
end     
     
switch aux_Y
    case 1
       aux_Y=0;
    case 2
       aux_Y=cos(cub_ang);
    case 3
       aux_Y=sin(cub_ang);  
end

switch aux_Z
    case 1
       aux_Z=0;
    case 2
       aux_Z=cos(cub_ang);
    case 3
       aux_Z=sin(cub_ang);
end

for w=1:tamanho_traj-1
    ang_ini=Trajetoria_funcao(w,1);
    ang_final=Trajetoria_funcao(w+1,1);
    cub_a_0=ang_ini;
    cub_a_1=0;
    cub_a_2=3*(ang_final-ang_ini)*(1/temp_desloc^2);
    cub_a_3=-2*(ang_final-ang_ini)*(1/temp_desloc^3);   
    cub_ang=cub_a_0+cub_a_1*t+cub_a_2*t^2+cub_a_3*t^3;
    
    traj_pos_funcao(w,1)= r_X*subs(aux_X)+ centro_X_func;
    traj_vel_funcao(w,1)=diff(traj_pos_funcao(w,1),t);
    traj_acel_funcao(w,1)=diff(traj_vel_funcao(w,1),t);
    
    traj_pos_funcao(w,2)=r_Y* subs(aux_Y) + centro_Y_func;
    traj_vel_funcao(w,2)=diff(traj_pos_funcao(w,2),t);
    traj_acel_funcao(w,2)=diff(traj_vel_funcao(w,2),t);
    
    traj_pos_funcao(w,3)= r_Z* subs(aux_Z) + centro_Z_func;
    traj_vel_funcao(w,3)=diff(traj_pos_funcao(w,3),t);
    traj_acel_funcao(w,3)=diff(traj_vel_funcao(w,3),t);
end

end