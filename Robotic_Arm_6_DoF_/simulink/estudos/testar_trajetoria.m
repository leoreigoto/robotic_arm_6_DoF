clc
clear all
close all
centro_z=907.64;
raio_traj=2050.15;
delta_t=6;
tempo_plot=1;
pos_final_traj1=-2040.39;
pos_final_traj2=-1822.22;
ang_final_traj1=acos(pos_final_traj1/raio_traj);
ang_final_traj2=acos(pos_final_traj2/raio_traj);
if ang_final_traj1<0
    ang_final_traj1=ang_final_traj1*(-1);
end
if ang_final_traj2<0
    ang_final_traj2=ang_final_traj2*(-1);
end
trajetoria1=[0;ang_final_traj1];
trajetoria2=[0;ang_final_traj2];
[pos1,vel1,acel1]=traj_cubica1D(trajetoria1,delta_t,3,2,1,raio_traj,raio_traj,raio_traj,0,0,centro_z);
pos_tempo1=substituir_tempo(pos1,delta_t,tempo_plot);
[pos2,vel2,acel2]=traj_cubica1D(trajetoria2,delta_t,3,2,1,raio_traj,raio_traj,raio_traj,0,0,centro_z);
pos_tempo2=substituir_tempo(pos2,delta_t,tempo_plot);
for t=1:size(pos_tempo1,1)
    hold on
    xlim([-1500 2500])
    ylim([-2500 2500])
    zlim([0 1500])
    plot3(0,0,0,'*b');
    plot3(pos_tempo1(t,1),pos_tempo1(t,2),pos_tempo1(t,3),'*r');
    pause(0.1)
end
pause(0.2)
figure()
for t=1:size(pos_tempo2,1)
    hold on
    xlim([-1500 2500])
    ylim([-2500 2500])
    zlim([0 1500])
    plot3(0,0,0,'*b');
    plot3(pos_tempo2(t,1),pos_tempo2(t,2),pos_tempo2(t,3),'*r');
    pause(0.1)
end