function [trajetoria_tempo]=substituir_tempo(trajetoria,tempo_final,delta_t)
pontos_traj=size(trajetoria,1);
indice_mat_posicoes=0;
trajetoria_tempo=zeros(floor(tempo_final/delta_t)+1,3); %o 0 só é plotado uma vez por isso +1, 
for j=1:pontos_traj                          %floor caso tempo_final/delta_t n�o seja inteiro
    for t=0:delta_t:tempo_final
         if t==tempo_final  %o tempo final do ponto da trajetoria coincide com o t0 do ponto seguinte
             if j~=pontos_traj %uma alternativa seria pular t0, entretanto para trajetorias com aceleracao inicial
                 continue      %como por exemplo a de 5ordem isso ocasionaria na perca de parte da reta nos pontos
             end               %em que a caneta estava no ar e tocou o papel para comecar uma nova reta                               
         end
         indice_mat_posicoes=indice_mat_posicoes+1;
        valor_temp=subs(trajetoria(j,:));
        trajetoria_tempo(indice_mat_posicoes,:)=valor_temp;
    end
end