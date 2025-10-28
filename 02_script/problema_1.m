clear all; clc; close all; 
fprintf("\n\nProblema 1 \n\n");

[t c1] = ode45(@(t, c) prob1(t, c, 2, 1, 1.25),[0 8],[1 0 0]);

plot_and_save(t, c1, "figura1_problema_1a.png");

[max_B1 pos_B1] = max(c1(:,2));
tpo_max_B1 = t(pos_B1);
fprintf("La concentracion maxima de B es %.4f mol/L y se alcanza a las %.4f horas.\n\n", max_B1, tpo_max_B1);

[t c2] = ode45(@(t, c) prob1(t, c, 2, 1, 1.5),[0 8],[1 0 0]);

plot_and_save(t, c2, "figura1_problema_1b.png");

[max_B2 pos_B2] = max(c2(:,2));
tpo_max_B2 = t(pos_B2);
fprintf("La concentracion maxima de B es %.4f mol/L y se alcanza a las %.4f horas.\n\n", max_B2, tpo_max_B2);

Nombres_Columnas = {'k2_1_25 ', 'k2_1_5' , 'delta_var'};
Nombres_Filas = {'valor maximo B (mol/L)', 'tiempo donde se alcanza el maximo (h)'};
T = table([max_B1 ; tpo_max_B1], [max_B2 ; tpo_max_B2 ] , [100*(max_B2 - max_B1)/max_B1 ; 100*(tpo_max_B2 - tpo_max_B1)/tpo_max_B1 ], 'VariableNames', Nombres_Columnas, 'RowNames', Nombres_Filas);
disp(T);

function dc = prob1(t,c, k1f, k1r, k2)

    % k1f 
    % k1r 
    % k2 

    dc = zeros(3,1);
    dc(1) = -k1f*c(1) + k1r*c(2);
    dc(2) = k1f*c(1) - (k1r + k2)*c(2);
    dc(3) = k2*c(2);
end

function plot_and_save(t, c, filename)
    fig = figure('Visible', 'off');
    plot(t,c(:,1),'k','LineWidth',1.5);
    hold on
    grid on
    plot(t,c(:,2),'k','LineWidth',1.5);
    plot(t,c(:,3),'k--','LineWidth',1.5);
    hold off
    axis([0 5 0 1]);
    xlabel('Tiempo (h)');
    ylabel('Concentracion (mol/L)');
    legend('A','B','C');
    saveas(fig, filename);
    close(fig);
end
