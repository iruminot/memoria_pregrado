clear all; clc; close all;
fprintf("\n\nProblema 2 \n\n");
options = optimset('TolX',1e-8,'TolFun',1e-8,'Display','off');
c = fsolve(@prob2,[1 1 1],optimset,0.5)

i = 1;
c0 = [1 1 1];
tau = 0.1:0.5:20;
for taui = tau
    c(i,:) = fsolve(@prob2,c0,options,taui);
    c0 = c(i,:);
    i = i + 1;
end

[max_B pos_B] = max(c(:,2));
tau_max_B = tau(pos_B);

fig = figure('Visible', 'off');
plot(tau,c(:,2),'k','LineWidth',1.5);
grid on
hold on 
plot(tau_max_B, max_B, 'ro', 'MarkerSize',8, 'MarkerFaceColor','r');
plot([tau_max_B tau_max_B], [min(c(:,2)) max_B], '--r', 'LineWidth',1.5);
xlabel('$\tau$ (h)', 'Interpreter', 'latex');
ylabel('Concentracion (mol/L)');
saveas(fig, "figura1_problema_2.png");
close(fig);
fprintf("La concentracion maxima de B es %.4f mol/L y se alcanza a las %.4f horas.\n\n", max_B, tau_max_B);

function ec = prob2(c,tau)
    k1 = 1;
    k2 = 0.1;
    ca0 = 1;

    ec(1) = (ca0 - c(1))/tau - 2*k1*c(1).^2;
    ec(1) = -c(2)/tau + k1*c(1).^2 -k2*c(2);
    ec(3) = -c(3)/tau + k2*c(2);
end
