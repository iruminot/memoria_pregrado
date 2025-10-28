clear all; clc; close all;

currentScriptPath = fileparts(mfilename('fullpath'));
functionsPath = fullfile(currentScriptPath, '..' ,'01_code', 'functions'); addpath(functionsPath);
dataPath = fullfile(currentScriptPath, '..' ,'03_data'); addpath(dataPath);
results_folder = fullfile(currentScriptPath, '..' ,'04_results');

iter_id = log_job_status('start', results_folder);

try

    D = 6/60*ones(1,3);
    sol = dde23(@modelo,D,@histoiainicial,[0 0 5]);

    fig = figure('Visible', 'off');
    subplot(2,1,1);
    plot(sol.x, sol.y(1,:),'k','LineWidth',1.5);
    hold on
    plot(sol.x, sol.y(2,:),'k','LineWidth',1.5);
    plot(sol.x, sol.y(3,:),'k--','LineWidth',1.5);
    hold off
    grid on
    title('Reactor 1');
    ylabel('Concentracion (mol/L)');
    legend('C_A','C_B','C_C');
    subplot(2,1,2);
    plot(sol.x, sol.y(4,:),'k','LineWidth',1.5);
    hold on
    plot(sol.x, sol.y(5,:),'k','LineWidth',1.5);
    plot(sol.x, sol.y(6,:),'k--','LineWidth',1.5);
    hold off
    grid on
    line([0.1 0.1],[0 1.5],'Color','k','LineStyle','--','LineWidth',1.5);
    title('Reactor 2');
    xlabel('Tiempo (h)');
    ylabel('Concentracion (mol/L)');

    print(fig,fullfile(results_folder, sprintf('%d_problema_3.pdf', iter_id)),'-dpdf');
    close(fig);

    %save(fullfile(results_folder, 'resultados_sim.mat'), 'sol');

    log_job_status('end', results_folder, iter_id);

catch ME
    error_time_str = datestr(now, 'yyyymmdd_HHMMSS');
    error_file = fullfile(results_folder, sprintf('%d_%s_status_error.txt', iter_id, error_time_str));
    
    fid = fopen(error_file, 'w');
    fprintf(fid, '--- EJECUCIÓN FALLIDA ---\n');
    fprintf(fid, 'Mensaje: %s\n', ME.message);
    fclose(fid);
    
    running_pattern_to_delete = fullfile(results_folder, sprintf('%d_*_status_running.txt', iter_id));
    running_files_to_delete = dir(running_pattern_to_delete);
    if ~isempty(running_files_to_delete)
        delete(fullfile(running_files_to_delete(1).folder, running_files_to_delete(1).name));
    end
    
    rethrow(ME);
end

function s = histoiainicial(t)
    if t < 0
        s = zeros(6,1);
    else
        s = [ 8.5 0 0 1 0 0];
    end
end

function dc = modelo(t,c,Z)
    V1 = 200; % L
    V2 = 300; % L
    F0 = 300; % L/h
    k11 = 4 ; k12 = 4; % 1/h
    k21 = 3.2; k22 = 3.2; % 1/h
    CAf = 0.5;

    CA1lag = Z(:,1);;
    CB1lag = Z(:,2);;
    CC1lag = Z(:,3);;

    dc = zeros(6,1);
    dc(1) = (F0/V1)*(CAf - c(1)) - (k11 + k21)*c(1);
    dc(2) = - (F0/V1)*c(2) + k11*c(1);
    dc(3) = - (F0/V1)*c(3) + k21*c(1);
    dc(4) = (F0/V2)*(CA1lag(1) - c(4)) - (k12 + k22)*c(4);
    dc(5) = (F0/V2)*(CB1lag(2) - c(5)) + k12*c(4);
    dc(6) = (F0/V2)*(CC1lag(3) - c(6)) + k22*c(4);
end