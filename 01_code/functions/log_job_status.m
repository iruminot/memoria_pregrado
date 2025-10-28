function [iter_id] = log_job_status(status, results_base_dir, current_iter_id)
% LOG_JOB_STATUS Escribe el estado de ejecución en archivos planos.
%
%   status: 'start' o 'end'
%   results_base_dir: Ruta a la carpeta principal de resultados (ej. '04_Results')
%   current_iter_id (Opcional): El ID de iteración a usar. Solo necesario para 'end'.
%
%   OUTPUT: iter_id: El número de iteración determinado automáticamente (solo para 'start').

    if nargin < 2
        error('Se requieren al menos dos argumentos: status y results_base_dir.');
    end

    % --- 1. Determinar el ID de Iteración ---
    if strcmpi(status, 'start')
        % ESCANEA ARCHIVOS PREVIOS PARA ENCONTRAR EL PRÓXIMO ID
        
        % Patrón para buscar archivos DONE anteriores (que terminan en .txt)
        done_pattern = fullfile(results_base_dir, '*_status_done.txt');
        archivos_done = dir(done_pattern);
        
        if isempty(archivos_done)
            % Si no hay archivos DONE, es la primera iteración
            iter_id = 1;
        else
            % Extraer los ID de iteración de los nombres de archivo
            % El ID es el número antes del primer '_'
            id_list = zeros(1, length(archivos_done));
            for k = 1:length(archivos_done)
                name = archivos_done(k).name;
                % Buscar el primer '_'
                idx = find(name == '_', 1, 'first');
                if ~isempty(idx)
                    try
                        % Convertir la parte inicial del nombre a número
                        id_list(k) = str2double(name(1:idx-1));
                    catch
                        % Ignorar archivos con nombres no estándar
                        id_list(k) = 0; 
                    end
                end
            end
            
            % El nuevo ID es el máximo encontrado más uno
            iter_id = max(id_list) + 1;
        end

    elseif strcmpi(status, 'end')
        % Para 'end', el ID debe ser proporcionado por el script principal (current_iter_id)
        if nargin < 3
            error('Para el estado ''end'', debe proporcionar el ID de iteración (current_iter_id).');
        end
        iter_id = current_iter_id;
    else
        warning('LOGGING:STATUS_INVALIDO', 'Estado de logging no reconocido: %s', status);
        iter_id = -1; % Devuelve un valor inválido
        return; 
    end


    % --- 2. Preparar la escritura del LOG ---

    % Definir nombres de archivos basados en el ID determinado
    current_time_str = datestr(now, 'yyyymmdd_HHMMSS');
    
    if strcmpi(status, 'start')
        log_filename = sprintf('%d_%s_status_running.txt', iter_id, current_time_str);
    else
        % Usamos el iter_id y una marca de tiempo NUEVA para el archivo DONE
        log_filename = sprintf('%d_%s_status_done.txt', iter_id, current_time_str);
    end

    log_file = fullfile(results_base_dir, log_filename);
    
    if ~exist(results_base_dir, 'dir')
        mkdir(results_base_dir);
    end

    % --- 3. Lógica de Escritura ---
    fid = fopen(log_file, 'w');
    if fid == -1
        error('No se pudo crear el archivo log: %s', log_file);
    end

    if strcmpi(status, 'start')
        fprintf(fid, 'TRABAJO INICIADO: %s\n', current_time_str);
        fprintf(fid, 'Iteración ID: %d\n', iter_id);
        disp(sprintf('\n[LOG] Iniciado (Iteración %d): %s\n', iter_id, current_time_str));
        
    elseif strcmpi(status, 'end')
        fprintf(fid, 'TRABAJO FINALIZADO CON ÉXITO: %s\n', current_time_str);
        
        % Eliminar el archivo RUNNING previo de esa misma iteración
        running_pattern_to_delete = fullfile(results_base_dir, sprintf('%d_*_status_running.txt', iter_id));
        running_files_to_delete = dir(running_pattern_to_delete);
        
        if ~isempty(running_files_to_delete)
            % Se asume que solo debería haber uno (el último iniciado)
            delete(fullfile(running_files_to_delete(1).folder, running_files_to_delete(1).name));
        end
        
        disp(sprintf('\n[LOG] Finalizado (Iteración %d): %s\n', iter_id, current_time_str));
    end
    
    fclose(fid);
end