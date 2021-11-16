% ----------------------------------------------------------------------- %
%                            SCVu 2021-2022                               %
%                          TRABAJO 1. GRUPO 04                            %
% ----------------------------------------------------------------------- % 
% Obtener funciones de transferencia longitudinales/laterales, diagramas de
% Bode y Nichols y respuestas temporales frente a escalón y rampa 
clear, clc, close all

% ------ DESCRIPCIÓN DEL CÓDIGO ------

% - GlobalHawk.m/Learjet24.m : contienen las derivadas de estabilidad de
% dichos aviones. Global Hawk es el avión del Trabajo 1 y Learjet 24 se usó
% para comparar con el libro Leyes de Control de Vuelo (FTs/respuestas
% adimensionales).

% - PLANE.M: crea una clase a partir de las subclases anteriores para
% modelos de aviones concretos que proyecta las derivadas en ejes adecuados 
% y llama al resto de funciones para operar sobre ellas (obtener funciones 
% de transferencia, factorizarlas, etc.). Es la parte más importante.
% Para chequear funciones de transferencia: plane.lon/lat.G.Gsalidaentrada,
% dentro están todos los parámetros de cada función de transferencia.
% NOTA: por defecto las K están en forma adimensional. Para obtenerlas en
% forma dimensional como K_dim, correr los puntos 1 de este script para
% longitudinal y lateral-direccional.
% NOTA 2: para obtener los parámetros de respuesta (overshoot, máximos,
% mínimos, settling time...) hay que ejecutar los apartados de respuesta a
% escalón y rampa (5 y 6 de longitudinal y lateral-direccional).

% - Dynamics.m: contiene los sistemas de ecuaciones dinámicos para cada
% canal por separado. Resuelve las ecuaciones por el método de Cramer y
% devuelve las funciones de transferencia separadas por canal, entrada y
% salida, factorizadas y con los datos de polos, ceros y ganancias.

% - Utils.m: recoge funciones de utilidad variada (conversión de unidades,
% ISA, factorización de funciones de transferencia, etc.).

%% ------ Creación del modelo dinámico y obtención de FTs ------
plane = Plane(GlobalHawk);

%% -------------------- CANAL LONGITUDINAL ---------------------
% -- 1. FTs adimensionales: contenidas en plane.lon.G --
% Añadir un factor de escala a cada FT para tener en cuenta las unidades.

% - Dimensiones -
% Si dim = false, todo adimensional y ángulos en radianes
% Si dim = true, velocidad en m/s y ángulos en grados
% Se ajustan los factores de escala (scales) y los títulos de los diagramas
% Para el trabajo, debe estar por defecto en true. Para comparar con el
% libro en false
% Si save = true, guarda gráficas en formato .eps en ./Graficas/Lon/

dim = true;
save = false;

varnames = ["udeltae"; "alphadeltae"; "thetadeltae"; "qdeltae"];
Gnames = strcat('G', varnames);

if dim == false
    scales = [1; 1; 1; 1];
    % - Parámetros diagramas de Bode/Nichols - 
    Magnames_latex  = ["$\left | G_{\hat{u} \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\alpha \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\theta \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\hat{q} \delta_{e}}(i\omega) \right |$"];
    Gnames_latex    = ["$G_{\hat{u} \delta_{e}}(i\omega)$";...
                       "$G_{\alpha \delta_{e}}(i\omega)$";...
                       "$G_{\theta \delta_{e}}(i\omega)$";...
                       "$G_{\hat{q} \delta_{e}}(i\omega)$"];
                
    % - Parámetros gráficas de respuesta -
    Respnames_latex = ["$\Delta \hat{u}$ [-]";...
                       "$\Delta \alpha$ [rad]";...
                       "$\Delta \theta$ [rad]";...
                       "$\Delta \hat{q}$ [-]"];
    angleunits      = '\hspace{0.75 mm}rad$';
else
	scales = [plane.FC.us/(180/pi);...  % Velocidad/ángulo
              1;...                     % Ángulo/ángulo
              1;...                     % Ángulo/ángulo
              1/plane.lon.t_lon];       % Vel.angular/ángulo (tiempo)

    for i=1:length(scales)
        plane.lon.G.(genvarname(Gnames(i))).K_dim = plane.lon.G.(genvarname(Gnames(i))).K * scales(i);
    end

    % - Parámetros diagramas de Bode/Nichols - 
    Magnames_latex  = ["$\left | G_{u \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\alpha \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\theta \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{q \delta_{e}}(i\omega) \right |$"];
    Gnames_latex    = ["$G_{u \delta_{e}}(i\omega)$";...
                       "$G_{\alpha \delta_{e}}(i\omega)$";...
                       "$G_{\theta \delta_{e}}(i\omega)$";...
                       "$G_{q \delta_{e}}(i\omega)$"];
                
    % - Parámetros gráficas de respuesta -
    Respnames_latex = ["$\Delta u$ [m/s]";...
                       "$\Delta \alpha$ [$^\circ$]";...
                       "$\Delta \theta$ [$^\circ$]";...
                       "$\Delta q$ [$^\circ$/s]"];
    angleunits      = '^{\circ}$';
end
    
%% -- 2. Diagramas de Bode --
% Las unidades vienen definidas en el apartado anterior (con dims)
lims = [1e-2 1e2];    % Límites de frecuencia del diagrama

bopts = bodeoptions;
	bopts.PhaseMatching = 'on';
	bopts.PhaseMatchingFreq = lims(1);
	bopts.PhaseMatchingValue = 0;
    
	bopts.XLabel.String = '$\omega$';
	bopts.XLabel.FontSize = 13;
	bopts.XLabel.Interpreter = 'latex';
    
	bopts.YLabel.FontSize = 13;
	bopts.YLabel.Interpreter = 'latex';
    
	bopts.Title.FontSize = 14;
	bopts.Title.Interpreter = 'latex';
    
for i = 1:length(varnames)
	figure
	bopts.YLabel.String = [Magnames_latex(i) "$\phi(\omega)$"];
    if save == false
	    bopts.Title.String = strcat("Bode plot ", Gnames_latex(i));
    else
        bopts.Title.String = "";
    end
	bodeplot(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), {lims(1), lims(2)}, bopts, 'k');
	grid on
    if save == true
        saveas(gcf,strcat("Graficas/Lon/Bode_G", varnames(i)),'eps')
    end
end
    
%% -- 3. Módulo de FT a las frecuencias de los modos propios --
% Se asumen con dimensiones (m/s vs º, º vs º, rad/s vs º)
% Frecuencias de los modos propios: corto período (SP) y fugoide (P)
wn = 1i*plane.lon.G.Gudeltae.poles.wn; % Contiene SP y P

modules.lon(:,1) = ["" strcat('|',Gnames,'|').'];
modules.lon(1,1:3) = ["" "w_SP" "w_P"];

for i = 1:length(varnames)
    for k = 1:length(wn)
        modules.lon(i+1,k+1) = 20*log10(abs(evalfr(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), wn(k))));
    end
end

%% -- 4. Diagramas de Nichols --
% Al igual que para el diagrama de Bode, se utilizan dimensiones de m/s y º
lims = [1e-2 1e2];

nopts = nicholsoptions;
	nopts.PhaseMatching = 'on';
	nopts.PhaseMatchingFreq = lims(1);
	nopts.PhaseMatchingValue = 0;
    
	nopts.XLabel.String = '$\phi(\omega)$';
	nopts.XLabel.FontSize = 13;
	nopts.XLabel.Interpreter = 'latex';
    
	nopts.YLabel.FontSize = 13;
	nopts.YLabel.Interpreter = 'latex';
    
	nopts.Title.FontSize = 14;
	nopts.Title.Interpreter = 'latex';
    
for i = 1:length(varnames)
	figure
	nopts.YLabel.String = Magnames_latex(i);
    if save == false
	    nopts.Title.String = strcat("Nichols plot ", Gnames_latex(i));
    else
        nopts.Title.String = "";
    end
	nicholsplot(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), {lims(1), lims(2)}, nopts, 'k');
	grid on
    if save == true
        saveas(gcf,strcat("Graficas/Lon/Nichols_G", varnames(i)),'epsc')
    end
end
    
%% -- 5. Respuesta a escalón en deltae --
% - Definir escalón -
tsimtr = [150 10 100 10];     % Tiempo total de simulación - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulación - estacionario [s]
amp    = 1;                   % Amplitud del escalón [º]
    
for i = 1:length(varnames)

    fig = figure('Position', [100, 100, 1120, 420]);
    
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = amp*ones(1,length(t));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==1
        hold on
        plot(t,u,'k-.')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = amp*ones(1,length(t));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==1
        hold on
        plot(t,u,'k-.')
        legend('Respuesta', '$\Delta \delta_e$', 'Interpreter', 'latex')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane.lon.G.(genvarname(Gnames(i))).step = stepinfo(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i));
    
    if save == false
        sgtitle(strcat('\underline{Respuesta step $\Delta\delta_e=', num2str(amp), angleunits, '}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Lon/Step_", varnames(i)),'eps')
    end
end

%% -- 6. Respuesta a rampa saturada en deltae --
% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulación - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulación - estacionario [s]
tramp  = 10;                  % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [º]
    
for i = 1:length(varnames)
    fig = figure('Position', [100, 100, 1120, 420]);
    
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==1
        hold on
        plot(t,u,'k-.')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==1
        hold on
        plot(t,u,'k-.')
        legend('Respuesta', '$\Delta \delta_e$', 'Interpreter', 'latex')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane.lon.G.(genvarname(Gnames(i))).ramp = lsiminfo(y, t_out);

    if save == false
        sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_e=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Lon/Ramp_", varnames(i)),'eps')
    end

end
%% ----------------------- Sacar Gs en consola ------------------------
% Para comparar con los valores del libro si procede
varnames = ["udeltae"; "alphadeltae"; "thetadeltae"; "qdeltae"];
for i = 1:length(varnames)
    %disp(strcat('------- G', varnames(i), ' -------'))
    fprintf('########## G%s ##########\n\n', varnames(i))
    pretty(vpa(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact_sym"))*scales(i), 4))
    fprintf('\n\n')
end

%% -------------------- CANAL LATERAL-DIRECCIONAL ---------------------
% -- 1. FTs adimensionales: contenidas en plane.lat.G --
% Añadir un factor de escala a cada FT para tener en cuenta las unidades.

% - Dimensiones -
% Si dim = false, todo adimensional y ángulos en radianes
% Si dim = true, velocidad en m/s y ángulos en grados
% Se ajustan los factores de escala (scales) y los títulos de los Bodes
% Para el trabajo debe estar por defecto en true. Para comparar con el
% libro en false.
% Si save = true, guarda gráficas en formato .eps en ./Graficas/Latdir/

dim = true;
save = false;

varnames = ["betadeltaa" "phideltaa" "rdeltaa" "pdeltaa" "psideltaa" "betadeltar" "phideltar" "rdeltar" "pdeltar" "psideltar"];
Gnames = strcat('G', varnames);

if dim == false
    scales = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
    % - Parámetros diagramas de Bode/Nichols - 
    Magnames_latex  = ["$\left | G_{\beta \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\phi \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\hat{r} \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\hat{p} \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\psi \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\beta \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{\phi \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{\hat{r} \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{\hat{p} \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{\psi \delta_{r}}(i\omega) \right |$"];       
    Gnames_latex    = ["$G_{\beta \delta_{a}}(i\omega)$";...
                       "$G_{\phi \delta_{a}}(i\omega)$";...
                       "$G_{\hat{r} \delta_{a}}(i\omega)$";...
                       "$G_{\hat{p} \delta_{a}}(i\omega)$";...
                       "$G_{\psi \delta_{a}}(i\omega)$";...
                       "$G_{\beta \delta_{r}}(i\omega)$";...
                       "$G_{\phi \delta_{r}}(i\omega)$";...
                       "$G_{\hat{r} \delta_{r}}(i\omega)$";...
                       "$G_{\hat{p} \delta_{r}}(i\omega)$";...
                       "$G_{\psi \delta_{r}}(i\omega)$"];
    
                
    % - Parámetros gráficas de respuesta -
    Respnames_latex = ["$\Delta \beta$ [rad]";...
                       "$\Delta \phi$ [rad]";...
                       "$\Delta \hat{r}$ [-]";...
                       "$\Delta \hat{p}$ [-]";...
                       "$\Delta \psi$ [rad]";...
                       "$\Delta \beta$ [rad]";...
                       "$\Delta \phi$ [rad]";...
                       "$\Delta \hat{r}$ [-]";...
                       "$\Delta \hat{p}$ [-]";...
                       "$\Delta \psi$ [rad]"];
    angleunits      = '\hspace{0.75 mm}rad$';
    
else
	scales = [1;...                     % Ángulo/Ángulo
              1;...                     % Ángulo/ángulo
              1/plane.lat.t_lat;...     % Vel.angular/ángulo (tiempo)
              1/plane.lat.t_lat;...     % Vel.angular/ángulo (tiempo)
              1;...                     % Ángulo/Ángulo
              1;...                     % Ángulo/Ángulo
              1;...                     % Ángulo/ángulo
              1/plane.lat.t_lat;...     % Vel.angular/ángulo (tiempo)
              1/plane.lat.t_lat;...     % Vel.angular/ángulo (tiempo)
              1];                       % Ángulo/Ángulo

    for i=1:length(scales)
        plane.lat.G.(genvarname(Gnames(i))).K_dim = plane.lat.G.(genvarname(Gnames(i))).K * scales(i);
    end

    % - Parámetros diagramas de Bode/Nichols - 
    Magnames_latex  = ["$\left | G_{\beta \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\phi \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{r \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{p \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\psi \delta_{a}}(i\omega) \right |$";...
                       "$\left | G_{\beta \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{\phi \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{r \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{p \delta_{r}}(i\omega) \right |$";...
                       "$\left | G_{\psi \delta_{r}}(i\omega) \right |$"];    
    Gnames_latex    = ["$G_{\beta \delta_{a}}(i\omega)$";...
                       "$G_{\phi \delta_{a}}(i\omega)$";...
                       "$G_{r \delta_{a}}(i\omega)$";...
                       "$G_{p \delta_{a}}(i\omega)$";...
                       "$G_{\psi \delta_{a}}(i\omega)$";...
                       "$G_{\beta \delta_{r}}(i\omega)$";...
                       "$G_{\phi \delta_{r}}(i\omega)$";...
                       "$G_{r \delta_{r}}(i\omega)$";...
                       "$G_{p \delta_{r}}(i\omega)$";...
                       "$G_{\psi \delta_{r}}(i\omega)$"];
                
    % - Parámetros gráficas de respuesta -
    Respnames_latex = ["$\Delta \beta$ [$^\circ$]";...
                       "$\Delta \phi$ [$^\circ$]";...
                       "$\Delta r$ [$^\circ/s$]";...
                       "$\Delta p$ [$^\circ/s$]";...
                       "$\Delta \psi$ [$^\circ$]";...
                       "$\Delta \beta$ [$^\circ$]";...
                       "$\Delta \phi$ [$^\circ$]";...
                       "$\Delta r$ [$^\circ/s$]";...
                       "$\Delta p$ [$^\circ/s$]";...
                       "$\Delta \psi$ [$^\circ$]"];
    angleunits      = '^{\circ}$';
end
    
%% -- 2. Diagramas de Bode --
% Las unidades vienen definidas en el apartado anterior (con dims)
lims = [1e-5 1e2];    % Límites de frecuencia del diagrama

bopts = bodeoptions;
	bopts.PhaseMatching = 'on';
	bopts.PhaseMatchingFreq = lims(1);
	bopts.PhaseMatchingValue = 0;
    
	bopts.XLabel.String = '$\omega$';
	bopts.XLabel.FontSize = 13;
	bopts.XLabel.Interpreter = 'latex';
    
	bopts.YLabel.FontSize = 13;
	bopts.YLabel.Interpreter = 'latex';
    
	bopts.Title.FontSize = 14;
	bopts.Title.Interpreter = 'latex';
    
for i = 1:length(varnames)
	figure
	bopts.YLabel.String = [Magnames_latex(i) "$\phi(\omega)$"];
    if save == false
	    bopts.Title.String = strcat("Bode plot ", Gnames_latex(i));
    else
        bopts.Title.String = "";
    end
	bodeplot(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), {lims(1), lims(2)}, bopts, 'k');
    grid on
    if save == true
        saveas(gcf,strcat("Graficas/Latdir/Bode_G", varnames(i)),'eps')
    end
end
    
%% -- 3. Módulo de FT a las frecuencias de los modos propios --
% Se asumen con dimensiones (m/s vs º, º vs º, rad/s vs º)
% Frecuencias de los modos propios: convergencia en balance (R, real), modo
% espiral (S, real) y balanceo holandés (DR, oscil.)
taun = plane.lat.G.Gbetadeltaa.poles.tau; % Contiene R, S
wn   = plane.lat.G.Gbetadeltaa.poles.wn;  % Contiene DR
modes = [1i*wn taun];
modules.lat(:,1) = ["" strcat('|',Gnames,'|')].';
modules.lat(1,1:4) = ["" "w_DR" "1/tauR" "1/tauS"];
for i = 1:length(varnames)
    for k = 1:length(modes)
        modules.lat(i+1,k+1) = 20*log10(abs(evalfr(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), modes(k))));
    end
end

%% -- 4. Diagramas de Nichols --
% Al igual que para el diagrama de Bode, se utilizan dimensiones de m/s y º
lims = [1e-5 1e2];

nopts = nicholsoptions;
	nopts.PhaseMatching = 'on';
	nopts.PhaseMatchingFreq = 0.01;
	nopts.PhaseMatchingValue = 0;
    
	nopts.XLabel.String = '$\omega$';
	nopts.XLabel.FontSize = 13;
	nopts.XLabel.Interpreter = 'latex';
    
	nopts.YLabel.FontSize = 13;
	nopts.YLabel.Interpreter = 'latex';
    
	nopts.Title.FontSize = 14;
	nopts.Title.Interpreter = 'latex';
    
for i = 1:length(varnames)
    figure
	nopts.YLabel.String = Magnames_latex(i);

    if save == false
	    nopts.Title.String = strcat("Nichols plot ", Gnames_latex(i));
    else
        nopts.Title.String = "";
    end

	nicholsplot(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), {lims(1), lims(2)}, nopts, 'k');
	grid on
    if save == true
        saveas(gcf,strcat("Graficas/Latdir/Nichols_G", varnames(i)),'epsc')
    end
end
    
%% -- 5. Respuesta a escalón en deltaa / deltar --
% - Definir escalón -
tsimtr = [150 10 100 10 10 150 10 100 10 10];        % Tiempo total de simulación - transitorio [s]
tsimst = [100000 100000 100000 100000 100000 100000 100000 100000 100000 100000];  % Tiempo total de simulación - estacionario [s]
amp    = 1;                                          % Amplitud del escalón [º]
    
for i = 1:5

    fig = figure('Position', [100, 100, 1120, 420]);
    
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = amp*ones(1,length(t));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==4
        hold on
        plot(t,u,'k-.')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = amp*ones(1,length(t));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==4
        hold on
        plot(t,u,'k-.')
        legend('Respuesta', '$\Delta \delta_a$', 'Interpreter', 'latex', 'Location', 'northeast')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane.lat.G.(genvarname(Gnames(i))).step = stepinfo(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i));
    
    if save == false
        sgtitle(strcat('\underline{Respuesta step $\Delta\delta_a=', num2str(amp), angleunits, '}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Latdir/Step_G", varnames(i)),'eps')
    end
end

for i = 6:10

    fig = figure('Position', [100, 100, 1120, 420]);
    
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = amp*ones(1,length(t));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    if i==6
        hold on
        plot(t,u,'k-.')
        ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(u)+(abs(max(y)-min(y)))*0.1])
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = amp*ones(1,length(t));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==6
        hold on
        plot(t,u,'k-.')
        legend('Respuesta', '$\Delta \delta_r$', 'Interpreter', 'latex', 'Location', 'east')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane.lat.G.(genvarname(Gnames(i))).step = stepinfo(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i));
    if save == false
        sgtitle(strcat('\underline{Respuesta step $\Delta\delta_r=', num2str(amp), angleunits, '}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Latdir/Step_G", varnames(i)),'eps')
    end
end

%% -- 6. Respuesta a rampa saturada en deltaa / deltar --
% Definir rampa
tsimtr = [150 50 30 50 100 50 50 50 100 100];     % Tiempo total de simulación - transitorio [s]
tsimst = [100000 100000 100000 100000 100000 100000 100000 100000 100000 100000];   % Tiempo total de simulación - estacionario [s]
tramp  = 10;                                          % Tiempo de rampeo [s]
amp    = 1;                                           % Amplitud de la rampa [º]


for i = 1:5
    fig = figure('Position', [100, 100, 1120, 420]);
    
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==4
        hold on
        plot(t,u,'k-.')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==4
        hold on
        plot(t,u,'k-.')
        legend('Respuesta', '$\Delta \delta_a$', 'Interpreter', 'latex', 'Location', 'northeast')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane.lat.G.(genvarname(Gnames(i))).ramp = lsiminfo(y, t_out);

    if save == false
        sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_a=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Latdir/Ramp_G", varnames(i)),'eps')
    end

end


for i = 6:10
    fig = figure('Position', [100, 100, 1120, 420]);
    
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    if i==6
        hold on
        plot(t,u,'k-.')
        ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(u)+(abs(max(y)-min(y)))*0.1])
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    if i==6
        hold on
        plot(t,u,'k-.')
        legend('Respuesta', '$\Delta \delta_r$', 'Interpreter', 'latex', 'Location', 'east')
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane.lat.G.(genvarname(Gnames(i))).ramp = lsiminfo(y, t_out);    
    if save == false
        sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_r=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Latdir/Ramp_G", varnames(i)),'eps')
    end
end

%% ----------------------- Sacar Gs en consola ------------------------
% Para comparar con los valores del libro si procede
varnames = ["betadeltaa" "phideltaa" "rdeltaa" "pdeltaa" "psideltaa" "betadeltar" "phideltar" "rdeltar" "pdeltar" "psideltar"];

for i = 1:length(varnames)
    fprintf('########## G%s ##########\n\n', varnames(i))
    pretty(vpa(plane.lat.G.(genvarname(Gnames(i))).(genvarname("Gfact_sym"))*scales(i), 4))
    fprintf('\n\n')
 end