% ----------------------------------------------------------------------- %
%                            SCVu 2021-2022                               %
%                          TRABAJO 2. GRUPO 04                            %
% ----------------------------------------------------------------------- %  
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
              1;...                     % �?ngulo/ángulo
              1;...                     % �?ngulo/ángulo
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

%% -- 2. Aplicaci�n de ley de control directo

% Esquema: delta --> K_DL --> G_actuador --> G_longitudinal --> variable
% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulación - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulación - estacionario [s]
tramp  = 5;                  % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [º]
delay = 0.1;                  % Delay del actuador
G_a = Utils.padeTF(delay);

% Se seleccionar� una ganancia de Direct-Link que proporcione un incremento
% de �ngulo de asiento de 1� por cada grado de deflexi�n del mando de cabeceo

K_dl = 1/plane.lon.G.Gthetadeltae.K; %Calculamos el K_DL para cada variable
    
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
    [y, t_out, ~] = lsim(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
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
    
    % Zoom en perturbaci�n

    axes('Position',[.30 .67 .15 .15])
    box on; % put box around new pair of axes
    indexOfInterest = (t < 4) & (t > 0); % range of t near perturbation
    plot(t(indexOfInterest),y(indexOfInterest),'k'); hold on % plot on new axes
    grid on;
    axis tight
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(plane.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
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
% Rise Time and Time Delay (Hablar con roman)

%% --


