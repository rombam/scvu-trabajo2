% ----------------------------------------------------------------------- %
%                            SCVu 2021-2022                               %
%                          TRABAJO 2. GRUPO 04                            %
% ----------------------------------------------------------------------- %  
% Diseño de un sistema de aumento de estabilidad y un autopiloto para el
% RPAS Global Hawk. Canal longitudinal.

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
plane_OL = Plane(GlobalHawk);

%% -------------------- CANAL LONGITUDINAL ---------------------
% -- 1. FTs adimensionales: contenidas en plane.lon.G --
% Añadir un factor de escala a cada FT para tener en cuenta las unidades.

% - Dimensiones -
% Si dim = false, todo adimensional y ángulos en radianes
% Si dim = true, velocidad en m/s y ángulos en grados
% Se ajustan los factores de escala (scales) y los tí­tulos de los diagramas
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
	scales = [plane_OL.FC.us/(180/pi);...  % Velocidad/ángulo
              1;...                     % Ángulo/ángulo
              1;...                     % Ángulo/Ángulo
              1/plane_OL.lon.t_lon];       % Vel.angular/Ángulo (tiempo)

    for i=1:length(scales)
        plane_OL.lon.G.(genvarname(Gnames(i))).K_dim = plane_OL.lon.G.(genvarname(Gnames(i))).K * scales(i);
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

%% -- 2. Aplicación de ley de control directo

% Esquema: delta --> K_DL --> G_actuador --> G_longitudinal --> variable
% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulación - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulación - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [º]
delay = 0.1;                  % Delay del actuador
G_a = Utils.padeTF(delay);

% Se seleccionará una ganancia de Direct-Link que proporcione un incremento
% de ángulo de asiento de 1º por cada grado de deflexión del mando de cabeceo

K_dl = 1/plane_OL.lon.G.Gthetadeltae.K; % Calculamos el K_DL para cada variable
    
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
    [y, t_out, ~] = lsim(plane_OL.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
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
    
    % Zoom en perturbación

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
    [y, t_out, ~] = lsim(plane_OL.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
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
    plane_OL.lon.G.(genvarname(Gnames(i))).ramp = lsiminfo(y, t_out);

    if save == false
        sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_e=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/Lon/Ramp_", varnames(i)),'eps')
    end

end
% Rise Time and Time Delay (Hablar con roman)

%% -- 3. Análisis de sensibilidad de la planta libre --
% Análisis del lugar de las raíces del modo de corto periodo
% Tarda mucho, ejecutar solo si es necesario.

plane_Aux = Plane(GlobalHawk);

% Factores de multiplicidad Cmalpha y Cmq
F_alpha  = [-0.5 0 0.5 1 2 3];
F_q      = [-1 0 1 2 4 5];
poles_SP = zeros(length(F_alpha),length(F_q));

% Root locus
figure
polesAUX = pole(plane_OL.lon.G.Galphadeltae.Gfact);
poles_SP_0 = polesAUX(1);
plot([poles_SP_0 conj(poles_SP_0)],'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
%set(gca,'linewidth',1)
xlabel('Re(s)','Interpreter','latex')
ylabel('Im(s)','Interpreter','latex')
xlim([-10 5])
ylim([-17 17])
ax=gca;
ax.FontSize = 20;
ax.YTickLabel.fontsize = 16;
ax.XTickLabel.fontsize = 16;
hold on
grid minor

for i=1:length(F_alpha)
    plane_Aux.model.Cm.alpha = plane_OL.model.Cm.alpha*F_alpha(i);
    for j=1:length(F_q)
        plane_Aux.model.Cm.q = plane_OL.model.Cm.q*F_q(j);
        plane_Aux = plane_Aux.recalc;
        polesAUX = pole(plane_Aux.lon.G.Galphadeltae.Gfact);
        poles_SP(i,j) = polesAUX(1);
        if F_alpha(i)==1 && F_q(j)==1
        else
            if imag(poles_SP(i,j)) ~= 0
                plot([poles_SP(i,j) conj(poles_SP(i,j))],'s','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w')
                hold on
            else
                plot([poles_SP(i,j) 2*real(poles_SP_0)-poles_SP(i,j)],[0 0],'s','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w')    
                hold on
            end
        end
    end
end
sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20])
quiver(-1.95, -4, 0, 3.5, 'k', 'AutoScale', 0, 'MaxHeadSize',0.7)
quiver(-1.95, 4, 0, -3.5, 'k', 'AutoScale', 0, 'MaxHeadSize',0.7)
quiver(-3.5, 0, -1.25, 0, 'k', 'AutoScale', 0, 'MaxHeadSize',10)
quiver(-0.5, 0, 1.4, 0, 'k', 'AutoScale', 0, 'MaxHeadSize',10)

