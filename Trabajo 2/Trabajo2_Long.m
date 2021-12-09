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
elseif dim == true
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

K_dl = 1/plane_OL.lon.G.Gthetadeltae.K; % Calculamos el K_DL para theta


fprintf("--------- Tiempos Característicos Lazo Abierto con K_dl + Gact + Gp -----------\n");

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
    
    % Rise Time and Time Delay Calculation and printing in console
    
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf(varnames(i) + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");

end

%% -- 3. Análisis de sensibilidad de la planta libre --
% Análisis del lugar de las raíces del modo de corto periodo
% Tarda mucho, ejecutar solo si es necesario.
clear
plane_OL  = Plane(GlobalHawk);
plane_Aux = Plane(GlobalHawk);

% Factores de multiplicidad Cmalpha y Cmq
F_alpha  = [-0.5 0 0.5 1 2 3 4];
F_q      = [-0.5 0 0.5 1 2 3 4];
poles_SP = zeros(length(F_alpha),length(F_q));

% Root locus
figure
polesAUX = pole(plane_OL.lon.G.Galphadeltae.Gfact);
poles_SP_0 = polesAUX(1);
plot([poles_SP_0 conj(poles_SP_0)],'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
xlabel('Re(s)','Interpreter','latex')
ylabel('Im(s)','Interpreter','latex')
xlim([-10 5])
ylim([-18 18])
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

%% -- 4. Barrido de ganancias de realimentación --
% WORK IN PROGRESS
% Análisis del lugar de las raíces del modo de corto periodo con
% realimentación en alpha y q.
% Análisis de la función de transferencia de deltae a theta en lazo
% cerrado.
clear
plane_OL = Plane(GlobalHawk);

% Funciones de transferencia y Kdl
Ga_deltae = Utils.forlagTF;
Gs_alpha = 1;
Gs_q = 1;
Gthetadeltae_pl = plane_OL.lon.G.Gthetadeltae.Gfact;
Galphadeltae_pl = plane_OL.lon.G.Galphadeltae.Gfact;
Gqdeltae_pl     = plane_OL.lon.G.Gqdeltae.Gfact;
K_DL = 1/plane_OL.lon.G.Gthetadeltae.K;

% Factores de multiplicidad Cmalpha y Cmq
F_alpha  = [-0.5 -0.2 0 0.5 1 2 3 4];
%F_alpha = 1;
F_q      = [-0.5 -0.2 0 0.5 1 2 3 4];
%F_q = 1;
Kalpha   = - (F_alpha-1)*plane_OL.lon.Cm.alpha/plane_OL.lon.Cm.deltae;
Kq       = - (F_q-1)*plane_OL.lon.Cm.q/plane_OL.lon.Cm.deltae*plane_OL.lon.t_lon;
poles_SP = zeros(length(Kalpha),length(Kq));

% Root locus. Ploteamos primero el corto periodo y el polo del actuador
% base
figure
polesAUX = pole(plane_OL.lon.G.Galphadeltae.Gfact);
poles_SP_0 = polesAUX(1);
plot([poles_SP_0 conj(poles_SP_0)],'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
hold on
plot(-10,0,'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
xlabel('Re(s)','Interpreter','latex')
ylabel('Im(s)','Interpreter','latex')
xlim([-20 3])
ylim([-14 14]) 
hold on
grid minor

for i=1:length(Kalpha)
    for j=1:length(Kq)
        % -- Ensamblaje de la FT en lazo cerrado --
        s = tf('s');
        Gthetadeltae_CL{i,j} = K_DL*(Ga_deltae*Gthetadeltae_pl)/(1+Ga_deltae*Kalpha(i)*Gs_alpha*Galphadeltae_pl +...
                                     Kq(j)*Gs_q*s*Gthetadeltae_pl);
        polesAUX = round(pole(Gthetadeltae_CL{i,j}),4);
        zerosAUX = round(zero(Gthetadeltae_CL{i,j}),4);
        % Quita a la lista de polos los ceros que sean iguales a estos.
        % Para eliminar los polos residuales de la planta libre.
        polesAUXclean = setdiff(polesAUX,zerosAUX);
        if Kalpha(i)==0 && Kq(j)==0
        else
            if Kalpha(i) < 0
                marker = 's';
            else
                marker = '^';
            end
            color = 'w';
            if Kq(j) < 0
                color = 'w';
            else
                color = [197/255 197/255 197/255];
            end
            for k=1:length(polesAUXclean)
                % -- Bucle en todos los polos --
                % Filtramos los polos cerca del origen (modo fugoide) a
                % mano. Mejorable.
                if abs(imag(polesAUXclean(k))) >= 0.3
                    plot(polesAUXclean(k),marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
                    hold on
                elseif abs(real(polesAUXclean(k))) >= 0.15
                    plot(polesAUXclean(k),0,marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
                    hold on
                elseif abs(real(polesAUXclean(k))) >= 10
                else
                end
            end
        end
    end
end
sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20 24])

%% -- 5. Diseño del SAS --
plane_cont = Plane(GlobalHawk);

% Ganancias de realimentación
Kalpha = 2;
Kq = -0.2;
% plane_cont.lon.Cont.Gudeltae = K_DL*(Ga_deltae*plane_cont.lon.G.Gudeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
%                                      Kq*Gs_q*s*plane_OL.lon.t_lon*plane_cont.lon.G.Gthetadeltae.Gfact);
% plane_cont.lon.Cont.Galphadeltae = K_DL*(Ga_deltae*plane_cont.lon.G.Galphadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
%                                      Kq*Gs_q*s*plane_OL.lon.t_lon*plane_cont.lon.G.Gthetadeltae.Gfact);
% plane_cont.lon.Cont.Gthetadeltae = K_DL*(Ga_deltae*plane_cont.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
%                                      Kq*Gs_q*s*plane_OL.lon.t_lon*plane_cont.lon.G.Gthetadeltae.Gfact);
% plane_cont.lon.Cont.Gqdeltae = K_DL*(Ga_deltae*plane_cont.lon.G.Gqdeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
%                                      Kq*Gs_q*s*plane_OL.lon.t_lon*plane_cont.lon.G.Gthetadeltae.Gfact);

plane_cont.lon.Cont.Gudeltae = (Ga_deltae*plane_cont.lon.G.Gudeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_cont.lon.G.Gqdeltae.Gfact);
plane_cont.lon.Cont.Galphadeltae = (Ga_deltae*plane_cont.lon.G.Galphadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_cont.lon.G.Gqdeltae.Gfact);
plane_cont.lon.Cont.Gthetadeltae = (Ga_deltae*plane_cont.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_cont.lon.G.Gqdeltae.Gfact);
plane_cont.lon.Cont.Gqdeltae = (Ga_deltae*plane_cont.lon.G.Gqdeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_cont.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_cont.lon.G.Gqdeltae.Gfact);

GSAS_theta = plane_cont.lon.Cont.Gthetadeltae;

                                 
%% -- 6.1. Diseño AP: Root Locus (MIRAR ESTO QUE CREO QUE ESTA MAL)

% Valores PID
k_p = [-0.05,-0.1,-0.15,-0.2,-0.5,-1,-2,-5];
k_i = [0,-0.01,-0.02,-0.05,-0.1,-0.2,-0.5,-1];
k_d = 0; % Lo ponemos a cero porque basicamente no afecta demasiado a la estabilidad del sistema

% base
figure
polesAUX = pole(plane_OL.lon.G.Galphadeltae.Gfact);
poles_SP_0 = polesAUX(1);
plot([poles_SP_0 conj(poles_SP_0)],'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
hold on
plot(-10,0,'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
xlabel('Re(s)','Interpreter','latex')
ylabel('Im(s)','Interpreter','latex')
xlim([-15 3])
ylim([-20 20])
hold on
grid minor

Gs_theta = 1;

for i=1:length(k_p)
    for j=1:length(k_i)             
        % -- Ensamblaje de la FT en lazo cerrado --
        PID_AP = pid(k_p(i),k_i(j),k_d);
        Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
        polesAUX = round(pole(Gap_theta),4);
        zerosAUX = round(zero(Gap_theta),4);
        % Quita a la lista de polos los ceros que sean iguales a estos.
        % Para eliminar los polos residuales de la planta libre.
        polesAUXclean = setdiff(polesAUX,zerosAUX);
        if k_p(i)==0 && k_i(j)==0
        else
            if k_p(i) < 0
                marker = 's';
            else
                marker = '^';
            end
            for k=1:length(polesAUXclean)
                % -- Bucle en todos los polos --
                % Filtramos los polos cerca del origen (modo fugoide) a
                % mano. Mejorable.
                if abs(imag(polesAUXclean(k))) >= 0.3
                    plot(polesAUXclean(k),marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w')
                    hold on
                elseif abs(real(polesAUXclean(k))) >= 0.15
                    plot(polesAUXclean(k),0,marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w')    
                    hold on
                else
                end
            end
        end
    end
end
sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20])




%% -- 6.2 Diseño de AP: PID

% Datos preliminares para respuesta en rampa
tsimtr = [150 50 100 100];    % Tiempo total de simulación - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulación - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [º]

syms = ["-","--",":","-.","-","--"];

% Barrido PID

% 1) Proporcional

lims = [1e-2 1e2];

% Valores PID
k_p = [-0.05,-0.1,-0.15,-0.2,-0.5,-1];
k_i = 0;
k_d = 0;
Gs_theta = Utils.padeTF();

% Margenes y fases

Gms_p = zeros(1,length(k_p));
Pms_p = zeros(1,length(k_p));
wp_p = zeros(1,length(k_p));

f_nichols = figure;
f_step = figure('Position', [100, 100, 600, 420]);

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
    
    nopts.YLabel.String = Magnames_latex(3);
    nopts.Title.String = strcat("Nichols plot ", Gnames_latex(3));

set(0, 'CurrentFigure', f_step)
    ax = gca; hold(ax,'on');
    ax.FontSize = 13; 
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    %ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])

for i=1:length(k_p)
    
    PID_AP = pid(k_p(i),k_i,k_d);
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    Gap_theta_ol = PID_AP*(GSAS_theta);
    
    % Nichols plotting
    set(0, 'CurrentFigure', f_nichols)
    nicholsplot(Gap_theta_ol, {lims(1), lims(2)},nopts); hold on
    
    % Step Plotting
    set(0, 'CurrentFigure', f_step)
    % Transitorio
    t = 0:0.05:tsimtr(3);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out, ~] = lsim(Gap_theta, u, t);
    hold(ax,'on');
    plot(t_out,y,'LineStyle',syms(i),'Color','k'); hold on
    if(i == 1)
        text(t(end)*0.75,y(end)*1.1,'$[k_{\delta_{e} \theta}]_P$' + " = " + num2str(k_p(i)),'FontSize',10,'Interpreter','latex');hold on
    else
        text(t(end)*0.75,y(end)*1.1,num2str(k_p(i)),'FontSize',10,'Interpreter','latex');hold on
    end
    ax = gca;
    ax.FontSize = 13; 
    title('Respuesta Rampa Saturada','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    

    % Sacamos los valores de margen de fase y ganancia
    [Gm,Pm,Wcg,Wcp] = margin(Gap_theta_ol);
    Gms_p(i) = Gm; Pms_p(i) = Pm; wp_p(i) = Wcp; 
    
    % Calculo de Tiempos característicos
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf("Kp = "+ k_p(i) + ", Kd = " + k_d + ", Ki = " + k_i + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");
    
end

set(0, 'CurrentFigure', f_step)
    grid minor
    box on
set(0, 'CurrentFigure', f_nichols)
    grid minor
    box on
    

    
% 2) Integrador

lims = [1e-2 1e2];

% Valores PID
k_p = -0.2;
k_i = [0,-0.01,-0.02,-0.05,-0.2];
k_d = 0;
Gs_theta = Utils.padeTF();

% Margenes y fases

Gms_i = zeros(1,length(k_i));
Pms_i = zeros(1,length(k_i));
wp_i = zeros(1,length(k_i));


f_nichols = figure;
f_step = figure('Position', [100, 100, 600, 420]);

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
    
    nopts.YLabel.String = Magnames_latex(3);
    nopts.Title.String = strcat("Nichols plot ", Gnames_latex(3));

set(0, 'CurrentFigure', f_step)
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    %ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
%     axes('Position',[.67 .73 .15 .15]); hold on
%     box on; % put box around new pair of axes

for i=1:length(k_i)
    
    PID_AP = pid(k_p,k_i(i),k_d);
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    Gap_theta_ol = PID_AP*(GSAS_theta);
    
    % Nichols plotting
    set(0, 'CurrentFigure', f_nichols)
    nicholsplot(Gap_theta_ol, {lims(1), lims(2)}, nopts, syms(i)); hold on
    
    % Step Plotting
    set(0, 'CurrentFigure', f_step)
    % Transitorio
    t = 0:0.05:tsimtr(4);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out, ~] = lsim(Gap_theta, u, t);
    plot(t_out, y, syms(i),'Color','k'); hold on
    ax = gca;
    ax.FontSize = 13; 
    title('Respuesta Rampa Saturada','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
        % Zoom en perturbación
%     ax2.indexOfInterest = (t < 5) & (t > 0); % range of t near perturbation
%     ax2.plot(t(indexOfInterest),y(indexOfInterest),'k'); hold on % plot on new axes
%     grid on;
%     axis tight
    
    % Calculo de margenes de ganancia y fase
    [Gm,Pm,Wcg,Wcp] = margin(Gap_theta_ol);
    Gms_i(i) = Gm; Pms_i(i) = Pm; wp_i(i) = Wcp; 
    
    % Calculo de Tiempos característicos
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf("Kp = "+ k_p + ", Kd = " + k_d + ", Ki = " + k_i(i) + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");
end

set(0, 'CurrentFigure', f_step)
    grid minor
    hold on
    box on
set(0, 'CurrentFigure', f_nichols)
    grid minor
    hold on
    box on
    

% 3) Derivativo

lims = [1e-2 1e2];

% Valores PID
k_p = -1;
k_i = -0.2;
k_d = [0,-0.05,-0.2];
Gs_theta = Utils.padeTF();

% Margenes y fases

Gms_d = zeros(1,length(k_i));
Pms_d = zeros(1,length(k_i));
wp_d = zeros(1,length(k_i));


f_nichols = figure;
f_step = figure('Position', [100, 100, 600, 420]);
f_step_zoom = figure('Position',[100, 100, 600, 420]);

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
    
    nopts.YLabel.String = Magnames_latex(3);
    nopts.Title.String = strcat("Nichols plot ", Gnames_latex(3));

    
set(0, 'CurrentFigure', f_step_zoom)
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')    

set(0, 'CurrentFigure', f_step)
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    %ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
%     axes('Position',[.67 .73 .15 .15]); hold on
%     box on; % put box around new pair of axes

for i=1:length(k_d)
    
    PID_AP = pid(k_p,k_i,k_d(i));
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    Gap_theta_ol = PID_AP*(GSAS_theta);
    
    % Nichols plotting
    set(0, 'CurrentFigure', f_nichols)
    nicholsplot(Gap_theta_ol, {lims(1), lims(2)}, nopts, syms(i)); hold on
    
    % Step Plotting
    set(0, 'CurrentFigure', f_step)
    % Transitorio
    t = 0:0.05:tsimtr(4);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out, ~] = lsim(Gap_theta, u, t);
    plot(t_out, y,syms(i),'Color','k'); hold on
    ax = gca;
    ax.FontSize = 13; 
    title('Respuesta Rampa Saturada','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    set(0,'CurrentFigure',f_step_zoom);
    t = 0:0.05:tsimtr(2);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out, ~] = lsim(Gap_theta, u, t);
    plot(t_out, y, syms(i)); hold on
    ax = gca;
    ax.FontSize = 13; 
    title('Respuesta Rampa Saturada','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(3),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    
    
    
%        Zoom en perturbación
%     indexOfInterest = (t < 5) & (t > 0); % range of t near perturbation
%     ax2.plot(t(indexOfInterest),y(indexOfInterest),'k'); hold on % plot on new axes
%     grid on;
%     axis tight
    
    %Calculo de margenes de ganancia y fase
    [Gm,Pm,Wcg,Wcp] = margin(Gap_theta_ol);
    Gms_d(i) = Gm; Pms_d(i) = Pm; wp_d(i) = Wcp; 
    
    % Calculo de Tiempos característicos
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf("Kp = "+ k_p + ", Kd = " + k_d(i) + ", Ki = " + k_i + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");
end

set(0, 'CurrentFigure', f_step)
    grid minor
    hold on
    box on
set(0, 'CurrentFigure', f_nichols)
    grid minor
    hold on
    box on


%% --6.3 Respuestas del resto de variables al AP


