% ----------------------------------------------------------------------- %
%                            SCVu 2021-2022                               %
%                          TRABAJO 2. GRUPO 04                            %
% ----------------------------------------------------------------------- %  
% DiseÃ±o de un sistema de aumento de estabilidad y un autopiloto para el
% RPAS Global Hawk. Canal longitudinal.

clear, clc, close all

% ------ DESCRIPCION DEL CODIGO ------

% - GlobalHawk.m/Learjet24.m : contienen las derivadas de estabilidad de
% dichos aviones. Global Hawk es el avion del Trabajo 1 y Learjet 24 se uso
% para comparar con el libro Leyes de Control de Vuelo (FTs/respuestas
% adimensionales).

% - PLANE.M: crea una clase a partir de las subclases anteriores para
% modelos de aviones concretos que proyecta las derivadas en ejes adecuados 
% y llama al resto de funciones para operar sobre ellas (obtener funciones 
% de transferencia, factorizarlas, etc.). Es la parte mas importante.
% Para chequear funciones de transferencia: plane.lon/lat.G.Gsalidaentrada,
% dentro estan todos los parametros de cada funcion de transferencia.
% NOTA: por defecto las K estan en forma adimensional. Para obtenerlas en
% forma dimensional como K_dim, correr los puntos 1 de este script para
% longitudinal y lateral-direccional.
% NOTA 2: para obtener los parametros de respuesta (overshoot, maximos,
% minimos, settling time...) hay que ejecutar los apartados de respuesta a
% escalon y rampa (5 y 6 de longitudinal y lateral-direccional).

% - Dynamics.m: contiene los sistemas de ecuaciones dinamicos para cada
% canal por separado. Resuelve las ecuaciones por el metodo de Cramer y
% devuelve las funciones de transferencia separadas por canal, entrada y
% salida, factorizadas y con los datos de polos, ceros y ganancias.

% - Utils.m: recoge funciones de utilidad variada (conversion de unidades,
% ISA, factorizacion de funciones de transferencia, etc.).

%% ------ Creacion del modelo dinamico y obtencion de FTs ------
plane_OL = Plane(GlobalHawk);

%% -------------------- CANAL LONGITUDINAL ---------------------
% -- 1. FTs adimensionales: contenidas en plane.lon.G --
% AÃ±adir un factor de escala a cada FT para tener en cuenta las unidades.

% - Dimensiones -
% Si dim = false, todo adimensional y angulos en radianes
% Si dim = true, velocidad en m/s y angulos en grados
% Se ajustan los factores de escala (scales) y los titulos de los diagramas
% Para el trabajo, debe estar por defecto en true. Para comparar con el
% libro en false
% Si save = true, guarda graficas en formato .eps en ./Graficas/Lon/

dim = true;
save = false;

varnames = ["udeltae"; "alphadeltae"; "thetadeltae"; "qdeltae"];
Gnames = strcat('G', varnames);

if dim == false
    scales = [1; 1; 1; 1];
    % - Parametros diagramas de Bode/Nichols - 
    Magnames_latex  = ["$\left | G_{\hat{u} \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\alpha \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\theta \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\hat{q} \delta_{e}}(i\omega) \right |$"];
    Gnames_latex    = ["$G_{\hat{u} \delta_{e}}(i\omega)$";...
                       "$G_{\alpha \delta_{e}}(i\omega)$";...
                       "$G_{\theta \delta_{e}}(i\omega)$";...
                       "$G_{\hat{q} \delta_{e}}(i\omega)$"];
                
    % - Parametros graficas de respuesta -
    Respnames_latex = ["$\Delta \hat{u}$ [-]";...
                       "$\Delta \alpha$ [rad]";...
                       "$\Delta \theta$ [rad]";...
                       "$\Delta \hat{q}$ [-]"];
    angleunits      = '\hspace{0.75 mm}rad$';
elseif dim == true
	scales = [plane_OL.FC.us/(180/pi);...  % Velocidad/angulo
              1;...                     % Angulo/angulo
              1;...                     % Angulo/angulo
              1/plane_OL.lon.t_lon];       % Vel.angular/angulo (tiempo)

    for i=1:length(scales)
        plane_OL.lon.G.(genvarname(Gnames(i))).K_dim = plane_OL.lon.G.(genvarname(Gnames(i))).K * scales(i);
    end

    % - Parametros diagramas de Bode/Nichols - 
    Magnames_latex  = ["$\left | G_{u \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\alpha \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{\theta \delta_{e}}(i\omega) \right |$";...
                       "$\left | G_{q \delta_{e}}(i\omega) \right |$"];
    Gnames_latex    = ["$G_{u \delta_{e}}(i\omega)$";...
                       "$G_{\alpha \delta_{e}}(i\omega)$";...
                       "$G_{\theta \delta_{e}}(i\omega)$";...
                       "$G_{q \delta_{e}}(i\omega)$"];
                
    % - Parametros graficas de respuesta -
    Respnames_latex = ["$\Delta u$ [m/s]";...
                       "$\Delta \alpha$ [$^\circ$]";...
                       "$\Delta \theta$ [$^\circ$]";...
                       "$\Delta q$ [$^\circ$/s]"];
    angleunits      = '^{\circ}$';
end

%% -- 2. Aplicacion de ley de control directo

% Esquema: delta --> K_DL --> G_actuador --> G_longitudinal --> variable
% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulacion - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulacion - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [deg]
G_a = Utils.forlagTF;

% Se seleccionara una ganancia de Direct-Link que proporcione un incremento
% de angulo de asiento de 1 grado por cada grado de deflexion del mando de cabeceo

K_dl = 1/plane_OL.lon.G.Gthetadeltae.K; % Calculamos el K_DL para theta


fprintf("--------- Tiempos Caracteristicos Lazo Abierto con K_dl + Gact + Gp -----------\n");

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
    set(gca,'TickLabelInterpreter','latex')
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Zoom en perturbacion
    if i ~=2
        axes('Position',[.30 .67 .15 .15])
    else
        axes('Position',[.30 .20 .15 .15])
    end
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
    set(gca,'TickLabelInterpreter','latex')
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane_OL.lon.G.(genvarname(Gnames(i))).ramp = lsiminfo(y, t_out);

    if save == false
        sgtitle("")
        %sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_e=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/SAS/Ramp_PL_", varnames(i)),'eps')
    end
    
    % Rise Time and Time Delay Calculation and printing in console
    
    [tr,td] = Utils.getTimesRD(y,t);
    TrTd_Kdl{i} = [tr,td];
    fprintf(varnames(i) + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");

end

%% -- 3. Analisis de sensibilidad de la planta libre --
% Analisis del lugar de las raices del modo de corto periodo
% Tarda mucho, ejecutar solo si es necesario.
clear
plane_OL  = Plane(GlobalHawk);
plane_Aux = Plane(GlobalHawk);

% Factores de multiplicidad Cmalpha y Cmq
F_alpha  = [-0.5 0 0.5 1 2 3 4];
F_q      = [-0.5 0 0.5 1 2 3 4];

% Root locus
figure('Position', [10 10 900 400])
polesAUX = pole(plane_OL.lon.G.Galphadeltae.Gfact);
poles_SP_0 = polesAUX(1);
plot([poles_SP_0 conj(poles_SP_0)],'x','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','w')
xlabel('Re(s)','Interpreter','latex')
ylabel('Im(s)','Interpreter','latex')
xlim([-0.3 0.5])
ylim([-1.5 1.5])
hold on
grid minor

for i=1:length(F_alpha)
    plane_Aux.model.Cm.alpha = plane_OL.model.Cm.alpha*F_alpha(i);
    for j=1:length(F_q)
        plane_Aux.model.Cm.q = plane_OL.model.Cm.q*F_q(j);
        plane_Aux = plane_Aux.recalc;
        polesAUX = pole(plane_Aux.lon.G.Galphadeltae.Gfact);
        poles_SP{i,j} = polesAUX;
        color = 'w';
        marker = 's';
        if F_alpha(i)==1 && F_q(j)==1
        else
            for k=1:length(polesAUX)
                if abs(imag(polesAUX(k))) > 0
                    plot(polesAUX(k),marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
                    hold on
        
                else
                    plot(polesAUX(k),0,marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
                    hold on
                end
            end
        end
    end
end
sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[1 4 8 12 16 20 24])

%% -- 4. Barrido de ganancias de realimentacion --
% Analisis del lugar de las raices del modo de corto periodo con
% realimentacion en alpha y q.
% Analisis de la funcion de transferencia de deltae a theta en lazo
% cerrado.

plane_OL = Plane(GlobalHawk);
save = true;

% Funciones de transferencia y Kdl
Ga_deltae = Utils.forlagTF;
Gs_alpha = 1;
Gs_q = 1;
Gthetadeltae_pl = plane_OL.lon.G.Gthetadeltae.Gfact;
Galphadeltae_pl = plane_OL.lon.G.Galphadeltae.Gfact;
Gqdeltae_pl     = plane_OL.lon.G.Gqdeltae.Gfact;
K_dl = 1/plane_OL.lon.G.Gthetadeltae.K;

% Factores de multiplicidad Cmalpha y Cmq
F_alpha  = [-0.5 -0.2 0 0.5 1 2 3 4];
% F_alpha = 1;
F_q      = [-0.5 -0.2 0 0.5 1 2 3 4];
% F_q = 1;
Kalpha   = - (F_alpha-1)*plane_OL.lon.Cm.alpha/plane_OL.lon.Cm.deltae;
Kq       = - (F_q-1)*plane_OL.lon.Cm.q/plane_OL.lon.Cm.deltae*plane_OL.lon.t_lon;
poles_SP = zeros(length(Kalpha),length(Kq));

% Root locus. Ploteamos primero el corto periodo, fugoide y polo del
% actuador base, sin controlar
figure('Position', [10 10 900 500])
polesAUX_0 = round(pole(plane_OL.lon.G.Gthetadeltae.Gfact),4);
zerosAUX_0 = round(zero(plane_OL.lon.G.Gthetadeltae.Gfact),4);
polesAUX_0 = sort(setdiff(polesAUX_0,zerosAUX_0));
for k=1:length(polesAUX_0)
    if abs(imag(polesAUX_0(k))) > 0
        plot(polesAUX_0(k),'x','MarkerSize',12,'MarkerEdgeColor','k')
        hold on       
    else
        plot(polesAUX_0(k),0,'x','MarkerSize',12,'MarkerEdgeColor','k')    
        hold on
    end
end
xlabel('Re(s)','Interpreter','latex','FontSize', 14)
ylabel('Im(s)','Interpreter','latex','FontSize', 14)
set(gca,'TickLabelInterpreter','latex')
xlim([-18 3])
ylim([-16 16]) 
hold on
grid minor

for i=1:length(Kalpha)
    for j=1:length(Kq)
        % -- Ensamblaje de la FT en lazo cerrado --
        Gthetadeltae_CL{i,j} = K_dl*(Ga_deltae*Gthetadeltae_pl)/(1+Ga_deltae*Kalpha(i)*Gs_alpha*Galphadeltae_pl +...
                                     Kq(j)*Gs_q*Gqdeltae_pl/plane_OL.lon.t_lon);
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
                if abs(imag(polesAUXclean(k))) > 0
                    plot(polesAUXclean(k),marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
                    hold on
                elseif real(polesAUXclean(k)) == -10
                else
                    plot(polesAUXclean(k),0,marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
                    hold on
                end
            end
        end
    end
end
sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20 24])

h = zeros(5, 1);
h(1) = plot(NaN,NaN,'x','MarkerEdgeColor','k');
h(2) = plot(NaN,NaN,'s','MarkerFaceColor','w','MarkerEdgeColor','k');
h(3) = plot(NaN,NaN,'s','MarkerFaceColor',[197/255 197/255 197/255],'MarkerEdgeColor','k');
h(4) = plot(NaN,NaN,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
h(5) = plot(NaN,NaN,'^','MarkerFaceColor',[197/255 197/255 197/255],'MarkerEdgeColor','k');
legend(h, 'Planta libre','$k_{\delta_e \alpha} < 0$, $k_{\delta_e q} < 0$','$k_{\delta_e \alpha} < 0$, $k_{\delta_e q} >= 0$',...
                           '$k_{\delta_e \alpha} >= 0$, $k_{\delta_e q} < 0$','$k_{\delta_e \alpha} >= 0$, $k_{\delta_e q} >= 0$', 'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);

% --- Diagramas de Nichols
lims = [1e-2 1e2];
nopts = nicholsoptions;
	nopts.PhaseMatching = 'on';
	nopts.PhaseMatchingFreq = lims(1);
	nopts.PhaseMatchingValue = 0;
    
	nopts.XLabel.String = '$\phi(\omega)$';
	nopts.XLabel.FontSize = 13;
	nopts.XLabel.Interpreter = 'latex';
    nopts.TickLabel.FontSize = 10;
    
	nopts.YLabel.FontSize = 13;
	nopts.YLabel.Interpreter = 'latex';

	nopts.Title.FontSize = 2;
	nopts.Title.Interpreter = 'latex';
    
    nopts.YLabel.String = '$\left | G_{OL}(i\omega) \right |$';
    nopts.Title.String = "";
figure
% - Barrido Kalpha -
Kalpha_nichols = [4.1 2 -2 -4.1];
Kq_nichols = 0;
lines = {'-' '--' '-.' ':'};
for i=1:length(Kalpha_nichols)
    for j=1:length(Kq_nichols)
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha_nichols(i)*Gs_alpha*plane_OL.lon.G.Galphadeltae.Gfact +...
                            Kq_nichols(j)*Gs_q*plane_OL.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    nicholsplot(Gthetadeltae_SAS,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
    end
end
grid on
legend('$k\_{\delta\_e \alpha} = 4.1$','$k\_{\delta\_e \alpha} = 2$','$k\_{\delta\_e \alpha} = -2$','$k\_{\delta\_e \alpha} = -4.1$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/SAS/Barrido_Nichols_Kalpha",'eps')
end

% - Barrido Kq -
Kq_nichols = [0.14 0.1 -0.1 -0.18];
Kalpha_nichols = 0;
figure
lines = {'-' '--' '-.' ':'};
for i=1:length(Kalpha_nichols)
    for j=1:length(Kq_nichols)
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha_nichols(i)*Gs_alpha*plane_OL.lon.G.Galphadeltae.Gfact +...
                            Kq_nichols(j)*Gs_q*plane_OL.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    nicholsplot(Gthetadeltae_SAS,strcat('k',lines{j}), {lims(1), lims(2)},nopts); hold on
    end
end
grid on
legend('$k\_{\delta\_e q} = 0.14$','$k\_{\delta\_e q} = 0.1$','$k\_{\delta\_e q} = -0.1$','$k\_{\delta\_e q} = -0.18$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/SAS/Barrido_Nichols_Kq",'eps')
end

% - Combinaciones -
Kq_nichols = [-0.1 -0.1 -0.3 -0.3];
Kalpha_nichols = [0.3 0.5 0.3 0.5];
figure
lines = {'-' '--' '-.' ':'};
for i=1:length(Kalpha_nichols)
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha_nichols(i)*Gs_alpha*plane_OL.lon.G.Galphadeltae.Gfact +...
                            Kq_nichols(i)*Gs_q*plane_OL.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    nicholsplot(Gthetadeltae_SAS,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('$\{k\_{\delta\_e \alpha},k\_{\delta\_e q}\}  = \{0.3, -0.1\}$','$\{k\_{\delta\_e \alpha},k\_{\delta\_e q}\}  = \{0.5, -0.1\}$',...
       '$\{k\_{\delta\_e \alpha},k\_{\delta\_e q}\}  = \{0.3, -0.2\}$','$\{k\_{\delta\_e \alpha},k\_{\delta\_e q}\}  = \{0.5, -0.2\}$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/SAS/Barrido_Nichols_Combinaciones",'eps')
end

%% -- 5. AnÃ¡lisis frente a tolerancias aerodinÃ¡micas --
save = true;
plane_OL  = Plane(GlobalHawk);
plane_Aux = Plane(GlobalHawk);

% Factores de multiplicidad Cmalpha y Cmq
F_alpha  = [0.8 1.2];
F_q  = [0.8 1.2];
F_deltae  = [0.8 1.2];

% Barrido (apartado anterior)
Falpha  = [-0.5 -0.2 0 0.5 1 2 3 4];
Fq      = [-0.5 -0.2 0 0.5 1 2 3 4];
K_alpha   = - (Falpha-1)*plane_OL.lon.Cm.alpha/plane_OL.lon.Cm.deltae;
K_q       = - (Fq-1)*plane_OL.lon.Cm.q/plane_OL.lon.Cm.deltae*plane_OL.lon.t_lon;

% Ganancias de realimentacion escogidas en el apartado anterior
Kalpha = 0.5;
Kq = -0.35;

% Funciones de transferencia y Kdl
Ga_deltae = Utils.forlagTF;
Gs_alpha = 1;
Gs_q = 1;
K_dl = 1/plane_OL.lon.G.Gthetadeltae.K;

% FunciÃ³n de transferencia closed-loop
Gthetadeltae_SAS_CL = K_dl*(Ga_deltae*plane_Aux.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                         Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);

% FunciÃ³n de transferencia open-loop
Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                         Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);

[GM_base, PM_base] = margin(Gthetadeltae_SAS);
GM_base_cumple = 'CUMPLE';
GM_base = mag2db(GM_base);
if GM_base < 6 || PM_base < 45
    GM_base_cumple = 'NO CUMPLE';
end
fprintf('----------- Robustez frente a tolerancias aerodinÃ¡micas ------------\n')
fprintf(strcat('Kalpha = ',num2str(Kalpha),'    Kq =',num2str(Kq),'\n'))
fprintf('--------- BASE ----------\n')
fprintf(strcat('GM = ',num2str(round(GM_base,4)),'  PM -20 = ',num2str(round(PM_base,4)),'\n'))
fprintf(strcat(GM_base_cumple,'\n'))

% Root locus
figure('Position', [10 10 900 400])
polesAUX_0 = round(pole(Gthetadeltae_SAS_CL),4);
zerosAUX_0 = round(zero(Gthetadeltae_SAS_CL),4);
polesAUX_0 = sort(setdiff(polesAUX_0,zerosAUX_0));
for k=1:length(polesAUX_0)
    if abs(imag(polesAUX_0(k))) > 0
        plot(polesAUX_0(k),'x','MarkerSize',12,'MarkerEdgeColor','k')
        hold on       
    else
        plot(polesAUX_0(k),0,'x','MarkerSize',12,'MarkerEdgeColor','k')    
        hold on
    end
end
% Obtenemos el amortiguamiento y frecuencia de cada modo base
count = 1;
tau_p = [];
wn_p = [];
xi_p = [];
while count <= length(polesAUX_0)
    if imag(polesAUX_0(count)) == 0
        if polesAUX_0(count) ~= 0
            tau_p(end+1) = -1/polesAUX_0(count);
            count = count + 1;
        else
            null_p = null_p + 1;
            count = count + 1;
        end
    else
        lambda_p = polesAUX_0(count);
        n_p = real(lambda_p);
        w_p = abs(imag(lambda_p));
        wn_p(end+1) = sqrt(n_p^2 + w_p^2);
        xi_p(end+1) = -n_p/sqrt(n_p^2 + w_p^2);
        count = count + 2;
    end
end
polesAUX_0_xiwn.tau = tau_p;
polesAUX_0_xiwn.wn = wn_p;
polesAUX_0_xiwn.xi = xi_p;
xlabel('Re(s)','Interpreter','latex','FontSize', 14)
ylabel('Im(s)','Interpreter','latex','FontSize', 14)
set(gca,'TickLabelInterpreter','latex')
xlim([-12 2])
ylim([-19.5 19.5])
hold on
grid minor
marker = ['v' '^'];

GM_alpha_cumple = 'CUMPLE';
for i=1:length(F_alpha)
    plane_Aux.model.Cm.alpha = plane_OL.model.Cm.alpha*F_alpha(i);
    plane_Aux = plane_Aux.recalc;
    Gthetadeltae_SAS_CL = K_dl*(Ga_deltae*plane_Aux.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                             Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                            Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    polesAUX = round(pole(Gthetadeltae_SAS_CL),4);
    zerosAUX = round(zero(Gthetadeltae_SAS_CL),4);
    polesAUX_alpha{i} = setdiff(polesAUX,zerosAUX);
    color = 'w';
    for k=1:length(polesAUX_alpha{i})
        if abs(imag(polesAUX_alpha{i}(k))) > 0
            plot(polesAUX_alpha{i}(k),marker(i),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
            hold on       
        else
            plot(polesAUX_alpha{i}(k),0,marker(i),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
            hold on
        end
    end
    [GM, PM] = margin(Gthetadeltae_SAS);
    GM_alpha(i) = mag2db(GM);
    PM_alpha(i) = PM;
    if GM_alpha(i) < 4.5 || PM_alpha(i) < 30
        GM_alpha_cumple = 'NO CUMPLE';
    end
end

GM_q_cumple = 'CUMPLE';
for i=1:length(F_q)
    plane_Aux.model.Cm.q = plane_OL.model.Cm.q*F_q(i);
    plane_Aux = plane_Aux.recalc;
    Gthetadeltae_SAS_CL = K_dl*(Ga_deltae*plane_Aux.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                             Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                             Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    polesAUX = round(pole(Gthetadeltae_SAS_CL),4);
    zerosAUX = round(zero(Gthetadeltae_SAS_CL),4);
    polesAUX_q{i} = setdiff(polesAUX,zerosAUX);
    color = [197/255 197/255 197/255];
    for k=1:length(polesAUX_q{i})
        if abs(imag(polesAUX_q{i}(k))) > 0
            plot(polesAUX_q{i}(k),marker(i),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
            hold on       
        else
            plot(polesAUX_q{i}(k),0,marker(i),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
            hold on
        end
    end
    [GM, PM] = margin(Gthetadeltae_SAS);
    GM_q(i) = mag2db(GM);
    PM_q(i) = PM;
    if GM_q(i) < 4.5 || PM_q(i) < 30
        GM_q_cumple = 'NO CUMPLE';
    end
end

GM_deltae_cumple = 'CUMPLE';
for i=1:length(F_deltae)
    plane_Aux.model.Cm.deltae = plane_OL.model.Cm.deltae*F_deltae(i);
    plane_Aux = plane_Aux.recalc;
    Gthetadeltae_SAS_CL = K_dl*(Ga_deltae*plane_Aux.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                             Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                            Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    polesAUX = round(pole(Gthetadeltae_SAS_CL),4);
    zerosAUX = round(zero(Gthetadeltae_SAS_CL),4);
    polesAUX_deltae{i} = setdiff(polesAUX,zerosAUX);
    color = 'k';
    for k=1:length(polesAUX_deltae{i})
        if abs(imag(polesAUX_deltae{i}(k))) > 0
            plot(polesAUX_deltae{i}(k),marker(i),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
            hold on       
        else
            plot(polesAUX_deltae{i}(k),0,marker(i),'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
            hold on
        end
    end
    [GM, PM] = margin(Gthetadeltae_SAS);
    GM_deltae(i) = mag2db(GM);
    PM_deltae(i) = PM;
    if GM_deltae(i) < 4.5 || PM_deltae(i) < 30
        GM_deltae_cumple = 'NO CUMPLE';
    end
end

for k=1:length(polesAUX_0)
    if abs(imag(polesAUX_0(k))) > 0
        plot(polesAUX_0(k),'x','MarkerSize',12,'MarkerEdgeColor','k')
        hold on       
    else
        plot(polesAUX_0(k),0,'x','MarkerSize',12,'MarkerEdgeColor','k')    
        hold on
    end
end

sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20 24])
h = zeros(6, 1);
h(1) = plot(NaN,NaN,'v','MarkerFaceColor','w','MarkerEdgeColor','k');
h(2) = plot(NaN,NaN,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
h(3) = plot(NaN,NaN,'v','MarkerFaceColor',[197/255 197/255 197/255],'MarkerEdgeColor','k');
h(4) = plot(NaN,NaN,'^','MarkerFaceColor',[197/255 197/255 197/255],'MarkerEdgeColor','k');
h(5) = plot(NaN,NaN,'v','MarkerFaceColor','k','MarkerEdgeColor','k');
h(6) = plot(NaN,NaN,'^','MarkerFaceColor','k','MarkerEdgeColor','k');
legend(h,'$0.8\cdot C_{m\alpha}$','$1.2\cdot C_{m\alpha}$', '$0.8\cdot C_{m\hat{q}}$','$1.2\cdot C_{m\hat{q}}$',...
                           '$0.8\cdot C_{m\delta_e}$','$1.2\cdot C_{m\delta_e}$', 'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);

% --- Chequeo de mÃ¡rgenes de estabilidad ---
fprintf('------ +-20 Cmalpha -------\n')
fprintf(strcat('GM -20 = ',num2str(round(GM_alpha(1),4)),'  PM -20 = ',num2str(round(PM_alpha(1),4)),'\n'))
fprintf(strcat('GM +20 = ',num2str(round(GM_alpha(2),4)),'  PM +20 = ',num2str(round(PM_alpha(2),4)),'\n'))
fprintf(strcat(GM_alpha_cumple,'\n'))
fprintf('------ +-20 Cmq -------\n')
fprintf(strcat('GM -20 = ',num2str(round(GM_q(1),4)),'  PM -20 = ',num2str(round(PM_q(1),4)),'\n'))
fprintf(strcat('GM +20 = ',num2str(round(GM_q(2),4)),'  PM +20 = ',num2str(round(PM_q(2),4)),'\n'))
fprintf(strcat(GM_q_cumple,'\n'))
fprintf('------ +-20 Cmdeltae -------\n')
fprintf(strcat('GM -20 = ',num2str(round(GM_deltae(1),4)),'  PM -20 = ',num2str(round(PM_deltae(1),4)),'\n'))
fprintf(strcat('GM +20 = ',num2str(round(GM_deltae(2),4)),'  PM +20 = ',num2str(round(PM_deltae(2),4)),'\n'))
fprintf(strcat(GM_deltae_cumple,'\n'))

% --- Diagramas de Nichols ---
lims = [1e-2 1e2];
nopts = nicholsoptions;
	nopts.PhaseMatching = 'on';
	nopts.PhaseMatchingFreq = lims(1);
	nopts.PhaseMatchingValue = 0;
    
	nopts.XLabel.String = '$\phi(\omega)$';
	nopts.XLabel.FontSize = 13;
	nopts.XLabel.Interpreter = 'latex';
    nopts.TickLabel.FontSize = 10;
    
	nopts.YLabel.FontSize = 13;
	nopts.YLabel.Interpreter = 'latex';

	nopts.Title.FontSize = 2;
	nopts.Title.Interpreter = 'latex';
    
    nopts.YLabel.String = '$\left | G_{OL}(i\omega) \right |$';
    nopts.Title.String = "";

% - Variacion Cmalpha -
figure
Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_OL.lon.G.Galphadeltae.Gfact +...
                         Kq*Gs_q*plane_OL.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
nicholsplot(Gthetadeltae_SAS,'k', {lims(1), lims(2)},nopts); hold on    
lines = {'--' '-.'};
for i=1:length(F_alpha)
    plane_Aux.model.Cm.alpha = plane_OL.model.Cm.alpha*F_alpha(i);
    plane_Aux = plane_Aux.recalc;
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                            Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    nicholsplot(Gthetadeltae_SAS,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('Base','$0.8\cdot C\_{m\_\alpha}$','$1.2\cdot C\_{m\_\alpha}$','Interpreter','latex','Location','south','FontSize',12)
if save == true
    saveas(gcf,"Graficas/SAS/Robustez_Nichols_Cmalpha",'eps')
end

% - Variacion Cmq -
figure
Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_OL.lon.G.Galphadeltae.Gfact +...
                         Kq*Gs_q*plane_OL.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
nicholsplot(Gthetadeltae_SAS,'k', {lims(1), lims(2)},nopts); hold on    
lines = {'--' '-.'};
for i=1:length(F_q)
    plane_Aux.model.Cm.q = plane_OL.model.Cm.q*F_q(i);
    plane_Aux = plane_Aux.recalc;
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                            Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    nicholsplot(Gthetadeltae_SAS,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('Base','$0.8\cdot C\_{m\_{\hat{q}}}$','$1.2\cdot C\_{m\_{\hat{q}}}$','Interpreter','latex','Location','south','FontSize',12)
if save == true
    saveas(gcf,"Graficas/SAS/Robustez_Nichols_Cmq",'eps')
end

% - Variacion Cmdeltae -
figure
Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_OL.lon.G.Galphadeltae.Gfact +...
                         Kq*Gs_q*plane_OL.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
nicholsplot(Gthetadeltae_SAS,'k', {lims(1), lims(2)},nopts); hold on    
lines = {'--' '-.'};
for i=1:length(F_deltae)
    plane_Aux.model.Cm.deltae = plane_OL.model.Cm.deltae*F_deltae(i);
    plane_Aux = plane_Aux.recalc;
    Gthetadeltae_SAS = K_dl*Ga_deltae*(Kalpha*Gs_alpha*plane_Aux.lon.G.Galphadeltae.Gfact +...
                            Kq*Gs_q*plane_Aux.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
    nicholsplot(Gthetadeltae_SAS,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('Base','$0.8\cdot C\_{m\_{\delta\_e}}$','$1.2\cdot C\_{m\_{\delta\_e}}$','Interpreter','latex','Location','south','FontSize',12)
if save == true
    saveas(gcf,"Graficas/SAS/Robustez_Nichols_Cmdeltae",'eps')
end

%% -- 5.1 DiseÃ±o del SAS --

save = true;
plane_OL  = Plane(GlobalHawk);
plane_SAS = Plane(GlobalHawk);
wn_FQ  = [sqrt(0.28*plane_OL.model.CL.alpha/plane_OL.FC.CLs) sqrt(3.6*plane_OL.model.CL.alpha/plane_OL.FC.CLs)];
xi_FQ = [0.35 1.30];

% Ganancias de realimentacion escogidas en el apartado anterior
Kalpha = 0.5;
Kq = -0.35;

% Funciones de transferencia y Kdl
Ga_deltae = Utils.forlagTF;
Gs_alpha = 1;
Gs_q = 1;
Gthetadeltae_pl = plane_OL.lon.G.Gthetadeltae.Gfact;
Galphadeltae_pl = plane_OL.lon.G.Galphadeltae.Gfact;
Gqdeltae_pl     = plane_OL.lon.G.Gqdeltae.Gfact;
K_dl = 1/plane_OL.lon.G.Gthetadeltae.K;

% Funciones de transferencia
plane_SAS.lon.Cont.Gudeltae.Gfact = K_dl*(Ga_deltae*plane_SAS.lon.G.Gudeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_SAS.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_SAS.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
plane_SAS.lon.Cont.Galphadeltae.Gfact = K_dl*(Ga_deltae*plane_SAS.lon.G.Galphadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_SAS.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_SAS.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
plane_SAS.lon.Cont.Gthetadeltae.Gfact = K_dl*(Ga_deltae*plane_SAS.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_SAS.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_SAS.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);
plane_SAS.lon.Cont.Gqdeltae.Gfact = K_dl*(Ga_deltae*plane_SAS.lon.G.Gqdeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_SAS.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_SAS.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);

% Guardamos las funciones de transferencia del SAS aparte                   
GSAS_u =  plane_SAS.lon.Cont.Gudeltae.Gfact/K_dl;
GSAS_alpha =  plane_SAS.lon.Cont.Galphadeltae.Gfact/K_dl;
GSAS_theta = plane_SAS.lon.Cont.Gthetadeltae.Gfact/K_dl;
GSAS_q =  plane_SAS.lon.Cont.Gqdeltae.Gfact/K_dl;

% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulacion - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulacion - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [deg]
G_a = Utils.forlagTF;

fprintf("--------- Tiempos Caracteristicos Lazo Cerrado con SAS en alpha y q -----------\n");

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
    [y, t_out, ~] = lsim(plane_SAS.lon.Cont.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    hold on
    [yorig, t_outorig, ~] = lsim(plane_SAS.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
    plot(t_outorig, yorig, 'k--')
    if i==1
        hold on
        plot(t,u,'k-.')
    end
    ax = gca;
    ax.FontSize = 13; 
    set(gca,'TickLabelInterpreter','latex')
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Zoom en perturbacion
    if i~=2
        axes('Position',[.30 .67 .15 .15])
    else
        axes('Position',[.30 .20 .15 .15])
    end
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
    [y, t_out, ~] = lsim(plane_SAS.lon.Cont.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i), u, t);
    plot(t_out, y, 'k')
    [yorig, t_outorig, ~] = lsim(plane_SAS.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
    plot(t_outorig, yorig, 'k--')
    if i==1
        hold on
        plot(t,u,'k-.')
        legend('SAS','Planta libre', '$\Delta \delta_e$', 'Interpreter', 'latex')
    elseif i==2
        legend('SAS','Planta libre', 'Interpreter', 'latex','Location','southeast')
    else
        legend('SAS','Planta libre', 'Interpreter', 'latex')
    end
    ax = gca;
    ax.FontSize = 13;
    set(gca,'TickLabelInterpreter','latex')
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    plane_SAS.lon.Cont.(genvarname(Gnames(i))).ramp = lsiminfo(y, t_out);

    if save == false
        sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_e=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/SAS/Ramp_SAS_", varnames(i)),'eps')
    end
    
    % Rise Time and Time Delay Calculation and printing in console
    
    [tr,td] = Utils.getTimesRD(y,t);
    TrTd_SAS{i} = [tr,td];
    fprintf(varnames(i) + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");

end
%%
figure
bode(plane_OL.lon.G.Gthetadeltae.Gfact*K_dl*Ga_deltae)
hold on
bode(GSAS_theta*K_dl)
hold on
legend('Kdl + Ga', 'Kdl + Ga + SAS')
grid on

%% -- 6.1. DiseÃ±o AP: Root Locus

% Valores PID
k_p_res = [0, -0.1, -0.2, -0.5, -1, -2];
k_i_res = [0, -0.01, -0.1, -0.2, -0.5, -1];
k_d_res = 0; % Lo ponemos a cero porque basicamente no afecta demasiado a la estabilidad del sistema
Gs_theta = Utils.padeTF();

% Valores SAS
Kalpha = 0.5;
Kq = -0.35;
GSAS_theta = (Ga_deltae*plane_SAS.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_SAS.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_SAS.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);

% Root locus
figure('Position', [10 10 900 400])
polesAUX_0 = round(pole(GSAS_theta),4);
zerosAUX_0 = round(zero(GSAS_theta),4);
polesAUX_0 = sort(setdiff(polesAUX_0,zerosAUX_0));
for k=1:length(polesAUX_0)
    if abs(imag(polesAUX_0(k))) > 0
        plot(polesAUX_0(k),'x','MarkerSize',12,'MarkerEdgeColor','k')
        hold on       
    else
        plot(polesAUX_0(k),0,'x','MarkerSize',12,'MarkerEdgeColor','k')    
        hold on
    end
end
xlabel('Re(s)','Interpreter','latex','FontSize', 14)
ylabel('Im(s)','Interpreter','latex','FontSize', 14)
set(gca,'TickLabelInterpreter','latex')
xlim([-12 2])
ylim([-19.5 19.5])
hold on
grid minor

for i=1:length(k_p_res)
    for j=1:length(k_i_res)
        % -- Ensamblaje de la FT en lazo cerrado --
        s = tf('s');
        PID_AP = pid(k_p_res(i),k_i_res(j),k_d_res);
        Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
        polesAUX = round(pole(Gap_theta),4);
        zerosAUX = round(zero(Gap_theta),4);
        % Quita a la lista de polos los ceros que sean iguales a estos.
        % Para eliminar los polos residuales de la planta libre.
        polesAUXclean = setdiff(polesAUX,zerosAUX);
        if k_p_res(i) < -0.5
            marker = 's';
        else
            marker = '^';
        end
        if k_i_res(j) < -0.5
            color = 'w';
        else
            color = [197/255 197/255 197/255];
        end
%         marker = 's';
%         color = 'w';
        for k=1:length(polesAUXclean)
            % -- Bucle en todos los polos --
            % Filtramos los polos cerca del origen (modo fugoide) a
            % mano. Mejorable.
            if abs(imag(polesAUXclean(k))) > 0
                plot(polesAUXclean(k),marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
                hold on
            elseif real(polesAUXclean(k)) == -10
            else
                plot(polesAUXclean(k),0,marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
                hold on
            end
        end
    end
end
sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20 24])
h = zeros(5, 1);
h(1) = plot(NaN,NaN,'x','MarkerEdgeColor','k');
h(2) = plot(NaN,NaN,'s','MarkerFaceColor','w','MarkerEdgeColor','k');
h(3) = plot(NaN,NaN,'s','MarkerFaceColor',[197/255 197/255 197/255],'MarkerEdgeColor','k');
h(4) = plot(NaN,NaN,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
h(5) = plot(NaN,NaN,'^','MarkerFaceColor',[197/255 197/255 197/255],'MarkerEdgeColor','k');
legend(h, 'Planta libre','$k_p < -0.5$, $k_i < -0.5$','$k_p < -0.5$, $k_i >= -0.5$',...
                           '$k_p >= -0.5$, $k_i < -0.5$','$k_p >= -0.5$, $k_i >= -0.5$', 'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);

%% -- 6.1.1 DiseÃ±o AP: Root Locus Kd
% Fijando Kp y Ki, variar Kd para ver la recuperaciÃ³n de amortiguamiento SP

% Valores PID
k_p_res = -1;
k_i_res = -0.2;
k_d_res = [-0.2 -0.1 -0.05]; % Lo ponemos a cero porque basicamente no afecta demasiado a la estabilidad del sistema
Gs_theta = Utils.padeTF();

% Valores SAS
Kalpha = 0.5;
Kq = -0.35;
GSAS_theta = (Ga_deltae*plane_SAS.lon.G.Gthetadeltae.Gfact)/(1+Ga_deltae*Kalpha*Gs_alpha*plane_SAS.lon.G.Galphadeltae.Gfact +...
                                     Kq*Gs_q*plane_SAS.lon.G.Gqdeltae.Gfact/plane_OL.lon.t_lon);

% Root locus sin AP
figure('Position', [10 10 900 400])
polesAUX_0 = round(pole(GSAS_theta),4);
zerosAUX_0 = round(zero(GSAS_theta),4);
polesAUX_0 = sort(setdiff(polesAUX_0,zerosAUX_0));
for k=1:length(polesAUX_0)
    if abs(imag(polesAUX_0(k))) > 0
        plot(polesAUX_0(k),'x','MarkerSize',12,'MarkerEdgeColor','k')
        hold on       
    else
        plot(polesAUX_0(k),0,'x','MarkerSize',12,'MarkerEdgeColor','k')    
        hold on
    end
end
xlabel('Re(s)','Interpreter','latex','FontSize', 14)
ylabel('Im(s)','Interpreter','latex','FontSize', 14)
set(gca,'TickLabelInterpreter','latex')
xlim([-12 2])
ylim([-19.5 19.5])
hold on
grid minor

% Root locus con AP base
PID_AP = pid(k_p_res,k_i_res,0);
Gap_theta_0 = feedback(GSAS_theta*PID_AP,Gs_theta);
polesAUX_0_KpKi = round(pole(Gap_theta_0),4);
zerosAUX_0_KpKi = round(zero(Gap_theta_0),4);
polesAUX_0_KpKi = sort(setdiff(polesAUX_0_KpKi,zerosAUX_0_KpKi));
for k=1:length(polesAUX_0_KpKi)
    if abs(imag(polesAUX_0_KpKi(k))) > 0
        plot(polesAUX_0_KpKi(k),'^','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k')
        hold on       
    else
        plot(polesAUX_0_KpKi(k),0,'^','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k')    
        hold on
    end
end

for i=1:length(k_d_res)
    PID_AP = pid(k_p_res,k_i_res,k_d_res(i));
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    polesAUX = round(pole(Gap_theta),4);
    zerosAUX = round(zero(Gap_theta),4);
    % Quita a la lista de polos los ceros que sean iguales a estos.
    % Para eliminar los polos residuales de la planta libre.
    polesAUXclean = setdiff(polesAUX,zerosAUX);
    marker = 's';
    color = 'w';
    for k=1:length(polesAUXclean)
        % -- Bucle en todos los polos --
        % Filtramos los polos cerca del origen (modo fugoide) a
        % mano. Mejorable.
        if abs(imag(polesAUXclean(k))) > 0
            plot(polesAUXclean(k),marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)
            hold on
        elseif real(polesAUXclean(k)) == -10
        else
            plot(polesAUXclean(k),0,marker,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color)    
            hold on
        end
    end
end

sgrid([0.1 0.2 0.3 0.4 0.5 0.6 0.8],[4 8 12 16 20 24])
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'x','MarkerEdgeColor','k');
h(2) = plot(NaN,NaN,'^','MarkerFaceColor','k','MarkerEdgeColor','k');
h(3) = plot(NaN,NaN,'s','MarkerFaceColor','w','MarkerEdgeColor','k');
legend(h, 'Sin AP','$k_d = 0$','$k_d < 0$', 'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);

%% -- 6.2 DiseÃ±o de AP: PID
save = false;
% Datos preliminares para respuesta en rampa
tsimtr = [150 50 100 100];    % Tiempo total de simulacion - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulacion - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [deg]

syms = ["-","--",":","-."];

lims = [1e-2 1e2];

% Valores PID
k_p_res = [-0.1, -0.2, -0.5, -1, -2];
k_i_res = [-0.01, -0.1, -0.2, -0.5, -1];
k_d_res = 0;
Gs_theta = Utils.padeTF();

% Margenes y fases

Gms_p = zeros(1,length(k_p_res));
Pms_p = zeros(1,length(k_p_res));
wp_p = zeros(1,length(k_p_res));

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
    
    nopts.YLabel.String = '$\left | G_{OL}(i\omega) \right |$';
    nopts.Title.String = '';


% - Barrido Kp -
k_p_nichols = [-0.5, -0.5, -1, -1];
figure
lines = {'-' '--' '-.' ':'};
for i=1:length(k_p_nichols)
    PID_AP = pid(k_p_nichols(i),0,0);
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;  
    nicholsplot(Gap_theta_ol,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('$k\_{p} = -0.1$','$k\_{p} = -0.2$','$k\_{p} = -0.5$','$k\_{p} = -1$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Nichols_Kp",'eps')
end

% - Barrido Ki -
k_p_nichols = -1;
k_i_nichols = [-0.01, -0.1, -0.2, -0.5];
figure
lines = {'-' '--' '-.' ':'};
for i=1:length(k_i_nichols)
    PID_AP = pid(k_p_nichols,k_i_nichols(i),0);
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;  
    nicholsplot(Gap_theta_ol,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('$k\_{i} = -0.01$','$k\_{i} = -0.1$','$k\_{i} = -0.2$','$k\_{i} = -0.5$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Nichols_Ki",'eps')
end

% - Barrido Kp-Ki -
k_p_nichols = [-0.5, -0.5, -1, -1];
k_i_nichols = [-0.1, -0.2, -0.1, -0.2];
figure
lines = {'-' '--' '-.' ':'};
for i=1:length(k_p_nichols)
    PID_AP = pid(k_p_nichols(i),k_i_nichols(i),0);
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;  
    nicholsplot(Gap_theta_ol,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('$\{k\_p,k\_i\}  = \{-0.5, -0.1\}$','$\{k\_p,k\_i\}  = \{-0.5, -0.2\}$',...
       '$\{k\_p,k\_i\}  = \{-1, -0.1\}$','$\{k\_p,k\_i\}  = \{-1, -0.2\}$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Nichols_KpKi",'eps')
end

% - Barrido Kd -
% Habiendo escogido kp = -1, ki = -0.2, barremos en kd para mejorar el
% amortiguamiento
figure
k_d_nichols = [-0.2 -0.1 -0.05];
lines = {'-' '--' '-.' ':'};
for i=1:3
    PID_AP = pid(-1,-0.2,k_d_nichols(i));
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;  
    nicholsplot(Gap_theta_ol,strcat('k',lines{i}), {lims(1), lims(2)},nopts); hold on
end
grid on
legend('$k\_{d} = -0.2$','$k\_{d} = -0.1$','$k\_{d} = -0.05$','Interpreter','latex','Location','northoutside','NumColumns',2,'FontSize',12)
if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Nichols_Kd",'eps')
end

% -- Chequeo de margenes --
k_p_check = -1;
k_i_check = -0.2;
k_d_check = -0.1;
PID_AP = pid(k_p_check,k_i_check,k_d_check);
Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;  
[GM_base, PM_base] = margin(Gap_theta_ol);
GM_base_cumple = 'CUMPLE';
GM_base = mag2db(GM_base);
if GM_base < 6 || PM_base < 45
    GM_base_cumple = 'NO CUMPLE';
end
fprintf('----------- Chequeo PID ------------\n')
fprintf(strcat('Kp = ',num2str(k_p_check),'    Ki =',num2str(k_i_check),'    Kd =',num2str(k_d_check),'\n'))
fprintf(strcat('GM = ',num2str(round(GM_base,4)),'  PM = ',num2str(round(PM_base,4)),'\n'))
fprintf(strcat(GM_base_cumple,'\n'))


% ------ Gráficas de respuesta -------

% 1) Proporcional 
figure
% Valores PID
k_p_res = [-0.1, -0.2, -0.5, -1];
k_i_res = 0;
k_d_res = 0;
Gs_theta = Utils.padeTF();

for i=1:length(k_p_res)

    PID_AP = pid(k_p_res(i),k_i_res,k_d_res);
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;
    
    % Step Plotting
    % Transitorio
    t = 0:0.05:tsimtr(3);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out] = lsim(Gap_theta, u, t);
    plot(t_out,y,'LineStyle',syms(i),'Color','k'); hold on
    set(gca,'TickLabelInterpreter','latex')
    ax = gca;
    ax.FontSize = 12; 
    ylabel(Respnames_latex(3),'FontSize', 14, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 14, 'Interpreter', 'latex')
    legend(strcat('$k_p=',num2str(k_p_res(1)),'$'),strcat('$k_p=',num2str(k_p_res(2)),'$'),strcat('$k_p=',num2str(k_p_res(3)),'$'),...
           strcat('$k_p=',num2str(k_p_res(4)),'$'),'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);
    grid minor
    hold on

    % Sacamos los valores de margen de fase y ganancia
    [Gm,Pm,Wcg,Wcp] = margin(Gap_theta_ol);
    Gms_p(i) = mag2db(Gm); Pms_p(i) = Pm; wp_p(i) = Wcp; 
    
    % Calculo de Tiempos caracteristicos
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf("------ Tiempos caracteristicos -------\n");
    fprintf("Kp = "+ k_p_res(i) + ", Kd = " + k_d_res + ", Ki = " + k_i_res + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");
    
end
if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Respuesta_Kp",'eps')
end

% 2) Integrador
figure
% Valores PID
k_p_res = -1;
k_i_res = [-0.01, -0.1, -0.2, -0.5];
k_d_res = 0;
Gs_theta = Utils.padeTF();

% Margenes y fases

Gms_i = zeros(1,length(k_i_res));
Pms_i = zeros(1,length(k_i_res));
wp_i = zeros(1,length(k_i_res));


for i=1:length(k_i_res)
    
    PID_AP = pid(k_p_res,k_i_res(i),k_d_res);
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;
    
    % Step Plotting
    % Transitorio
    t = 0:0.05:tsimtr(4);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out, ~] = lsim(Gap_theta, u, t);
    plot(t_out, y, syms(i),'Color','k'); hold on
    set(gca,'TickLabelInterpreter','latex')
    ax = gca;
    ax.FontSize = 12; 
    ylabel(Respnames_latex(3),'FontSize', 14, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 14, 'Interpreter', 'latex')
    legend(strcat('$k_i=',num2str(k_i_res(1)),'$'),strcat('$k_i=',num2str(k_i_res(2)),'$'),strcat('$k_i=',num2str(k_i_res(3)),'$'),...
           strcat('$k_i=',num2str(k_i_res(4)),'$'),'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);
    grid minor
    hold on
    % Calculo de margenes de ganancia y fase
    [Gm,Pm,Wcg,Wcp] = margin(Gap_theta_ol);
    Gms_i(i) = mag2db(Gm); Pms_i(i) = Pm; wp_i(i) = Wcp; 
    
    % Calculo de Tiempos caracteristicos
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf("------ Tiempos caracteristicos -------\n");
    fprintf("Kp = "+ k_p_res + ", Kd = " + k_d_res + ", Ki = " + k_i_res(i) + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");
end

if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Respuesta_Ki",'eps')
end

% 3) Derivativo
figure
% Valores PID
k_p_res = -1;
k_i_res = -0.2;
k_d_res = [-0.05 -0.1 -0.2];
Gs_theta = Utils.padeTF();

% Margenes y fases

Gms_d = zeros(1,length(k_d_res));
Pms_d = zeros(1,length(k_d_res));
wp_d = zeros(1,length(k_d_res));


for i=1:length(k_d_res)
    
    PID_AP = pid(k_p_res,k_i_res,k_d_res(i));
    Gap_theta = feedback(GSAS_theta*PID_AP,Gs_theta);
    Gap_theta_ol = PID_AP*(GSAS_theta)*Gs_theta;
    
    % Step Plotting
    % Transitorio
    t = 0:0.05:tsimtr(4);
    u = max(0,min(amp/tramp*(t),amp));
    [y, t_out, ~] = lsim(Gap_theta, u, t);
    plot(t_out, y, syms(i),'Color','k'); hold on
    set(gca,'TickLabelInterpreter','latex')
    ax = gca;
    ax.FontSize = 12; 
    ylabel(Respnames_latex(3),'FontSize', 14, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 14, 'Interpreter', 'latex')
    legend(strcat('$k_d=',num2str(k_d_res(1)),'$'),strcat('$k_d=',num2str(k_d_res(2)),'$'),strcat('$k_d=',num2str(k_d_res(3)),'$'),...
           'Interpreter', 'latex', 'Location','northoutside', 'NumColumns', 3, 'Fontsize', 12);
    grid minor
    hold on
    % Calculo de margenes de ganancia y fase
    [Gm,Pm,Wcg,Wcp] = margin(Gap_theta_ol);
    Gms_d(i) = mag2db(Gm); Pms_d(i) = Pm; wp_d(i) = Wcp; 
   
end

if save == true
    saveas(gcf,"Graficas/AP/AP_Barrido_Respuesta_Kd",'eps')
end

%% --6.3 Respuestas del resto de variables al AP

k_p_res = -1;
k_i_res = -0.2;
k_d_res = -0.1;

% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulacion - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulacion - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [deg]
delay = 0.1;                  % Delay del actuador

PID_AP = pid(k_p_res,k_i_res,k_d_res);
Gs_theta = Utils.padeTF();
load('Gaps.mat');
GSASs = [Gap_u,Gap_alpha,GSAS_theta,Gap_q];

fprintf("--------- Tiempos Caracteristicos AP para cada variable -----------\n");

for i = 1:length(varnames)
    fig = figure('Position', [100, 100, 1120, 420]);
    
    if(i == 3)
        Gap = feedback(GSASs(i)*PID_AP,Gs_theta);
    else
        Gap = GSASs(i);
    end
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(Gap*scales(i), u, t);
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
    
    % Zoom en perturbacion
    
    if(i == 4 || i == 1)
        axes('Position',[.30 .67 .15 .15])
        box on; % put box around new pair of axes
        indexOfInterest = (t < 4) & (t > 0); % range of t near perturbation
        plot(t(indexOfInterest),y(indexOfInterest),'k'); hold on % plot on new axes
        grid on;
        axis tight
    end
    if(i == 3 || i == 2)
        axes('Position',[.3 .2 .15 .15])
        box on; % put box around new pair of axes
        indexOfInterest = (t < 4) & (t > 0); % range of t near perturbation
        plot(t(indexOfInterest),y(indexOfInterest),'k'); hold on % plot on new axes
        grid on;
        axis tight
    end
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    [y, t_out, ~] = lsim(Gap*scales(i), u, t);
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

    if (save == false)
        sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_e=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/AP/Ramp_", varnames(i),"AP_ONLY"),'eps')
    end
    
    % Rise Time and Time Delay Calculation and printing in console
    
    [tr,td] = Utils.getTimesRD(y,t);
    fprintf(varnames(i) + "\n");
    fprintf("Rise Time: " + num2str(tr) + " s\n");
    fprintf("Time Delay: " + num2str(td) + " s\n");

end

%% --7. Comparación Final SAS vs No SAS vs AP

k_p_res = -1;
k_i_res = -0.2;
k_d_res = -0.1;

syms = ["-","--","-."];

% Definir rampa
tsimtr = [150 50 150 100];    % Tiempo total de simulacion - transitorio [s]
tsimst = [500 500 500 500];   % Tiempo total de simulacion - estacionario [s]
tramp  = 5;                   % Tiempo de rampeo [s]
amp    = 1;                   % Amplitud de la rampa [deg]
delay = 0.1;                  % Delay del actuador

PID_AP = pid(k_p_res,k_i_res,k_d_res);
Gs_theta = Utils.padeTF();
load('Gaps.mat');
Gaps = [Gap_u,Gap_alpha,Gap_theta,Gap_q];
GSASs = [GSAS_u,GSAS_alpha,GSAS_theta,GSAS_q];

fprintf("--------- Tiempos Caracteristicos AP para cada variable en el sistema PLANTA LIBRE vs SAS vs AUTOPILOT -----------\n");

for i = 1:length(varnames)
    fig = figure('Position', [100, 100, 1120, 420]);
    
    if(i == 3)
        Gap = feedback(Gaps(i)*PID_AP,Gs_theta);
    else
        Gap = Gaps(i);
    end
    % Transitorio
    subplot(1, 2, 1)
    hold on
    t = 0:0.05:tsimtr(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
   if(i == 3)
        [y_ap, t_out_ap, ~] = lsim(Gaps(i)*scales(i), u, t);
        [y_SAS,t_out_sas,~] = lsim(GSASs(i)*K_dl*scales(i),u,t);
        [y,t_out, ~] = lsim(plane_OL.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
    else
        [y_ap, t_out_ap, ~] = lsim(Gaps(i)*scales(i), u, t);
        [y_SAS,t_out_sas,~] = lsim(GSASs(i)*K_dl*scales(i),u,t);
        [y,t_out, ~] = lsim(plane_OL.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
    end
    plot(t_out, y, syms(3),'Color','k'); hold on
    plot(t_out_sas, y_SAS, syms(2),'Color','k'); hold on
    plot(t_out_ap, y_ap, syms(1),'Color','k'); hold on
    set(gca,'TickLabelInterpreter','latex')
    ax = gca;
    ax.FontSize = 13; 
    title('Transitorio','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    
    % Zoom en perturbacion
    
    if(i == 4 || i == 3 || i == 1)
        axes('Position',[.30 .67 .15 .15])
        box on; % put box around new pair of axes
        indexOfInterest = (t < 4) & (t > 0); % range of t near perturbation
        plot(t(indexOfInterest),y(indexOfInterest),syms(3),'Color','k'); hold on % plot on new axes
        plot(t(indexOfInterest),y_SAS(indexOfInterest),syms(2),'Color','k'); hold on 
        plot(t(indexOfInterest),y_ap(indexOfInterest),syms(1),'Color','k'); hold on 
        grid on;
        axis tight
    end
    if(i == 2)
        axes('Position',[.3 .2 .15 .15])
        box on; % put box around new pair of axes
        indexOfInterest = (t < 4) & (t > 0); % range of t near perturbation
        plot(t(indexOfInterest),y(indexOfInterest),syms(3),'Color','k'); hold on % plot on new axes
        plot(t(indexOfInterest),y_SAS(indexOfInterest),syms(2),'Color','k'); hold on 
        plot(t(indexOfInterest),y_ap(indexOfInterest),syms(1),'Color','k'); hold on 
        grid on;
        axis tight
    end
    
    % Estacionario
    subplot(1, 2, 2)
    t = 0:0.05:tsimst(i);
    u = max(0,min(amp/tramp*(t),amp));
    grid minor
    hold on
    box on
    if(i == 3)
        [y_ap, t_out_ap, ~] = lsim(Gaps(i)*scales(i), u, t);
        [y_SAS,t_out_sas,~] = lsim(GSASs(i)*K_dl*scales(i),u,t);
        [y,t_out, ~] = lsim(plane_OL.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
    else
        [y_ap, t_out_ap, ~] = lsim(Gaps(i)*scales(i), u, t);
        [y_SAS,t_out_sas,~] = lsim(GSASs(i)*K_dl*scales(i),u,t);
        [y,t_out, ~] = lsim(plane_OL.lon.G.(genvarname(Gnames(i))).(genvarname("Gfact"))*scales(i)*G_a*K_dl, u, t);
    end
    plot(t_out, y, syms(3),'Color','k'); hold on
    plot(t_out_sas, y_SAS, syms(2),'Color','k'); hold on
    plot(t_out_ap, y_ap, syms(1),'Color','k'); hold on
    set(gca,'TickLabelInterpreter','latex')
    ax = gca;
    ax.FontSize = 13; 
    title('Estacionario','FontSize', 15, 'Interpreter', 'latex')
    ylabel(Respnames_latex(i),'FontSize', 15, 'Interpreter', 'latex')
    xlabel('$t$ [s]','FontSize', 15, 'Interpreter', 'latex')
    ylim([min(y)-(abs(max(y)-min(y)))*0.1 max(y)+(abs(max(y)-min(y)))*0.1])
    legend("Planta libre","SAS","Autopiloto",'Interpreter','latex');

    if (save == false)
        %sgtitle(strcat('\underline{Respuesta rampa $\Delta\delta_e=', num2str(amp), angleunits, ', ', ' $t_{ramp}=', num2str(tramp), 's$}'), 'Interpreter', 'latex')
        sgtitle("")
    else
        sgtitle("")
        saveas(gcf,strcat("Graficas/AP/Ramp_", varnames(i),"SASvsNoSASvsAP"),'eps')
    end
    
    % Rise Time and Time Delay Calculation and printing in console
    
    [tr_p,td_p] = Utils.getTimesRD(y,t);
    [tr_ap,td_ap] = Utils.getTimesRD(y_ap,t);
    [tr_SAS,td_SAS] = Utils.getTimesRD(y_SAS,t);
    fprintf(" ---- " + varnames(i)+ " ----\n");
    fprintf("         PLANTA LIBRE | SAS | AP \n");
    fprintf("RISE TIME: " + num2str(tr_p) + " s   | " +num2str(tr_SAS) + " s  |" + num2str(tr_ap) + " s \n");
    fprintf("DELAY TIME: " + num2str(td_p) + " s   | " +num2str(td_SAS) + " s  |" + num2str(td_ap) + " s \n");
    
end


