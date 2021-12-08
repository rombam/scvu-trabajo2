%% Simulink Validation
close all, clear all;

plane_OL = Plane(GlobalHawk);

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


%% Obtencion de parametros para simulink 

delay = 0.1;
G_a_deltae = Utils.padeTF(delay); [G_a_deltae_num,G_a_deltae_den ] = tfdata(G_a_deltae);
K_dl = 1/plane_OL.lon.G.Gthetadeltae.K_dim;
Gs_alpha = Utils.padeTF(); [Gs_alpha_num,Gs_alpha_den] = tfdata(Gs_alpha);
Gs_q = Utils.padeTF(); [Gs_q_num,Gs_q_den] = tfdata(Gs_q);
Gs_theta = Utils.padeTF(); [Gs_theta_num,Gs_theta_den] = tfdata(Gs_theta);

Gs_alpha_num = cell2mat(Gs_alpha_num); Gs_alpha_den = cell2mat(Gs_alpha_den);
Gs_q_num = cell2mat(Gs_q_num); Gs_q_den = cell2mat(Gs_q_den);
Gs_theta_num = cell2mat(Gs_theta_num); Gs_theta_den = cell2mat(Gs_theta_den);
G_a_deltae_num = cell2mat(G_a_deltae_num); G_a_deltae_den = cell2mat(G_a_deltae_den);

K_deltae_q = 0;
K_deltae_alpha = 0;

[Gp_u_deltae_num,Gp_u_deltae_den]=tfdata(zpk(zero(plane_OL.lon.G.(genvarname(Gnames(1))).(genvarname("Gfact"))*scales(1)),pole(plane_OL.lon.G.(genvarname(Gnames(1))).(genvarname("Gfact"))*scales(1)),plane_OL.lon.G.(genvarname(Gnames(1))).K_dim));
[Gp_alpha_deltae_num,Gp_alpha_deltae_den]=tfdata(zpk(zero(plane_OL.lon.G.(genvarname(Gnames(2))).(genvarname("Gfact"))*scales(2)),pole(plane_OL.lon.G.(genvarname(Gnames(2))).(genvarname("Gfact"))*scales(2)),plane_OL.lon.G.(genvarname(Gnames(2))).K_dim));
[Gp_theta_deltae_num,Gp_theta_deltae_den]=tfdata(zpk(zero(plane_OL.lon.G.(genvarname(Gnames(3))).(genvarname("Gfact"))*scales(3)),pole(plane_OL.lon.G.(genvarname(Gnames(3))).(genvarname("Gfact"))*scales(3)),plane_OL.lon.G.(genvarname(Gnames(3))).K_dim));
[Gp_q_deltae_num,Gp_q_deltae_den]=tfdata(zpk(zero(plane_OL.lon.G.(genvarname(Gnames(4))).(genvarname("Gfact"))*scales(4)),pole(plane_OL.lon.G.(genvarname(Gnames(4))).(genvarname("Gfact"))*scales(4)),plane_OL.lon.G.(genvarname(Gnames(4))).K_dim));

Gp_u_deltae_num = cell2mat(Gp_u_deltae_num); Gp_u_deltae_den = cell2mat(Gp_u_deltae_den);
Gp_alpha_deltae_num = cell2mat(Gp_alpha_deltae_num); Gp_alpha_deltae_den = cell2mat(Gp_alpha_deltae_den);
Gp_theta_deltae_num = cell2mat(Gp_theta_deltae_num); Gp_theta_deltae_den = cell2mat(Gp_theta_deltae_den);
Gp_q_deltae_num = cell2mat(Gp_q_deltae_num); Gp_q_deltae_den = cell2mat(Gp_q_deltae_den);
