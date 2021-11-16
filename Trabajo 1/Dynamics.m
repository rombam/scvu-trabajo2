classdef Dynamics
    % -- Generates longitudinal and lateral-directional TFs for a plane --
 
    methods (Static)    
        %% ------- LONGITUDINAL CHANNEL --------
        function G = LongDyn(plane)
            
            % -- Calculate all longitudinal transfer functions --
            % Assumes control by elevator deflection (deltae) and thrust
            % (deltaT).
            
            % - Symbolic variables -
            clear s
            syms s
    
            % - State matrix -
            Alon = sym(zeros(3));
            Alon(1,1) = 2*plane.lon.mu_lon*plane.lon.t_lon*s - plane.lon.CX.u;
            Alon(1,2) = -plane.lon.CX.alpha;
            Alon(1,3) = -plane.FC.CZs;
            Alon(2,1) = -(plane.lon.CZ.u + 2*plane.FC.CZs);
            Alon(2,2) = (2*plane.lon.mu_lon - plane.lon.CZ.alphadot)*plane.lon.t_lon*s...
                        - plane.lon.CZ.alpha;
            Alon(2,3) = -(2*plane.lon.mu_lon + plane.lon.CZ.q)*plane.lon.t_lon*s;
            Alon(3,1) = -plane.lon.Cm.u;
            Alon(3,2) = -plane.lon.Cm.alphadot*plane.lon.t_lon*s - plane.lon.Cm.alpha;
            Alon(3,3) = plane.lon.Iyb_nd*plane.lon.t_lon^2*s^2 - plane.lon.Cm.q*plane.lon.t_lon*s;
    
            % - Control matrix -
            Blon = sym(zeros(3,2));
            Blon(1,1) = plane.lon.CX.deltae;
            Blon(1,2) = plane.lon.CX.T;
            Blon(2,1) = plane.lon.CZ.deltae;
            Blon(2,2) = plane.lon.CZ.T;
            Blon(3,1) = plane.lon.Cm.deltae;
            Blon(3,2) = plane.lon.Cm.T;
    
            % - Characteristic equation -
            Det = vpa(det(Alon,'Algorithm','minor-expansion'));
    
            % -- TRANSFER FUNCTIONS --
            names = ["udeltae" "uT"; "alphadeltae" "alphaT"; "thetadeltae" "thetaT"; "qdeltae" "qdeltaT"];
            [x, u] = size(names);
            for i = 1:u
                for j = 1:(x-1)
                Asub = Alon;
                Asub(:,j) = Blon(:,i);
                Num = vpa(det(Asub,'Algorithm','minor-expansion'));
        
                G.(genvarname(strcat('G',names(j,i)))) = Utils.factorize_G(Utils.syms2tf(Num/Det));
                end
                Num = Num*s;
                % Pitch rate transfer function (qhat = Dtheta/Dt)
                G.(genvarname(strcat('G',names(j+1,i)))) = Utils.factorize_G(Utils.syms2tf(Num/Det)*plane.lon.t_lon);
            end
        end

        %% ------- LATERAL-DIRECTIONAL CHANNEL --------
        function G = LatDyn(plane)
            
            % -- Calculate all lateral-directional transfer functions --
            % Assumes control by aileron deflection (deltaa) and rudder
            % deflection (deltar).
            
            % - Symbolic variables -
            clear s
            syms s
    
            % - State matrix -
            Alat = sym(zeros(3));
            Alat(1,1) = 2*plane.lat.mu_lat*plane.lat.t_lat*s - plane.lat.CY.beta;
            Alat(1,2) = plane.FC.CZs - plane.lat.CY.p*plane.lat.t_lat*s;
            Alat(1,3) = 2*plane.lat.mu_lat - plane.lat.CY.r;
            Alat(2,1) = -plane.lat.Cl.beta;
            Alat(2,2) = plane.lat.Ixb_nd*plane.lat.t_lat^2*s^2 - plane.lat.Cl.p*plane.lat.t_lat*s;
            Alat(2,3) = -plane.lat.Ixzb_nd*plane.lat.t_lat*s - plane.lat.Cl.r;
            Alat(3,1) = -plane.lat.Cn.beta;
            Alat(3,2) = -plane.lat.Ixzb_nd*plane.lat.t_lat^2*s^2 - plane.lat.Cn.p*plane.lat.t_lat*s;
            Alat(3,3) = plane.lat.Izb_nd*plane.lat.t_lat*s - plane.lat.Cn.r;
    
            % - Control matrix -
            Blat = sym(zeros(3,2));
            Blat(1,1) = plane.lat.CY.deltaa;    % Generally, 0
            Blat(1,2) = plane.lat.CY.deltar;
            Blat(2,1) = plane.lat.Cl.deltaa;
            Blat(2,2) = plane.lat.Cl.deltar;
            Blat(3,1) = plane.lat.Cn.deltaa;
            Blat(3,2) = plane.lat.Cn.deltar;
    
            % - Characteristic equation -
            Det = vpa(det(Alat,'Algorithm','minor-expansion'));
            %Det = vpa(det(Alat));
    
            % -- TRANSFER FUNCTIONS --
            names = ["betadeltaa" "betadeltar"; "phideltaa" "phideltar"; "rdeltaa" "rdeltar"; "pdeltaa" "pdeltar"; "psideltaa" "psideltar"];
            [x, u] = size(names);
            for i = 1:u
                for j = 1:(x-2)
                Asub = Alat;
                Asub(:,j) = Blat(:,i);
                Num(j,i) = vpa(det(Asub,'Algorithm','minor-expansion'));
        
                G.(genvarname(strcat('G',names(j,i)))) = Utils.factorize_G(Utils.syms2tf(Num(j,i)/Det));
                end
            end
            % p transfer functions
            Num(4,1) = Num(2,1)*s*plane.lat.t_lat;
            Num(4,2) = Num(2,2)*s*plane.lat.t_lat;
            G.(genvarname(strcat('G',names(4,1)))) = Utils.factorize_G(Utils.syms2tf(Num(4,1)/Det));
            G.(genvarname(strcat('G',names(4,2)))) = Utils.factorize_G(Utils.syms2tf(Num(4,2)/Det));

            % Psi transfer functions
            Num(5,1) = Num(3,1)/s/plane.lat.t_lat;
            Num(5,2) = Num(3,2)/s/plane.lat.t_lat;
            G.(genvarname(strcat('G',names(5,1)))) = Utils.factorize_G(Utils.syms2tf(Num(5,1)/Det));
            G.(genvarname(strcat('G',names(5,2)))) = Utils.factorize_G(Utils.syms2tf(Num(5,2)/Det));
        end

    end
end

