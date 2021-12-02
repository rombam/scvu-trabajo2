classdef Utils
    
    % -- Contains useful functions (converting units, etc) --
    % Also contains ISA functions, scaling functions,TF manipulation utils
    % and sensor/actuator modelling TFs.
    
    methods (Static)    
        %% ---- Unit conversion ----
        
        function m_value = ft_to_m(ft_value)
            % Converts from ft to m 
            m_value = ft_value/3.2808;
        end
        
        function ms_value = kt_to_ms(kt_value)
            % Converts from knots to m/s
            ms_value = kt_value*0.5144;
        end
        
        function deg_value = rad_to_deg(rad_value)
            % Converts from radians to degrees
            deg_value = rad_value * 180/pi;
        end
        
        function rad_value = deg_to_rad(deg_value)
            % Converts from radians to degrees
            rad_value = deg_value * pi/180;
        end
        
        function kg_value = lb_to_kg(lb_value)
           % Converts from lbs to kgs
           kg_value = lb_value*0.453592;
        end
        
        function kgm2_value = slgft2_tokgm2(slgft2_value)
            % Converts from slugs*ft^2 to kg*m^2
            kgm2_value = slgft2_value*1.35581795;
        end

        function pa_value = lbft2_to_pa(lbft2_value)
            % Converts from lb/ft^2 to Pa
            pa_value = lbft2_value*47.880172;
        end
        
        %% ---- Atmosphere ----
        
        function [T, P, rho] = ISA(h)
            % ISA model. h: meters.
            R = 287.04;    % [m^2/(K*s^2)]
            P_0 = 101325;  % [Pa]
            T_0 = 288.15;  % [K]
            g_0 = 9.80665; % [m/s^2]
            
            if h < 11000
                T = T_0 - 6.5*h/1000;
                P = P_0*(1 - 0.0065*h/T_0)^5.2561;
            else
                T = T_0 - 6.5*11000/1000;
                P_11 = P_0*(1 - 0.0065*11000/T_0)^5.2561;
                P = P_11*exp(-g_0/(R*T)*(h-11000));
            end
            rho = P/(R*T);           
        end
        
        %% ---- Get dimension scales ----
        
        function scales = get_scales(plane)
            % Returns characteristic values
            [~, ~, rho] = Utils.ISA(plane.FC.hs);
            % Longitudinal
            scales.lon.length = plane.c/2;
            scales.lon.area = plane.Sw;
            scales.lon.speed = plane.FC.us;
            scales.lon.mass = plane.c*plane.Sw*rho/2;
            scales.lon.time = plane.c/(2*plane.FC.us);
            scales.lon.press = 0.5*rho*plane.FC.us^2;
            scales.lon.force = 0.5*rho*plane.FC.us^2*plane.Sw;
            scales.lon.moment = 0.5*rho*plane.FC.us^2*plane.Sw*plane.c;
            scales.lon.inertia = rho*plane.Sw*(plane.c/2)^3;
            % Lateral-directional
            scales.lat.length = plane.b/2;
            scales.lat.area = plane.Sw;
            scales.lat.speed = plane.FC.us;
            scales.lat.mass = plane.b*plane.Sw*rho/2;
            scales.lat.time = plane.b/(2*plane.FC.us);
            scales.lat.press = 0.5*rho*plane.FC.us^2;
            scales.lat.force = 0.5*rho*plane.FC.us^2*plane.Sw;
            scales.lat.moment = 0.5*rho*plane.FC.us^2*plane.Sw*plane.b;
            scales.lat.inertia = rho*plane.Sw*(plane.b/2)^3;            
        end
        
        %% -- Transfer function utils --

        function TF_from_sym = syms2tf(G)
            % Returns a transfer function from a symbolic expression
            % Crystal Nassouri 2009
            [symNum, symDen] = numden(G);   %Get num and den of Symbolic TF
            num = sym2poly(symNum);         %Convert Symbolic num to polynomial
            den = sym2poly(symDen);         %Convert Symbolic den to polynomial
            TF_from_sym = tf(num, den);
        end
        
        function sym_from_TF = tf2syms(G)
            % Returns a symbolic expression from a transfer function
            syms s
            [num, den] = tfdata(G,'v');
            sym_from_TF = poly2sym(num,s)/poly2sym(den,s);
        end
        
        function Gf = factorize_G(G)
            % Returns a factorized version of a transfer function and the
            % associated characteristic times, natural frequencies and
            % damping factors
            
            syms s
            % -- Poles (can be real or complex) --
            poles = pole(G);
            count = 1;
            tau_p = [];
            wn_p = [];
            xi_p = [];
            null_p = 0;
            while count <= length(poles)
                if imag(poles(count)) == 0
                    if poles(count) ~= 0
                        tau_p(end+1) = -1/poles(count);
                        count = count + 1;
                    else
                        null_p = null_p + 1;
                        count = count + 1;
                    end
                else
                    lambda_p = poles(count);
                    n_p = real(lambda_p);
                    w_p = abs(imag(lambda_p));
                    wn_p(end+1) = sqrt(n_p^2 + w_p^2);
                    xi_p(end+1) = -n_p/sqrt(n_p^2 + w_p^2);
                    count = count + 2;
                end
            end
            Denfact = 1;
            for i = 1:length(tau_p)
                Denfact = Denfact * (tau_p(i)*s + 1);
            end
            for i = 1:length(wn_p)
                Denfact = Denfact * ((s/wn_p(i))^2 + 2*xi_p(i)*(s/wn_p(i)) + 1);
            end   
            for i = 1:null_p
                Denfact = Denfact * s;
            end
            
            % -- Zeros (can be real or complex) --
            zeros  = zero(G);
            count  = 1;
            tau_z  = [];
            wn_z   = [];
            xi_z   = [];
            null_z = 0;
            while count <= length(zeros)
                if imag(zeros(count)) == 0
                    if zeros(count) ~= 0
                        tau_z(end+1) = -1/zeros(count);
                        count = count + 1;
                    else
                        null_z = null_z + 1;
                        count = count + 1;
                    end
                else
                    lambda_z = zeros(count);
                    n_z = real(lambda_z);
                    w_z = abs(imag(lambda_z));
                    wn_z(end+1) = sqrt(n_z^2 + w_z^2);
                    xi_z(end+1) = -n_z/sqrt(n_z^2 + w_z^2);
                    count = count + 2;
                end
            end

            % -- Static gain -- (removing s associated with null zeros and poles)
            K = vpa(limit(Utils.tf2syms(G)*s^null_p/s^null_z));
            
            Numfact = 1;
            for i = 1:length(tau_z)
                Numfact = Numfact * (tau_z(i)*s + 1);
            end
            for i = 1:length(wn_z)
                Numfact = Numfact * ((s/wn_z(i))^2 + 2*xi_z(i)*(s/wn_z(i)) + 1);
            end
            for i = 1:null_z
                Numfact = Numfact * s;
            end
            
            % -- Transfer function --
            Gf.Gfact_sym    = K*Numfact/Denfact;
            Gf.Gfact        = Utils.syms2tf(Gf.Gfact_sym);
            Gf.K            = double(K);
            Gf.zeros.tau    = tau_z;
            Gf.zeros.wn     = wn_z;
            Gf.zeros.xi     = xi_z;
            Gf.poles.tau    = tau_p;
            Gf.poles.wn     = wn_p;
            Gf.poles.xi     = xi_p;            
            
        end

        % ----- Sensor and actuator modelling -----
        function pade_TF = padeTF(lag,order)
            % Returns a Padé approximation given a lag and order
            % Used for sensor modelling
            if ~exist('order','var')
                % If Padé order is not specified, default to 3
                order = 1;
            end
            s = tf('s');
            baseTF = exp(-lag*s);    
            pade_TF = pade(baseTF,order);
        end
        
        function forlag_TF = forlagTF(wb)
            % Returns a first order lag transfer function
            % Used for actuator modelling
            s = tf('s');   
            forlag_TF = 1/(s/wb + 1);

        end
        
        function K_dl =getKdl(G)
            
            K_dl = 1/G.K;
        end
        
        function [t_r,t_d] = getTimesRD(y,t,amp)
            
        if ~exist('amp','var')
           amp = y(end);
        end
        dy = [];
        j = 0;
        t1=0;t2=0;
        [ymax,imax] = max(y);
        if amp + 0.01 < 0
            [ymax,imax] = min(y);
        end
        for i = 1:imax
            %Time delay calc
            if i <= imax
                dy = [dy,(y(i+1) - y(i))/(t(i + 1) - t(i))];
            end     
            % Rise time calc
            if abs(y(i)) >= abs(amp*0.1)
                if j == 0
                   t1 = t(i);
                   j = 1;
                end
                if abs(y(i)) <= abs(amp*0.9)
                    t2 = t(i);
                end
            end
            t_r = t2 - t1;

        end
        [dy_max,i_td] = max(dy);
        
        t_d = t(i_td) - (y(i_td)/dy_max);
        end
    end
end

