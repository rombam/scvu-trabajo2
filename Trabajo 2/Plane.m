classdef Plane
    
    % ------- Plane class --------
    % Contains geometry, mass, stability derivatives and TFs of a given plane
    % Transforms dimensional quantities to non dimensional and projects
    % stability derivatives on body axes. Calculates TFs from derivatives
    % and dynamic equations.
    
    properties
        
        % -- Name [str] -- 
        model_name
        
        % -- Geometric [float] --
        Sw        {mustBeNumeric}
        c         {mustBeNumeric}
        b         {mustBeNumeric}
        
        % -- Mass [float] --
        m         {mustBeNumeric}
        Ixb       {mustBeNumeric}
        Iyb       {mustBeNumeric}
        Izb       {mustBeNumeric}
        Ixzb      {mustBeNumeric}
        
        % -- Model [float] --
        model

        % -- Flight condition [struct] -- 
        FC

        % -- Channels --
        % Each one contains necessary non-dimensional magnitudes,
        % stability derivatives and transfer functions
        lon
        lat
        
            % - Longitudinal channel -
            %   Body-axes long. stability derivatives
            % C_Xu, C_Xalpha, C_Xdeltae
            % C_Zu, C_Zalpha, C_Zalphadot, C_Zq, C_Zdeltae
            % C_mu, C_malpha, C_malphadot, C_mq, C_mTu,
            % C_mdeltae
            %   Non dimensional time, mass, inertia
            % t_lon, mu_lon, Iyb_nd
            %   Transfer functions
            % G
        
            % - Lateral-directional channel -
            %   Body-axes lat.-dir. stability derivatives
            % C_Ybeta, C_Yp, C_Yr, C_Ydeltaa, C_Ydeltar
            % C_lbeta, C_lp, C_lr, C_ldeltaa, C_ldeltar
            % C_nbeta, C_np, C_nr, C_nTb, C_ndeltaa, 
            % C_ndeltar
            %   Non dimensional time, mass, inertia
            % t_lat, mu_lat, Ixb_nd, Izb_nd, Ixzb_nd
            %   Transfer functions
            % G        
    end
    
    methods
        function obj = Plane(model)
            % -- Load flight condition and derivatives into class --
            % Includes non dimensional characteristics
            scales = Utils.get_scales(model);
            
            % - Geometric and mass properties -
            obj.model_name    = model.name;
            obj.Sw            = model.Sw;
            obj.c             = model.c;
            obj.b             = model.b;
            obj.m             = model.m;
            obj.Ixb           = model.Ixb;
            obj.Iyb           = model.Iyb;
            obj.Izb           = model.Izb;
            obj.Ixzb          = model.Ixzb;
            
            % - Non dimensional properties -
            obj.lon.t_lon     = scales.lon.time;
            obj.lat.t_lat     = scales.lat.time;
            obj.lon.mu_lon    = model.m/scales.lon.mass;
            obj.lat.mu_lat    = model.m/scales.lat.mass;
            obj.lat.Ixb_nd    = model.Ixb/scales.lat.inertia;
            obj.lon.Iyb_nd    = model.Iyb/scales.lon.inertia;
            obj.lat.Izb_nd    = model.Izb/scales.lat.inertia;
            obj.lat.Ixzb_nd   = model.Ixzb/scales.lat.inertia;
            
            % - Flight condition -
            obj.FC     = model.FC;
            obj.FC.CXs = model.FC.CTs*cos(model.FC.epss) + model.FC.CLs*model.FC.alphabs - model.FC.CDs;
            obj.FC.CZs = -model.FC.CTs*sin(model.FC.epss) - model.FC.CLs - model.FC.CDs*model.FC.alphabs;
            obj.FC.Cms = 0;
            
            % - Flight derivatives. By default, dimensionless -
            % Store raw inputs
            obj.model.CL = model.CL;
            obj.model.CD = model.CD;
            obj.model.CT = model.CT;
            obj.model.CY = model.CY;
            obj.model.Cl = model.Cl;
            obj.model.Cm = model.Cm;
            obj.model.Cn = model.Cn;
            
            % Longitudinal
            obj.lon.CX.u = obj.model.CT.u*cos(obj.FC.epss) - obj.model.CD.u;
            obj.lon.CX.alpha = obj.model.CT.alpha*cos(obj.FC.epss) + obj.FC.CLs - obj.model.CD.alpha;
            obj.lon.CX.deltae = -obj.model.CD.deltae;
            obj.lon.CX.T = obj.FC.CTs*cos(obj.FC.epss);
            obj.lon.CZ.u = -obj.model.CT.u*sin(obj.FC.epss) - obj.model.CL.u;
            obj.lon.CZ.alpha = -obj.model.CT.alpha*sin(obj.FC.epss) -obj.model.CL.alpha - obj.FC.CDs;
            obj.lon.CZ.alphadot = -obj.model.CL.alphadot;
            obj.lon.CZ.q = -obj.model.CL.q;
            obj.lon.CZ.deltae = -obj.model.CL.deltae;
            obj.lon.CZ.T = -obj.FC.CTs*sin(obj.FC.epss);
            obj.lon.Cm = obj.model.Cm;
            obj.lon.Cm.u = obj.model.Cm.u + obj.model.Cm.Tu;
            
            % Lateral-directional
            obj.lat.CY = model.CY;
            obj.lat.Cl = model.Cl;
            obj.lat.Cn = model.Cn;
            
            % - Transfer functions -
            % Longitudinal
            obj.lon.G = Dynamics.LongDyn(obj);
            
            % Lateral
            obj.lat.G = Dynamics.LatDyn(obj);
        end      
        
        function obj = recalc(obj)
            % -- Recalculate body axes derivatives and TFs --
            % Useful if stability derivatives have been manually modified
            
            % Longitudinal channel
            obj.lon.CX.u = obj.model.CT.u*cos(obj.FC.epss)/cos(obj.FC.alphabs) - obj.model.CD.u ...
                       + obj.model.CL.u*obj.FC.alphabs;
            obj.lon.CX.alpha = obj.FC.CLs - obj.model.CD.alpha;
            obj.lon.CX.deltae = -obj.model.CD.deltae;
            obj.lon.CX.T = obj.FC.CTs*cos(obj.FC.epss);
            obj.lon.CZ.u = -obj.model.CT.u*sin(obj.FC.epss)/cos(obj.FC.alphabs) - obj.model.CL.u ...
                       - obj.model.CD.u*obj.FC.alphabs;
            obj.lon.CZ.alpha = -obj.model.CL.alpha - obj.FC.CDs;
            obj.lon.CZ.alphadot = -obj.model.CL.alphadot;
            obj.lon.CZ.q = -obj.model.CL.q;
            obj.lon.CZ.deltae = -obj.model.CL.deltae;
            obj.lon.CZ.T = -obj.FC.CTs*sin(obj.FC.epss)/cos(obj.FC.alphabs);
            obj.lon.Cm = obj.model.Cm;
            
            obj.lon.G = Dynamics.LongDyn(obj);

            % Lateral-directional channel
            obj.lat.CY = obj.model.CY;
            obj.lat.Cl = obj.model.Cl;
            obj.lat.Cn = obj.model.Cn;

            % - Transfer functions -
            % Longitudinal
            obj.lon.G = Dynamics.LongDyn(obj);
            
            % Lateral
            obj.lat.G = Dynamics.LatDyn(obj);
        end  
       

    end
end

