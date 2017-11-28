function [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION] = ST(P_e,options,display)
% ST Steam power plants modelisation
% ST(P_e,options,display) compute the thermodynamics states for a Steam
% power plant (combustion, exchanger, cycle) turbine based on several 
% inputs (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [W]
% OPTIONS is a structure containing :
%   -options.nsout     [-] : Number of feed-heating 
%   -options.reheat    [-] : Number of reheating
%   -options.T_max     [°C] : Maximum steam temperature
%   -options.T_cond_out[°C] : Condenseur cold outlet temperature
%   -options.p3_hp     [bar] : Maximum pressure
%   -options.drumFlag  [K] : if =1 then drum if =0 => no drum. 
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax     [°C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.p_3       [-] : High pressure after last reheating
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (titre)
%   -options.T0        [°C] : Reference temperature
%   -options.TpinchSub [°C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [°C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[°C] : Temperature pinch at condenser 
%   -options.Tdrum     [°C] : drum temperature
%   -option.eta_SiC    [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_gen, Steam generator energy efficiency
%   -eta(6) : eta_gex, Steam generator exergy efficiency
%   -eta(7) : eta_combex, Combustion exergy efficiency
%   -eta(8) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(9) : eta_transex, Heat exchanger overall exergy efficiency
% Xmassflow is a vector with each feedheating massflow (respect to figure 
%           2.33, page 91 "Thermal Power Plants" English version).
%           Xmassflow(1) = mass flow at 6_1 etc...
% DATEN is a vector with : 
%   -daten(1) : perte_gen [W]
%   -daten(2) : perte_mec [W]
%   -daten(3) : perte_cond [W]
% DATEX is a vector with :
%   -datex(1) : perte_mec    [W]
%   -datex(2) : perte_totex  [W]
%   -datex(3) : perte_rotex  [W]
%   -datex(4) : perte_combex [W]
%   -datex(5) : perte_condex [W]
%   -datex(6) : perte_chemex [W]
%   -datex(7) : perte_transex[W]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , ...       , T_4;  [°C]
%        p_1       , p_2       , ...       , p_4;  [bar]
%        h_1       , h_2       , ...       , h_4;  [kJ/kg]
%        s_1       , s_2       , ...       , s_4;  [kJ/kg/K]
%        e_1       , e_2       , ...       , e_4;  [kJ/kg]
%        x_1       , x_2       , ...       , x_4;};[-]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_v, water massflow at 2 [kg/s]
%   -massflow(3) = m_c, combustible massflow [kg/s] 
%   -massflow(4) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
%   -combustion.fum  : is a vector of the exhaust gas composition :
%       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s]

if nargin<3
    display = 1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 250e6; % [W] Puissance énergétique de l'installation
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% n_sout limite a 20  %%
%% n_reheat limite a 5 %%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Description des etats
% etat 1  : entr�e pompe alimentaire        (sortie dernier resurchauffeur)
% etat 2  : entr�e chaudi�re                (sortie pompe alimentaire)
% etat 21 : vapeur sature chaudiere  
% etat 22 : vapeur surchaufee chaudiere
% etat 3  : sortie chaudiere avant resurchauffe
% etat 4  : sortie HP pour resurchauffe
% etat 5  : sortie chaudiere apres resurchauffe
% etat 6  : sortie BP                       (entree condenseur)
% etats sout_turb(1->nsout+nresur) : soutirages turbines (etats 6_I, 6_II...)
% etat 7  : sortie condenseur               (entree pompe) 
% etats sout_resur(1->nsout+nresur) : sortie resurchauffeurs R_I, R_II,...R_nsout
% etat 8  : sortie pompe apres condenseur   (entr�e resurchauffeur R0)
% etat 9  : sortie resurchauffeur R0        
% etats sout_princ(1->n_sout+nresur) : sortie principale resurchauffeurs R_I,
% R_II,...R_nsout (etats 9_I, 9_II...)






%% Input initialisation
%   -options.T0        [°C] : Reference temperature
if isfield(options,'T_0')           
    T_0 = options.T_0;    
else
    T_0 = 288.15;  % [éC] 
end
%   -options.nsout     [-] : Number of feed-heating
if isfield(options,'nsout')           
    data.nsout = options.nsout;    
else
    data.nsout = 5;  % [éC] 
end
%   -options.reheat    [-] : Number of reheating
if isfield(options,'reheat')           
    data.reheat = options.reheat;    
else
    data.reheat = 1;  % [éC] 
end
%   -options.T_max     [°C] : Maximum steam temperature
if isfield(options,'T_max')           
    data.T_max = options.T_max;    
else
    data.T_max = 525;  % [éC] 
end
%   -options.T_cond_out[°C] : Condenseur cold outlet temperature
if isfield(options,'T_cond_out')           
    data.T_cond_out = options.T_cond_out;    
else
    data.T_cond_out = 15;  % [éC] 
end
%   -options.p3_hp     [bar] : Maximum pressure
if isfield(options,'p3_hp')           
    data.p3_hp = options.p3_hp;    
else
    data.p3_hp = 15;  % [éC] 
end
%   -options.drumFlag  [K] : if =1 then drum if =0 => no drum. 
if isfield(options,'drumFlag')           
    data.drumFlag = options.drumFlag;    
else
    data.drumFlag = 1;  % [éC] 
end
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
if isfield(options,'eta_mec')           
    data.eta_mec = options.eta_mec;    
else
    data.eta_mec = 0.9;  % [éC] 
end
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax     [°C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
if isfield(options,'comb')           
    data.comb = options.comb;    
else
    data.comb.Tmax   = 1500;  % [éC] 
    data.comb.lambda =  0.6;  % [-] 
    data.comb.x      = 0.05;
    data.comb.y      =  1.2;
end
%   -options.p_3       [-] : High pressure after last reheating
if isfield(options,'p_3')           
    data.p_3 = options.p_3;    
else
    data.p_3 = 30;   
end
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (titre)
if isfield(options,'x4')           
    data.x4 = options.x4;    
else
    data.x4 = 0.5;   
end
%   -options.TpinchSub [°C] : Temperature pinch at the subcooler
if isfield(options,'TpinchSub')           
    data.TpinchSub = options.TpinchSub;    
else
    data.TpinchSub = 5;   
end
%   -options.TpinchEx  [°C] : Temperature pinch at a heat exchanger
if isfield(options,'TpinchEx')           
    data.TpinchEx = options.TpinchEx;    
else
    data.TpinchEx = 5;   
end
%   -options.TpinchCond[°C] : Temperature pinch at condenser 
if isfield(options,'TpinchCond')           
    data.TpinchCond = options.TpinchCond;    
else
    data.TpinchCond = 5;   
end
%   -options.Tdrum     [°C] : drum temperature
if isfield(options,'Tdrum')           
    data.Tdrum = options.Tdrum;    
else
    data.Tdrum = 20;   
end
%   -option.eta_SiC    [-] : Isotrenpic efficiency for compression
if isfield(options,'eta_SiC')           
    data.eta_SiC = options.eta_SiC;    
else
    data.eta_SiC = 0.9;   
end
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine
if isfield(options,'eta_SiT')           
    data.eta_SiT = options.eta_SiT;    
else
    data.eta_SiT(1) = 0.9;
    data.eta_SiT(2) = 0.9;
    data.eta_SiT(3) = 0.9;
end
%% Intermidiary Initialisation
% Reference State
    data.T0 =25;
    data.p_ref=1;
    data.h_ref=XSteam('h_pT',data.p_ref,data.T0);
    data.s_ref=XSteam('s_pT',data.p_ref,data.T0);
% Compression

    data.TurbLP_safe = 5;
    data.TurbHP_p_out = 50;
    data.TurbIP_p_out =  10;
    data.TurbHP_comp = data.p3_hp/data.TurbHP_p_out;
    data.TurbIP_comp = data.TurbHP_p_out/data.TurbIP_p_out;
    data.TurbLP_comp = 6;
    
    
    data.Sout_Pe_ratio = 2;
% Steam Generator
    data.SG_ploss = 0.1;
    
    data.v_eau = 1/1000; %(volume massique eau)
    
% Exergy
    function [ex] = exergy(h,h_ref,s,s_ref,T0)
        ex = (h-h_ref+T0*(s-s_ref));
    end
%% Simulation

%% INIT : 
% SORTIE CHAUDIERE AVANT (possible) RESURCHAUFFE ---> OK
    data.result(3).p  = data.p3_hp;
    data.result(3).T  = data.T_max;
    data.result(3).h  = XSteam('h_pT',data.result(3).p,data.result(3).T);
    data.result(3).s  = XSteam('s_pT',data.result(3).p,data.result(3).T);
    data.result(3).v  = XSteam('v_pT',data.result(3).p,data.result(3).T);
    data.result(3).x  = titre_corr(XSteam('x_ph',data.result(3).p,data.result(3).h),data.result(3).h,data.result(3).p);
    data.result(3).ex = exergy(data.result(3).h,data.h_ref,data.result(3).s,data.s_ref,data.T0);
    
% SORTIE CONDENSEUR ---> OK
    data.result(7).T  = T_cond_out;
    data.result(7).p  = XSteam('psat_T',data.result(7).T);
    data.result(7).h  = XSteam('hL_T',data.result(7).T);
    data.result(7).s  = XSteam('sL_T',data.result(7).T);
    data.result(7).x  = titre_corr(XSteam('x_ph',data.result(7).p,data.result(7).h),data.result(7).h,data.result(7).p);
    %VERIFIER QUE BIEN TITRE NUL -> liquide sature
    data.result(7).v  = XSteam('vL_p',data.result(7).p);
    data.result(7).ex = exergy(data.result(7).h, data.h_ref, data.result(7).s, data.s_ref, data.T0);
 
if data.reheat>0
    %% SORTIE CHAUDIERE APRES RESURCHAUFFE ---> OK
     data.result(5).T  = data.T_max;
     data.result(5).p  = data.p_3; %haute pression apres derniere resurchauffe
     data.result(5).h  = XSteam('h_pT', data.result(5).p, data.result(5).T);
     data.result(5).s  = XSteam('s_pT', data.result(5).p, data.result(5).T);
     data.result(5).x  = titre_corr(XSteam('x_ph', data.result(5).p, data.result(5).T),data.result(5).h,data.result(5).p);
     data.result(5).v  = XSteam('v_pT', data.result(5).p,data.result(5).T);
     data.result(5).ex = exergy(data.result(5).h, data.h_ref, data.result(5).s, data.s_ref, data.T0);

    %% SORTIE HP Turbine (pour resurchauffe) ---> OK
     data.result(4).p = data.result(5).p;
     s_4_is = data.result(3).s;
     h_4_is = XSteam('h_ps', data.result(4).p, s_4_is);
     data.result(4).h  = data.result(3).h + data.eta_SiT(1)*(h_4_is-data.result(3).h);
     data.result(4).s  = XSteam('s_ph', data.result(4).p, data.result(4).h);
     data.result(4).T  = XSteam('T_ph', data.result(4).p, data.result(4).h);
     data.result(4).x  = titre_corr(XSteam('x_ph',data.result(4).p,data.result(4).h),data.result(4).h,data.result(4).p);
     data.result(4).v  = XSteam('v_pT', data.result(4).p,data.result(4).T);
     data.result(4).ex = exergy(data.result(4).h, data.h_ref, data.result(4).s, data.s_ref, data.T0);
end

% 
%     %% IP Turbine ----> !!!!!!!!! PAS OK
%     data.result(32).p = data.result(311).p/data.TurbIP_comp;
%     s_32_is = data.result(311).s;
%     h_32_is = XSteam('h_ps', data.result(32).p, s_32_is);
%     data.result(32).h  = data.result(311).h + data.eta_SiT(2)*(h_32_is-data.result(311).h);
%     data.result(32).T  = XSteam('T_ph', data.result(32).p, data.result(32).h);
%     data.result(32).x  = XSteam('x_ph', data.result(32).p, data.result(32).h);
%     data.result(32).s  = XSteam('s_ph', data.result(32).p, data.result(32).h);
%     data.result(32).v  = XSteam('v_ph',data.result(32).p,data.result(32).h);
%     data.result(32).ex = exergy(data.result(32).h, data.h_ref, data.result(32).s, data.s_ref, data.T0);

%% LP Turbine : SORTIE ----> OK
% etat isentropique : hyp : temperature fin de detente idem condenseur
    data.result(6).T = data.T_cond_out;
    data.result(6).p = XSteam('psat_t', data.result(6).T);
    s_6_is = data.result(50).s;
    h_6_is = XSteam('h_ps', data.result(6).p, s_6_is);
% etat reel
    data.result(6).h  = data.result(5).h + data.eta_SiT(3)*(h_6_is-data.result(5).h);
    data.result(6).T  = XSteam('T_ph',data.result(6).p,data.result(6).h);
    data.result(6).x  = titre_corr(XSteam('x_ph', data.result(6).p, data.result(6).h),data.result(6).h,data.result(6).p);
    data.result(6).s  = XSteam('s_ph', data.result(6).p, data.result(6).h);
    data.result(6).v  = XSteam('v_ph',data.result(6).p,data.result(6).h);
    data.result(6).ex = exergy(data.result(6).h, data.h_ref, data.result(6).s, data.s_ref, data.T0);
    
%% Soutirages - etats soutire des turbines
data.Sout_turbine = zeros(data.nsout,5); % p,T,h,s,x
if data.reheat>0    
    delta_h_sout = (data.result(5).h-data.result(6).h)/(data.nsout+1+data.reheat);
end
for i=1:data.nsout 
    data.Sout_turbine(i,3) = data.result(5).h-delta_h_sout*i;
    h_4i_is                = data.result(5).h - (data.result(50).h-data.Sout_turbine(i,3))/data.eta_SiT(3);
    s_4i_is                = data.result(5).s;
    data.Sout_turbine(i,1) = XSteam('p_hs',h_4i_is,s_4i_is); %idem que p_6i_is
    data.Sout_turbine(i,2) = XSteam('T_ph',data.Sout_turbine(i,1),data.Sout_turbine(i,3));
    data.Sout_turbine(i,4) = XSteam('s_ph',data.Sout_turbine(i,1),data.Sout_turbine(i,3));
    data.Sout_turbine(i,5) = titre_corr(XSteam('x_ph',data.Sout_turbine(i,1),data.Sout_turbine(i,3)),data.Sout_turbine(i,3),data.Sout_turbine(i,1));
    
    if i==data.nsout && data.reheat>0 % dernier soutirage = etat avant resurchauffe
        data.Sout_turbine(i,1) = data.result(4).p;
        data.Sout_turbine(i,2) = data.result(4).T;
        data.Sout_turbine(i,3) = data.result(4).h;
        data.Sout_turbine(i,4) = data.result(4).s;
        data.Sout_turbine(i,5) = data.result(4).x;
    end
end

% Etat vapeur sature lors du soutirage (pour plot)
data.nsout_sat = zeros(data.nsout,5);
for i=1:data.nsout
    data.nsout_sat(i,1) = data.Sout_turbine(i,1);
    data.nsout_sat(i,2) = XSteam('Tsat_p',data.nsout_sat(i,1));
    data.nsout_sat(i,3) = XSteam('hV_p',data.nsout_sat(i,1));
    data.nsout_sat(i,4) = XSteam('sV_p',data.nsout_sat(i,1));
    data.nsout_sat(i,5) = titre_corr(XSteam('x_ph',data.nsout_sat(i,1),data.nsout_sat(i,3)),data.nsout_sat(i,3),data.nsout_sat(i,1));
end

 
%% SORTIE RECHAUFFEURS R_I, R_II,...R_nsout
% amene le fluide chaud de maniere isobare jusque etat de liquide sature
% liquide sature a la pression de soutirage
data.Sout_resur = zeros(data.nsout,5); % p,T,h,s,x 
for i=1:data.nsout
    data.Sout_resur(i,1) = data.Sout_turbine(i,1); %isobare
    data.Sout_resur(i,2) = XSteam('Tsat_p',data.Sout_resur(i,1));
    data.Sout_resur(i,3) = XSteam('hL_p',data.Sout_resur(i,1));
    data.Sout_resur(i,4) = XSteam('sL_p',data.Sout_resur(i,1));
    data.Sout_resur(i,5) = titre_corr(XSteam('x_ph',data.Sout_resur(i,1),data.Sout_resur(i,3)),data.Sout_resur(i,3),data.Sout_resur(i,1));
end

%% SORTIE POMPE d'extraction Pe, (ENTREE R0) ---> OK
data.result(8).p = 10; % = p_Pe A FIXER AU PREALABLE :
% imposer Pe telle que la pression en 9_nsout (ou avant le degazeur) soit suffisamment �lev�e pour
% qu'a la temp�rature de 9_nsout on soit toujours a l'�tat liquide et ne
% pas avoir de vapeur (donc prendre Pe> Psat(T_9_nsout)
h_8_is           = data.result(7).h + data.v_eau*(data.result(8).p-data.result(7).p)*(10^5)/(10^3);
data.result(8).h = data.result(7).h - data.eta_SiC*(data.result(7).h-h_8_is);
data.result(8).T = data.result(7).T + ((data.result(8).h-data.result(7).h)-0.09*(data.result(8).p-data.result(7).p))/4.18;
data.result(8).s = XSteam('s_pT',data.result(8).p,data.result(8).T);
data.result(8).x =  titre_corr(XSteam('x_ph', data.result(8).p, data.result(8).h),data.result(8).h,data.result(8).p);
data.result(8).v  = XSteam('v_ph',data.result(8).p,data.result(8).h);
%VERIFIER QUE BIEN UN TITRE INEXISTANT    

%% SORTIE principale ECHANGEURS 
% echangeurs non parfaits : tpinch
data.Sout_princ = zeros(data.nsout,5); % p,T,h,s,x 
for i=1:data.nsout
    data.Sout_princ(i,1) = data.result(8).p;
    data.Sout_princ(i,2) = data.Sout_resur(i,2) - data.TpinchEx;
    data.Sout_princ(i,3) = XSteam('h_pT',data.Sout_princ(i,1),data.Sout_princ(i,2));
    data.Sout_princ(i,4) = XSteam('s_pT',data.Sout_princ(i,1),data.Sout_princ(i,2));
    data.Sout_princ(i,5) = titre_corr(XSteam('x_ph',data.Sout_princ(i,1),data.Sout_princ(i,3)),data.Sout_princ(i,3),data.Sout_princ(i,1));
end

%% Steam Generator

% Etat 2
    data.result(2).p = data.result(4).p;
    h20_s = data.result(1).h + data.v_eau*((data.result(2).p-data.result(1).p)*10^5); %delta_h = w_m car isentropique adiabatique, q=0, w_f=0)
    data.result(2).h = data.result(1).h - data.eta_SiC*(data.result(1).h-h2_s);
    data.result(2).T = data.result(4).T ; %difference de temperature negligeable, decommenter ligne suivante pour le prouver...
                        %...+(1/1000*10e5/10e3-0.006*(data.result(20).p-data.result(1).p))/4.18;
                        %%pour etre plus precis
    data.result(2).s = XSteam('s_ph',data.result(2).p,data.result(2).h); %OLD mais a changer..
    data.result(2).v = XSteam('v_pT',data.result(2).p,data.result(2).T);
    data.result(2).x = 0/0;
    data.result(2).ex = exergy(data.result(2).h, data.h_ref,data.result(2).s,data.s_ref,data.T0);
    

% Etat 2'
    data.result(21).p = data.result(3).p;
     %data.result(21).p = data.result(3).p*(1+data.SG_ploss/2);
    data.result(21).T = XSteam('Tsat_p',data.result(21).p);
    data.result(21).h = XSteam('hL_p', data.result(21).p);
    data.result(21).s = XSteam('sL_p',data.result(21).p);
    data.result(21).v = XSteam('vL_p',data.result(21).p);
    data.result(21).x = 0;
    data.result(21).ex = exergy(data.result(21).h, data.h_ref, data.result(21).s, data.s_ref, data.T0);

% Etat 2''    
    data.result(22).p = data.result(3).p;
    data.result(22).T = XSteam('Tsat_p',data.result(21).p);
    data.result(22).h = XSteam('hV_p', data.result(21).p);
    data.result(22).s = XSteam('sV_p',data.result(21).p);
    data.result(22).v = XSteam('vV_p',data.result(21).p);
    data.result(22).x = 1;
    data.result(22).ex = exergy(data.result(22).h, data.h_ref, data.result(22).s, data.s_ref, data.T0);
    
    
%% Reheating
% if data.reheat>0
%     for i=1:data.reheat
%         data.result(50+i).T = data.result(30).T;
%         data.result(50+i).p = data.TurbHP_p_out;
%         data.result(50+i).h = XSteam('h_pT', data.result(50+i).p, data.result(50+i).T);
%         data.result(50+i).s = XSteam('s_pT', data.result(50+i).p, data.result(50+i).T);
%         data.result(50+i).x = XSteam('x_ph', data.result(50+i).p, data.result(50+i).h);
%         data.result(50+i).ex = exergy(data.result(50+i).h, data.h_ref, data.result(50+i).s, data.s_ref, data.T0);
%     end
% end

% %% Calcule du nombre de desurchauffes :
% count_desurch = 0;
% for i=1:n_sout
%     while (data.result(60+n_sout).T < data.result(60+n_sout-i)) %tant que T_HP < T_soutire
%         count_desurch=count_desurch+1;
%     end
% end


% %% Fraction soutirages (avant degazeur)
% %resoudre systeme Ax = b avec X les soutirages
% A_bf_deg = zeros(n_sout, n_sout); %A_before_degazeur
% b_bf = zeros(n_sout,1);
% %remplissage matrice A_bf
% for i=1:n_sout
%    for j=1:n_sout
%        if i==1
%            A_bf_deg(i,j) = -(data.result(90).h-data.result(80).h)+(data.result(100).h-data.result(120+i).h); %juste avant le condenseur
%        elseif i<j
%            A_bf_deg(i,j) = (data.result(70+i).h-data.result(120+i+1).h)-(data.result(90+i).h-data.result(90+i-1).h);
%        elseif i==j
%            A_bf_deg(i,j) = (data.result(70+i).h-data.result(60+i).h)-(data.result(90+i).h-data.result(90+i-1).h);
%        else
%            A_bf_deg(i,j) = -(data.result(90+i).h-data.result(90+i-1)).h;
%        end
%    end
% end
% %remplissage termes independants b_bf
% for i=1:n_sout
%     if i==1
%         b_bf(i,1) = data.result(90).h-data.result(80).h;
%     end
%     b_bf(i,1) = data.result(90+i).h-data.result(90+i-1).h;
% end
% 
% %% Fraction soutirages (apres degazeur)
% %resoudre systeme Ax = b avec X les soutirages
% A_af_deg = zeros(n_sout, n_sout); %A_after_degazeur
% b_af = zeros(n_sout,1);
% %remplissage matrice A_bf
% for i=1:n_sout
%    for j=1:n_sout
%        if i==n_sout
%            A_af_deg(i,j) = -(data.result(1).h-data.result(90+i-1).h); %juste avant la pompe chaudiere
%        elseif i<j
%            A_af_deg(i,j) =  (data.result(70+i).h-data.result(120+i+1).h)-(data.result(90+i).h-data.result(90+i-1).h);
%        elseif i==j
%            A_af_deg(i,j) =  (data.result(70+i).h-data.result(60+i).h)-(data.result(90+i).h-data.result(90+i-1).h);
%        else
%            A_af_deg(i,j) = -(data.result(90+i).h-data.result(90+i-1).h);
%        end
%    end
% end
% %remplissage termes independants b_bf
% for i=1:n_sout
%     if i==1
%         b_af(i,1) = data.result(90).h-data.result(80).h;
%     end
%     b_af(i,1) = data.result(90+i).h-data.result(90+i-1).h;
% end
   
%% Titre XSteam (correction)
function [x_corr] = titre_corr(x_current,h_current,p)
    if x_current == 0 && h_current<XSteam('hL_p',p)
        x_corr = nan; 
    elseif x_current == 1 && h_current>XSteam('hV_p',p)
        x_corr = nan;
    else 
        x_corr = x_current;
    end    
end


 
%% Display
    if display == true
        figure
        T=linspace(0,375,1000);
    
        for i=1:length(T)
            S1(i)=XSteam('sL_T',T(i));
            S2(i)=XSteam('sV_T',T(i));
        end
    
        S=[S1,S2];
        %S=S1;
        T=[T,T];
    
        
        xlabel('S [kJ/kg K]')
        ylabel('T [degC]')
        subplot(2,2,1);
        plot(S,T); hold on
        %hold on;
        title('Rankine-Hirn Cycle')
        %xlabel('s[kJ/kgK]')
        %ylabel('T[K]')
        plot(data.result(1).s,data.result(1).T,'.','MarkerSize',15)
        text(data.result(1).s,data.result(1).T,'1')
        plot(data.result(20).s,data.result(20).T,'.','MarkerSize',15)
        text(data.result(20).s,data.result(20).T,'20')
        plot(data.result(21).s,data.result(21).T,'.','MarkerSize',15)
        text(data.result(21).s,data.result(21).T,'21')
        plot(data.result(22).s,data.result(22).T,'.','MarkerSize',15)
        text(data.result(22).s,data.result(22).T,'22')
        plot(data.result(3).s,data.result(3).T,'.','MarkerSize',15)
        text(data.result(3).s,data.result(3).T,'3')
        plot(data.result(31).s,data.result(31).T,'.','MarkerSize',15)
        text(data.result(31).s,data.result(31).T,'31')
        plot(data.result(311).s,data.result(311).T,'.','MarkerSize',15)
        text(data.result(311).s,data.result(311).T,'311')
        plot(data.result(32).s,data.result(32).T,'.','MarkerSize',15)
        text(data.result(32).s,data.result(32).T,'32')
        plot(data.result(4).s,data.result(4).T,'.','MarkerSize',15)
        text(data.result(4).s,data.result(4).T,'4')

        plot([data.result(1).s data.result(20).s], [data.result(1).T data.result(20).T])
        %splot([data.result(20).s data.result(3).s], [data.result(20).T data.result(3).T] )
        plot([data.result(21).s data.result(22).s], [data.result(21).T data.result(22).T] )
        plot([data.result(22).s data.result(3).s], [data.result(22).T data.result(3).T] )
        plot([data.result(3).s data.result(31).s], [data.result(3).T data.result(31).T] )
        plot([data.result(31).s data.result(311).s], [data.result(31).T data.result(311).T] )
        plot([data.result(311).s data.result(32).s], [data.result(311).T data.result(32).T] )
        plot([data.result(32).s data.result(4).s], [data.result(32).T data.result(4).T] )
        plot([data.result(4).s data.result(1).s], [data.result(4).T data.result(1).T] )
        


        hold off;
        
        
        
        subplot(2,2,2);
        hold on;
        title('Rankine-Hirn Cycle')
        xlabel('v[m^3/kg]')
        ylabel('p[Pa]')
        
        plot(data.result(1).v,data.result(1).p,'.','MarkerSize',15)
        text(data.result(1).v,data.result(1).p,'1')
        plot(data.result(20).v,data.result(20).p,'.','MarkerSize',15)
        text(data.result(20).v,data.result(20).p,'20')
        plot(data.result(21).v,data.result(21).p,'.','MarkerSize',15)
        text(data.result(21).v,data.result(21).p,'21')
        plot(data.result(22).v,data.result(22).p,'.','MarkerSize',15)
        text(data.result(22).v,data.result(22).p,'22')
        plot(data.result(3).v,data.result(3).p,'.','MarkerSize',15)
        text(data.result(3).v,data.result(3).p,'3')
        plot(data.result(31).v,data.result(31).p,'.','MarkerSize',15)
        text(data.result(31).v,data.result(31).p,'31')
        plot(data.result(32).v,data.result(32).p,'.','MarkerSize',15)
        text(data.result(32).v,data.result(32).p,'32')
        plot(data.result(4).v,data.result(4).p,'.','MarkerSize',15)
        text(data.result(4).v,data.result(4).p,'4')

        plot([data.result(1).v data.result(20).v], [data.result(1).p data.result(20).p] )
        plot([data.result(20).v data.result(21).v], [data.result(20).p data.result(21).p] )
        plot([data.result(21).v data.result(22).v], [data.result(21).p data.result(22).p] )
        plot([data.result(22).v data.result(3).v], [data.result(22).p data.result(3).p] )
        plot([data.result(3).v data.result(31).v], [data.result(3).p data.result(31).p] )
        plot([data.result(31).v data.result(32).v], [data.result(31).p data.result(32).p] )
        plot([data.result(32).v data.result(4).v], [data.result(32).p data.result(4).p] )
        plot([data.result(4).v data.result(1).v], [data.result(4).p data.result(1).p] )
        hold off;
        grid off;


subplot(2,2,3);
title('ENERGY')
labels = {'pump','steam generator','turbine','condensor'};
sum_energ = [data.result(1).h data.result(20).h+data.result(21).h+data.result(22).h data.result(3).h data.result(4).h];
pie(sum_energ,labels);

subplot(2,2,4);
title('EXERGY')
labels = {'pump','steam generator','turbine','condensor'};
sum_exerg = [data.result(1).ex data.result(20).ex+data.result(21).ex+data.result(22).ex data.result(3).ex data.result(4).ex];
pie(sum_exerg,labels);
    end
    
    
    
end
