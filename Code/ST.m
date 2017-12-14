function [ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG] = ST(P_e,options,display)
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
%   -options.T_max     [??C] : Maximum steam temperature
%   -options.T_cond_out[??C] : Condenseur cold outlet temperature
%   -options.p3_hp     [bar] : Maximum pressure
%   -options.drumFlag  [K] : if =1 then drum if =0 => no drum. 
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax     [??C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.p_3       [-] : High pressure after last reheating
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (titre)
%   -options.T0        [??C] : Reference temperature
%   -options.TpinchSub [??C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [??C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[??C] : Temperature pinch at condenser 
%   -options.Tdrum     [??C] : drum temperature
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
% dat = {T_1       , T_2       , ...       , T_4;  [??C]
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
%
% FIG is a vector of all the figure you plot. Before each figure, define a
% figure environment such as:  

%  "FIG(1) = figure;
%  plot(...);
%  [...]
%   FIG(2) = figure;
%  pie(...);
%  [...]"
%  Your vector FIG will contain all the figure plot during the run of this
%  code (whatever the size of FIG).
%
if nargin<3
    display = 1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 300e6; % [W] Puissance effective de l'installation
        end
    end
end

%% Input initialisation

%nombre de sous intervalle pour etats 2' et 2'' si surcritique
N_discr = 100; 

% Number of feed-heating [-]
if isfield(options,'nsout')           
    n_sout = options.nsout;    
else
    n_sout = 0;  %[-] 
end

% Number of reheating [-]
if isfield(options,'reheat')   
    n_reheat = options.reheat;
    if n_reheat<n_sout || n_reheat==0
        n_reheat = n_sout;
    end
else
    n_reheat = 0;  %[-] 
end

% Drum : if =1 then drum if = 0 => no drum. 
if isfield(options,'drumFlag')           
    drumFlag = options.drumFlag;    
elseif n_sout==0
    drumFlag = 0;    %[-] 
else
    drumFlag = 0;    %[-] 
end

% Drum temperature
if isfield(options,'Tdrum')           
    T_drum = options.Tdrum;    
else
    T_drum = 115; %[?C]   
end

% Maximum steam temperature [?C]
if isfield(options,'T_max')           
    T_max = options.T_max;    
else
    T_max = 525;%565;  %[?C] 
end

% Maximum pressure [bar]
if isfield(options,'p3_hp')           
    p3_hp = options.p3_hp;    
else
    p3_hp = 200;     %[bar] (enonce : 250)
end

% High pressure after last reheating
if isfield(options,'p_3')           
    p_3 = options.p_3;    
else
    p_3 = 70; %[bar]  
end

% Condenseur cold outlet temperature [?C]
if isfield(options,'T_cond_out')           
    Tcond_out = options.T_cond_out;    
else
    Tcond_out = 33; %[?C]
end

% Vapor ratio [gaseous/liquid] (titre)
if isfield(options,'x4')           
    x4 = options.x4;    
else
    x4 = 0.88; %[-]   
end

% Temperature pinch at the subcooler
if isfield(options,'TpinchSub')           
    TPinchSub = options.TpinchSub;    
else
    TPinchSub = 4; %[?C]  
end

% Temperature pinch at a heat exchanger
if isfield(options,'TpinchEx')           
    TPinchEx = options.TpinchEx;    
else
    TPinchEx = 4; %[?C]  
end

% Combustion data : 
if isfield(options,'comb')           
    comb = options.comb;    
else
    comb.Tmax   =  1400;  % [?C] 
    comb.lambda =  1.05;  % [-] 
    comb.x      =  0;
    comb.y      =  4;     %CH4
end

% Mecanic efficiency of shafts bearings
if isfield(options,'eta_mec')           
    eta_mec = options.eta_mec;    
else
    eta_mec = 0.98; %[-] 
end

% Isotrenpic efficiency for compression
if isfield(options,'eta_SiC')           
    eta_SiC = options.eta_SiC;    
else
    eta_SiC = 0.8; %[-]  
end

% Isotrenpic efficiency for Turbine
if isfield(options,'eta_SiT')           
    eta_SiT = options.eta_SiT;    
else %A CHANGER 
    eta_SiT(1) = 0.88; %[-]
    eta_SiT(2) = 0.88; %[-]
    eta_SiT(3) = 0.88; %[-]
end

%% Intermidiary Initialisation

% Reference State :eau etat liquide a 15 degres
if isfield(options,'T_0')           
    data.T0 = options.T_0;    
else
    data.T0 = 15;  %[?C] 
end
h_ref = XSteam('hL_p',XSteam('psat_T',data.T0)); %[kj/kg]
s_ref = XSteam('sL_p',XSteam('psat_T',data.T0)); %[kj/kgK]
v_eau = 1/1000; %(volume massique eau)
    
%% Auxiliary functions    
% Exergy
function [ex] = exergy(h,h_ref,s,s_ref,T0)
    ex = (h-h_ref)-(T0+273.15)*(s-s_ref);
end

% Titre XSteam (correction)
function [x_corr] = titre_corr(x_current,h_current,p)
    if x_current == 0 && h_current<XSteam('hL_p',p)
        x_corr = nan; 
    elseif x_current == 1 && h_current>XSteam('hV_p',p)
        x_corr = nan;
    else 
        x_corr = x_current;
    end    
end


%% Simulation

%% INIT : 
% SORTIE CHAUDIERE APRES(possible) RESURCHAUFFE avant turbine LP
% (vapeur surchaufee)
if n_reheat==0
    etat_3(1,1)  = p_3;
    etat_3(1,2)  = T_max;
    etat_3(1,3)  = XSteam('h_pT',etat_3(1,1),etat_3(1,2));
    etat_3(1,4)  = XSteam('s_pT',etat_3(1,1),etat_3(1,2));
    etat_3(1,5)  = titre_corr(XSteam('x_ph',etat_3(1,1),etat_3(1,3)),etat_3(1,3),etat_3(1,1));
    etat_3(1,6) = exergy(etat_3(1,3),h_ref,etat_3(1,4),s_ref,data.T0);
end

%% Etats 4i et 5i avant et apres resurchauffe i     
if n_reheat>0
    etat_3(1,1)  = p3_hp;
    etat_3(1,2)  = T_max;
    etat_3(1,3)  = XSteam('h_pT',etat_3(1,1),etat_3(1,2));
    etat_3(1,4)  = XSteam('s_pT',etat_3(1,1),etat_3(1,2));
    etat_3(1,5)  = titre_corr(XSteam('x_ph',etat_3(1,1),etat_3(1,3)),etat_3(1,3),etat_3(1,1));
    etat_3(1,6) = exergy(etat_3(1,3),h_ref,etat_3(1,4),s_ref,data.T0);
    
    % SORTIES CHAUDIERE APRES RESURCHAUFFE ---> ??
    etat_5i = zeros(n_reheat,6); % p,T,h,s,x,e
    %p_reheat = linspace(etat_3(1,1),p_3,n_reheat+1);
%     delta_p = (p_3-etat_3(1,1))/(n_reheat);
        etat_5i(n_reheat,1)=p_3;
        etat_5i(n_reheat,2)=T_max;
        etat_5i(n_reheat,3)=XSteam('h_pT',etat_5i(n_reheat,1),etat_5i(n_reheat,2));
        etat_5i(n_reheat,4)=XSteam('s_pT',etat_5i(n_reheat,1),etat_5i(n_reheat,2));
        etat_5i(n_reheat,5)=titre_corr(XSteam('x_ph', etat_5i(n_reheat,1), etat_5i(n_reheat,3)),etat_5i(n_reheat,3),etat_5i(n_reheat,1));
        etat_5i(n_reheat,6)=exergy(etat_5i(n_reheat,3),h_ref,etat_5i(n_reheat,4),s_ref,data.T0);
    for i=1:n_reheat-1
        
%         if i==1
%             etat_5i(i,1) = delta_p*i+etat_3(1,1);   
%             etat_5i(i,2) = T_max;
%             etat_5i(i,3) = XSteam('h_pT', etat_5i(i,1), etat_5i(i,2));
%             etat_5i(i,4) = XSteam('s_pT', etat_5i(i,1), etat_5i(i,2));
%         
%         else 
%         end
%         etat_5i(i,1) = delta_p*i+etat_3(1,1);   
%         etat_5i(i,2) = T_max;
%         etat_5i(i,3) = XSteam('h_pT', etat_5i(i,1), etat_5i(i,2));
%         etat_5i(i,4) = XSteam('s_pT', etat_5i(i,1), etat_5i(i,2));
%         etat_5i(i,5) = titre_corr(XSteam('x_ph', etat_5i(i,1), etat_5i(i,3)),etat_5i(i,3),etat_5i(i,1));
%         etat_5i(i,6) = exergy(etat_5i(i,3), h_ref, etat_5i(i,4), s_ref, data.T0); 
           

           delta_S = (etat_5i(n_reheat,4)-etat_3(1,4))/n_reheat;
           etat_5i(i,4)  = delta_S*i+etat_3(1,4);
           etat_5i(i,2)=T_max;
           f_s = @(p) XSteam('T_ps',p,etat_5i(i,4))-etat_5i(i,2);
           etat_5i(i,1) = fsolve(f_s,0.1,optimoptions('fsolve','Display','off'));
          
           etat_5i(i,3) = XSteam('h_pT',etat_5i(i,1),etat_5i(n_reheat,2));
           etat_5i(i,5) = titre_corr(XSteam('x_ph', etat_5i(i,1), etat_5i(i,3)),etat_5i(i,3),etat_5i(i,1));
           etat_5i(i,6) = exergy(etat_5i(i,3),h_ref,etat_5i(i,4),s_ref,data.T0);
    
    end
    %%BAUMANN : quand on gagne 1% d'humidite, on perd 1% de rendement
    %if (n_reheat(i,4)<XSteam('sV_p',n_reheat(i,1)))
    %    s_L = XSteam('sL_p',n_reheat(i,1));
    %    x = (n_reheat(i,4)-s_L)/(XSteam('sV_p',reheat(i,1))-s_L);
    %    eta_SiT(1) = eta_SiT(1) - (1-x);    
    %end

    %p3 =etat_3(1,1)
    % SORTIES HP Turbine (pour resurchauffe) ---> ??
    etat_4i = zeros(n_reheat,6); % p,T,h,s,x,e
    for i=1:n_reheat 
        etat_4i(i,1) = etat_5i(i,1);
        if i==1
            s_4_is          = etat_3(1,4);
        else
            s_4_is          = etat_5i(i-1,4);
        end
        h_4_is              = XSteam('h_ps', etat_4i(i,1), s_4_is);
        %VERIFIER SI PAS COMMENCER GENRE AVEC UN IF i==1 pour calculer le
        %h_3 puis prendre a chaque fois les valeurs de 4i(i-1) pour les h
        %entre chaque detente
        if i==1
            etat_4i(i,3) = etat_3(1,3) + eta_SiT(1)*(h_4_is-etat_3(1,3));
        else
            etat_4i(i,3) = etat_5i(i-1,3) + eta_SiT(1)*(h_4_is-etat_5i(i-1,3));
        end
        etat_4i(i,2) = XSteam('T_ph', etat_4i(i,1), etat_4i(i,3));
        etat_4i(i,4) = XSteam('s_ph', etat_4i(i,1), etat_4i(i,3));
        etat_4i(i,5) = titre_corr(XSteam('x_ph',etat_4i(i,1), etat_4i(i,3)),etat_4i(i,3), etat_4i(i,1));
        etat_4i(i,6) = exergy(etat_4i(i,3), h_ref, etat_4i(i,4), s_ref, data.T0);
    end
    reheat_T=[];
    reheat_s=[];
    for i=1:n_reheat
        reheat_T = [reheat_T,[etat_4i(i,2)' etat_5i(i,2)']];
        reheat_s = [reheat_s,[etat_4i(i,4)' etat_5i(i,4)']];
    end
end %end if reheat
   
%% LP Turbine : SORTIE ----> ?? TRAITER CAS OU x4 est donne!
% etat isentropique : hyp : temperature fin de detente idem condenseur
% en fin de detente, le fluide est en partie liquide, en partie vapeur    
    etat_6(1,2) = Tcond_out;
    etat_6(1,1) = XSteam('psat_t', etat_6(1,2));
% etat isentropique    
    if n_reheat==0
        s_6_is = etat_3(1,4); 
        h_6_is = XSteam('h_ps', etat_6(1,1), s_6_is);
        etat_6(1,3)  = etat_3(1,3) + eta_SiT(3)*(h_6_is-etat_3(1,3));
    else
        s_6_is = etat_5i(n_reheat,4); %derniere resurchauffe
        h_6_is = XSteam('h_ps', etat_6(1,1), s_6_is);
        etat_6(1,3)  = etat_5i(n_reheat,3) + eta_SiT(3)*(h_6_is-etat_5i(n_reheat,3));
    end

% etat reel
    etat_6(1,4)  = XSteam('s_ph', etat_6(1,1), etat_6(1,3));
    etat_6(1,5)  = titre_corr(XSteam('x_ph', etat_6(1,1), etat_6(1,3)),etat_6(1,3),etat_6(1,1));
    etat_6(1,6) = exergy(etat_6(1,3), h_ref, etat_6(1,4), s_ref, data.T0);

% SORTIE CONDENSEUR ---> OK
% liquide sature en sortie de condenseur
    etat_7(1,2)  = Tcond_out;
    etat_7(1,1)  = XSteam('psat_T',etat_7(1,2));
    etat_7(1,3)  = XSteam('hL_T',etat_7(1,2));
    etat_7(1,4)  = XSteam('sL_T',etat_7(1,2));
    etat_7(1,5)  = titre_corr(XSteam('x_ph',etat_7(1,1),etat_7(1,3)),etat_7(1,3),etat_7(1,1));
    %VERIFIER QUE BIEN TITRE NUL -> liquide sature
    etat_7(1,6) = exergy(etat_7(1,3), h_ref, etat_7(1,4), s_ref, data.T0);
 
%% Soutirages - etats soutire des turbines 6i
if n_sout>0  

if n_reheat>0    
    etat_6i = zeros(n_sout,6); % p,T,h,s,x,e
    delta_h_sout = (etat_5i(n_reheat,3)-etat_6(1,3))/(n_sout+1-n_reheat);
else
    etat_6i = zeros(n_sout,6); % p,T,h,s,x,e
    delta_h_sout = (etat_3(1,3)-etat_6(1,3))/(n_sout+1);
end
for i=1:n_sout 
    if n_reheat>0 && i>(n_sout-n_reheat) % de HP
        etat_6i(i,1) = etat_4i(n_sout-i+1,1); %p
        etat_6i(i,2) = etat_4i(n_sout-i+1,2); %T
        etat_6i(i,3) = etat_4i(n_sout-i+1,3); %h
        etat_6i(i,4) = etat_4i(n_sout-i+1,4); %s
        etat_6i(i,5) = etat_4i(n_sout-i+1,5); %x
        etat_6i(i,6) = etat_4i(n_sout-i+1,6); %ex
    else % de LP
        if n_reheat==0
            etat_6i(i,3) = etat_6(1,3) + delta_h_sout*i;
            h_6i_is           = etat_3(1,3) - (etat_3(1,3)-etat_6i(i,3))/eta_SiT(3);
            s_6i_is           = etat_3(1,4);
        else
            etat_6i(i,3) = etat_6(1,3) + delta_h_sout*i;
            h_6i_is           = etat_5i(n_reheat,3) - (etat_5i(n_reheat,3)-etat_6i(i,3))/eta_SiT(3);
            s_6i_is           = etat_5i(n_reheat,4);
        end
        
        etat_6i(i,1) = XSteam('p_hs',h_6i_is,s_6i_is); %idem que p_6i_is
        etat_6i(i,2) = XSteam('T_ph',etat_6i(i,1),etat_6i(i,3));
        etat_6i(i,4) = XSteam('s_ph',etat_6i(i,1),etat_6i(i,3));
        etat_6i(i,5) = titre_corr(XSteam('x_ph',etat_6i(i,1),etat_6i(i,3)),etat_6i(i,3),etat_6i(i,1));
        etat_6i(i,6) = exergy(etat_6i(i,3), h_ref, etat_6i(i,4), s_ref, data.T0);
        
        %%BAUMANN : quand on gagne 1% d'humidite, on perd 1% de rendement
        %if (etat_6i(i,4)<XSteam('sV_p',etat_6i(i,1)))
        %    s_L = XSteam('sL_p',etat_6i(i,1));
        %    x = (etat_6i(i,4)-s_L)/(XSteam('sV_p',etat_6i(i,1))-s_L);
        %    eta_SiT(3) = eta_SiT(3) - (1-x);
        %end
    end
end

%% SORTIE RECHAUFFEURS R_I, R_II,...R_nsout etats 7i
% amene le fluide chaud de maniere isobare jusque etat de liquide sature
% liquide sature a la pression de soutirage
etat_7i      = zeros(n_sout,6); % p,T,h,s,x,e 

for i=1:n_sout
    etat_7i(i,1) = etat_6i(i,1); %isobare
    etat_7i(i,2) = XSteam('Tsat_p',etat_7i(i,1));
    etat_7i(i,3) = XSteam('hL_p',etat_7i(i,1));
    etat_7i(i,4) = XSteam('sL_p',etat_7i(i,1));
    etat_7i(i,5) = titre_corr(XSteam('x_ph',etat_7i(i,1),etat_7i(i,3)),etat_7i(i,3),etat_7i(i,1));
    etat_7i(i,6) = exergy(etat_7i(i,3), h_ref, etat_7i(i,4), s_ref, data.T0);
end

%trier les temperatures dans l'ordre croissant afin de disposer 
%les echangeurs dans le bon ordre 
[~,new_index] = sort(etat_7i(:,2));  

%% SORTIE principale ECHANGEURS 9_i
% Echangeurs non parfaits : tpinch
etat_9i      = zeros(n_sout,6); % p,T,h,s,x 
etat_9i(:,2) = etat_7i(new_index(:),2) - TPinchEx;
if drumFlag==1
    % imposer Pe telle que la pression en 9_nsout (ou avant le degazeur) soit suffisamment ?lev?e pour
    % qu'a la temperature de 9_nsout on soit toujours a l'etat liquide et ne
    % pas avoir de vapeur (donc prendre Pe> Psat(T_9_nsout)
    
    % Position degazeur
    drum_index = 2; %
    while etat_9i(drum_index-1,2) < T_drum  
        drum_index=drum_index+1;
    end
    etat_8(1,1)               = XSteam('psat_T',etat_9i(drum_index,2));
    etat_9i(1:drum_index-1,1) = etat_8(1,1);
    %CHANGER le +1 en un + interessant (10? 35?)
    etat_9i(drum_index:end,1) = XSteam('psat_T',max(etat_9i(:,2)))+35;
    h_9s_af_drum = etat_7i(new_index(drum_index),3) + v_eau*(etat_9i(drum_index,1)-etat_7i(new_index(drum_index),1))*(1e5)/(1e3);
    etat_9i(drum_index,3) = etat_7i(new_index(drum_index),3) + (h_9s_af_drum-etat_7i(new_index(drum_index),3))/eta_SiC;
    %changer le 0.04 (pas vraiment necessaire temperature quasi constante
    etat_9i(drum_index,2)     = etat_7i(new_index(drum_index),2) + (etat_9i(drum_index,3)-etat_7i(new_index(drum_index),3)-0.04*(etat_9i(drum_index,1)-etat_7i(new_index(drum_index),1)))/4.18;
else
    etat_8(1,1)  = XSteam('psat_T',max(etat_9i(:,2)))+35; %A CHANGER UNIQUEMENT TEST
    etat_9i(:,1) = etat_8(1,1);
end
for i=1:n_sout
    if etat_9i(i,1) == XSteam('psat_T',etat_9i(i,2))
        etat_9i(i,3) = XSteam('hL_p',etat_9i(i,1));
        etat_9i(i,4) = XSteam('sL_p',etat_9i(i,1));
    else
        etat_9i(i,3) = XSteam('h_pT',etat_9i(i,1),etat_9i(i,2));
        etat_9i(i,4) = XSteam('s_pT',etat_9i(i,1),etat_9i(i,2));
    end
    etat_9i(i,5) = titre_corr(XSteam('x_ph',etat_9i(i,1),etat_9i(i,3)),etat_9i(i,3),etat_9i(i,1));
    etat_9i(i,6) = exergy(etat_9i(i,3),h_ref,etat_9i(i,4),s_ref, data.T0);
end

%% SORTIE POMPE d'extraction Pe, (ENTREE R0), etat 8
h_8_is    = etat_7(1,3) + v_eau*(etat_8(1,1)-etat_7(1,1))*(1e5)/(1e3);
etat_8(1,3)  = etat_7(1,3) - (etat_7(1,3)-h_8_is)/eta_SiC;
etat_8(1,2)  = etat_7(1,2) + ((etat_8(1,3)-etat_7(1,3))-0.09*(etat_8(1,1)-etat_7(1,1)))/4.18;
etat_8(1,4)  = XSteam('s_pT',etat_8(1,1),etat_8(1,2));
etat_8(1,5)  = titre_corr(XSteam('x_ph', etat_8(1,1), etat_8(1,3)),etat_8(1,3),etat_8(1,1));
etat_8(1,6) = exergy(etat_8(1,3),h_ref,etat_8(1,4),s_ref,data.T0);    

%% Etat 10 sortie R0 sous-refroidisseur des soutirages
etat_10(1,2)  = etat_8(1,2) + TPinchSub;
etat_10(1,1)  = etat_7i(1,1);
etat_10(1,3)  = XSteam('h_pT',etat_10(1,1),etat_10(1,2));
etat_10(1,4)  = XSteam('s_pT',etat_10(1,1),etat_10(1,2));
etat_10(1,5)  = titre_corr(XSteam('x_ph',etat_10(1,1),etat_10(1,3)),etat_10(1,3),etat_10(1,1));
etat_10(1,6) = exergy(etat_10(1,3),h_ref,etat_10(1,4),s_ref,data.T0);

%% Etat 11 apres detente de 10
etat_11(1,1)  = etat_7(1,1);
etat_11(1,3)  = etat_10(1,3);
etat_11(1,2)  = XSteam('T_ph',etat_11(1,1),etat_11(1,3));
etat_11(1,4)  = XSteam('s_ph',etat_11(1,1),etat_11(1,3));
etat_11(1,5)  = titre_corr(XSteam('x_ph',etat_11(1,1),etat_11(1,3)),etat_11(1,3),etat_11(1,1)); 
etat_11(1,6) = exergy(etat_11(1,3),h_ref,etat_11(1,4),s_ref,data.T0);

%% Etat 1 = Etat dernier soutirage
etat_1(1,1)  = etat_9i(n_sout,1);
etat_1(1,2)  = etat_9i(n_sout,2);
etat_1(1,3)  = etat_9i(n_sout,3);
etat_1(1,4)  = etat_9i(n_sout,4);
etat_1(1,5)  = etat_9i(n_sout,5);
etat_1(1,6) = etat_9i(n_sout,6);

%% SOUTIRAGES : 
etat_9(1,3) = etat_8(1,3); %valeur de depart
error            = 1;
% Fraction soutirages (pas de drum)
%resoudre systeme Ax = b avec X les soutirages
%remplissage termes independants
b = zeros(n_sout,1);
A = zeros(n_sout,n_sout);
while error>0.01
if drumFlag==0
    for i=1:n_sout
        if i==1
            b(i) = etat_9i(i,3) - etat_9(1,3);  
        else
            b(i) = etat_9i(i,3) - etat_9i(i-1,3);
        end
        for j=1:n_sout
            if i==1
                if j>1
                    A(i,j) = -(etat_9i(i,3)-etat_9(1,3))+(etat_7i(new_index(i+1),3)-etat_7i(new_index(i),3));
                else 
                    A(i,j) = -(etat_9i(i,3)-etat_9(1,3))+(etat_6i(new_index(i),3)-etat_7i(new_index(i),3)); 
                end
            elseif i<j
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3))+(etat_7i(new_index(i+1),3)-etat_7i(new_index(i),3));
            elseif i==j
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3))+(etat_6i(new_index(i),3)-etat_7i(new_index(i),3));
            else
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3));
            end
        end
    end
end

% Fraction soutirages (avec drum degazeur)
%resoudre systeme Ax = b avec X les soutirages

if drumFlag==1
for i=1:n_sout
    if i==1
        b(i) = etat_9i(i,3) - etat_9(1,3);
    elseif i==drum_index
        b(i) = etat_7i(new_index(i),3) - etat_9i(i-1,3);
    else
        b(i) = etat_9i(i,3) - etat_9i(i-1,3);
    end
    for j=1:n_sout
       if ((i<drum_index && j<drum_index) || (i>drum_index && j>drum_index))
            if i==1
                if j>1
                    A(i,j) = -(etat_9i(i,3)-etat_9(1,3))+(etat_7i(new_index(i+1),3)-etat_7i(new_index(i),3));
                else 
                    A(i,j) = -(etat_9i(i,3)-etat_9(1,3))+(etat_6i(new_index(i),3)-etat_7i(new_index(i),3)); 
                end
            elseif i<j
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3))+(etat_7i(new_index(i+1),3)-etat_7i(new_index(i),3));
            elseif i==j
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3))+(etat_6i(new_index(i),3)-etat_7i(new_index(i),3));
            else
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3));
            end
       elseif i==drum_index 
            if i>j
                A(i,j) = -(etat_7i(new_index(i),3) - etat_9i(i-1,3));
            elseif i==j
                A(i,j) =  (etat_6i(new_index(i),3) - etat_7i(new_index(i),3));
            else 
                A(i,j) =  (etat_7i(new_index(i+1),3)-etat_7i(new_index(i),3));
            end
       elseif i>drum_index && j<=drum_index
                A(i,j) = -(etat_9i(i,3)-etat_9i(i-1,3));
       end
    end
end
end %end drumflag1
A=A;
b=b;
%Resolution du systeme
X = A\b;
X=X';

%Calcul du nouveau h90
if drumFlag==1
    sum_X   = sum(X(1:drum_index-1));
    h90_new = etat_8(1,3) + (etat_7i(new_index(1),3)-etat_10(1,3))*(sum_X/(1+sum_X));
else
    sum_X   = sum(X(1:end));
    h90_new = etat_8(1,3) + (etat_7i(new_index(1),3)-etat_10(1,3))*(sum_X/(1+sum_X));
end
error       = abs(h90_new-etat_9(1,3));
etat_9(1,3) = h90_new;
end %end while
%X_sorted  = zeros(1,n_sout);
X_sorted  = X(new_index); 
X(:)      = X_sorted(:);

%% Etat 9_0, sortie du R0
etat_9(1,1)  = etat_8(1,1);
etat_9(1,2)  = XSteam('T_ph',etat_9(1,1),etat_9(1,3));
etat_9(1,4)  = XSteam('s_ph',etat_9(1,1),etat_9(1,3));
etat_9(1,5)  = titre_corr(XSteam('x_ph',etat_9(1,1),etat_9(1,3)),etat_9(1,3),etat_9(1,1));
etat_9(1,6) = exergy(etat_9(1,3),h_ref,etat_9(1,4),s_ref,data.T0);

else
%% Etat 1, si pas de soutirages   
    etat_1(1,1)  = etat_7(1,1);
    etat_1(1,2)  = etat_7(1,2);
    etat_1(1,3)  = etat_7(1,3);
    etat_1(1,4)  = etat_7(1,4);
    etat_1(1,5)  = etat_7(1,5);
    etat_1(1,6) = etat_7(1,6);
end %end n_sout>0

%% Steam Generator
% Etat 2
    etat_2(1,1)  = etat_3(1,1); 
    h2_s         = etat_1(1,3) + v_eau*((etat_2(1,1)-etat_1(1,1)))*(1e5)/(1e3); %delta_h = w_m car isentropique adiabatique, q=0, w_f=0)
    etat_2(1,3)  = etat_1(1,3) - (etat_1(1,3)-h2_s)/eta_SiC;
    etat_2(1,2)  = etat_1(1,2)+((etat_2(1,3)-etat_1(1,3))-0.001*(etat_2(1,1)-etat_1(1,1)))/4.18 ; %difference de temperature negligeable CHANGER
    etat_2(1,4)  = XSteam('s_pT',etat_2(1,1),etat_2(1,2));
    etat_2(1,5)  = titre_corr(XSteam('x_ph',etat_2(1,1),etat_2(1,3)),etat_2(1,3),etat_2(1,1));
    etat_2(1,6) = exergy(etat_2(1,3), h_ref,etat_2(1,4),s_ref,data.T0);
T2_sur_test = XSteam('Tsat_p',etat_2(1,1)); 
if isfinite(T2_sur_test)
    supercritique=0;
else
    supercritique=1;
    T2_to_T3 = linspace(etat_2(1,2),etat_3(1,2),N_discr);
    s2_to_s3 = zeros(size(T2_to_T3));
    h2_to_h3 = zeros(size(T2_to_T3));
    for i=1:N_discr
        s2_to_s3(i)=XSteam('s_pT',etat_2(1,1),T2_to_T3(i));
        h2_to_h3(i)=XSteam('h_pT',etat_2(1,1),T2_to_T3(i));
    end
end
% Etat 2'
    etat_21(1,1)  = etat_3(1,1);
    etat_21(1,2)  = XSteam('Tsat_p',etat_21(:,1)');
    etat_21(1,3)  = XSteam('hL_p', etat_21(:,1)');
    etat_21(1,4)  = XSteam('sL_p',etat_21(:,1)');
    etat_21(1,5)  = titre_corr(XSteam('x_ph',etat_21(:,1)',etat_21(:,3)'),etat_21(:,3)',etat_21(:,1)');
    etat_21(1,6) = exergy(etat_21(:,3)', h_ref, etat_21(:,4)', s_ref, data.T0);

% Etat 2''    
    etat_22(1,1) = etat_3(1,1);
    etat_22(1,2) = XSteam('Tsat_p',etat_21(:,1)');
    etat_22(1,3) = XSteam('hV_p', etat_21(:,1)');
    etat_22(1,4) = XSteam('sV_p',etat_21(:,1)');
    etat_22(1,5) = titre_corr(XSteam('x_ph',etat_22(1,1),etat_22(1,3)),etat_22(1,3),etat_22(1,1));
    etat_22(1,6) = exergy(etat_22(1,3), h_ref, etat_22(1,4), s_ref, data.T0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Travail moteur : turbines et Chaudiere
%CHAUDIERE
Q_II     = 0;
Q_II_ex  = 0;
if n_reheat>0
    for i=1:n_reheat+1
        if i==1
            Q_II     = Q_II + (etat_3(1,3)-etat_2(1,3));
            Q_II_ex  = Q_II_ex + (etat_3(1,6)-etat_2(1,6));
        elseif i==n_reheat+1
           % SUMPLUS1 = 1+sum(X(1:n_sout-1+1))
            if n_sout>0
                Q_II     = Q_II    + (etat_5i(n_reheat,3)-etat_4i(n_reheat,3))*(1+sum(X(1:n_sout-i+1)))/(1+sum(X));
                Q_II_ex  = Q_II_ex + (etat_5i(n_reheat,6)-etat_4i(n_reheat,6))*(1+sum(X(1:n_sout-i+1)))/(1+sum(X));
            else
                Q_II     = Q_II    + (etat_5i(n_reheat,3)-etat_4i(n_reheat,3));
                Q_II_ex  = Q_II_ex + (etat_5i(n_reheat,6)-etat_4i(n_reheat,6));
            end
        else
            if n_sout>0
                Q_II     = Q_II    + (etat_5i(i-1,3)-etat_4i(i,3))*(1+sum(X(1:n_sout-i+1)))/(1+sum(X));
                Q_II_ex  = Q_II_ex + (etat_5i(i,6)-etat_4i(i,6))*(1+sum(X(1:n_sout-i+1)))/(1+sum(X));
            else
                Q_II     = Q_II    + (etat_5i(i-1,3)-etat_4i(i,3));
                Q_II_ex  = Q_II_ex + (etat_5i(i,6)-etat_4i(i,6));
            end
        end
    end
else %RANKINE
    Q_II     = Q_II + (etat_3(1,3)-etat_2(1,3));
    Q_II_ex  = Q_II + (etat_3(1,6)-etat_2(1,6));
end

Q_II        = Q_II*1e3;            %[J]
Q_II_ex     = Q_II_ex*1e3;         %[J]

%TRAVAIL MOTEUR
wm_HP    = 0;
wm_HP_ex = 0;
wm_LP    = 0;
wm_LP_ex = 0;
if n_sout>0
    if n_reheat==0 
            for i=1:n_sout
                if i==n_sout 
                    wm_LP    = wm_LP    + (etat_3(1,3)-etat_6i(i,3))*(1+sum(X(1:i-1)))/(1+sum(X));
                    wm_LP_ex = wm_LP_ex + (etat_3(1,6)-etat_6i(i,6))*(1+sum(X(1:i-1)))/(1+sum(X));
                elseif i==1
                    wm_LP    = wm_LP    + (etat_6i(i,3)-etat_6(1,3))*1/(1+sum(X));
                    wm_LP_ex = wm_LP_ex + (etat_6i(i,6)-etat_6(1,6))*1/(1+sum(X));
                else
                    wm_LP    = wm_LP    + (etat_6i(i+1,3)-etat_6i(i,3))*(1+sum(X(1:i-1)))/(1+sum(X));
                    wm_LP_ex = wm_LP_ex + (etat_6i(i+1,6)-etat_6i(i,6))*(1+sum(X(1:i-1)))/(1+sum(X));
                end
            end
    else 
        for i=1:n_sout+1
            if i>n_sout+1-n_reheat
                if i==n_sout+1 %ICI
                    wm_HP    = wm_HP     + (etat_3(1,3)-etat_4i(1,3))*(1+sum(X(1:end-1)))/(1+sum(X));
                    wm_HP_ex = wm_HP_ex  + (etat_3(1,6)-etat_4i(1,6))*(1+sum(X(1:end-1)))/(1+sum(X));
                else
                    wm_HP    = wm_HP     + (etat_5i(n_sout-i+1,3)-etat_4i(n_sout-i+2,3))*(1+sum(X(1:i-2)))/(1+sum(X));
                    wm_HP_ex = wm_HP_ex  + (etat_5i(n_sout-i+1,6)-etat_4i(n_sout-i+2,6))*(1+sum(X(1:i-2)))/(1+sum(X));
                end
            elseif i==1 
                    wm_LP    = wm_LP    + (etat_6i(i,3)-etat_6(1,3))*1/(1+sum(X));
                    wm_LP_ex = wm_LP_ex + (etat_6i(i,6)-etat_6(1,6))*1/(1+sum(X));
            elseif i==n_sout+1-n_reheat 
                    wm_LP    = wm_LP    + (etat_5i(n_reheat,3)-etat_6i(i-1,3))*(1+sum(X(1:i-1)))/(1+sum(X));
                    wm_LP_ex = wm_LP_ex + (etat_5i(n_reheat,6)-etat_6i(i-1,6))*(1+sum(X(1:i-1)))/(1+sum(X));
            else 
                wm_LP = wm_LP       + (etat_6i(i,3)-etat_6i(i-1,3))*(1+sum(X(1:i-1)))/(1+sum(X));
                wm_LP_ex = wm_LP_ex + (etat_6i(i,6)-etat_6i(i-1,6))*(1+sum(X(1:i-1)))/(1+sum(X));
            end
        end
    end
    wm_HP    = wm_HP+wm_LP;       %[kJ]
    wm_HP_ex = wm_HP_ex+wm_LP_ex; %[kJ]    
elseif n_sout==0 && n_reheat>0 %OOOOK
	for i=1:n_reheat+1
        if i==1
            wm_HP    = wm_HP    + (etat_3(i,3)-etat_4i(i,3));
            wm_HP_ex = wm_HP_ex + (etat_3(i,6)-etat_4i(i,6));
        elseif i==n_reheat+1
            wm_HP    = wm_HP     + (etat_5i(n_reheat,3)-etat_6(1,3));
            wm_HP_ex = wm_HP_ex  + (etat_5i(n_reheat,6)-etat_6(1,6));
        else
            wm_HP    = wm_HP    + (etat_5i(i-1,3)-etat_4i(i,3));
            wm_HP_ex = wm_HP_ex + (etat_5i(i-1,6)-etat_4i(i,6));
        end
    end   
else %Rankine Hirn
     wm_HP    = etat_3(1,3)-etat_6(1,3);     %[kJ]  
     wm_HP_ex = etat_3(1,6)-etat_6(1,6);   %[kJ]
end
wm_turb_eff = (wm_HP)*eta_mec; %[kJ]
wm_turb     = (wm_HP);             %[kJ]
wm_turb_ex  = (wm_HP_ex);           %[kJ]
%Travail moteur effectif pompes
 if n_sout>0
     if drumFlag==0
         wm_pump    = (etat_2(1,3)-etat_1(1,3)) ... %Pa
                    + (etat_8(1,3)-etat_7(1,3));    %Pe 
         wm_pump_ex = (etat_2(1,6)-etat_1(1,6)) ... %Pa
                    + (etat_8(1,6)-etat_7(1,6));    %Pe     
     else
         wm_pump    = (etat_2(1,3)-etat_1(1,3)) ... %Pa
                    + (etat_8(1,3)-etat_7(1,3))*(1+sum(X(1:drum_index-1)))/(1+sum(X)) ... %Pe 
                    + (etat_9i(drum_index,3)-etat_7i(drum_index,3)); %Pb
         wm_pump_ex = (etat_2(1,6)-etat_1(1,6)) ... %Pa
                    + (etat_8(1,6)-etat_7(1,6))*(1+sum(X(1:drum_index-1)))/(1+sum(X)) ... %Pe 
                    + (etat_9i(drum_index,6)-etat_7i(drum_index,6)); %Pb
     end
 else
         wm_pump    = (etat_2(1,3)-etat_1(1,3));    %Pa
         wm_pump_ex = (etat_2(1,6)-etat_1(1,6));   %Pa
 end
 wm_pump_eff = wm_pump/eta_mec;             %[kJ]
 
%travail total
wm_tot_eff  = (wm_turb_eff-wm_pump_eff)*1e3;    %[J]
wm_tot      = (wm_turb-wm_pump)*1e3;            %[J]
wm_tot_ex   = (wm_turb_ex-wm_pump_ex)*1e3;      %[J]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %% Combustion
MM_CO2   = (12+2*16)*1e-3; %[kg/mol]
MM_H2O   = (2*1+16)*1e-3;  %[kg/mol]
MM_O2    = (2*16)*1e-3;    %[kg/mol]
MM_N2    = (2*14)*1e-3;    %[kg/mol]
y        = comb.y;
x        = comb.x;
lambda   =  comb.lambda; 
%lambda = 1.05; %test
w_stoech = 1+(y-2*x)/4;
w        = lambda*w_stoech;
m_a1     = (32+3.76*28)*(1+(y-2*x)/4)/(12+y+16*x);  %[kg_air/kg_fuel]
%m_a1   = 17.2; %test

n_CO2  = 1;                        %[mol]
n_H2O  = y/2;                      %[mol]
n_O2   = w-n_CO2-n_H2O/2-x/2;      %[mol]
n_N2   = 3.76*w;                   %[mol]

n_f  = n_CO2+n_H2O + n_O2 + n_N2;                              %[mol]
MM_f = (n_CO2*MM_CO2+n_H2O*MM_H2O+n_N2*MM_N2+n_O2*MM_O2); %[kg/mol]
MM_fuel = (12+y*1+x*16)*1e-3;                                 %[kg/mol]

%PCI du carburant [kJ/kmol]
PCI = 393400 + 102250*y-x/(1+y/2)*(111000+102250*y); %[kJ/kmol]%approximation
PCI = PCI/MM_fuel;  %[J/kg_fuel]
%PCI = 50150e3; %[J/kg_CH4] %test

%REVOIR POUR TOUS LES CAS!!!
%exergie du carburant
%e_c = PCI*f : (p.211 livre1)
%  f=...1.04...(gaz naturel, charbon) 
%  f=...1.06...(produits petroliers)
e_c        = 1.04*PCI;  %[J]
Tmax_comb  = comb.Tmax;
T_exhaust  = 120;
T_a        = 15;
n_comb_f   = [n_CO2 n_H2O n_O2 n_N2];
T0         = data.T0;

h_f_m_echap = abs(h_fm(T0+273.15,T_exhaust+273.15,n_comb_f));     %[J/(kg)]
h_f_m       = abs(h_fm(T0+273.15,Tmax_comb+273.15,n_comb_f));    %[J/(kg)]
cp_exhaust  = h_f_m_echap/(T_exhaust-T0);                         %[J/(kg*K)]
%cp_exhaust = cp_fm(T0+273.15,T_fechap+273.15,n_comb_f,'no_bar') %[J/(kg*K)]
cp_f_m      =  abs(h_f_m)/(Tmax_comb-T0);                        %[J/(kg*K)]
%cp_f  = cp_fm(T0+273.15,Tmax_comb+273.15,n_comb_f,'no_bar');    %[J/(kg*K)]
cp_f_m_bar  = cp_fm(T0+273.15,Tmax_comb+273.15,n_comb_f,'bar');  %[J/(kg*K)]
cp_a_m      = cp_am(T0+273.15,T_a+273.15) ;                      %[J/(kg*K)]


%exergie de fumee :  (1.72) p.31
%h_a et h_c, enthalpies sensibles de l'air et combustible, negligees
e_f       = PCI/(lambda*m_a1+1)-cp_f_m_bar*(T0+273.15)*log(1+PCI/((lambda*m_a1+1)*cp_f_m*(T0+273.15))); %[J]
e_exhaust = h_f_m_echap-(T0+273.15)*cp_exhaust*log((T_exhaust+273.15)/(T0+273.15));                      %[J]

if T0 == T_a
    dh_a = 0;
    ds_a = 0;
else
    dh_a = cp_am(T0+273.15,T_a+273.15) * (T_a - T0) ;
    ds_a = cp_aml(T0+273.15,Ta+273.15) * log((T_a+273.15)/(T0+273.15)); 
end
e_a = dh_a - T0 * ds_a; %[J]
e_cr =  0 ; % exergie du combustible associee a son enthalpie sensible : negligeable
e_r = (e_a * (lambda*m_a1)+ e_cr)/(lambda*m_a1+1); 

%fraction energetique perdu a la cheminee
eps_exhaust = ((lambda*m_a1+1)*cp_exhaust*(T_exhaust+273.15)-lambda*m_a1*cp_a_m*(T_a+273.15))/PCI;
%deperditions parietales (1%)
eps_p       = 0.01;
%rendement energetique GV
eta_gen     = 1-eps_p-eps_exhaust; 
%PROBLEME : definition du cp_f utilisee dans la fonction cp_fm ne semble 
%           etonnement pas donner de resultat correct avec le livre, 
%           passage donc avec les delta d'enthalpie semble resoudre 
%           le "probleme"

COMBUSTION.LHV    = PCI; 
COMBUSTION.e_c    = e_c;
COMBUSTION.lambda = lambda;
COMBUSTION.Cp_g   = cp_f_m;  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%debit total de vapeur d'eau
debit_m_v     = P_e/(wm_tot_eff);               %[kg/s]
%debit total de combustible
debit_m_fuel  = debit_m_v*Q_II/(eta_gen*PCI);   %[kg/s]
%debit total d'air 
debit_m_air   = debit_m_fuel*lambda*m_a1;      %[kg/s]
%debit total des fumees
debit_m_f     = debit_m_fuel*(1+lambda*m_a1);    %[kg/s]
%debits composition des fumees
debit_m_f_O2  = debit_m_f*(n_O2*MM_O2)/(n_f*MM_f);   %[kg/s]
debit_m_f_N2  = debit_m_f*(n_N2*MM_N2)/(n_f*MM_f);   %[kg/s]
debit_m_f_CO2 = debit_m_f*(n_CO2*MM_CO2)/(n_f*MM_f); %[kg/s]
debit_m_f_H2O = debit_m_f*(n_H2O*MM_H2O)/(n_f*MM_f); %[kg/s]
%debit soutires
if n_sout>0
XMASSFLOW = (debit_m_v*(X/(1+sum(X))))' ;
else
   XMASSFLOW = 0;
end%[kg/s]
MASSFLOW(1) = debit_m_air;  %[kg/s]
MASSFLOW(2) = debit_m_v;    %[kg/s]
MASSFLOW(3) = debit_m_fuel; %[kg/s]
MASSFLOW(4) = debit_m_f;    %[kg/s]

%composition fumees
COMBUSTION.fum = zeros(1,4);
COMBUSTION.fum(1) = debit_m_f_O2;
COMBUSTION.fum(2) = debit_m_f_N2;
COMBUSTION.fum(3) = debit_m_f_CO2;
COMBUSTION.fum(4) = debit_m_f_H2O;
fum1 = COMBUSTION.fum(1);
fum2 = COMBUSTION.fum(2);
fum3 = COMBUSTION.fum(3);
fum4 = COMBUSTION.fum(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rendements energitque et exergetiques

%COMMENTER
eta_cyclen  = wm_tot/(Q_II);                                   %---->OK
eta_toten   = eta_mec*eta_cyclen*eta_gen;                      %---->OK
eta_gen;                                                       %---->~OK!!
eta_cyclex  = wm_tot/(Q_II_ex);                                %---->OK
eta_totex   = P_e/(debit_m_fuel*e_c);                          %---->OK
eta_gex     = (debit_m_v*Q_II_ex)/(debit_m_fuel*e_c);          %---->OK
eta_combex  = debit_m_f*e_f/(debit_m_fuel*e_c);                %---->OK
eta_chemex  = (e_f-e_exhaust)/(e_f-e_r);                       %---->OK
eta_transex = (debit_m_v*Q_II_ex)/(debit_m_f*(e_f-e_exhaust)); %---->OK
eta_rotex   = wm_tot/wm_tot_ex;                                %---->OK

ETA(1) = eta_cyclen;
ETA(2) = eta_toten;
ETA(3) = eta_cyclex;
ETA(4) = eta_totex;
ETA(5) = eta_gen;
ETA(6) = eta_gex;
ETA(7) = eta_combex;
ETA(8) = eta_chemex;
ETA(9) = eta_transex;

%% Pertes energetiques et exergetiques


%flux d'energie primaire
P_toten   = debit_m_fuel*PCI;               %[W] ---->OK
P_gen     = (1-eta_gen)*(debit_m_v*Q_II);   %[W] ---->OK
P_mecen   = P_e*(1-eta_mec);           %[W] ---->OK
P_conden_old  = P_toten*eta_gen*(1-eta_cyclen); %[W] ---->OK

if n_sout>0
P_conden   = (debit_m_v-sum(XMASSFLOW))*(etat_6(1,3)-etat_7(1,3))*1e3;
else
P_conden   =  debit_m_v*(etat_6(1,3)-etat_7(1,3))*1e3 ;  
end

DATEN(1) = P_gen;    %[W]
DATEN(2) = P_mecen;  %[W]
DATEN(3) = P_conden; %[W]

%flux d'exergie primaire
P_totex   = debit_m_fuel*e_c;                              %[W] ---->OK
P_mecex   = P_e*(1-eta_mec);                          %[W] ---->OK
P_rotex   = wm_tot_ex*(1-eta_rotex)*debit_m_v;             %[W] ---->OK
P_combex  = P_totex*(1-eta_combex);                        %[W] ---->OK
P_chemex  = P_totex*eta_combex*(1-eta_chemex);             %[W] ---->OK
P_transex = P_totex*eta_combex*eta_chemex*(1-eta_transex); %[W] ---->OK 
P_condex  = (debit_m_v-sum(XMASSFLOW))*(etat_6(1,6)-etat_7(1,6))*1e3; 
%P_condex  = P_totex*eta_gex*(1-eta_cyclex)-P_rotex;        %[W] ---->OK      

DATEX(1) = P_mecex;   %[W]
DATEX(2) = P_totex;   %[W]
DATEX(3) = P_rotex;   %[W]
DATEX(4) = P_combex;  %[W]
DATEX(5) = P_condex;  %[W]
DATEX(6) = P_chemex;  %[W]
DATEX(7) = P_transex; %[W]

% Etat vapeur sature lors du soutirage (pour plot)
if n_sout>0
    etat_6i_sat = zeros(n_sout,5);
    for i=1:n_sout
        etat_6i_sat(i,1) = etat_6i(i,1);
        etat_6i_sat(i,2) = XSteam('Tsat_p',etat_6i_sat(i,1));
        etat_6i_sat(i,3) = XSteam('hV_p',etat_6i_sat(i,1));
        etat_6i_sat(i,4) = XSteam('sV_p',etat_6i_sat(i,1));
        etat_6i_sat(i,5) = titre_corr(XSteam('x_ph',etat_6i_sat(i,1),etat_6i_sat(i,3)),etat_6i_sat(i,3),etat_6i_sat(i,1));
    end      
end

%%Remplissage des DATA

if n_reheat>0 && n_sout>0
    p = [etat_1(1,1) etat_2(1,1) etat_21(:,1)' etat_22(1,1) etat_3(1,1) ...
        etat_4i(:,1)' etat_5i(:,1)' fliplr(etat_6i(:,1)')...
        etat_6(1,1) etat_7(1,1) etat_7i(:,1)' etat_8(1,1) etat_9(1,1) etat_9i(:,1)'...
        etat_10(1,1) etat_11(1,1)];
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) ...
        etat_4i(:,2)' etat_5i(:,2)' fliplr(etat_6i(:,2)')...
        etat_6(1,2) etat_7(1,2) etat_7i(:,2)' etat_8(1,2) etat_9(1,2) etat_9i(:,2)'...
        etat_10(1,2) etat_11(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) ...
        etat_4i(:,3)' etat_5i(:,3)' fliplr(etat_6i(:,3)')...
        etat_6(1,3) etat_7(1,3) etat_7i(:,3)' etat_8(1,3) etat_9(1,3) etat_9i(:,3)'...
        etat_10(1,3) etat_11(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) ...
        etat_4i(:,4)' etat_5i(:,4)' fliplr(etat_6i(:,4)')...
        etat_6(1,4) etat_7(1,4) etat_7i(:,4)' etat_8(1,4) etat_9(1,4) etat_9i(:,4)'...
        etat_10(1,4) etat_11(1,4)];
    x = [etat_1(1,5) etat_2(1,5) etat_21(:,5)' etat_22(1,5) etat_3(1,5) ...
        etat_4i(:,5)' etat_5i(:,5)' fliplr(etat_6i(:,5)')...
        etat_6(1,5) etat_7(1,5) etat_7i(:,5)' etat_8(1,5) etat_9(1,5) etat_9i(:,5)'...
        etat_10(1,5) etat_11(1,5)];
    e = [etat_1(1,6) etat_2(1,6) etat_21(:,6)' etat_22(1,6) etat_3(1,6) ...
        etat_4i(:,6)' etat_5i(:,6)' fliplr(etat_6i(:,6)')...
        etat_6(1,6) etat_7(1,6) etat_7i(:,6)' etat_8(1,6) etat_9(1,6) etat_9i(:,6)'...
        etat_10(1,6) etat_11(1,6)];
elseif n_reheat==0 && n_sout>0
    p = [etat_1(1,1) etat_2(1,1) etat_21(:,1)' etat_22(1,1) etat_3(1,1) ...
        fliplr(etat_6i(:,1)') etat_6(1,1) etat_7(1,1) etat_7i(:,1)' etat_8(1,1) ...
        etat_9(1,1) etat_9i(:,1)' etat_10(1,1) etat_11(1,1)];
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) ...
        fliplr(etat_6i(:,2)') etat_6(1,2) etat_7(1,2) etat_7i(:,2)' etat_8(1,2)...
        etat_9(1,2)  etat_9i(:,2)' etat_10(1,2) etat_11(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) ...
        fliplr(etat_6i(:,3)') etat_6(1,3) etat_7(1,3) etat_7i(:,3)' etat_8(1,3) etat_9(1,3)...
        etat_9i(:,3)' etat_10(1,3) etat_11(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) ...
        fliplr(etat_6i(:,4)') etat_6(1,4) etat_7(1,4) etat_7i(:,4)' etat_8(1,4) ...
        etat_9(1,4) etat_9i(:,4)' etat_10(1,4) etat_11(1,4)];
    x = [etat_1(1,5) etat_2(1,5) etat_21(:,5)' etat_22(1,5) etat_3(1,5) ...
        fliplr(etat_6i(:,5)') etat_6(1,5) etat_7(1,5) etat_7i(:,5)' etat_8(1,5)...
        etat_9(1,5) etat_9i(:,5)' etat_10(1,5) etat_11(1,5)];
    e = [etat_1(1,6) etat_2(1,6) etat_21(:,6)' etat_22(1,6) etat_3(1,6) ...
        fliplr(etat_6i(:,6)') etat_6(1,6) etat_7(1,6) etat_7i(:,6)' etat_8(1,6) ...
        etat_9(1,6) etat_9i(:,6)' etat_10(1,6) etat_11(1,6)];
elseif n_reheat>0 && n_sout==0
    if supercritique == 1
        etat_21(:,:) = NaN; 
        etat_22(:,:) = NaN; 
    end
    p = [etat_1(1,1) etat_2(1,1) etat_21(:,1)' etat_22(1,1) etat_3(1,1) ...
        etat_4i(:,1)' etat_5i(:,1)' etat_6(1,1)];
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) ...
        etat_4i(:,2)' etat_5i(:,2)' etat_6(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) ...
        etat_4i(:,3)' etat_5i(:,3)' etat_6(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) ...
        etat_4i(:,4)' etat_5i(:,4)' etat_6(1,4)];
    x = [etat_1(1,5) etat_2(1,5) etat_21(:,5)' etat_22(1,5) etat_3(1,5) ...
        etat_4i(:,5)' etat_5i(:,5)' etat_6(1,5)];
    e = [etat_1(1,6) etat_2(1,6) etat_21(:,6)' etat_22(1,6) etat_3(1,6) ...
        etat_4i(:,6)' etat_5i(:,6)' etat_6(1,6)];
else
    p = [etat_1(1,1) etat_2(1,1) etat_21(:,1)' etat_22(1,1) etat_3(1,1) etat_6(1,1)];
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) etat_6(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) etat_6(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) etat_6(1,4)];
    x = [etat_1(1,5) etat_2(1,5) etat_21(:,5)' etat_22(1,5) etat_3(1,5) etat_6(1,5)];
    e = [etat_1(1,6) etat_2(1,6) etat_21(:,6)' etat_22(1,6) etat_3(1,6) etat_6(1,6)];
end

DAT=[T;p;h;s;e;x];

if display==1

% correction des etats 2' et 2'' s'ils n'existent pas, pour les plots
if supercritique == 1
	for i=1:6
        etat_22(1,i) = NaN;
    end
	for i=1:N_discr
        etat_21(i,1) = NaN;
        etat_21(i,2) = T2_to_T3(i);
        etat_21(i,3) = h2_to_h3(i);
        etat_21(i,4) = s2_to_s3(i);
        etat_21(i,5) = NaN;
        etat_21(i,6) = NaN;
    end
end  
if n_reheat>0 && n_sout==0
        reheat_T = [];
        reheat_h = [];
        reheat_s = [];    
    for i=1:n_reheat
        reheat_T = [reheat_T,[etat_4i(i,2)' etat_5i(i,2)']];
        reheat_h = [reheat_h,[etat_4i(i,3)' etat_5i(i,3)']];
        reheat_s = [reheat_s,[etat_4i(i,4)' etat_5i(i,4)']];
    end
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) ...
        reheat_T etat_6(1,2) etat_7(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) ...
        reheat_h etat_6(1,3) etat_7(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) ...
        reheat_s etat_6(1,4) etat_7(1,4)];
elseif n_reheat==0 && n_sout>0
    T6i = fliplr(etat_6i(:,2)');
    h6i = fliplr(etat_6i(:,3)');
    s6i = fliplr(etat_6i(:,4)');
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) ...
         T6i(:)' etat_6(1,2) etat_7(1,2) etat_7i(:,2)' etat_8(1,2)...
        etat_9(1,2) etat_9i(:,2)'];%etat_10(1,2) etat_11(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) ...
         h6i(:)' etat_6(1,3) etat_7(1,3) etat_7i(:,3)' etat_8(1,3) etat_9(1,3) etat_9i(:,3)'];%...
        %etat_10(1,3) etat_11(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) ...
         s6i(:)' etat_6(1,4) etat_7(1,4) etat_7i(:,4)' etat_8(1,4) ...
        etat_9(1,4) etat_9i(:,4)'];% etat_10(1,4) etat_11(1,4)];
elseif n_reheat>0 && n_sout>0
        reheat_T = []; reheat_h = []; reheat_s = []; 
    for i=1:n_reheat
        reheat_T = [reheat_T,[etat_4i(i,2)' etat_5i(i,2)']];
        reheat_h = [reheat_h,[etat_4i(i,3)' etat_5i(i,3)']];
        reheat_s = [reheat_s,[etat_4i(i,4)' etat_5i(i,4)']];
    end
    T6i = fliplr(etat_6i(1:(n_sout-n_reheat),2)');
    h6i = fliplr(etat_6i(1:(n_sout-n_reheat),3)');
    s6i = fliplr(etat_6i(1:(n_sout-n_reheat),4)');
    
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) ...
        reheat_T T6i(:)' etat_6(1,2) etat_7(1,2) etat_8(1,2) etat_9(1,2) etat_9i(:,2)'...
        etat_10(1,2) etat_11(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) ...
        reheat_h h6i(:)' etat_6(1,3) etat_7(1,3) etat_8(1,3) etat_9(1,3) etat_9i(:,3)'...
        etat_10(1,3) etat_11(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) ...
        reheat_s s6i(:)' etat_6(1,4) etat_7(1,4)  etat_8(1,4) etat_9(1,4) etat_9i(:,4)'...
        etat_10(1,4) etat_11(1,4)];
else
    T = [etat_1(1,2) etat_2(1,2) etat_21(:,2)' etat_22(1,2) etat_3(1,2) etat_6(1,2) etat_7(1,2)];
    h = [etat_1(1,3) etat_2(1,3) etat_21(:,3)' etat_22(1,3) etat_3(1,3) etat_6(1,3) etat_7(1,3)];
    s = [etat_1(1,4) etat_2(1,4) etat_21(:,4)' etat_22(1,4) etat_3(1,4) etat_6(1,4) etat_7(1,4)];
end

FIG(1) = figure;
T_int=linspace(0,400,1000);
        for i=1:length(T_int)
            S_c_L(i)=XSteam('sL_T',T_int(i));
            S_c_V(i)=XSteam('sV_T',T_int(i));
        end
        S_c=[S_c_L,S_c_V];
        T_c=[T_int,T_int];
        plot(S_c,T_c)
        xlabel('S [kJ/kg K]')
        ylabel('T [degC]')
        hold on
plot(s,T,'k')
hold on
for i=1:n_sout
	if etat_6i(i,5)<1
        plot([etat_6i(i,4),etat_7i(i,4)],[etat_6i(i,2),etat_7i(i,2)],'--r');    
        hold on;
    else
        plot([etat_6i(i,4),etat_6i_sat(i,4),etat_7i(i,4)],[etat_6i(i,2),etat_6i_sat(i,2),etat_7i(i,2)],'--r');
        hold on;
	end
end
plot(DAT(4,:),DAT(1,:),'.','markersize',20)

FIG(2) = figure;
T_int=linspace(0,400,1000);
        for i=1:length(T_int)
            S_c_L(i)=XSteam('sL_T',T_int(i));
            S_c_V(i)=XSteam('sV_T',T_int(i));
            h_c_L(i)=XSteam('hL_T',T_int(i));
            h_c_V(i)=XSteam('hV_T',T_int(i));
        end
        S_c=[S_c_L,S_c_V];
        H_c=[h_c_L,h_c_V];
        plot(S_c,H_c)
        xlabel('S [kJ/kg K]')
        ylabel('H [degC]')
        hold on
plot(s,h,'k')
hold on
for i=1:n_sout
	if etat_6i(i,5)<1
        plot([etat_6i(i,4),etat_7i(i,4)],[etat_6i(i,3),etat_7i(i,3)],'--r');    
        hold on;
    else
        plot([etat_6i(i,4),etat_6i_sat(i,4),etat_7i(i,4)],[etat_6i(i,3),etat_6i_sat(i,3),etat_7i(i,3)],'--r');
        hold on;
	end
end
plot(DAT(4,:),DAT(3,:),'.','markersize',20)

FIG(3) = figure;
pie([P_e,DATEN])
legend(sprintf('Puissance effective: %10.1fMW',P_e*1e-6),...
       sprintf('Pertes au generateur de vapeur: %10.1fMW',P_gen*1e-6),...
       sprintf('Pertes mecaniques: %10.1fMW',P_mecen*1e-6),...
       sprintf('Pertes au condenseur: %10.1fMW',P_conden*1e-6),'Location','northeastoutside')
title(sprintf('Puissance energetique primaire: %10.1fMW',P_toten*1e-6))

FIG(4) = figure;
pie([P_e,DATEX(1),DATEX(3:end)])
title(sprintf('Puissance exergetique primaire: %10.1fMW',P_totex*1e-6))
legend(sprintf('Puissance effective: %10.1fMW',P_e*1e-6),...
       sprintf('Irreversibilites turbine et pompes: %10.1fMW',P_rotex*1e-6),...
       sprintf('Irreversibilite combustion: %10.1fMW',P_combex*1e-6),...
       sprintf('Pertes condenseur: %10.1fMW',P_condex*1e-6),...
       sprintf('Pertes cheminee: %10.1fMW',P_chemex*1e-6),...
       sprintf('Pertes transfert de chaleur generateur de vapeur: %10.1fMW',P_transex*1e-6),'Location','northeastoutside')

else
    FIG=[];
end
end

