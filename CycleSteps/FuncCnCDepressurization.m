function derivatives = FuncCnCDepressurization(~, state_vars, params, isotherm_params)
%#codegen
%CnC_Depressurization: Calculate the change in state variables for CnC depres step
%   When coupled with MATLAB ode solvers, (ode15s), solve the system of PDEs
%   for the PSA step. The system of PDEs is solved by using the finite volume
%   method (FVM) where the spatial domain is discretized into N volumes. 
%   Each of these volumes is assumed to have a uniform value for the state 
%   variables. Weighted Essentially Non-oscillatory (WENO) is used to improve
%   calculations by calculating the value of the state variables at the walls
%   of the finite volumes. 
%   
%   Input:
%       ~: This is the time variable. Since time is not used to calculate
%          the derivatives but a place holder is required for the ode solver,
%          this is used
%   
%       state_vars: This is the dimensionless state variables at time being
%          evaluated. The ordering is as followed [P_1, P_2,...,P_N+2, y_1,
%          ...,y_N+2, q_CO2_1,...,q_CO2_N+2, q_N2_1,...,q_N2_N+2, T_1,...,
%          T_N+2]
%   
%       params: All parameters needed for the simulation. Provided in the 
%          ProcessInputParameters function file
%   
%       isotherm_params: All parameters for the isotherm. Provided in the 
%          ProcessInputParameters function file
%   
%   Output:
%       derivatives: These are the temporal derivatives of the state
%       variables. The ordering is the same as the state variables
%   
%% Retrieve process parameters
    N				=	params(1)	;
    deltaU_1	    =	params(2)	;
    deltaU_2	    =	params(3)	;
    ro_s			=	params(4)	;
    T_0				=	params(5)	;
    epsilon			=	params(6)	;
    r_p				=	params(7)	;
    mu				=	params(8)	;
    R				=	params(9)	;
    v_0				=	params(10)	;
    q_s0			=	params(11)	;
    C_pg			=	params(12)	;
    C_pa			=	params(13)	;
    C_ps			=	params(14)	;
    D_m				=	params(15)	;
    K_z				=	params(16)	;
    P_0				=	params(17)	;
    L				=	params(18)	;
    MW_CO2			=	params(19)	;
    MW_N2			=	params(20)	;
    k_1_LDF		    =	params(21)	;
    k_2_LDF		    =	params(22)	;
    tau				=	params(24)	;
    P_l				=	params(25)	;
%   
%% Initialize state variables
    P  = zeros(N+2, 1) ;
    y  = zeros(N+2, 1) ;
    x1 = zeros(N+2, 1) ;
    x2 = zeros(N+2, 1) ;
    T  = zeros(N+2, 1) ;
    
    P(1:N+2)  = state_vars(1:N+2)               ;
    y(1:N+2)  = max(state_vars(N+3:2*N+4), 0)   ;
    x1(1:N+2) = max(state_vars(2*N+5:3*N+6), 0) ;
    x2(1:N+2) = state_vars(3*N+7:4*N+8)         ;
    T(1:N+2)  = state_vars(4*N+9:5*N+10)        ;
%   
%% Initialize all variables used in the function
    % Temporal derivatives
    derivatives = zeros(5*N+10, 1) ;
    dPdt        = zeros(N+2, 1)    ;
    dPdt1       = zeros(N+2, 1)    ;
    dPdt2       = zeros(N+2, 1)    ;
    dPdt3       = zeros(N+2, 1)    ;
    dydt        = zeros(N+2, 1)    ;
    dydt1       = zeros(N+2, 1)    ;
    dydt2       = zeros(N+2, 1)    ;
    dydt3       = zeros(N+2, 1)    ;
    dx1dt       = zeros(N+2, 1)    ;
    dx2dt       = zeros(N+2, 1)    ;
    dTdt        = zeros(N+2, 1)    ;
    dTdt1       = zeros(N+2, 1)    ;
    dTdt2       = zeros(N+2, 1)    ;
    dTdt3       = zeros(N+2, 1)    ;
    % Spatial derivatives
    dpdz        = zeros(N+2, 1)    ;
    dpdzh       = zeros(N+1, 1)    ;
    dydz        = zeros(N+2, 1)    ;
    d2ydz2      = zeros(N+2, 1)    ;
    dTdz        = zeros(N+2, 1)    ;
    d2Tdz2      = zeros(N+2, 1)    ;
%   
%%  
%   The following estimations are parameters that are used. The axial
%   dispersion coefficient is calculated using the following equation
%   from Ruthvan "Pressure Swing Adsorption"
%   
%%  
%   $$ D_{L}= 0.7 D_{m} + r_{p} \cdot v_{0} $$
%   
%%  
%   The Concentration of gas is calculated using the ideal gas law
%   
%%  
%   $$ C_{g}=\frac{\bar{P}}{R\bar{T}} \cdot \frac{P_{0}}{T_{0}} $$
%   
%% Calculate all parameters used
    dz   = 1/N                                ;
    D_l  = 0.7*D_m + v_0*r_p                  ;
    Pe   = v_0*L/D_l                          ;
    phi  = R*T_0*q_s0*(1-epsilon)/epsilon/P_0 ;
    ro_g = P(1:N+2).*P_0/R./T(1:N+2)/T_0      ;
%   
%% Boundary Conditions
%   The boundary condtions for the CnC depressurization, corresponding to the
%   heavy product end of the column being opened and the the light product 
%   end of the column being closed. For the heavy product end, the following 
%   boundary condtions are used:
%   
%%  
%   $$ \frac{\partial \bar{P}}{\partial \tau} = \lambda(\bar{P_L}-\bar{P}) \quad | Z=0^- $$
%   
%   $$ \frac{\partial y}{\partial \tau} = 0 \quad | Z=0^- $$
%   
%   $$ \frac{\partial \bar{T}}{\partial \tau} = 0 \quad | Z=0^- $$
%   
%%  
    y(1) = y(2)     ;
    T(1) = T(2)     ;
    if P(1) > P(2)
        P(1) = P(2) ;
    end
%   
%%  
%   For the light product end, the following boundary condtions are used:
%   
%%  
%   $$ \frac{\partial \bar{P}}{\partial \tau} = 0 \quad | Z=1^+ $$
%   
%   $$ \frac{\partial y}{\partial \tau} = 0 \quad | Z=1^+ $$
%   
%   $$ \frac{\partial \bar{T}}{\partial \tau} = 0 \quad | Z=1^+ $$
%   
%%  
    y(N+2) = y(N+1) ;
    T(N+2) = T(N+1) ;
    P(N+2) = P(N+1) ;
%   
%% Spatial Derivative Calculations
%   
%   *1st Derivatives*
%   
%   For the first derivative, the value at each volume is estimated based on
%   the values at the walls of the volumes. These values are estimated using
%   the Weighted Essentially NonOscillatory (WENO) scheme.
%   
%%  
%   $$ \frac{\partial f_j}{\partial Z}=\frac{f_{j+0.5}-f_{j-0.5}}{\Delta Z} $$
%   
%%  
%   Pressure: at the center of the volumes and the walls
    
    Ph          = WENO(P, 'downwind')    ;
    
    dpdz(2:N+1) = (Ph(2:N+1)-Ph(1:N))/dz ;
    dpdzh(2:N)  = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(1)    =  2*(P(2)-P(1))/dz      ;
    dpdzh(N+1)  =  2*(P(N+2)-P(N+1))/dz  ;
%   
%%  
%   Mole Fraction: at the center of the volumes
    
    yh          = WENO(y, 'downwind')    ;
    
    dydz(2:N+1) = (yh(2:N+1)-yh(1:N))/dz ;
%   
%%  
%   Temperature: at the center of the volumes
    
    Th          = WENO(T, 'downwind')    ;
    
    dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;
%   
%%  
%   *2nd Derivatives*
%   
%   The second derivatives are calculated based off of the values of the
%   nodes. It is only necessary to calculate the 2nd derivative of the
%   temperature and mole fraction. It is noted that at the ends of the
%   column, no diffusion/conduction is occuring, so these are set to 0 in
%   the value of the derivative.
%   
%%  
%   $$ \frac{{\partial}^2 f_{j}}{\partial {Z}^2} = \frac{f_{j+1}+f_{j-1}-
%       2f_{j}}{{\Delta Z}^2} $$
%   
%%  
%   Mole Fraction
    d2ydz2(3:N) = (y(4:N+1)+y(2:N-1)-2*y(3:N))/dz/dz ;
    d2ydz2(2)   = (y(3)-y(2))/dz/dz                  ;
    d2ydz2(N+1) = (y(N)-y(N+1))/dz/dz                ;
%   
%%  
%   Temperature
    d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;
%   
%% Velocity Calculations
%   Calculates the interstitial velocity of the gas at the walls of the
%   volumes based off of the pressure gradients. NOTE: bear in mind that
%   the velocity in the Ergun's equation is the superficial velocity, and
%   not the interstitial, that's why in the denominator in the viscous term
%   it is found epsilon^2 and not epsilon^3, the same for kinetic terms,
%   where it is found epsilon and not epsilon^3. Superficial velocity is 
%   equal to interstitial velocity times the void fraction
%   
%%  
%   $$ U = v \cdot \varepsilon $$
%   
%   $$ - \frac{\partial \bar{P}}{\partial Z}\frac{P_0}{L} =
%   \frac{150\mu(1-\varepsilon)^2}{4r_{p}\varepsilon^2}\bar{v}v_0 +
%   \frac{1.75(1-\varepsilon)}{2r_{p}\varepsilon}
%   (\sum_{i}y_{i}MW_{i}C_{g})\bar{v}^2{v_{0}}^2 $$
%   
%%  
    ro_gh          = (P_0/R/T_0)*Ph(1:N+1)./Th(1:N+1)               ;
    
    viscous_term   = 150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2         ;
    kinetic_term_h = (ro_gh.*(MW_N2+(MW_CO2-MW_N2).*yh)).*(1.75*(1-epsilon)...
                     /2/r_p/epsilon)                                         ;
                    
    % Velocities at walls of volumes 
    vh = -sign(dpdzh).*(-viscous_term+(abs(viscous_term^2+4*kinetic_term_h... 
               .*abs(dpdzh)*P_0/L)).^(.5))/2./kinetic_term_h/v_0             ;
%   
%% Temporal Derivatives
%   
%%  
%   *1) Adsorbed Mass Balance (molar loading for components 1 and 2)*
%       Linear Driving Force (LDF) is used to calculate the uptake rate of
%       the gas. The LDF mass transfer coefficients are assumed to be
%       constant over the pressure and temperature range of importance.
%   
%%  
%   $$ \frac{\partial x_{i}}{\partial \tau}q_{s0}= k_{i}({q_{i}}^*-q_{i}) $$
%   
%%  
%   1.1) Calculate equilibrium molar loading
    q   = Isotherm(y, P*P_0, T*T_0, isotherm_params) ;
    q_1 = q(:, 1)*ro_s                               ;
    q_2 = q(:, 2)*ro_s                               ;
%   
%%  
%   1.2) Calculate the LDF parameter
    k_1 = k_1_LDF*L/v_0 ;
    k_2 = k_2_LDF*L/v_0 ;
%   
%%  
%   1.3) Calculate the temporal derivative 
    dx1dt(2:N+1) = k_1*(q_1(2:N+1)/q_s0 - x1(2:N+1)) ;
    dx2dt(2:N+1) = k_2*(q_2(2:N+1)/q_s0 - x2(2:N+1)) ;
%   
%%  
%   *2) Column Energy Balance (Column Temperature)*
%       For the energy balance in the column, it is assumed that:
%       * There is thermal equilibrium between the solid and gas phase
%       * Conduction occurs through both the solid and gas phase 
%       * Advection only occurs in the gas phase 
%   
%%  
%   $$ \big[ \varepsilon C_{g}C_{p,g}+(1-\varepsilon)(C_{p,s}\rho_{s}+
%       C_{p,a}q_{s0})\big]\frac{\partial \bar{T}}{\partial \tau}= 
%       \frac{K_z}{v_0L} \frac{\partial^2 \bar{T}}{\partial Z^2}
%     - \varepsilon C_{g}C_{p,g} \bar{v}\frac{\partial \bar{T}}{\partial Z}
%     + \sum_{i} (1-\varepsilon)(-\Delta H_{i})\frac{q_{s0}}{T_{0}}\frac{
%       \partial \bar{x_{i}}}{\partial \tau} $$
%   
%%  
    % [J/m^3/K]
    sink_term = ((1-epsilon)*(ro_s*C_ps+q_s0*C_pa)+(epsilon.*ro_g(2:N+1).*C_pg)) ;
%   
%%  
%   2.1) Calculate the temperature change due to conduction in the column 
%        (both solid and gas phase)
%   
%%  
%   $$ \frac{K_z}{v_0L} \frac{\partial^2 \bar{T}}{\partial Z^2} $$
%   
%%  
    transfer_term = K_z./v_0./L                             ;
    dTdt1(2:N+1)  = transfer_term.*d2Tdz2(2:N+1)./sink_term ;
%   
%%  
%   2.2) Calculate the temperature change due to advection
%   
%%  
%   $$ -\varepsilon C_{g}C_{p,g} \bar{v}\frac{\partial \bar{T}}{\partial Z} $$
%   
%%  
    PvT          = Ph(1:N+1).*vh(1:N+1)./Th(1:N+1) ;
    Pv           = Ph(1:N+1).*vh(1:N+1)            ;
    dTdt2(2:N+1) = -epsilon.*C_pg.*P_0./R./T_0.*((Pv(2:N+1)-Pv(1:N))- ... 
                    T(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz./sink_term      ;
%   
%%  
%   2.3) Calculate the temperature change due to adsorption/desoprtion enthalpy
%   
%%  
%   $$ \sum_{i} (1-\varepsilon)(-\Delta
%   H_{i})\frac{q_{s0}}{T_{0}}\frac{\partial \bar{x_{i}}}{\partial \tau} $$
%   
%   $$ \Delta H_{i} = \Delta U_i - R\bar{T}T_{0} $$
%   
%%  
    generation_term_1 = (1-epsilon).*q_s0.*(-(deltaU_1-R*T(2:N+1)*T_0))./T_0 ;
    generation_term_2 = (1-epsilon).*q_s0.*(-(deltaU_2-R*T(2:N+1)*T_0))./T_0 ;
    
    dTdt3(2:N+1)      = (generation_term_1.*dx1dt(2:N+1)+...
                         generation_term_2.*dx2dt(2:N+1))./sink_term  ;
%   
%%  
%   2.4) Total sum of all temperature derivatives 
    dTdt(2:N+1) = dTdt1(2:N+1) + dTdt2(2:N+1) + dTdt3(2:N+1) ;
%   
%%  
%   *3) Total mass balance*
%   
%%  
%   $$ \frac{\partial \bar{P}}{\partial \tau} = -\frac{\partial( \bar{v}
%      \bar{P}/\bar{T})}{\partial Z} - \Psi \bar{T} \sum_{i}\frac{\partial 
%      \bar{x_{i}}}{\partial \tau} + \frac{\bar{P}}{\bar{T}}\frac{\partial 
%      \bar{T}}{\partial \tau} $$
%   
%%  
%   3.1) Calculate the change in pressure due to advection
    dPdt1(2:N+1) = -T(2:N+1).*(PvT(2:N+1)-PvT(1:N))./dz  ;
%   
%%  
%   3.2) Calculate the change in pressure due to adsorption/desorption
    dPdt2(2:N+1) = -phi*T(2:N+1).*(dx1dt(2:N+1)+dx2dt(2:N+1)) ;
%   
%%  
%   3.3) Calculate the change in pressure due to temperature changes
    dPdt3(2:N+1) = P(2:N+1).*dTdt(2:N+1)./T(2:N+1)            ;
%   
%%  
%   3.4) Total sum of all presure changes
    dPdt(2:N+1) = dPdt1(2:N+1) + dPdt2(2:N+1) + dPdt3(2:N+1)  ;
%   
%%  
%   *4) Component Mass Balance (Based on Mole Fraction)*
%   
%%  
%   $$ \frac{\partial y}{\partial \tau} = \frac{1}{Pe} \big(\frac{{
%      \partial}^2 y}{\partial {Z}^2}+\frac{1}{\bar{P}}\frac{\partial
%      \bar{P}}{\partial Z}\frac{\partial y}{\partial Z}-\frac{1}{\bar{T}}
%      \frac{\partial \bar{T}}{\partial Z}\frac{\partial y}{\partial Z}
%      \big)-\bar{v}\frac{\partial y}{\partial Z}+\frac{\Psi \bar{T}}{
%      \bar{P}} \big((y-1)\frac{\partial \bar{x_{1}}}{\partial \tau}+y
%      \frac{\partial\bar{x_{2}}}{\partial \tau}\big) $$
%   
%%  
%   4.1) Calculate the change in mole fraction due to diffusion
%   
%%  
%   $$ \frac{1}{Pe} \big(\frac{{\partial}^2 y}{\partial {Z}^2}+\frac{1}{
%      \bar{P}}\frac{\partial \bar{P}}{\partial Z}\frac{\partial y}{
%      \partial Z}-\frac{1}{\bar{T}}\frac{\partial \bar{T}}{\partial Z}
%      \frac{\partial y}{\partial Z} \big) $$
%   
%%  
    dydt1(2:N+1) = (1/Pe)*(d2ydz2(2:N+1)+(dydz(2:N+1).*dpdz(2:N+1)./P(2:N+1))... 
                  -(dydz(2:N+1).*dTdz(2:N+1)./T(2:N+1)))                        ;
%   
%%  
%   4.2) Calculate the change in mole fraction due to advection
%   
%%  
%   $$ -\bar{v}\frac{\partial y}{\partial Z} $$
%   
%%  
    ypvt         = yh(1:N+1).*Ph(1:N+1).*vh(1:N+1)./Th(1:N+1)        ;
    dydt2(2:N+1) = -(T(2:N+1)./P(2:N+1)).*((ypvt(2:N+1)-ypvt(1:N))... 
                   -y(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz             ;
%   
%%  
%   4.3) Calculate the change in mole fraction due to adsorption/desorption
%   
%%  
% $$ \frac{\Psi \bar{T}}{\bar{P}} \big((y-1)\frac{\partial \bar{x_{1}}}{
%    \partial \tau}+y\frac{\partial \bar{x_{2}}}{\partial \tau}\big) $$
%   
%%  
    dydt3(2:N+1) = (phi*T(2:N+1)./P(2:N+1)).*((y(2:N+1)-1).*dx1dt(2:N+1)... 
                  + y(2:N+1).*dx2dt(2:N+1))                                ;
%   
%%  
%   4.4) Total sum of all mole fraction changes
    dydt(2:N+1) = dydt1(2:N+1) + dydt2(2:N+1) + dydt3(2:N+1) ;
%   
%%  Boundary Derivatives
    %dPdt(1)    = tau*(P_l/P_0-P(1))       ;
	dPdt(1)    = tau*L/v_0*(P_l/P_0-P(1)) ;
    dPdt(N+2)  = dPdt(N+1)                ;
    dydt(1)    = dydt(2)                  ;
    dydt(N+2)  = dydt(N+1)                ;
    dx1dt(1)   = 0                        ;
    dx2dt(1)   = 0                        ;
    dx1dt(N+2) = 0                        ;
    dx2dt(N+2) = 0                        ;
    dTdt(1)    = dTdt(2)                  ;
    dTdt(N+2)  = dTdt(N+1)                ;
%   
%%  Export derivatives to output
    derivatives(1:N+2)        = dPdt(1:N+2)  ;
    derivatives(N+3:2*N+4)    = dydt(1:N+2)  ;
    derivatives(2*N+5:3*N+6)  = dx1dt(1:N+2) ;
    derivatives(3*N+7:4*N+8)  = dx2dt(1:N+2) ;
    derivatives(4*N+9:5*N+10) = dTdt(1:N+2)  ;
%   
end 