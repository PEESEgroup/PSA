function flux_w = WENO(flux_c, FlowDir)
%WENO: Apply Weighted Essentially NonOscillatory scheme
%   Receive fluxes at finite volume centers and direction of flow (upwind -
%   downwind), the flow direction is given by the velocity or pressure drop
%   If v>=0 or dP<0, then an upwind scheme is applied. If v<0 or dP>0, then
%   an downwind scheme is applied. The function return the fluxes at finite
%   volume walls. upwind means that the computation is done from left to 
%   right on the domain, or in terms of the PSA column "co-current flow 
%   regarding the feed inlet". downwind means that the computation is done 
%   from right to left on the domain, or in terms of the PSA column "counter
%   -current flow regarding the feed inlet"
%   
%   Input:
%       flux_c : flux at the finite volume centers
%       FlowDir: direction of flow. OPTIONS: upwind and downwind
%   
%   Output:
%       flux_w : flux at the edges or walls of the finite volumes
%   
%%  
%   For co-current flow (upwind) the fluxes at the walls of finite volumes
%   are calculated as follows:
%   
%%  
%   $$ f_{j+0.5}=\frac{\alpha_{0,j}}{\alpha_{0,j}+\alpha_{1,j}} 
%   \Big[\frac{1}{2}(f_{j}+f_{j+1})\Big] + \frac{\alpha_{1,j}}{\alpha_{0,j}+
%   \alpha_{1,j}}\Big[\frac{3}{2}f_{j}-\frac{1}{2}f_{j-1}\Big] $$
%   
%   $$ \alpha_{0,j}= \frac{2/3}{(f_{j+1}-f_{j}+\delta)^4} $$
%   
%   $$ \alpha_{1,j}= \frac{1/3}{(f_{j}-f_{j-1}+\delta)^4} $$
%   
%%  
%   For counter-current flow (downwind) the fluxes at the walls of finite 
%   volumes are calculated as follows:
%   
%%  
%   $$ f_{j+0.5}=\frac{\alpha_{0,j}}{\alpha_{0,j}+\alpha_{1,j}} 
%   \Big[\frac{1}{2}(f_{j}+f_{j+1})\Big] + \frac{\alpha_{1,j}}{\alpha_{0,j}+
%   \alpha_{1,j}}\Big[\frac{3}{2}f_{j+1}-\frac{1}{2}f_{j+2}\Big] $$
%   
%   $$ \alpha_{0,j}= \frac{2/3}{(f_{j}-f_{j+1}+\delta)^4} $$
%   
%   $$ \alpha_{1,j}= \frac{1/3}{(f_{j+1}-f_{j+2}+\delta)^4} $$
%   
%%  
    oo     = 10^-10              ;
    [N, m] = size(flux_c)        ;
    N      = N-2                 ;
    flux_w = zeros(N+1, m)       ;
    alpha0 = zeros(size(flux_c)) ;
    alpha1 = zeros(size(flux_c)) ;
    
    % Fluxes at boundaries of the domain
    flux_w(1, :)   = flux_c(1, :)   ;
    flux_w(N+1, :) = flux_c(N+2, :) ;
    
    if strcmpi(FlowDir, 'upwind') == 1
        
        alpha0(2:N, :) =(2/3)./((flux_c(3:N+1, :)-flux_c(2:N, :)+oo).^4) ;
        alpha1(3:N, :) =(1/3)./((flux_c(3:N, :)-flux_c(2:N-1, :)+oo).^4) ;
        alpha1(2, :)   =(1/3)./((2*(flux_c(2, :)-flux_c(1, :))+oo).^4)   ;
        
        flux_w(3:N, :) = (alpha0(3:N, :)./(alpha0(3:N, :)+alpha1(3:N, :)))...      
                       .*((flux_c(3:N, :)+flux_c(4:N+1, :))./2)+(alpha1(3:N, :)... 
                       ./(alpha0(3:N, :)+alpha1(3:N, :))).*(1.5*flux_c(3:N, :)...  
                        -.5*flux_c(2:N-1, :))                                     ;
        
        flux_w(2, :)   = (alpha0(2, :)./(alpha0(2, :)+alpha1(2, :)))...    
                       .*((flux_c(2, :)+flux_c(3, :))./2)+(alpha1(2, :)... 
                       ./(alpha0(2, :)+alpha1(2, :))).*(2*flux_c(2, :)...  
                         -flux_c(1, :))                                   ;
        
    elseif strcmpi(FlowDir, 'downwind') == 1
        
        alpha0(2:N, :)   = (2/3)./((flux_c(2:N, :)-flux_c(3:N+1, :)+oo).^4)   ;
        alpha1(2:N-1, :) = (1/3)./((flux_c(3:N, :)-flux_c(4:N+1, :)+oo).^4)   ;
        alpha1(N, :)     = (1/3)./((2*(flux_c(N+1, :)-flux_c(N+2, :))+oo).^4) ;
        
        flux_w(2:N-1, :) = (alpha0(2:N-1, :)./(alpha0(2:N-1, :)+alpha1(2:N-1, :)))...   
                         .*((flux_c(2:N-1, :)+flux_c(3:N, :))./2)+(alpha1(2:N-1, :)...  
                         ./(alpha0(2:N-1, :)+alpha1(2:N-1, :))).*(1.5*flux_c(3:N, :)... 
                          -.5*flux_c(4:N+1, :))                                        ;
                      
        flux_w(N, :)     = (alpha0(N, :)./(alpha0(N, :)+alpha1(N, :)))...      
                         .*((flux_c(N, :)+flux_c(N+1, :))./2)+(alpha1(N, :)... 
                         ./(alpha0(N, :)+alpha1(N, :))).*(2*flux_c(N+1, :)...  
                           -flux_c(N+2, :))                                   ;
    else
        error('Please specify the direction of flow. OPTIONS: upwind and downwind')
    end 
%   
end 