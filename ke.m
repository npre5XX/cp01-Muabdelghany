% NPRE 555 Computer Project-1 "Monte Carlo Code of Neutron Transport"
function keff=ke
% clc
% clear all

npart=1000;             % number of particles

for j=1:npart

x(j)=100*rand();        % the particle is born in a random position in x and y
y(j)=rand();

lost=0;                 % parameter specifiy whether the particle is lost from the slab or not? 
history(j)=0;           % counter for the number of histories before the particle is absorbed or leaked out

    while lost == 0
                                         % for the left half of the slab
    if x(j)<=50                          
        seg_a=0.12;                      % cross-setion for absrption
        seg_s=0.05;                      % corss-section for scattering
        nwseg_f=0.15;                    % nuw* cross-section of fission
        kinf = nwseg_f/seg_a;            % K_infinity of neutrons produced per fission
        seg_t = seg_a + seg_s;
                                         % for the right half of the slab
    else
        seg_a=0.10;                      % cross-setion for absrption
        seg_s=0.05;                      % corss-section for scattering
        nwseg_f=0.12;                    % nuw* cross-section of fission
        kinf = nwseg_f/seg_a;            % K_infinity of neutrons produced per fission
        seg_t = seg_a + seg_s;
    end 
    
        r= -(1/seg_t)*log(1-rand());     % sampling the traveled distance by (r)
        theta=2*pi*rand();               % generating a random number for theta between [0,2pi]

        dx=r*cos(theta);
        dy=r*sin(theta);
                                        % updating the traveled distance for the particle (j)
        x(j)=x(j)+dx;
        y(j)=y(j)+dy;

        history(j)=history(j)+1;
        
            if x(j)<=50                          
                seg_a=0.12;                      % cross-setion for absrption
                seg_s=0.05;                      
                nwseg_f=0.15;                   
                kinf = nwseg_f/seg_a;           
                seg_t = seg_a + seg_s;
                                             % for the right half of the slab
            else
                seg_a=0.10;                      
                seg_s=0.05;                      
                nwseg_f=0.12;                    
                kinf = nwseg_f/seg_a;            
                seg_t = seg_a + seg_s;
            end 
    
        if x(j)>0 && x(j)<100 && rand()< seg_a/seg_t   % No leakage and absorbtion occured
            xn=1+fix(x(j));                            % defining a discretization for x (delta(x)= 1 Cm) 
            epsi(j,xn)=0;
            k(j)=fix(kinf+rand());
            epsi(j,xn)=k(j);                           % calculating the flux due to each fission in this space mesh
            lost=1;
        elseif  x(j)<=0 || x(j)>=100                   % if the neutron is leaked from the system 
                lost=1;
        end
        
% for small number of particles we can plot their position, but for large
% number it will take a very long time
%         plot(x,y,'r.')
%         hold on
%         xlim([0 100])
    
    end

end

alpha=0;            % alpha is the integrated flux over the volume of the slab
% calculating the flux at each space mesh; phi(xn) & the integrated flux alpha
for xn=1:100
    phi(xn)=0;
for j=1:npart
    phi(xn)=phi(xn)+epsi(j,xn);
end
    alpha=alpha+phi(xn);
end

% plotting the flux with x
figure(2)
z=linspace(1,100);
plot(z,phi(z));

keff=alpha/npart        % calculating the multiplication factor (k_effective)