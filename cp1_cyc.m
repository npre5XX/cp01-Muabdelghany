% NPRE 555 Computer Project-1 "Monte Carlo Code of Neutron Transport"
clc
clear all
cyc= 5;         % number of cycles
for i=1:cyc
    m(i)=ke;          % function that returns the value of K_effective for each cycle
end
m                     % K_effective array
keffective= mean(m)   % Average multiplication factor over all the cycles