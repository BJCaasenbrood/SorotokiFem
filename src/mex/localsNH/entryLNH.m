eNode = [7, 10, 9, 8];
eDof = [13, 14, 19, 20, 17, 18, 15, 16];
dV   = 0;
Rb   = 1;
Dim  = 2;
W = ones(4,1)*0.5;
Utmp = zeros(44,1);
dUtmp = zeros(44,1);
Rho = 9.7000e-10;
Zeta = 0.1;
Grav = [0;-9810];
Mu = 0.0038;
Lambda = 0.0058;


for ii = 1:1e4
tic;
[Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,Tke,Re,Ue] = ...
        LocalsNHFast(eNode,eDof,dV,Rb,...
        Dim,...   % Fem.Dim
        Node0,... % Node0
        N,...     % shpfnc for Fem.ShapeFnc{nn}
        dNdxi,... % derive shpfnc for Fem.ShapeFnc{nn}
        W,...     % weights shapefnc
        Utmp,...  % current displacements
        dUtmp,... % current velocities
        Rho,...   % density
        Zeta,...  % dampings
        Grav,...  % gravity
        Mu,...     
        Lambda...
        );
T(ii) = toc;
end

disp(['Average elapsed time is ', num2str(mean(T),6), '(s)']);
    
function y = Node0    
y = [-0.0000    2.9815
   -0.0000    0.0000
    4.0839    0.0000
    4.0839    2.9815
    8.1596    0.0000
    8.1148    3.0000
   12.2196    0.0000
   12.2196    2.9365
   16.2581    3.0000
   16.2581    0.0000
   20.2713    3.0000
   20.2713    0.0000
   24.2580    0.0000
   24.2578    3.0000
   28.2195    3.0000
   28.2192    0.0000
   32.1595    0.0000
   32.1592    3.0000
   36.0838    3.0000
   36.0837    0.0000
   40.0000    3.0000
   40.0000    0.0000];
end

function y = N
y = zeros(4,1,4);
y(:,:,1) = [0.4167; 0.4167; 0.0833; 0.0833];
y(:,:,2) = [0.0833; 0.4167; 0.4167; 0.0833];
y(:,:,3) = [0.0833; 0.0833; 0.4167; 0.4167];
y(:,:,4) = [0.4167; 0.0833; 0.0833; 0.4167];
end

function y = dNdxi
y = zeros(4,2,4);
y(:,:,1) = [-0.1667   -0.6667
    0.6667    0.1667
   -0.1667    0.3333
   -0.3333    0.1667];
y(:,:,2) = [   -0.1667   -0.3333
    0.6667   -0.1667
   -0.1667    0.6667
   -0.3333   -0.1667];
y(:,:,3) = [
    0.1667   -0.3333
    0.3333   -0.1667
    0.1667    0.6667
   -0.6667   -0.1667];
y(:,:,4) = [    0.1667   -0.6667
    0.3333    0.1667
    0.1667    0.3333
   -0.6667    0.1667];
end