clear;
clc;

fprintf('\n####################################################################\n')
fprintf(' Equations of motion in classical physics ')
fprintf('\n####################################################################\n')
fprintf(' Kayode Olumoyin ')
fprintf(' COMS 7100 ')
fprintf(' Project 3 ')
fprintf(' Part 2 ')
fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long e % using the highest matlab precision

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 2
% Planetary Orbits of Mercury,Venus,Earth,Mars,Jujpiter,Saturn,Uranus,Neptune
% Using only the gravitational force of the sun

G = 6.67428e-11;            % Gravitational Constant
N = 1e6;

set(0,'DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sun data

Msun = 1.98855e30;
pdsun = [0 0 0];
pvsun = [0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mercury data

Mmercury = 0.3301e24;       % kg    % Mass of sun
pdmercury = [-1.946172639275932e10 -6.691327522588462e10 -3.679854414596553e9];
pvmercury = [3.699499188030234e4 -1.116441595690670e4 -4.307627980092748e3];

% Pdmercury = 46e9;           % m    % Perihelion distance to Sun
% Pvmercury = 58.980e3;       % m/s  % Perihelion spped to Sun
Tmercury = 87.97*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% venus data

Mvenus = 4.8676e24;         % kg
pdvenus = [-1.074564940489116e11 -4.885015029930510e9 6.135634314000621e9];
pvvenus = [1.381906047920155e3 -3.514029517606325e4 -5.600423209496953e2];

% Pdvenus = 107.48e9;         % m
% Pvvenus = 35.26e3;          % m/s
Tvenus = 224.701*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% earth data

Mearth =  5.9726e24;        % kg
pdearth = [-2.649903375682292e10 1.446972967792532e11 -6.112214304122390e5];
pvearth = [-2.979426006719171e4 -5.469294956024861e3 0.1817900213849671e0];

% Pdearth = 147.09e9;         % m
% Pvearth = 30.29e3;          % m/s
Tearth = 365.256*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mars data

Mmars = 0.64174e24;         % kg
pdmars = [2.080481406481886e11 -2.007052475139214e9 -5.156288876678156e9];
pvmars = [1.162672383641686e3 2.629606454658444e4 5.222970037729194e2];

% Pdmars = 206.62e9;          % m
% Pvmars = 26.50e3;           % m/s
Tmars = 686.980*24*3600;    % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jupiter data

Mjupiter = 1898.3e24;       % kg
pdjupiter = [5.985676196098850e11 4.396046864373181e11 -1.522689720543665e10];
pvjupiter = [-7.909860396044607e3 1.115621744689623e4 1.308653756345028e2];

% Pdjupiter = 740.52e9;       % m
% Pvjupiter = 13.72e3;        % m/s
Tjupiter = 4332.59*24*3600; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saturn data

Msaturn = 568.36e24;        % kg
pdsaturn = [9.583853585628320e11 9.828562834129111e11 -5.521298021048592e10];
pvsaturn = [-7.431212962376309e3 6.736755751121786e3 1.777382936720429e2];

% Pdsaturn = 1352.55e9;       % m
% Pvsaturn = 10.18e3;         % m/s
Tsaturn = 10759.22*24*3600; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uranus data

Muranus = 86.816e24;        % kg
pduranus = [2.158975041926949e12 -2.054625328046016e12 -3.562533932149468e10];
pvuranus = [4.637273675667855e3 4.627598915337914e3 -4.291762022516144e1];

% Pduranus = 2741.30e9;       % m
% Pvuranus = 7.11e3;          % m/s
Turanus = 30685.4*24*3600;  % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neptune data

Mneptune = 102.42e24;       % kg
pdneptune = [2.515046330944007e12 -3.738714537156832e12 1.903228410412154e10];
pvneptune = [4.465277067705161e3 3.075982230271975e3 -1.662486247561596e2];

%Pdneptune = 4444.45e9;      % m
%Pvneptune = 5.50e3;         % m/s
Tneptune = 60189*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mass = [Msun Mmercury Mvenus Mearth Mmars Mjupiter Msaturn Muranus Mneptune];
pos = [pdsun pdmercury pdvenus pdearth pdmars pdjupiter pdsaturn pduranus pdneptune];
sp = [pvsun pvmercury pvvenus pvearth pvmars pvjupiter pvsaturn pvuranus pvneptune];
position = reshape(pos,[3,length(pos)/3]);
speed = reshape(sp,[3,length(sp)/3]);
% %%%%%%%%%%%%%%%%%%%%% euler method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n euler method 3d \n')

fprintf('\n Mercury \n')
Y_e_m = euler_3d(G,Mass(2),position(:,2),Mass,position,speed,5200380000,57782,90000);
%plot3(Y_e_m(:,1), Y_e_m(:,2), Y_e_m(:,3),'g')

% %%%%%%%%%%%%%%%%%%%%% cromer method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n cromer method 3d \n')

fprintf('\n Mercury \n')
Y_c_m = cromer_3d(G,Mass(2),position(:,2),Mass,position,speed(:,2),5200380000,57782,90000);
Y_c_v = cromer_3d(G,Mass(3),position(:,3),Mass,position,speed(:,3),5200380000,57782,90000);
Y_c_e = cromer_3d(G,Mass(4),position(:,4),Mass,position,speed(:,4),5200380000,57782,90000);
Y_c_ma = cromer_3d(G,Mass(5),position(:,5),Mass,position,speed(:,5),5200380000,57782,90000);
Y_c_j = cromer_3d(G,Mass(6),position(:,6),Mass,position,speed(:,6),5200380000,57782,90000);
Y_c_s = cromer_3d(G,Mass(7),position(:,7),Mass,position,speed(:,7),5200380000,57782,90000);
Y_c_u = cromer_3d(G,Mass(8),position(:,8),Mass,position,speed(:,8),5200380000,57782,90000);
Y_c_n = cromer_3d(G,Mass(9),position(:,9),Mass,position,speed(:,9),5200380000,57782,90000);

figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot3(Y_c_m(:,1), Y_c_m(:,2), Y_c_m(:,3),'g')
hold on
plot3(Y_c_v(:,1), Y_c_v(:,2), Y_c_v(:,3),'m')
plot3(Y_c_e(:,1), Y_c_e(:,2), Y_c_e(:,3),'b')
plot3(Y_c_ma(:,1), Y_c_ma(:,2), Y_c_ma(:,3),'r')
plot3(Y_c_j(:,1), Y_c_j(:,2), Y_c_j(:,3),'k')
plot3(Y_c_s(:,1), Y_c_s(:,2), Y_c_s(:,3),'y')
plot3(Y_c_u(:,1), Y_c_u(:,2), Y_c_u(:,3),'Color',color)
plot3(Y_c_n(:,1), Y_c_n(:,2), Y_c_n(:,3),'c')
hold off

toc
%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = cromer_3d(G, mass_p, pos_p, Mass, position, speed_p, T, N, dt)
rr0 = pos_p;
vv0 = speed_p; 
r = zeros(3,N);
v = zeros(3,N);
r(:,1) = rr0;
v(:,1) = vv0;

for i = 1:N-1
    f = f_attr(G, mass_p, pos_p, Mass, position);
    v(:,i+1) = v(:,i) + f*r(:,i)*dt;
    r(:,i+1) = r(:,i) + v(:,i+1)*dt;
    
end
    
Y = [r(1,:)' r(2,:)' r(3,:)' ];
end

function Y = euler_3d(G, mass_p, pos_p, Mass, position, speed, T, N, dt)
rr0 = position(:,2);
vv0 = speed(:,2); 
r = zeros(3,N);
v = zeros(3,N);
r(:,1) = rr0;
v(:,1) = vv0;

for i = 1:N-1
    r(:,i+1) = r(:,i) + v(:,i)*dt;
    f = f_attr(G, mass_p, pos_p, Mass, position);
    v(:,i+1) = v(:,i) + f*r(:,i)*dt;
end
    
Y = [r(1,:)' r(2,:)' r(3,:)' ];
end

% Total force of attraction
function f = f_attr(G, mass_p, pos_p, Mass, position)
r1 = pos_p - position(:,1);
if norm(r1) == 0
    f1 = 0;
else
    f1 = (1/mass_p)*(-G*mass_p*Mass(1)/norm(r1)^3);
end
r2 = pos_p - position(:,2);
if norm(r2) == 0
    f2 = 0;
else
    f2 = (1/mass_p)*(-G*mass_p*Mass(2)/norm(r2)^3);
end
r3 = pos_p - position(:,3);
if norm(r3) == 0
    f3 = 0;
else
    f3 = (1/mass_p)*(-G*mass_p*Mass(3)/norm(r3)^3);
end
r4 = pos_p - position(:,4);
if norm(r4) == 0
    f4 = 0;
else
    f4 = (1/mass_p)*(-G*mass_p*Mass(4)/norm(r4)^3);
end
r5 = pos_p - position(:,5);
if norm(r5) == 0
    f5 = 0;
else
    f5 = (1/mass_p)*(-G*mass_p*Mass(5)/norm(r5)^3);
end
r6 = pos_p - position(:,6);
if norm(r6) == 0
    f6 = 0;
else
    f6 = (1/mass_p)*(-G*mass_p*Mass(6)/norm(r6)^3);
end
r7 = pos_p - position(:,7);
if norm(r7) == 0
    f7 = 0;
else
    f7 = (1/mass_p)*(-G*mass_p*Mass(7)/norm(r7)^3);
end
r8 = pos_p - position(:,8);
if norm(r8) == 0
    f8 = 0;
else
    f8 = (1/mass_p)*(-G*mass_p*Mass(8)/norm(r8)^3);
end
r9 = pos_p - position(:,9);
if norm(r9) == 0
    f9 = 0;
else
    f9 = (1/mass_p)*(-G*mass_p*Mass(9)/norm(r9)^3);
end
f = f1+f2+f3+f4+f5+f6+f7+f8+f9;
end



















