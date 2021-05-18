clear;
clc;

fprintf('\n####################################################################\n')
fprintf(' Equations of motion in classical physics ')
fprintf('\n####################################################################\n')
fprintf(' Kayode Olumoyin ')
fprintf(' COMS 7100 ')
fprintf(' Project 3 ')
fprintf(' Part 1 ')
fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long e % using the highest matlab precision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART 1
% Planetary Orbits of Mercury,Venus,Earth,Mars,Jujpiter,Saturn,Uranus,Neptune
% Using only the gravitational force of the sun

G = 6.67428e-11;            % Gravitational Constant
N = 1e6;

set(0,'DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sun data

Msun = 1.98855e30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mercury data

Mmercury = 0.3301e24;       % kg    % Mass of sun
Pdmercury = 46e9;           % m    % Perihelion distance to Sun
Pvmercury = 58.980e3;       % m/s  % Perihelion spped to Sun
Tmercury = 87.97*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% venus data

Mvenus = 4.8676e24;         % kg
Pdvenus = 107.48e9;         % m
Pvvenus = 35.26e3;          % m/s
Tvenus = 224.701*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% earth data

Mearth =  5.9726e24;        % kg
Pdearth = 147.09e9;         % m
Pvearth = 30.29e3;          % m/s
Tearth = 365.256*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mars data

Mmars = 0.64174e24;         % kg
Pdmars = 206.62e9;          % m
Pvmars = 26.50e3;           % m/s
Tmars = 686.980*24*3600;    % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jupiter data

Mjupiter = 1898.3e24;       % kg
Pdjupiter = 740.52e9;       % m
Pvjupiter = 13.72e3;        % m/s
Tjupiter = 4332.59*24*3600; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saturn data

Msaturn = 568.36e24;        % kg
Pdsaturn = 1352.55e9;       % m
Pvsaturn = 10.18e3;         % m/s
Tsaturn = 10759.22*24*3600; % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uranus data

Muranus = 86.816e24;        % kg
Pduranus = 2741.30e9;       % m
Pvuranus = 7.11e3;          % m/s
Turanus = 30685.4*24*3600;  % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neptune data

Mneptune = 102.42e24;       % kg
Pdneptune = 4444.45e9;      % m
Pvneptune = 5.50e3;         % m/s
Tneptune = 60189*24*3600;   % s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%% euler method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n euler method \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mercury
fprintf('\n Mercury \n')
Y_e_m = euler_2d(G,Msun,Pdmercury,Pvmercury,Tmercury,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tmercury/N)/86400, Y_e_m(1,5)/1000, Y_e_m(1,6)/1000, Y_e_m(1,7)]);
T_e_m = length(Y_e_m(:,1));
for kk = 1:T_e_m/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tmercury/N)*kk*5000/86400, Y_e_m(kk*5000,5)/1000, Y_e_m(kk*5000,6)/1000, Y_e_m(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_m
    if abs(Y_e_m(jj,7)-180) < 0.0001 && Y_e_m(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_m(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_m(jj,5)/1000,[Y_e_m(jj,3)/1000 Y_e_m(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_m(jj,7));
        fprintf(' Time (days) = %f\n',(Tmercury/N)*jj/86400);
        aphe_m = Y_e_m(jj,6);
    end
end
fprintf(' \n other selected parameters from simulation \n ')

side_m = (Tmercury/N)*T_e_m/86400;
a_m = (aphe_m+Pdmercury)/2000;
c_m = a_m - (Pdmercury/1000);
e_m = c_m/a_m;
mean_m = mean(Y_e_m(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_m);
fprintf(' Sidereal orbit period (days) = %f \n',side_m);
fprintf(' orbit eccentricity = %f \n',e_m);
fprintf(' mean orbital spped = %f km/s \n',mean_m);

subplot(2,4,1)
plot(Y_e_m(:,1), Y_e_m(:,2),'g')
%title(' Mercury Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Mercury Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Venus
fprintf('\n Venus \n')
Y_e_v = euler_2d(G,Msun,Pdvenus,Pvvenus,Tvenus,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tvenus/N)/86400, Y_e_v(1,5)/1000, Y_e_v(1,6)/1000, Y_e_v(1,7)]);
T_e_v = length(Y_e_v(:,1));
for kk = 1:T_e_v/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tvenus/N)*kk*5000/86400, Y_e_v(kk*5000,5)/1000, Y_e_v(kk*5000,6)/1000, Y_e_v(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_v
    if abs(Y_e_v(jj,7)-180) < 0.0005 && Y_e_v(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_v(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_v(jj,5)/1000,[Y_e_v(jj,3)/1000 Y_e_v(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_v(jj,7));
        fprintf(' Time (days) = %f\n',(Tvenus/N)*jj/86400);
        aphe_v = Y_e_v(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_v = (Tvenus/N)*T_e_v/86400;
a_v = (aphe_v+Pdvenus)/2000;
c_v = a_v - (Pdvenus/1000);
e_v = c_v/a_v;
mean_v = mean(Y_e_v(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_v);
fprintf(' Sidereal orbit period (days) = %f \n',side_v);
fprintf(' orbit eccentricity = %f \n',e_v);
fprintf(' mean orbital spped = %f km/s \n',mean_v);

subplot(2,4,2)
plot(Y_e_v(:,1), Y_e_v(:,2),'m')
%title(' Venus Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Venus Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth
fprintf('\n Earth \n')
Y_e_e = euler_2d(G,Msun,Pdearth,Pvearth,Tearth,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tearth/N)/86400, Y_e_e(1,5)/1000, Y_e_e(1,6)/1000, Y_e_e(1,7)]);
T_e_e = length(Y_e_e(:,1));
for kk = 1:T_e_e/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tearth/N)*kk*5000/86400, Y_e_e(kk*5000,5)/1000, Y_e_e(kk*5000,6)/1000, Y_e_e(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_e
    if abs(Y_e_e(jj,7)-180) < 0.0005 && Y_e_e(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_e(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_e(jj,5)/1000,[Y_e_e(jj,3)/1000 Y_e_e(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_e(jj,7));
        fprintf(' Time (days) = %f\n',(Tearth/N)*jj/86400);
        aphe_e = Y_e_e(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_e = (Tearth/N)*T_e_e/86400;
a_e = (aphe_e+Pdearth)/2000;
c_e = a_e - (Pdearth/1000);
e_e = c_e/a_e;
mean_e = mean(Y_e_e(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_e);
fprintf(' Sidereal orbit period (days) = %f \n',side_e);
fprintf(' orbit eccentricity = %f \n',e_e);
fprintf(' mean orbital spped = %f km/s \n',mean_e);

subplot(2,4,3)
plot(Y_e_e(:,1), Y_e_e(:,2),'b')
%title(' Earth Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Earth Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mars
fprintf('\n Mars \n')
Y_e_ma = euler_2d(G,Msun,Pdmars,Pvmars,Tmars,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tmars/N)/86400, Y_e_ma(1,5)/1000, Y_e_ma(1,6)/1000, Y_e_ma(1,7)]);
T_e_ma = length(Y_e_ma(:,1));
for kk = 1:T_e_ma/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tmars/N)*kk*5000/86400, Y_e_ma(kk*5000,5)/1000, Y_e_ma(kk*5000,6)/1000, Y_e_ma(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_ma
    if abs(Y_e_ma(jj,7)-180) < 0.0002 && Y_e_ma(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_ma(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_ma(jj,5)/1000,[Y_e_ma(jj,3)/1000 Y_e_ma(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_ma(jj,7));
        fprintf(' Time (days) = %f\n',(Tmars/N)*jj/86400);
        aphe_ma = Y_e_ma(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_ma = (Tmars/N)*T_e_ma/86400;
a_ma = (aphe_ma+Pdmars)/2000;
c_ma = a_ma - (Pdmars/1000);
e_ma = c_ma/a_ma;
mean_ma = mean(Y_e_ma(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_ma);
fprintf(' Sidereal orbit period (days) = %f \n',side_ma);
fprintf(' orbit eccentricity = %f \n',e_ma);
fprintf(' mean orbital spped = %f km/s \n',mean_ma);

subplot(2,4,4)
plot(Y_e_ma(:,1), Y_e_ma(:,2),'r')
%title(' Mars Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Mars Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jupiter
fprintf('\n Jupiter \n')
Y_e_j = euler_2d(G,Msun,Pdjupiter,Pvjupiter,Tjupiter,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',1, ...
        [(Tjupiter/N)/86400, Y_e_j(1,5)/1000, Y_e_j(1,6)/1000, Y_e_j(1,7)]);
T_e_j = length(Y_e_j(:,1));
for kk = 1:T_e_j/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tjupiter/N)*kk*5000/86400, Y_e_j(kk*5000,5)/1000, Y_e_j(kk*5000,6)/1000, Y_e_j(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_j
    if abs(Y_e_j(jj,7)-180) < 0.0002 && Y_e_j(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_j(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_j(jj,5)/1000,[Y_e_j(jj,3)/1000 Y_e_j(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_j(jj,7));
        fprintf(' Time (days) = %f\n',(Tjupiter/N)*jj/86400);
        aphe_j = Y_e_j(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_j = (Tjupiter/N)*T_e_j/86400;
a_j = (aphe_j+Pdjupiter)/2000;
c_j = a_j - (Pdjupiter/1000);
e_j = c_j/a_j;
mean_j = mean(Y_e_j(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_j);
fprintf(' Sidereal orbit period (days) = %f \n',side_j);
fprintf(' orbit eccentricity = %f \n',e_j);
fprintf(' mean orbital spped = %f km/s \n',mean_j);

subplot(2,4,5)
plot(Y_e_j(:,1), Y_e_j(:,2),'k')
%title(' Jupiter Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Jupiter Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saturn
fprintf('\n Saturn \n')
Y_e_s = euler_2d(G,Msun,Pdsaturn,Pvsaturn,Tsaturn,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tsaturn/N)/86400, Y_e_s(1,5)/1000, Y_e_s(1,6)/1000, Y_e_s(1,7)]);
T_e_s = length(Y_e_s(:,1));
for kk = 1:T_e_s/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tsaturn/N)*kk*5000/86400, Y_e_s(kk*5000,5)/1000, Y_e_s(kk*5000,6)/1000, Y_e_s(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_s
    if abs(Y_e_s(jj,7)-180) < 0.0002 && Y_e_s(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_s(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_s(jj,5)/1000,[Y_e_s(jj,3)/1000 Y_e_s(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_s(jj,7));
        fprintf(' Time (days) = %f\n',(Tsaturn/N)*jj/86400);
        aphe_s = Y_e_s(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_s = (Tsaturn/N)*T_e_s/86400;
a_s = (aphe_s+Pdsaturn)/2000;
c_s = a_s - (Pdsaturn/1000);
e_s = c_s/a_s;
mean_s = mean(Y_e_s(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_s);
fprintf(' Sidereal orbit period (days) = %f \n',side_s);
fprintf(' orbit eccentricity = %f \n',e_s);
fprintf(' mean orbital spped = %f km/s \n',mean_s);

subplot(2,4,6)
plot(Y_e_s(:,1), Y_e_s(:,2),'y')
%title(' Saturn Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Saturn Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uranus
fprintf('\n Uranus \n')
Y_e_u = euler_2d(G,Msun,Pduranus,Pvuranus,Turanus,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Turanus/N)/86400, Y_e_u(1,5)/1000, Y_e_u(1,6)/1000, Y_e_u(1,7)]);
T_e_u = length(Y_e_u(:,1));
for kk = 1:T_e_u/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Turanus/N)*kk*5000/86400, Y_e_u(kk*5000,5)/1000, Y_e_u(kk*5000,6)/1000, Y_e_u(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_u
    if abs(Y_e_u(jj,7)-180) < 0.0005 && Y_e_u(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_u(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_u(jj,5)/1000,[Y_e_u(jj,3)/1000 Y_e_u(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_u(jj,7));
        fprintf(' Time (days) = %f\n',(Turanus/N)*jj/86400);
        aphe_u = Y_e_u(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_u = (Turanus/N)*T_e_u/86400;
a_u = (aphe_u+Pduranus)/2000;
c_u = a_u - (Pduranus/1000);
e_u = c_u/a_u;
mean_u = mean(Y_e_u(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_u);
fprintf(' Sidereal orbit period (days) = %f \n',side_u);
fprintf(' orbit eccentricity = %f \n',e_u);
fprintf(' mean orbital spped = %f km/s \n',mean_u);

subplot(2,4,7)
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_e_u(:,1), Y_e_u(:,2),'Color',color)
%title(' Uranus Euler method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Uranus Euler method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neptune
fprintf('\n Neptune \n')
Y_e_n = euler_2d(G,Msun,Pdneptune,Pvneptune,Tneptune,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f   \n',1, ...
        [(Tneptune/N)/86400, Y_e_n(1,5)/1000, Y_e_n(1,6)/1000, Y_e_n(1,7)]);
T_e_n = length(Y_e_n(:,1));
for kk = 1:T_e_n/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tneptune/N)*kk*5000/86400, Y_e_n(kk*5000,5)/1000, Y_e_n(kk*5000,6)/1000, Y_e_n(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_e_n
    if abs(Y_e_n(jj,7)-180) < 0.0005 && Y_e_n(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_e_n(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_e_n(jj,5)/1000,[Y_e_n(jj,3)/1000 Y_e_n(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_e_n(jj,7));
        fprintf(' Time (days) = %f\n',(Tneptune/N)*jj/86400);
        aphe_n = Y_e_n(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_n = (Tneptune/N)*T_e_n/86400;
a_n = (aphe_n+Pdneptune)/2000;
c_n = a_n - (Pdneptune/1000);
e_n = c_n/a_n;
mean_n = mean(Y_e_n(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_n);
fprintf(' Sidereal orbit period (days) = %f \n',side_n);
fprintf(' orbit eccentricity = %f \n',e_n);
fprintf(' mean orbital spped = %f km/s \n',mean_n);

subplot(2,4,8)
plot(Y_e_n(:,1), Y_e_n(:,2),'c')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Neptune Euler method'])
ylabel('position in y-axis')

sgtitle('Planet orbits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all terrestrial planets Eulers method
figure
plot(Y_e_m(:,1), Y_e_m(:,2),'g')
hold on
plot(Y_e_v(:,1), Y_e_v(:,2),'m')
plot(Y_e_e(:,1), Y_e_e(:,2),'b')
plot(Y_e_ma(:,1), Y_e_ma(:,2),'r')
hold off
title(' all terrestrial planets orbit using Euler method ')
legend({'Mercury','Venus','Earth','Mars'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e11 5e11])
ylim([-5e11 5e11])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all Jovian planets Eulers method
figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_e_j(:,1), Y_e_j(:,2),'k')
hold on
plot(Y_e_s(:,1), Y_e_s(:,2),'y')
plot(Y_e_u(:,1), Y_e_u(:,2),'Color',color)
plot(Y_e_n(:,1), Y_e_n(:,2),'c')
hold off
title(' all jovian planets orbit using Euler method ')
legend({'Jupiter','Saturn','Uranus','Neptune'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all planets Eulers method
figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_e_m(:,1), Y_e_m(:,2),'g')
hold on
plot(Y_e_v(:,1), Y_e_v(:,2),'m')
plot(Y_e_e(:,1), Y_e_e(:,2),'b')
plot(Y_e_ma(:,1), Y_e_ma(:,2),'r')
plot(Y_e_j(:,1), Y_e_j(:,2),'k')
plot(Y_e_s(:,1), Y_e_s(:,2),'y')
plot(Y_e_u(:,1), Y_e_u(:,2),'Color',color)
plot(Y_e_n(:,1), Y_e_n(:,2),'c')
hold off
title(' all planets orbit using Euler method ')
legend({'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%% cromer method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n cromer method \n')

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mercury
fprintf('\n Mercury \n')
Y_c_m = cromer_2d(G,Msun,Pdmercury,Pvmercury,Tmercury,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tmercury/N)/86400, Y_c_m(1,5)/1000, Y_c_m(1,6)/1000, Y_c_m(1,7)]);
T_c_m = length(Y_c_m(:,1));
for kk = 1:T_c_m/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tmercury/N)*kk*5000/86400, Y_c_m(kk*5000,5)/1000, Y_c_m(kk*5000,6)/1000, Y_c_m(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_m
    if abs(Y_c_m(jj,7)-180) < 0.0003 && Y_c_m(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_m(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_m(jj,5)/1000,[Y_c_m(jj,3)/1000 Y_c_m(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_m(jj,7));
        fprintf(' Time (days) = %f\n',(Tmercury/N)*jj/86400);
        aphe_c_m = Y_c_m(jj,6);
    end
end
fprintf(' \n other selected parameters from simulation \n ')

side_c_m = (Tmercury/N)*T_c_m/86400;
a_c_m = (aphe_c_m + Pdmercury)/2000;
c_c_m = a_c_m - (Pdmercury/1000);
e_c_m = c_c_m/a_c_m;
mean_c_m = mean(Y_c_m(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_m);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_m);
fprintf(' orbit eccentricity = %f \n',e_c_m);
fprintf(' mean orbital spped = %f km/s \n',mean_c_m);

figure
subplot(2,4,1)
plot(Y_c_m(:,1), Y_c_m(:,2),'g')
%title(' Mercury Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Mercury Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Venus
fprintf('\n Venus \n')
Y_c_v = cromer_2d(G,Msun,Pdvenus,Pvvenus,Tvenus,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tvenus/N)/86400, Y_c_v(1,5)/1000, Y_c_v(1,6)/1000, Y_c_v(1,7)]);
T_c_v = length(Y_c_v(:,1));
for kk = 1:T_c_v/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tvenus/N)*kk*5000/86400, Y_c_v(kk*5000,5)/1000, Y_c_v(kk*5000,6)/1000, Y_c_v(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_v
    if abs(Y_c_v(jj,7)-180) < 0.0003 && Y_c_v(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_v(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_v(jj,5)/1000,[Y_c_v(jj,3)/1000 Y_c_v(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_v(jj,7));
        fprintf(' Time (days) = %f\n',(Tvenus/N)*jj/86400);
        aphe_c_v = Y_c_v(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_v = (Tvenus/N)*T_c_v/86400;
a_c_v = (aphe_c_v+Pdvenus)/2000;
c_c_v = a_c_v - (Pdvenus/1000);
e_c_v = c_c_v/a_c_v;
mean_c_v = mean(Y_c_v(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_v);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_v);
fprintf(' orbit eccentricity = %f \n',e_c_v);
fprintf(' mean orbital spped = %f km/s \n',mean_c_v);

subplot(2,4,2)
plot(Y_c_v(:,1), Y_c_v(:,2),'m')
%title(' Venus Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Venus Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth
fprintf('\n Earth \n')
Y_c_e = cromer_2d(G,Msun,Pdearth,Pvearth,Tearth,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tearth/N)/86400, Y_c_e(1,5)/1000, Y_c_e(1,6)/1000, Y_c_e(1,7)]);
T_c_e = length(Y_c_e(:,1));
for kk = 1:T_c_e/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tearth/N)*kk*5000/86400, Y_c_e(kk*5000,5)/1000, Y_c_e(kk*5000,6)/1000, Y_c_e(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_e
    if abs(Y_c_e(jj,7)-180) < 0.0005 && Y_c_e(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_e(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_e(jj,5)/1000,[Y_c_e(jj,3)/1000 Y_c_e(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_e(jj,7));
        fprintf(' Time (days) = %f\n',(Tearth/N)*jj/86400);
        aphe_c_e = Y_c_e(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_e = (Tearth/N)*T_c_e/86400;
a_c_e = (aphe_c_e+Pdearth)/2000;
c_c_e = a_c_e - (Pdearth/1000);
e_c_e = c_c_e/a_c_e;
mean_c_e = mean(Y_c_e(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_e);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_e);
fprintf(' orbit eccentricity = %f \n',e_c_e);
fprintf(' mean orbital spped = %f km/s \n',mean_c_e);

subplot(2,4,3)
plot(Y_c_e(:,1), Y_c_e(:,2),'b')
%title(' Earth Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Earth Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mars
fprintf('\n Mars \n')
Y_c_ma = cromer_2d(G,Msun,Pdmars,Pvmars,Tmars,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tmars/N)/86400, Y_c_ma(1,5)/1000, Y_c_ma(1,6)/1000, Y_c_ma(1,7)]);
T_c_ma = length(Y_c_ma(:,1));
for kk = 1:T_c_ma/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tmars/N)*kk*5000/86400, Y_c_ma(kk*5000,5)/1000, Y_c_ma(kk*5000,6)/1000, Y_c_ma(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_ma
    if abs(Y_c_ma(jj,7)-180) < 0.0004 && Y_c_ma(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_ma(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_ma(jj,5)/1000,[Y_c_ma(jj,3)/1000 Y_c_ma(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_ma(jj,7));
        fprintf(' Time (days) = %f\n',(Tmars/N)*jj/86400);
        aphe_c_ma = Y_c_ma(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_ma = (Tmars/N)*T_c_ma/86400;
a_c_ma = (aphe_c_ma+Pdmars)/2000;
c_c_ma = a_c_ma - (Pdmars/1000);
e_c_ma = c_c_ma/a_c_ma;
mean_c_ma = mean(Y_c_ma(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_ma);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_ma);
fprintf(' orbit eccentricity = %f \n',e_c_ma);
fprintf(' mean orbital spped = %f km/s \n',mean_c_ma);

subplot(2,4,4)
plot(Y_c_ma(:,1), Y_c_ma(:,2),'r')
%title(' Mars Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Mars Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jupiter
fprintf('\n Jupiter \n')
Y_c_j = cromer_2d(G,Msun,Pdjupiter,Pvjupiter,Tjupiter,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',1, ...
        [(Tjupiter/N)/86400, Y_c_j(1,5)/1000, Y_c_j(1,6)/1000, Y_c_j(1,7)]);
T_c_j = length(Y_c_j(:,1));
for kk = 1:T_c_j/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tjupiter/N)*kk*5000/86400, Y_c_j(kk*5000,5)/1000, Y_c_j(kk*5000,6)/1000, Y_c_j(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_j
    if abs(Y_c_j(jj,7)-180) < 0.00035 && Y_c_j(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_j(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_j(jj,5)/1000,[Y_c_j(jj,3)/1000 Y_c_j(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_j(jj,7));
        fprintf(' Time (days) = %f\n',(Tjupiter/N)*jj/86400);
        aphe_c_j = Y_c_j(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_j = (Tjupiter/N)*T_c_j/86400;
a_c_j = (aphe_c_j+Pdjupiter)/2000;
c_c_j = a_c_j - (Pdjupiter/1000);
e_c_j = c_c_j/a_c_j;
mean_c_j = mean(Y_c_j(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_j);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_j);
fprintf(' orbit eccentricity = %f \n',e_c_j);
fprintf(' mean orbital spped = %f km/s \n',mean_c_j);

subplot(2,4,5)
plot(Y_c_j(:,1), Y_c_j(:,2),'k')
%title(' Jupiter Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Jupiter Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saturn
fprintf('\n Saturn \n')
Y_c_s = cromer_2d(G,Msun,Pdsaturn,Pvsaturn,Tsaturn,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tsaturn/N)/86400, Y_c_s(1,5)/1000, Y_c_s(1,6)/1000, Y_c_s(1,7)]);
T_c_s = length(Y_c_s(:,1));
for kk = 1:T_c_s/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tsaturn/N)*kk*5000/86400, Y_c_s(kk*5000,5)/1000, Y_c_s(kk*5000,6)/1000, Y_c_s(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_s
    if abs(Y_c_s(jj,7)-180) < 0.0005 && Y_c_s(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_s(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_s(jj,5)/1000,[Y_c_s(jj,3)/1000 Y_c_s(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_s(jj,7));
        fprintf(' Time (days) = %f\n',(Tsaturn/N)*jj/86400);
        aphe_c_s = Y_c_s(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_s = (Tsaturn/N)*T_c_s/86400;
a_c_s = (aphe_c_s+Pdsaturn)/2000;
c_c_s = a_c_s - (Pdsaturn/1000);
e_c_s = c_c_s/a_c_s;
mean_c_s = mean(Y_c_s(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_s);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_s);
fprintf(' orbit eccentricity = %f \n',e_c_s);
fprintf(' mean orbital spped = %f km/s \n',mean_c_s);

subplot(2,4,6)
plot(Y_c_s(:,1), Y_c_s(:,2),'y')
%title(' Saturn Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Saturn Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uranus
fprintf('\n Uranus \n')
Y_c_u = cromer_2d(G,Msun,Pduranus,Pvuranus,Turanus,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Turanus/N)/86400, Y_c_u(1,5)/1000, Y_c_u(1,6)/1000, Y_c_u(1,7)]);
T_c_u = length(Y_c_u(:,1));
for kk = 1:T_c_u/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Turanus/N)*kk*5000/86400, Y_c_u(kk*5000,5)/1000, Y_c_u(kk*5000,6)/1000, Y_c_u(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_u
    if abs(Y_c_u(jj,7)-180) < 0.0004 && Y_c_u(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_u(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_u(jj,5)/1000,[Y_c_u(jj,3)/1000 Y_c_u(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_u(jj,7));
        fprintf(' Time (days) = %f\n',(Turanus/N)*jj/86400);
        aphe_c_u = Y_c_u(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_u = (Turanus/N)*T_c_u/86400;
a_c_u = (aphe_c_u+Pduranus)/2000;
c_c_u = a_c_u - (Pduranus/1000);
e_c_u = c_c_u/a_c_u;
mean_c_u = mean(Y_c_u(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_u);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_u);
fprintf(' orbit eccentricity = %f \n',e_c_u);
fprintf(' mean orbital spped = %f km/s \n',mean_c_u);

subplot(2,4,7)
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_c_u(:,1), Y_c_u(:,2),'Color',color)
%title(' Uranus Cromer method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Uranus Cromer method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neptune
fprintf('\n Neptune \n')
Y_c_n = cromer_2d(G,Msun,Pdneptune,Pvneptune,Tneptune,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f   \n',1, ...
        [(Tneptune/N)/86400, Y_c_n(1,5)/1000, Y_c_n(1,6)/1000, Y_c_n(1,7)]);
T_c_n = length(Y_c_n(:,1));
for kk = 1:T_c_n/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tneptune/N)*kk*5000/86400, Y_c_n(kk*5000,5)/1000, Y_c_n(kk*5000,6)/1000, Y_c_n(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_c_n
    if abs(Y_c_n(jj,7)-180) < 0.0005 && Y_c_n(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_c_n(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_c_n(jj,5)/1000,[Y_c_n(jj,3)/1000 Y_c_n(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_c_n(jj,7));
        fprintf(' Time (days) = %f\n',(Tneptune/N)*jj/86400);
        aphe_c_n = Y_c_n(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_c_n = (Tneptune/N)*T_c_n/86400;
a_c_n = (aphe_c_n+Pdneptune)/2000;
c_c_n = a_c_n - (Pdneptune/1000);
e_c_n = c_c_n/a_c_n;
mean_c_n = mean(Y_c_n(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_c_n);
fprintf(' Sidereal orbit period (days) = %f \n',side_c_n);
fprintf(' orbit eccentricity = %f \n',e_c_n);
fprintf(' mean orbital spped = %f km/s \n',mean_c_n);

subplot(2,4,8)
plot(Y_c_n(:,1), Y_c_n(:,2),'c')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Neptune Cromer method'])
ylabel('position in y-axis')

sgtitle('Planet orbits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all terrestrial planets Cromers method
figure
plot(Y_c_m(:,1), Y_c_m(:,2),'g')
hold on
plot(Y_c_v(:,1), Y_c_v(:,2),'m')
plot(Y_c_e(:,1), Y_c_e(:,2),'b')
plot(Y_c_ma(:,1), Y_c_ma(:,2),'r')
hold off
title(' all terrestrial planets orbit using Cromer method ')
legend({'Mercury','Venus','Earth','Mars'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e11 5e11])
ylim([-5e11 5e11])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all Jovian planets Cromers method
figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_c_j(:,1), Y_c_j(:,2),'k')
hold on
plot(Y_c_s(:,1), Y_c_s(:,2),'y')
plot(Y_c_u(:,1), Y_c_u(:,2),'Color',color)
plot(Y_c_n(:,1), Y_c_n(:,2),'c')
hold off
title(' all jovian planets orbit using Cromer method ')
legend({'Jupiter','Saturn','Uranus','Neptune'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all planets Cromers method
figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_c_m(:,1), Y_c_m(:,2),'g')
hold on
plot(Y_c_v(:,1), Y_c_v(:,2),'m')
plot(Y_c_e(:,1), Y_c_e(:,2),'b')
plot(Y_c_ma(:,1), Y_c_ma(:,2),'r')
plot(Y_c_j(:,1), Y_c_j(:,2),'k')
plot(Y_c_s(:,1), Y_c_s(:,2),'y')
plot(Y_c_u(:,1), Y_c_u(:,2),'Color',color)
plot(Y_c_n(:,1), Y_c_n(:,2),'c')
hold off
title(' all planets orbit using Cromer method ')
legend({'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% rk2 method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n rk2 method \n')

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mercury
fprintf('\n Mercury \n')
Y_r_m = rk2_2d(G,Msun,Pdmercury,Pvmercury,Tmercury,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tmercury/N)/86400, Y_r_m(1,5)/1000, Y_r_m(1,6)/1000, Y_r_m(1,7)]);
T_r_m = length(Y_c_m(:,1));
for kk = 1:T_r_m/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tmercury/N)*kk*5000/86400, Y_r_m(kk*5000,5)/1000, Y_r_m(kk*5000,6)/1000, Y_r_m(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_m
    if abs(Y_r_m(jj,7)-180) < 0.0003 && Y_r_m(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_m(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_m(jj,5)/1000,[Y_r_m(jj,3)/1000 Y_r_m(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_m(jj,7));
        fprintf(' Time (days) = %f\n',(Tmercury/N)*jj/86400);
        aphe_r_m = Y_r_m(jj,6);
    end
end
fprintf(' \n other selected parameters from simulation \n ')

side_r_m = (Tmercury/N)*T_r_m/86400;
a_r_m = (aphe_r_m + Pdmercury)/2000;
c_r_m = a_r_m - (Pdmercury/1000);
e_r_m = c_r_m/a_r_m;
mean_r_m = mean(Y_r_m(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_m);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_m);
fprintf(' orbit eccentricity = %f \n',e_r_m);
fprintf(' mean orbital spped = %f km/s \n',mean_r_m);

figure
subplot(2,4,1)
plot(Y_r_m(:,1), Y_r_m(:,2),'g')
%title(' Mercury rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Mercury rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Venus
fprintf('\n Venus \n')
Y_r_v = rk2_2d(G,Msun,Pdvenus,Pvvenus,Tvenus,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tvenus/N)/86400, Y_r_v(1,5)/1000, Y_r_v(1,6)/1000, Y_r_v(1,7)]);
T_r_v = length(Y_r_v(:,1));
for kk = 1:T_r_v/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tvenus/N)*kk*5000/86400, Y_r_v(kk*5000,5)/1000, Y_r_v(kk*5000,6)/1000, Y_r_v(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_v
    if abs(Y_r_v(jj,7)-180) < 0.0003 && Y_r_v(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_v(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_v(jj,5)/1000,[Y_r_v(jj,3)/1000 Y_r_v(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_v(jj,7));
        fprintf(' Time (days) = %f\n',(Tvenus/N)*jj/86400);
        aphe_r_v = Y_r_v(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_v = (Tvenus/N)*T_r_v/86400;
a_r_v = (aphe_r_v+Pdvenus)/2000;
c_r_v = a_r_v - (Pdvenus/1000);
e_r_v = c_r_v/a_r_v;
mean_r_v = mean(Y_r_v(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_v);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_v);
fprintf(' orbit eccentricity = %f \n',e_r_v);
fprintf(' mean orbital spped = %f km/s \n',mean_r_v);

subplot(2,4,2)
plot(Y_r_v(:,1), Y_r_v(:,2),'m')
%title(' Venus rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Venus rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth
fprintf('\n Earth \n')
Y_r_e = rk2_2d(G,Msun,Pdearth,Pvearth,Tearth,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tearth/N)/86400, Y_r_e(1,5)/1000, Y_r_e(1,6)/1000, Y_r_e(1,7)]);
T_r_e = length(Y_r_e(:,1));
for kk = 1:T_r_e/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tearth/N)*kk*5000/86400, Y_r_e(kk*5000,5)/1000, Y_r_e(kk*5000,6)/1000, Y_r_e(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_e
    if abs(Y_r_e(jj,7)-180) < 0.0005 && Y_r_e(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_e(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_e(jj,5)/1000,[Y_r_e(jj,3)/1000 Y_r_e(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_e(jj,7));
        fprintf(' Time (days) = %f\n',(Tearth/N)*jj/86400);
        aphe_r_e = Y_r_e(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_e = (Tearth/N)*T_r_e/86400;
a_r_e = (aphe_r_e+Pdearth)/2000;
c_r_e = a_r_e - (Pdearth/1000);
e_r_e = c_r_e/a_r_e;
mean_r_e = mean(Y_r_e(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_e);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_e);
fprintf(' orbit eccentricity = %f \n',e_r_e);
fprintf(' mean orbital spped = %f km/s \n',mean_r_e);

subplot(2,4,3)
plot(Y_r_e(:,1), Y_r_e(:,2),'b')
%title(' Earth rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Earth rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mars
fprintf('\n Mars \n')
Y_r_ma = rk2_2d(G,Msun,Pdmars,Pvmars,Tmars,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tmars/N)/86400, Y_r_ma(1,5)/1000, Y_r_ma(1,6)/1000, Y_r_ma(1,7)]);
T_r_ma = length(Y_r_ma(:,1));
for kk = 1:T_r_ma/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tmars/N)*kk*5000/86400, Y_r_ma(kk*5000,5)/1000, Y_r_ma(kk*5000,6)/1000, Y_r_ma(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_ma
    if abs(Y_r_ma(jj,7)-180) < 0.0004 && Y_r_ma(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_ma(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_ma(jj,5)/1000,[Y_r_ma(jj,3)/1000 Y_r_ma(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_ma(jj,7));
        fprintf(' Time (days) = %f\n',(Tmars/N)*jj/86400);
        aphe_r_ma = Y_r_ma(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_ma = (Tmars/N)*T_r_ma/86400;
a_r_ma = (aphe_r_ma+Pdmars)/2000;
c_r_ma = a_r_ma - (Pdmars/1000);
e_r_ma = c_r_ma/a_r_ma;
mean_r_ma = mean(Y_r_ma(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_ma);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_ma);
fprintf(' orbit eccentricity = %f \n',e_r_ma);
fprintf(' mean orbital spped = %f km/s \n',mean_r_ma);

subplot(2,4,4)
plot(Y_r_ma(:,1), Y_r_ma(:,2),'r')
%title(' Mars rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Mars rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jupiter
fprintf('\n Jupiter \n')
Y_r_j = rk2_2d(G,Msun,Pdjupiter,Pvjupiter,Tjupiter,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',1, ...
        [(Tjupiter/N)/86400, Y_r_j(1,5)/1000, Y_r_j(1,6)/1000, Y_r_j(1,7)]);
T_r_j = length(Y_r_j(:,1));
for kk = 1:T_r_j/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tjupiter/N)*kk*5000/86400, Y_r_j(kk*5000,5)/1000, Y_r_j(kk*5000,6)/1000, Y_r_j(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_j
    if abs(Y_r_j(jj,7)-180) < 0.00035 && Y_r_j(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_j(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_j(jj,5)/1000,[Y_r_j(jj,3)/1000 Y_r_j(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_j(jj,7));
        fprintf(' Time (days) = %f\n',(Tjupiter/N)*jj/86400);
        aphe_r_j = Y_r_j(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_j = (Tjupiter/N)*T_r_j/86400;
a_r_j = (aphe_r_j+Pdjupiter)/2000;
c_r_j = a_r_j - (Pdjupiter/1000);
e_r_j = c_r_j/a_r_j;
mean_r_j = mean(Y_r_j(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_j);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_j);
fprintf(' orbit eccentricity = %f \n',e_r_j);
fprintf(' mean orbital spped = %f km/s \n',mean_r_j);

subplot(2,4,5)
plot(Y_r_j(:,1), Y_r_j(:,2),'k')
%title(' Jupiter rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Jupiter rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saturn
fprintf('\n Saturn \n')
Y_r_s = rk2_2d(G,Msun,Pdsaturn,Pvsaturn,Tsaturn,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Tsaturn/N)/86400, Y_r_s(1,5)/1000, Y_r_s(1,6)/1000, Y_r_s(1,7)]);
T_r_s = length(Y_r_s(:,1));
for kk = 1:T_r_s/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',kk*5000, ...
        [(Tsaturn/N)*kk*5000/86400, Y_r_s(kk*5000,5)/1000, Y_r_s(kk*5000,6)/1000, Y_r_s(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_s
    if abs(Y_r_s(jj,7)-180) < 0.0004 && Y_r_s(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_s(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_s(jj,5)/1000,[Y_r_s(jj,3)/1000 Y_r_s(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_s(jj,7));
        fprintf(' Time (days) = %f\n',(Tsaturn/N)*jj/86400);
        aphe_r_s = Y_r_s(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_s = (Tsaturn/N)*T_r_s/86400;
a_r_s = (aphe_r_s+Pdsaturn)/2000;
c_r_s = a_r_s - (Pdsaturn/1000);
e_r_s = c_r_s/a_r_s;
mean_r_s = mean(Y_r_s(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_s);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_s);
fprintf(' orbit eccentricity = %f \n',e_r_s);
fprintf(' mean orbital spped = %f km/s \n',mean_r_s);

subplot(2,4,6)
plot(Y_r_s(:,1), Y_r_s(:,2),'y')
%title(' Saturn rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Saturn rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uranus
fprintf('\n Uranus \n')
Y_r_u = rk2_2d(G,Msun,Pduranus,Pvuranus,Turanus,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f  \n',1, ...
        [(Turanus/N)/86400, Y_r_u(1,5)/1000, Y_r_u(1,6)/1000, Y_r_u(1,7)]);
T_r_u = length(Y_r_u(:,1));
for kk = 1:T_r_u/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Turanus/N)*kk*5000/86400, Y_r_u(kk*5000,5)/1000, Y_r_u(kk*5000,6)/1000, Y_r_u(kk*5000,7) ]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_u
    if abs(Y_r_u(jj,7)-180) < 0.0004 && Y_r_u(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_u(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_u(jj,5)/1000,[Y_r_u(jj,3)/1000 Y_r_u(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_u(jj,7));
        fprintf(' Time (days) = %f\n',(Turanus/N)*jj/86400);
        aphe_r_u = Y_r_u(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_u = (Turanus/N)*T_r_u/86400;
a_r_u = (aphe_r_u+Pduranus)/2000;
c_r_u = a_r_u - (Pduranus/1000);
e_r_u = c_r_u/a_r_u;
mean_r_u = mean(Y_r_u(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_u);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_u);
fprintf(' orbit eccentricity = %f \n',e_r_u);
fprintf(' mean orbital spped = %f km/s \n',mean_r_u);

subplot(2,4,7)
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_r_u(:,1), Y_r_u(:,2),'Color',color)
%title(' Uranus rk2 method ')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Uranus rk2 method'])
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neptune
fprintf('\n Neptune \n')
Y_r_n = rk2_2d(G,Msun,Pdneptune,Pvneptune,Tneptune,N);
fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f   \n',1, ...
        [(Tneptune/N)/86400, Y_r_n(1,5)/1000, Y_r_n(1,6)/1000, Y_r_n(1,7)]);
T_r_n = length(Y_r_n(:,1));
for kk = 1:T_r_n/5000
    fprintf(' step=%i, time[day]=%f, v[km/s]=%f, dist2Sun[km]=%e, theta=%f \n',kk*5000, ...
        [(Tneptune/N)*kk*5000/86400, Y_r_n(kk*5000,5)/1000, Y_r_n(kk*5000,6)/1000, Y_r_n(kk*5000,7)]);
end

fprintf('\n Aphelion parameters calculated \n');

for jj = 1:T_r_n
    if abs(Y_r_n(jj,7)-180) < 0.0005 && Y_r_n(jj,7)>180
        fprintf(' Distance to the Sun = %e km\n',Y_r_n(jj,6)/1000);
        fprintf(' Speed = %f km/s, [%f %f] \n',Y_r_n(jj,5)/1000,[Y_r_n(jj,3)/1000 Y_r_n(jj,4)/1000]);
        fprintf(' True Anomaly = %f degrees\n',Y_r_n(jj,7));
        fprintf(' Time (days) = %f\n',(Tneptune/N)*jj/86400);
        aphe_r_n = Y_r_n(jj,6);
    end
end

fprintf(' \n other selected parameters from simulation \n ')

side_r_n = (Tneptune/N)*T_r_n/86400;
a_r_n = (aphe_r_n+Pdneptune)/2000;
c_r_n = a_r_n - (Pdneptune/1000);
e_r_n = c_r_n/a_r_n;
mean_r_n = mean(Y_r_n(:,5))/1000;

fprintf(' Semimajor axis = %e km\n',a_r_n);
fprintf(' Sidereal orbit period (days) = %f \n',side_r_n);
fprintf(' orbit eccentricity = %f \n',e_r_n);
fprintf(' mean orbital spped = %f km/s \n',mean_r_n);

subplot(2,4,8)
plot(Y_r_n(:,1), Y_r_n(:,2),'c')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel(['position in x-axis',newline,'\bf Neptune rk2 method'])
ylabel('position in y-axis')

sgtitle('Planet orbits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all terrestrial planets rk2 method
figure
plot(Y_r_m(:,1), Y_r_m(:,2),'g')
hold on
plot(Y_r_v(:,1), Y_r_v(:,2),'m')
plot(Y_r_e(:,1), Y_r_e(:,2),'b')
plot(Y_r_ma(:,1), Y_r_ma(:,2),'r')
hold off
title(' all terrestrial planets orbit using rk2 method ')
legend({'Mercury','Venus','Earth','Mars'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e11 5e11])
ylim([-5e11 5e11])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all Jovian planets rk2 method
figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_r_j(:,1), Y_r_j(:,2),'k')
hold on
plot(Y_r_s(:,1), Y_r_s(:,2),'y')
plot(Y_r_u(:,1), Y_r_u(:,2),'Color',color)
plot(Y_r_n(:,1), Y_r_n(:,2),'c')
hold off
title(' all jovian planets orbit using rk2 method ')
legend({'Jupiter','Saturn','Uranus','Neptune'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all planets rk2 method
figure
str = '#FF0000';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
plot(Y_r_m(:,1), Y_r_m(:,2),'g')
hold on
plot(Y_r_v(:,1), Y_r_v(:,2),'m')
plot(Y_r_e(:,1), Y_r_e(:,2),'b')
plot(Y_r_ma(:,1), Y_r_ma(:,2),'r')
plot(Y_r_j(:,1), Y_r_j(:,2),'k')
plot(Y_r_s(:,1), Y_r_s(:,2),'y')
plot(Y_r_u(:,1), Y_r_u(:,2),'Color',color)
plot(Y_r_n(:,1), Y_r_n(:,2),'c')
hold off
title(' all planets orbit using rk2 method ')
legend({'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune'},'Location','bestoutside')
pbaspect([1 1 1])
xlim([-5e12 5e12])
ylim([-5e12 5e12])
xlabel('position in x-axis')
ylabel('position in y-axis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n####################################################################\n')

toc
%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = euler_2d(G,Ms,r0,v0,T,N)
rr0 = zeros(2,1);
vv0 = zeros(2,1);
rr0(1,1) = r0;
vv0(2,1) = v0;
dt = T/N;
V = zeros(N,1);
R = zeros(N,1);
theta = zeros(N,1);

r = zeros(2,N);
v = zeros(2,N);

r(:,1) = rr0;
v(:,1) = vv0;
theta(1) = atand(rr0(2,1)/rr0(1,1));
R(1) = norm(rr0);
V(1) = norm(vv0);
for i = 1:N-1
    r(:,i+1) = r(:,i) + v(:,i)*dt;
    if (r(2,i+1)/r(1,i+1))>0 && r(1,i+1)>0
        theta(i+1) = atand(r(2,i+1)/r(1,i+1));
    elseif (r(2,i+1)/r(1,i+1))<0 && r(1,i+1)<0
        theta(i+1) = 180 + atand(r(2,i+1)/r(1,i+1)); 
    elseif (r(2,i+1)/r(1,i+1))>0 && r(1,i+1)<0
        theta(i+1) = 180 + atand(r(2,i+1)/r(1,i+1)); 
    elseif (r(2,i+1)/r(1,i+1))<0 && r(1,i+1)>0
        theta(i+1) = 360 + atand(r(2,i+1)/r(1,i+1));
    end
    v(:,i+1) = v(:,i) - (G*Ms/R(i)^3)*r(:,i)*dt;
    R(i+1) = norm(r(:,i+1));
    V(i+1) = norm(v(:,i+1));
end

Y = [r(1,:)' r(2,:)' v(1,:)' v(2,:)' V R theta];
end

function Y = cromer_2d(G,Ms,r0,v0,T,N)
rr0 = zeros(2,1);
vv0 = zeros(2,1);
rr0(1,1) = r0;
vv0(2,1) = v0;
dt = T/N;
V = zeros(N,1);
R = zeros(N,1);
theta = zeros(N,1);

r = zeros(2,N);
v = zeros(2,N);

r(:,1) = rr0;
v(:,1) = vv0;
theta(1) = atand(rr0(2,1)/rr0(1,1));
R(1) = norm(rr0);
V(1) = norm(vv0);
for i = 1:N-1
    v(:,i+1) = v(:,i) - (G*Ms/R(i)^3)*r(:,i)*dt;
    r(:,i+1) = r(:,i) + v(:,i+1)*dt;
    if (r(2,i+1)/r(1,i+1))>0 && r(1,i+1)>0
        theta(i+1) = atand(r(2,i+1)/r(1,i+1));
    elseif (r(2,i+1)/r(1,i+1))<0 && r(1,i+1)<0
        theta(i+1) = 180 + atand(r(2,i+1)/r(1,i+1)); 
    elseif (r(2,i+1)/r(1,i+1))>0 && r(1,i+1)<0
        theta(i+1) = 180 + atand(r(2,i+1)/r(1,i+1)); 
    elseif (r(2,i+1)/r(1,i+1))<0 && r(1,i+1)>0
        theta(i+1) = 360 + atand(r(2,i+1)/r(1,i+1));
    end
    R(i+1) = norm(r(:,i+1));
    V(i+1) = norm(v(:,i+1));
end

Y = [r(1,:)' r(2,:)' v(1,:)' v(2,:)' V R theta];
end

function Y = rk2_2d(G,Ms,r0,v0,T,N)
rr0 = zeros(2,1);
vv0 = zeros(2,1);
rr0(1,1) = r0;
vv0(2,1) = v0;
dt = T/N;
V = zeros(N,1);
R = zeros(N,1);
theta = zeros(N,1);

r = zeros(2,N);
v = zeros(2,N);

r(:,1) = rr0;
v(:,1) = vv0;
theta(1) = atand(rr0(2,1)/rr0(1,1));
R(1) = norm(rr0);
V(1) = norm(vv0);
for i = 1:N-1
    r_old = r(:,i) + v(:,i)*dt/2;
    v_old = v(:,i) - (G*Ms/R(i)^3)*r(:,i)*dt/2;
    r(:,i+1) = r(:,i) + v_old*dt;
    if (r(2,i+1)/r(1,i+1))>0 && r(1,i+1)>0
        theta(i+1) = atand(r(2,i+1)/r(1,i+1));
    elseif (r(2,i+1)/r(1,i+1))<0 && r(1,i+1)<0
        theta(i+1) = 180 + atand(r(2,i+1)/r(1,i+1)); 
    elseif (r(2,i+1)/r(1,i+1))>0 && r(1,i+1)<0
        theta(i+1) = 180 + atand(r(2,i+1)/r(1,i+1)); 
    elseif (r(2,i+1)/r(1,i+1))<0 && r(1,i+1)>0
        theta(i+1) = 360 + atand(r(2,i+1)/r(1,i+1));
    end
    R(i+1) = norm(r(:,i+1));
    v(:,i+1) = v(:,i) - (G*Ms/norm(r_old)^3)*r_old*dt;
    V(i+1) = norm(v(:,i+1));
end

Y = [r(1,:)' r(2,:)' v(1,:)' v(2,:)' V R theta];
end

  
























