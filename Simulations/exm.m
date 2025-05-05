T = @(theta) [ ...
        cos(theta)^2,          sin(theta)^2,          2*sin(theta)*cos(theta); ...
        sin(theta)^2,          cos(theta)^2,         -2*sin(theta)*cos(theta); ...
       -sin(theta)*cos(theta), sin(theta)*cos(theta), cos(theta)^2 - sin(theta)^2 ...
    ];

syms t q Nxy
Q = [100 10 0; 10 20 0; 0 0 20]
qbar1 = T(deg2rad(45))^-1 * Q * (T(deg2rad(45))^-1)'

qbar2 = T(deg2rad(-45))^-1 * Q * (T(deg2rad(-45))^-1)'

T(deg2rad(45))
T(deg2rad(-45))

s45 = [Nxy/(5*t);Nxy/(5*t); Nxy/(4*t) ]
sn45 = [Nxy/(-5*t);Nxy/(-5*t); Nxy/(4*t) ]

sprime45 = T(deg2rad(45)) * s45
sprimen45 = T(deg2rad(-45)) * sn45