function [r_e,r_v,v_e,v_v] = EarthVenusPosVel(t0,dt)

%INPUTS
%   t0 - Initial time (JD)
%   dt - Cruise duration (days)
%
%OUTPUTS
%   r_e - 3x1 col vector - Heliocentric position of Earth (in AU) at t0
%   r_v - 3x1 col vector - Heliocentric position of Venus (in AU) at t0+dt
%   v_e - 3x1 col vector - Heliocentric-Inertial velocity of Earth 
%                          (in AU/day) at t0
%   v_v - 3x1 col vector - Heliocentric-Inertial velocity of Venus 
%                          (in AU/day) at t0+dt

%Earth and Venus orbital elements
a_e = 1.000373836656026E+00;        %sma (AU)
a_v = 7.233304294142561E-01;
e_e = 1.712127710968187E-02;        %eccentricity
e_v = 6.801883435239245E-03;
I_e = 2.777040607882003E-03*pi/180; %Inclination (rad)
I_v = 3.394490723804590E+00*pi/180; 
w_e = 3.043573249748720E+02*pi/180; %arg. of periapsis (rad)
w_v = 5.476413547432202E+01*pi/180;
O_e = 1.596967974767415E+02*pi/180; %long. of ascending node (rad)
O_v = 7.662845580893872E+01*pi/180;
t_p_e = 2458853.731945450883;       % time of periapsis passage (JD)
t_p_v = 2458029.748282369226; 
mu_e = 2.9591309705483544E-04;      %G(m_sun + m_planet) (AU^3/day^2)
mu_v = 2.9591293263082414E-04;

[r_e,v_e] = rfromorbes(mu_e,t0,t_p_e,a_e,e_e,O_e,I_e,w_e);
[r_v,v_v] = rfromorbes(mu_v,t0+dt,t_p_v,a_v,e_v,O_v,I_v,w_v);

    function [ri,vi] = rfromorbes(mu,t0,tp,a,e,Omega,I,omega)
        n = sqrt(mu/(a^3));
        t = t0 - tp;
        M = n*t;
        M = mod(M,2*pi);
        
        %Newton-Raphson to find Eccentric anomaly
        tol = eps(2*pi);
        if M/(1-e) <sqrt(6*(1-e)/e)
            E0=M/(1-e);
        else
            E0=(6*M./e).^(1/3);
        end
        E=E0-(M-E0+e*sin(E0))/(e*cos(E0)-1);
        while abs(E-E0)>tol
            E0=E;
            E=E0-(M-E0+e*sin(E0))/(e*cos(E0)-1);
        end
        
        %orbital radius and velocity in the perifocal frame:
        rp = [a.*(cos(E)-e); a*sqrt(1-e.^2).*sin(E); 0];
        rmag = a*(1 - e*cos(E));
        vp = a*n/rmag*[-a*sin(E); a*sqrt(1-e.^2).*cos(E); 0];

        
        PcI = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1]*...
              [1 0 0; 0 cos(I) sin(I); 0 -sin(I) cos(I)]*...
              [cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1];
        ri = PcI.'*rp;
        vi = PcI.'*vp;
    end

%ensure outputs are column vectors
r_e = r_e(:);
r_v = r_v(:);
v_e = v_e(:);
v_v = v_v(:);
end