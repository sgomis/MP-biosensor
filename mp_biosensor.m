
% Title: Model for reagentless biosensor based on the dynamics of an inverted pendulum
% Author: Surath Gomis 
% Research Article: Reagentless Biomolecular Analysis Using a Molecular Pendulum 
% Updated: April 19, 2020

%% Generate potential field 
clear;clf;
domain = ones(5e3,5e3); %larger domains take too much data
dy = 1e-11; %each index is 0.01nm
dy_domain = domain*dy; %domain size in real units

C = 0.0137; %137 mM salt for 1x PBS %150 mM for pbs or 15mM for 0.1x pbs; or 0.5mM total conc adding pos and neg ions
Na = 6.022e23;
n = C*Na/0.001; %convert mol/L to #/m^3
z = 1; %charge of ion (assumes single ion solution, so we have to estimate)
er = 80; %permittivity of water
e0 = 8.85e-12; %permittivity of free space
q = 1.602e-19; %electronic charge
k = 1.38e-23; %boltzmann constant
T = 300; %temp
kappa_full = (2*n*(z^2)*(q^2)/(er*e0*k*T))^0.5; %full formulation
% kappa = 1e2*(3.29e7)*z*C^0.5; %simplified expression for aqueous sol, converted to /m
V0 = 0.5; %applied potential

V = zeros(length(domain),1);
for i = 1:1:length(domain)
    V(i) = (4*k*T/z/q)*atanh(tanh(z*q*V0/4/k/T)*exp(-kappa_full*i*dy)); %solve for voltage in domain, pulled from Electrochemical Methods - Fundamentals and Applications (Bard 2001)
end 


figure(1)
plot(1e9*(dy:dy:length(domain)*dy),V,'k','LineWidth',5)
xlabel('distance from electrode [nm]','fontsize',18)
ylabel('magnitude of potential field [V]','fontsize',18)
% xlim([0 20])
set(gca,'FontSize',24,'LineWidth',3,'FontName','Arial','units','inches','position',[1.5 1.5 5 5])
xlim([0 50])
%% Derive electric field
Ey = V*0;
%Calulate Ey field using the central difference method
for i = 2:length(V)-1 
        Ey(i) = -(V(i+1)-V(i-1))/(2*dy);
end

%set first and last values missed in central diff 
Ey(1) = Ey(2);
Ey(length(Ey)) = Ey(length(Ey)-1);

figure(2)
clf
plot(1e9*(dy:dy:length(domain)*dy),Ey,'k','LineWidth',5)
xlabel('distance from electrode [nm]','fontsize',18)
ylabel('magnitude of electric field [V/m]','fontsize',18)
% xlim([0 20])
set(gca,'FontSize',24,'LineWidth',3,'FontName','Arial','units','inches','position',[1.5 1.5 5 5])
xlim([0 50])
% axis([0 5 0 1e9])
%%
kd = 1.66054e-24; %kg per kilodalton
bp = 0.65*kd; %mass of 1 base pair
elec = 1.602e-19; %electronic charge
num_dna = 24; %dna linker length
q_multiplier = 100; %linear multiplier for dna charge to speed up ODE solving
q_dna = num_dna*2*q_multiplier; %dna linker charge 
L = num_dna*(0.34e-9) + 6e-9; %length of DNA to middle of Ab (Ab = 10nm, 150kD)
e0 = 8.854e-12; %permittivity of free space
er = 1; %medium permittivity, linearly scales output, set to 1 to speed up computation
eta = 8.9e-4; %dynamic viscosity of dionized water [Pa*s]
kc = 8.99e9; %coulombs constant
dist = 8e-9; %distance between dna linkers, i.e. probe density 
qr_i = 0; %value to decay charge if strand hits surface before middle strand
ql_i = 0; %value to decay charge if strand hits surface before middle strand

%Experimental parameters of system - e.g. CTI protein + Ab
R = 0.34e-9 + 5e-9 + 2.01e-9;
q0 = -((q_dna-3)/er)*elec; 
qr = -((q_dna-3)/er)*elec; 
ql = -((q_dna-3)/er)*elec; 
m = num_dna*bp + 150*kd + 24*kd;

%% Generate Gaussian distribution to sample starting angles (via inverse transform sampling)
figure(3)
clf
a = [-90:0.001:90];
% r = normpdf(a,30,40);
r = normpdf(a,0,15); %starting angles centered around 0, sd of 15

cdf = cumtrapz(r);
cdf = cdf/max(cdf);
uniform = rand(1000000,1);
Y = interp1(cdf,a,uniform);
histogram(Y)
% 
% clf
% a_elec = [-10e-18:1e-20:10e-18];
% r_elec = normpdf(a_elec,-4e-18,1e-17);
% plot(r_elec);
% 
% cdf_elec = cumtrapz(r_elec);
% cdf_elec = cdf_elec/max(cdf_elec);
% plot(cdf_elec);
% 
% uniform_elec = rand(1000000,1);
% Y_elec = interp1(cdf_elec,a_elec,uniform_elec);
% histogram(Y_elec)
    
%% Perform 4th-order Runge Kutta to solve EOM 

clearvars TTF    

for sample = 1:5000 %solve for 5000 dna strands %actual numDNA per electrode = 3.5e-4cm^2/66nm^2 = 5.3e8

    clearvars x0 y0 y0_cell y0_cell_old theta0 omega0 xr yr yr_cell yr_cell_old thetar omegar xl yl yl_cell yl_cell_old thetal omegal; 
    
    %reset charges
    qr_i = 0;
    ql_i = 0;
    
    %set the time interval
    T0 = 0;
    Tf = 1e-5; 
    T = 1e-12;
    t = T0:T:Tf;

    %preallocate matrices for the variables to solve for
    numPoints = length(t);
       
    %define initial conditions
    theta0_init = Y(sample)*pi/180;
    omega0_init = 0;
    thetar_init = Y(sample)*pi/180;
    omegar_init = 0;
    thetal_init = Y(sample)*pi/180;
    omegal_init = 0;
    
    %include initial conditions
    theta0 = theta0_init;
    omega0 = omega0_init;
    thetar = thetar_init;
    omegar = omegar_init;
    thetal = thetal_init;
    omegal = omegal_init;

    %pre-allocate arrays
    x0 = zeros(length(numPoints)-1,1);   
    vx0 = zeros(length(numPoints)-1,1);
    y0 = zeros(length(numPoints)-1,1); 
    vy0 = zeros(length(numPoints)-1,1);

    xr = zeros(length(numPoints)-1,1);  
    vxr = zeros(length(numPoints)-1,1);
    yr = zeros(length(numPoints)-1,1);
    vyr = zeros(length(numPoints)-1,1);

    xl = zeros(length(numPoints)-1,1);  
    vxl = zeros(length(numPoints)-1,1);
    yl = zeros(length(numPoints)-1,1);
    vyl = zeros(length(numPoints)-1,1);
    
    %iterate to calculate all k values to evaluate the next point for all four
    %variables (3 thetas and 3 omegas) we are solving for, using the RK4 method 
    for i = 1:1:numPoints-1

        %set current thetas and omegas
        theta0_t = theta0;
        omega0_t = omega0;
        thetar_t = thetar;
        omegar_t = omegar;
        thetal_t = thetal;
        omegal_t = omegal;

        %calculate the (x,y), (vx,vy) true coordinates
        x0(i) = L*sin(theta0);   
        vx0(i) = L*cos(theta0)*omega0;
        y0(i) = L*cos(theta0); 
        vy0(i) = -L*sin(theta0)*omega0;
        
        xr(i) = L*sin(thetar) + dist;   
        vxr(i) = L*cos(thetar)*omegar;
        yr(i) = L*cos(thetar); 
        vyr(i) = -L*sin(thetar)*omegar;
        
        xl(i) = L*sin(thetal) - dist;   
        vxl(i) = L*cos(thetal)*omegal;
        yl(i) = L*cos(thetal); 
        vyl(i) = -L*sin(thetal)*omegal;
   
        %condition to stop simulation of the middle DNA hits the electrode
        if y0(i) <= 0
            time_steps = i;
            time_real = i*T*80*q_multiplier; %correct for medium permittivity (er=80) and dna charge multiplier,
                                             %which were used to speed computation 
            disp('The protein hit the electrode');
            fprintf('Time steps = %d\n',i)
            fprintf('Time [us] = %d\n',time_real*1e6)
            fprintf('Sample = %d\n',sample)
            break;
        end

        %linearly interpolate (x_cell,y_cell) from (x,y) to get domain
        %values - used to maintain right and left strand final positions if they
        %hit the electrode first
        y0_cell = round(1 + (y0(i)-0)*((length(domain)-1))/(dy*((length(domain)-1)-0)));
        yr_cell = round(1 + (yr(i)-0)*((length(domain)-1))/(dy*((length(domain)-1)-0)));
        yl_cell = round(1 + (yl(i)-0)*((length(domain)-1))/(dy*((length(domain)-1)-0)));
        
        %slowly decay right strand charge over time if strand hits before middle strand
        if yr(i) <= 0 
            yr_cell = yr_cell_old;
            yr(i) = yr(i-1);
            xr(i) = xr(i-1);
            qr = (-num_dna+qr_i)*elec;
            qr_i = qr_i + 0; %parameter to tweak if this case is met (can happen when using large starting angles)
        end
        
        %slowly decay left strand charge over time if strand hits before middle strand
        if yl(i) <= 0 
            yl_cell = yl_cell_old;
            yl(i) = yl(i-1);
            xl(i) = xl(i-1);
            ql = (-num_dna+ql_i)*elec;
            ql_i = ql_i + 0; %parameter to tweak if this case is met (can happen when using large starting angles)
        end
        
        %Lorentz force 
        F0_e = -q0*Ey(y0_cell)*sin(theta0);
        Fr_e = -qr*Ey(yr_cell)*sin(thetar);
        Fl_e = -ql*Ey(yl_cell)*sin(thetal);
        
        %Stoke drag force
        F0_d = -6*pi*eta*R*L*omega0;
        Fr_d = -6*pi*eta*R*L*omegar;
        Fl_d = -6*pi*eta*R*L*omegal;

        %Coulomb force
        F0_c = kc*q0*qr*(L*sin(theta0-thetar)-dist*cos(theta0)) / (2*L^2 + dist^2 + 2*dist*L*(sin(thetar)-sin(theta0)) - 2*(L^2)*cos(thetar-theta0))^1.5 ...
            +  kc*q0*ql*(L*sin(theta0-thetal)+dist*cos(theta0)) / (2*L^2 + dist^2 - 2*dist*L*(sin(thetal)-sin(theta0)) - 2*(L^2)*cos(thetal-theta0))^1.5; %avg probe density: 1 molecule/30nm^2
        
        Fr_c = kc*qr*ql*(L*sin(thetar-thetal)-dist*cos(thetar)) / (2*L^2 + dist^2 + 2*dist*L*(sin(thetal)-sin(thetar)) - 2*(L^2)*cos(thetal-thetar))^1.5 ...
            +  kc*qr*q0*(L*sin(thetar-theta0)+dist*cos(thetar)) / (2*L^2 + dist^2 - 2*dist*L*(sin(theta0)-sin(thetar)) - 2*(L^2)*cos(theta0-thetar))^1.5;
        
        Fl_c = kc*ql*q0*(L*sin(thetal-theta0)-dist*cos(thetal)) / (2*L^2 + dist^2 + 2*dist*L*(sin(theta0)-sin(thetal)) - 2*(L^2)*cos(theta0-thetal))^1.5 ...
            +  kc*ql*qr*(L*sin(thetal-thetar)+dist*cos(thetal)) / (2*L^2 + dist^2 - 2*dist*L*(sin(thetar)-sin(thetal)) - 2*(L^2)*cos(thetar-thetal))^1.5;
        
        
        %find first ks (derivative functions for omegas are in separate
        %function)
        k1_theta0 = omega0_t;
        k1_omega0 = (1/m/L)*(F0_d + F0_e + F0_c);
        k1_thetar = omegar_t;
        k1_omegar = (1/m/L)*(Fr_d + Fr_e + Fr_c);
        k1_thetal = omegal_t;
        k1_omegal = (1/m/L)*(Fl_d + Fl_e + Fl_c);

        %solve for next ks using previous ks, with time multiplying as per RK4
        theta0_t = theta0 + k1_theta0*T/2;
        omega0_t = omega0 + k1_omega0*T/2;
        thetar_t = thetar + k1_thetar*T/2;
        omegar_t = omegar + k1_omegar*T/2;
        thetal_t = thetal + k1_thetal*T/2;
        omegal_t = omegal + k1_omegal*T/2;

        %evaluate for the second ks
        k2_theta0 = omega0_t;
        k2_omega0 = (1/m/L)*(F0_d + F0_e + F0_c);
        k2_thetar = omegar_t;
        k2_omegar = (1/m/L)*(Fr_d + Fr_e + Fr_c);
        k2_thetal = omegal_t;
        k2_omegal = (1/m/L)*(Fl_d + Fl_e + Fl_c);

        %again, use previous ks to find arguments for next ks
        theta0_t = theta0 + k2_theta0*T/2;
        omega0_t = omega0 + k2_omega0*T/2;
        thetar_t = thetar + k2_thetar*T/2;
        omegar_t = omegar + k2_omegar*T/2;
        thetal_t = thetal + k2_thetal*T/2;
        omegal_t = omegal + k2_omegal*T/2;

        %evalute for third ks
        k3_theta0 = omega0_t;
        k3_omega0 = (1/m/L)*(F0_d + F0_e + F0_c);
        k3_thetar = omegar_t;
        k3_omegar = (1/m/L)*(Fr_d + Fr_e + Fr_c);
        k3_thetal = omegal_t;
        k3_omegal = (1/m/L)*(Fl_d + Fl_e + Fl_c);

        %perform last evaluation of parameters for next ks
        theta0_t = theta0 + k3_theta0*T;
        omega0_t = omega0 + k3_omega0*T;
        thetar_t = thetar + k3_thetar*T;
        omegar_t = omegar + k3_omegar*T;
        thetal_t = thetal + k3_thetal*T;
        omegal_t = omegal + k3_omegal*T;

        %evaluate for fourth ks
        k4_theta0 = omega0_t;
        k4_omega0 = (1/m/L)*(F0_d + F0_e + F0_c);
        k4_thetar = omegar_t;
        k4_omegar = (1/m/L)*(Fr_d + Fr_e + Fr_c);
        k4_thetal = omegal_t;
        k4_omegal = (1/m/L)*(Fl_d + Fl_e + Fl_c);

        %evalute next thetas and omegas points using all four ks, as per RK4
        theta0 = theta0 + 1/6*(k1_theta0 + 2*k2_theta0 + 2*k3_theta0 + k4_theta0)*T;
        omega0 = omega0 + 1/6*(k1_omega0 + 2*k2_omega0 + 2*k3_omega0 + k4_omega0)*T;   
        thetar = thetar + 1/6*(k1_thetar + 2*k2_thetar + 2*k3_thetar + k4_thetar)*T;
        omegar = omegar + 1/6*(k1_omegar + 2*k2_omegar + 2*k3_omegar + k4_omegar)*T;   
        thetal = thetal + 1/6*(k1_thetal + 2*k2_thetal + 2*k3_thetal + k4_thetal)*T;
        omegal = omegal + 1/6*(k1_omegal + 2*k2_omegal + 2*k3_omegal + k4_omegal)*T;   
        
        %save last cell value if needed to store l and r last positions
        yl_cell_old = yl_cell;
        yr_cell_old = yr_cell;

    end
    
    %save time to fall if simulating multiple samples
    TTF(sample) = time_real;
    TTF_steps = time_steps;
end

%% Plot motion of DNA strands (for single test case; set sample = 1 in previous block)
figure(4)
for i = 1:100:length(y0)-1 %choose how fast you plot by changing middle value
    clf
    hold all
    plot([0,x0(i)],[0,y0(i)],'LineWidth',3);
    plot([x0(i)],[y0(i)],'-o','MarkerSize',15,'markerfacecolor','r','LineWidth',3);
    plot([dist,xr(i)],[0,yr(i)],'LineWidth',3);
    plot([xr(i)],[yr(i)],'-o','MarkerSize',15,'markerfacecolor','r','LineWidth',3);
    plot([-dist,xl(i)],[0,yl(i)],'LineWidth',3);
    plot([xl(i)],[yl(i)],'-o','MarkerSize',15,'markerfacecolor','r','LineWidth',3);
    axis([-1e-8-dist 1e-8+dist 0 2e-8])
    set(gca,'DataAspectRatio',[1 1 1])
    drawnow %limitrate
    i
end
drawnow
