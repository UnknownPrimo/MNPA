%% ELEC 4700 
%% PA 7 - MNA Building

clf;
    
% Variables      
% Resistors
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000; 
Cap = 0.25;
L = 0.2;
vs = 1;
a = 100;

% Conductance
G1 = 1/R1;  
G2 = 1/R2;  
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;

% Part (a) creating G and C matrices
% G matrix
G = [G1 -G1 0 0 0 1 0 0;
    -G1 G1+G2 0 0 0 0 1 0;
    0 0 G3 0 0 0 -1 0 ;
    0 0 0 G4 -G4 0 0 1 ;
    0 0 0 -G4 G4+G0 0 0 0;
    1 0 0 0 0 0 0 0;
    0 1 -1 0 0 0 0 0;
    0 0 -a/R3 1 0 0 0 0;] ;

% C matrix
C = [Cap -Cap 0 0 0 0 0 0 ;
    -Cap Cap 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 -1 0;
    0 0 0 0 0 0 0 0;];

% F matrix
F = zeros(8,1);

% Part (b). For the DC case sweep the input voltage V1 from -10V to 10V 
% and plot VO and the voltage at V3.

VI = zeros(100, 1);
V01 = zeros(100, 1);
V3 = zeros(100, 1);

for Vo = -10:0.1:10
    %DC Sweep
    V2 = G\F;
    VI(v) = Vo;
    V01(v) = V2(5);
    V3(v) = V2(3);
    vs = vs + 1;
end 

figure(1)
hold on;
plot(VI, V01);
plot(VI, V3);
legend('VO','V3');
title('Plot of DC sweep at V0 and V3');
ylabel('Voltage (V)');
xlabel('Vin (V)');

% Part (c). For the AC case plot VO as a function of omega 
% also plot the gain VO/V1 in dB.

% Calculating and plotting AC sweep and Gain(dB)
omega = logspace(1,2,1000); 

for Vo = 1:length(omega)
    % AC sweep
    Vac = (G+1j*omega(Vo)*C)\F;
    % Gain   
    Gain = 20*log10(abs(Vac(5))/10);   
end
% Plot of Angular frequency (omega) as a function of VO   
figure(2)
plot(Vo,Vac);
xlabel('Angular Frequency (rad/s)')
ylabel('V0 (V)')
title('Plot of V0 - AC analysis')

% Plot of Gain in dB as a function of VO  
figure(3)
semilogx(Vo, Gain);
xlabel('Angular Frequency (rad/s)')
ylabel('Gain (dB)')
title('Plot of Gain - AC analysis ');






