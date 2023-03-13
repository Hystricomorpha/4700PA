clear;
close all;

R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R5 = 1000;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G5 = 1/R5;

C = 0.25;

L = 0.2;

a = 100;

G = [1       0       0       0   0       0 0 0;
     G1      -G1     0       0   0       0 1 0;
     -G1     G1+G2   0       0   0       1 0 0;
     0       0       G3      0   0       1 0 0;
     0       0       0       G4  -G4     0 0 1;
     0       0       0       -G4 G4+G5   0 0 0;
     0       1       -1      0   0       0 0 0;
     0       0       a*G3    -1  0       0 0 0]; % G matrix

C = [1  0   0 0 0 0  0 0;
     C  -C  0 0 0 0  0 0;
     -C C   0 0 0 0  0 0;
     0  0   0 0 0 0  0 0;
     0  0   0 0 0 0  0 0;
     0  0   0 0 0 0  0 0;
     0  0   0 0 0 -L 0 0;
     0  0   0 0 0 0  0 0]; % C matrix

% Matrix in the form:
% x(s) = [V1; V2; V3; V4; V5; I_L; I_S1; I_S2];

F = [1; 0; 0; 0; 0; 0; 0; 0];
 
% DC Case (no C matrix):
% Sweep the input voltage V1 from -10V to 10V
% Plot VO (V5) and the voltage at V3
i = 1;
for V1 = -10:10
    
    % A*x = b => x = A^-1*b
    x = G\(F*V1);

    V3(i) = x(3);
    V5(i) = x(5);

    i = i + 1;
end

Vin = -10:10;
subplot(2, 3, 1);
plot(Vin, V3); hold on;
plot(Vin, V5);
xlabel('V_{in} (V)');
ylabel('Voltage at Node (V)');
legend('V_3', 'V_{out}');
hold off;

% AC case:
% plot VO as a function of ω 
% also plot the gain VO/V1 in dB
i = 1;
for w1 = 0:100

    % Setup of A matrix A*V = F(w)
    A = G + 1j*w1*C;

    x = A\F;

    V_out(i) = x(5);
    Gain(i) = x(5)/x(1);

    i = i + 1;
end

w = 0:100;
subplot(2, 2, 2);
plot(w, V_out); 
xlabel('\omega (Hz)');
ylabel('V_{out} (V)');

subplot(2, 2, 3);
plot(w, Gain);
xlabel('\omega (Hz)');
ylabel('Gain');

% AC case:
% plot the gain as function of random perturbations 
% on C using a normal distribution with 
% std = .05 at ω = π. Do a histogram of the gain.

%Random distribution for C
C_rand = 0.05*randn(1000) + 0.25;

for i = 1:size(C_rand)
    w = 3*pi;
    
    %Creating C matrix with random numbers
    C_random = [0          0           0 0 0 0 0 0;
            C_rand(i)   -C_rand(i)  0 0 0 0 0 0;
            C_rand(i)    -C_rand(i) 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 L;
            0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0];

    % Making A matrix with random C values
    A = G + 1j*w*C_random;
    
    x = A\F;
    
    Gain_C_random(i) = x(5)/x(1);

end

subplot(2, 2, 4);
histogram(abs(Gain_C_random));
xlabel('Gain');
ylabel('Number of times');

