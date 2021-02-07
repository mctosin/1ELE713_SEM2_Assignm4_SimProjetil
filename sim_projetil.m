%%Carrega e define variáveis
clear all
close all
load Evoo.mat; %Arquivo com os dados adquiridos de empuxo no tempo
tsim = 32; %Tempo de simulação
f = 1000; %pontos adquiridos por segundo (Hz)
tplot = [0:(1/f):tsim]'; %Cria vetor com 32000 elementos onde tplot(n+1)-tplot(n)=1/f

h_londrina = 1; 
g = 9.8; %gravidade'
vsom = 340;
Is = 101; %impulso específico
tq = 4.0; %segundos
mi = 12;  %massa inicial do foguete

%definir o vetor empuxo
E(1:tsim*f) = 0;
E = E';
E(1:tq*f) = Evou(401:4400);

% Definindo o vetor de massa
I = cumsum(E)/f; % vetor impulso
It = sum(E)/f; % impulso total
mp = It/(g*Is); % massa de propelente
m = mi - mp*(I/It); %vetor de massa do foguete

%define o vetor força peso;
P = m * g;

A = (pi/4)*(0.1^2); %área da seção transvesal principal do foguete

%To = %temperatura ambiente (518.4 R = 15 C)
%b = %gradiente da temperatura com a altitude (0.003566R/ft)
%R = %constante dos gases para o ar (53.33089 ft/R)
%rho_0 = 0.002378; %densidade inicial do ar (slug/ft^3) em 744mmHg e 15C

%valores iniciais das 
rho_0 = 0.002378 * 515.4; %densidade inicial do ar (kg/m^3)
hmar = h_londrina;
v = zeros(tsim*f+1,1);
h = zeros(tsim*f+1,1);
Cd = zeros(tsim*f+1,1);
Fd = zeros(tsim*f+1,1);
Fr = zeros(tsim*f+1,1);
Acc = zeros(tsim*f+1,1);
v(1) = 0;
h(1) = 0;
Cd0 = 0.23;

for i = 1 : tsim*f
    % rho(i) =  rho_0*(1-(b/To)*(hmar+h(i)))^((1-b*R)/(b*R));
    rho(i) = rho_0 * (1 - (0.00000688/0.3048) * (hmar+h(i))).^ 4.256;
    if v(i) < 238 % v<238
        Cd(i) = Cd0;
    elseif v(i) < 442
        Cd(i) = 0.32308*log(v(i)/vsom)+0.3452 ; %238<=v<442
    else
        break; %v>=442 modelo não funciona
    end
    Fd(i) = (1/2)*Cd(i)*A*rho(i)*(v(i)^2);
    Fr(i) = E(i) - P(i) - Fd(i);
    Acc(i) = Fr(i) / m(i);
    v(i+1) = v(i) + (1/f) * Acc(i);
    h(i+1) = h(i) + (1/f) * v(i);
end

%plots
subplot(4,1,1)
plot(tplot(1:tsim*f), Acc(1:tsim*f)/g, 'b')
ylabel('Aceleração (g)')
subplot(4,1,2)
plot(tplot(1:tsim*f),   v(1:tsim*f), 'b')
ylabel('Velocidade (m/s)')
subplot(4,1,3)
plot(tplot(1:tsim*f),   h(1:tsim*f), 'b')
ylabel('Altura (m)')
subplot(4,1,4);
plot(tplot(1:tsim*f),   -Fd(1:tsim*f), 'b')
ylabel('Drag (kgm/s?2)')
xlabel('Tempo (s)')


