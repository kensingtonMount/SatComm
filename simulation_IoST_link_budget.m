Pt = 10;
Ptdb = 10*log10(Pt);
fc = [60 120 180 300 600 1000]* 10^9 ;
lambda = (3*10^8)./fc;
B = 10^7;
T0 = 1500;
Pn = 1.381*10^(-23)*B*T0;
Pndb = 10*log10(Pn);
Prdb = 10 + Pndb;
d = 10:20:700;
Gdb = zeros(1,length(d));
G = zeros(1,length(d));
Bwrad = zeros(1,length(d));
Bwdegree= zeros(1,length(d));
Length= zeros(1,length(d));

for n = 1:6
for k = 1:1:length(d)
Gdb(n,k) = -1*(Ptdb - Prdb + 20*log10(lambda(n)/(4*pi*d(k)*1000)))/2;
G(n,k)= 10^(Gdb(n,k)/10);
Bwrad(n,k)=sqrt((4*pi)/G(n,k));
Bwdegree(n,k)= (Bwrad(n,k)*180)/pi;
Length(n,k)=2*sind(Bwdegree(n,k)/2)*d(k);
end
end

% figure(1);
% for l = 1:6
% hold on;
% plot(d,Gdb(l,:));
% end
% % title('Fc = 0.3 Thz,Bandwidth = 1 MHz,Temperature = 173K,SNR required=10dB');
% xlabel('Transmission Distance [km]');
% ylabel('Path Gain [dB]');
% figure(2);
% for l = 1:6
% hold on;
% plot(d,Bwdegree(l,:));
% end
% % title('Fc = 0.3 Thz,Bandwidth = 1 MHz,Temperature = 173K,SNR required=10dB');
% xlabel('Transmission Distance [km]');
% ylabel('Beamwidth in degrees');



% figure(3);
% plot(d,Length);
% title('Fc = 0.3 Thz,Bandwidth = 1 MHz,Temperature = 173K,SNR required=10dB');
% xlabel('Distance in kilometers');
% ylabel('Arc length in kilometers');

%%GSL to calculate data rate
Pt = 50;
Ptdb = 10*log10(Pt);
fc = [3 12 30 60 120 180 300]* 10^9 ;
lambda = (3*10^8)./fc;
filename = 'bandwidth.xlsx';
A = xlsread(filename,1);
B = A(1:7,7:9);
B_wv = A(10:16,7:9);
%B = 10^10;
T0 = 300; %25 celcius
T_c = 25;
P = 101300.0; %pressure

G_ground = 50;
% Prdb = 10 + Pndb;
d = [500 700 900];
% Gdb = zeros(1,length(d));
% G = zeros(1,length(d));
G_sat = [7.7305 19.7717 27.7305 33.7511 39.7717 43.2935 47.7305];%[15 18 21 Gdb(1:4,25)'];
Bwrad = zeros(1,length(d));
Bwdegree= zeros(1,length(d));
Length= zeros(1,length(d));

for n = 1:length(fc)
    for k = 1:1:length(d)
        Pn = 1.381*10^(-23)*B(n,k)*T0;
        Pndb = 10*log10(Pn);
        %Gdb(n,k) = -1*(Ptdb - Prdb + 20*log10(lambda(n)/(4*pi*d(k)*1000)))/2;
        SNR_0(n,k) = Ptdb + G_ground + G_sat(n) - Pndb + 20*log10(lambda(n)/(4*pi*d(k)*1000))-gaspl(d(k),fc(n), T_c,P,0); %in dry air
        
        Pn_wv = 1.381*10^(-23)*B_wv(n,k)*T0;
        Pndb_wv = 10*log10(Pn_wv);
        SNR_wv(n,k) = Ptdb + G_ground + G_sat(n) - Pndb_wv + 20*log10(lambda(n)/(4*pi*d(k)*1000))-gaspl(d(k),fc(n), T_c,P,7.5); %in water vapor 7.5 g/m^3
        %SNR_wv1(n,k) = Ptdb + G_ground + G_sat(n) - Pndb + 20*log10(lambda(n)/(4*pi*d(k)*1000))-gaspl(d(k),fc(n), T_c,P,12); %in water vapor 12 g/m^3
        % G(n,k)= 10^(Gdb(n,k)/10);
        % Bwrad(n,k)=sqrt((4*pi)/G(n,k));
        % Bwdegree(n,k)= (Bwrad(n,k)*180)/pi;
        % Length(n,k)=2*sind(Bwdegree(n,k)/2)*d(k);
        capacity_dry(n,k) = B(n,k)*log2(1+SNR_0(n,k));
        capacity_wv(n,k) = B_wv(n,k)*log2(1+SNR_wv(n,k));
        %capacity_wv1 = B*log2(1+SNR_wv1);
    end
end

figure;
for m = 1:3
plot(fc/1e9,capacity_dry(:,m)/1e9,'-o');
hold on;
plot(fc/1e9,capacity_wv(:,m)/1e9,'-^');
%plot(fc/1e9,capacity_wv1/1e9,'--d');

end
xlabel('Frequency [GHz]');
ylabel('Data Rate [Gbps]');




