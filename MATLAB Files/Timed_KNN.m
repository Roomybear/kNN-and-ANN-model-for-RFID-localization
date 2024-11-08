
dataName = 'config1_tx31_rx75';
dataName2 = 'config1_tx31_rx75DAQ';

rawdata = SpeedwayImport(dataName);
rawdata2 = DAQ_Import(dataName2);

rows = 13;
columns = (1:total_seconds);
total_seconds = rawdata(end,4); %% sets total_seconds equal to last second of data in data set
%% creates several cells the size of total_seconds
second_Data = {1, total_seconds}; 
timed_error = {1, total_seconds};
second_Error = {4,total_seconds};
y0_Coords = {1:total_seconds};
targetY_Coords = {1:total_seconds}; 
transformed_Data = {1:total_seconds}; %% array that will contain Fourier transform of Y coordinates over time
%csvOut(1,:) =1:13;

for q = 1:total_seconds
    S = find(rawdata(:,4)== q); %% organizes array by seconds
    sec_q = rawdata(S,:);
    second_Data{q} = sec_q; %% stores single second into an array
    current_timed_error = test_KNN(sec_q,rawdata2(q,2));%%sends 1 second of data to test_KNN function to be calculated
    timed_error{1,q} = current_timed_error; %%each second of calculations for 4 target tags are put into their own array
    y0_Coords = [];
end

csv_out = zeros(total_seconds,14); %storage for matrix
for i = 1:total_seconds
    for j = 1:14
        csv_out(i,j) = timed_error{i}(j);
    end
end
zeroed_DAQ = (1:total_seconds);
for k = 1:total_seconds
    zeroed_DAQ(1,k) = rawdata2(k,2);
end

daq = fft(zeroed_DAQ, total_seconds);
x = 1:total_seconds;
y = (1:total_seconds);

y1 = (1:total_seconds);
y2 = (1:total_seconds);
y3 = (1:total_seconds);
y4 = (1:total_seconds);

y5 = (1:total_seconds);
y6 = (1:total_seconds);
y7 = (1:total_seconds);
y8 = (1:total_seconds);

writematrix(csv_out, "Tag4_3.csv")
%{
for k = 1:total_seconds
    y(1,k) = zeroed_DAQ(1,k)+29;
    
    y1(1,k) = timed_error{1,k}(1,1); %% array for each coordinate of the first target tag
    y2(1,k) = timed_error{1,k}(2,1); %% array for each coordinate of the second target tag
    y3(1,k) = timed_error{1,k}(3,1); %% array for each coordinate of the third target tag
    y4(1,k) = timed_error{1,k}(4,1); %% array for each coordinate of the fourth target tag
    
end

figure
plot(x,y,x,y4); %% figure compares Y coordinates of tags vs position of actuator
title('Target Coordinates compared to Actuator position')
xlabel('Time(s)');
ylabel('Y coordinates(cm)');

Y = fft(y2,total_seconds);
Y_mag = abs(Y);
%% takes the average of all Y-coordinates too see accuracy
averageY1 = mean(y1);
averageY2 = mean(y2);
averageY3 = mean(y3);
averageY4 = mean(y4);



signal1 = -y2;
signal2 = daq;

Fs = 100;
T = 1/Fs;
L = total_seconds;
t = (0:L-1)*T;
Y = fft(rawdata(:,2));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

plot(singal1, signal2)
title('Single-Side Amplitude Spectrum of Actuator')
xlabel('f(Hz)')
ylabel('|P1(f)|')

%{
fftSignal1 = fft(signal1);
fftSignal1 = fftshift(fftSignal1);
fftSignal2 = fft(signal2);
fftSignal2 = fftshift(fftSignal2);
%%f = Fs/2*linspace(-1,1,Fs);

figure
plot(f,abs(fftSignal1));
title('Magnitude fft of target tag')
xlabel('Frequency(Hz)');
ylabel('Magnitude');

figure
plot(f,abs(fftSignal2));
title('magnitutde fft of DAQ')
xlabel('Frequency(Hz)');
ylabel('Magnitude');

%}




    %}