dataName = 'config1_tx30_rx65';
dataName2 = 'config1_tx30_rx65DAQ';

rawdata = SpeedwayImport(dataName);
rawdata2 = DAQ_Import(dataName2);

total_seconds = rawdata(end,4);
second_Data = {1, total_seconds};
timed_error = {1, total_seconds};
second_Error = {4,total_seconds};
y0_Coords = {1:total_seconds};
targetY_Coords = {1:total_seconds};

for q = 1:total_seconds
    S = find(rawdata(:,4)== q); %% organizes array by seconds
    sec_q = rawdata(S,:);
    second_Data{q} = sec_q; %% stores single second into an array
    current_timed_error = test_KNN(sec_q,rawdata2(q,2));%%sends 1 second of data to test_KNN function to be calculated
    timed_error{1,q} = current_timed_error; %%each second of calculations for 4 target tags are put into their own array
    y0_Coords = [];
    
end


x = 1:total_seconds;
y = (1:total_seconds);
y1 = (1:total_seconds);
y2 = (1:total_seconds);
y3 = (1:total_seconds);
y4 = (1:total_seconds);

for k = 1:total_seconds
    y(1,k) = rawdata2(k,2);
    y1(1,k) = timed_error{1,k}(1,1); %% array for each coordinate of the first target tag
    y2(1,k) = timed_error{1,k}(2,1); %% array for each coordinate of the second target tag
    y3(1,k) = timed_error{1,k}(3,1); %% array for each coordinate of the third target tag
    y4(1,k) = timed_error{1,k}(4,1); %% array for each coordinate of the fourth target tag
end

figure
plot(x,y,x,y4);
