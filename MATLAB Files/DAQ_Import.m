function[daq_Info] = DAQ_Import(filename)

[num, txt, raw] = xlsread(filename);
%{
for p = 4
    for l = 1:5
        raw(:,p) = [];
    end
end

for p = 1:4
    for l = 1
        raw(l,:) = [];
    end
end
%}

distance = string(raw(:,3)); % Puts distance column into its own array
distance = str2double(distance); % Changes distances to a double
distance(:,1) = -(distance(:,1) * 2.54); %Converts distance from inches to centimeters

time = string(raw(:,1)); %puts time into its own array
time = str2double(time); %makes strings into doubles that can be used in for loop

for k = 1:length(raw) %% loop changes time to numbers instead of seconds (Ex. time 12:01:01 would now be just 1)
    time(k,1) = k;
end

daq_Info = [time, distance]; % outputs final array that contains converted distance and time

end
