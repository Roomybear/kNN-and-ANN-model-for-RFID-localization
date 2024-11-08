%% Calls data file from excel and seperates into three types
dataName = 'config1_tx28_rx65';
[num, txt, raw] = xlsread('config1_tx28_rx65');

%%deletes lines of data if needed
%%raw(:,7) = [];
%%raw(:,6) = [};
txt(:,5) = [];

%% Array is created that truncates unneaded information
goodArray = erase(txt, '{"antennaPort":');
goodArray = erase(goodArray, 'epc:"');
goodArray = erase(goodArray, 'firstSeenTimestamp:"');
goodArray = erase(goodArray, 'peakRssi:');
goodArray = erase(goodArray, '"');
goodArray = erase(goodArray, 'Z');
goodArray = erase(goodArray, 'ntennaPort:');
goodArray(end,:) = [];

%% Gives titles to table
%%goodArray(1,1:4) = {'Antenna', 'EPC', 'TimeStamp', 'Peak RSSI'};

%% creates start array
b = goodArray(:,2); 

%% Creates a usable double file for EPC's
truncEpc = b;
truncEpc = string(truncEpc);
truncEpc = str2double(truncEpc);


%% creates final data set of same type and puts into one table
d = (goodArray(:,4));
d = string(d);
d = strrep(d,' ','');
d = str2double(d);
%%d = str2num(d);

e = (goodArray(:,1));
e = string(e);
e = str2double(e);

time = goodArray(:,3);
time = string(time);
time = eraseBetween(time,1,11);
time = strrep(time,':','');
time = str2double(time);

finalSet = [e,d,truncEpc,time];
%%finalSet = sortrows(finalSet);


%% n = number of reference tags
n = 16;

%% m = number of reader antennas
m = 1;

%% targetStrength = signal strength vector for ith target tag
%%targetStrength = [];
%% referenceStrength = reference tag strength by reader
%%referenceStrength = [];

%% eDistance = Euclidean distance from target tag to reference tag
%%inside = symsum((targetStrength - referenceStrength), m, 1, 4);
%%eDistance = sqrt((inside)^2);

%% eOrder = sorted euclidean distance in ascending order
%%eOrder = sort(eDistance);

%% weight = weighting factor for each selected reference tag
%%denom = symsum((1/(eOrder)^2),k,1,16);
%% weight = (1/(eOrder)^2)/denom;

%% coordinates[] = estimated coordinates of target tag
%%estcoord = (x:y);
%%coordinates = symsum(weight,k,1,16)*estcoord;

%% p = transmission power


%% lHat = approximate distances from reader to references
%% lHat = []
%% l = real distance between target and reader
%% cooooooords???

%% epsolon position of ith target
%% epsolon = symsum(((lHat-l)/lHat)^2), m, 0, 4)



