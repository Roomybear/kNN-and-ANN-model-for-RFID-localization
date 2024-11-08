function[exportData] = test_KNN(secondData,DAQ_Data)

%dataName = 'config2_tx32_rx65';
%dataName2 = 'config2_tx32_rx65DAQ';

%%rawdata = SpeedwayImport(dataName);
%%rawdata2 = DAQ_Import(dataName2);

%total_seconds = rawdata(end,4); %% sets total_seconds equal to last second of data in data set
%% creates several cells the size of total_seconds
%secondData = {1, total_seconds}; 
%timed_error = {1, total_seconds};
%second_Error = {4,total_seconds};
%y0_Coords = {1:total_seconds};
%targetY_Coords = {1:total_seconds}; 

rawdata = secondData;
rawdata2 = DAQ_Data;

n =4; %Number of RF Readers
m = 28; %Number of RFID reference tags
u = 4; %Number of RFID tracking tags
g = 4; %Number of Nearest Neighbors
total_Tags = 32; %Number of Total tags being read

reference_tag_nums = [1:16 21:32];
tracking__tag_nums = [17:20];

tag_cell = {1, 32};
ant_cell = {1, 32};
avg_RSSI = zeros(n,32);



for k = 1:32 %create a cell for each reference tag + tracking tag
         I = find(rawdata(:,3) == k); %finds the index for tag 'k'
         tag_k = rawdata(I,:); % This is the data for the current tag 
         tag_cell{k} = tag_k;  % Place the current tag data into the array
    
     for p = 1:n  %Iterate through each RFID antenna
             I = find(tag_k(:,1) == p); %search through tag k for reads by antenna p and find the index
             ant_p = tag_k(I,:); %for tag k read by antenna p, make an array
             ant_cell{p,k} = ant_p; %for tag k read by antenna p, store array in a cell
             avg_RSSI(p,k) = mean(ant_p(:,2));
     end
end


%%%
% tag_cell{k} is the data for reads of tag k
% ex. tag_cell{7} is the data for tag 7, read by any antenna
%%%
% ant_cell{m ,n} is the data for reads of tag n by antenna m
% ex. ant_cell{3 ,7} is the data for tag 7, read by antenna 3
%%%
% avg_RSSI(p,k) is the RSSI of tag k, read by antenna p 
% ex. avg_RSSI(3 ,7) is the RSSI of tag 7, read by antenna 3

E = zeros(u, length(reference_tag_nums)); %u rows for each tracking tag. 
tempE = []; %Create a temporary storage of E. Since there may be NaN values we will do it this way

for q = 1:length(tracking__tag_nums) %number of tracking tags
    for k = 1:length(reference_tag_nums) %number of reference tags

        for p = 1:n %Number of RFID antennas

            theta_i = avg_RSSI( p,reference_tag_nums(k) ); %Take the RSSI of current reference tag
            S_i     = avg_RSSI( p,tracking__tag_nums(q) ); %Take the RSSI of the current tracking tag
            
            if isnan(theta_i) || isnan(S_i) % if either is not a number, move to next iteration of p (next RFID antenna)
                thisE = 0;
            else
                thisE = (theta_i - S_i); 
            end
            tempE = [tempE, thisE]; %Add current 
            
            
        end

        E(q, reference_tag_nums(k)) = sqrt( sum( tempE.^2 ));
        tempE = [];
        

    end

end

E(:, 17: 18: 19: 20) = [];

%% Coordinates (x,y) of reference tag placement
rtag_coordsX = [0,46,91,140,0,46,91,140,0,46,91,140,0,46,91,140,0,50,0,100,100,150,50,0,50,150,150,100]; %% x coordinates of reference tags
rtag_coordsY = [0,0,0,0,35,35,35,35,78,78,78,78,110,110,110,110,10,9,8.5,9.5,8,9,9,9,9,9,9.5,8]; %% y coordinates of reference tags (is there an easier way to plot coords?)
rtag_coordsZ = [150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,150,0,0,50,50,0,0,100,100,50,100,50,100]; %% z coordinates of reference tags


%% Adds X and Y coordinates of reference tags too E for each antenna read

E(5,:) = rtag_coordsX;
E(6,:) = rtag_coordsY;
E(7,:) = rtag_coordsZ;
%% creates a cell array to add coordinates to tags
coordinate_double = {1:length(n)};

%% creates multiple cells with coordinates and sorts them according to tag E distance
for k = 1:n
   
    coordinate_double{k,:} = [(E(k,:)) ; (E(5,:)) ; (E(6,:)) ; (E(7,:))];  %% Goes through each row of distances and adds coordnate values to each %% if using 3-D add (E(7,:)
    tempC = coordinate_double{k,:}(1,:);        %% creates a temp storage for coordinate_double values
    tempC(tempC == 0) = inf;    %% finds values that are equal to 0 and changes them to inf so they will not be calculated
    coordinate_double{k,:}(1,:) = tempC; %% new coordinate double contains no 0 values
    
    coordinate_double{k,:} = sortrows(coordinate_double{k,:}.',1).'; %% Sorts with respective of distance rows.

        
    for p = 5
        
        for q = 5:length(reference_tag_nums)
            coordinate_double{k,:}(:,5) = []; % Deletes everything but 4 smallest E values with respective coordinates
        end
    end
    tempC = [];
end


W = zeros(g:g); % Zeros array for W

thisW = zeros(n:g); % Zeros Array for calculated numerator

for v = 1:g
    
    for t = 1:g
        current_array = coordinate_double{t,:};
        for q = 1:g
        weight_calc = current_array(1,q);
        
        thisW(t,q) = (1/(weight_calc).^2); %% Numerator calculation for weight       
      
        end 
    
     W(v,t) = (thisW(v,t))/(sum(thisW(v,:))); %%weight calculation 
     
    end 
    coordinate_double{v}(1,:) = W(v,:);
end


targetX = (zeros(1,4))'; % storage for target x coordinate calculation
targetY = (zeros(1,4))'; % storage for target y coordinate calculation
targetZ = (zeros(1,4))'; % storage for target z coordinate calculation

for q = 1:g
    est_coord_array = coordinate_double{q,:}; % Final estimated coordinate array storage
    targetX(q,1) = sum((est_coord_array(1,:) .* est_coord_array(2,:))); % x coordinate calculation
    targetY(q,1) = sum((est_coord_array(1,:) .* est_coord_array(3,:))); % y coordinate calculation
    targetZ(q,1) = sum((est_coord_array(1,:) .* est_coord_array(4,:))); % z coordinate calculation
end


% Real target x coordinate and y coordinate
x0 = [64,64,51,46]';
y0 = [23,35,30,30]';
z0 = [85,85,84,77]';


y0 = y0 + rawdata2; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Creates zero array for estimation error calculation
error = (zeros(4,1));

errorX = (zeros(4,1));
errorY = (zeros(4,1));
errorZ = (zeros(4,1));

for w = 1:4

error(w,1) = sqrt(((targetX(w,1) - x0(w,1))^2 + (targetY(w,1) - y0(w,1))^2  + (targetZ(w,1) - z0(w,1))^2));

errorX(w,1) = targetX(w,1) - x0(w,1); %all of these give the coordinate distance away from actual
errorY(w,1) = targetY(w,1) - y0(w,1);
errorZ(w,1) = targetZ(w,1) - z0(w,1);

%% if using 3 dimensions, use error(w,1) = sqrt(((targetX(w,1) - x0(w,1))^2 + (targetY(w,1) - y0(w,1))^2  + (targetZ(w,1) - z0(w,1))^2));
%% if using 2 dimensions, use error(w,1) = sqrt(((targetX(w,1) - x0(w,1))^2 + (targetY(w,1) - y0(w,1))^2));
%% start making an array for values to input into neural network

exportData = (zeros(1,14));
%%include targetX[0], targetY[0], targetZ[0], est_coord_array [0:], x0[0], y0[0],
%%z0[0], errorX[0], errorY[0] errorZ[0]
%exportData(1,:) = (targetX(1), targetY(1), targetZ(1), est_coord_array(1,:), x0(1), y0(1), z0(1), errorX(1), errorY(1), errorZ(1));
exportData(1) = targetX(4);
exportData(2) = targetY(4);

exportData(3) = targetZ(4);
exportData(4:7) = est_coord_array(4,:);
exportData(8) = x0(4);
exportData(9) = y0(4);
exportData(10) = z0(4);
exportData(11) = errorX(4);
exportData(12) = errorY(4);
exportData(13) = errorZ(4);
exportData(14) = error(4);


end




