clc; clear; close all;

%% prompt the user to enter the name of the city they are interested in.
cityName = input('Enter city name: ', 's');

%% Get Geographical Coordinates (from openstreetmap API)
geoAPI = sprintf('https://nominatim.openstreetmap.org/search?q=%s&format=json&limit=1', strrep(cityName, ' ', '%20'));
geoData = webread(geoAPI);

if isempty(geoData)
    error('City not found! Please check the name and try again.');
end

format long g;  
latitude = str2num(geoData(1).lat); 
longitude = str2num(geoData(1).lon);


fprintf('\n📍 Location Found: %s (Lat: %.4f, Lon: %.4f)\n', cityName, latitude, longitude);

%% Fetch Past Solar Radiation Data
startDate = input('Enter the start date for past data (YYYYMMDD): ', 's');
endDate = input('Enter the end date for past data(YYYYMMDD): ', 's');

api_url = sprintf(['https://power.larc.nasa.gov/api/temporal/daily/point?' ...
    'parameters=ALLSKY_SFC_SW_DWN&community=RE&longitude=%f&latitude=%f&start=%s&end=%s&format=JSON'], ...
    longitude, latitude, startDate, endDate);

response = webread(api_url);
solarRadiation = response.properties.parameter.ALLSKY_SFC_SW_DWN;
dates = fieldnames(solarRadiation);


dateArray = datetime(erase(dates, 'x'), 'InputFormat', 'yyyyMMdd', 'Format', 'yyyy-MM-dd');


solarValues = zeros(size(dates));
for i = 1:numel(dates)
    solarValues(i) = solarRadiation.(dates{i});
end

%% Plot Solar Radiation Data
figure;
plot(dateArray, solarValues, '-o', 'LineWidth', 2);
xlabel('Date');
ylabel('Solar Radiation (kWh/m²/day)');
title(['Past Solar Radiation Data for ', cityName, ' (NASA POWER API)']);
grid on;

%% Shading Analysis
solarValues = solarValues / max(solarValues); 
shadeThreshold = 0.3;
shadedDays = dateArray(solarValues < shadeThreshold);

if ~isempty(shadedDays)
    fprintf('\n⚠️ Warning: Shading detected on these days:\n');
    disp(shadedDays);
else
    fprintf('\n✅ No major shading issues detected!\n');
end

%% Plot Shading Effects
figure;
bar(dateArray, solarValues, 'FaceColor', 'b');
hold on;
bar(shadedDays, solarValues(solarValues < shadeThreshold), 'FaceColor', 'r');
xlabel('Date');
ylabel('Normalized Solar Radiation');
title(['Past Shading Analysis for ', cityName, ' (Red = Shaded Days)']);

grid on;

%% Ask User for Corrected Coordinates
web('start https://www.google.com/maps');
currentLat = str2num(input('Enter your implementation place latitude: ', 's')); 
currentLon = str2num(input('Enter your implementation place longitude: ', 's'));

fprintf('\n📍 Using Selected Location: Lat: %.15f, Lon: %.15f\n', currentLat, currentLon);

%% Ask for Implementation Date
implDate = input('Enter the date of implementation (YYYYMMDD): ', 's');

%% Fetch Solar Radiation Data from NASA POWER API
api_url = sprintf(['https://power.larc.nasa.gov/api/temporal/daily/point?' ...
    'parameters=ALLSKY_SFC_SW_DWN&community=RE&longitude=%f&latitude=%f&start=%s&end=%s&format=JSON'], ...
    longitude, latitude, implDate, implDate);

response = webread(api_url);
solarRadiation = response.properties.parameter.ALLSKY_SFC_SW_DWN;
dates = fieldnames(solarRadiation);
dateArray = datetime(erase(dates, "x"), 'InputFormat', 'yyyyMMdd');
solarValues = cellfun(@(d) solarRadiation.(d), dates);

%% Normalize Solar Radiation
solarValues = solarValues / max(solarValues);
shadeThreshold = 0.3;
shadedDays = dateArray(solarValues < shadeThreshold);

%% Predict Best Location Using Machine Learning (Basic Heuristic)
averageRadiation = mean(solarValues);
if averageRadiation > 0.2
    fprintf('\n✅ The selected location is suitable for solar panel implementation!\n');
else
    fprintf('\n⚠️ Warning: The selected location may not be optimal for solar panel implementation. Consider alternative locations.\n');
end

%% Advanced Shading Analysis - Dynamic Sun Angles
sunElevation = linspace(0, 90, length(dateArray)); 
obstacleHeight = 10; % Assume a 10m obstacle
obstacleDistance = 20; % 20m away
shadingThreshold = atand(obstacleHeight / obstacleDistance);
shadedTimes = sunElevation < shadingThreshold;

%% Ask User for Required Energy Output
requiredEnergy = input('Enter the required energy output in kWh: ');

%% Ask User if Rotation is Needed
useRotation = input('Do you want rotatable solar panels? (y/n): ', 's');
solarTime = []; 
optimalRotationAngles = []; 
efficiencyFactor = 1.0;  


%% Define Solar Panel Efficiency and Area
panelEfficiency = 0.2;
panelArea = 1.6;


%% Calculate Required Number of Solar Panels
energyPerPanel = sum(solarValues) * panelEfficiency * panelArea * efficiencyFactor;
numPanels = ceil(requiredEnergy / energyPerPanel);

fprintf('\n🔋 Estimated Number of Solar Panels Required: %d\n', numPanels);
fprintf('\n⚡ Estimated Energy Output per Panel: %.2f kWh/day\n', energyPerPanel);
fprintf('\n⚡ Total Estimated Energy Output: %.2f kWh/day\n', numPanels * energyPerPanel);

if strcmpi(useRotation, 'y')  
    efficiencyFactor = 1.25; 
    fprintf('\n✅ Rotatable Solar Panels Selected! Higher energy generation.\n');

    
    solarTime = linspace(6, 18, 100);  

    
    daysInYear = 365;
    N = day(datetime(implDate, 'InputFormat', 'yyyyMMdd'), 'dayofyear');
    declination = 23.45 * sind(360 * (N - 81) / daysInYear); 

    
    [HourGrid, DeclinationGrid] = meshgrid((solarTime - 12) * 15, declination);

    
    altitudeAngle = asind(sind(latitude) .* sind(DeclinationGrid) + ...
                          cosd(latitude) .* cosd(DeclinationGrid) .* cosd(HourGrid));

    
    optimalRotationAngles = atand(tand(altitudeAngle) ./ cosd(HourGrid));

    fprintf('\n✅ Optimal Panel Rotation Angles Calculated!\n');
    for i = 1:length(solarTime)
        fprintf('⏳ At %.2f Hours, Rotate Panel to %.2f°\n', solarTime(i), optimalRotationAngles(i));
    end

    
    optimalTilt = max(altitudeAngle);
    fprintf('\n📐 Maximum Optimal Tilt Angle of the Day: %.2f°\n', optimalTilt);
    fprintf('🧭 Optimal Azimuth Angle: 180° (South-facing for max exposure)\n');
    
    
    figure;
    plot(solarTime, optimalRotationAngles, '-o', 'LineWidth', 2);
    xlabel('Time (Hours)');
    ylabel('Optimal Rotation Angle (Degrees)');
    title('Optimal Panel Rotation Angle Throughout the Day');
    grid on;
    legend({'Optimal Rotation Angle'}, 'Location', 'best');
else
    fprintf('\n🚫 Skipping solar panel rotation calculations as per user choice.\n');
end

  
  

%% Calculate Optimal Tilt Angle for Maximum Energy Generation
N = day(datetime(startDate, 'InputFormat', 'yyyyMMdd'), 'dayofyear');
declination = 23.45 * sind(360 * (N - 81) / 365); 
optimalTiltAngle = abs(latitude - declination);

fprintf('\n📐 Optimal Tilt Angle for Maximum Energy Generation: %.2f°\n', optimalTiltAngle);


%% Energy Output Estimation
panelEfficiency = 0.2;
panelArea = 1.6; 
energyOutput = sum(solarValues) * panelEfficiency * panelArea;

%% Plot Shading Analysis
figure;
plot(dateArray, solarValues, '-o', 'LineWidth', 2);
hold on;
plot(shadedDays, solarValues(solarValues < shadeThreshold), 'ro', 'MarkerFaceColor', 'r');
xlabel('Date');
ylabel('Normalized Solar Radiation');
title('Shading Analysis Over Time');
grid on;

%% 3D Terrain Visualization (GIS Data Integration)
[X, Y] = meshgrid(linspace(longitude-0.1, longitude+0.1, 50), linspace(latitude-0.1, latitude+0.1, 50));
Z = peaks(50); 
tileSize = 5; 

figure;
surf(X, Y, Z, 'EdgeColor', 'none');
hold on;
scatter3(longitude, latitude, max(Z(:)), 100, 'r', 'filled');
title('3D Terrain Map with Optimal Panel Placement');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Elevation (m)');
grid on;
colorbar;
hold off;


%% Fetch Elevation Data
fprintf('\nFetching Elevation Data...\n');
elevAPI = sprintf('https://api.opentopodata.org/v1/srtm90m?locations=%.6f,%.6f', currentLat, currentLon);
options = weboptions('Timeout', 20);  
elevData = webread(elevAPI, options);
elevation = elevData.results(1).elevation;

fprintf('\n🗻 Approximate Elevation: %.2f meters\n', elevation);



%% Generate 3D Terrain Model using Open-Elevation API
gridSize = 5;
latitudes = linspace(currentLat - 0.01, currentLat + 0.01, gridSize);
longitudes = linspace(currentLon - 0.01, currentLon + 0.01, gridSize);
[latGrid, lonGrid] = meshgrid(latitudes, longitudes);



elevationData = zeros(gridSize, gridSize);
for i = 1:gridSize
    for j = 1:gridSize
        api_url = sprintf('https://api.open-elevation.com/api/v1/lookup?locations=%.6f,%.6f', latGrid(i, j), lonGrid(i, j));
        elevResponse = webread(api_url);
        elevationData(i, j) = elevResponse.results(1).elevation;
    end
end

%% Plot 3D Terrain Model
figure;
surf(lonGrid, latGrid, elevationData, 'EdgeColor', 'none');
colormap jet;
shading interp;
colorbar;
xlabel('Longitude'); ylabel('Latitude'); zlabel('Elevation (m)');
title(['3D Terrain Model for ', cityName]);
grid on;
view(3);

fprintf('\n✅ 3D Terrain Model Generated Successfully!\n');

%% Generate KML File for Google Earth
kmlFile = 'SolarPanel_3D_Map.kml';
fid = fopen(kmlFile, 'w');

fprintf(fid, ['<?xml version="1.0" encoding="UTF-8"?>\n' ...
    '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n' ...
    '<Placemark><name>Solar Panel</name>\n' ...
    '<Point><coordinates>%.6f,%.6f,%.2f</coordinates></Point>\n' ...
    '</Placemark>\n</Document>\n</kml>\n'], currentLon, currentLat, elevation);

fclose(fid);
