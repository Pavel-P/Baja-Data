%% Import Data

addpath('LineCurvature')
data = dlmread('C:\Users\Pavel\Desktop\rc_3.csv', ',', 100,0);

interval_start = 3720;
interval_end = 4040;

time = data(interval_start:interval_end,1);
lat = data(interval_start:interval_end,2);
lon = data(interval_start:interval_end,3);
speed = data(interval_start:interval_end,4);

%% Use moving average filter to reduce noise in gps coords

plot(lat,lon, 'ro');
waitforbuttonpress;

% moving average is equivalent to convolving with rectangle
window_length = 11;
avg_lon = conv(lon, ones(1, window_length)/window_length, 'full');
avg_lat = conv(lat, ones(1, window_length)/window_length, 'full');

% truncate intervals for other plotted data to fit filtered signals
avg_lon = avg_lon(window_length:end-window_length+1);
avg_lat = avg_lat(window_length:end-window_length+1);
avg_time = time(floor(window_length/2):floor(end-window_length/2));
plot(avg_lat, avg_lon, 'bo');

%% Calculate speed from filtered coordinates

% Standard lat/lon distance calculation
R = 6371e3;
lat_r = deg2rad(avg_lat);
lon_r = deg2rad(avg_lon);

del_lat = diff(lat_r);
del_lon = diff(lon_r);

a = sin(del_lat / 2).^2 + cos(lat_r(1:end-1)).* cos(lat_r(2:end)).* sin(del_lon/2).^2;
c = 2*atan2(sqrt(a), sqrt(1-a));

distance = R * c;
avg_speed = distance./diff(avg_time/1000);

hold on;
% avg_speed vector is 1 shorter than other avg data due to diff call
plot(avg_time(2:end), avg_speed * 2.23694, 'r-');
plot(avg_time(2:end), speed(floor(window_length/2 + 1):floor(end - window_length/2)),'b-');
hold off;
%% Use low-pass filter on original speeds

plot(time, speed, 'b');
waitforbuttonpress;

% use multiplication with a rectangle in freq domain as LPF
S = fft(speed);
threshold = 100;
rect_window = cat(1, ones(threshold + 1, 1), zeros(length(S) - 2*threshold - 1, 1), ones(threshold, 1));
lpf_speed = real(ifft(rect_window .* S, length(S)));

hold on;
plot(time, lpf_speed);
hold off;
waitforbuttonpress;

del_s = diff(lpf_speed * 0.44704); %converting speed from mph to m/s
del_t = diff(time/1000);
lpf_accel = del_s./del_t;

% accel vector is 1 shorter than other data due to diff call
plot(time(2:end), lpf_accel(1:end));


%% Plot normals and curvature from original and filtered coordinates

N = LineNormals2D([lat lon]);
k = LineCurvature2D([lat lon]);
quiver(lat, lon, k.*N(:,1), k.*N(:,2));

waitforbuttonpress;

avg_N = LineNormals2D([avg_lat avg_lon]);
avg_k = LineCurvature2D([avg_lat avg_lon]);

quiver(avg_lat, avg_lon, avg_k.*avg_N(:,1), avg_k.*avg_N(:,2));

%% compare original/filtered curvature over time

hold on
plot(time, k);
plot(avg_time, avg_k);
hold off

waitforbuttonpress;

window_length_k = 11;
avg_avg_k = conv(avg_k, ones(1, window_length_k)/window_length_k, 'full');
avg_avg_k = avg_avg_k(window_length:end-window_length+1);

% use this truncation each time you average, identical to earlier section
avg_avg_time = avg_time(floor(window_length_k/2):floor(end-window_length_k/2));
avg_avg_lat = avg_lat(floor(window_length_k/2):floor(end-window_length_k/2));
avg_avg_lon = avg_lon(floor(window_length_k/2):floor(end-window_length_k/2));

threshold = 9500;

plot(avg_avg_time, avg_avg_k);
hold on
plot(avg_time, avg_k);
% everything above or below these lines is a turn (for now)
plot(avg_time, threshold*ones(length(avg_time), 1), 'g-');
plot(avg_time, -threshold*ones(length(avg_time), 1), 'g-');
hold off

%% Plot moving average of filtered acceleration vs position and vs time

window_length_accel = 11;
avg_lpf_accel = conv(lpf_accel, ones(1, window_length_accel)/window_length_accel, 'full');
avg_lpf_accel = avg_lpf_accel(window_length_accel:end-window_length_accel+1);

avg_accel_time = time(floor(window_length_accel/2):floor(end-window_length_accel/2));
avg_accel_lat = lat(floor(window_length_accel/2):floor(end-window_length_accel/2));
avg_accel_lon = lon(floor(window_length_accel/2):floor(end-window_length_accel/2));

figure
scatter(avg_accel_lat(2:end),avg_accel_lon(2:end),6,avg_lpf_accel)
colormap(jet);
colorbar;
title('Clean Accel vs Track Location')
xlabel('Latitude')
ylabel('Longitude')
zlabel('Accel')

waitforbuttonpress;

plot(avg_accel_time(2:end), avg_lpf_accel, 'r');
hold on
plot(time(2:end), lpf_accel,'b');
hold off

%% plot track with turns color-coded
threshold = 9500;
l = logical(abs(avg_avg_k) > threshold);
figure(2);
hold on;
plot(avg_avg_lat(l), avg_avg_lon(l), 'bo');
plot(avg_avg_lat(~l), avg_avg_lon(~l), 'ro');
quiver(avg_lat, avg_lon, avg_k.*avg_N(:,1), avg_k.*avg_N(:,2));
hold off;