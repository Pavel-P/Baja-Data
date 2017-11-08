%% Import Data

data = dlmread('C:\Users\Pavel\Desktop\rc_3.csv', ',', 100,0);

i_start = 3720;
i_end = 4040;

time = data(i_start:i_end,1);
lat = data(i_start:i_end,2);
lon = data(i_start:i_end,3);
speed = data(i_start:i_end,4);

%% Use moving average filter to reduce noise in gps coords

plot(lat,lon, 'ro');
waitforbuttonpress;

window = 11;
avg_lon = conv(lon, ones(1, window)/window, 'full');
avg_lat = conv(lat, ones(1, window)/window, 'full');

avg_lon = avg_lon(window:end-window+1);
avg_lat = avg_lat(window:end-window+1);
valid_time = time(floor(window/2):floor(end-window/2));
plot(avg_lat, avg_lon, 'bo');

%% Calculate speed from filtered coordinates

R = 6371e3;
lat_r = deg2rad(avg_lat);
lon_r = deg2rad(avg_lon);

del_lat = diff(lat_r);
del_lon = diff(lon_r);

a = sin(del_lat / 2).^2 + cos(lat_r(1:end-1)).* cos(lat_r(2:end)).* sin(del_lon/2).^2;
c = 2*atan2(sqrt(a), sqrt(1-a));

d = R * c;
avg_speed = d./diff(valid_time/1000);

hold on;
plot(valid_time(2:end), avg_speed * 2.23694, 'r-');
plot(valid_time(2:end), speed(floor(window/2 + 1):floor(end - window/2)),'b-');
hold off;
%% Use FFT zeroing to filter original values

plot(time, speed, 'b');
waitforbuttonpress;

S = fft(speed);
cutoff = 100;
N = length(S);
s = real(ifft(cat(1, S(1:cutoff + 1), zeros(N - 2*cutoff - 1, 1), S(end - cutoff + 1:end)), N));

hold on;
plot(time, s);
hold off;
waitforbuttonpress;

del_s = diff(s * 0.44704);
del_t = diff(time/1000);
accel = del_s./del_t;
plot(time(2:end), accel(1:end));


%% Calculate normals from original and filtered coordinates

N_raw = LineNormals2D([lat lon]);
k_raw = LineCurvature2D([lat lon]);
quiver(lat, lon, k_raw.*N_raw(:,1), k_raw.*N_raw(:,2));

waitforbuttonpress;

N_avg = LineNormals2D([avg_lat avg_lon]);
k_avg = LineCurvature2D([avg_lat avg_lon]);

quiver(avg_lat, avg_lon, k_avg.*N_avg(:,1), k_avg.*N_avg(:,2));

%% compare curvatures over time

hold on
plot(time, k_raw);
plot(valid_time, k_avg);
hold off

waitforbuttonpress;

window = 11;
k_filter = conv(k_avg, ones(1, window)/window, 'full');
k_filter = k_filter(window:end-window+1);

time_filter = valid_time(floor(window/2):floor(end-window/2));
lat_filter = avg_lat(floor(window/2):floor(end-window/2));
lon_filter = avg_lon(floor(window/2):floor(end-window/2));

cutoff = 9500;

plot(time_filter, k_filter);
hold on
plot(valid_time, k_avg);
plot(valid_time, cutoff*ones(length(valid_time), 1), 'g-');
plot(valid_time, -cutoff*ones(length(valid_time), 1), 'g-');
hold off

%% plot curve meeting threshold
cutoff = 9500;
l = logical(abs(k_filter) > cutoff);
figure(2);
hold on;
plot(lat_filter(l), lon_filter(l), 'bo');
plot(lat_filter(~l), lon_filter(~l), 'ro');
quiver(avg_lat, avg_lon, k_avg.*N_avg(:,1), k_avg.*N_avg(:,2));
hold off;