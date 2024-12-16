% Define filter parameters
[b, a] = butter(4, 0.15);  % 4th-order Butterworth filter
THR = 6;
% Apply zero-phase filtering
filtered_data = filtfilt(b, a, data{1}.Values.Data);
for i=1:length(filtered_data)
    if(abs(filtered_data(i)) < THR)
        filtered_data(i) = 0;
    end
end
% Use the filtered data in the Parameter Estimator app
plot(filtered_data);
