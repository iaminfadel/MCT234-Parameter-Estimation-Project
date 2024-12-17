% Filter parameters
[b, a] = butter(4, 0.15);
THR = 6;
% Apply filter
filtered_data = filtfilt(b, a, data{1}.Values.Data);
for i=1:length(filtered_data)
    if(abs(filtered_data(i)) < THR)
        filtered_data(i) = 0;
    end
end
plot(filtered_data);
