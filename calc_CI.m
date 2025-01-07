function CI = calc_CI(data, alpha)


n = length(data); % Number of observations in your sample

df = n - 1; % Degrees of freedom

sampleMean = mean(data); % Sample mean

sampleStd = std(data); % Sample standard deviation

% Calculate the t critical value for the 2-tailed test
tCritical = tinv(1 - alpha/2, df);

% Calculate the margin of error
marginError = tCritical * (sampleStd / sqrt(n));

% Calculate the 90% confidence interval
CI = [sampleMean - marginError, sampleMean + marginError];

end