function [corr,MSE] = calcModelLoss(model,data)

% each row = observation

% flip model matrix if unequal # columns 
if size(model,2) ~= size(data,2)
model = model';
end

for i = 1:size(model,1)
    % correlation coefficient

temp = corrcoef(model(i,:),data);
corr(i,1) = temp(1,2);

% mean sq error
temp = (model(i,:)-data).^2;
MSE(i,1) = mean(temp);

% st deviation of mse
MSE(i,2) = std(temp);
end

end