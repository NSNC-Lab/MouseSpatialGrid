


% %Take apart p
% variableNames = fieldnames(p);
% 
% % Loop through each string and create variables
% for zz = 1:length(variableNames)
%     eval(['global ' variableNames{zz} [' = p.' variableNames{zz}]]);
% end

params = load('params.mat','p');
p = params.p;

fields = fieldnames(p);
                
for zz = 1:numel(fields)
    fieldName = fields{zz};
    value = p.(fieldName);
    assignin('base', fieldName, value);
end