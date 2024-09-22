function csv_file = runJuliaCode(matfile)

% 06-12-2024 - Jio Nocon
% Calls Julia from command line/terminal to run Jake's code
% See README for more details


global exe

if ispc
    % For Windows, I installed an executable version of Julia, which the
    % exe variable points to
    % exe = "C:\Users\noconjio\AppData\Local\Programs\Julia-1.7.3\bin\julia";
    exe_str = exe;
    exe_str = split(exe_str,filesep); exe_str = join(exe_str,'\\'); exe_str = convertCharsToStrings(exe_str{1});
    

    folder = split(cd,filesep);  folder = join(folder,'\\'); folder = convertCharsToStrings(folder{1});

    command = sprintf('"%s" -e "cd(raw\\"%s\\"); ARGS = [\\"%s\\"]; include(\\"discreteEstimator.jl\\")"', exe_str{1}, folder{1}, matfile);
elseif ismac
    % For Macs, this exe variable should be replaced with the launcher for
    % julia, found using "echo $PATH"

    % Note: installed Julia using the command line (see website for
    % installation details)
    % exe = '/Users/jionocon/.juliaup/bin/julialauncher';
    command = sprintf('%s discreteEstimator.jl "%s"',exe,matfile);
end

[status,cmdout] = system(command);

csv_file = strrep(matfile,'distances.mat','MIs.csv');

end