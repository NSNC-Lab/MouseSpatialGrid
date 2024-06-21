function sim_job7

SimIDs=[7]
addpath C:\Users\sanke\Desktop\Github\dynasim
addpath C
addpath \Users\sanke\Desktop\Github\dynasim\functions;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\Coverage;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\Coverage\Coverage;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\Coverage\Coverage\lib;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\Coverage\Coverage\lib\absolutepath;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\Coverage\Coverage\src;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\DataHash_20160618;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\classes;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies\catstruct;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies\kakearney-boundedline-pkg-32f2a1f;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies\kakearney-boundedline-pkg-32f2a1f\boundedline;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies\kakearney-boundedline-pkg-32f2a1f\catuneven;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies\kakearney-boundedline-pkg-32f2a1f\singlepatch;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\dependencies\subplot_grid;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\functions;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\junk;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\library_plots;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD\templates;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MDD_SM;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\MapN;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\Octave;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\altmany-export_fig-2763b78;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\compareFigFiles;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\distinguishable_colors;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\ear;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\github-sync-matlab;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\isfunction;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\m2html;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\m2html\templates;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\m2html\templates\3frames;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\m2html\templates\blue;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\m2html\templates\brain;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\m2html\templates\frame;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\mathworks_exchange;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\mathworks_exchange\uigetvar_putvar;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\mtit;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\subplot_grid;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\dependencies\tight_subplot;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\internal;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\internal\unit-test;C
addpath \Users\sanke\Desktop\Github\dynasim\functions\internal\xp_libraries;
addpath C:\Users\sanke\Desktop\Github\dynasim\models
if strcmp(reportUI,'matlab')
  c=parcluster;
  saveAsProfile(c, sprintf('local_%s', getenv('JOB_ID'))); 
  parpool(c,4); 
else
  disp('   For GNU Octave users: Do not expect any speed up by using DynaSims "parfor_flag". In GNU Octave, parfor loops currently default to regular for loops.');
end
parfor s=1:length(SimIDs)
	SimID=SimIDs(s);
	try
		studyinfoFile = load(fullfile('$HOME\batchdirs\4-channel-PV-inputs',sprintf('studyinfo_%g.mat',SimID)),'studyinfo');
		studyinfo = studyinfoFile.studyinfo;
		[valid,message]=dsCheckHostPaths(studyinfo);
		if ~valid
		  lasterr(message);
		  dsUpdateStudy(studyinfo.study_dir,'sim_id',tSimID,'status','failed'); 
				end
		siminfo = studyinfo.simulations(SimID);
		options = rmfield(siminfo.simulator_options,{'modifications','studyinfo','analysis_functions','plot_functions','sim_id'});
		keyvals = dsOptions2Keyval(options);
		fprintf('-----------------------------------------------------\n');
		fprintf('Processing simulation %g (%g of %g in this job)...\n',SimID,s,length(SimIDs));
		fprintf('-----------------------------------------------------\n');
		data = dsSimulate(studyinfo.base_model,'modifications',siminfo.modifications,'studyinfo',studyinfo,'sim_id',SimID,keyvals{:});
		nFn = length(siminfo.result_functions);
		for iFn = 1:nFn
			try
				fprintf('Result fn: %s (%i/%i)\n', func2str(siminfo.result_functions{iFn}), iFn, nFn)
				dsAnalyze(data, siminfo.result_functions{iFn}, 'result_file',siminfo.result_files{iFn}, 'save_results_flag',1, siminfo.result_options{iFn}{:}, 'in_sim_flag',1);
			catch fnErr
				displayError(fnErr);
			end
		end
	catch err
		displayError(err);
	end
end
if strcmp(reportUI,'matlab')
  delete(gcp)
end
exit
