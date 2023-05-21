function getStats2Overlay()
    
    %get projectdir (must run externally as func)
    %/Applications/MATLAB_R2018a.app/bin/matlab -nojvm -nodesktop -r "getStats2Overlay()"
    b=mfilename;
    disp(b);
    tmp = strsplit(b,'/');
    projpath = strjoin(tmp(1:length(tmp)-2),filesep)
    
    projpath='/Users/jamesroe/Dropbox/GitHub/PopAsym/';
    outdir=[projpath '/results/heritability'];
    labelname=[projpath '/annotInfo/vtx.csv'];
    addpath('/Applications/freesurfer/matlab')

    
    % read label
    a = fopen(labelname);
    vtx = textscan(a,'%s');
    lab = str2double(vtx{1}) + 1;
    
	
    % get all raw data - open first image as mgh - template
	[vol, M, mr_parms, volsz] = load_mgh([projpath '/annotInfo/lh.cortex.mgh']);
	vol(:) = 0;


	%get solutions 
	cd(outdir);
    maps=dir('map*.csv');

	for i = 1:length(maps)
        [~, name,ext] = fileparts(maps(i).name); 
        filename = [outdir filesep name, '.mgh'];
        %if ~exist(filename)
            disp(name)
            fid = fopen(maps(i).name); 
            A = csvread(maps(i).name);
            vol = A;
            save_mgh(vol, filename, M, mr_parms);
        %end
    end
  
end
