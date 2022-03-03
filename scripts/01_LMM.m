function PopAsym_LMM(wdir, adir, sdir, odir, Y, Qdec, subj, surf, mask, cohort, f) 

addpath(genpath('/cluster/projects/p23/tools/mri/freesurfer/current/matlab'))


%inputs
disp(wdir)
disp(adir)
disp(sdir)
disp(odir)
disp(Y)
disp(Qdec)
disp(subj)
disp(surf)
disp(mask)
disp(cohort)


if ~exist(odir); disp(['making outdir: ' odir]); mkdir(odir); end     
cd(odir);


%load data
disp('loading concatenated data')
[Y, mri] = fs_read_Y(Y);
disp('loading sphere')
lhsphere = fs_read_surf([sdir '/' subj '/surf/' surf]);
disp('loading label')
lhcortex = fs_read_label([sdir '/' subj '/label/' mask]);


% sort data
disp('loading Qdec')
Qdec = fReadQdec([adir '/' Qdec]);
fieldNames = Qdec(1, :);
if strcmp('IXI_55minus',cohort)
    disp('covas are age age*hemi sex scanner')
    nbetacols=7
else
    disp('covars are age age*hemi sex')
    nbetacols=5
end
disp(Qdec(1:5,1:nbetacols))


disp('prepping Qdec')
sID = Qdec(2:end,1);   % (grabs the subjects' IDs)
Qdec = rmQdecCol(Qdec,1);  % (removes the subjects'ID column)
M = Qdec2num(Qdec);  % (converts to a numeric matrix)
[M,Y,ni] = sortData(M,1,Y,sID);  % (sorts the data)
%csvwrite('sID.csv',sID);
%csvwrite('M.csv',M);



% matrix and contrasts 
% ....................%
X = [ones(length(M),1) M]; %add intercept


%load([analysisdir filesep 'Xmat.mat'])
% X:    B1 --> Intercept
%       B2 --> Age
%       B3 --> Hemi
%       B4 --> Age*Hemi
%       B5 --> Sex
%       B6 --> Scanner1 (LCBC & IXI)
%       B7 --> Scanner2 (LCBC & IXI)

C = [];
C.age = [[0 1] zeros(1,nbetacols-2)]; % age
C.hemi = [[0 0 1] zeros(1,nbetacols-3)]; %hemi
C.ageHemi = [[0 0 0 1] zeros(1,nbetacols-4)]; % ageHemi
C.sex = [zeros(1,4) 1 zeros(1,nbetacols-5)]; % sex
disp(C)


Cnames = fieldnames(C);
disp(Cnames)
%csvwrite([odir filesep 'Xmat.csv'],X);
%save([odir filesep 'Xmat.mat'],'X');


%Starting values at each location for the linear mixed-effects iterative  estimation.
disp('starting vals')
X(1:5,1:end)
Y(1:5,1:5)
[lhTh0, lhRe] = lme_mass_fit_EMinit(X, [1], Y,  ni, lhcortex, 5);

% covariance estimates are segmented into homogeneous regions. 
disp('cov estimates')
[lhRgs, lhRgMeans] = lme_mass_RgGrow(lhsphere, lhRe, lhTh0, lhcortex, 2,95); 


% Region-wise linear mixed-effects estimation.
lhstats = lme_mass_fit_Rgw(X, [1], Y, ni, lhTh0, lhRgs, lhsphere, [odir filesep 'stats']);

% if loading in precomputed
%lhstats = load([odir filesep 'stats.mat']);
%lhstats = lhstats.stats;

%h0 testing contrasts
for ii = 1:numel(Cnames)
    i = Cnames{ii};
    disp('running contrasts')
    disp(i)

    CM.C = C.(i)
    F_lhstats = lme_mass_F(lhstats, CM);
    
    % write stats
    fs_write_fstats(F_lhstats, mri,[odir filesep i '.sig.mgh'], 'sig');
    fs_write_fstats(F_lhstats, mri,[odir filesep i '.F.mgh'], 'fval');
    
    if (i=="hemi")
        disp('saving degrees of freedom')
        %template
        [vol, DD, mr_parms, volsz]=load_mgh([odir filesep 'age.sig.mgh']);
        % write stats
        tmp=F_lhstats.df()';
        vol=tmp(:,2);
        save_mgh(vol, [odir filesep i '.df.mgh'], DD, mr_parms);
    end
    
    % FDR-correction
    do_saveFDR = 1
    if do_saveFDR == 1  
        
        levels={'0.05','0.01','0.005','0.001'};
        loglev={'13','2','23','3'};
        
        for l = 1:4
            
            level=levels{l};
            logl=loglev{l};
            
            disp(['running FDR ' level])
            mri1 = mri;
            mri1.volsz(4) = 1; 
            [dvtx, spval, pth] = lme_mass_FDR2(F_lhstats.pval, F_lhstats.sgn, lhcortex, level,0);

            if (i=="hemi") && (l==1)
                disp(['writing unsigned sigmap'])
                spvallog = -log10(spval);
                fs_write_Y(spvallog, mri1, [odir filesep i '.spvallog_abs.mgh']);
            end
            disp('writing FDR-corrected siglevel')
            pthlog = [i ',' num2str(pth)]
            dlmwrite([odir filesep 'FDR_' i logl '.txt'],pthlog,'')
            
        end
    end
end

do_saveBeta = 1
nBetas = 1:(nbetacols-2); % select Betas to save
if do_saveBeta == 1
    disp('saving Bs')
    Betadir = [odir filesep 'Beta'];
    if ~exist(Betadir,'dir'); mkdir(Betadir); end
    mri1 = mri;
    mri1.volsz(4) = 1; 
    for j = nBetas
        nv=length(lhstats);
        Beta = zeros(1,nv);
        for i=1:nv % loop vertices
           if ~isempty(lhstats(i).Bhat)
              Beta(i) = lhstats(i).Bhat(j);
           end
        end
        disp(['saving beta' num2str(j)])
        fs_write_Y(Beta,mri1,[Betadir filesep 'Beta' num2str(j) '.mgh']);
    end
    save([odir filesep 'matlab.mat']);
end

