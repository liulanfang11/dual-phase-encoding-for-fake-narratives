% Code used in Vidaurre et al. (2017) PNAS
% from PNAS2017_fullpipeline_llf
% Detailed documentation and further examples can be found in:
% https://github.com/OHBA-analysis/HMM-MAR
% This pipeline must be adapted to your particular configuration of files. 
%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS AND FILE NAMES
clear
addpath(genpath('.\MovieBrainDynamics-master\scripts\hmm_new')) % add toolbox
mydir='.\HMM_results\';
J=9; % the number of networks
repetitions = 10; % to run it multiple times (keeping all the results)

load('.\sub_tc_pool.mat')
N=size(sub_tc,3); % number of sujbects

TRs=300;
for K=2:20
    DirOut = [mydir 'state',num2str(K,'%02d'),'\iterations\'];
    mkdir(DirOut)
    TR = 1;  
    use_stochastic = 0; % set to 1 if you have loads of data

    f = cell(N,1); T = cell(N,1);
    % each .mat file contains the data (ICA components) for a given subject, 
    % in a matrix X of dimension (4800time points by 50 ICA components). 
    % T{j} contains the lengths of each session (in time points)
    for j=1:N
        f{j} = sub_tc(:,:,j);
        T{j} = TRs;
    end

    options = struct();
    options.K = K; % number of states 
    options.order = 0; % no autoregressive components
    options.zeromean = 0; % model the mean
    options.covtype = 'full'; % full covariance matrix
    options.Fs = 1/TR; 
    options.verbose = 1;
    options.standardise = 1;
    options.inittype = 'HMM-MAR';
    options.cyc = 500;
    options.initcyc = 10;
    options.initrep = 3;

    % stochastic options
    if use_stochastic
        options.BIGNbatch = round(N/30);
        options.BIGtol = 1e-7;
        options.BIGcyc = 500;
        options.BIGundertol_tostop = 5;
        options.BIGforgetrate = 0.7;
        options.BIGbase_weights = 0.9;
    end

    % We run the HMM multiple times
    for  r = 1:repetitions

         disp(['RUN ' num2str(r)]);
         [hmm, Gamma, ~, vpath, ~, ~, fe] = hmmmar(f,T,options);
         save([DirOut 'HMMrun_rep_', num2str(r,'%02d'), '.mat'],'Gamma','vpath',...
            'hmm','T','J','K','N','fe');

        % calculate summary measures for this HMM and save those
        % too to the disk:
         mean_em=zeros(J,K); % mean emissions
         for k=1:K
             mean_em(:,k)=getMean(hmm,k);
         end
         prob=hmm.P; % transition probabilities

         do = 1; % do = Flag to replace 0s with NaNs (in all summary measures)
                    % do=0 keeps all misisng values as zeros; do=1 sets them
                    % to Nans; this is important while comparing between groups!!
        [FO,dign,avg_life]=summary_measures(vpath,N,do);
        maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM 
                % is capturing dynamics or grand between-subject 
                % differences (see Wiki)           
         fe_final = hmmfe(f,T,hmm);
        save([DirOut 'Summary_measures_rep_',num2str(r,'%02d') '.mat'],...
         'mean_em','prob','FO','dign','avg_life','maxFO','fe_final');
       
    end
      
        
    cd(DirOut)

end

    
   