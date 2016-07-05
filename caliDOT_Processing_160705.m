function caliData
%% CALI CART DOT PROCESSING SCRIPT
% 1. Preliminaries
% 2. Flags
% 3. Load/Prepare Data
% 4. Data Quality Visualizations
% 4. Filter, etc.
% 5. Image

%% 1. Preliminaries
disp('>> DOT PRELIMINARIES...')
% Get Images to process
[files,path] = uigetfile('*.mag','Select runs to reconstruct',...
    'MultiSelect','on');

if ~iscell(files)
    foo = files;
    clear files
    files{1} = foo;
end

% Load anatomy
[MNI,infoA] = MRItoDOT([],'mni152nl_T1_on_333','4dfp');
t4=eye(4);

% Load grid

% Load A

here = pwd;
cd(path);

%% 2. Flags
disp('>> SETTING FLAGS...')

flags.tag='_MNI_Amat_160705';
flags.cap='cali';

flags.bthresh=0.075;        % Noisy channel threshold
flags.HP0_co=0.003;
flags.omega_hp=0.02;        % Highpass cutoff: REST=0.009; TASK=0.02
flags.omega_lp1=1;          % Lowpass filter 1 cutoff
flags.omega_lp2=0.5;
flags.gi='nn123';
flags.lambda1=1e-2;         % Tikhonov regularization (resolution)
flags.lambda2=1e-1;         % Spatially-dependent regularization
flags.gbox=5;               % Smoothing Gaussian range
flags.gsigma=1.2;           % Smoothing Gaussian width

flags.t4=t4;
flags.gsr=1;                % Global (Superficial) Signal Regression?
flags.resample=1;
flags.freqout=1;
flags.rstol=1e-5;

flags.infoA=infoA;

load(flags.a,'iRad')

%% 3. Load/Prepare Data
for f = 1:length(files)
    disp('>> LOADING DATA...')
    fname = files{1}(1:end-4);
    
    % Load the .mag file and read the header
    fid = fopen([fname '.mag'], 'rb');
    [ndet nsrc ncolor] = read_header(fid);
    data=fread(fid, 'double');
    fclose(fid);
    
    % Reshape the data into standard src x det x color x time
    T=numel(data)/(ndet*nsrc*ncolor);
    data=reshape(data,ncolor,nsrc,ndet,T);
    data=data(:,:,:,2:end);
    [Ns,Nd,Nc,Nt]=size(data);
    
    % Grab only nn12345
    data=reshape(data,Ns*Nd,Nc,Nt);
    data=data(info.Rad.nn12345,:,:);
    [Nm,Nc,Nt]=size(data);
    
    
    %% Apply HP0 before proceeding
    disp('>> HP0...')
    data=reshape(data,[],Nt);
    for j=1:(Nm*Nc)
        foo=squeeze(data(j,:));
        foob=mean(foo);
        foo=-log(foo./foob);
        foo=highpass(foo,flags.HP0_co,info.framerate);
        data(j,:)=foob.*exp(-foo);
    end
    data=reshape(data,Nm,Nc,Nt);
    
    %% Filter Data, Find Noisy Measurements, Logmean
    % logmean
    nlrdata=logmean(data);
    
    % good measurements
    for c=1:Nc
        info.goodindex.(['c',num2str(c)])=...
            intersect(find(std(squeeze(nlrdata(:,c,t0:tF)),1,2)...
            <=flags.bthresh),Rad.([flags.gi]));
    end
    
    % hp
    nlrdata=highpass(nlrdata,flags.omega_hp,info.framerate);
    % lp1
    nlrdata=lowpass(nlrdata,flags.omega_lp1,info.framerate);
    % ssr
    hem=mean(squeeze(nlrdata(intersect(info.goodindex,iRad.nn1),1,:)),1);
    nlrdata=regcorr(nlrdata,hem);
    % lp2
    nlrdata=lowpass(nlrdata,flags.omega_lp2,info.framerate);
    
    %% Imaging
    [cortexmua,cortexhb,dim]=DotImage3D_ATE(nlrdata,flags,info,A);
    Map_DOT_to_4dfp(cortexhb,dim,infoA,t4,datafile);
end

cd(here)
end

%% read_header()
function [ndet nsrc ncolor fmt timesecs] = read_header(fid)
disp('>> READING HEADER...')
fmt = fread(fid, 1, 'uint16');
nheader_bytes = fread(fid, 1, 'uint16');
ndet = fread(fid, 1, 'uint16');
nsrc = fread(fid, 1, 'uint16');
ncolor = fread(fid, 1, 'uint16',2); % the extra bit here has a meaning, but I don't remember it
timesecs = fread(fid, 1, 'uint32');
end

%% logmean()
function [data2]=logmean(data)

[data Sin Sout]=datacondition(data,1);
data2=zeros(Sout,'double');

for n=1:Sout(1)
    switch class(data)
        case 'double'
            tmp=squeeze(data(n,:));
        otherwise
            tmp=double(squeeze(data(n,:)));
    end
    tmp2=-log(tmp./mean(tmp));
    data2(n,:)=tmp2;
    clear tmp tmp2;
end
data2=reshape(data2,Sin);
end

function [cortexmua,cortexhb,dim]=DotImage3D(nlrdata,flags,info,A)
% Nvox=numel(dim.Good_Vox);
[Nm,Nc,Nt]=size(nlrdata);
Nvox=size(A,3);

% Invert
disp('>Inverting A')
tic
for c=1:Nc
    disp(['size A = ',num2str(size(...
        single(squeeze(A(c,info.reconindex.(['c',num2str(c)]),:)))))])
    
    iA.(['c',num2str(c)])=single(TikhonovInvert...
        (single(squeeze(A(c,info.reconindex.(['c',num2str(c)]),:))),...
        flags.lambda1,flags.lambda2));
end
toc

disp('>Smoothing iA')
tic
for c=1:Nc
    iA.(['c',num2str(c)])=single(smoothA3D(iA.(['c',num2str(c)]),...
        flags.gbox,flags.gsigma,dim));
end
toc

% Perform imaging
disp('>Performing Imaging')
cortexmua=zeros(Nvox,Nc,Nt,'single');
for c=1:Nc
    cortexmua(:,c,:)=iA.(['c',num2str(c)])*...
        squeeze(nlrdata(info.reconindex.(['c',num2str(c)]),c,:));
end
cortexmua=cortexmua./100; % AD-HOC UNITS FIX.

% Load spectroscopy matrix
disp('>Performing Spectroscopy')
load(['E_',flags.e])
[Cin Cout]=size(E);
iE=inv(E);

% Initialize Outputs
cortexhb=zeros(Nvox,Nc,Nt,'single');
for h=1:Cout
    tseq=zeros(Nvox,Nt);
    for c=1:Cin
        tseq=tseq+squeeze(iE(h,c))*squeeze(cortexmua(:,c,:));
    end
    cortexhb(:,h,:)=tseq;
end
clear tseq
cortexhb=cortexhb.*1000;
end

function Map_DOT_to_4dfp(cortexhb,dim,infoA,t4,filename)
%
% This function maps DOT imaged data (cortexhb)to a 4dfp volume defined 
% by the infoA file.

Hb={'HbO','HbR','HbT'};
Nt=size(cortexhb,3);

for H=1:3
    foob=zeros(dim.nVx*dim.nVy*dim.nVz,Nt);
    
    switch H
        case {1,2}
            for j=1:Nt
                foob(dim.Good_Vox,j)=squeeze(cortexhb(:,H,j));
            end
        case 3
            for j=1:Nt
                foob(dim.Good_Vox,j)=squeeze(cortexhb(:,1,j))+...
                    squeeze(cortexhb(:,2,j));
            end
    end
    foob=reshape(foob,dim.nVx,dim.nVy,dim.nVz,Nt);
    foob = Affine3d(foob,dim,infoA,t4);
    
    infoHb=infoA;
    infoHb.nVt=Nt;
    MATLABto4dfp(foob,infoHb,[filename,'_',Hb{H}]);

end
end