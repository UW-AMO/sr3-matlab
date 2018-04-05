
% can change system type (only default is implemented, currently)

systemType = 'default';

% create bin folder if it doesn't exist

if ~exist('./bin', 'dir')
  mkdir('./bin');
end

addpath([pwd,'/src']);
addpath([pwd,'/test']);
addpath([pwd,'/demo']);
addpath([pwd,'/bin']);

% compile mex binaries

if (strcmp(systemType,'default'))

    fortheaderdir = fullfile(matlabroot,'extern','examples','refbook');
    fortfile = fullfile(matlabroot,'extern','examples','refbook','fort.c');

    mex('-output','bin/xgeqp3_m','-v','-largeArrayDims',['-I' fortheaderdir],'src/xgeqp3_m.c',fortfile,'-lmwlapack')
    mex('-output','bin/xormqr_m','-v','-largeArrayDims',['-I' fortheaderdir],'src/xormqr_m.c',fortfile,'-lmwlapack')
end
