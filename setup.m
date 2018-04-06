
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

if (~exist('xgeqp3_m','file') || ~exist('xormqr_m','file') ...
        || (exist('recompile','var') && recompile))
    if (strcmp(systemType,'default'))

        % this seems to work on most systems
        
        fortheaderdir = fullfile(matlabroot,'extern','examples','refbook');
        fortfile = fullfile(matlabroot,'extern','examples','refbook','fort.c');
        
        if (~exist('xgeqp3_m','file') || recompile)
            mex('-output','bin/xgeqp3_m','-v','-largeArrayDims',['-I' fortheaderdir],'src/xgeqp3_m.c',fortfile,'-lmwlapack')
        end
        if (~exist('xormqr_m','file') || recompile)
            mex('-output','bin/xormqr_m','-v','-largeArrayDims',['-I' fortheaderdir],'src/xormqr_m.c',fortfile,'-lmwlapack')
        end
    end
end