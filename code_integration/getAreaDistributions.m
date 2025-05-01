function areaStruct = getAreaDistributions(baseDir, cellType, alphaRad, nbins)
% getAreaDistributions
%   baseDir   – top folder that holds all the JBB0.05_JOOx_JBO0.10 subfolders
%   cellType  – which integer in Types_*.dat to use (0 = blue, 1 = orange)
%   alphaRad  – alpha radius for alphaShape (e.g. 1.5)
%   nbins     – #bins for optional histogram (set [] to skip hist)
%
% RETURNS: areaStruct, a struct array with fields
%          .folder   – subfolder name
%          .areas    – vector of region areas
%          .hist     – probability histogram (length = nbins) if nbins>0
%
%  Saves all results to AreaDistributions.mat

if nargin<2, cellType = 0;  end          % default: blue
if nargin<3, alphaRad  = 1.5; end
if nargin<4, nbins     = 20;  end        % 20‑bin histogram by default

subfolders = dir(fullfile(baseDir,'JBB*'));      % all param folders
subfolders = subfolders([subfolders.isdir]);

areaStruct = struct('folder',{},'areas',{},'hist',{});

for k = 1:numel(subfolders)
    foldName = subfolders(k).name;
    dataDir  = fullfile(baseDir,foldName,'ParamSweep_1_Output');
    posFile  = fullfile(dataDir,'Pos_0500000.dat');
    typFile  = fullfile(dataDir,'Types_0500000.dat');
    
    if ~isfile(posFile) || ~isfile(typFile)
        warning('Missing files in %s – skipped.',foldName);
        continue
    end
    
    % --- load complex positions
    txt = fileread(posFile);
    toks = strsplit(txt,',');
    posC = str2double(toks).';          % complex column
    X = real(posC);  Y = imag(posC);
    
    % --- load types and keep desired type
    types = load(typFile);
    keep  = (types == cellType);
    if ~any(keep)
        warning('No cells of requested type in %s',foldName);
        continue
    end
    X = X(keep);   Y = Y(keep);
    
    % --- alphaShape & connected components
    shp = alphaShape(X,Y,alphaRad);
    nReg = numRegions(shp);
    A = zeros(nReg,1);
    for j = 1:nReg
        A(j) = area(shp,j);
    end
    
    % --- build histogram if requested
    if ~isempty(nbins) && nbins>0
        edges = linspace(0,max(A)+eps,nbins+1);
        counts = histcounts(A,edges,'Normalization','probability');
    else
        counts = [];
    end
    
    % --- store
    areaStruct(end+1).folder = foldName;
    areaStruct(end).areas    = A;
    areaStruct(end).hist     = counts;
end
fprintf('Total folders processed: %d\n', numel(subfolders));
fprintf('Area distributions collected: %d\n', numel(areaStruct));

save AreaDistributions.mat areaStruct
fprintf('Saved %d area‑distributions to AreaDistributions.mat\n',numel(areaStruct));
end
