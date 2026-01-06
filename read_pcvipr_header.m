function PCVIPR_HEADER = read_pcvipr_header(headerFile)
% READ_PCVIPR_HEADER
% Reads and parses a PCVIPR header text file
%
% OUTPUT:
%   PCVIPR_HEADER.raw      - raw key/value header fields
%   PCVIPR_HEADER.nframes  - number of reconstructed frames
%   PCVIPR_HEADER.timeres  - temporal resolution (ms)
%   PCVIPR_HEADER.res      - spatial resolution (mm), isotropic
%   PCVIPR_HEADER.matrix  - [Nx Ny Nz]
%   PCVIPR_HEADER.VENC    - velocity encoding (cm/s)

if nargin < 1 || isempty(headerFile)
    headerFile = 'pcvipr_header.txt';
end

fid = fopen(headerFile,'r');
assert(fid>0,'Could not open PCVIPR header file');

formatSpec = '%s%s%[^\n\r]';
dataArray  = textscan(fid, formatSpec, ...
    'Delimiter',' ', ...
    'MultipleDelimsAsOne',true, ...
    'ReturnOnError',false);
fclose(fid);

% Convert value column to numeric where possible
vals = cellfun(@str2num, dataArray{2}, 'UniformOutput', false);
keys = dataArray{1};

PCVIPR_HEADER.raw = cell2struct(vals(:), keys(:), 1);

% ==========================================================
% Derived fields
% ==========================================================
PCVIPR_HEADER.nframes = PCVIPR_HEADER.raw.frames;
PCVIPR_HEADER.timeres = PCVIPR_HEADER.raw.timeres;

PCVIPR_HEADER.res = nonzeros(abs([ ...
    PCVIPR_HEADER.raw.ix, ...
    PCVIPR_HEADER.raw.iy, ...
    PCVIPR_HEADER.raw.iz ]));

PCVIPR_HEADER.matrix = [ ...
    PCVIPR_HEADER.raw.matrixx, ...
    PCVIPR_HEADER.raw.matrixy, ...
    PCVIPR_HEADER.raw.matrixz ];

PCVIPR_HEADER.VENC = PCVIPR_HEADER.raw.VENC;

end