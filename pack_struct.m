function S = pack_struct(varargin)
% PACK_STRUCT  Pack input variables into a struct using their variable names
%
% Example:
%   S = pack_struct(a, b, c)

S = struct();

for k = 1:nargin
    fname = inputname(k);
    if isempty(fname)
        error('All inputs must be named variables.');
    end
    S.(fname) = varargin{k};
end

end