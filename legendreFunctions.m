function P = legendreFunctions(theta, max_deg)
% m-function wrapper for parsing the inputs
if (nargin < 2)
	error('not enough input arguments');
elseif (nargin > 2)
	error('too many input arguments');
end

if (~isreal(theta) || ~isscalar(theta))
	error('theta muste be a real scalar');
elseif (~isreal(max_deg) || ~isscalar(max_deg))
	error('max_deg must be a real scalar');
end

P = legendreFunctionsmx(theta, max_deg);
