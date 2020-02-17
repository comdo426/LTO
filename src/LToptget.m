function o = LToptget(options,name,default,flag)
%LTOPTGET - gets the LT options
%
%  Syntax:
%     o = LTOPTGET(options,name,default,flag)
%
%  Example:
%		isFeasible = LToptget(Option.LTO, 'FeasibilitySolver');
%		isOptimize = LToptget(Option.LTO, 'Optimizer');
%		isMesh = LToptget(Option.LTO, 'MeshRefinement');
%
%	See also: LTOPTSET
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 17-Feb-2020

if (nargin == 4) && isequal(char(flag),'fast')
   o = getknownfield(options,name,default);
   return
end

if nargin < 2
  error(message('MATLAB:LTodeget:NotEnoughInputs'));
end
if nargin < 3
  default = [];
end
if isstring(name) && isscalar(name)
  name = char(name);
end
if ~isempty(options) && ~isa(options,'struct')
  error(message('MATLAB:LTodeget:Arg1NotLTODESETstruct'));
end

if isempty(options)
  o = default;
  return;
end

Names = [
    'EngineType       '
    'FeasibilitySolver'
    'Optimizer        '
    'InitialStep      '
    'MeshRefinement   '
    'NaturalDynamics  '
    'Ephemeris        '
     ];

names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error(message('MATLAB:LTodeget:InvalidPropName', name));
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    matches = deblank(Names(j(1),:));
    for k = j(2:length(j))'
      matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
    end
    error(message('MATLAB:LTodeget:AmbiguousPropName',name,matches));
  end
end

if any(strcmp(fieldnames(options),deblank(Names(j,:))))
  o = options.(deblank(Names(j,:)));
  if isempty(o)
    o = default;
  end
else
  o = default;
end

% --------------------------------------------------------------------------
function v = getknownfield(s, f, d)
%GETKNOWNFIELD  Get field f from struct s, or else yield default d.

if isfield(s,f)   % s could be empty.
  v = s.(f);
  if isempty(v)
    v = d;
  end
else
  v = d;
end

