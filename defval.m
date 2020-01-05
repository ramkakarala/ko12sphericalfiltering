function defval(name,value)
% DEFVAL(name,value)
%
% Assigns a default value to the named variable
%
% INPUT:
% 
% name    A string, enclosed in single quotes, with a variable name
% value   The value, whatever it is, that you want the variable to have 
%
% OUTPUT:
%
%      None. The variables appear as if by magic into your workspace or
%      will be available inside your function.
%
% NOTE: 
%
% This won't work for an unassigned structure variable.
%
% Last modified by fjsimons-at-alum.mit.edu, 2/16/2009

if ~ischar(name),
  error(sprintf(['The first argument of DEFVAL',...
		'has to be a string with a variable name']));
end

% Always do it
si=1;
% If it exists...
if evalin('caller',[ 'exist(''' name ''')']);
  % ... and it's empty, do it; but don't do it if it's non empty
  si=evalin('caller',[ 'isempty(' name ')']);
end
% Do it or not
if si,
  assignin('caller',name,value);
  na=dbstack;
end

  
