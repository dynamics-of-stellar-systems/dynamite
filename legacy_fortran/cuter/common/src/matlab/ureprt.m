function calls=ureprt
%   ( Last modified on 23 Dec 2000 at 17:29:50 )
%
% UREPRT determines statistics on the current run,
%        for unconstrained or bound-constrained problems.
%        
%
%   [time,calls]=ureprt returns user/system time in time
%                     and the number of calls to f, g, H in calls.
%
calls=utools('ureprt');

