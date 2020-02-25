function [mc, mcOpt] = getmc(System, State, Spacecraft, Option)
%GETNX - gets the number of constraints for the fsolve/newtonRaphson/fmincon
%
%  Syntax:
%     [mc, mcOpt] = GETNX(State)
%
%  Description:
%     gets the number of constraints for the each iPhase. 
%
%  Outputs:
%     mx - column vector of numbers of constraints for each phase
%		mxOpt - column vector of numbers of constraints for each phase(optimizer)
%
%  See also: SETPROBLEMFSOLVE, SETPROBLEMFMINCON, GETNX
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(State);
mc = nan(nPhase,1);
mcOpt = nan(nPhase,1);

for iPhase = 1:nPhase
	mLnrAddCon = 0; % number of linear additional conditions
	mNlnrAddCon = 0; % number of nonlinear additional conditions
	
	[~, nSegment, ~, ~, ~, mDefect, ~] = getPhaseStateInfo(State{iPhase});
	mThrust = nSegment;
	mCont = 7*(nSegment-1);
	
	if ~isempty(Option.AddCon{iPhase, 1})
		mnAddCon = size(Option.AddCon);
		nAddCon = mnAddCon(1, 2);
		
		for iAddCon = 1:nAddCon
			Con = Option.AddCon{iPhase, iAddCon};
			isNlnr = strcmp(Con{1}, 'nlnr'); % is nonlinear
			c = str2func(strcat('get',Con{3},'ConstraintNo')); % define function
			if isNlnr
				m = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase},	Con);
				mNlnrAddCon = mNlnrAddCon + m;
			else
				m = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase}, ...
					Con);
				mLnrAddCon = mLnrAddCon + m;
			end
		end
	end
	mc(iPhase) = mDefect + mCont + mThrust + mLnrAddCon + mNlnrAddCon; 
	mcOpt(iPhase) = mDefect + mCont + mNlnrAddCon; % linear constraints are separate
end
end