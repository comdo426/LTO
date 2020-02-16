function [mc, mcOpt] = getmc(System, State, Spacecraft, Option)

nPhase = length(State);
mc = nan(nPhase,1);
mcOpt = nan(nPhase,1);

for iPhase = 1:nPhase
	mLnrAddCon = 0;
	mNlnrAddCon = 0;
	
	[~, mThrust, ~, ~, ~, mDefect, ~] = getPhaseStateInfo(State{iPhase});
	
	if ~isempty(Option.AddCon{iPhase, 1})
		mnAddCon = size(Option.AddCon);
		nAddCon = mnAddCon(1, 2);
		
		for iAddCon = 1:nAddCon
			Con = Option.AddCon{iPhase, iAddCon};
			isNlnr = strcmp(Con{1}, 'nlnr'); % is nonlinear
			c = str2func(strcat('get',Con{3},'ConstraintNo'));
			if isNlnr
				m = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase}, ...,
					Con);
				mNlnrAddCon = mNlnrAddCon + m;
			else
				mLnrAddCon = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase}, ...
					Con);
				mLnrAddCon = mLnrAddCon + m;
			end
		end
	end
	mc(iPhase) = mDefect + mThrust + mLnrAddCon + mNlnrAddCon;
	mcOpt(iPhase) = mDefect + mNlnrAddCon;
end
end