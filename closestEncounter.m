function closestEncounter(System, State)

nPhase = length(System);

for iPhase = 1:nPhase
   
   lstar = System{iPhase}.parameter.lstar;
   isCR3BP = strcmp(System{iPhase}.dynamics{1}, 'CR3BP');
   if isCR3BP
      mu = System{iPhase}.parameter.mu;
      primary1Position = [-mu, 0, 0];
      primary2Position = [1-mu, 0, 0];
      [stateMat, ~, ~] = getStateControlMat(State{iPhase});
      [~, ~, nNode, ~, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
      primary1Distance = nan(nNode, 1);
      primary2Distance = nan(nNode, 1);
      for iNode = 1:nNode
         position = stateMat(iNode,1:3);
         primary1Distance(iNode) = lstar*norm(position-primary1Position);
         primary2Distance(iNode) = lstar*norm(position-primary2Position);
      end
      
      [minDistance1 minNode1] = min(primary1Distance);
      [minDistance2 minNode2] = min(primary2Distance);
      minSegment1 = ceil(minNode1/3);
      minSegment2 = ceil(minNode2/3);
      fprintf('Phase no. %d, Closest encounter: %0.5e(Primary 1) at Segment no. %d, %0.5e(Primary2) at Segment no. %d\n', ...
         iPhase, minDistance1, minSegment1, minDistance2, minSegment2);
   else
      fprintf('Other dynamics not supported yet')
   end
   
end

end