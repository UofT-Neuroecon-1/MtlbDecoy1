function [ X_1decoy, X_2decoy , TargetAndAltX1, TargetAndAltX2 ] = AddDecoy( X_indif , attrVals, attrSign )
%ADDDECOY Adds one or 2 decoys to the indif set
target_and_alt = randsample(1:2,2);
step_sizes = [attrVals{1}(2)-attrVals{1}(1) attrVals{2}(2)-attrVals{2}(1)];
step_range = [numel(attrVals{1}) numel(attrVals{2})];
%% Compute maximum number of step sizes
% dist is how much better the target is (+ is better)
K = numel(attrSign);
dist = attrSign .* (X_indif(target_and_alt(1),:) - X_indif(target_and_alt(2),:))
better = dist > 0;
max_step_sizes = floor(step_range / 2);
for k = 1:K
    if better(k)
        max_step_sizes(k) = dist(k) / step_sizes(k);
    elseif attrSign(k) > 0
        max_step_sizes(k) = X_indif(target_and_alt(1),k)/step_sizes(k);
    end
end
step_limit = min(max_step_sizes) - 1;


%% Add decoy
decoys = [X_indif(target_and_alt(1),:);X_indif(target_and_alt(1),:)];
try
    num_step = randsample(floor(step_limit/2):step_limit,1);
catch
    sca;
    error('Call the supervisor');
end
decoys(1,:) = decoys(1,:) - attrSign .* (num_step+1) .* step_sizes; %Remove num_step+1 stepsize from decoy attributes
decoys(2,:) = decoys(2,:) - attrSign .* num_step .* step_sizes; %Remove num_step stepsize from decoy attribute
X_full = [X_indif(1,:);X_indif(2,:);decoys];
X_full = max(X_full,0);

% decoy_list = randsample(3:4,2);
x_list_single = randsample(1:3,3);
x_list_double = randsample(1:4,4);

X_1decoy = X_full(x_list_single,:);
X_2decoy = X_full(x_list_double,:);

TargetAndAltX1 = [find(target_and_alt(1)==x_list_single) find(target_and_alt(2)==x_list_single) ];
TargetAndAltX2 = [find(target_and_alt(1)==x_list_double) find(target_and_alt(2)==x_list_double) ];

end

