function [ X_optim ] = PathToMin(objective,X0, attrVals )
%PATHTOMAX find the NYX path to maximum on a grid.
K=size(attrVals,1);
startvalue = objective(X0);
num_dir = numel(X0);
X_walker = X0';
current_dir = 1;
last_move = num_dir;
moving = true;
which_k = @(elt) (mod(elt,K)==0)*K + mod(elt,K);
% find the number of levels for each k
num_levels = 1:K;
for k=1:K
    num_levels(k) = numel(attrVals{k});
end
%% returns highest ascent direction
    function [best_dir,move_forward,is_min,X_walker,startvalue,current_k] = find_dir(X_walker,startvalue)
        value_dir = inf(num_dir,2);
        prop_level = zeros(num_dir,2);
        current_k = 1;
        parfor direction=1:num_dir
            current_k = (mod(direction,K)==0)*K + mod(direction,K);
            level_index = find(X_walker(direction)==attrVals{current_k});
            if level_index < num_levels(current_k)
                X_prop = X_walker;
                prop_level(direction,2) = level_index + 1;
                X_prop(direction) = attrVals{current_k}(prop_level(direction,2));
                propvalue = objective(X_prop');
                if propvalue < startvalue
                    value_dir(direction,2) = propvalue;       
                end     
            end
        end
        parfor direction=1:num_dir
            current_k = (mod(direction,K)==0)*K + mod(direction,K);
            level_index = find(X_walker(direction)==attrVals{current_k});
            if level_index > 1
                X_prop = X_walker;
                prop_level(direction,1) = level_index - 1;
                X_prop(direction) = attrVals{current_k}(prop_level(direction,1));
                propvalue = objective(X_prop');
                if propvalue < startvalue
                    value_dir(direction,1) = propvalue;
                end     
            end
        end
        [best_dir,forward] = find(min(value_dir(:))==value_dir,1);
        startvalue = value_dir(best_dir,forward);
        move_forward = forward == 2;
        is_min = isinf(min(min(value_dir)));
        if ~is_min
            current_k = (mod(best_dir,K)==0)*K + mod(best_dir,K);
            X_walker(best_dir) = attrVals{current_k}(prop_level(best_dir,forward));
        end
    end
%% 
is_min = false;
while ~is_min
    [current_dir,move_forward,is_min,X_walker,startvalue,current_k] = find_dir(X_walker,startvalue);
    level_index = find(X_walker(current_dir)==attrVals{current_k});
    % try to move forward
    if move_forward && ~is_min
        move_success = true;
        while move_success
           X_prop = X_walker;
           prop_level = level_index + 1;
           if prop_level <= num_levels(current_k)
               X_prop(current_dir) = attrVals{current_k}(prop_level);
               propvalue = objective(X_prop');
               if propvalue < startvalue
                   move_forward = true;
                   X_walker = X_prop;
                   startvalue = propvalue;
                   last_move = current_dir;
                   level_index = prop_level;
               else
                   move_success = false;
               end
           else
               move_success = false;
           end
        end
    elseif ~move_forward && ~is_min
    % try to move backward if forward did not work
        move_success = true;
        while move_success
           X_prop = X_walker;
           prop_level = level_index - 1;
           if prop_level >= 1
               X_prop(current_dir) = attrVals{current_k}(prop_level);
               propvalue = objective(X_prop');
               if propvalue < startvalue
                   X_walker = X_prop;
                   startvalue = propvalue;
                   last_move = current_dir;
                   level_index = prop_level;
               else
                   move_success = false;
               end
           else
               move_success = false;
           end
        end
    end
end
X_optim=X_walker';

end

