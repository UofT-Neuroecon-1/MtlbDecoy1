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
%% 
while moving
    current_k = which_k(current_dir);
    level_index = find(X_walker(current_dir)==attrVals{current_k});
    % try to move forward
    move_forward = false;
    if level_index < num_levels(current_k)
        move_success = true;
        while move_success
           X_prop = X_walker;
           prop_level = level_index + 1;
           if prop_level <= num_levels(current_k)
               X_prop(current_dir) = attrVals{current_k}(prop_level);
               propvalue = objective(X_prop');
               if propvalue < startvalue
%                    if last_move ~= current_dir
%                        fprintf("%d",current_dir);
%                    end
%                    fprintf("+");
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
    end
    % try to move backward if forward did not work
    if level_index > 1 && ~move_forward
        move_success = true;
        while move_success
           X_prop = X_walker;
           prop_level = level_index - 1;
           if prop_level >= 1
               X_prop(current_dir) = attrVals{current_k}(prop_level);
               propvalue = objective(X_prop');
               if propvalue < startvalue
%                    if last_move ~= current_dir
%                        fprintf("%d",current_dir);
%                    end
%                    fprintf("-");
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
%     if last_move == current_dir
%         fprintf("0\n");
%     end
    %find next dimension
    if current_dir == num_dir
       current_dir= 1;
    else
       current_dir = current_dir + 1; 
    end
    %check if end loop
    moving = current_dir ~= last_move;   
end
X_optim=X_walker';

end

