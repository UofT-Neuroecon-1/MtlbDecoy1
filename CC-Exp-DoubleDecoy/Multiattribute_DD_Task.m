
%%%%%%%%%%% TERNARY CHOICE TASK

%%  Screen parameters
screens = Screen('Screens');
screenNumber = max(screens);
[mainwin, screenrect] =  Screen(screenNumber, 'OpenWindow');
Screen('FillRect', mainwin, bgcolor);
center = [screenrect(3)/2 screenrect(4)/2];

%% Compute text positions
hShiftFromCenter = 500;
vShiftFromCenter = 250;
rectWidth = 500;
rectHeight = 500;
leftRect = [center - [hShiftFromCenter+(rectWidth/2) (rectHeight/2)], center - [hShiftFromCenter-(rectWidth/2) (-rectHeight/2)]];
rightRect = [center + [hShiftFromCenter+(rectWidth/2) (rectHeight/2)], center + [hShiftFromCenter-(rectWidth/2) (-rectHeight/2)]];
topRect = [ center - [(rectWidth/2) , vShiftFromCenter+(rectHeight/2)] , center + [(rectWidth/2) , -vShiftFromCenter+(rectHeight/2)] ];
botRect = [ center + [(rectWidth/2) , vShiftFromCenter+(rectHeight/2)] , center - [(rectWidth/2) , -vShiftFromCenter+(rectHeight/2)] ];
centerRect = [center - [(rectWidth/2) (rectHeight/2)], center + [(rectWidth/2) (rectHeight/2)]];

Screen('TextSize', mainwin, 30);
[nx, ny, bbox] = DrawFormattedText(mainwin,'Loading', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
Screen('Flip',mainwin);

% Savefile variables
computername = getenv('computername');
filenamesave = ['data' filesep 'DD-' computername '-' datestr(datetime('now'),'yyyy-mm-dd-HH.MM.SS') '-' num2str(subid)  '.mat'];


%% Comprehension check
run('Multiattribute_comprehension_check.m');

%%
timerVal = tic;
timeRecords = struct;
timeRecords.show = zeros(opt_num_quest,1);
timeRecords.answer = zeros(opt_num_quest,1);

%%
for obs=1:num_non_check+num_consistency_check;
    
    %% Compute Question
    if obs <= opt_num_quest
        %% Optimal Questions
        X_optim = OptimDesign( Particles, num_option_list(obs), attrVals, attrSign, 'Shannon' );
        Xs{obs,1} = reshape( X_optim,[param.K size(X_optim,2)/param.K] )';
    elseif obs <= opt_num_quest + num_double_decoy
        %% Compute Decoys
        X_indif = FindIndif( Particles, attrVals, attrSign );
        X_indif = reshape( X_indif,[param.K size(X_indif,2)/param.K] )';
        XIndifs{obs,1} = X_indif;
        [ X_1decoy, X_2decoy , TargetAndAltX1, TargetAndAltX2 ] = AddDecoy(X_indif, attrVals, attrSign );
        Xs{obs,1} = X_1decoy; %Single decoy
        Xs{list_double_decoy(obs-opt_num_quest),1} = X_2decoy; %Double decoy
        TargetAndAltX(obs,:) = TargetAndAltX1;
        TargetAndAltX(list_double_decoy(obs-opt_num_quest),:) = TargetAndAltX2;
        num_option_list = [num_option_list;J];
        %% debug
        list_double_decoy(obs-opt_num_quest)
        Xs{list_double_decoy(obs-opt_num_quest),1}
        TargetAndAltX(list_double_decoy(obs-opt_num_quest),:)
        TargetAndAltX2
    elseif obs <= opt_num_quest + 2*num_double_decoy
        J = size(Xs{obs,1},1);
    else
        %% Consistency check
        redo_obs = list_const_check(obs-num_non_check);
        J = size(Xs{redo_obs,1},1);
        num_option_list = [num_option_list;J];
        option_permuted = randperm(J);
        Xs{obs,1} =  Xs{redo_obs,1}(option_permuted,:);
        expected_choice = sum((1:J) .* (ChoiceList(redo_obs)==option_permuted));
        TargetAndAltX(obs,1) = find(TargetAndAltX(redo_obs,1)==option_permuted);
        TargetAndAltX(obs,2) = find(TargetAndAltX(redo_obs,2)==option_permuted);
    end
    
    J = size(Xs{obs,1},1);
    K = size(Xs{obs,1},2);
    num_option_list(obs) = J;
    
    % create text
    X_text = cell(J,1);
    for j=1:J
        X_text{j} = '';
        for k=1:K
            X_text{j} = char(strcat(X_text{j},attrNames{k},{' : '},num2str(Xs{obs}(j,k),' %.2f'),'\n'));
        end
    end
    
    %% show alternatives
    Screen('TextSize', mainwin, 20);
    textLeft = [X_text{1} '\nPress Left'];
    textRight = [X_text{2} '\nPress Right'];
    [nx, ny, bbox] = DrawFormattedText(mainwin,textLeft, 'center', 'center', 0, [], [], [], [1.5], [],leftRect);
    [nx, ny, bbox] = DrawFormattedText(mainwin,textRight, 'center', 'center', 0, [], [], [], [1.5], [],rightRect);
    if J > 2
        textTop = [X_text{3} '\nPress Up'];
        [nx, ny, bbox] = DrawFormattedText(mainwin,textTop, 'center', 'center', 0, [], [], [], [1.5], [],topRect);
    end
    if J > 3
        textBot = [X_text{4} '\nPress Down'];
        [nx, ny, bbox] = DrawFormattedText(mainwin,textBot, 'center', 'center', 0, [], [], [], [1.5], [],botRect);
    end

    %%
%     if Debug==1
%         Screen('FrameRect', mainwin ,[0 0 255],[leftRect;rightRect;topRect]',1);
%         Screen('FillOval', mainwin ,[0 0 255],[center-[2 2] center+[2 2]]);
%         Screen('DrawText',mainwin,'Debug is on',10,10,textcolor);
%     end
    
    Screen('Flip',mainwin);
    timeRecords.show(obs) = toc(timerVal);
    respToBeMade = true;
    
    
    while respToBeMade
        
        [keyIsDown,secs, keyCode] = KbCheck;
        
        
        if  keyCode(LeftKey) || keyCode(LeftNKey)
            choice = 1;
            respToBeMade = false;
        elseif keyCode(UpKey) || keyCode(UpNKey)
            if J > 2
                choice = 3;
                respToBeMade = false;
            end
        elseif keyCode(DownKey) || keyCode(DownNKey)
            if J > 3
                choice = 4;
                respToBeMade = false;
            end
        elseif keyCode(RightKey) || keyCode(RightNKey)
            choice = 2;
            respToBeMade = false;
        elseif keyCode(escKey)
            sca;
            error('leaving experiment before completion.')
        end
        
        
    end
    timeRecords.answer(obs) = toc(timerVal);
    ChoiceList(obs) = choice;
    % Consistency check
    if obs > num_non_check
        ConsistencyCheck(obs-num_non_check) = expected_choice == choice;
    end
    
    Screen('TextSize', mainwin, 25);
    [nx, ny, bbox] = DrawFormattedText(mainwin,'Loading', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
    Screen('Flip',mainwin);
    %
    [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
    %% Save Data
    NumDecoychoice = sum(ChoiceList(1:obs) == TargetAndAltX(1:obs,:))
    save(filenamesave,'Xs','ChoiceList','Particles','timeRecords','subid', 'subage', 'gender','ConsistencyCheck','list_const_check','NumDecoychoice','TargetAndAltX');
    
    
end




%Thank you screen
Screen('TextSize', mainwin, 18);
DrawFormattedText(mainwin,'Thank you !\nThe experiment is now over, you can go to the experimenter for the debriefing.', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
Screen('Flip',mainwin);
    
KbStrokeWait;
sca;

%Send file through FTP
if ~Debug
    ts = ftp('ftp.remidaviet.com','MatlabData@davserv.net','t6ubJ3Mn1qQm7gC9AAJb');
    mput(ts,filenamesave);
    close(ts);
end

