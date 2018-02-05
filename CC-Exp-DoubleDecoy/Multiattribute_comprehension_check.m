Screen('TextSize', mainwin, 18);
DrawFormattedText(mainwin,'This Lab is considered as an exam room.\n Please turn your phone off and put it in your bag. \nPlease remain silent at all time.\nPress Enter to continue.', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
Screen('Flip',mainwin);
KbStrokeWait;

Screen('TextSize', mainwin, 18);
DrawFormattedText(mainwin,'You are going to see 3 example questions.\n They will not be recorded and you can take as much time as you want.\nPress Enter to begin.', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
Screen('Flip',mainwin);
KbStrokeWait;

%% Run Comprehension check example
for num_option= 2:4
    X_optim = OptimDesign( Particles, num_option, attrVals, attrSign, 'Shannon' ); 
    X_check = reshape( X_optim,[param.K size(X_optim,2)/param.K] )';
    J = size(X_check,1);
    K = size(X_check,2);

    % create text
    X_text = cell(J,1);
    for j=1:J
        X_text{j} = '';
        for k=1:K
            X_text{j} = char(strcat(X_text{j},attrNames{k},{' : '},num2str(X_check(j,k),' %.2f'),'\n'));
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
    respToBeMade = true;
    
    
    while respToBeMade
        
        [keyIsDown,secs, keyCode] = KbCheck;
        
        
        if  keyCode(LeftKey) || keyCode(LeftNKey)
            choice = 1;
            respToBeMade = false;
            SelectedText = textLeft;
            selectedTextPos = leftRect;
        elseif keyCode(UpKey)  || keyCode(UpNKey)
            if J > 2
                choice = 3;
                respToBeMade = false;
                SelectedText = textTop;
                selectedTextPos = topRect;
            end
        elseif keyCode(DownKey) || keyCode(DownNKey)
            if J > 3
                choice = 4;
                respToBeMade = false;
                SelectedText = textBot;
                selectedTextPos = botRect;
            end
        elseif keyCode(RightKey) || keyCode(RightNKey)
            choice = 2;
            respToBeMade = false;
            SelectedText = textRight;
            selectedTextPos = rightRect;
        elseif keyCode(escKey)
            sca;
            error('leaving experiment before completion.')
        end
        
        
    end
    % Feedback
    Screen('TextSize', mainwin, 18);
    DrawFormattedText(mainwin,'You selected the following option.\nPress Enter to continue.', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
    [nx, ny, bbox] = DrawFormattedText(mainwin,SelectedText, 'center', 'center', 0, [], [], [], [1.5], [],selectedTextPos);
    Screen('Flip',mainwin);
    KbStrokeWait;
end

Screen('TextSize', mainwin, 18);
DrawFormattedText(mainwin,'The experiment is now ready to begin.\nPress Enter when you are ready.', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
Screen('Flip',mainwin);
KbStrokeWait;