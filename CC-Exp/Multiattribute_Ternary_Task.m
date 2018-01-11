
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



%%
timerVal = tic;
timeRecords = struct;
timeRecords.show = zeros(opt_num_quest,1);
timeRecords.answer = zeros(opt_num_quest,1);

%%
for obs=1:opt_num_quest ;
    
    X_optim = OptimDesign( Particles, num_option_list(obs), attrVals, attrSign, 'ADO' );
    Xs{obs,1} = reshape( X_optim,[param.K size(X_optim,2)/param.K] )';
    J = size(Xs{obs,1},1);
    K = size(Xs{obs,1},2);
    
    % create text
    X_text = cell(J,1);
    for j=1:J
        X_text{j} = '';
        for k=1:K
            X_text{j} = char(strcat(X_text{j},attrNames{k},{' : '},num2str(Xs{obs,1}(j,k),' %.2f'),'\n'));
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
    if Debug==1
        Screen('FrameRect', mainwin ,[0 0 255],[leftRect;rightRect;topRect]',1);
        Screen('FillOval', mainwin ,[0 0 255],[center-[2 2] center+[2 2]]);
        Screen('DrawText',mainwin,'Debug is on',10,10,textcolor);
    end
    
    Screen('Flip',mainwin);
    timeRecords.show(obs) = toc(timerVal);
    respToBeMade = true;
    
    
    while respToBeMade
        
        [keyIsDown,secs, keyCode] = KbCheck;
        
        
        if  keyCode(LeftKey)
            choice = 1;
            respToBeMade = false;
        elseif keyCode(UpKey)
            if J > 2
                choice = 3;
                respToBeMade = false;
            end
        elseif keyCode(DownKey)
            if J > 3
                choice = 4;
                respToBeMade = false;
            end
        elseif keyCode(RightKey)
            choice = 2;
            respToBeMade = false;
        elseif keyCode(escKey)
            sca;
            error('leaving experiment before completion.')
        end
        
        
    end
    timeRecords.answer(obs) = toc(timerVal);
    ChoiceList = [ChoiceList; choice];
    
    Screen('TextSize', mainwin, 25);
    [nx, ny, bbox] = DrawFormattedText(mainwin,'Loading', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
    Screen('Flip',mainwin);
    %
    [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
    BF = exp(Particles{2}.log_marg_like - Particles{1}.log_marg_like)
    
    
    
end

%Save Data
filenamesave = ['data' filesep 'Optim-' num2str(subid) '-' datestr(datetime('now'),'yyyy-mm-dd-HH.MM.SS') '.mat'];
save(filenamesave,'Xs','ChoiceList','Particles','BF','timeRecords','subid', 'subage', 'gender');
%Send file through FTP
ts = ftp('ftp.remidaviet.com','MatlabData@davserv.net','t6ubJ3Mn1qQm7gC9AAJb');
mput(ts,filenamesave);
close(ts);


%Thank you screen
Screen('TextSize', mainwin, 18);
DrawFormattedText(mainwin,'Thank you !\nThe experiment is now over, you can go to the experimenter and collect your participation fees.', 'center', 'center', 0, [], [], [], [1.5], [],centerRect);
Screen('Flip',mainwin);
    
KbStrokeWait;
sca;

