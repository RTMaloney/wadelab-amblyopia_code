function Batch_run_MID_thresholds_2IFC (SubjCode, CohOrCont)

% Sort out the different conditions.
% We want 7 amplitudes, 9 frequencies, and 2 cues: CD or IOVD
% So, 7 * 9 * 2 = 126
% Then *2 again for contrast/coherence.
%
% We include the flag 'CohOrCont' to decide whether to collect contrast or coherence (done separately).
%
% R Maloney, 16 January 2016
% 
% NOTE: this is the modified version with reduced conditions, for running on the Viewpixx2/PC in Room B116.
% R Maloney, September 2016

%Just check that the coherence/contrast flag has been entered correctly, fix if not
if strcmp(CohOrCont, 'coh') || strcmp(CohOrCont, 'cont')
    %that's correct, so do nothing
elseif strcmp(CohOrCont, 'Coh')
    CohOrCont = 'coh';
elseif strcmp(CohOrCont, 'Cont')
    CohOrCont = 'cont';
else
    %any other combination will be wrong!
    error ('you didn''t enter the coherence/contrast flag correctly (coh or cont)!')
end

% If it doesn't exist, then we need to set up all the variables.
NumAmplitudes = 5; %7;
NumFrequencies = 5; %9;
NumCues = 2;
% all conditions in a single repetition (ie a single run for each condition) ignoring threshold type:
NumRuns = NumAmplitudes * NumFrequencies * NumCues;
RunRepeats = 3; %how many repeat runs of each condition to perform

%Set up cue type
Cues = [repmat({'CD'}, NumAmplitudes * NumFrequencies, 1); ...
    repmat({'IOVD'}, NumAmplitudes * NumFrequencies, 1)];

%Set up signal type: not necessary.
% CohOrCont = [repmat({'coh'}, NumAmplitudes * NumFrequencies, 1); ...
%     repmat({'cont'}, NumAmplitudes * NumFrequencies, 1)];
% CohOrCont = [CohOrCont; CohOrCont];

% Set up frequencies. Multiplying the result by 100 & then dividing by 100 gives the values to 2 dec points
% (important for when the values are used in the data file names).
frequencies = round(logspace(log10(0.5), log10(8), NumFrequencies) * 100) / 100;
frequencies = repmat(reshape(repmat(frequencies,NumFrequencies,1), 1, NumRuns/2)',2,1);

% Now do the same for amplitudes:
CDamplitudes = round(logspace(log10(1.67), log10(16.67), NumAmplitudes) * 100) / 100;
%Let's just fix a few of the overlapping values so they are the same in the 2 cues:
%CDamplitudes(5) = 7.75; % was 7.74 => not necessary now
CDamplitudes = repmat(CDamplitudes',NumFrequencies,1);

IOVDamplitudes = round(logspace(log10(1.67), log10(167), NumAmplitudes) * 100) / 100;
%Let's just fix a few of the overlapping values so they are the same in the 2 cues:
IOVDamplitudes(3) = 16.67;  %was 16.7
IOVDamplitudes = repmat(IOVDamplitudes',NumFrequencies,1);

% Put both cue amplitudes together now in a single vector:
amplitudes = [CDamplitudes; IOVDamplitudes];

% Determine whether this subject has already started testing or not
RunOrderFileName = fullfile('data', [SubjCode, '_MID_thresholds_RunOrder_', CohOrCont, '.mat']);

% load this subject's run order if it does not exist.
if exist(RunOrderFileName, 'file') == 2
    load (RunOrderFileName) %load the existing order file for that subject
else
    ExpIncrement = 0;
    % Work out the random order for this particular subject.
    % We repeat it 3 times for each of the 3 run repeats
    RunOrder = [randperm(NumRuns), randperm(NumRuns), randperm(NumRuns)];
    save(RunOrderFileName, 'RunOrder', 'ExpIncrement')
end


while ExpIncrement < NumRuns * RunRepeats
    
    R = ExpIncrement + 1;
    
    % Determine which Run repetition we're up to, depending on how far through testing we are.
    % We will finish a single repeat for each condition before moving on to the next repeat (of 3).
    if R <= NumRuns  %1st repeat
        runN = 1;
    elseif R > NumRuns && R <= NumRuns*2 %2nd repeat
        runN = 2;
    elseif R > NumRuns*2 && R <= NumRuns*3 % 3rd repeat
        runN = 3;
    end
    
    %--------------------------------%
    %        Display stimulus:
    %--------------------------------%
    % Now we can take the current stimulus and display it continuously on the screen so the subject knows what to look for
    % Most of the code below is just copied directly from the MID threshold paradigm (& many others that display the stimuli...)
    
    % Define some of the display parameters:
    PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)
    
    % If you're using the PROpixx or Viewpixx
    UsingVP = true;
    useHardwareStereo = false;
    
    % Flag to draw the annulus
    DrawAnnulus = true;
    
    %Choose the screen: it is usually the max screen no. available.
    %Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
    %So if we're on that machine, need to -1 from the screen number:
    [~, CompName] = system('hostname'); %find out the computer name
    if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... %normal strcmp not working here, can't figure out why...
            && (length(Screen('Screens'))>1) %and there is more than 1 display connected...
        WhichScreen = max( Screen( 'Screens' ) )-1;
    else
        WhichScreen = max( Screen( 'Screens' ) ); %should be right for any other machine!
    end
    
    screenRect = Screen('Rect',WhichScreen); %get the screen resolution.
    centreX = screenRect(3)/2;
    centreY = screenRect(4)/2;
    RefreshRate = Screen('NominalFrameRate', WhichScreen);
    
    %%%%-------------------------------%%%%
    %       define response keys
    %%%%-------------------------------%%%%
    
    %Unify the keyboard names in case we run this on a mac:
    KbName('UnifyKeyNames')
    % Finally, press a q at any time to quit the program.
    RespQuit = KbName('q'); % q to quit.
    
    if ismac
        % Also get the Gamepad indices:
        GamepadIndex = GetGamepadIndices;
    end
    
    %%%%-------------------------------%%%%
    %           Set up stimulus
    %%%%-------------------------------%%%%
    
    % Timing:
    % the effective frame rate PER EYE: since each eye is stimulated on successive video frames
    % Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
    PeyeFR = RefreshRate/2; % Per eye f.r. is always half the absolute f.r.
    
    % Contrast:
    pdc = 1; % the peak dot contrast
    
    % Define the dot texture, a square-shaped sheet of dots.
    % Make the texture the same size as the height of the screen
    % (well really we are using many little dot textures, but 'imsize' defines the region they inhabit)
    % Note, this should really be the same for both CD and IOVD...
    imsize = screenRect(4);
    
    % specify dot size:
    dot_sigma_in_degrees = 0.05;            % 0.064; %size of SD of the dot profile in degs/vis angle
    dot_sigma = dot_sigma_in_degrees * PPD; % sigma in pixels
    dotsize = round(dot_sigma * 10);        % make the dots some multiple of sigma
    % NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
    % It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.
    
    % make the dot profile:
    x = (1-dotsize)/2:(dotsize-1)/2;
    %x = x/PPD; %rescale x into vis angle
    [x,y] = meshgrid(x,x);
    y = -y;
    [a,r] = cart2pol(x,y);
    
    % This gives us a white-peaked dot (+ve contrast polarity)
    env = exp(-r.^2/(2*dot_sigma.^2)); % gaussian window
    env = env./max(max(abs(env)));     % normalize peak to +/- 1
    env2 = -env;                       % make a negative version to give us a black dot (-ve contrast polarity)
    
    % set up the raised cosine annular window.
    % specify parameters for the annulus:
    inrad = PPD * 1;     % inner radius of annulus (in pixels), for fixation spot
    outrad = PPD * 12/2; % outer radius of annulus (in pixels)
    % define extent of spatial raised cosine at edge of aperture (in pixels)
    cos_smooth = dotsize; % make it one dot size wide
    % This should plonk the window in the middle of the matrix, which is what we want
    imsize2 = imsize*2;   % double the texture size
    x0 = (imsize2+1)/2;
    y0 = (imsize2+1)/2;
    J = ones(imsize2);
    for (ii=1:imsize2)
        for (jj=1:imsize2)
            r2 = (ii-x0)^2 + (jj-y0)^2;
            if (r2 > outrad^2)
                J(ii,jj) = 0;
            elseif (r2 < inrad^2)
                J(ii,jj) = 0;
            elseif (r2 > (outrad - cos_smooth)^2)
                J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
            elseif (r2 < (inrad + cos_smooth)^2)
                J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
            end
        end
    end
    
    %%%%-------------------------------%%%%
    %       Set up fixation
    %%%%-------------------------------%%%%
    
    % Set up the fixation cross or spot:
    % This is drawn directly to the screen using Screen('FillRect')
    crossWidth = 2;
    crossHeight = 10;
    fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);
    
    ForcedQuit = false; % this is a flag for the exit function to indicate whether the program was aborted
    
    try % Start a try/catch statement, in case something goes awry with the PTB functions
        
        %----------------------------
        %     Set up the screen
        %----------------------------
        
        % initialization of the display
        AssertOpenGL;
        % Open PTB onscreen window: We request a 32 bit per colour component
        % floating point framebuffer if it supports alpha-blending. Otherwise
        % the system shall fall back to a 16 bit per colour component framebuffer:
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        % required for gamma correction through the PsychImaging pipeline:
        PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
        % Set the color range to be normalised between 0 and 1 (rather than 0-255):
        PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
        
        % Open an on screen (grey) window and configure the imaging pipeline
        % Info about the 'blueline' mechanism for synching to the 3D glasses:
        % There seems to be a blueline generation bug on some OpenGL systems.
        % SetStereoBlueLineSyncParameters(windowPtr, windowRect(4)) corrects the
        % bug on some systems, but breaks on other systems.
        % We'll just disable automatic blueline, and manually draw our own bluelines!
        if useHardwareStereo
            [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5, [], [], [], 1); %flag of 1 engages stereomode
            SetStereoBlueLineSyncParameters(win, windowRect(4)+10);
        else
            [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
        end
        
        % Initialise the Vpixx device:
        if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
            
            PsychImaging('PrepareConfiguration');
            PsychImaging('AddTask', 'General', 'UseDataPixx');
            Datapixx('Open');
            % The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
            Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
            Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
            Datapixx('EnableVideoStereoBlueline');
            Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
            % Prepare for the M16 mode for 16 bit luminance output precision: important when varying contrast
            % PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
            
            % Now also prepare the overlay. This is to allow drawing of the blue lines (I think)?
            %overlay = PsychImaging('GetOverlayWindow', win);
            
            if Datapixx('IsViewpixx3D') % If it's the Viewpixx3D
                
                % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
                % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
                % The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
                % The values are the K channel (all RGB) gamma function fits as measured through the goggles, 25/8/16.
                % The spectral measurements were made using the Jaz spectrometer
                K_gamma_left = 3.7537;
                K_gamma_right = 3.7282;
                % We'll use the average of the right and left gamma values
                PsychColorCorrection('SetEncodingGamma', win, 1/mean([K_gamma_left, K_gamma_right]));
                
                Datapixx('EnableVideoLcd3D60Hz');
                subjectData.DisplayType = 'Viewpixx3D'; % set aside the device type for reference
                Datapixx('RegWr');
                
            elseif Datapixx('IsPropixx') % if it's the Propixx DLP projector
                
                subjectData.DisplayType = 'PROpixx';        % set aside the device type for reference
                Datapixx('SetPropixxDlpSequenceProgram',0); % set to normal RGB video processing for driving the LEDs & DLP MMDs
                %Datapixx('RegWr');
                Datapixx('RegWrRd');                        % seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
            end
        end
        % Define the 'blue line' parameters
        blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
        blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
        blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
        blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
        % No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
        
        % Query the screen refresh rate:
        ifi = Screen('GetFlipInterval',win); % in sec
        
        % Set the alpha-blending:
        % We want a linear superposition of the dots should they overlap:
        % Just like the Gabors in GarboriumDemo.m (see there for further info).
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
        % We also want alpha-blending for smooth (anti-aliased) dots...
        % not sure how this will conflict with the above command
        % about the linear superposition of dots... but it doesn't seem to cause problems
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        % Make the textures for the dots, 1 for black and white
        Dots(1) = Screen('MakeTexture',win,env,[],[],2);  % white dot
        Dots(2) = Screen('MakeTexture',win,env2,[],[],2); % black dot
        
        % Generate the annulus texture:
        AnnImages = 0.5.*ones(imsize2,imsize2,2); % specify RGB matrix of annulus, grey
        AnnImages(:,:,2) = 1.*J;                  % specify Alpha channel of annulus
        annulus = Screen('MakeTexture',win,AnnImages,[],[],2);
        
        % Preallocate array with destination rectangles:
        % This also defines initial dot locations
        % for the very first drawn stimulus frame:
        texrect = Screen('Rect', Dots(1));
        
        % If motion coherence < 1.0, we will need to adjust the dot x positions of the newly-determined noise dots so they are in the centre of the screen.
        % If we don't do this they will be shifted off to the left of screen. This is because imsize is smaller than the full width (ie x distance) of the screen
        AdjustXBy = (screenRect(3) - screenRect(4))/2; % shift X dot positions by this amount to centre them, since image size <screenRect(3)
        
        %---------------------------------------------------------%
        %           Load pre-generated dot positions
        %---------------------------------------------------------%
        
        % Determine and load the sequence parameters for the current interval of the current trial.
        % Note that trial_order is entered as a scalar to define the current interval where 0=noise, 1=signal
        % SetUpDotPositions is a nested sub-function defined below
        tic
        [DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
            (Cues(RunOrder(R)), 1, 1, 0, CohOrCont, amplitudes(RunOrder(R)), frequencies(RunOrder(R)), AdjustXBy, PPD, Dots, texrect, imsize);
        toc
        
        %%%%---------------------------------------------------------------%%%%
        %                     Present the stimuli!
        %%%%---------------------------------------------------------------%%%%
        
        % Display the welcome screen and wait for the user to begin.
        Screen('TextFont',win, 'Arial');
        Screen('TextSize',win, 24);
        vbl = Screen('Flip',win); %sync vbl to start time
        F = 1; % F increments continuously across frames.
        ButtonPressed = 0;
        WaitSecs(0.2) 
        KbCheck(); % take a quick KbCheck to load it now & flush any stored events
        
        while ~ButtonPressed % Loop until any button is pressed, either keyboard or gamepad
            
            % Select left-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 0);   % flag of 0= left eye
                DrawFormattedText(win,['Welcome ' SubjCode , ...
                    '. \n \nWhich presentation contains the smooth 3D motion below?', ...
                    ['\n \nThis is ' char(Cues(RunOrder(R))) ' at ' num2str(frequencies(RunOrder(R))) ' Hz, to ' num2str(amplitudes(RunOrder(R))) ' arcmin.'], ...
                    '\n \nPress '' q '' on the keyboard to quit at any time.', ...
                    '\n \nPress any button/key to begin.'], ...
                    windowRect(4)/4, windowRect(3)/4, 0);
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw left eye stimulus:
            %%%%------------------------------------------------%%%%
            
            % Draw dots:
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); % turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('DrawTextures',win,DotsIdxL, ...             % texture index
                [], ...
                dstRectsL(:,:,mod(F,FramesFullCycle)+1), ...    % select appropriate dot positions for current frame
                [],[],pdc)                                      % final argument is contrast
            
            % Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            % Draw the black fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
            
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            % Draw blue lines :
            Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
            Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
            
            DrawFormattedText(win,['Welcome ' SubjCode , ...
                '. \n \nWhich presentation contains the smooth 3D motion below?', ...
                ['\n \nThis is ' char(Cues(RunOrder(R))) ' at ' num2str(frequencies(RunOrder(R))) ' Hz, to ' num2str(amplitudes(RunOrder(R))) ' arcmin.'], ...
                '\n \nPress '' q '' on the keyboard to quit at any time.', ...
                '\n \nPress any button/key to begin.'], ...
                windowRect(4)/8, windowRect(3)/12, 0);
            
            % Select right-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 1);
                DrawFormattedText(win,['Welcome ' SubjCode , ...
                    '. \n \nWhich presentation contains the smooth 3D motion below?', ...
                    ['\n \nThis is ' char(Cues(RunOrder(R))) ' at ' num2str(frequencies(RunOrder(R))) ' Hz, to ' num2str(amplitudes(RunOrder(R))) ' arcmin.'], ...
                    '\n \nPress '' q '' on the keyboard to quit at any time.', ...
                    '\n \nPress any button/key to begin.'], ...
                    windowRect(4)/8, windowRect(3)/12, 0);
            else
                Screen('DrawingFinished', win);
                [vbl , ~ , ~, ~] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); % update display on next refresh (& provide deadline)
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw right eye stimulus:
            %%%%------------------------------------------------%%%%
            
            % Draw dots:
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); % turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('DrawTextures',win,DotsIdxR, ...             % texture index
                [], ...
                dstRectsR(:,:,mod(F,FramesFullCycle)+1), ...    % select appropriate dot positions for current frame
                [],[], pdc)                                     % final argument is contrast
            
            % Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            % Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            %Draw the fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
            
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            % Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
            Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
            
            DrawFormattedText(win,['Welcome ' SubjCode , ...
                '. \n \nWhich presentation contains the smooth 3D motion below?', ...
                ['\n \nThis is ' char(Cues(RunOrder(R))) ' at ' num2str(frequencies(RunOrder(R))) ' Hz, to ' num2str(amplitudes(RunOrder(R))) ' arcmin.'], ...
                '\n \nPress '' q '' on the keyboard to quit at any time.', ...
                '\n \nPress any button/key to begin.'], ...
                windowRect(4)/8, windowRect(3)/12, 0);
            
            Screen('DrawingFinished', win);
            
            [vbl , ~ , ~, ~] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
            
            % Increment frame counters on a PER EYE basis:
            F = F+1;   % For stimulus position index
            
            % Check for button or key presses to resume experiment
            if ismac
                %  Check both the keyboard and the gamepad;
                % if 'q' is pressed, abort
                [KeyIsDown, ~, keyCode] = KbCheck();
                if KeyIsDown % a key has been pressed
                    if keyCode(RespQuit)
                        ForcedQuit = true
                        ExitGracefully(UsingVP, ForcedQuit)
                        %if any other button on the keyboard has been pressed
                    else
                        ButtonPressed = 1;
                    end
                else %if not, check the response pad
                    [ButtonPressed, ~, ~] = KbCheck(GamepadIndex);
                end
            else % if we're not on the mac:
                [KeyIsDown, ~, keyCode] = KbCheck();
                if KeyIsDown
                    if keyCode(RespQuit)
                        ForcedQuit = true
                        ExitGracefully(UsingVP, ForcedQuit)
                        %if any other button on the keyboard has been pressed
                    else
                        ButtonPressed = 1;
                    end
                end
            end
            
        end % end of stimulus demo loop
        
    catch MException
        
        rethrow (MException)
        psychrethrow(psychlasterror)
        ExitGracefully (UsingVP, ForcedQuit)
        error('Error!')
        
    end % End of try/catch statement
    
    Screen('CloseAll') % close down the demo screen between runs
    
    % Button has been pressed, end intro demo & begin experiment
    MID_Cue_thresholds_2IFC_M16 (SubjCode, char(Cues(RunOrder(R))), CohOrCont, frequencies(RunOrder(R)), amplitudes(RunOrder(R)), runN)
    %Increment the test number, in case we need to start again
    ExpIncrement = ExpIncrement + 1
    % And save it:
    save(RunOrderFileName, 'RunOrder', 'ExpIncrement')
    
end

end %End of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitGracefully (UsingVP, ForcedQuit)
%...need to shut everything down here...

% turn off the prioritisation:
Priority( 0 ); % restore priority

if UsingVP        % close down the ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
    if Datapixx('IsViewpixx3D')
        Datapixx('DisableVideoLcd3D60Hz');
    end
    Datapixx('RegWr');
    %Datapixx('Close'); %closing it here might cause it to crash?
end

% Close down the screen:
Screen('CloseAll')

Datapixx('Close'); % closing the Datapixx here (after closing the screen) might stop it from crashing

% Bring back the mouse cursor:
ShowCursor();

% announce to cmd window if the program was aborted by the user
if ForcedQuit
    error('You quit the program!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
    (Cue, trial_order, C, Trial, CohOrCont, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize)

% This is where the pre-generated dot positions are loaded and the disparity added.
% (see GenerateIOVDdots.m & GenerateCDdots.m)
% These are picked at random (with replacement) from those already in the folder.
% Note that trial_order is entered as a scalar to define the order of the current interval, where 0 = noise, 1=signal

% Indicate how many existing sets of dot stimuli are in the subfolder 'stimuli' for loading later.
StimInFolder = 30;

% Define the max. change in disparity (horiz shift of dots) for the 2 cues:
% the numerator is in arcmin (but converted to pixels) akin to the amplitude of the sine wave
IOVD_disparity = amplitude/60 * PPD;
CD_disparity = amplitude/60 * PPD;

% Convert frequency into a string so we can load the appropriate stimulus:
IOVDstrfreq = num2str(frequency);
IOVDstrfreq(IOVDstrfreq == '.') = '_'; % replace decimal point with underscore
CDstrfreq = num2str(frequency);
CDstrfreq(CDstrfreq == '.') = '_';

% Draw a random value to load the stimuli with:
RandDraw = ceil(rand*StimInFolder);

%%%%-------------------------------------------------------------------------%%%%
%                     Shift dot trajectories: add disparity
%%%%-------------------------------------------------------------------------%%%%

% *** ---------------------------------- ***
%               Set up CD
% *** ---------------------------------- ***
if strcmp(Cue, 'CD')
    
    FileNameCD = fullfile('stimuli',['CD_dots_', CDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
    CD_stim = load(FileNameCD); % includes parameters file and dot positions matrices
    % Just cut the parameters back to a smaller struct:
    CDparams = CD_stim.DotParams;
    
    % *** set up dot coherence ***
    % Now we'll determine dot coherence according to the staircase
    if strcmp(CohOrCont, 'coh')                % If we're measuring coherence thresholds
        NumDots = round(CDparams.NumDots * C); % modify number of coherent dots by value by current staircase setting
        
        % Now that we know the number of coherent dots, set up some random indices that will determine which dots are 'randomised'
        % this is necessary to ensure that dots of both colours are made noise dots equally.
        % We also want these indices to be consistent across frames in a trial, so they are set here.
        % Remember that we want to index the NOISE dots to change (ie=1), not the signal dots,
        % so it may seem a bit counterintuitive that noise =1, coherent =0
        NoiseDotIndices = Shuffle([true(CDparams.NumDots-NumDots,1); false(NumDots,1)]); %The 1s in here should sum to IOVDparams.NumDots-NumDots: the number of noise dots.
        
    else
        NumDots = CDparams.NumDots;            % otherwise, all dots are coherent.
    end
    
    % Define the number of frames in a full cycle:
    FramesFullCycle = CDparams.FramesFullCycle;
    
    % Assign the pre-generated dot position matrices.
    % First, we duplicate the left eye dot positions to make the right eye dot positions.
    % At this stage all dots are identical for the two eyes (disparity has not been added yet).
    dot_posL = CD_stim.dot_posL;
    dot_posR = CD_stim.dot_posL;
    
    % set aside the disparity
    CDparams.disparity = CD_disparity;
    % For reference also store the disparity in Arcmin
    CDparams.DisparityInArcmin = (CD_disparity/PPD)*60;
    
    % set up a few extra parameters of the sine waves for each cue:
    % the period, in sec
    CDparams.period = 1/CDparams.frequency;
    
    % the angular frequency of the sine waves:
    CDparams.angFreq = 2 * pi * CDparams.frequency;
    
    % Length of the sine wave:
    % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    % end of the cycle and the beginning of the next cycle.
    % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    CDparams.t = linspace(CDparams.period/CDparams.FramesFullCycle, CDparams.period, CDparams.FramesFullCycle);
    
    % Now make one full cycle of the sine wave, no matter the frequency
    % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    CDparams.SineWave =  CDparams.disparity * sin(CDparams.angFreq * CDparams.t); %disparity * sin(angFreq * t)
    
    % assign dstRects: destination rect matrices
    dstRectsL = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
    dstRectsR = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
    
    % Shift dot trajectory: add disparity
    for fr = 1:CDparams.FramesFullCycle
        
        % determine dot coordinates: remember, y position does not change: horiz disparity only
        % Update dot position according to sinsoidal trajectory on each (per-eye) frame
        % Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ... % the size of the dot texture
            dot_posL(:,1,fr) - CDparams.SineWave(fr), ...  % the x positions
            dot_posL(:,2,fr))';                            % the y positions
        % Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
            dot_posR(:,1,fr) + CDparams.SineWave(fr), ...
            dot_posR(:,2,fr))';
        
        % Now, randomise the dot positions of the noise dots when coherence <1.0
        if (C < 1) && (strcmp(CohOrCont, 'coh'))
            %Left eye:
            dstRectsL(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(CDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ... % random x positions
                rand(CDparams.NumDots-NumDots,1) * imsize)';               % random y positions
            %right eye:
            dstRectsR(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(CDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ...
                rand(CDparams.NumDots-NumDots,1) * imsize)';
        end
        
    end
    
    % Now, if it's a NOISE interval, we simply suffle the frames temporally to give us the CD control stimulus.
    % This gives us a stimulus with the same distribution of disparities but no consistent/smooth
    % change in disparity over time
    % This should be all we need to change for the CD noise stimulus control
    if trial_order == 0
        ShuffledFrames = randperm(CDparams.FramesFullCycle);
        dstRectsL = dstRectsL(:,:,ShuffledFrames);
        dstRectsR = dstRectsR(:,:,ShuffledFrames);
    end
    % If it is a SIGNAL interval we don't need to change anything
    
    % Now we have the dot position indices (the dstRects), define the dot texture indices.
    % These are simply half white (Dots(1) and half black (Dots(2)) in each eye.
    DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
    DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
    
    % *** ---------------------------------- ***
    %               Set up IOVD
    % *** ---------------------------------- ***
elseif strcmp(Cue, 'IOVD')
    
    FileNameIOVD = fullfile('stimuli',['IOVD_dots_', IOVDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
    IOVD_stim = load(FileNameIOVD); % includes parameters file and dot positions matrices
    % Just cut the parameters back to a smaller struct:
    IOVDparams = IOVD_stim.DotParams;
    
    % *** set up dot coherence ***
    % Now we'll determine dot coherence according to the staircase
    if strcmp(CohOrCont, 'coh')                 % If we're measuring coherence thresholds
        NumDots = round(IOVDparams.NumDots * C); % modify number of coherent dots by value by current staircase setting
        
        % Now that we know the number of coherent dots, set up some random indices that will determine which dots are 'randomised'
        % this is necessary to ensure that dots of both colours are made noise dots equally.
        % We also want these indices to be consistent across frames in a trial, so they are set here.
        % Remember that we want to index the NOISE dots to change (ie=1), not the signal dots,
        % so it may seem a bit counterintuitive that noise =1, coherent =0
        NoiseDotIndices = Shuffle([true(IOVDparams.NumDots-NumDots,1); false(NumDots,1)]); %The 1s in here should sum to IOVDparams.NumDots-NumDots: the number of noise dots.
        
    else
        NumDots = IOVDparams.NumDots;            % otherwise, all dots are coherent.
    end
    
    % Define the number of frames in a full cycle:
    FramesFullCycle = IOVDparams.FramesFullCycle;
    
    % set aside the disparity
    IOVDparams.disparity = IOVD_disparity;
    % For reference also store the disparity in Arcmin
    IOVDparams.DisparityInArcmin = (IOVD_disparity/PPD)*60;
    
    % set up a few extra parameters of the sine waves for each cue:
    % the period, in sec
    IOVDparams.period = 1/IOVDparams.frequency;
    
    % the angular frequency of the sine waves:
    IOVDparams.angFreq = 2 * pi * IOVDparams.frequency;
    
    % Length of the sine wave:
    % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    % end of the cycle and the beginning of the next cycle.
    % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    IOVDparams.t = linspace(IOVDparams.period/IOVDparams.FramesFullCycle, IOVDparams.period, IOVDparams.FramesFullCycle);
    
    % Now make one full cycle of the sine wave, no matter the frequency
    % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    IOVDparams.SineWave =  IOVDparams.disparity * sin(IOVDparams.angFreq * IOVDparams.t); %disparity * sin(angFreq * t)
    
    % assign dstRects: destination rect matrices
    dstRectsL = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
    dstRectsR = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
    
    % Shift dot trajectory: add disparity.
    % Remember that with IOVD the dot positions for each eye are NOT identical (before the horizontal shifts are added)
    % as they are with CD, that is why we have 2 matrices for the 2 eyes
    % Dot lifetime is determined by how many dot positions are re-randomised with each frame (see GenerateIOVDdots.m)
    
    % Determine a value to modulate the sine wave by.
    % In this way the sinusoidal disparity is either subtracted or added to the dot positions with each frame.
    % This way we can make the noise/control stimulus with a very simple switch
    % The two columns in SineModulate are for all dots in each eye: col 1 for LE, col 2 for RE.
    
    if trial_order == 0 % for a noise trial
        % Make vector of 1s & -1s
        SineModulate = ones(IOVDparams.NumDots,2);
        SineModulate(1:2:end,:) = -1; % this bit will make stimuli the IOVD control: both eyes have dots moving both directions
        
        %**** need random 1/2 to be black, 1/2 white; because color is assigned in alternate manner
        % motion direction cannot also.
        % If motion direction is assigned as above, in an alternate manner, the control has dots of opposite colours moving in
        % opposite directions. We can't have that.
        % So, re-assign SineModulate. Now it is a 2 * num_dots vector, upper half are -1s, lower half +1s
        % Remember that there is no relationship between dot position on the screen & its location in the matrix, so making the SineModulate
        % matrix in this way won't matter; the motion directions will still be assigned randomly in the IOVD control.
        % And because color is assigned in an alternate manner, that means equal numbers of black/white dots will
        % be going both monocular directions.
        
        SineModulate = ones(IOVDparams.NumDots,2);
        SineModulate(1:end/2,:) = -1;
        
        % *** Not doing the below any more ...
        %(Also, define indices to select a random sample of NumDots/2 from each eye.
        %All these dots will be presented to BOTH eyes to ensure the same amount of dots per eye.
        %These indices should stay constant across frames to ensure continuous motion of the dots (so they are assigned here).
        %Equal numbers of dots are sampled from each eye.
        %To do this, we simply re-assign the random sample of dots from the L eye to the R eye.
        %The resulting R eye matrix then has 50% of dots from the L & R eyes (because 50% of dots in the R eye matrix have not been re-assigned).
        %Then we simply make LeyeMatrix = ReyeMatrix so that the same pattern is displayed to both eyes.
        %LERandDotSample = Shuffle([true(IOVDparams.NumDots/2,1); false(IOVDparams.NumDots/2,1)]); )... ***
        
    elseif trial_order % when trial_order = 1 or -1 (or any value other than 0 really)
        % subtract left eye, add right eye for MID: signal interval
        SineModulate = ones(IOVDparams.NumDots,2);
        SineModulate(:,1) = -1; %-1 * sin for the left eye only in actual IOVD; +1 * sin for right eye
    end
    
    for fr = 1:IOVDparams.FramesFullCycle
        
        % determine dot coordinates: remember, y position does not change: horiz disparity only
        % Update dot position according to sinsoidal trajectory on each (per-eye) frame
        
        % Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ...                                      % the size of the dot texture
            IOVD_stim.DotPosBinoc{1}(:,1,fr) + SineModulate(:,1) * IOVDparams.SineWave(fr), ... % x positions
            IOVD_stim.DotPosBinoc{1}(:,2,fr))';                                                 % y positions
        
        %Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
            IOVD_stim.DotPosBinoc{2}(:,1,fr) + SineModulate(:,2) * IOVDparams.SineWave(fr), ...
            IOVD_stim.DotPosBinoc{2}(:,2,fr))';
        
        % *** Don't do the below any more:
        % AS usual we present 2 sets of dots to L/R eyes,
        % But the direction of motion is alternately left or right within a single eye, & this is determined by
        % 'SineModulate' parameter....
        
        % If it's a noise interval, we take a random sample of NumDots dots from both eyes,
        % these dots are displayed to BOTH eyes simultaneously
        %         if trial_order == 0
        %             %if a noise interval, L & R are identical
        %             dstRectsR(:,LERandDotSample,fr) = dstRectsL(:,LERandDotSample,fr); %reassign 50% of dots from L eye into R eye
        %             dstRectsL(:,:,fr) = dstRectsR(:,:,fr); %now make LE pattern = RE pattern
        %         end % ***
        
        %Now, randomise the dot positions of the noise dots when coherence <1.0
        %L & R eye matrices both exist by this point in time
        %The random indices for the NOISE dots were determined above, according to IOVDparams.NumDots-NumDots
        if (C < 1) && (strcmp(CohOrCont, 'coh'))
            
            %Left eye:
            dstRectsL(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(IOVDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ... % random x positions
                rand(IOVDparams.NumDots-NumDots,1) * imsize)';               % random y positions
            %right eye:
            dstRectsR(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(IOVDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ...
                rand(IOVDparams.NumDots-NumDots,1) * imsize)';
        end
        
    end
    
    % Now we have the dot position indices (the dstRects), define the dot texture indices.
    % We use separate vectors for L & R eye just in case you want to split the eyes by colour.
    % *** This is where the dot color must be assigned in an alternate fashion for L & R eye
    % to make (mostly) sure that neighbouring dots across the eyes are of a different color
    
    % Previously we only assigned dot color in alternating fashion for IOVD, not the control.
    % This was because we previously used the same set of dots in the 2 eyes for the control,
    % so if color was assigned in an alternating fashion then the dots would flicker between the 2 colours as the
    % dots were presented in alternating fashion to the 2 eyes.
    % Now, because we use a different set of dots, it is ok for the dots to be assigned in alternating manner
    % across the eyes. This also maintains the method used to reduce binocular correlations in the new IOVD algorithm, applying it to the control also.
    
    %Not doing this now:
    %     if trial_order == 0
    %         DotsIdxL = [repmat(Dots(1),1,IOVDparams.NumDots/2), repmat(Dots(2),1,IOVDparams.NumDots/2)];
    %         DotsIdxR = [repmat(Dots(1),1,IOVDparams.NumDots/2), repmat(Dots(2),1,IOVDparams.NumDots/2)];
    %     else
    
    %Left eye: WBWBWB ....
    DotsIdxL = [repmat(Dots(1),1,IOVDparams.NumDots/2); repmat(Dots(2),1,IOVDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
    DotsIdxL = reshape(DotsIdxL,1,IOVDparams.NumDots);
    %Right eye: BWBWBW....
    DotsIdxR = [repmat(Dots(2),1,IOVDparams.NumDots/2); repmat(Dots(1),1,IOVDparams.NumDots/2)];
    DotsIdxR = reshape(DotsIdxR,1,IOVDparams.NumDots);
    %end
    
    %Or to split by colour:
    %DotsIdxL = repmat(Dots(2),1,IOVDparams.NumDots); %white dots only
    %DotsIdxR = repmat(Dots(2),1,IOVDparams.NumDots); %black dots only
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
