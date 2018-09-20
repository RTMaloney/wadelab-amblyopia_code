function MID_Full_Cue_thresholds_2IFC_M16_ambl (SubjCode, CohOrCont, frequency, amplitude, runN)

% Frequency should be in Hz
% Amplitude should be in ARCMIN!
%
% Run the 2IFC paradigm to measure either contrast or coherence thresholds for the FULL cue condition,
% containing both the IOVD and CD cues.
% Loads pre-generated stimuli of appropriate frequency.
%   CohOrCont = String; 'coh' for coherence, 'cont' for contrast thresholds.
%   SubjCode = Unique subject identifier, as a string, eg 'RM'
%   runN = current run for this freq/ampl/threshold type.
%
% Allows for use of a USB-HID gamepad to receive subject responses (more ergonomic than keyboard).
% Note that the Gamepad functions do not work on Windows, but it will fall back on
% expecting the keyboard if run on a PC (ie for testing/development).
% Keyboard responses can also still be made if running on the mac; can use either.
% Also put fixation_cross subfunction at the bottom
%
% Version _M16: modified to make use of the DATAPixx M16 (16-bit monochrome) video mode using the PsychToolbox imaging pipeline.
%
% R Maloney, September 2016

% Define some of the display parameters:
PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)

% If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false;

% Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

% Draw the fixation lock/rings
DrawRings = true;

% Number of trials per staircase
num_trials = 30;

% Flag to give audio feedback?
GiveFeedback = false;

%%%%---------------------------------------%%%%
%       first, create a unique file name
%%%%---------------------------------------%%%%

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

% Set up file name
% Enter the frequency and amplitude into file name
strfreq = num2str(frequency);
strfreq(strfreq == '.') = '_'; % replace decimal point with underscore
strampl = num2str(amplitude);
strampl(strampl == '.') = '_';
fName = fullfile('data', [SubjCode, '_FULLcue_', CohOrCont, '_', ...
    strfreq, 'Hz_', strampl, 'arcmin_', num2str(runN), '.mat']);

%Check that it doesn't already exist
if exist( fName, 'file' )
    userResponse = input( 'WARNING: File exists. Overwrite? Enter y or n: ', 's' );
    if ~strcmp( userResponse, 'y' )
        error('Aborting function' );
    end
end

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

if ~ismac
    jheapcl; %clear the java heap space.
end

%Set aside some of the initial parameters:
parameters.Cue = 'FULL';
if strcmp(CohOrCont, 'cont')
    parameters.SignalVaried = 'contrast';
elseif strcmp(CohOrCont, 'coh')
    parameters.SignalVaried = 'coherence';
end

%%%%-------------------------------%%%%
%       define response keys
%%%%-------------------------------%%%%

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

% Subjects can press '1' for signal in interval 1: either on the keyboard or the 1 on the num keypad are fine.
Interval1 = [KbName('1!'), KbName('1')];

% Or press '2' for signal in interval 2.
Interval2 = [KbName('2@'), KbName('2')];
% Finally, press a q at any time to quit the program.
RespQuit = KbName('q'); % q to quit.

% Add the keycodes for the appropriate buttons on the gamepad, if we're
% using the mac setup.
% Now, at time of writing we are using the 'Afterglow gamepad for PS3' device.
% We want to use the L2 and R2 buttons on the top left & right of the pad,
% for the 1st and 2nd intervals, respectively.
% Because the gamepad registers as a HID device, like a keyboard, we can
% define the button responses in the same way as we would a keyboard response.
% On the Afterglow gamepad, L2 = button 7; R2 = button 8.
% They are equivalent to KbName('d') & KbName('e'), on the keyboard, respectively.
% This of course means that if a subject presses 'd' or 'e' on the keyboard
% (as well as 1 or 2), that will register as a response, so we will just
% have to warn against that.

if ismac
    Interval1 = [Interval1, KbName('d')];
    Interval2 = [Interval2, KbName('e')];
    % Also get the Gamepad indices:
    GamepadIndex = GetGamepadIndices;
end


%%%%-------------------------------%%%%
%           Set up stimulus
%%%%-------------------------------%%%%

% Timing:
TrialDuration = 2; % duration of the stimulus, IN SEC

% the effective frame rate PER EYE: since each eye is stimulated on successive video frames
% Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
PeyeFR = RefreshRate/2; % Per eye f.r. is always half the absolute f.r.

% stimuli for each eye are indexed with the same number: they increment on a per-eye frame rate basis
% But PTB refresh 'flips' of course happen at the machine frame rate: there are 2 flips within a single stimulus loop: one for each eye.
% PeyeFrames is the number of (per-eye) frames needed to achieve the trial duration.
% This is an important number and will determine how many points in the sine wave are sampled (& whether the points are wrapped around).
% Note that 'PeyeFrames' is determined by the stimulus duration and may not be the same as 'FramesFullCycle'
PerEyeFrames = round(PeyeFR*TrialDuration); % also round to integer number of frames

% Contrast:
pdc = 1; % the peak dot contrast

% define raised cosine temporal window.
% the value of 'cont' can be fed in for each frame to determine the contrast. It can be multiplied by
% whatever the current contrast value is if contrast thresholds are being measured
cont = pdc.*ones(1,PerEyeFrames);                                        % initialize to peak dot contrast (defined above)
win_length = round(PerEyeFrames/4);                                      % define window length
cont(1:win_length) = pdc.*0.5.*(1-cos(pi*(1:win_length)./win_length));   % ramp up
cont(end:-1:end-win_length+1) = cont(1:win_length);                      % ... and down
cont = [0 cont 0 0];                                                     % Add a few extra zeros in there

% Define the dot texture, a square-shaped sheet of dots.
% Make the texture the same size as the height of the screen
% (well really we are using many little dot textures, but 'imsize' defines the region they inhabit)
% Note, this should really be the same for both CD and IOVD...
imsize = screenRect(4);

% specify dot size:
dot_sigma_in_degrees = 0.075;            % 0.064; %size of SD of the dot profile in degs/vis angle
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

% Set aside some more of the stimulus parameters:
parameters.TrialDurationInSec = TrialDuration;
%parameters.PeyeFR = PeyeFR;
parameters.PerEyeFrames = PerEyeFrames;
parameters.PeakDotContrast = pdc;
parameters.DotSigmaInDeg = dot_sigma_in_degrees;
parameters.OuterStimRadInDeg = outrad/PPD;
parameters.InnerStimRadInDeg = inrad/PPD;


%%%%-------------------------------%%%%
%       Set up fixation
%%%%-------------------------------%%%%

% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

% Make the fixation lock rings:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.5;                % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;        % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/4; % 1/4 of a degree thick
ringWidthOuter = ringRadiusOuter - PPD/3; % 1/3 of a degree thick

% Make the rings. Both are in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
% Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.
ring_alphaOuter = ((r>ringWidthOuter+1) & (r<ringRadiusOuter-1));

%%%%--------------------------------------------%%%%
%      Set up the staircase for the Psi procedure
%%%%--------------------------------------------%%%%

% Use Psi functions
num_psi = 2; % 2 interleaved staircases

xLevels = linspace( .01, 1, 100 ); % the signal coherence/contrast levels, steps of 1%
aLevels = xLevels;
bLevels = linspace( 0.5, 10, 20 );

for ii = 1:num_psi
    stair{ii} = Psi( 'xLevels', xLevels, ...
        'alphaLevels', aLevels, ...
        'betaLevels', bLevels, ...
        'psychFunction', 'quick'); % default lapse rate for Quick function is 0.04; giving a threshold of 80.3% correct
    stair{ii} = PsiStep( stair{ii} );
end

%sigCoh = stair{1}.xLevels( stair{1}.currentXIndex ); %the starting coherence/contrast level

% ----------------------------------
% Determine order of the psi staircases:
num_stim = num_trials*num_psi;                       % total number of trials over whole run.
ResponseArray = zeros(num_stim,6);                   % store the data in here: the psi/trial, the stim level,
% the 1st interval (0=noise, 1=signal),
% the response (1 or 0 for correct/wrong), the alphaHat value,
% And the alphaHat error value.
psi_order = Shuffle(ceil([1:num_stim]./num_trials)); % randomise presentation of the 2 staircases
%psi_order = ceil([1:num_stim]./num_trials);         % non-randomised
Responses = zeros(num_stim,1);                       % all responses coded as zeros (ie wrong) by default

% Determine the order for each trial: is the signal in the first or 2nd interval?
% In the order, 1st:2nd interval, 0 = noise, 1 = signal
trial_order = round(rand(num_stim,1));
trial_order = [trial_order, ~trial_order];

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
    %PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    PsychImaging('AddTask', 'General', 'FloatingPoint32Bit'); % Just force it to use 32 bit precision (seems to work)
    % Initialise the M16 mode for 16 bit luminance output precision
    PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay'); 

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
    
    % Now also prepare the overlay. This is to allow drawing of the blue lines whilst in M16 mode.
    overlay = PsychImaging('GetOverlayWindow', win);
    
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
            
            %Datapixx('EnableVideoLcd3D60Hz');
            Datapixx('DisableVideoLcd3D60Hz'); %=> weirdly, this seems to cause less crosstalk (see DBaker)
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
    
    % Set up the overlay CLUT giving us blue only:
    OverlayClut = repmat([0 0 1], [256,1]);   % By default, all overlays are transparent
    % We need to make the top 2 rows black (0,0,0). This seems to fix the clut and give us 
    % black for row 1, blue for row 256.
    % It seems there is an 'off by 1' error in there somewhere.
    % From the demo DatapixxM16Demo:
    % "% ! ON WINDOWS, DrawFormattedText scales the color by 255/256, therefore
    % the color is off by 1 for the upper half of the CLUT"
    OverlayClut(1:2,:) = 0; %Make top 2 rows black.

    % load the overlay CLUT: it is always loaded into the overlay when M16 mode is on:
    Datapixx('SetVideoClut', OverlayClut);
    
    HideCursor;
    % raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
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
    
    % Generate the Inner ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaInner;
    fixationRingTextureInner = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Generate the Outer ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaOuter;
    fixationRingTextureOuter = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Preallocate array with destination rectangles:
    % This also defines initial dot locations
    % for the very first drawn stimulus frame:
    texrect = Screen('Rect', Dots(1));
    
    % If motion coherence < 1.0, we will need to adjust the dot x positions of the newly-determined noise dots so they are in the centre of the screen.
    % If we don't do this they will be shifted off to the left of screen. This is because imsize is smaller than the full width (ie x distance) of the screen
    AdjustXBy = (screenRect(3) - screenRect(4))/2; % shift X dot positions by this amount to centre them, since image size <screenRect(3)
    
    % Pick and load a stimulus first-off, just to set aside memory for the variables but also to
    % store the cue stimulus parameters now (rather than on each trial).
    % trial_order must be entered as a scalar. Here we enter -1 as a flag to store the parameters & indicate it's not a trial yet
    C = 1;
    Trial = 0;
    [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
        (-1, C, Trial, CohOrCont, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize, parameters);
    
    % Display the welcome screen and wait for the user to begin.
    
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 24);
    if useHardwareStereo
        Screen('SelectStereoDrawBuffer', win, 0);   % flag of 0= left eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 1 '' if the smooth 3D motion is in the 1st presentation interval,', ...
            '\n \nor press '' 2 '' if the smooth 3D motion is in the 2nd presentation interval.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        
        Screen('SelectStereoDrawBuffer', win, 1);   % flag of 1= right eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 1 '' if the smooth 3D motion is in the 1st presentation interval,', ...
            '\n \nor press '' 2 '' if the smooth 3D motion is in the 2nd presentation interval.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); % , [], [], [], 1);
    else
        
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 1 '' if the smooth 3D motion is in the 1st presentation interval,', ...
            '\n \nor press '' 2 '' if the smooth 3D motion is in the 2nd presentation interval.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); %, [], [], [], 1);
    end
    
    WaitSecs(0.2)
    KbCheck(); % take a quick KbCheck to load it now & flush any stored events
    
    % Wait for user response to continue...
    ButtonPressed = 0;
    while ~ButtonPressed
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
    end
    
    % Blank the screen and wait 1 secs before beginning.
    Screen('Flip', win);
    WaitSecs(1);
    missedFrames = 0;
    
    %%%%-----------------------------------------------------------%%%%
    %               Perform tasks across trials
    %%%%-----------------------------------------------------------%%%%
    
    for Trial = 1:num_stim  % loop across trials
        
        %---------------------------------------------------------%
        %  Determine current stim based on current Psi staircase
        %---------------------------------------------------------%
        
        % Select which psi to present (random order determined above)...
        trident = psi_order(Trial)
        Trial                       % print out trial number
        stair{trident}.alphaHat     % print out current threshold estimate (NaN on first trial)
        
        % Determine the signal coherence/contrast for the next trial (and print out)
        C = stair{trident}.xLevels( stair{trident}.currentXIndex );
        
        % If we're measuring contrast thresholds, modifying contrast is very easy. We'll do that now.
        % Coherence is slightly tricker & depends on the number of dots so that is done later.
        if strcmp(CohOrCont, 'cont')  % If we're measuring contrast thresholds
            TrialContrast = cont * C; % modify contrast value by current staircase setting
        else
            TrialContrast = cont;     % otherwise, leave it at full contrast, as set above.
        end
        
        % Draw the fixation rings/cross for 1 sec to indicate the start of a trial and prepare the subject:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
        Screen('DrawTexture',win, fixationRingTextureInner);
        Screen('DrawTexture',win, fixationRingTextureOuter);
        % Draw the fixation cross:
        Screen('FillRect',win,[0 0 0],fixationCross);
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        Screen('DrawingFinished', win);
        Screen('Flip', win);
        WaitSecs(1);
        
        for Interval = 1:2            % Loop across the 2 intervals.
            
            %---------------------------------------------------------%
            %           Load pre-generated dot positions
            %---------------------------------------------------------%
            
            % Determine and load the sequence parameters for the current interval of the current trial.
            % Note that trial_order is entered as a scalar to define the current interval where 0=noise, 1=signal
            % SetUpDotPositions is a nested sub-function defined below
            tic
            %C=1 % for 100% coh every time
            [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
                (trial_order(Trial, Interval), C, Trial, CohOrCont, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize, parameters);
            toc
            
            %%%%---------------------------------------------------------------%%%%
            %                     Present the stimuli!
            %%%%---------------------------------------------------------------%%%%
            
            if ~ismac
                jheapcl; %clear the java heap space.
            end
            
            % Draw the rings again for 1 frame so that there is no 'blank'
            % frame presented in between, & hopefully stop the 'flickering' of the goggles between intervals
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);
            % Draw the fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            Screen('DrawingFinished', win);
            vbl = Screen('Flip',win); %sync vbl to start time
            
            f0 = 1; % This value increments across PER EYE frames (ie. every 2nd video frame): for contrast
            % Also select a random start point within the cycle to begin stim presentation at (motion wraps for the trial duration)
            % It will also increment on a PER EYE frame basis, but with a random starting point
            % Even though f0 and F could be the same thing, F starts at a random time point in the motion cycle, and indexes the motion
            % no matter how many frames are sampled (using mod(F,FramesFullCycle).
            % f0 only indexes the dot contrast values (for the temporal ramp) so must increment 1:PerEyeFrames
            F = round(rand*PerEyeFrames);
            
            while f0 <= PerEyeFrames % loop across stimulus frames (per eye)
                
                % Select left-eye image buffer for drawing:
                if useHardwareStereo
                    Screen('SelectStereoDrawBuffer', win, 0);
                end
                %%%%------------------------------------------------%%%%
                %               Draw left eye stimulus:
                %%%%------------------------------------------------%%%%
                
                % Draw dots:
                Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); % turn on alpha blending for lin superposition: see GarboriumDemo.m
                Screen('DrawTextures',win,DotsIdxL(:,mod(F,FramesFullCycle)+1), ...             % texture index
                    [], ...
                    dstRectsL(:,:,mod(F,FramesFullCycle)+1), ...    % select appropriate dot positions for current frame
                    [],[],TrialContrast(f0))                        % final argument is contrast
                
                % Superimpose the annulus:
                if DrawAnnulus
                    Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                    Screen('DrawTexture',win,annulus);
                    Screen('Blendfunction', win, GL_ONE, GL_ZERO);
                end
                
                % Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
                if DrawRings
                    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
                    Screen('DrawTexture',win, fixationRingTextureInner);
                    Screen('DrawTexture',win, fixationRingTextureOuter);
                    % Draw the black fixation cross:
                    Screen('FillRect',win,[0 0 0],fixationCross);
                end
                
                Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
                Screen('BlendFunction', win, GL_ONE, GL_ZERO);
                
                % Draw blue lines in the overlay channel :
                if UsingVP
                    Screen('FillRect', overlay, 0); % always need to make the full screen overlay transparent first
                    Screen('FillRect', overlay, 256, blueRectLeftOn); % Should be blue
                    Screen('FillRect', overlay, 1, blueRectLeftOff); % should be black
                end
                
                % Select right-eye image buffer for drawing:
                if useHardwareStereo
                    Screen('SelectStereoDrawBuffer', win, 1);
                else
                    Screen('DrawingFinished', win);
                    [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); % update display on next refresh (& provide deadline)
                    
                    if missed > 0
                        missedFrames = missedFrames + 1;
                    end
                end
                
                %%%%------------------------------------------------%%%%
                %               Draw right eye stimulus:
                %%%%------------------------------------------------%%%%
                
                % Draw dots:
                Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); % turn on alpha blending for lin superposition: see GarboriumDemo.m
                Screen('DrawTextures',win,DotsIdxR(:,mod(F,FramesFullCycle)+1), ...             % texture index
                    [], ...
                    dstRectsR(:,:,mod(F,FramesFullCycle)+1), ...    % select appropriate dot positions for current frame
                    [],[], TrialContrast(f0))                        % final argument is contrast
                
                % Superimpose the annulus:
                if DrawAnnulus
                    Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                    Screen('DrawTexture',win,annulus);
                    Screen('Blendfunction', win, GL_ONE, GL_ZERO);
                end
                
                % Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
                if DrawRings
                    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
                    Screen('DrawTexture',win, fixationRingTextureInner);
                    Screen('DrawTexture',win, fixationRingTextureOuter);
                    %Draw the fixation cross:
                    Screen('FillRect',win,[0 0 0],fixationCross);
                end
                
                Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
                Screen('BlendFunction', win, GL_ONE, GL_ZERO);
                
                % Draw blue lines in the overlay channel:
                if UsingVP
                    Screen('FillRect', overlay, 0); % always need to make the full screen overlay transparent first
                    Screen('FillRect', overlay, 256, blueRectRightOn); % should be blue
                    Screen('FillRect', overlay, 1, blueRectRightOff);  % should be black 
                end
                  
                Screen('DrawingFinished', win);
                
                [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
                
                % keep record of any missed frames:
                if missed > 0
                    missedFrames = missedFrames + 1;
                end
                
                % Increment frame counters on a PER EYE basis:
                f0 = f0+1; % For contrast index
                F = F+1;   % For stimulus position index
                
            end % end of single trial stimulus loop
            
            % Keep fixation rings/lock on screen between the intervals for 1 s
            if Interval == 1
                Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
                Screen('DrawTexture',win, fixationRingTextureInner);
                Screen('DrawTexture',win, fixationRingTextureOuter);
                %Draw the fixation cross:
                Screen('FillRect',win,[0 0 0],fixationCross);
                Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
                Screen('BlendFunction', win, GL_ONE, GL_ZERO);
                Screen('DrawingFinished', win);
                Screen('Flip', win);
                WaitSecs(0.5);
            end
            
        end % End of loop across the 2 intervals.
        
        % Blank the screen while we await a response
        Screen('Flip', win);
        
        %%%%------------------------------------------------%%%%
        %               Check for response
        %%%%------------------------------------------------%%%%
        
        % Wait for a valid button press. If that response is correct, then indicate the response
        % as such. If not correct, response is coded as incorrect by default.
        % The loop will continue to wait until a valid button press is made
        % Responses can be made on the keyboard, or, if you're on the mac,
        % on the gamepad.
        
        % Wait for response if we haven't had one already ...
        validResp = false;
        % Continues to check for responses until  vald one is made.
        while ~validResp
            
            % Check for keyboard response
            [keyIsDown, ~, keyCode] = KbCheck();
            
            % Only use gamepad if you're on the mac:
            if ismac
                if ~keyIsDown % And only check if a keyboard response hasn't been made
                    [keyIsDown, ~, keyCode] = KbCheck(GamepadIndex); % Check for gamepad response
                end
            end
            
            if keyIsDown %If they've pressed a button on keyboard or gamepad,
                % determine if it is valid and proceed
                
                %If they respond 'interval 1' and the signal is in interval 1: **CORRECT**
                if any(keyCode(Interval1)) ...
                        && (trial_order(Trial, 1))      % Signal in interval 1
                    Responses(Trial) = 1;               % code response as correct
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(2000,0.4,0.05); Beeper(2000,0.4,0.05);
                    end
                    
                    % If they respond 'interval 2' and the signal is in interval 2: **CORRECT**
                elseif any(keyCode(Interval2)) ...
                        && ~(trial_order(Trial, 1))     % Signal not in interval 1
                    Responses(Trial) = 1;               % code response as correct
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(2000,0.4,0.05); Beeper(2000,0.4,0.05);
                    end
                    
                    % If they respond 'interval 1' and the signal is in interval 2: **INCORRECT**
                elseif any(keyCode(Interval1)) ...
                        && ~(trial_order(Trial, 1))     % Signal not in interval 1
                    % No need to recode response (coded as incorrect by default).
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(400,0.4,0.15); Beeper(400,0.4,0.5);
                    end
                    
                    % If they respond 'interval 2' and the signal is in interval 1: **INCORRECT**
                elseif any(keyCode(Interval2)) ...
                        && (trial_order(Trial, 1))      % Signal in interval 1
                    % No need to recode response (coded as incorrect by default).
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(400,0.4,0.15); Beeper(400,0.4,0.5);
                    end
                    
                elseif keyCode(RespQuit)
                    % break out of program if 'q' is pressed
                    ForcedQuit = true;
                    ExitGracefully(UsingVP, ForcedQuit)
                    
                end
            end
        end
        
        % write stimulus & response details to ResponseArray before updating staircase...
        % In the order the psi/trial, the stim level,
        % the 1st interval (0=noise, 1=signal) designation,
        % the response (1 or 0 for correct/wrong), the alphaHat value
        % and the error associated with that alphaHat value
        ResponseArray(Trial,:) = [trident, C, trial_order(Trial, 1), Responses(Trial), stair{trident}.alphaHat, stair{trident}.alphaHatSD];
        
        %---------------------------------------------------------%
        %                   Update the staircase:
        %---------------------------------------------------------%
        
        % update the staircase for the next trial, based on the subject's response
        
        % *** NOTE *** the staircases are updated at the end of the final trial AFTER the value
        % has been placed in 'ResponseArray', meaning that the very last value stored in the staircase struct arrays
        % (StairCase{trident}.alphaHat is the *TRUE* final threshold estimate, not the final value in the ResponseArray matrix
        
        
        stair{trident} = PsiUpdate( stair{trident}, Responses(Trial));
        stair{trident} = PsiStep( stair{trident} );
        
        
    end %end of loop across trials
    
    %%%%-----------------------------------------------------------%%%%
    %                       Save the results
    %%%%-----------------------------------------------------------%%%%
    
    % Store the data and some info about the run
    results.ResponseArray = ResponseArray;
    results.SubjectCode = SubjCode;
    results.RunNumber = runN;
    results.StairCase = stair;
    % For quick reference, set aside the frequency & amplitude parameters
    parameters.MotionFrequencyInHz = frequency;
    parameters.DisparityInArcMin = amplitude;
    
    % Save the data and parameters information!
    save(fName, 'parameters', 'results')
    
    missedFrames % print out number of missed frames during experiment
    
    % Finally, plot the staircases.
    figure
    if strcmp(parameters.SignalVaried, 'coherence')
        Ylab = {'Coherence (%)'};
    elseif strcmp(parameters.SignalVaried, 'contrast')
        Ylab = {'Contrast (%)'};
    end
    Titles = {'Psi 1', 'Psi 2'};
    
    % All cont/coherence and alpha hat values *100 to convert to percentage
    % Remember, ResponseArray is in the order
    % the psi/trial, the stim level,
    % the 1st interval (0=noise, 1=signal),
    % the response (1 or 0 for correct/wrong), the alphaHat value,
    % and the alphaHat SDs.
    for X = 1:max(results.ResponseArray(:,1)) % equivalent to num_psi
        
        ix = (results.ResponseArray(:,1)==X);                                     % index for values of the current staircase
        ixCorr = (results.ResponseArray(:,1)==X & results.ResponseArray(:,4)==1); % index for correctly-answered trials for current staircase
        subplot(1,2,X)
        % Plot incorrect responses in white first:
        plot(results.ResponseArray(ix,2) .* 100, 'ko-', 'MarkerFaceColor', 'w', 'MarkerSize', 10)
        set(gca, 'box', 'off', 'TickDir', 'out', 'FontSize', 14)
        % Now determine correct responses and plot on top:
        corrResp = nan(length(results.ResponseArray),1);
        corrResp(ixCorr) = results.ResponseArray(ixCorr,2);
        hold on
        plot(corrResp(ix) .* 100, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
        % Now plot alpha hat values (threshold estimates) for each trial in current staircase:
        %plot(results.ResponseArray(ix,5) .* 100, 'r.-')
        % Also plot the Alpha Hat standard deviations:
        errorbar(results.ResponseArray(ix,5) .* 100, results.ResponseArray(ix,6) .* 100, 'r.-')
        xlabel('Trial')
        ylabel(Ylab)
        title(Titles{X})
        ylim([0 100])
        
    end
    suptitle(['Observer ', results.SubjectCode, ', ', ...
        parameters.Cue, ' Psi staircases, ', num2str(parameters.MotionFrequencyInHz), ...
        ' Hz, ', num2str(amplitude), ' arcmin, Run ', num2str(runN)])
    legend('Incorrect', 'Correct', 'alpha hat')
    
    % Close everything down & exit:
    ExitGracefully (UsingVP, ForcedQuit)
    
catch MException
    
    rethrow (MException)
    psychrethrow(psychlasterror)
    ExitGracefully (UsingVP, ForcedQuit)
    error('Error!')
    
end % End of try/catch statement
end % End of main function

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
function rect = fixation_cross(width,height,centrex,centrey)
% $Id: fixation_cross.m 4 2007-06-23 11:13:19Z damienm $
%Make a small fixation cross to display on screen.
%Width and height are in pixels.
%centrex & centrey give the x & y screen coordinates where you want the cross centred.

rect = zeros(4,2);

width = width/2;
height = height/2;

rect(1,1) = -height;
rect(2,1) = width;
rect(3,1) = height;
rect(4,1) = -width;

rect(1,2) = -width;
rect(2,2) = height;
rect(3,2) = width;
rect(4,2) = -height;


rect(1:2:4,:) = rect(1:2:4,:) + centrex;
rect(2:2:4,:) = rect(2:2:4,:) + centrey;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
    (trial_order, C, Trial, CohOrCont, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize, parameters)
% This is where the pre-generated dot positions are loaded and the disparity added.
% (see GenerateFULLCUEdots.m)
% These are picked at random (with replacement) from those already in the folder.
% Note that trial_order is entered as a scalar to define the order of the current interval, where 0 = noise, 1=signal

disp(trial_order)

% Indicate how many existing sets of dot stimuli are in the subfolder 'stimuli' for loading later.
StimInFolder = 30;

% Define the max. change in disparity (horiz shift of dots):
% the numerator is in arcmin (but converted to pixels) akin to the amplitude of the sine wave
FULL_disparity = amplitude/60 * PPD;

% Convert frequency into a string so we can load the appropriate stimulus:
FULLstrfreq = num2str(frequency);
FULLstrfreq(FULLstrfreq == '.') = '_'; % replace decimal point with underscore

% Draw a random value to load the stimuli with:
RandDraw = ceil(rand*StimInFolder);

%%%%-------------------------------------------------------------------------%%%%
%                     Shift dot trajectories: add disparity
%%%%-------------------------------------------------------------------------%%%%

% *** ---------------------------------- ***
%               Set up FULL cue
% *** ---------------------------------- ***

FileNameFULL = fullfile('stimuli',['FULLCUE_dots_', FULLstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
FULL_stim = load(FileNameFULL); % includes parameters file and dot positions matrices
% Just cut the parameters back to a smaller struct:
FULLparams = FULL_stim.DotParams;

% *** set up dot coherence ***
% Now we'll determine dot coherence according to the staircase
if strcmp(CohOrCont, 'coh')                % If we're measuring coherence thresholds
    NumDots = round(FULLparams.NumDots * C); % modify number of coherent dots by value by current staircase setting
    
    % Now that we know the number of coherent dots, set up some random indices that will determine which dots are 'randomised'
    % this is necessary to ensure that dots of both colours are made noise dots equally.
    % We also want these indices to be consistent across frames in a trial, so they are set here.
    % Remember that we want to index the NOISE dots to change (ie=1), not the signal dots,
    % so it may seem a bit counterintuitive that noise =1, coherent =0
    NoiseDotIndices = Shuffle([true(FULLparams.NumDots-NumDots,1); false(NumDots,1)]); %The 1s in here should sum to FULLparams.NumDots-NumDots: the number of noise dots.
    
else
    NumDots = FULLparams.NumDots;            % otherwise, all dots are coherent.
end

% Define the number of frames in a full cycle:
FramesFullCycle = FULLparams.FramesFullCycle;

% Assign the pre-generated dot position matrices.
% First, we duplicate the left eye dot positions to make the right eye dot positions.
% At this stage all dots are identical for the two eyes (disparity has not been added yet).
dot_posL = FULL_stim.dot_posL;
dot_posR = FULL_stim.dot_posL;

% set aside the disparity
FULLparams.disparity = FULL_disparity;
% For reference also store the disparity in Arcmin
FULLparams.DisparityInArcmin = (FULL_disparity/PPD)*60;

% set up a few extra parameters of the sine waves for each cue:
% the period, in sec
FULLparams.period = 1/FULLparams.frequency;

% the angular frequency of the sine waves:
FULLparams.angFreq = 2 * pi * FULLparams.frequency;

% Length of the sine wave:
% The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
% end of the cycle and the beginning of the next cycle.
% So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
FULLparams.t = linspace(FULLparams.period/FULLparams.FramesFullCycle, FULLparams.period, FULLparams.FramesFullCycle);

% Now make one full cycle of the sine wave, no matter the frequency
% Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
% Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
FULLparams.SineWave =  FULLparams.disparity * sin(FULLparams.angFreq * FULLparams.t); %disparity * sin(angFreq * t)

% assign dstRects: destination rect matrices
dstRectsL = zeros(4, FULLparams.NumDots,FULLparams.FramesFullCycle);
dstRectsR = zeros(4, FULLparams.NumDots,FULLparams.FramesFullCycle);

% *** index the sine wave, depending on whether it's signal or noise trial ***
% Now we need to set up an index of the sine wave so each dot has disparity added at a different phase for the NOISE intervals.
% We must index each point of the sine wave, which is the same as the frames in a full cycle.
% Each dot should have the disparity added smoothly, although it will be at a random phase for each dot,
% meaning the overall average disparity will be the same
% (this bit added 2/11/16)

% Remember that 'trial_order' can have values of -1, 0 or 1.
if trial_order == 0 % For a noise interval:
    for nd = 1:FULLparams.NumDots
        SineWvIdx(nd, :) = circshift(1:FULLparams.FramesFullCycle, [0 ceil(FULLparams.FramesFullCycle*rand)]); %random phase for each dot
    end
    % Now we have a row of indices for every single dot. These indices are for the disparity sine wave, each at a random phase.
else  % For a signal interval
    SineWvIdx = repmat(1:FULLparams.FramesFullCycle, FULLparams.NumDots, 1); %simply replicate for each dot.
end

% Shift dot trajectory: add disparity
for fr = 1:FULLparams.FramesFullCycle
    
    % determine dot coordinates: remember, y position does not change: horiz disparity only
    % Update dot position according to sinsoidal trajectory on each (per-eye) frame
    % Left eye: -sin
    dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ... % the size of the dot texture
        dot_posL(:,1,fr) - FULLparams.SineWave(SineWvIdx(:,fr))', ...  % the x positions
        dot_posL(:,2,fr))';                            % the y positions
    % Right eye: +sin
    dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
        dot_posR(:,1,fr) + FULLparams.SineWave(SineWvIdx(:,fr))', ...
        dot_posR(:,2,fr))';
    
    % Now, randomise the dot positions of the noise dots when coherence <1.0
    if (C < 1) && (strcmp(CohOrCont, 'coh'))
        %Left eye:
        dstRectsL(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
            rand(FULLparams.NumDots-NumDots,1) * imsize + AdjustXBy, ... % random x positions
            rand(FULLparams.NumDots-NumDots,1) * imsize)';               % random y positions
        %right eye:
        dstRectsR(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
            rand(FULLparams.NumDots-NumDots,1) * imsize + AdjustXBy, ...
            rand(FULLparams.NumDots-NumDots,1) * imsize)';
    end
    
end

%Left eye:
DotsIdxL = [repmat(Dots(1),1,FULLparams.NumDots/2); repmat(Dots(2),1,FULLparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
DotsIdxL = repmat(reshape(DotsIdxL,1,FULLparams.NumDots)', 1, FULLparams.FramesFullCycle); % Replicate for each frame
%Right eye:
DotsIdxR = [repmat(Dots(1),1,FULLparams.NumDots/2); repmat(Dots(2),1,FULLparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
DotsIdxR = repmat(reshape(DotsIdxR,1,FULLparams.NumDots)', 1, FULLparams.FramesFullCycle); % Replicate for each frame

% Set aside details of the stimulus parameters.
% Most of these are already stored in the pre-generated stimulus files,
% but it is probably worthwhile to have them here again in a single output file
if trial_order == 1
    parameters.StimulusNumberEachTrial{Trial,1} = FileNameFULL; % Store the filename of the randomly-drawn stimulus for SIGNAL
elseif trial_order == 0
    parameters.StimulusNumberEachTrial{Trial,2} = FileNameFULL; % Store the filename of the randomly-drawn stimulus for NOISE
elseif trial_order == -1
    % -1 indicates that it is during set up, so set aside parameters
    parameters.CDparameters = FULLparams;                       % Store the CD parameters once, during set-up
end


end