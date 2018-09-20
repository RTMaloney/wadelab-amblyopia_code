function GenerateFULLCUEdots (TrialsToSave, dot_dens_per_deg2, SpatialJitter, frequency)

% This function will set up and save 'TrailsToSave' sets of random dot positions for a random dot
% stereogram (RDS).
% The RDS contains both Changing disparities over time (CD) and interocular velocity difference (IOVD) cues to motion-in-depth.
% These can then be loaded and used in an experiment as needed.
% The random dot x and y positions are separated by the minimum separation 'SpatialJitter'
% ie pseudo-random positioning to prevent dot overlap/clustering
% Dots will have the same lifetime as the IOVD dots (ie 3 frames per eye, = 50 ms).
%
% Dots for each separate trial are saved as individual files. Inside each is a struct called DotParams that gives the parameters
% that went into determining the stimuli.
% It also contains a 3D matrix called dot_posL. This contains the initial dot positions for every frame in a full cycle.
% This is for a single eye: because the 2 eyes' stimuli will be identical, we need only make 1 version. Later, prior to presentation,
% it can be duplicated and lateral shifts of the dots in different directions are added to each eye, in opposite directions.
% Inside dot_posL are 3D matrices of dimensions: number of dots rows * 2 columns (x&y positions) * number of frames
% The 3rd dimension, number of frames, is determined by the frequency of the sine wave of the motion. A full cycle is generated (even if it isn't eventually used)
% of as many (per-eye) frames as necessary. Thus a lower frequency will require more frames (longer period), more space in the matrix and a longer time to generate.
% Making a full cycle allows the stimulus to be wrapped & re-presented over and over if necessary (or presented at a random phase).
%
% Drawing the dots may take some time (depending on number of trials, dot density, frequency and spatial jitter),
% so please be patient and consider making all stimuli well in advance of actual testing.
% Parameters:
%   TrialsToSave: the number of unique sets of dot stimuli you want (ie number of unique trials).
%   dot_dens_per_deg2: Dot density should be in dots per square degree
%   SpatialJitter: minimum separation of dots, should be in deg of visual angle
%   frequency: the temp frequency of the motion in depth, in Hz.
%
% R Maloney, September 2016

%Determine some of the parameters needed for the stimuli:
PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)
DotLifetimeFrames = 3; %Length of dot lifetime in frames (remember, this is the per-eye frame rate)

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
RefreshRate = Screen('NominalFrameRate', WhichScreen);

% Define the dot texture, a square-shaped sheet of dots.
%Make the texture the same size as the height of the screen
imsize = screenRect(4);

%compute number of dots in the total available area
num_dots = round(dot_dens_per_deg2 * (imsize/PPD)^2); %this many dots in the full dot field

%Just check whether num_dots is odd or even: (important for when contrast polarity is assigned)
%If odd, -1 to make it even
if mod(num_dots,2) %if odd, mod = 1
    num_dots = num_dots-1;
end

%Define the minimum spacing between the PEAKS of the dot positions (in pixels):
DotParams.SpatialJitterDeg = SpatialJitter; % Just store the spatial jitter in degrees value before converting to pixels.
SpatialJitter = round(SpatialJitter * PPD); %+ dotsize/2;
%NOTE: add dot radius (ie dotsize/2) ensures the EDGES of the dots are separated by the minimum, but there may not be space in the matrix!

%If the frequency of the sine wave is specified in Hz (ie cycles/sec)
%then other units of time must also be in SECONDS!!!
period = 1/frequency; %in sec
%the effective frame rate PER EYE: since each eye is stimulated on successive video frames
%Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
PeyeFR = RefreshRate/2; %Per eye f.r. is always half the absolute f.r.

%The no. of frames for a full cycle:
%Remember that each eye samples the motion vector at the same point, so this is in terms of per-eye frame rate
%(otherwise the different eyes would sample successive points on the trajectory/sine wave: one eye would lead the other)
FrmsFullCyc = round(PeyeFR*period); %must have integer no. of frames

%Store some of the parameters of the stimuli along with each stimulus set:
DotParams.FramesFullCycle = FrmsFullCyc; %No. of frames in a full cycle of the motion in depth.
DotParams.frequency = frequency; %in Hz, freq of the sinusoidal change in disparity
DotParams.PerEyeFR = PeyeFR; %in Hz
DotParams.SpatialJitterPix = SpatialJitter; %in pixels
DotParams.NumDots = num_dots; %number of dots in matrix
DotParams.ImageSize = imsize; %size of image matrix
DotParams.TrialsInBatch = TrialsToSave; %total number of trials in this batch
DotParams.DotDensityDeg2 = dot_dens_per_deg2; %the dot density in dots/deg^2
DotParams.DotLifetimeMS = (1000/PeyeFR)*DotLifetimeFrames;
DotParams.PixelsPerDeg = PPD; %In case it's run on different resolution

%Set up dot indices to determine dot lifetime:
%We need to re-randomise 1/3 of all dots on every frame
%This ensures all dots have a lifetime of 3 frames (or 50 ms, with a per-eye frame rate of 60 Hz).
%So we set up the indices for (roughly) every 1/3 of the dots, leaving no dot unturned
%This is the same for every trial so only needs to be set once.
DotThirdIndices = round(linspace(1,num_dots,4));
DotThirdIndices = [DotThirdIndices(1:3)', DotThirdIndices(2:4)'];
DotThirdIndices(2:3,1) = DotThirdIndices(2:3,1)+1;

for T = 1:TrialsToSave %loop across all trials
    
    disp(['Generating trial ' num2str(T) ' of ' num2str(TrialsToSave)])
    tic
    
    %%%%-------------------------------------------------------------------------%%%%
    %           Determine the dot positions for the very first frame
    %%%%-------------------------------------------------------------------------%%%%
    
    CurrDot=1;
    dot_pos = zeros(num_dots,2); %assign dot location matrix
    while CurrDot <= num_dots
        
        if CurrDot == 1 %set random coordinates for very first dot
            dot_pos(CurrDot,:) = imsize.*rand(1,2);
            CurrDot = CurrDot+1;
        else
            %set the next dot's random position
            dot_pos(CurrDot,:) = imsize.*rand(1,2);
            %find the smallest distance (in pixels) between the current dot and any other existing dot
            idx = 1:CurrDot-1; %index each existing dot except the current one
            d = min((dot_pos(idx,1)-dot_pos(CurrDot,1)).^2 + (dot_pos(idx,2)-dot_pos(CurrDot,2)).^2);
            d = sqrt(d);
            
            %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
            %This will continue until (at least) the minimum distance is met
            if d < SpatialJitter
                dot_pos(CurrDot,:) = imsize.*rand(1,2);
            else %if that minimum distance is met, move on to the next dot
                CurrDot = CurrDot+1;
            end
        end
    end
    
    %%%%-------------------------------------------------------------------------%%%%
    %       replicate dot_pos by the number of frames needed (to allocate memory)
    %%%%-------------------------------------------------------------------------%%%%
    
    % Only need to replicate once, because at this stage, both eyes' images are identical
    % (disparity/lateral shift has not been added)
    
    dot_pos = repmat(dot_pos,1,1,FrmsFullCyc); %for the left
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Determine dot lifetime
    %%%%-------------------------------------------------------------------------%%%%
    
    %Shift the position of a random 1/3 of dots. This effectively ends the lifetime of the dot and begins it anew somewhere else
    
    %set dot positions for each frame in a full cycle
    for m = 1:FrmsFullCyc
        
        %Make the dots the same position as the previous frame.
        %1/3 will then be moved.
        if m~=1
            dot_pos(:,:,m) = dot_pos(:,:,m-1);
        end
        
        Curr3rd = mod(m,3)+1; %tells us whether this is a 1st, 2nd or 3rd frame, & determines which 3rd of dots to change
        CurrRows = DotThirdIndices(Curr3rd,:);
        dot_pos(CurrRows(1):CurrRows(2),:,m) = nan; %make the third we want to change NaNs to put them out of the calculations
        CurrDot = CurrRows(1);
        
        while CurrDot <= CurrRows(2) %go through all dots in this 3rd
            
            %set the next dot's random position
            dot_pos(CurrDot,:,m) = imsize.*rand(1,2);
            %find the smallest distance (in pixels) between the current dot and any other existing dot
            %Index all the existing dots except the current one to do this
            %This means excluding all the dots currently unassigned (defined as NaNs).
            CurrentEmpty = find(isnan(dot_pos(:,1,m)))'; %all rows currently empty
            idx = setdiff(1:num_dots, [CurrDot, CurrentEmpty]);
            
            %Find the smallest distance between the current dot and any other dot
            d = min((dot_pos(idx,1,m)-dot_pos(CurrDot,1,m)).^2 + (dot_pos(idx,2,m)-dot_pos(CurrDot,2,m)).^2);
            d = sqrt(d);
            
            %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
            %This will continue until (at least) the minimum distance is met
            if d < SpatialJitter
                dot_pos(CurrDot,:,m) = imsize.*rand(1,2);
            else %if that minimum distance is met, move on to the next dot
                CurrDot = CurrDot+1;
            end
            
        end %end of loop across the dots in the current 3rd
        
    end %end of loop across all frames in a full cycle
    
    %adjust the x positions so they are in the centre of the screen. Do this for all frames in one go:
    %If we don't do this they will be shifted off to the left of screen
    AdjustXBy = (screenRect(3) - screenRect(4))/2; %shift X dot positions by this amount to centre them, since image size <screenRect(3)
    dot_pos(:,1,:) = dot_pos(:,1,:) + AdjustXBy;
    dot_posL = dot_pos; % Save these positions as the Left eye. It can be duplicated to produce the right eye just before the disparity is added; prior to presentation.
    
    %Save the result:
    DotParams.StimulusNumber = T;   % store the number of the current set of dots (also in file name)
    strfreq = num2str(frequency);   % generate a string of the frequency
    strfreq(strfreq == '.') = '_';  % replace the decimal pt (if any) with underscore so it doesn't mess up the file name
    fname = fullfile('stimuli', ['FULLCUE_dots_', strfreq, '_Hz_stim_', num2str(T), '.mat']);    save(fname, 'dot_posL', 'DotParams');
    toc
    
end %end of loop across trials/unique stimulus sets

