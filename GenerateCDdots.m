function GenerateCDdots (TrialsToSave, dot_dens_per_deg2, SpatialJitter, frequency)

% This function will set up and save 'TrailsToSave' sets of random dot positions for a dynamic random dot 
% stimulus for isolating the changing disparity (CD) cue to motion in depth. 
% These can then be loaded and used in an experiment as needed. 
% The random dot x and y positions are separated by the minimum separation 'SpatialJitter'
% ie pseudo-random positioning to prevent dot overlap/clustering
%
% Dots for each separate trial are saved as individual files. Inside each is a struct called DotParams that gives the parameters
% that went into determining the stimuli.
% They are in the form of 3D matrices of dimensions: number of dots rows * 2 columns (x&y positions) * number of frames.
% The 3rd dimension, number of frames, is determined by the frequency of the sine wave of the motion. A full cycle is generated (even if it isn't eventually used)
% of as many (per-eye) frames as necessary. Thus a lower frequency will require more frames (longer period), more space in the matrix and a longer time to generate.
% Making a full cycle allows the stimulus to be wrapped & re-presented over and over if necessary (or presented at a random phase).
%
% Only the left eye is generated because at this stage it is idential to the right eye, so it is simply duplicated.
% They will differ when the disparity is added later (prior to presentation). 
% Drawing the dots may take some time (depending on number of trials, dot density, frequency and spatial jitter),
% so please be patient and consider making all stimuli well in advance of actual testing.
% Parameters:
%   TrialsToSave: the number of unique sets of dot stimuli you want (ie number of unique trials). 
%   dot_dens_per_deg2: Dot density should be in dots per square degree
%   SpatialJitter: minimum separation of dots, should be in deg of visual angle
%   frequency: the temp frequency of the motion in depth, in Hz.
%R Maloney, August 2015

%Determine some of the parameters needed for the stimuli:
PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)

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
imsize = screenRect(4) ;

%compute number of dots in the total available area
num_dots = round(dot_dens_per_deg2 * (imsize/PPD)^2); %this many dots in the full dot field

%Just check whether num_dots is odd or even: (important for when contrast polarity is assigned)
%If odd, -1 to make it even
if mod(num_dots,2) %if odd, mod = 1
    num_dots = num_dots-1;
end

%Define the minimum spacing between the PEAKS of the dot positions (in pixels):
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
 
for T = 1:TrialsToSave
    
    disp(['Generating trial ' num2str(T) ' of ' num2str(TrialsToSave)])
    %assign space to dot location matrices for all frames
    %just do the left eye for now; they are currently identical (disparity hasn't been added)
    dot_posL = zeros(num_dots,2,FrmsFullCyc);
    
    %set dot positions for each frame in a full cycle
    tic
    for m = 1:FrmsFullCyc
        
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
        
        %Set aside the dot position matrix for each frame (left eye only)
        dot_posL(:,:,m) = dot_pos;
        
    end
    
    %adjust the x positions so they are in the centre of the screen. Do this for all frames in one go:
    %If we don't do this they will be shifted off to the left of screen
    AdjustXBy = (screenRect(3) - screenRect(4))/2; %shift X dot positions by this amount to centre them, since image size <screenRect(3)
    dot_posL(:,1,:) = dot_posL(:,1,:) + AdjustXBy;
    
    DotParams.StimulusNumber = T;   % store the number of the current set of dots (also in file name)
    strfreq = num2str(frequency);   % generate a string of the frequency 
    strfreq(strfreq == '.') = '_';  % replace the decimal pt (if any) with underscore so it doesn't mess up the file name
    fname = fullfile('stimuli', ['CD_dots_', strfreq, '_Hz_stim_', num2str(T), '.mat']);
    save(fname, 'dot_posL', 'DotParams');
    toc
end

%All done.
%Bye bye



