function GenerateIOVDdots_2 (TrialsToSave, dot_dens_per_deg2, SpatialJitter, frequency)

% This function will set up and save 'TrailsToSave' sets of random dot positions for a random dot
% stimulus for isolating the interocular velocity difference (IOVD) cue to motion in depth.
% These can then be loaded and used in an experiment as needed.
% The random dot x and y positions are separated by the minimum separation 'SpatialJitter'
% ie pseudo-random positioning to prevent dot overlap/clustering
%
% Dots for each separate trial are saved as individual files. Inside each is a struct called DotParams that gives the parameters
% that went into determining the stimuli.
% Also is a cell array called 'DotPosBinoc'. DotPosBinoc{1} gives the left eye stimulus, DotPosBinoc{2} gives the right eye stimulus.
% Currently these are completely unique sets of dots. Lateral shifts of the dots in different directions will be added later (prior to presentation).
% Inside the two cells of DotPosBinoc are 3D matrices of dimensions: number of dots rows * 2 columns (x&y positions) * number of frames
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
%   Version _2 of this code was developed to remove as many spurious binocular matches from the 2 eyes' images as possible,
% giving us a 'purer' IOVD stimulus.
% Assign dots in L & R eye within different strips for IOVD.
% Dot positions are assigned in adjacent horizontal strips in the L & R eye to prevent overlap
% between the 2 eyes.
% The horizontal strips are 2 dot widths wide.
% Dot positions can be anywhere within this strip for each eye.
% Dot sigma (SD of the Gaussian envelope), and hence dot size/width, is now important for determining dot positions.
% After placing the dot positions within the strips, the smallest distances between dots is determined
% in order to deal with the dots that might fall close to the borders.
% The dots in the R eye are sorted according to their smallest distance to any dot in the L eye.
% Once a L eye dot has been determined to be closest to one R eye dot, it is removed from the
% distance computations that follow for the next R eye dot.
% That way a L eye dot is never assigned more than once, which could be possible if more than one
% L eye dot has their closest distance to a single R eye dot.
% After this sorting, the color (ie Black/white) of those dots must be assigned in alternating sequence:
% that way if a L eye dot is near the border of the strips and is close to a
% R eye dot position, then at least it will be assigned the opposite polarity
% (and hence be less likely to be paired binocularly). But the dot color assignment must happen
% prior to the presentation of the dots (so subsequent to this algorithm that just assigns dot positions).
%
%R Maloney, August 2015
%Version _2, Nov 2015
%
% Updated March 2016 to account for dot colour assignment bug in Right eye sorting. See documentation in MID_Dots_IOVD_3.m
% Now includes a new matrix (MinDist) that appropriately assigns dot colour for both eyes, ensuring they are opposite for the 
% dots that are closest in position across the eyes.
%
% Modified April 2016 to include flag for larger dots depending on whether dots are coloured or achromatic.
% Colour dots with the larger dot sigma are saved under a different file name so that we can run the MID pilot and the MID colour experiments 
% side by side (& hence we can have co-existing stimuli).
% Remember that dot_sigma is relevant here because it determines the width of the strips in the 2 eyes
% NOTE: we should not need to do this for the CD dots, since dot_sigma is modified before they are drawn, and there are no
% strips needed for the CD stimuli. 

% The flag for the colour stimuli:
MakeColour = false;

% Just provide a warning:
if MakeColour
    disp('NOTE: generating colour IOVD stimuli')
end

%Do you want to plot a figure showing the resulting dot positions?
PlotResult = false;

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
SpatialJitter = round(SpatialJitter * PPD); %+ dotsize/2;
%NOTE: add dot radius (ie dotsize/2) ensures the EDGES of the dots are separated by the minimum, but there may not be space in the matrix!

%Determine dot size/width.
%Because dot size determines the width of the L/R eye strips, we now need dot size information.
%This means the dot positions (for IOVD at least) may need to be re-drawn even if all that changes is the dot size.
% *** Work out the dot strips. We will make the dot strips 2 dot widths wide.
% Positions can vary +/-1 dot width from the midpoint of each strip.
% The dot sigma (SD of the Gaussian envelope) must be known here.

% Note that dot_sigma is different for coloured and achromatic stimuli.
if MakeColour
    dot_sigma_in_degrees = 0.2;
else
    dot_sigma_in_degrees = 0.05;
end

% Determine the dot width.
% At full width, half maximum (FWHM), there are approx 2.355 sigma under the Gaussian.
% So dot_sigma * 2.355 will give us the dot width at half the max height.
% Multiply by PPD to get result in pixels.
% Multiply by 2 to give us approx the full dot width
Dot_width = dot_sigma_in_degrees * 2.355 * PPD * 2; %The approx dot width.
StripWidth = 2*Dot_width; %So the strips are (approx) 2 dot widths wide.

%Compute the midpoints of the strips:
midpts = Dot_width:StripWidth:imsize;

%Take left and right eye strip midpoints:
%use a cell array because these vectors may have different lengths
y_pos_midpts{1} = midpts(1:2:end); %Left eye strip midpoints
y_pos_midpts{2} = midpts(2:2:end); %Right eye strip midpoints

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
DotParams.LReyeStripWidthInDeg = StripWidth/PPD; %the L/R eye strip width, in deg

%Set up dot indices to determine dot lifetime:
%We need to re-randomise 1/3 of all dots on every frame
%This ensures all dots have a lifetime of 3 frames (or 50 ms, with a per-eye frame rate of 60 Hz).
%So we set up the indices for (roughly) every 1/3 of the dots, leaving no dot unturned
%This is the same for every trial so only needs to be set once.
DotThirdIndices = round(linspace(1,num_dots,4));
DotThirdIndices = [DotThirdIndices(1:3)', DotThirdIndices(2:4)'];
DotThirdIndices(2:3,1) = DotThirdIndices(2:3,1)+1;

% Dot positions are assigned by first randomly selecting one of the strip midpoints for that eye.
% Once the random midpoint is selected, the random value r is generated and added.
% It is somewhere between -1 dot width & +1 dot widths, given by 'r = (2*rand-1) * Dot_width'
% The dots are still checked & re-randomised if they fail the minimum separation 'spatial jitter'.

for T = 1:TrialsToSave %loop across all trials
    
    disp(['Generating trial ' num2str(T) ' of ' num2str(TrialsToSave)])
    tic
    
    %%%%-------------------------------------------------------------------------%%%%
    %           Determine the dot positions for the very first frame
    %%%%-------------------------------------------------------------------------%%%%
    
    for EachEye = 1:2 %loop across both eyes, L then R
        
        CurrDot=1;
        dot_pos = zeros(num_dots,2); %assign dot location matrix
        y_pos = y_pos_midpts{EachEye}; %select the midpoints, left, then right
        
        while CurrDot <= num_dots
            
            if CurrDot == 1 %set random coordinates for very first dot
                r = (2*rand-1) * Dot_width; % *** gives random val betw -1 dot width & 1 dot width
                %Select a random strip midpoint for this eye, and randomly select 'y' position within the strip
                dot_pos(CurrDot,:) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %was imsize.*rand(1,2);
                CurrDot = CurrDot+1;
            else
                %set the next dot's random position
                r = (2*rand-1) * Dot_width; % *** gives random val betw -1 dot width & 1 dot width
                dot_pos(CurrDot,:) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %imsize.*rand(1,2);
                %find the smallest distance (in pixels) between the current dot and any other existing dot
                idx = 1:CurrDot-1; %index each existing dot except the current one
                d = min((dot_pos(idx,1)-dot_pos(CurrDot,1)).^2 + (dot_pos(idx,2)-dot_pos(CurrDot,2)).^2);
                d = sqrt(d);
                
                %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
                %This will continue until (at least) the minimum distance is met
                if d < SpatialJitter
                    r = (2*rand-1) * Dot_width;
                    dot_pos(CurrDot,:) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %imsize.*rand(1,2);
                else %if that minimum distance is met, move on to the next dot
                    CurrDot = CurrDot+1;
                end
            end
        end
        
        %Set aside dot positions for each eye:
        DotPosBinoc{EachEye} = dot_pos;
        
    end
    
    %%%%-------------------------------------------------------------------------%%%%
    %       replicate dot_pos by the number of frames needed (to allocate memory)
    %%%%-------------------------------------------------------------------------%%%%
    
    %Note that, previously, we were just swapping the x & y positions just generated for the L eye to produce
    %an essentially different set of positions for the R eye (& hence saving computation time), while still satisfying the 'spatial jitter' requirement.
    %We can't do that any more (because dots are now assigned in alternating strips for L & R).
    %The previous method did not guarantee the dots wouldn't overlap in L/R eyes, & they may have coincided along the diagonal anyway
    %DotPosBinoc{2} = repmat([dot_pos(:,2), dot_pos(:,1)],1,1,FrmsFullCyc); %for the right: swap x & y coordinates
    
    DotPosBinoc{1} = repmat(DotPosBinoc{1},1,1,FrmsFullCyc);
    DotPosBinoc{2} = repmat(DotPosBinoc{2},1,1,FrmsFullCyc);
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Determine dot lifetime
    %%%%-------------------------------------------------------------------------%%%%
    
    %Shift the position of a random 1/3 of dots. This effectively ends the lifetime of the dot and begins it anew somewhere else
    %Repeat for the 2 eyes, since both eyes need different (uncorrelated) dots for IOVD
    for TwoEyes = 1:2
        
        dot_pos = DotPosBinoc{TwoEyes}; %re-assign dot_pos, for either the left (1) or right (2) eye
        y_pos = y_pos_midpts{TwoEyes}; %select the midpoints, left, then right
        
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
                r = (2*rand-1) * Dot_width; % *** gives random val betw -1 dot width & 1 dot width
                %Select a random strip midpoint for this eye, and randomly select 'y' position within the strip
                dot_pos(CurrDot,:,m) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r];
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
                    r = (2*rand-1) * Dot_width;
                    dot_pos(CurrDot,:,m) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %imsize.*rand(1,2);
                else %if that minimum distance is met, move on to the next dot
                    CurrDot = CurrDot+1;
                end
                
            end %end of loop across the dots in the current 3rd
            
        end %end of loop across all frames in a full cycle
        
        %adjust the x positions so they are in the centre of the screen. Do this for all frames in one go:
        %If we don't do this they will be shifted off to the left of screen
        AdjustXBy = (screenRect(3) - screenRect(4))/2; %shift X dot positions by this amount to centre them, since image size <screenRect(3)
        dot_pos(:,1,:) = dot_pos(:,1,:) + AdjustXBy;
        
        DotPosBinoc{TwoEyes} = dot_pos; %set aside the newly-assigned dots for this eye
        
    end %end of loop across eyes
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Sort all dots by smallest distances between L/R eyes
    %%%%-------------------------------------------------------------------------%%%%
    
    %Remember this just sorts the order by which the dot positions are listed in the matrices so
    %that the L & R eye dots are paired according to their smallest distance to each other: this means we can
    %assign color to L & R eye in an alternate manner (eg BWBWBW...; WBWBW) meaning that any dots appearing close to one
    %another across the 2 eyes *should* have the opposite color/contrast polarity (*though there may still be the odd occasion
    %where this does not happen).
    
    % *** %
    % On any given frame 1/3 of dots change positions from their previous location.
    % This ensures there are the same amount of transient events across frames (even the 1st).
    % we also need to make sure these changes are tracked so a dot doesn't accidentally change colour during its lifetime
    % when dot colour is assigned.
    % Because of the sorting below, position in the MATRIX now matters when assigning dot colour, but 
    % the indices we retain in 'MinDist' allow us to keep track of that.
    % *** %
    
    %Open a new figure if you want to plot the resulting positions: each frame will be drawn over
    %the top of the preceding one on a single figure
    if PlotResult, figure, end 
    
    % container for all min distance INDICES for each dot pair; for each frame
    % 1/3 of these indices change with each frame, ensuring each dot pairing lasts no longer than the 3 frame dot lifetime
    MinDist = nan(num_dots,FrmsFullCyc); 
    for F = 1:FrmsFullCyc %Need to do this for every frame
        
        Curr3rd = mod(F,3)+1; %tells us whether this is a 1st, 2nd or 3rd frame, & determines which 3rd of dots to change
        %NotCurr3rd = setdiff([1:3],Curr3rd)
        CurrRows = DotThirdIndices(Curr3rd,:);
        
        if F == 1
            %For the very first frame, we want to sort all of the dots
            dot_posL(:,:) = DotPosBinoc{1}(:,:,F);
            dot_posR(:,:) = DotPosBinoc{2}(:,:,F);
            D = pdist2(dot_posL(:,:), dot_posR(:,:));
            CurrRows = [1, num_dots];
        else
            % Make all values nan to take the 2/3 of dots we don't need to sort out
            dot_posL(:,:) = nan;
            dot_posR(:,:) = nan;
            %Then place the 1/3 of dot coordinates on this frame that we DO need to sort, back in
            dot_posL(CurrRows(1):CurrRows(2),:) = DotPosBinoc{1}(CurrRows(1):CurrRows(2),:,F);
            dot_posR(CurrRows(1):CurrRows(2),:) = DotPosBinoc{2}(CurrRows(1):CurrRows(2),:,F);
            %Now compute the distance matrix, excluding the 2/3 that do not need sorting
            D = pdist2(dot_posL(:,:), dot_posR(:,:));
            %Make MinDist the same as the previous frame (then we re-arrange 1/3)
            MinDist(:,F) = MinDist(:,F-1);
        end
        
        % Compute the distance matrix (unsorted):
        % Gives the distances between each pairwise combination of dots in each eye
        %D = pdist2(DotPosBinoc{1}(:,:,F), DotPosBinoc{2}(:,:,F));
        
        % Go through each ROW in the distance matrix, D, which correspond to each ROW in the L eye matrix, DotPosBinoc{1}
        % Once a dot in the R eye (DotPosBinoc{2}) has been assigned to a Left eye dot, it is removed from further distance computations.
        % So we would take the next smallest distance, then the next smallest, and so on.
        RE_ind_rep = true(num_dots,1); %set logical indices for RE dots (all dots for first iteration).
        
        for ix = CurrRows(1):CurrRows(2) %Loop across all dots in the current third
            
            %find the column with the smallest distance
            %MinDist indexes the COLUMNS of D and the ROWS of DotPosBinoc{2}
            %for each item 1..n of DotPosBinoc{1}
            MinDist(ix,F) = find(D(ix,:) == nanmin(D(ix,RE_ind_rep))); %Find the INDEX of the smallest distance between LE dot ix and all (unassigned/remaining) RE dots
            
            % Now update the RE index, to remove the dot with the smallest distance at this current (i'th)
            % iteration from further computations of the smallest distance for subsequent LE dots
            RE_ind_rep(MinDist(ix,F)) = false;
            
        end %end of loop across dots for current frame
        
        %So now we have the dots sorted according to their smallest possible distances.
        %(without replacement).
        
        % **** At this point, we previously would sort the dots in the dot position matrix according to the smallest distances.
        % However, we no longer need to do that. Because we record the indices of the dot pairings/sortings in MinDist,
        % we can just use that to sort the dot colour assignment matrices instead. Dot position matrices need not change.
        % Sorting the dot position matrix for the R eye only was part of the reason we encountered the dot colour flickering bug.        
        % *** So no longer doing this below *** :
        %DotPosBinoc{2}(:,:,F) = DotPosBinoc{2}(MinDist,:,F);
        
        %If you want to plot the resulting dot positions:
        % The figure below should reveal 2 things.
        % First, dots in the Left eye (rep. as dots) and Right eye (rep. as rings)
        % should only appear in alternating strips.
        % Second, for dots near the strip borders, whenever a L & R eye dot come close, they should have opposite color
        % (Ie if rings and dots are close, they should not be the same colour).
        % *** NOTE *** : this will only work now if you first sort the dot positions using the 'MinDist' indices
        
        %Plot left eye, alternating colours according to x position:
        if PlotResult
            plot(DotPosBinoc{1}(1:2:end,1,F), DotPosBinoc{1}(1:2:end,2,F), 'r.')
            hold on
            plot(DotPosBinoc{1}(2:2:end,1,F), DotPosBinoc{1}(2:2:end,2,F), 'b.')
            %Now plot right eye, in opposite colour alternation:
            %plot circles to see where dots appear in nearby locations
            plot(DotPosBinoc{2}(2:2:end,1,F), DotPosBinoc{2}(2:2:end,2,F), 'ro')
            plot(DotPosBinoc{2}(1:2:end,1,F), DotPosBinoc{2}(1:2:end,2,F), 'bo')
            axis square
            axis([1+AdjustXBy imsize+AdjustXBy 1 imsize])
            %Draw strip borders on if you wish:
            for ii = 1:StripWidth:imsize, refline(0, ii), end
            hold off
            pause
        end
        
    end %end of sorting loop across frames
    
    %Save the result:
    DotParams.StimulusNumber = T;   % store the number of the current set of dots (also in file name)
    strfreq = num2str(frequency);   % generate a string of the frequency
    strfreq(strfreq == '.') = '_';  % replace the decimal pt (if any) with underscore so it doesn't mess up the file name
    
    % Provide a different filename for coloured or achromatic stimuli:
    if MakeColour
        fname = fullfile('stimuli', ['IOVD_dots_colour_', strfreq, '_Hz_stim_', num2str(T), '.mat']);
    else
        fname = fullfile('stimuli', ['IOVD_dots_', strfreq, '_Hz_stim_', num2str(T), '.mat']);
    end
    % Note that MinDist is also now saved (14/3/16).
    save(fname, 'DotPosBinoc', 'DotParams', 'MinDist');
    toc
    
end %end of loop across trials/unique stimulus sets

