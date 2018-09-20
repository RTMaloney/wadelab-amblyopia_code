function Batch_run_FULL_Cue_thresholds_2IFC (SubjCode, CohOrCont)

% Batch file for running the FULL cue condition for MID thresholds.
% There is no motion preview at the start of each run, because I think it is unnecessary.
% We use 5 amplitudes * 5 frequencies for the FULL cue condition. The amplitudes are the same (range) as IOVD.
% R Maloney, Oct 2016

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
%NumCues = 2; we do not differentiate by 'cue' now.
% all conditions in a single repetition (ie a single run for each condition) ignoring threshold type:
NumRuns = NumAmplitudes * NumFrequencies;
RunRepeats = 3; %how many repeat runs of each condition to perform

% Set up frequencies. Multiplying the result by 100 & then dividing by 100 gives the values to 2 dec points
% (important for when the values are used in the data file names).
frequencies = round(logspace(log10(0.5), log10(8), NumFrequencies) * 100) / 100;
frequencies = reshape(repmat(frequencies,NumAmplitudes,1), NumRuns,1); % Each frequency by each amplitude

% Now do the same for amplitudes:
amplitudes = round(logspace(log10(1.67), log10(167), NumAmplitudes) * 100) / 100; % NOTE: these are the same as the IOVD amplitudes.
%Let's just fix a few of the overlapping values so they are the same in the 2 cues:
amplitudes(3) = 16.67;  %was 16.7
amplitudes = repmat(amplitudes',NumFrequencies,1);

% Determine whether this subject has already started testing or not
RunOrderFileName = fullfile('data', [SubjCode, '_MID_thresholds_FULLcue_RunOrder_', CohOrCont, '.mat']);

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
    
    % Execute the current run!
    MID_Full_Cue_thresholds_2IFC_M16 (SubjCode, CohOrCont, frequencies(RunOrder(R)), amplitudes(RunOrder(R)), runN)
    %Increment the test number, in case we need to start again
    ExpIncrement = ExpIncrement + 1
    % And save it:
    save(RunOrderFileName, 'RunOrder', 'ExpIncrement')
    
end