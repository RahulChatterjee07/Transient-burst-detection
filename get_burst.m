% function for extracting burst
% April 10, 2021, Rahul Chatterjee
% burst_base, burst_base_moved etc are organized data cell of {Irial}, each cell is a [channel*time] matirx
% burst_base_pmov etc are organized data cell of {Irial*channel}, each cell
% is a {#burst (e.g, 1st burst, 2nd burst ...)} cell


function [burst_base_pmov, burst_base_moving, burst_base_pmoved, ...
    burst_p90_pmov, burst_p90_moving, burst_p90_moved]...
    = get_burst(burst_base, burst_base_moved, burst_90, burst_90_moved, fs, CH)
    % ---- extract envelope to determine threshold ---- %
    for iTrial = 1:50
        for iChannel = 1:length(CH)
            [trial_base_pmov(iTrial,iChannel,:),~] = smoothdata(envelope(burst_base{iTrial}(CH(iChannel),0.26*fs+1:1.26*fs)),'gaussian');
            [trial_base_moving(iTrial,iChannel,:),~] = smoothdata(envelope(burst_base{iTrial}(CH(iChannel),1.76*fs+1:4.76*fs)),'gaussian');
            [trial_base_moved(iTrial,iChannel,:),~] = smoothdata(envelope(burst_base_moved{iTrial}(CH(iChannel),1.56*fs+1:2.06*fs)),'gaussian');
            [trial_90_pmov(iTrial,iChannel,:),~] = smoothdata(envelope(burst_90{iTrial}(CH(iChannel),0.26*fs+1:1.26*fs)),'gaussian');
            [trial_90_moving(iTrial,iChannel,:),~] = smoothdata(envelope(burst_90{iTrial}(CH(iChannel),1.76*fs+1:4.76*fs)),'gaussian');
            [trial_90_moved(iTrial,iChannel,:),~] = smoothdata(envelope(burst_90_moved{iTrial}(CH(iChannel),1.56*fs+1:2.06*fs)),'gaussian');
        end
    end
    env.base_pmov = squeeze(mean(trial_base_pmov,1)); % mean envelope over 50 trials, 50*8 matrix
    env.p90_pmov = squeeze(mean(trial_90_pmov,1));
    env.base_moving = squeeze(mean(trial_base_moving,1));
    env.p90_moving = squeeze(mean(trial_90_moving,1));
    env.base_moved = squeeze(mean(trial_base_moved,1));
    env.p90_moved = squeeze(mean(trial_90_moved,1));
   
    for iTrial = 1:50
        for iChannel = 1:length(CH)
            Ediff.base_pmov(iTrial,iChannel,:) = trial_base_pmov(iTrial,iChannel,:) - 0.75 * mean(env.base_pmov(iChannel,:));
            Ediff.base_moving(iTrial,iChannel,:) = trial_base_moving(iTrial,iChannel,:) - 0.75 * mean(env.base_moving(iChannel,:));
            Ediff.base_moved(iTrial,iChannel,:) = trial_base_moved(iTrial,iChannel,:) - 0.75 * mean(env.base_moved(iChannel,:));
            Ediff.p90_pmov(iTrial,iChannel,:) = trial_90_pmov(iTrial,iChannel,:) - 0.75 * mean(env.p90_pmov(iChannel,:));
            Ediff.p90_moving(iTrial,iChannel,:) = trial_90_moving(iTrial,iChannel,:) - 0.75 * mean(env.p90_moving(iChannel,:));
            Ediff.p90_moved(iTrial,iChannel,:) = trial_90_moved(iTrial,iChannel,:) - 0.75 * mean(env.p90_moved(iChannel,:));
            
            idx.base_pmov{iTrial,1}{iChannel} = find(Ediff.base_pmov(iTrial,iChannel,:) > 0);
            idx.base_moving{iTrial,1}{iChannel} = find(Ediff.base_moving(iTrial,iChannel,:) > 0);
            idx.base_moved{iTrial,1}{iChannel} = find(Ediff.base_moved(iTrial,iChannel,:) > 0);
            I.base_pmov{iTrial,1}{iChannel} = find(diff(idx.base_pmov{iTrial,1}{iChannel}) ~= 1);
            if isempty(I.base_pmov{iTrial,1}{iChannel})
                I.base_pmov{iTrial,1}{iChannel} = length(idx.base_pmov{iTrial,1}{iChannel});
            end
            I.base_moving{iTrial,1}{iChannel} = find(diff(idx.base_moving{iTrial}{iChannel}) ~= 1);
            if isempty(I.base_moving{iTrial,1}{iChannel})
                I.base_moving{iTrial,1}{iChannel} = length(idx.base_moving{iTrial,1}{iChannel});
            end   
            I.base_moved{iTrial,1}{iChannel} = find(diff(idx.base_moved{iTrial}{iChannel}) ~= 1);
            if isempty(I.base_moved{iTrial,1}{iChannel})
                I.base_moved{iTrial,1}{iChannel} = length(idx.base_moved{iTrial,1}{iChannel});
            end   
            idx.p90_pmov{iTrial,1}{iChannel} = find(Ediff.p90_pmov(iTrial,iChannel,:) > 0);
            idx.p90_moving{iTrial,1}{iChannel} = find(Ediff.p90_moving(iTrial,iChannel,:) > 0);
            idx.p90_moved{iTrial,1}{iChannel} = find(Ediff.p90_moved(iTrial,iChannel,:) > 0);
            I.p90_pmov{iTrial,1}{iChannel} = find(diff(idx.p90_pmov{iTrial,1}{iChannel}) ~= 1);
            if isempty(I.p90_pmov{iTrial,1}{iChannel})
                I.p90_pmov{iTrial,1}{iChannel} = length(idx.p90_pmov{iTrial,1}{iChannel});
            end
            I.p90_moving{iTrial,1}{iChannel} = find(diff(idx.p90_moving{iTrial}{iChannel}) ~= 1);
            if isempty(I.p90_moving{iTrial,1}{iChannel})
                I.p90_moving{iTrial,1}{iChannel} = length(idx.p90_moving{iTrial,1}{iChannel});
            end
            I.p90_moved{iTrial,1}{iChannel} = find(diff(idx.p90_moved{iTrial}{iChannel}) ~= 1);
            if isempty(I.p90_moved{iTrial,1}{iChannel})
                I.p90_moved{iTrial,1}{iChannel} = length(idx.p90_moved{iTrial,1}{iChannel});
            end
        end
    end
    
    for iTrial = 1:50
        for iChannel = 1:length(CH)
            startpoint = 1;    % sets start index at the first value in your array
            for i = 1:size(I.base_pmov{iTrial}{iChannel},1)
                End_Idx = I.base_pmov{iTrial}{iChannel}(i);   %set end index
                seq.base_pmov{iTrial,1}{iChannel}{i} = idx.base_pmov{iTrial}{iChannel}(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                startpoint = End_Idx+1;   %update start index for the next consecutive sequence
            end
            
            startpoint = 1;    % sets start index at the first value in your array
            for i = 1:size(I.base_moving{iTrial}{iChannel},1)
                End_Idx = I.base_moving{iTrial}{iChannel}(i);   %set end index
                seq.base_moving{iTrial,1}{iChannel}{i} = idx.base_moving{iTrial}{iChannel}(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                startpoint = End_Idx+1;   %update start index for the next consecutive sequence
            end
            
            startpoint = 1;    % sets start index at the first value in your array
            for i = 1:size(I.base_moved{iTrial}{iChannel},1)
                End_Idx = I.base_moved{iTrial}{iChannel}(i);   %set end index
                seq.base_moved{iTrial,1}{iChannel}{i} = idx.base_moved{iTrial}{iChannel}(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                startpoint = End_Idx+1;   %update start index for the next consecutive sequence
            end
            
            startpoint = 1;    % sets start index at the first value in your array
            for i = 1:size(I.p90_pmov{iTrial}{iChannel},1)
                End_Idx = I.p90_pmov{iTrial}{iChannel}(i);   %set end index
                seq.p90_pmov{iTrial,1}{iChannel}{i} = idx.p90_pmov{iTrial}{iChannel}(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                startpoint = End_Idx+1;   %update start index for the next consecutive sequence
            end
            
            startpoint = 1;    % sets start index at the first value in your array
            for i = 1:size(I.p90_moving{iTrial}{iChannel},1)
                End_Idx = I.p90_moving{iTrial}{iChannel}(i);   %set end index
                seq.p90_moving{iTrial,1}{iChannel}{i} = idx.p90_moving{iTrial}{iChannel}(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                startpoint = End_Idx+1;   %update start index for the next consecutive sequence
            end
            
            startpoint = 1;    % sets start index at the first value in your array
            for i = 1:size(I.p90_moved{iTrial}{iChannel},1)
                End_Idx = I.p90_moved{iTrial}{iChannel}(i);   %set end index
                seq.p90_moved{iTrial,1}{iChannel}{i} = idx.p90_moved{iTrial}{iChannel}(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
                startpoint = End_Idx+1;   %update start index for the next consecutive sequence
            end
        end
    end
    
    for iTrial = 1:50
        for iChannel = 1:1:length(CH)    
            
            startpoint = 1;            
            if size(seq.base_pmov{iTrial}{iChannel},2)
                for BST = 1:size(seq.base_pmov{iTrial}{iChannel},2)
                    if size(seq.base_pmov{iTrial}{iChannel}{BST},1) > fs/30
                        seq_burst.base_pmov{iTrial,iChannel}{startpoint,1} = seq.base_pmov{iTrial}{iChannel}{BST};
                        startpoint = startpoint + 1;
                    end
                end
            end
            
            startpoint = 1;            
            if size(seq.base_moving{iTrial}{iChannel},2)
                for BST = 1:size(seq.base_moving{iTrial}{iChannel},2)
                    if size(seq.base_moving{iTrial}{iChannel}{BST},1) > fs/30
                        seq_burst.base_moving{iTrial,iChannel}{startpoint,1} = seq.base_moving{iTrial}{iChannel}{BST};
                        startpoint = startpoint + 1;
                    end
                end
            end
            
            startpoint = 1;            
            if size(seq.base_moved{iTrial}{iChannel},2)
                for BST = 1:size(seq.base_moved{iTrial}{iChannel},2)
                    if size(seq.base_moved{iTrial}{iChannel}{BST},1) > fs/30
                        seq_burst.base_moved{iTrial,iChannel}{startpoint,1} = seq.base_moved{iTrial}{iChannel}{BST};
                        startpoint = startpoint + 1;
                    end
                end
            end
            
            startpoint = 1;            
            if size(seq.p90_pmov{iTrial}{iChannel},2)
                for BST = 1:size(seq.p90_pmov{iTrial}{iChannel},2)
                    if size(seq.p90_pmov{iTrial}{iChannel}{BST},1) > fs/30
                        seq_burst.p90_pmov{iTrial,iChannel}{startpoint,1} = seq.p90_pmov{iTrial}{iChannel}{BST};
                        startpoint = startpoint + 1;
                    end
                end
            end
            
            startpoint = 1;            
            if size(seq.p90_moving{iTrial}{iChannel},2)
                for BST = 1:size(seq.p90_moving{iTrial}{iChannel},2)
                    if size(seq.p90_moving{iTrial}{iChannel}{BST},1) > fs/30
                        seq_burst.p90_moving{iTrial,iChannel}{startpoint,1} = seq.p90_moving{iTrial}{iChannel}{BST};
                        startpoint = startpoint + 1;
                    end
                end
            end              
            
            startpoint = 1;            
            if size(seq.p90_moved{iTrial}{iChannel},2)
                for BST = 1:size(seq.p90_moved{iTrial}{iChannel},2)
                    if size(seq.p90_moved{iTrial}{iChannel}{BST},1) > fs/30
                        seq_burst.p90_moved{iTrial,iChannel}{startpoint,1} = seq.p90_moved{iTrial}{iChannel}{BST};
                        startpoint = startpoint + 1;
                    end
                end
            end              
        end
    end
    
    for iTrial = 1:50
        for iChannel = 1:length(CH)
            for j = 1:size(seq_burst.base_pmov{iTrial,iChannel},1)
                burst_base_pmov{iTrial,iChannel}{j,1} = squeeze(trial_base_pmov(iTrial,iChannel,seq_burst.base_pmov{iTrial,iChannel}{j}(1) : seq_burst.base_pmov{iTrial,iChannel}{j}(end)));
                burst_dur_base_pmov{iTrial,iChannel}{j,1} = size(seq_burst.base_pmov{iTrial,iChannel}{j},1);
            end
            for j = 1:size(seq_burst.base_moving{iTrial,iChannel},1)
                burst_base_moving{iTrial,iChannel}{j,1} = squeeze(trial_base_moving(iTrial,iChannel,seq_burst.base_moving{iTrial,iChannel}{j}(1) : seq_burst.base_moving{iTrial,iChannel}{j}(end)));
                burst_dur_base_moving{iTrial,iChannel}{j,1} = size(seq_burst.base_moving{iTrial,iChannel}{j},1);
            end
            for j = 1:size(seq_burst.base_moved{iTrial,iChannel},1)
                burst_base_pmoved{iTrial,iChannel}{j,1} = squeeze(trial_base_moved(iTrial,iChannel,seq_burst.base_moved{iTrial,iChannel}{j}(1) : seq_burst.base_moved{iTrial,iChannel}{j}(end)));
                burst_dur_base_moved{iTrial,iChannel}{j,1} = size(seq_burst.base_moved{iTrial,iChannel}{j},1);
            end
            for j = 1:size(seq_burst.p90_pmov{iTrial,iChannel},1)
                burst_p90_pmov{iTrial,iChannel}{j,1} = squeeze(trial_90_pmov(iTrial,iChannel,seq_burst.p90_pmov{iTrial,iChannel}{j}(1) : seq_burst.p90_pmov{iTrial,iChannel}{j}(end)));
                burst_dur_p90_pmov{iTrial,iChannel}{j,1} = size(seq_burst.p90_pmov{iTrial,iChannel}{j},1);
            end
            for j = 1:size(seq_burst.p90_moving{iTrial,iChannel},1)
                burst_p90_moving{iTrial,iChannel}{j,1} = squeeze(trial_90_moving(iTrial,iChannel,seq_burst.p90_moving{iTrial,iChannel}{j}(1) : seq_burst.p90_moving{iTrial,iChannel}{j}(end)));
                burst_dur_p90_moving{iTrial,iChannel}{j,1} = size(seq_burst.p90_moving{iTrial,iChannel}{j},1);
            end
            for j = 1:size(seq_burst.p90_moved{iTrial,iChannel},1)
                burst_p90_moved{iTrial,iChannel}{j,1} = squeeze(trial_90_moved(iTrial,iChannel,seq_burst.p90_moved{iTrial,iChannel}{j}(1) : seq_burst.p90_moved{iTrial,iChannel}{j}(end)));
                burst_dur_p90_moved{iTrial,iChannel}{j,1} = size(seq_burst.p90_moved{iTrial,iChannel}{j},1);
            end
        end
    end    
end