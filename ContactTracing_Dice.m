% ContactTracing_Dice.m - Contact Tracing Generation Script
% This script is executed when someone enters quarantine (Q)
% It traces their contacts and moves them to traced compartments

% Determine contact tracing effectiveness based on current time
if currentTime > input.SDtime
    input.ctredH = 0.5;     % Reduced effectiveness post-social distancing
    input.ctredL = 0.5;
else
    input.ctredH = 1;       % Full effectiveness pre-social distancing
    input.ctredL = 1;
end

% Calculate number of contacts to trace based on the type of case
if (aa==7) || (aa==9) || (aa==14) || (aa==20) || (aa==21)
    % High-risk individual or high-risk traced individual enters Q
    etaLL = round(input.CTrateL * input.cLL * input.t_suv);
    etaHL = round(input.CTrateL * input.cHL * input.t_suv);
    etaLH = round(input.CTrateH * input.cLH * 4.04);
    etaHH = round(input.CTrateH * input.cHH * 4.04);
    eta = [etaHH, etaLH];
    eta = min([sum(eta), CH + CL + E1H + E2H + E1L + E2L + I1H + I2H + I1L + I2L]);
else
    % Low-risk individual or low-risk traced individual enters Q  
    etaLL = round(input.CTrateL * input.cLL * 4.04);
    etaHL = round(input.CTrateL * input.cHL * 4.04);
    etaLH = round(input.CTrateH * input.cLH * input.t_suv);
    etaHH = round(input.CTrateH * input.cHH * input.t_suv);
    eta = [etaHL, etaLL];
    eta = min([sum(eta), CH + CL + E1H + E2H + E1L + E2L + I1H + I2H + I1L + I2L]);
end

% Perform contact tracing for each contact
for i = 1:eta
    % Create probability distribution for different contact types
    tempunit = [input.ctredH*CH, input.ctredL*CL, input.ctredH*E1H, ...
                input.ctredH*E2H, input.ctredL*E1L, input.ctredL*E2L, ...
                I1H, I2H, I1L, I2L];
    
    if sum(tempunit) > 0
        dice_CT = cumsum(tempunit)/sum(tempunit);
        r1 = rand;
        
        if r1 < dice_CT(1)
            % Contact high-risk susceptible/contact (CH)
            if rand < 1 && CH > 0
                CH = CH - 1;
                CHt = CHt + 1;
                DT_CHtQ = [DT_CHtQ, currentTime + unifrnd(input.t1, input.t2)];
                if ~isempty(DT_CHtoSH)
                    DT_CHtoSH(randsample(1:length(DT_CHtoSH),1)) = [];
                end
            end
            
        elseif r1 < dice_CT(2)
            % Contact low-risk susceptible/contact (CL)
            if rand < 1 && CL > 0
                CL = CL - 1;
                CLt = CLt + 1;
                DT_CLtQ = [DT_CLtQ, currentTime + unifrnd(input.t1, input.t2)];
                if ~isempty(DT_CLtoSL)
                    DT_CLtoSL(randsample(1:length(DT_CLtoSL),1)) = [];
                end
            end
            
        elseif r1 < dice_CT(3)
            % Contact high-risk exposed E1H
            if rand < 1 && E1H > 0
                E1H = E1H - 1;
                E1Ht = E1Ht + 1;
                DT_E1HIt = [DT_E1HIt, currentTime + unifrnd(input.t1, input.t2)];
                if ~isempty(DT_EHtoI1H)
                    DT_EHtoI1H(randsample(1:length(DT_EHtoI1H),1)) = [];
                end
            end
            
        elseif r1 < dice_CT(4)
            % Contact high-risk exposed E2H
            if rand < 1 && E2H > 0
                E2H = E2H - 1;
                E2Ht = E2Ht + 1;
                DT_E2HIt = [DT_E2HIt, currentTime + unifrnd(input.t1, input.t2)];
                if ~isempty(DT_EHtoI2H)
                    DT_EHtoI2H(randsample(1:length(DT_EHtoI2H),1)) = [];
                end
            end
            
        elseif r1 < dice_CT(5)
            % Contact low-risk exposed E1L
            if rand < 1 && E1L > 0
                E1L = E1L - 1;
                E1Lt = E1Lt + 1;
                DT_E1LIt = [DT_E1LIt, currentTime + unifrnd(input.t1, input.t2)];
                if ~isempty(DT_ELtoI1L)
                    DT_ELtoI1L(randsample(1:length(DT_ELtoI1L),1)) = [];
                end
            end
            
        elseif r1 < dice_CT(6)
            % Contact low-risk exposed E2L
            if rand < 1 && E2L > 0
                E2L = E2L - 1;
                E2Lt = E2Lt + 1;
                DT_E2LIt = [DT_E2LIt, currentTime + unifrnd(input.t1, input.t2)];
                if ~isempty(DT_ELtoI2L)
                    DT_ELtoI2L(randsample(1:length(DT_ELtoI2L),1)) = [];
                end
            end
            
        elseif r1 < dice_CT(7)
            % Contact high-risk infectious I1H
            if rand < 1 && I1H > 0
                I1H = I1H - 1;
                I1Ht = I1Ht + 1;
                if ~isempty(DT_I1Hrep)
                    trackIDX = randsample(1:length(DT_I1Hrep),1);
                    DT_I1HtQ = [DT_I1HtQ, min(DT_I1Hrep(trackIDX), currentTime + unifrnd(input.t1, input.t2))];
                    DT_I1Hrep(trackIDX) = [];
                end
            end
            
        elseif r1 < dice_CT(8)
            % Contact high-risk infectious I2H
            if rand < 1 && I2H > 0
                I2H = I2H - 1;
                I2Ht = I2Ht + 1;
                if ~isempty(DT_I2HtoR)
                    trackIDX = randsample(1:length(DT_I2HtoR),1);
                    DT_I2HtQ = [DT_I2HtQ, min(DT_I2HtoR(trackIDX), currentTime + unifrnd(input.t1, input.t2))];
                    DT_I2HtoR(trackIDX) = [];
                end
            end
            
        elseif r1 < dice_CT(9)
            % Contact low-risk infectious I1L
            if rand < 1 && I1L > 0
                I1L = I1L - 1;
                I1Lt = I1Lt + 1;
                if ~isempty(DT_I1Lrep)
                    trackIDX = randsample(1:length(DT_I1Lrep),1);
                    DT_I1LtQ = [DT_I1LtQ, min(DT_I1Lrep(trackIDX), currentTime + unifrnd(input.t1, input.t2))];
                    DT_I1Lrep(trackIDX) = [];
                end
            end
            
        else
            % Contact low-risk infectious I2L
            if rand < 1 && I2L > 0
                I2L = I2L - 1;
                I2Lt = I2Lt + 1;
                if ~isempty(DT_I2LtoR)
                    trackIDX = randsample(1:length(DT_I2LtoR),1);
                    DT_I2LtQ = [DT_I2LtQ, min(DT_I2LtoR(trackIDX), currentTime + unifrnd(input.t1, input.t2))];
                    DT_I2LtoR(trackIDX) = [];
                end
            end
        end
    end
end