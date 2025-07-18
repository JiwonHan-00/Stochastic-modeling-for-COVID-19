function results = ModiGillespie_algorithm(input)
% Modified Gillespie Algorithm for COVID-19 Contact Tracing Simulation
% Two-group model: High-risk and Low-risk populations

%% Initialize Compartments
% Susceptible compartments
SH = input.SH0;     % Susceptible high-risk
SL = input.SL0;     % Susceptible low-risk

% Contact compartments (under surveillance)
CH = 0; CHt = 0;    % Contacts high-risk (normal/traced)
CL = 0; CLt = 0;    % Contacts low-risk (normal/traced)
QHC = 0; QLC = 0;   % Quarantined contacts

% Exposed compartments
E1H = input.E1H0; E2H = input.E2H0; E1Ht = 0; E2Ht = 0;  % High-risk exposed
E1L = input.E1L0; E2L = input.E2L0; E1Lt = 0; E2Lt = 0;   % Low-risk exposed

% Infectious compartments
I1H = input.I1H0; I2H = input.I2H0; I1Ht = 0; I2Ht = 0;  % High-risk infectious
I1L = input.I1L0; I2L = input.I2L0; I1Lt = 0; I2Lt = 0;   % Low-risk infectious

% Removed compartments
Q = 0;              % Quarantined
R = 0;              % Recovered

%% Initialize Time Arrays
TList = zeros(input.Iter_max, 1);
TList(1) = 0;
currentTime = 0;
reactionlist = zeros(input.Iter_max, 1);

%% Initialize Delay Event Arrays
% Quarantine release delays
DT_QHCtoSH = []; DT_QLCtoSL = [];

% Disease progression delays
DT_EHtoI1H = []; DT_EHtoI2H = []; DT_ELtoI1L = []; DT_ELtoI2L = [];

% Isolation delays
DT_I1HtQ = []; DT_I1LtQ = []; DT_I2HtQ = []; DT_I2LtQ = [];

% Recovery delays
DT_I2HtoR = []; DT_I2LtoR = []; DT_QtoR = [];

% Self-reporting delays
DT_I1Hrep = []; DT_I1Lrep = [];

% Contact tracing delays
DT_CHtQ = []; DT_CHtoSH = []; DT_CLtQ = []; DT_CLtoSL = [];
DT_E1HIt = []; DT_E2HIt = []; DT_E1LIt = []; DT_E2LIt = [];

%% Initialize Delayed Events for Initial Conditions
% High-risk group initial delays
if input.E1H0 > 0
    for i = 1:input.E1H0
        DT_EHtoI1H = [DT_EHtoI1H, currentTime + gamrnd(input.l1, input.l2)];
    end
end
if input.E2H0 > 0
    for i = 1:input.E2H0
        DT_EHtoI2H = [DT_EHtoI2H, currentTime + gamrnd(input.l1, input.l2)];
    end
end
if input.I1H0 > 0
    for i = 1:input.I1H0
        DT_I1Hrep = [DT_I1Hrep, currentTime + lognrnd(input.aH1, input.aH2)];
    end
end
if input.I2H0 > 0
    for i = 1:input.I2H0
        DT_I2HtoR = [DT_I2HtoR, currentTime + unifrnd(input.i1, input.i2)];
    end
end

% Low-risk group initial delays
if input.E1L0 > 0
    for i = 1:input.E1L0
        DT_ELtoI1L = [DT_ELtoI1L, currentTime + gamrnd(input.l1, input.l2)];
    end
end
if input.E2L0 > 0
    for i = 1:input.E2L0
        DT_ELtoI2L = [DT_ELtoI2L, currentTime + gamrnd(input.l1, input.l2)];
    end
end
if input.I1L0 > 0
    for i = 1:input.I1L0
        DT_I1Lrep = [DT_I1Lrep, currentTime + lognrnd(input.aL1, input.aL2)];
    end
end
if input.I2L0 > 0
    for i = 1:input.I2L0
        DT_I2LtoR = [DT_I2LtoR, currentTime + unifrnd(input.i1, input.i2)];
    end
end

%% Main Simulation Loop
time_index = 0;
for ii = 2:input.Iter_max
    
    % Display progress occasionally
    if time_index < currentTime
        if mod(time_index, 5) == 0
            disp(['time: ', num2str(time_index), ', S_H: ', num2str(SH), ...
                  ', S_L: ', num2str(SL), ', # active: ', ...
                  num2str(E1H+E2H+E1L+E2L+I1H+I2H+I1L+I2L)]);
        end
        time_index = time_index + 1;
    end
    
    % Update parameters based on social distancing
    if currentTime > input.SDtime
        % Post-social distancing
        input.cHH = input.cHH_post;
        input.cHL = input.cHL_post;
        input.cLH = input.cLH_post;
        input.cLL = input.cLL_post;
        input.CTrateH = input.CTrateH_post;
        input.CTrateL = input.CTrateL_post;
    else
        % Pre-social distancing
        input.cHH = input.cHH_pre;
        input.cHL = input.cHL_pre;
        input.cLH = input.cLH_pre;
        input.cLL = input.cLL_pre;
        input.CTrateH = input.CTrateH_pre;
        input.CTrateL = input.CTrateL_pre;
    end
    
    % Calculate current population sizes
    NH = SH + CH + CHt + E1H + E2H + E1Ht + E2Ht + I1H + I2H + I1Ht + I2Ht;
    NL = SL + CL + CLt + E1L + E1Lt + E2L + E2Lt + I1L + I2L + I1Lt + I2Lt;
    total_pop = NH + NL + R;
    
    % Calculate transmission rates
    infectious_H = I1H + I1Ht + I2H + I2Ht;
    infectious_L = I1L + I1Lt + I2L + I2Lt;
    
    reaction_contSH = SH * (input.cHH * infectious_H + input.cHL * infectious_L) / total_pop;
    reaction_contCH = CH * (input.cHH * infectious_H + input.cHL * infectious_L) / total_pop;
    reaction_contCHt = CHt * (input.cHH * infectious_H + input.cHL * infectious_L) / total_pop;
    reaction_contSL = SL * (input.cLH * infectious_H + input.cLL * infectious_L) / total_pop;
    reaction_contCL = CL * (input.cLH * infectious_H + input.cLL * infectious_L) / total_pop;
    reaction_contCLt = CLt * (input.cLH * infectious_H + input.cLL * infectious_L) / total_pop;
    
    reactiondice = [reaction_contSH, reaction_contCH, reaction_contCHt, ...
                   reaction_contSL, reaction_contCL, reaction_contCLt];
    
    % Calculate time to next reaction
    total_rate = sum(reactiondice);
    if total_rate > 0
        Tau = -1/total_rate * log(rand);
    else
        Tau = inf;
    end
    
    % Sort all delayed events
    DT_QHCtoSH = sort(DT_QHCtoSH); DT_QLCtoSL = sort(DT_QLCtoSL);
    DT_EHtoI1H = sort(DT_EHtoI1H); DT_EHtoI2H = sort(DT_EHtoI2H);
    DT_ELtoI1L = sort(DT_ELtoI1L); DT_ELtoI2L = sort(DT_ELtoI2L);
    DT_I1HtQ = sort(DT_I1HtQ); DT_I1LtQ = sort(DT_I1LtQ);
    DT_I2HtQ = sort(DT_I2HtQ); DT_I2LtQ = sort(DT_I2LtQ);
    DT_I2HtoR = sort(DT_I2HtoR); DT_I2LtoR = sort(DT_I2LtoR);
    DT_QtoR = sort(DT_QtoR); DT_I1Hrep = sort(DT_I1Hrep); DT_I1Lrep = sort(DT_I1Lrep);
    DT_CHtQ = sort(DT_CHtQ); DT_CHtoSH = sort(DT_CHtoSH);
    DT_CLtQ = sort(DT_CLtQ); DT_CLtoSL = sort(DT_CLtoSL);
    DT_E1HIt = sort(DT_E1HIt); DT_E2HIt = sort(DT_E2HIt);
    DT_E1LIt = sort(DT_E1LIt); DT_E2LIt = sort(DT_E2LIt);
    
    % Find next delayed event
    totstackTime = [DT_QHCtoSH, DT_QLCtoSL, DT_EHtoI1H, DT_EHtoI2H, DT_ELtoI1L, DT_ELtoI2L, ...
                   DT_I1HtQ, DT_I1LtQ, DT_I2HtQ, DT_I2LtQ, DT_I2HtoR, DT_I2LtoR, DT_QtoR, ...
                   DT_I1Hrep, DT_I1Lrep, DT_CHtQ, DT_CHtoSH, DT_CLtQ, DT_CLtoSL, ...
                   DT_E1HIt, DT_E2HIt, DT_E1LIt, DT_E2LIt];
    
    if ~isempty(totstackTime)
        minStack = min(totstackTime);
    else
        minStack = inf;
    end
    
    % Update time to next event
    if currentTime + Tau < minStack
        % Regular transmission event occurs first
        currentTime = currentTime + Tau;
        reactiondice = cumsum(reactiondice) / sum(reactiondice);
        rd = rand;
        
        if rd < reactiondice(1)
            % SH transmission
            SH = SH - 1;
            reactionlist(ii) = 24;
            if rand < input.pHH
                if rand < input.rhoH
                    E1H = E1H + 1;
                    reactionlist(ii) = 25;
                    DT_EHtoI1H = [DT_EHtoI1H, currentTime + gamrnd(input.l1, input.l2)];
                else
                    E2H = E2H + 1;
                    reactionlist(ii) = 26;
                    DT_EHtoI2H = [DT_EHtoI2H, currentTime + gamrnd(input.l1, input.l2)];
                end
            else
                CH = CH + 1;
                DT_CHtoSH = [DT_CHtoSH, currentTime + input.t_suv];
            end
            
        elseif rd < reactiondice(2)
            % CH transmission
            reactionlist(ii) = 27;
            if ~isempty(DT_CHtoSH)
                DT_CHtoSH = DT_CHtoSH(2:end);
            end
            if rand < input.pHH
                CH = CH - 1;
                if rand < input.rhoH
                    E1H = E1H + 1;
                    reactionlist(ii) = 28;
                    DT_EHtoI1H = [DT_EHtoI1H, currentTime + gamrnd(input.l1, input.l2)];
                else
                    E2H = E2H + 1;
                    reactionlist(ii) = 29;
                    DT_EHtoI2H = [DT_EHtoI2H, currentTime + gamrnd(input.l1, input.l2)];
                end
            else
                DT_CHtoSH = [DT_CHtoSH, currentTime + input.t_suv];
            end
            
        elseif rd < reactiondice(3)
            % CHt transmission
            reactionlist(ii) = 30;
            if ~isempty(DT_CHtQ)
                DT_CHtQ = DT_CHtQ(2:end);
            end
            if rand < input.pHH
                CHt = CHt - 1;
                if rand < input.rhoH
                    E1H = E1H + 1;
                    reactionlist(ii) = 31;
                    DT_EHtoI1H = [DT_EHtoI1H, currentTime + gamrnd(input.l1, input.l2)];
                else
                    E2H = E2H + 1;
                    reactionlist(ii) = 32;
                    DT_EHtoI2H = [DT_EHtoI2H, currentTime + gamrnd(input.l1, input.l2)];
                end
            else
                DT_CHtQ = [DT_CHtQ, currentTime + unifrnd(input.t1, input.t2)];
            end
            
        elseif rd < reactiondice(4)
            % SL transmission
            SL = SL - 1;
            reactionlist(ii) = 33;
            if rand < input.pLL
                if rand < input.rhoL
                    E1L = E1L + 1;
                    reactionlist(ii) = 34;
                    DT_ELtoI1L = [DT_ELtoI1L, currentTime + gamrnd(input.l1, input.l2)];
                else
                    E2L = E2L + 1;
                    reactionlist(ii) = 35;
                    DT_ELtoI2L = [DT_ELtoI2L, currentTime + gamrnd(input.l1, input.l2)];
                end
            else
                CL = CL + 1;
                DT_CLtoSL = [DT_CLtoSL, currentTime + input.t_suv];
            end
            
        elseif rd < reactiondice(5)
            % CL transmission
            reactionlist(ii) = 36;
            if ~isempty(DT_CLtoSL)
                DT_CLtoSL = DT_CLtoSL(2:end);
            end
            if rand < input.pLL
                CL = CL - 1;
                if rand < input.rhoL
                    E1L = E1L + 1;
                    reactionlist(ii) = 37;
                    DT_ELtoI1L = [DT_ELtoI1L, currentTime + gamrnd(input.l1, input.l2)];
                else
                    E2L = E2L + 1;
                    reactionlist(ii) = 38;
                    DT_ELtoI2L = [DT_ELtoI2L, currentTime + gamrnd(input.l1, input.l2)];
                end
            else
                DT_CLtoSL = [DT_CLtoSL, currentTime + input.t_suv];
            end
            
        else
            % CLt transmission
            reactionlist(ii) = 39;
            if ~isempty(DT_CLtQ)
                DT_CLtQ = DT_CLtQ(2:end);
            end
            if rand < input.pLL
                CLt = CLt - 1;
                if rand < input.rhoL
                    E1L = E1L + 1;
                    reactionlist(ii) = 40;
                    DT_ELtoI1L = [DT_ELtoI1L, currentTime + gamrnd(input.l1, input.l2)];
                else
                    E2L = E2L + 1;
                    reactionlist(ii) = 41;
                    DT_ELtoI2L = [DT_ELtoI2L, currentTime + gamrnd(input.l1, input.l2)];
                end
            else
                DT_CLtQ = [DT_CLtQ, currentTime + unifrnd(input.t1, input.t2)];
            end
        end
        
    else
        % Delayed event occurs first
        currentTime = minStack;
        
        % Determine which delayed event occurs
        minstack_array = [
            get_min_or_inf(DT_QHCtoSH), get_min_or_inf(DT_QLCtoSL), ...
            get_min_or_inf(DT_EHtoI1H), get_min_or_inf(DT_EHtoI2H), ...
            get_min_or_inf(DT_ELtoI1L), get_min_or_inf(DT_ELtoI2L), ...
            get_min_or_inf(DT_I1HtQ), get_min_or_inf(DT_I1LtQ), ...
            get_min_or_inf(DT_I2HtQ), get_min_or_inf(DT_I2LtQ), ...
            get_min_or_inf(DT_I2HtoR), get_min_or_inf(DT_I2LtoR), ...
            get_min_or_inf(DT_QtoR), get_min_or_inf(DT_I1Hrep), ...
            get_min_or_inf(DT_I1Lrep), get_min_or_inf(DT_CHtQ), ...
            get_min_or_inf(DT_CHtoSH), get_min_or_inf(DT_CLtQ), ...
            get_min_or_inf(DT_CLtoSL), get_min_or_inf(DT_E1HIt), ...
            get_min_or_inf(DT_E2HIt), get_min_or_inf(DT_E1LIt), ...
            get_min_or_inf(DT_E2LIt)
        ];
        
        [~, aa] = min(minstack_array);
        
        % Execute the delayed event
        if aa == 1  % QHC to SH
            QHC = QHC - 1; SH = SH + 1;
            DT_QHCtoSH = DT_QHCtoSH(2:end);
            reactionlist(ii) = 1;
        elseif aa == 2  % QLC to SL
            QLC = QLC - 1; SL = SL + 1;
            DT_QLCtoSL = DT_QLCtoSL(2:end);
            reactionlist(ii) = 2;
        elseif aa == 3  % E1H to I1H
            E1H = E1H - 1; I1H = I1H + 1;
            DT_EHtoI1H = DT_EHtoI1H(2:end);
            DT_I1Hrep = [DT_I1Hrep, currentTime + lognrnd(input.aH1, input.aH2)];
            reactionlist(ii) = 3;
        elseif aa == 4  % E2H to I2H
            E2H = E2H - 1; I2H = I2H + 1;
            DT_EHtoI2H = DT_EHtoI2H(2:end);
            DT_I2HtoR = [DT_I2HtoR, currentTime + unifrnd(input.i1, input.i2)];
            reactionlist(ii) = 4;
        elseif aa == 5  % E1L to I1L
            E1L = E1L - 1; I1L = I1L + 1;
            DT_ELtoI1L = DT_ELtoI1L(2:end);
            DT_I1Lrep = [DT_I1Lrep, currentTime + lognrnd(input.aL1, input.aL2)];
            reactionlist(ii) = 5;
        elseif aa == 6  % E2L to I2L
            E2L = E2L - 1; I2L = I2L + 1;
            DT_ELtoI2L = DT_ELtoI2L(2:end);
            DT_I2LtoR = [DT_I2LtoR, currentTime + unifrnd(input.i1, input.i2)];
            reactionlist(ii) = 6;
        elseif aa == 7  % I1Ht to Q
            I1Ht = I1Ht - 1; Q = Q + 1;
            DT_I1HtQ = DT_I1HtQ(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 7;
        elseif aa == 8  % I1Lt to Q
            I1Lt = I1Lt - 1; Q = Q + 1;
            DT_I1LtQ = DT_I1LtQ(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 8;
        elseif aa == 9  % I2Ht to Q
            I2Ht = I2Ht - 1; Q = Q + 1;
            DT_I2HtQ = DT_I2HtQ(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 9;
        elseif aa == 10  % I2Lt to Q
            I2Lt = I2Lt - 1; Q = Q + 1;
            DT_I2LtQ = DT_I2LtQ(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 10;
        elseif aa == 11  % I2H to R
            I2H = I2H - 1; R = R + 1;
            DT_I2HtoR = DT_I2HtoR(2:end);
            reactionlist(ii) = 11;
        elseif aa == 12  % I2L to R
            I2L = I2L - 1; R = R + 1;
            DT_I2LtoR = DT_I2LtoR(2:end);
            reactionlist(ii) = 12;
        elseif aa == 13  % Q to R
            Q = Q - 1; R = R + 1;
            DT_QtoR = DT_QtoR(2:end);
            reactionlist(ii) = 13;
        elseif aa == 14  % I1H self-report
            I1H = I1H - 1; Q = Q + 1;
            DT_I1Hrep = DT_I1Hrep(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 14;
        elseif aa == 15  % I1L self-report
            I1L = I1L - 1; Q = Q + 1;
            DT_I1Lrep = DT_I1Lrep(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 15;
        elseif aa == 16  % CHt to QHC
            CHt = CHt - 1; QHC = QHC + 1;
            DT_CHtQ = DT_CHtQ(2:end);
            DT_QHCtoSH = [DT_QHCtoSH, currentTime + input.t_isoC];
            reactionlist(ii) = 16;
        elseif aa == 17  % CH to SH
            CH = CH - 1; SH = SH + 1;
            DT_CHtoSH = DT_CHtoSH(2:end);
            reactionlist(ii) = 17;
        elseif aa == 18  % CLt to QLC
            CLt = CLt - 1; QLC = QLC + 1;
            DT_CLtQ = DT_CLtQ(2:end);
            DT_QLCtoSL = [DT_QLCtoSL, currentTime + input.t_isoC];
            reactionlist(ii) = 18;
        elseif aa == 19  % CL to SL
            CL = CL - 1; SL = SL + 1;
            DT_CLtoSL = DT_CLtoSL(2:end);
            reactionlist(ii) = 19;
        elseif aa == 20  % E1Ht to Q
            E1Ht = E1Ht - 1; Q = Q + 1;
            DT_E1HIt = DT_E1HIt(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 20;
        elseif aa == 21  % E2Ht to Q
            E2Ht = E2Ht - 1; Q = Q + 1;
            DT_E2HIt = DT_E2HIt(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 21;
        elseif aa == 22  % E1Lt to Q
            E1Lt = E1Lt - 1; Q = Q + 1;
            DT_E1LIt = DT_E1LIt(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 22;
        else  % aa == 23, E2Lt to Q
            E2Lt = E2Lt - 1; Q = Q + 1;
            DT_E2LIt = DT_E2LIt(2:end);
            DT_QtoR = [DT_QtoR, currentTime + input.t_iso];
            % Call contact tracing
            ContactTracing_Dice;
            reactionlist(ii) = 23;
        end
    end
    
    % Update time list
    TList(ii) = currentTime;
    
    % Check stopping conditions
    if E1H == 0 && E2H == 0 && E1Ht == 0 && E2Ht == 0 && E1L == 0 && E2L == 0 && ...
       E1Lt == 0 && E2Lt == 0 && I1Ht == 0 && I2Ht == 0 && I1H == 0 && I2H == 0 && ...
       I1Lt == 0 && I2Lt == 0 && I1L == 0 && I2L == 0
        break
    end
    if SH == 0 && SL == 0
        break
    end
    if currentTime > 34
        break
    end
end

%% Return Results
results.tspan = TList(1:ii);
results.reactionlist = reactionlist(1:ii);

end

function result = get_min_or_inf(array)
% Helper function to get minimum or infinity if empty
if isempty(array)
    result = inf;
else
    result = min(array);
end
end