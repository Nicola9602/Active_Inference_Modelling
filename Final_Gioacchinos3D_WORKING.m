%% Gioacchino task
% Description of the task to be modelled. Nouns of graspable objects are
% used as stimuli. Participants had to decide if each noun referred to a
% natural or artificat, by performing either a precision or a power reach
% to grasp movement. Response grasp could be compatible or incompatible
% with the grasp typically used to manipulate the objects. They were
% instructed about the grasp at the beginning of the trial for example, preicison with natural objects and
% power with artificial, this didn't change within trial and 50% prob of one
% vs the other.

% There is a trede off between accuracy and velocity, so that one participant
% can choose to be accurate but slower in incongruent trials or quick but 
% inaccurate in incongruent trials 

% Insert habit in the sense that if I see a Big-N I have the habit to grasp
% with a Pw but i might have to overwrite that if I am in the other context
%

% Moreover, the words can suggest a power or precision
% grip (olive vs hammer), hence one can find oneself in a Compatible
% situation (Artifact-Power + Word-power) or in an incompatible
% (Artifact-power + Word-precision) and the same for natural stimuli.

% Instruction about whether artifact or natural should call for a precision
% vs power grip do not change in-trial, the numbers of trial is 512, with 8
% nouns natural + precision, 8 natural + power, 8 artifact + precision, 8
% artifact + power each repeated 16 times

% from now on Pw = power, Pr = precision, A = artifact, N = natural

clear
clear all 
rng('shuffle') 
%% 
% Simulation options after model building below:

% If Sim = 1, simulate single trial. This will reproduce fig. 8. (Although
            % note that, for this and the following simulations, results 
            % will vary each time due to random sampling)

% If Sim = 2, simulate multiple trials where the left context is active 
            % (D{1} = [1 0]'). This will reproduce fig. 10.
             
% If Sim = 3, simulate reversal learning, where the left context is active 
            % (D{1} = [1 0]') in early trials and then reverses in later 
            % trials (D{1} = [0 1]'). This will reproduce fig. 11.
            
% If Sim = 4, run parameter estimation on simulated data with reversal
            % learning. This will reproduce the top panel of fig. 17.
            
% If Sim = 5, run parameter estimation on simulated data with reversal
            % learning from multiple participants under different models
            % (i.e., different parameter values) and perform model comparison. 
            % This will reproduce the bottom panel of fig. 17. This option
            % will also save two structures that include results of model
            % comparison, model fitting, parameter recoverability analyses,
            % and inputs needed for group (PEB) analyses.

Sim = 2;

%% 
%==========================================================================

% D prior beliefs about intial states 
%--------------------------------------------------------------------------

% State that is told to the participant at the beginning of the trial. 

D{1} = [1 0]'; % Context state {'Pw-A + Pr-N','Pw-N + Pr-A'} 

% Time in trial 

D{2} = [1 0 0 0]'; % {'Start', 'Stim', 'Pw', 'Pr'}

D{3} = [0.25, 0.25, 0.25, 0.25]'; % Stimuli state {'A-Big', 'N-Big', 'A-Small', 'N-Small'} 


d{1} = [1 0]'; % The agent knows in which context it is 
 
d{2} = [1 0 0 0]'; %The agent knows he will start from the start state 

d{3} = D{3}*200;

% A probabilistic likelihood mapping
%--------------------------------------------------------------------------


%---- A{1} STIMULI OBSERVATION
% I start by specifying that all behavioural state lead to a null

for i = 1:4
    for k = 1:4
        A{1}(:,:,i,k) = [1 1; % Null
                   0 0; % A-Big
                   0 0; % N-Big
                   0 0; % A-Small
                   0 0]; % N-Small
    end
end 


% Instead Stim behavioral state leads to the observation of one of the four
% stimuli with equal probability in both contexts.  

%{'Pw-A + Pr-N','Pw-N + Pr-A'} 

Pstim = 0.25;

for i = 2
    for k = 1
        A{1}(:,:,i,k) = [0 0; % Null
                   Pstim Pstim; % A-Big
                   0 0; % N-Big
                   0 0; % A-Small
                   0 0]; % N-Small
    end 
end

for i = 2
    for k = 2
        A{1}(:,:,i,k) = [0 0; % Null
                   0 0; % A-Big
                   Pstim Pstim; % N-Big
                   0 0; % A-Small
                   0 0]; % N-Small
    end 
end

for i = 2
    for k = 3
        A{1}(:,:,i,k) = [0 0; % Null
                   0 0; % A-Big
                   0 0; % N-Big
                   Pstim Pstim; % A-Small
                   0 0]; % N-Small
    end 
end

for i = 2
    for k = 4
        A{1}(:,:,i,k) = [0 0; % Null
                   0 0; % A-Big
                   0 0; % N-Big
                   0 0; % A-Small
                   Pstim Pstim]; % N-Small
    end 
end

% A{2} --- Grasp Outomes observations correct incorrect 

for i = 1:4
    for k = 1:4
        A{2}(:,:,i,k) = [1 1;  % Null 
                   0 0;  % Correct
                   0 0]; % Incorrect
    end
end

% Pw
for i = 3
    for k = [1,3]
        A{2}(:,:,i,k) = [0 0;  % Null 
                   1 0;  % Correct
                   0 1]; % Incorrect
    end
end

for i = 3
    for k = [2,4]
        A{2}(:,:,i,k) = [0 0;  % Null 
                   0 1;  % Correct
                   1 0]; % Incorrect
    end
end

%Pr
for i = 4
    for k = [1,3]
        A{2}(:,:,i,k) = [0 0;  % Null 
                   0 1;  % Correct
                   1 0]; % Incorrect
    end
end

for i = 4
    for k = [2,4]
        A{2}(:,:,i,k) = [0 0;  % Null 
                   1 0;  % Correct
                   0 1]; % Incorrect
    end
end


for i = 1:4
    for k = 1:4
        A{3}(i,:,i,k) = [1 1];

    end
end


%--------------------------------------------------------------------------
% B Transition probabilities
%--------------------------------------------------------------------------
% Here I need one B for each D I have, therefore, three Bs. Here columns are 
% time t and rows are time t+1. Importantly in U/V you can only have the
% numbers you have in the third dimensions of B{2} (in this case)

% Context is stable within a trial

B{1}(:,:,1) = [1 0;
               0 1];

% Time in trial {'Start', 'Stim', 'Grasp'} the agent deterministically 
% transitions through trial sequence 

% If I set these Bs that it can move anywhere the multi trial simulation is
% really bad only one straight line 


B{2}(:,:,1) = [1 1 1 1; % Start
               0 0 0 0; % Stim
               0 0 0 0; % Pw
               0 0 0 0]; % Pr
B{2}(:,:,2) = [0 0 0 0; % Start
               1 1 1 1; % Stim
               0 0 0 0; % Pw
               0 0 0 0]; % Pr
B{2}(:,:,3) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 1 1 1; % Pw
               0 0 0 0]; % Pr
B{2}(:,:,4) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 0 0 0; % Pw
               0 1 1 1]; % Pr
%}

% When I set this, the multi trial simulation works way better, but still
% is not 100% convincing. There is something i am missing. 

%{
B{2}(:,:,1) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 1 1 1; % Pw
               0 0 0 0]; % Pr

B{2}(:,:,2) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 0 0 0; % Pw
               0 1 1 1]; % Pr


B{2}(:,:,3) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 1 1 1; % Pw
               0 0 0 0]; % Pr

B{2}(:,:,4) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 0 0 0; % Pw
               0 1 1 1]; % Pr
%}

% Stimuli state, stimuli is stable within trial

B{3}(:,:,1) = [1 0 0 0; % A-B
               0 1 0 0; % N-B
               0 0 1 0; % A-S
               0 0 0 1]; % N-S


b{1} = B{1}*200;
b{2} = B{2};
b{3} = B{3}*200;


%--------------------------------------------------------------------------
% C preferences
%--------------------------------------------------------------------------
% Here I specify preferences over each set of outcomes (C), with one matrix
% per outcome modality (numbers of A's). Here columns indicate time points, and rows
% indicate observations (same order as in the corresponding A matrices).
% Therefore columns ALWAYS = T

% The agent starts in start state (1), sees a stimuli (2), performs a grasp
% (3) and goes back to the start state (4)

T = 3;

No = [size(A{1}, 1) size(A{2},1) size(A{3}, 1)];
% Size gets the first dimension (Row) of A, hence the No matrix is a 5 x 3

% I have 3 A matrices so I will define 3 Cs.

C{1} = zeros(No(1), T); % Stimuli
C{2} = zeros(No(2), T); % Correct incorrect
C{3} = zeros(No(3), T); % Observed behaviour state factor

% Each C maps to its corrsponding A. The agent has no preference about
% being in the start state, seeing the stimuli, performing the grasp. 
% In the second case the agent has prefrence about his grasp, he prefers to
% perform the right grasps depending on the context. He prefers observing a
% grasp that is coherent with the context which I set as Cp (context
% preference)


C{2}(:,:) = [0 0  0  ;  % Null
             0 0  0.1;  % Correct
             0 0  0  ];  % Incorrect
                         



% Allowable policies: U or V. 
%==========================================================================

T ;  % number of timesteps defined above in C
Pi = 4; % number of policies
Nf = 3; % number of factors
% U = ones(1,Pi,Nf); % Context and stimuli are not controllable
% V = ones(T-1, Pi, Nf)

V = ones(T-1,Pi, Nf);

% 1. From the start state observe a stimulus
% 2. From the stimulus perform a Pw
% 3. From the stimulus perform a Pr
% 4. From grasp go back to start 

% U(:,:,1) is the context state hence is not controllab.e U(:,:,3) is the 
% stimulus state also not controllable. It can only control the grasp state 
% U(:,:,3).
% {'Null','Pw-Big','Pr-Small', 'Pw-Small', 'Pr-Big'}


V(:,:,2) = [1 2 3 4;
            1 2 3 4];


%V(:,:,2) = [1 2 3 4;
%            3 4 1 2];


%V(:,:,2) = [1 2 2 3 4;
%            2 3 4 1 1];


        
             
E = [0 0 0 0]';          

% MDP structure 
%==========================================================================

mdp.T = T; % Number of time steps

mdp.V = V; % Allowable policies

mdp.A = A; % state-outcome mapping 

mdp.B = B; % Transition probabilities

mdp.C = C; % preferred states 

mdp.D = D; % prior over initial states

mdp.E = E;


%mdp.a = a; % enable learning over likelihoods probabilities 

mdp.b = b; % enable learning over state transitions 

mdp.d = d; % enable learning priors over intial states


mdp.eta = 0.5; % Learning rate (0-1)

mdp.omega = 1; % Forgettin rate no forgetting (0-1)

mdp.beta = 1; % Lower values increase the influence of habits (E).

mdp.alpha = 32; % Action precision controls the level of randomness.

mdp.erp = 1; % Degree of belief resetting at each timestep 

mdp.tau = 12; % Times constant for evidence accumulation 

%--------------------------------------------------------------------------
% Labelling and plotting

% Adding some lables for subsequent plotting 

mdp.Aname = {'Stimuli', 'Correct/Incorrect', 'Observed Behaviour'};
mdp.Bname = {'Context', 'Time in Trial', 'Stimuli'};

% One label factor for each D 

label.factor{1} = 'Context State'; label.name{1} = {'Pw-A + Pr-N','Pw-N + Pr-A'};
label.factor{2} = 'Time in Trial'; label.name{2} = {'Start', 'Stim', 'Pw', 'Pr'};
label.factor{3} = 'Stimuli'; label.name{3} = {'A-Big', 'N-Big', 'A-Small', 'N-Small'};

% One label modeality for each A \
label.modality{1} = 'Stimuli'; label.outcome{1} = {'Null', 'A-Big', 'N-Big', 'A-Small', 'N-Small'};
label.modality{2} = 'Feedback'; label.outcome{2} = {'Null', 'Correct', 'Incorrect'};
label.modality{3} = 'Observed Behaviour'; label.outcome{3} = {'Start', 'Stim', 'Pw', 'Pr'};


% One label action for each B, usualyl not the first when it says that the
% context does not change. 

label.action{2} =  {'Pw', 'Pr'};
mdp.label = label;

% Check whether all matrix dimensions are correct 

mdp = spm_MDP_check(mdp);

if Sim == 1
%% Single trial simulation

% Run active inference estimantion 

MDP = spm_MDP_VB_X_tutorial(mdp);

% We can then use standard plotting routines to visualize simulated neural 
% responses

spm_figure('GetWin','Figure 1'); clf    % display behavior
spm_MDP_VB_LFP(MDP); 

%  and to show posterior beliefs and behavior:

spm_figure('GetWin','MDP'); clf    % display behavior
spm_MDP_VB_trial(MDP); 



elseif Sim == 2

N = 30; % number of trials

MDP = mdp;

[MDP(1:N)] = deal(MDP);

MDP = spm_MDP_VB_X_tutorial(MDP);

% We can again visualize simulated neural responses

spm_figure('GetWin','Figure 3'); clf    % display behavior
spm_MDP_VB_game_tutoria_G_Copy(MDP); 

end