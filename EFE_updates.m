
%% Expected free energy precisioin updates 
% Replication of code and text taken from  Step-by-Step Tutorial on Active Inference Modelling and its 
% Application to Empirical Data (Smith, Friston & Whyte, 2021) for personal
% understanding 

% By: Chinchella Nicola 


%--------------------------------------------------------------------------
% PRIOR DISTRIBUTION OVER POLICIES
% π0 = σ(ln E − γ G) --> encodes the prior distribution over policies

pi_0 = exp(log(E) - gamma*G) / sum(exp(log(E) -gamma*G));

% The initial distribution pi_O is made up of the learned prior over
% policies encoded in the E vector (habits - or reflecting the number of times
% a policy has been chosen before) and the expected free energy (G) of each
% allowable policy. 

%--------------------------------------------------------------------------
% POSTERIOR DISTRIBUTION OVER POLICIES
% π = σ(ln E − F − γ G) --> encodes the posterior distribution over
% policies

pi_posterior = exp(log(E) - gamma* G - F) / sum(exp(log(E) - gamma* G - F)); 

% The posterior distribution is determined by E, G, and the variational
% free energy (F) under each policy after making a new observation. The
% influence of G is modulated by gamma which encodes prior confidence in
% beliefs about G

%--------------------------------------------------------------------------
% CALCULATION OF PREDICTION ERROR 
% Gerror --> (π − π0) · (−G)

G_error = (pi_posterior - pi_0)' * -G;

% (pi_posterior - pi_0) --> captures how the new observation changes the
% agent's belief about policies 

% -G --> serves as a desirability or value function. The goal is to align
% the posterior belief with policies that minimize free energy. To
% understnad bettwe why we use -G try changing it and see what happens with
% the beta_update. Weighting by -G ensures that: 

    % 1. Positive deviations of pi_posterior from pi_0 (i.e., when posterior
    % belief increases for a policy) are rewarded if that policy has low G
    % (i.e., is more desirable)
    
    % 2. Negative deviations of pi_posterior from pi_0 (i.e., when
    % posterior belief decreases for a policy) are penalized if that policy
    % has high G (i.e., it is less desirable)

%--------------------------------------------------------------------------
% Beta (β) and beta_prior (β0) is an hyperparameter on the expected free energy precision term
% gamma (γ)

% Gamma (γ) controls the precision of G, based on the agent's confidence
% in its estimates of expected free energy. This confidence changes when
% new observations are consistent or inconsistent with G. More
% specifically, gamma modulates the influence of G on policy selection bsed
% on G_error. The difference between pi_posterior and pi_0 reflectes the
% extent to which new observations (scored by F) make policies more or less
% likely.

% If the vector encoding the posterior over policies increseas in magnitude
% in comparison to the prior, and still points in the same direction, the
% difference vector between the posterior and the prior will point in the
% same direction as the -G vector (i.e., less than 90% angle), hence the
% value of gamma will increase, thereby increasing the impact of G on
% policy selection. In constrast, if the difference vector between the
% posdterior and the prior does not point in the same direction as the -G
% vector, (i.e., greater than 90% angle) gamma will decrease and thereby
% reduce the impact of G on policy selection, as the agent's confidence in
% its estimates of expected free energy has decreased. 

% β ← 1/γ

beta_prior = 1/gamma;

beta_posterior = besta_prior; % At the beginning we can initialize the as the same 

% βupdate ← β − β0 + Gerror

beta_update = beta_posterior - beta_prior + G_error;

% β ← β − βupdate/ψ

beta_posterior = beta_posterior - beta_update/psi; % psi is a step size to promote stable
                                                   % convergence

% γ ← 1/β
gamma = 1/beta_posterior; % Update expected free energy precision

    
    
    
% When a new observation is inconsistent with prior beliefs about policies
% (pi_0 based on G), the agent asssigns a lower expected precision (gamma)
% to G when arriving at posteriors over policies