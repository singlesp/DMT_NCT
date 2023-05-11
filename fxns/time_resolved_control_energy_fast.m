function [global_CE, regional_CE] = time_resolved_control_energy_fast(Anorm, T, TS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this version uses the controlability gramian to quickly calculate
    % control energies and can only be used in scenarios where the identity
    % matrix is used as the control strategy. Often will be scaled from the
    % slower calculations, the amount of which depends on how the B matrix 
    % is constructed in those cases (with 1's or 2's along the diagonal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Anorm - normalized adjacency matrix
    %T - time-horizon
    %TS - time-series overwhich to compute time-resolved control energy
    %TS should be size nparc x volumes
    % S. Parker Singleton, 2023

    %final and initial states are adjacent volumes
    x0 = TS(:,1:size(TS,2)-1);
    xf = TS(:,2:size(TS,2));
    
    % compute controlability gramian, used for simplifying calculation
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T


    % Loop over each transition
    [global_CE,regional_CE] = MIN_CONTROL_ENERGY(Anorm,WcI,x0,xf,T,false);

end
