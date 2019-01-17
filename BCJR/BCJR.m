function [mapLLRs,output] = BCJR(trellis,coded)
% =========================================================================
% Max-Log BCJR (forward-backward algorithm)
% -------------------------------------------------------------------------
%   Max-log soft-in soft-out decoder [Bahl et.al.]. The decoder has two
%   outputs: extrinsic a posteriori LLRs for coded bits and sliced output
%   estimates (bits). Use the mex version for significant speedups!
%
%   [mapLLRs,output] = BCJR(trellis,coded)
%
%     trellis : trellis structure
%     coded   : input LLRs
%     output  : output bit estimates 
%     mapLLRs : output LLRs (intrinsic output LLR values)
%
% -------------------------------------------------------------------------
% Revision history:
%   - may-02-2018  v1.2  cs: cleanup for release
%   - may-29-2006  v1.1  cs: initial version 
% -------------------------------------------------------------------------
% (c) 2018 Sven Jacobsson and Christoph Studer 
% e-mail: sven.jacobsson@ericsson.com and studer@cornell.edu 
%=========================================================================

  warning('Use MEX version of BCJR decoder for significant speedups!');
 
  % max-log LLR scaling factor
  epsilon = 1.0; % set to 1.0 if you don't know what this is doing
  
  % --  extract trellis information
  numInputSymbols  = trellis.numInputSymbols;
  numOutputSymbols = trellis.numOutputSymbols;
  numStates        = trellis.numStates;
  nextStates       = trellis.nextStates;
  outputs          = trellis.outputs;

  % -- required parameters
  numInputBits  = log2(numInputSymbols);
  numOutputBits = log2( numOutputSymbols);
  numLLRs       = length(coded);
  numBranches   = numLLRs/numOutputBits;  
  bitTable      = de2bi((0:2^numOutputBits-1),numOutputBits);
  
  % -- INPUT METRICS : precompute gammas
  gamma = zeros(numOutputSymbols,numBranches+1);
  for l=1:numBranches
    for o=1:numOutputSymbols
      %% SLOW: Bits = 2*(de2bi(o-1,numOutputBits)-0.5);
      Bits = 2*(bitTable(o,:)-0.5);      
      for b=1:numOutputBits
        % MSB & LSB exchanged
        gamma(o,l) = gamma(o,l) + (1/2)*Bits(b)*coded(1+(l)*numOutputBits-(b));
      end
    end
  end
    
  % -- FORWARD ITERATION : compute alpha metrics  
  alpha = -ones(numStates,numBranches+1)*inf;
  alpha(1,1) = 0; % force initial state  
  for l=1:numBranches
    % -- calculate next alphas
    for s=1:numStates        
      % get current metric (alpha)
      curMetric = alpha(s,l);  
      for i=1:numInputSymbols
        % get gamma value for current branch
        curGamma = gamma(1+outputs(s,i),l);
        % update alpha for next state (max-log)
        alpha(1+nextStates(s,i),l+1) = max(alpha(1+nextStates(s,i),l+1),curMetric+curGamma);        
      end                    
    end
  end 
  
  % -- BACKWARD ITERATION : compute beta metrics  
  beta = -ones(numStates,numBranches+1)*inf;
  beta(1,numBranches+1) = 0; % force initial state  
  for l=numBranches:-1:1
    % -- calculate next alphas
    for s=1:numStates              
      for i=1:numInputSymbols 
        % get current metric (beta)
        curMetric = beta(1+nextStates(s,i),l+1);  
        % get gamma value for current branch
        curGamma = gamma(1+outputs(s,i),l);
        % update beta for next state (max-log)
        beta(s,l) = max(beta(s,l),curMetric+curGamma);
      end                    
    end
  end

  % -- CALCULATE OUTPUTS  
  px0llr = -ones(1,numLLRs)*inf; 
  px1llr = -ones(1,numLLRs)*inf;   
  px0out = -ones(1,numInputBits*numBranches)*inf;
  px1out = -ones(1,numInputBits*numBranches)*inf;
  for l=1:numBranches     
    for s=1:numStates      
      for i=1:numInputSymbols
        % -- sigma (output) metric
        curMetric = alpha(s,l) + beta(1+nextStates(s,i),l+1) + gamma(1+outputs(s,i),l);        
        % -- calculate estimated output bits
        if i==1
           px0out(1,l) = max(px0out(1,l),curMetric);
        else
           px1out(1,l) = max(px1out(1,l),curMetric); 
        end
        % -- calculate extrinsic information
        %% SLOW: Bits = de2bi(outputs(s,i),numOutputBits);
        Bits = bitTable(outputs(s,i)+1,:);   
        for b=[1:numOutputBits]
          % *** BUG *** (just kidding)
          if Bits(1,b)==0
            px0llr(1,1+l*numOutputBits-b) = max(px0llr(1,1+l*numOutputBits-b),curMetric);
          else
            px1llr(1,1+l*numOutputBits-b) = max(px1llr(1,1+l*numOutputBits-b),curMetric);
          end
        end     
      end
    end
  end

  % return LLR metrics and output estimates
  output  = 0.5*(sign(px1out-px0out)+1);
  mapLLRs = (px1llr-px0llr)*epsilon;
  
return