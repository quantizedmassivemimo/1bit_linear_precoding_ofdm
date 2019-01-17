function res = precoder_linear_ofdm_sim(varargin)
% =========================================================================
% Simulator for "Linear precoding with low-resolution DACs for massive 
% MU-MIMO-OFDM downlink"
% -------------------------------------------------------------------------
% Revision history:
%   - jan-17-2019  v0.1   sj: simplified/commented code for GitHub
% -------------------------------------------------------------------------
% (c) 2019 Sven Jacobsson and Christoph Studer 
% e-mail: sven.jacobsson@ericsson.com and studer@cornell.edu 
% -------------------------------------------------------------------------
% If you this simulator or parts of it, then you must cite our paper:
% -- S. Jacobsson, G. Durisi, M. Coldrey, and C. Studer, “Linear precoding 
% with low-resolution DACs for massive MU-MIMO-OFDM downlink,” IEEE Trans.
% Wireless Commun., Jun., to appear.
%=========================================================================

    % -- set up default/custom parameters
    if isempty(varargin)

        % set default simulation parameters
        disp('using default simulation settings and parameters...');
        par.runId = 1; % simulation ID (used to reproduce results)
        par.save = false; % save results (true,false)
        par.plot = true; % plot results (true,false)
        par.B = 32; % number of BS antenns
        par.U = 8; % number of single-antenna UEs
        par.S = 300; % number of occupied subcarriers: '72', '144', '300', '600', '1200'
        par.T = 10; % number of channels taps
        par.L = 2; % number of DAC levels per I and Q dimension
        par.approximation = 'rounding'; % distortion model: 'none', 'diagonal', 'rounding'
        par.sampling_rate_multiplier = 1; % sampling rate multiplier
        par.nrof_channels = 1e1; % number of channel realizations
        par.nrof_ofdm_symbols = 1e1; % number of OFDM symbols per channel realization
        par.precoder = {'ZFQ', 'ZFI'}; % select precoding scheme(s) to be evaluated
        par.SNRdB_list = -10:3:20; % list of SNR [dB] values to be simulated
        par.relerr = 0; % relative channel estimate error
        par.mod = 'QPSK'; % modulation type: 'QPSK', '8PSK', '16QAM','64QAM'
        par.code.rate = 1; % code rate: '1/2','3/4','2/3','5/6', '1'

    else

        % use custom simulation parameters
        disp('use custom simulation settings and parameters...')
        par = varargin{1}; % load custom simulation parameters

    end

    % -- initialization
    
    % use runId random seed (enables reproducibility)
    rng(par.runId);
    
    % set unique filename
    par.simName = ['PRE_',num2str(par.U),'x',num2str(par.B),'_',num2str(par.mod),'_',num2str(par.runId),' - ',datestr(clock,0)];  % simulation name (used for saving results)
    
    % -- plotting properties
    if par.plot == 1
        
        % number of rows/columns in constellation/spectrum plots
        if length(par.precoder)==1
            nd = 1;
        elseif length(par.precoder)<=6
            nd = 2;
        elseif length(par.precoder)<=9
            nd = 3;
        elseif length(par.precoder)<=12
            nd = 4;
        end
        md = ceil(length(par.precoder)/nd);
        
        % marker properties
        marker_style = {'--o','--s','--v','--+','--<','-->','--x','--^','--*','--d','--h','--p'};
        marker_color = [0.000, 0.447, 0.741;0.850, 0.325, 0.098; 0.929, 0.694, 0.125; 0.494, 0.184, 0.556;
                    0.466, 0.674, 0.188; 0.301, 0.745, 0.933; 0.635, 0.078, 0.184];
                
    end

    % set up constellation alphabets
    switch (par.mod)
        case 'QPSK'
            par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
        case '8PSK'
            par.symbols = [...
                exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
                exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
                exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
                exp(1i*2*pi/8*4), exp(1i*2*pi/8*5)];
        case '16QAM'
            par.symbols = [...
                -3-3i,-3-1i,-3+3i,-3+1i, ...
                -1-3i,-1-1i,-1+3i,-1+1i, ...
                +3-3i,+3-1i,+3+3i,+3+1i, ...
                +1-3i,+1-1i,+1+3i,+1+1i ];
        case '64QAM'
            par.symbols = [...
                -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    end

    % normalize symbol energy
    par.symbols = par.symbols/sqrt(mean(abs(par.symbols).^2));

    % subcarrier allocation and size of discrete Fourier transform
    switch par.S
        case 72
            par.N = 128*par.sampling_rate_multiplier;
            par.submap = [par.N-par.S/2+1:par.N, 2:ceil(par.S/2)+1];
        case 144
            par.N = 256*par.sampling_rate_multiplier;
            par.submap = [par.N-par.S/2+1:par.N, 2:ceil(par.S/2)+1];
        case 300
            par.N = 512*par.sampling_rate_multiplier;
            par.submap = [par.N-par.S/2+1:par.N, 2:ceil(par.S/2)+1];
        case 600
            par.N = 1024*par.sampling_rate_multiplier;
            par.submap = [par.N-par.S/2+1:par.N, 2:ceil(par.S/2)+1];
        case 900
            par.N = 1536*par.sampling_rate_multiplier;
            par.submap = [par.N-par.S/2+1:par.N, 2:ceil(par.S/2)+1];
        case 1200
            par.N = 2048*par.sampling_rate_multiplier;
            par.submap = [par.N-par.S/2+1:par.N, 2:ceil(par.S/2)+1];
    end
    
    % sampling rate
    par.samplingrate = 15e3*par.N;

    % precompute bit labels
    par.card = length(par.symbols); % cardinality
    par.bps = log2(par.card); % number of bits per symbol
    par.bits = de2bi(0:par.card-1,par.bps,'left-msb'); % symbols-to-bits
    
    % -- code parameters 
    if par.code.rate < 1
        
        addpath('BCJR'); % add path to decoder
        
        par.code.constraintlength = 7; % constraint Length
        par.code.generator = [133  171]; % generator polynomial      
        par.code.trellis = poly2trellis(par.code.constraintlength, par.code.generator); % trellis
        par.code.n = par.nrof_ofdm_symbols*par.bps*par.S; % length of code word
        par.code.m = floor(par.code.n*par.code.rate); % length of message
        for u=1:par.U
            par.code.interleaverperm(u,:) = randperm(par.code.n); % random interleaver
        end

        % puncturing pattern
        switch (par.code.rate) 
            case 1/2
              par.code.puncturing.period = 1;
              par.code.puncturing.pattern = 1;
            case 3/4
              par.code.puncturing.period = 18;
              par.code.puncturing.pattern = [1 2 3 6 7 8 9 12 13 14 15 18];
            case 5/6
              par.code.puncturing.period = 10;
              par.code.puncturing.pattern = [1 2 3 6 7 10];
            otherwise  
              error('coding rate not supported')    
        end
        
        % extend puncturing pattern 
        nrofcodedbits = par.code.n; % bits per terminal
        patlen = length(par.code.puncturing.pattern); 
        for i=1:ceil(nrofcodedbits/patlen)       
            par.code.puncturing.index(patlen*(i-1)+1:patlen*i) = ...
              (i-1)*par.code.puncturing.period+par.code.puncturing.pattern;   
        end
        par.code.puncturing.index(:,nrofcodedbits+1:end) = [];
        
    end
    
    % quantizer parameters
    par.clip_prb = 1e-3; % clipping probability (needs to be adjusted with the number of levels)
    par.clip_lvl = sqrt(par.S/par.N)/sqrt(2*par.B) * qfuncinv(par.clip_prb/2); % clipping level
    par.lsb = 2*par.clip_lvl/par.L; % least significant bit    
    par.labels = par.lsb *((0:par.L-1) - (par.L-1)/2); % uniform quantization labels
    par.thresholds = [-10^100, bsxfun(@minus, par.labels(:,2:end), par.lsb/2), 10^100];	% uniform quantization thresholds
    par.alpha = sqrt(2*par.B/(par.S/par.N)*sum(par.labels.^2.* ... 
        (normcdf(par.thresholds(2:end)*sqrt(2*par.B)/sqrt(par.S/par.N)) ...
        -normcdf(par.thresholds(1:end-1)*sqrt(2*par.B)/sqrt(par.S/par.N)))))^-1; % normalization constant (approximate, but accurate)
    par.labels = par.alpha*par.labels; % normalize quantization labels
    
    % clipping and quantization
    par.clipper = @(x) max(min(x,par.clip_lvl-par.lsb/1e5),-(par.clip_lvl-par.lsb/1e5)); % clipper
    if mod(par.L,2) == 0
        par.quantizer = @(x) par.alpha * (par.lsb*floor(par.clipper(x)/par.lsb) + par.lsb/2); % midrise quantizer (without clipping)
    else
        par.quantizer = @(x) par.alpha * par.lsb*floor(par.clipper(x)/par.lsb + 1/2); % midtread quantizer (without clipping)
    end
    par.quantizer = @(x) par.quantizer(par.clipper(real(x))) + 1i*par.quantizer(par.clipper(imag(x))); % quantizer
    
    % initialize result arrays
    [res.BER_uncoded, res.BER_coded] = deal(zeros(length(par.precoder),length(par.SNRdB_list)));
    [res.BER_uncoded_approx, res.sumrate_approx] = deal(zeros(length(par.precoder),length(par.SNRdB_list)));
    [res.TxAvgPower, res.RxAvgPower, res.TxMaxPower, res.RxMaxPower] = deal(zeros(length(par.precoder),1));
    pow_xf = zeros(length(par.precoder), par.N); 
    pow_xf_emp = zeros(length(par.precoder), par.N);
    pow_yf = zeros(length(par.precoder), par.N);
    pow_yf_emp = zeros(length(par.precoder), par.N);

    % save detected symbols for later viewing (if not too many of them)
    if par.nrof_ofdm_symbols * par.nrof_channels * par.S <= 1e5
        shat_list = nan(par.U, par.S*par.nrof_ofdm_symbols,par.nrof_channels,length(par.precoder));
    end

    % track simulation time
    time_elapsed = 0; tic;
    
    % -- start of simulation
    
    fprintf(' running numerical simulation. \n');

    for cc = 1:par.nrof_channels
        
        % time-domain channel matrix
        Ht = sqrt(0.5/par.T)*(randn(par.U,par.B,par.T) + 1i*randn(par.U,par.B,par.T));
        
        % channel estimation error
        if par.relerr > 0
            Ht_est = sqrt(1-par.relerr)*Ht + sqrt(par.relerr/2)*(randn(par.U,par.B,par.T) + 1i*randn(par.U,par.B,par.T));
        else
            Ht_est = Ht;
        end
        
        % account for sampling rate multiplier
        if par.T > 1 
            Ht = permute(upsample(permute(Ht,[3,2,1]),round(par.sampling_rate_multiplier)),[3,2,1]);
            if par.relerr > 0
                Ht_est = permute(upsample(permute(Ht_est,[3,2,1]),round(par.sampling_rate_multiplier)),[3,2,1]);
            else
                Ht_est = Ht;
            end
        end
        
        % frequency-domain channel matrix
        if par.T > 1
            Hf = fft(Ht,par.N,3);
            Hf_est = fft(Ht_est,par.N,3); 
        else
            Hf = nan(par.U,par.B,par.N);
            Hf_est = nan(par.U,par.B,par.N);
            for n = 1:par.N
                Hf(:,:,n) = Ht;
                Hf_est(:,:,n) = Ht_est;
            end
        end

        % generate noise vector
        nt = sqrt(0.5)*(randn(par.U,par.N,par.nrof_ofdm_symbols)+1i*randn(par.U,par.N,par.nrof_ofdm_symbols)); % time domain
        nf = sqrt(1/par.N)*fft(nt,par.N,2); % frequency domain
        
        % information bits and coded bits
        if par.code.rate < 1
            par.code.padbits = ceil(par.U*par.code.n*par.code.rate)-floor(par.U*par.code.n*par.code.rate);
            b_data = round(rand(par.U,ceil(par.code.n*par.code.rate)));
            b_data(1:par.U,size(b_data,2)+1-(par.code.constraintlength-1)-par.code.padbits:size(b_data,2)) = 0;
            b_encoded = nan(par.U, par.code.n); % encoded bits
            for u = 1:par.U
                b_trellis = convenc(b_data(u,:),par.code.trellis); % encoding
                b_punctured = b_trellis(par.code.puncturing.index); % puncturing
                b_encoded(u,:) = b_punctured(par.code.interleaverperm(u,:)); % interleaving
            end
        else
            b_encoded = randi([0 1], par.U, par.bps*par.S*par.nrof_ofdm_symbols); % data bits
        end
        
        % map bits to subcarriers
        s = zeros(par.U,par.N,par.nrof_ofdm_symbols);
        for u = 1:par.U
            idx = reshape(bi2de(reshape(b_encoded(u,:),par.bps,[])','left-msb')+1,par.S,[]); % index
            s(u,par.submap,:) = par.symbols(idx); 
        end
         
        % algorithm loop
        for pp = 1:length(par.precoder) 

            % ZF precoding
            zf = zeros(par.B,par.N,par.nrof_ofdm_symbols); % precoded frequency-domain vector
            Pf = zeros(par.B,par.U,par.N); % precoding matrix
            [zf(:,par.submap,:), beta, Pf(:,:,par.submap)] = ZF(par, s(:,par.submap,:), Hf_est(:,:,par.submap));
            zt = sqrt(par.N)*ifft(zf, par.N, 2); % transform to time domain
            
            % quantization
            if ismember(par.precoder{pp}, 'ZFQ')
                xt = par.quantizer(zt); % quantize time-domain signal
                bussgang_gain = par.alpha*par.lsb*sqrt(par.B/(par.S/par.N)/pi)* ...
                    sum(exp(-par.B/(par.S/par.N)*par.lsb^2*((1:par.L-1)-par.L/2).^2)); % Bussgang gain 
            else
                xt = zt; % do nothing
                bussgang_gain = 1; % Bussgang gain 
            end
            
            if ~strcmpi(par.approximation, 'none')

                % input covariance matrix
                Czf = zeros(par.B,par.B,par.N); 
                for k = par.submap
                    Czf(:,:,k) = Pf(:,:,k)*Pf(:,:,k)'; 
                end
                Czt = fft(Czf, par.N, 3)/par.N;

                % linear gain matrix and output covariance matrix
                if ismember(par.precoder{pp}, 'ZFQ')
                   [G, Cxt, Cdt] = bussgang_decomposition(par, Czt);
                else
                    G = eye(par.B);
                    Cxt = Czt;
                    Cdt = zeros(size(Czt));
                end
                Cxf = par.N*ifft(Cxt, par.N, 3);
                Cdf = par.N*ifft(Cdt, par.N, 3);

                % signal-to-interference-noise-and-distortion ratio (SINDR)
                sindr = nan(par.U,par.S,length(par.SNRdB_list));
                for k = 1:par.S 
                    S = repmat(diag(abs(Hf(:,:,par.submap(k))*G*Pf(:,:,par.submap(k))).^2),1,length(par.SNRdB_list)); % signal power
                    I = repmat(sum(abs(Hf(:,:,par.submap(k))*G*Pf(:,:,par.submap(k))).^2,2),1,length(par.SNRdB_list)) - S; % interference power
                    D = repmat(real(diag(Hf(:,:,par.submap(k))*Cdf(:,:,par.submap(k))*Hf(:,:,par.submap(k))')),1,length(par.SNRdB_list)); % distortion power
                    N = repmat(10.^(-par.SNRdB_list/10),par.U,1); % noise power
                    sindr(:,k,:) = S ./ (I + D + N);
                end

                % compute uncoded BER and achievable rate (analytical)
                res.BER_uncoded_approx(pp,:) = res.BER_uncoded_approx(pp,:) + shiftdim(mean(mean(qfunc(sqrt(sindr)),1),2),1)/par.nrof_channels;
                res.sumrate_approx(pp,:) = res.sumrate_approx(pp,:) + shiftdim(sum(sum(log2(1 + sindr),1),2),1)/par.nrof_channels/par.S;

                % power spectral density (analytical)
                for n = 1:par.N
                    pow_xf(pp,n) = pow_xf(pp,n) + mean(diag(Cxf(:,:,n)))/par.nrof_channels;
                    pow_yf(pp,n) = pow_yf(pp,n) + mean(diag(Hf(:,:,n)*Cxf(:,:,n)*Hf(:,:,n)'))/par.nrof_channels;
                end
                
            end
            
            for ss = 1:length(par.SNRdB_list) % SNR loop

                % noise variance 
                N0 = 10.^(-par.SNRdB_list(ss)/10); 
                
                % convert transmit signal to frequency domain
                xf = 1/sqrt(par.N)*fft(xt, par.N, 2);

                % transmit signal over wireless channel
                Hxf = nan(par.U,par.N,par.nrof_ofdm_symbols);
                for n = 1:par.N
                    Hxf(:,n,:) = Hf(:,:,n)*squeeze(xf(:,n,:));
                end
                
                % add noise
                yf = Hxf + sqrt(N0)*nf; 

                % extract transmitted/received power
                res.TxMaxPower(pp) = max(res.TxMaxPower(pp), max(squeeze(sum(sum(abs(xf).^2,1),2)))/par.S);
                res.TxAvgPower(pp) = res.TxAvgPower(pp) + sum(sum(squeeze(sum(abs(xf).^2))))/par.S/par.nrof_ofdm_symbols/par.nrof_channels/length(par.SNRdB_list);

                % estimated symbols (remove guard carriers and compensate for gain loss)
                shat = beta * yf(:,par.submap,:) / bussgang_gain; % = yf(:,par.submap,:) ./ sqrt(mean(abs(yf(:,par.submap,:)).^2,2) - N0);
                 
                % symbol detection and decoding
                if par.code.rate == 1 % nearest-neighbor
                    b_detected = nan(size(b_encoded));
                    for u = 1:par.U
                        [~,idxhat] = min(abs(reshape(shat(u,:,:),[],1)*ones(1,length(par.symbols))-ones(par.S*par.nrof_ofdm_symbols,1)*par.symbols).^2,[],2);
                        b_detected(u,:) = reshape(par.bits(idxhat,:)',1,[]);
                    end
                else % log-max BCJR
                    llr = detector(par,reshape(shat,par.U,par.S*par.nrof_ofdm_symbols),N0);
                    b_detected = (llr>0);
                    [~,b_decoded] = decoder(par,llr);
                end
                
                % compute error metrics
                res.BER_uncoded(pp,ss) = res.BER_uncoded(pp,ss) + sum(sum(b_encoded~=b_detected))/par.U/par.bps/par.S/par.nrof_ofdm_symbols/par.nrof_channels;
                if par.code.rate < 1
                    res.BER_coded(pp,ss) = res.BER_coded(pp,ss) + sum(sum(b_data~=b_decoded))/par.U/par.code.m/par.nrof_channels;
                end
                
            end % end of SNR loop
            
            % save data for later viewing
            if par.nrof_ofdm_symbols * par.nrof_channels * par.S <= 1e5
                shat_list(:,:,cc,pp) = reshape(shat,par.U,[]);
            end
            
            % power spectral density (simulated)
            pow_xf_emp(pp,:) = pow_xf_emp(pp,:) + mean(mean(abs(xf).^2,3),1)/par.nrof_channels;
            pow_yf_emp(pp,:) = pow_yf_emp(pp,:) + mean(mean(abs(Hxf).^2,3),1)/par.nrof_channels;

        end % end of algorithm loop

         % keep track of simulation time
        if toc>10
            time=toc;
            time_elapsed = time_elapsed + time;
            fprintf('\t estimated remaining simulation time: %3.0f min. \n',time_elapsed*(par.nrof_channels/cc-1)/60);
            tic
        end

    end % end of channels loop

    fprintf('\t numerical simulation finished after %.2f seconds. \n', time_elapsed);

    % -- end of simulation

    if par.plot
        
       close all;
        
        % plot tranmitted and received power spectrum
        ylim_min = -50; ylim_max = ceil(10*log10(max([max(abs(pow_xf_emp(pp,:))), max(abs(pow_yf_emp(pp,:)))])))+10; % limits
        fig_txspec = figure; set(fig_txspec,'name','Tx Spectrum','numbertitle','off'); % transmitted spectrum
        for pp = 1:length(par.precoder)
            subplot(md,nd,pp); hold all;
            plot(15e-3*(-par.N/2:par.N/2-1),10*log10(eps+fftshift(pow_xf_emp(pp,:)+eps)),'-x','color',marker_color(pp,:));
            if ~strcmpi(par.approximation, 'none')
                plot(15e-3*(-par.N/2:par.N/2-1),10*log10(fftshift(abs(pow_xf(pp,:))+eps)),'k-');
                set(legend('Simulated', 'Analytical'), 'fontsize',14, 'location', 'southeast');
            else
                set(legend('Simulated'), 'fontsize',14, 'location', 'southeast');
            end
            xlim(15e-3*[-par.N/2, par.N/2-1]); ylim([ylim_min, ylim_max]);
            xlabel('Frequency [MHz]','fontsize',14); ylabel('PSD [dB]','fontsize',14); box on; grid on;
            title(par.precoder{pp},'fontsize',12);
        end
        fig_rxspec = figure; set(fig_rxspec,'name','Rx Spectrum','numbertitle','off'); % received spectrum
        for pp = 1:length(par.precoder)
            subplot(md,nd,pp); hold all;
            plot(15e-3*(-par.N/2:par.N/2-1),10*log10(eps+fftshift(pow_yf_emp(pp,:)+eps)),'-x','color',marker_color(pp,:));
            if ~strcmpi(par.approximation, 'none')
                plot(15e-3*(-par.N/2:par.N/2-1),10*log10(fftshift(abs(pow_yf(pp,:))+eps)),'k-');
                set(legend('Simulated', 'Analytical'), 'fontsize',14, 'location', 'southeast');
            else
                set(legend('Simulated'), 'fontsize',14, 'location', 'southeast');
            end
            xlim(15e-3*[-par.N/2, par.N/2-1]); ylim([ylim_min, ylim_max]);
            xlabel('Frequency [MHz]','fontsize',14); ylabel('PSD [dB]','fontsize',14); box on; grid on;
            title(par.precoder{pp},'fontsize',12);
        end
        
        % plot constellation
        if par.nrof_ofdm_symbols * par.nrof_channels * par.S <= 1e5
            fig_const = figure; set(fig_const,'name','Const.','numbertitle','off');
            for pp = 1:length(par.precoder)
                subplot(md,nd,pp); hold all;
                plot(reshape(shat_list(:,:,:,pp),1,[]),'*', 'color', marker_color(pp,:),'markersize',7);
                plot(par.symbols, 'ko','MarkerSize',7);
                axis(max(reshape(abs(shat_list(:,:,:,pp)),1,[]))*[-1 1 -1 1]); 
                axis square; box on;
                title(par.precoder{pp},'fontsize',12);
                xlabel(['P_{avg}= ',num2str(pow2db(res.TxAvgPower(pp)),'%0.2f'),' dB  and  P_{max}= ',num2str(pow2db(res.TxMaxPower(pp)),'%0.2f'),' dB'],'fontsize',12);
            end

        end

        % plot uncoded BER
        fig_uncodedber = figure; set(fig_uncodedber,'name','Uncoded BER','numbertitle','off');
        for pp=1:length(par.precoder) % simulated BER
            semilogy(par.SNRdB_list,res.BER_uncoded(pp,:),marker_style{pp},'color',marker_color(pp,:),'LineWidth',2); hold on;
        end
        if strcmpi(par.mod, 'QPSK') && ~strcmpi(par.approximation, 'none')
            for pp=1:length(par.precoder)
                semilogy(par.SNRdB_list,res.BER_uncoded_approx(pp,:),'-','color',marker_color(pp,:),'LineWidth',2); hold on;
            end
        end
        grid on; box on;
        xlabel('SNR [dB]','FontSize',12)
        ylabel('uncoded BER','FontSize',12);
        if length(par.SNRdB_list) > 1
            axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1]);
        end
        legend(par.precoder,'FontSize',12,'location','southwest')
        set(gca,'FontSize',12);
        
        % plot coded BER
        if par.code.rate < 1
            fig_codedber = figure; set(fig_codedber,'name','Coded BER','numbertitle','off');
            for pp=1:length(par.precoder) % simulated BER
                semilogy(par.SNRdB_list,res.BER_coded(pp,:),marker_style{pp},'color',marker_color(pp,:),'LineWidth',2); hold on;
            end
            grid on; box on;
            xlabel('SNR [dB]','FontSize',12)
            ylabel('coded BER','FontSize',12);
            if length(par.SNRdB_list) > 1
                axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-6 1]);
            end
            legend(par.precoder,'FontSize',12,'location','southwest')
            set(gca,'FontSize',12);
        end

    end
        
    % save final results
    if par.save
        save(par.simName,'par','res');
    end
    
    if nargin == 0
       keyboard; 
    end
    
end

function [xf, beta, Pf] = ZF(par, s, Hf)
    
    % initialize vectors
    xf = zeros(par.B,par.S,par.nrof_ofdm_symbols); 
    Pf = nan(par.B,par.U,par.S);
    beta = 0;
    
    % precoding
    for k = 1:par.S
        Pf(:,:,k) = Hf(:,:,k)'/(Hf(:,:,k)*Hf(:,:,k)'); % precoding matrix
        xf(:,k,:) = Pf(:,:,k)*reshape(s(:,k,:),par.U,[]); % precoded vector
        beta = beta + trace(Pf(:,:,k)*Pf(:,:,k)')/par.S; % precoding factor (squared) 
    end
    beta = sqrt(beta); 
    
    % scale output
    xf = xf/beta;
    Pf = Pf/beta; 
    
end

function [G, Cxt, Cdt] = bussgang_decomposition(par, Czt)
    
    % diagonal of the zero-lag covariance matrix
    Dzt = real(diag(diag(Czt(:,:,1)))); 

    % Bussgang matrix
    if strcmpi(par.approximation, 'rounding') || strcmpi(par.approximation, 'diagonal')
        G = zeros(par.B,par.B); 
        for i = 1:par.L-1
            G = G + par.alpha*par.lsb/sqrt(pi)*Dzt^-.5 * diag(exp(-par.lsb^2*(i-par.L/2)^2*diag(Dzt^-1)));
        end
    end
    
    switch par.approximation
        
        case 'diagonal'
                
            % output covariance matrix
            Cxt = zeros(par.B,par.B,par.N);
            for i = 1:par.L
                Cxt(:,:,1) = Cxt(:,:,1) + diag(2 * par.labels(i)^2 * (qfunc(sqrt(2)*par.thresholds(i)./sqrt(diag(Dzt))) - qfunc(sqrt(2)*par.thresholds(i+1)./sqrt(diag(Dzt)))));
            end
            Cxt(:,:,1) = Cxt(:,:,1) + G*(Czt(:,:,1) - Dzt)*G';
            for n = 2:par.N
                Cxt(:,:,n) = G*Czt(:,:,n)*G';
            end
                
        case 'rounding'
            
            % output covariance matrix
            Cxt = zeros(size(Czt));
            if par.L == 2
                
                % arcsine law
                for n = 1:par.N
                    Cxt(:,:,n) =  2/pi * par.S/(par.N*par.B) * (asin(Dzt^-.5*real(Czt(:,:,n))*Dzt^-.5) ...
                        + 1i*asin(Dzt^-.5*imag(Czt(:,:,n))*Dzt^-.5)); 
                end  
                
            elseif par.L > 2
                
                Cet = zeros(size(Czt));
                
                for n = 1:par.N

                    % limit the infinite sum
                    lmax = 30; lvec = combvec(1:lmax,1:lmax);

                    % quantization error covariance matrix
                    for l = 1:lmax^2
                        cv = 2*pi^2/par.lsb^2*(lvec(1,l).^2*real(Dzt*ones(par.B,par.B))/2 + lvec(2,l).^2*real(ones(par.B,par.B)*Dzt)/2);
                        cr = 4*pi^2/par.lsb^2*lvec(1,l).*lvec(2,l)*real(Czt(:,:,n))/2;
                        ci = 4*pi^2/par.lsb^2*lvec(1,l).*lvec(2,l)*imag(Czt(:,:,n))/2;
                        Cet(:,:,n) = Cet(:,:,n) + 2* par.lsb^2/(2*pi^2)  .* ...
                            cos(pi*par.L).^(lvec(1,l)+lvec(2,l)) ./ (lvec(1,l).*lvec(2,l)) .* ...
                            (exp(-cv+cr) - exp(-cv-cr) + 1i*exp(-cv+ci) - 1i*exp(-cv-ci));
                    end

                end
                 
                 % output covariance matrix
                 for n = 1:par.N
                    Cxt(:,:,n) = par.alpha^2*(Cet(:,:,n) + Czt(:,:,n)*(G/par.alpha-eye(par.B)) + (G/par.alpha-eye(par.B))*Czt(:,:,n) + Czt(:,:,n));
                 end

            end
            
    end
    
    % distortion covariance matrix
    Cdt = zeros(size(Czt));
    for n = 1:par.N
        Cdt(:,:,n) = Cxt(:,:,n) - G*Czt(:,:,n)*G';
    end

end

% -- auxilliary functions

% soft-output data detection using the max-log approximation.
function llr = detector(par,shat,N0)

    % initialization 
    llr = nan(par.U, par.code.n);
    bin_array = sign(par.bits-.5);
   
    % perform detection for each terminal
    for u = 1:par.U
        for k = 1:par.S*par.nrof_ofdm_symbols
        
            % compute distance metric
            metrics = abs(shat(u,k)-par.symbols).^2/N0;

            % compute max-log LLRs
            for b = 1:par.bps
                
                pmin = inf; 
                nmin = inf;
                
                for z = 1:2^par.bps         
                    if bin_array(z,b)==1
                        pmin = min(pmin,metrics(z));
                    else
                        nmin = min(nmin,metrics(z));          
                    end                                
                end
                
                % log-likelihood ratio
                llr(u,(k-1)*par.bps+b) = nmin-pmin;
                
            end 

        end

    end  
  
end

% deinterleaving, depuncturing, and decoding using the BCJR algorithm
function [llr_out,bhat] = decoder(par,llr_in)

    % initialization
    llr_temp = nan(1,size(llr_in,2));
    llr_out = nan(par.U, 2*par.code.m);
    bhat = nan(par.U, par.code.m);

    for u = 1:par.U

        % deinterleaving
        llr_temp(1,par.code.interleaverperm(u,:)) = llr_in(u,:);   

        % depuncturing
        llr_punctured = zeros(1,2*length(llr_temp(1,:))*par.code.rate);    
        llr_punctured(par.code.puncturing.index) = llr_temp;
        
        % decoding
        [llr_out(u,:),bhat(u,:)] = BCJR(par.code.trellis,llr_punctured); 

    end
  
end






