function [nextBegin, nextorigin, nextsr] = PhaseVocoderMATLABExample ( start, now , origin, sr )

% This implements a phase vocoder for time-stretching/compression
% Put the file name of the input file in "infile"
%     the file to be created is "outfile"
% The amount of time-stretch is determined by the ratio of the hopin
% and hopout variables. For example, hopin=242 and hopout=161.3333
% (integers are not required) increases the tempo by 
% hopin/hopout = 1.5. To slow down a comparable amount,
% choose hopin = 161.3333, hopout = 242.
% 5/2005 Bill Sethares

    time=0;                                    % total time to process
    hopin=121;                                 % hop length for input 
    hopout=121;                                % hop length for output

    all2pi=2*pi*(0:100);                       % all multiples of 2 pi (used in PV-style freq search) 
    max_peak=50;                               % parameters for peak finding: number of peaks
    eps_peak=0.005;                            % height of peaks
    nfft=2^12; nfft2=nfft/2;                   % fft length 
    win=hanning(nfft)';                        % windows and windowing variables   
    %[y,sr]=audioread(infile);                    % read song file
    
    
    SPB = 60/(109*4);
    beginSecond = start / 1000; % Define where you want to start playing
    nowSecond = now / 1000;
    beatDiff =  round( (nowSecond - beginSecond) / SPB ) ;
    beginIndex = floor(beginSecond*sr); % Find beginning index of where to sample
    nowIndex = floor(nowSecond*sr);
    [y, ~] = audioread('input.mp3', [beginIndex nowIndex]); % Extract out the signal from the starting point to the end of the file
   % sound(y, sr); % Play the sound
    hopin = hopin * (nowSecond - beginSecond)/ (SPB * beatDiff);



    siz= size(y);           % length of song in samples
    stmo=min(siz); leny=max(siz);              % stereo (stmo=2) or mono (stmo=1)
    if siz(2)<siz(1), y=y'; end
    if time==0, time = 100000; end
    tt=zeros(stmo,ceil(hopout/hopin)*min(leny,time*sr)); % place for output 
    lenseg=floor((min(leny,time*sr)-nfft)/hopin); % number of nfft segments to process

    ssf=sr*(0:nfft2)/nfft;                     % frequency vector
    phold=zeros(stmo,nfft2+1); phadvance=zeros(stmo,nfft2+1);
    outbeat=zeros(stmo,nfft); pold1=[]; pold2=[];
    dtin=hopin/sr;                             % time advances dt per hop for input
    dtout=hopout/sr;                           % time advances dt per hop for output
    for k=1:lenseg-1                           % main loop - process each beat separately
        if k/1000==round(k/1000), disp(k), end % for display so I know where we are
        indin=round(((k-1)*hopin+1):((k-1)*hopin+nfft));
        for sk=1:stmo                          % do L/R channels separately
            s=win.*y(sk,indin);                % get this frame and take FFT
            ffts=fft(s);
            mag=abs(ffts(1:nfft2+1)); 
            ph=angle(ffts(1:nfft2+1));

            % find peaks to define spectral mapping

            peaks=findPeaks4(mag, max_peak, eps_peak, ssf);
            [dummy,inds]=sort(mag(peaks(:,2)));
            peaksort=peaks(inds,:);
            pc=peaksort(:,2);
            
            bestf=zeros(size(pc));
            for tk=1:length(pc)                % estimate frequency using PV strategy  
                dtheta=(ph(pc(tk))-phold(sk,pc(tk)))+all2pi;
                fest=dtheta./(2*pi*dtin);      % see pvanalysis.m for same idea
                [er,indf]=min(abs(ssf(pc(tk))-fest));
                bestf(tk)=fest(indf);          % find best freq estimate for each row
            end

            % generate output mag and phase
 
            magout=mag; phout=ph;
            for tk=1:length(pc)
                fdes=bestf(tk);                           % reconstruct with original frequency
                freqind=(peaksort(tk,1):peaksort(tk,3));  % indices of the surrounding bins
                
                % specify magnitude and phase of each partial
                magout(freqind)=mag(freqind);
                phadvance(sk,peaksort(tk,2))=phadvance(sk,peaksort(tk,2))+2*pi*fdes*dtout;
                pizero=pi*ones(1,length(freqind));
                pcent=peaksort(tk,2)-peaksort(tk,1)+1;
                indpc=(2-mod(pcent,2)):2:length(freqind);
                pizero(indpc)=zeros(1,length(indpc));
                phout(freqind)=phadvance(sk,peaksort(tk,2))+pizero;
            end

            % reconstruct time signal (stretched or compressed)

            compl=magout.*exp(sqrt(-1)*phout);
            compl(nfft2+1)=ffts(nfft2+1);
            compl=[compl,fliplr(conj(compl(2:(nfft2))))];
            wave=real(ifft(compl));
            outbeat(sk,:)=wave;
            phold(sk,:)=ph; 
            
        end % end stereo
        indout=round(((k-1)*hopout+1):((k-1)*hopout+nfft));
        tt(:,indout)=tt(:,indout)+outbeat;
    end

    tmp = ~any(tt,1);
    tmp2 = ~any(tt,2);
    tt=0.8*tt/max(max(abs(tt)));
    

   if hopout / hopin >= 1
        tmp = ones( ceil(hopout/hopin)*min(leny,time*sr) + length(origin), stmo );
        tmp(:) = nan;
        tmp(beginIndex:ceil(hopout/hopin)*min(leny,time*sr)+beginIndex-1, :) = tt';
        locations = find(isnan(tmp));
        tmp(locations) = origin;
        origin = tmp;
        rowsOfZero = ~any(origin,2);
        rowsOfZero(1:beginIndex - 1) = 0;
        rowsOfZero(ceil(hopout/hopin)*min(leny,time*sr)+beginIndex : length(rowsOfZero)) = 0;
        origin( rowsOfZero, : ) = [];  %rows
    else
        origin(beginIndex : nowIndex, :) = tt';
        rowsOfZero = ~any(origin,2);
        rowsOfZero(1:beginIndex - 1) = 0;
        rowsOfZero(nowIndex + 1 : length(rowsOfZero)) = 0;
        origin( rowsOfZero, : ) = [];  %rows
    end


    %{  
    tt( ~any(tt,2), : ) = [];  %rows
    tt( :, ~any(tt,1) ) = []; 
    [rtt,ctt]=size(tt); if rtt==2, tt=tt'; end
    audiowrite(outfile, y, sr, 'BitsPerSample',32 );
    fclose('all');
    %}
    
   
   % origin( :, ~any(origin,1) ) = []; 
    [rtt,ctt]=size(origin); if rtt==2, origin=origin'; end
    %audiowrite(outfile, origin, sr, 'BitsPerSample',32 );
    %fclose('all');
    
    
    nextBegin = (beginSecond + (SPB * beatDiff) ) * 1000;
    nextorigin = origin;
    nextsr = sr;
end