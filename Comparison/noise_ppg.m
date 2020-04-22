function [pulse_noise] = noise_ppg(sig,noiseMode,snrdb,freq)
N=length(sig);
ts=1/freq; 
t = ts*((0:N-1)); % create time vector

switch noiseMode
    case 'white'
        pulse_noise = add_awgn_noise(sig,snrdb);
    case 'pink'
        Ps = 10*log10(std(sig).^2);  % signal power     dBV^2
        Pn = Ps - snrdb;           % noise power       dBV^2
        Pn = 10^(Pn/10);         
        sigma = sqrt(Pn);        
        n = sigma*pinknoise(1,N);
        pulse_noise = sig + n;
    case 'mov1' % movement artifact 1
        artefakt1=1.5*t; 
        pulse_noise=artefakt1(1:N)+sig;
    case 'mov2' % movement artifact 2
        artefakt2=1.5*t;
        sig= 0.5*sig;
        pulse_noise=sig+artefakt2(1:N);
    otherwise
        errordlg('Unknown noise mode selected','Input Error','modal');
    return;
end
end

