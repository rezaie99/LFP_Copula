data = load("POm_data.mat");
for N = 1:9
    disp('Load Data 1')
    downSample = 500;

    LFP = data.POm_data{N}.LFP;
    LFP = LFP(1:downSample*floor(size(LFP)/downSample));
    LFP = reshape(LFP,floor(size(LFP,1)/downSample),downSample);
    LFP = sum(LFP,2)/downSample;
    LFP = LFP(1:100*floor(size(LFP)/100));

    spike_serie = (data.POm_data{N}.filteredResponse>500).*1;
    spike_serie = spike_serie(1:downSample*floor(size(spike_serie)/downSample));
    spike_serie = reshape(spike_serie,floor(size(spike_serie,1)/downSample),downSample);
    spike_serie = sum(spike_serie,2);
    spike_serie = spike_serie(1:100*floor(size(spike_serie,1)/100));

    LFP = reshape(LFP, [floor(size(LFP,1)/100),100]);
    spike_serie = reshape(spike_serie, [floor(size(spike_serie,1)/100),100]);
    fprintf('Analysis for N= %d\n',N)
    [gcm,gcer] = demo(LFP,spike_serie,N);
    display(gcm);
    display(gcer);
end

data = load("POm_L5_data.mat");
for N = 1:14
    disp('Load Data 2')
    downSample = 500;

    LFP = data.POm_L5_data{N}.LFP;
    LFP = LFP(1:downSample*floor(size(LFP)/downSample));
    LFP = reshape(LFP,floor(size(LFP,1)/downSample),downSample);
    LFP = sum(LFP,2)/downSample;
    LFP = LFP(1:100*floor(size(LFP)/100));

    spike_serie = (data.POm_L5_data{N}.filteredResponse>500).*1;
    spike_serie = spike_serie(1:downSample*floor(size(spike_serie)/downSample));
    spike_serie = reshape(spike_serie,floor(size(spike_serie,1)/downSample),downSample);
    spike_serie = sum(spike_serie,2);
    spike_serie = spike_serie(1:100*floor(size(spike_serie,1)/100));

    LFP = reshape(LFP, [floor(size(LFP,1)/100),100]);
    spike_serie = reshape(spike_serie, [floor(size(spike_serie,1)/100),100]);
    fprintf('Analysis for N= %d\n',N)
    [gcm,gcer] = demo(LFP,spike_serie,N);
    display(gcm);
    display(gcer);
end