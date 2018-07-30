N=1
%data = load("T13_17_28_00_pom.mat");
%data = load("T16_16_14_00_pom.mat");
%data = load("T18_15_43_00_pom.mat");
%data = load("T10_18_01_00_vpm.mat");
%data = load("T14_16_04_00_vpm.mat");
data = load("T14_16_20_00_vpm.mat");
disp('T14_16_20_00_vpm.mat')
downSampleSpike = 20;
downSampleLFP = 1;

LFP = data.LFP;
LFP = LFP(1:downSampleLFP*floor(size(LFP,2)/downSampleLFP));
LFP = reshape(LFP,floor(size(LFP,2)/downSampleLFP),downSampleLFP);
LFP = sum(LFP,2)/downSampleLFP;
LFP = LFP(1:100*floor(size(LFP)/100));

spike_serie = (data.filteredResponse).*0;
spike_serie(data.spikes) = 1;
spike_serie = spike_serie(1:downSampleSpike*floor(size(spike_serie,2)/downSampleSpike));
spike_serie = reshape(spike_serie,floor(size(spike_serie,2)/downSampleSpike),downSampleSpike);
spike_serie = sum(spike_serie,2);
spike_serie = spike_serie(1:100*floor(size(spike_serie,1)/100));

LFP = reshape(LFP, [floor(size(LFP,1)/100),100]);
spike_serie = reshape(spike_serie, [floor(size(spike_serie,1)/100),100]);

minSize = min(size(LFP,1),size(spike_serie,1));
LFP = LFP(1:minSize,:);
spike_serie = spike_serie(1:minSize,:);
fprintf('Analysis for N= %d\n',N)
[gcm,gcer] = demo(LFP,spike_serie,N);
display(gcm);
display(gcer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whisker_serie = (data.Triggers.whisker).*1;

whisker_serie = whisker_serie(1:downSampleSpike*floor(size(whisker_serie,2)/downSampleSpike));
whisker_serie = reshape(whisker_serie,downSampleSpike,floor(size(whisker_serie,2)/downSampleSpike));
whisker_serie = sum(whisker_serie,1);
whisker_serie = whisker_serie(1:100*floor(size(whisker_serie,2)/100));

whisker_serie = reshape(whisker_serie, [floor(size(whisker_serie,2)/100),100]);
whisker_serie = whisker_serie(1:minSize,:);
fprintf('Analysis for N= %d\n',N)
[gcm,gcer] = demo(LFP,whisker_serie,N);
display(gcm);
display(gcer);