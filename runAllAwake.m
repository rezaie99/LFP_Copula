fileList = dir('awake/POm/');

downSampleSpike = 20;
downSampleLFP = 1;

for i = 1:length(fileList)
    if length(fileList(i).name)>4
        if all(fileList(i).name(end-3:end) == '.mat')
            data = load([fileList(i).folder '/' fileList(i).name]);
            display(fileList(i).name)

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

            whisker = data.Triggers.whisker*1;
            whisker = whisker(1:downSampleSpike*floor(size(whisker,2)/downSampleSpike));
            whisker = reshape(whisker,floor(size(whisker,2)/downSampleSpike),downSampleSpike);
            whisker = sum(whisker,2);
            whisker = whisker(1:100*floor(size(spike_serie,1)/100));
            
            
            LFP = reshape(LFP, [floor(size(LFP,1)/100),100]);
            spike_serie = reshape(spike_serie, [floor(size(spike_serie,1)/100),100]);
            whisker = reshape(whisker, [floor(size(whisker,1)/100),100]);

            minSize = min(size(LFP,1),size(spike_serie,1));
            LFP = LFP(1:minSize,:);
            spike_serie = spike_serie(1:minSize,:);
            fprintf('Analysis for N= %d\n',i)
            try
                mask = whisker>-10;
                [gcm,gcer,gc12_all,gc21_all] = demo(LFP,spike_serie,1,mask);
                display(gcm);
                display(gcer);
                fileList(i).gcm12 = gcm(1);
                fileList(i).gcm21 = gcm(2);
                fileList(i).gcmall12 = gc12_all;
                fileList(i).gcmall21 = gc21_all;
                fileList(i).gcer12 = gcer(1);
                fileList(i).gcer21 = gcer(2);
            catch
                continue
            end
        end
    end
end

fileListPOm = fileList;

fileList = dir('awake/VPm/');

downSampleSpike = 20;
downSampleLFP = 1;

for i = 1:length(fileList)
    if length(fileList(i).name)>4
        if all(fileList(i).name(end-3:end) == '.mat')
            data = load([fileList(i).folder '/' fileList(i).name]);
            display(fileList(i).name)

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
            
            whisker = data.Triggers.whisker*1;
            whisker = whisker(1:downSampleSpike*floor(size(whisker,2)/downSampleSpike));
            whisker = reshape(whisker,floor(size(whisker,2)/downSampleSpike),downSampleSpike);
            whisker = sum(whisker,2);
            whisker = whisker(1:100*floor(size(spike_serie,1)/100));


            LFP = reshape(LFP, [floor(size(LFP,1)/100),100]);
            spike_serie = reshape(spike_serie, [floor(size(spike_serie,1)/100),100]);
            whisker = reshape(whisker, [floor(size(whisker,1)/100),100]);

            minSize = min(size(LFP,1),size(spike_serie,1));
            LFP = LFP(1:minSize,:);
            spike_serie = spike_serie(1:minSize,:);
            fprintf('Analysis for N= %d\n',i)
            try
                mask = whisker>-10;
                [gcm,gcer,gc12_all,gc21_all] = demo(LFP,spike_serie,1,mask);
                display(gcm);
                display(gcer);
                fileList(i).gcm12 = gcm(1);
                fileList(i).gcm21 = gcm(2);
                fileList(i).gcmall12 = gc12_all;
                fileList(i).gcmall21 = gc21_all;
                fileList(i).gcer12 = gcer(1);
                fileList(i).gcer21 = gcer(2);
            catch
                continue
            end
        end
    end
end
fileListVPm = fileList;