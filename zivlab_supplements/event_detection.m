function [eventsMatrix,onsetWidthMatrix]=event_detection(traces,eventParams)
% This function performs event detection on the df/f0 traces

thresh=eventParams.eventTh;
tau=eventParams.tau;
% grade_thresh=eventParams.grade;
fs=eventParams.fs;
M=size(traces,1);
num_cells=size(traces,2);
dt_min=0.15;
dt_samp_min=round(dt_min*fs);

smooth_normalized_traces=traces./repmat(median(abs(traces)),M,1);

events=zeros(size(traces));
events(smooth_normalized_traces>=thresh)=smooth_normalized_traces(smooth_normalized_traces>=thresh);

events_logical=events>0;
diff_events_logical=diff(events_logical);

eventsMatrix=zeros(size(traces));
onsetWidthMatrix=zeros(size(traces));
relativeAmpMatrix=zeros(size(traces));
multiEventsMatrix=zeros(size(traces));
eventlengthMatrix=zeros(size(traces));
areaUnderCurve=zeros(size(traces));
numEventsMatrix=zeros(size(traces));
h = waitbar(0,'Performing event detection');
for n=1:num_cells
    waitbar((n-1)/num_cells,h,['Event detection - detecting events in cell number ' num2str(n) '/' num2str(num_cells)])
    events_beg=find(diff_events_logical(:,n)==1);
    events_end=find(diff_events_logical(:,n)==-1);
    if length(events_beg)>length(events_end)
        events_beg(end)=[];
    elseif length(events_beg)<length(events_end)
        events_end(1)=[];
    end
    num_events_cell=length(events_beg);
    for k=1:num_events_cell
        num_simple_events_temp=0;
        is_max_amp_true_event=0;
        multiEventsMatrix(events_beg(k)+1,n)=1;
        eventlengthMatrix(events_beg(k)+1,n)=events_end(k)-events_beg(k);
        last_detected_amp=0;
        temp_event=events(events_beg(k)+1:events_end(k),n);
        temp_trace=traces(events_beg(k)+1:events_end(k),n);
        areaUnderCurve(events_beg(k)+1,n)=sum(temp_trace)/fs;
        if events_beg(k)-round(fs)>0
            expanded_event=smooth_normalized_traces(events_beg(k)-round(fs):events_end(k),n);
            expanded_trace=traces(events_beg(k)-round(fs):events_end(k),n);
        else
            expanded_event=smooth_normalized_traces(1:events_end(k),n);
            expanded_trace=traces(1:events_end(k),n);
        end
        pos=sum(diff(temp_event)>0);
        neg=sum(diff(temp_event)<0);
        % Specific change for 2-photon acquired at less than 3Hz
        if neg==0
            if fs<3
                neg=1;
            end
        end
        if (neg/(fs*tau)>1 & length(temp_event)>2 & fs>3) || (neg/(fs*tau)>1 & length(temp_event)>0 & fs<3)
            [amp_max, ind_max]=max(temp_event);
            if events_beg(k)+ind_max-dt_samp_min>0
                if sum(eventsMatrix(events_beg(k)+ind_max-dt_samp_min:events_beg(k)+ind_max,n))==0
                    eventsMatrix(events_beg(k)+ind_max,n)=traces(events_beg(k)+ind_max,n);
                    num_simple_events_temp=num_simple_events_temp+1;
                    is_max_amp_true_event=1;
                    if events_beg(k)-round(fs)>0
                        expanded_peak=flip(expanded_event(1:ind_max+round(fs)));
                        expanded_peak_trace=flip(expanded_trace(1:ind_max+round(fs)));
                    else
                        expanded_peak=flip(expanded_event(1:ind_max+events_beg(k)));
                        expanded_peak_trace=flip(expanded_trace(1:ind_max+events_beg(k)));
                    end
                    last_min=find(diff(expanded_peak)>0,1,'first');
                    if isempty(last_min)
                        last_min=length(expanded_peak);
                        last_min_amp=expanded_peak_trace(end);
                    else
                        last_min_amp=expanded_peak_trace(last_min);
                    end
                    relativeAmpMatrix(events_beg(k)+ind_max,n)=traces(events_beg(k)+ind_max,n)-last_min_amp;
                    onsetWidthMatrix(events_beg(k)+ind_max,n)=last_min;
                end
            end
            if length(temp_event)>2
                [~, ind]=findpeaks(temp_event);
                if length(ind)>1
                    for m=1:length(ind)
                        if ind(m)>ind_max
                            last_detected_amp=amp_max;
                        end
                        if m==1 % the first peak:
                            this_peak=temp_event(ind(m):end);
                            if events_beg(k)-round(fs)>0
                                expanded_peak=flip(expanded_event(1:ind(m)+round(fs)));
                                expanded_peak_trace=flip(expanded_trace(1:ind(m)+round(fs)));
                            else
                                expanded_peak=flip(expanded_event(1:ind(m)+events_beg(k)));
                                expanded_peak_trace=flip(expanded_trace(1:ind(m)+events_beg(k)));
                            end
                            next_min=find(diff(this_peak)>0,1,'first');
                            next_min_ind=ind(m)+next_min-1;
                            last_min=find(diff(expanded_peak)>0,1,'first');
                            if isempty(last_min)
                                last_min=length(expanded_peak);
                                last_min_amp=expanded_peak_trace(end);
                            else
                                last_min_amp=expanded_peak_trace(last_min);
                            end
                            amp_change=temp_event(ind(m))-temp_event(next_min_ind);
                            if amp_change>thresh/2
                                if sum(eventsMatrix(events_beg(k)+ind(m)-dt_samp_min:events_beg(k)+ind(m),n))==0
                                    eventsMatrix(events_beg(k)+ind(m),n)=traces(events_beg(k)+ind(m),n);
                                    num_simple_events_temp=num_simple_events_temp+1;
                                    last_detected_amp=temp_event(ind(m));
                                    onsetWidthMatrix(events_beg(k)+ind(m),n)=last_min;
                                    relativeAmpMatrix(events_beg(k)+ind(m),n)=traces(events_beg(k)+ind(m),n)-last_min_amp;
                                end
                            end
                            
                        else % the rest of the peaks:
                            this_peak=flip(temp_event(1:ind(m))');
                            this_trace=flip(temp_trace(1:ind(m))');
                            last_min=find(diff(this_peak)>0,1,'first');
                            last_min_amp=this_trace(last_min);
                            last_min_ind=1+length(this_peak)-last_min;
                            amp_change=temp_event(ind(m))-temp_event(last_min_ind);
                            if amp_change>thresh | (temp_event(ind(m))-last_detected_amp)>thresh
                                if sum(eventsMatrix(events_beg(k)+ind(m)-dt_samp_min:events_beg(k)+ind(m),n))==0
                                    if m==length(ind) % the last peak
                                        eventsMatrix(events_beg(k)+ind(m),n)=traces(events_beg(k)+ind(m),n);
                                        num_simple_events_temp=num_simple_events_temp+1;
                                        last_detected_amp=temp_event(ind(m));
                                        onsetWidthMatrix(events_beg(k)+ind(m),n)=(ind(m)-last_min_ind);
                                        relativeAmpMatrix(events_beg(k)+ind(m),n)=traces(events_beg(k)+ind(m),n)-last_min_amp;
                                    else
                                        this_peak=temp_event(ind(m):end);
                                        next_min=find(diff(this_peak)>0,1,'first');
                                        next_min_ind=ind(m)+next_min;
                                        amp_change=temp_event(ind(m))-temp_event(next_min_ind);
                                        if amp_change>thresh/2
                                            eventsMatrix(events_beg(k)+ind(m),n)=traces(events_beg(k)+ind(m),n);
                                            num_simple_events_temp=num_simple_events_temp+1;
                                            last_detected_amp=temp_event(ind(m));
                                            onsetWidthMatrix(events_beg(k)+ind(m),n)=(ind(m)-last_min_ind);
                                            relativeAmpMatrix(events_beg(k)+ind(m),n)=1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        numEventsMatrix(events_beg(k)+1,n)=num_simple_events_temp;
        if num_simple_events_temp==1 & is_max_amp_true_event==1
            relativeAmpMatrix(events_beg(k)+ind_max,n)=eventsMatrix(events_beg(k)+ind_max,n);
        end
    end
end
close(h);

temp_binary_activity=zeros(size(eventsMatrix));
temp_binary_activity(eventsMatrix>0)=1;

end



