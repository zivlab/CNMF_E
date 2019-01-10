function [idxOfGoodICs] = trace_sorting(traces,thresh,tau,fs)
% This function gives a grade for each trace and eliminates those with
% lower grades

M=size(traces,1);
num_cells=size(traces,2);
dt=1/fs;
t=0:dt:(M-1)*dt;

fpass=0.5/tau;
[N,Hd]=LPF(fs,fpass);
smooth_traces=filter(Hd,traces);
smooth_normalized_traces=smooth_traces./repmat(median(abs(smooth_traces)),M,1);

events=zeros(size(traces));
events(smooth_normalized_traces>=thresh)=smooth_normalized_traces(smooth_normalized_traces>=thresh);

grade=zeros(1,num_cells);
for n=1:num_cells
    pos=sum(diff(events(:,n))>0); 
    neg=sum(diff(events(:,n))<0); 
    grade(n)=log10(pos)*neg/pos;
end

idxOfGoodICs=find(grade>thresh);

end

