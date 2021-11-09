% Firstly set up the parameters for the evolutionary matrix
% BLOUSUMno=30;
% ResultFile=['Distance' int2str(BLOUSUMno) '.mat'];
% 
% Prepare the main data file

BLOSUMno=30;
ResultFile = ['Distance_' int2str(BLOSUMno) '_2021.mat'];
[Headers, SeqV] = fastaread('all_viral_interacting_2021.fasta');
% just to check if there are duplicates 
% check = {};
% for a=1:1:length(Headers)
%     test = strsplit(Headers{a},"|");
%     entry = test{2};
%     check{a} = entry;
% end
% filter out the ones that do not appear in ProteinB_Unique
true_member={};
for a=1:1:length(Headers)
    check = strsplit(Headers{a},"|");
    entry = check{2};
    if ismember(entry,ProteinB_Unique)
        true_member{end+1} = Headers{a};
    end
end

% store the sequences for true_member
true_SeqV={};
for i=1:length(Headers)
    check = strsplit(Headers{i},"|");
    entry = check{2};
    if ismember(entry,ProteinB_Unique)
        true_SeqV{end+1} = SeqV{i};
    end
end
% global alignment 
VpCount = length(true_SeqV);
ScoringMatrix = ['BLOSUM' num2str(BLOSUMno)];
ScoreVps = zeros(VpCount,VpCount);
parfor a = 1:1:VpCount
    for b = 1:1:VpCount
        ScoreVps(a,b) = nwalign(true_SeqV{1,a},true_SeqV{1,b},'ScoringMatrix',ScoringMatrix);
    end
end

% examine the matrix 
imagesc(ScoreVps);
colorbar;

% having checked the plot, looks like 233 is an outlier 
% below remove the outlier 
% ScoreVps(:,233) = 0;
% ScoreVps(233,:) = 0;
% 
Weights=zeros(VpCount,VpCount);
for r=1:VpCount
    Row=ScoreVps(r,:);
    %Normalize on [0-1] scale, then reverse it (1-x)
    Weights(r,:)=1-((Row-min(Row))/(max(Row)-min(Row)));
end

save(ResultFile,'VpCount','BLOSUMno','Weights','true_member','true_SeqV')