% Firstly set up the parameters for the evolutionary matrix
% BLOUSUMno=30;
% ResultFile=['Distance' int2str(BLOUSUMno) '.mat'];
% 
% Prepare the main data file
opts = detectImportOptions('9606_2021_cleaned');
opts.Delimiter=",";
data = readtable('9606_2021_cleaned',opts);
ProteinA_Unique = unique(data.ProteinA);
ProteinB_Unique = unique(data.ProteinB);
ProteinA = data.ProteinA;
ProteinB = data.ProteinB;
GeneA = data.GeneA;
GeneB = data.GeneB;
SubfamilyA = data.FamilyA;
SubfamilyB = data.FamilyB;
Score = data.Score;
PMID = data.PMID;
Sub_Unique = unique(SubfamilyB);
Tax_Unique = unique(data.TaxonB);
TaxonA = unique(data.TaxonA);
TaxonB = unique(data.TaxonB);
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

% global alignment 
VpCount = length(SeqV);
ScoringMatrix = ['BLOSUM' num2str(BLOSUMno)];
ScoreVps = zeros(VpCount,VpCount);
parfor a = 1:1:VpCount
    for b = 1:1:VpCount
        ScoreVps(a,b) = nwalign(SeqV{1,a},SeqV{1,b},'ScoringMatrix',ScoringMatrix);
    end
end

% examine the matrix 
imagesc(ScoreVps);
colorbar;

% having checked the plot, looks like 233 is an outlier 
% below remove the outlier 
ScoreVps(:,233) = 0;
ScoreVps(233,:) = 0;

Weights=zeros(VpCount,VpCount);
for r=1:VpCount
    Row=ScoreVps(r,:);
    %Normalize on [0-1] scale, then reverse it (1-x)
    Weights(r,:)=1-((Row-min(Row))/(max(Row)-min(Row)));
end

save(ResultFile,'VpCount','BLOUSUMno','Weights')