clear;

[junk1,chr,start,stop] = textread('human_probes.txt','%s %s %d %d');

uchr = unique(chr);


filenames = textread('cgmapfilesturk.txt','%s');

methmat = zeros(length(chr),length(filenames));
countsmat = zeros(length(chr),length(filenames));
maxcountsmat = zeros(length(chr),length(filenames));

for fidx=1:length(filenames)
    %fidx
    filenames{fidx}
    
    filepath = sprintf('%s',filenames{fidx});
    
    str_chr = 'chr';
    
    clear chrp;
    clear pos;
    clear countsC;
    clear countsTot;
    
    [chrp nuc pos CG CGP meth countsC countsTot] = textread(filepath,'%s %c %d %s %s %f %d %d');
    
    for i=1:length(uchr)
        idx{i} = find(strcmp(append(str_chr,uchr{i}),chrp));   % Find all methylation values for sites in all chromosomes. 
    end

    for i=1:length(chr)

        idx2 = strcmp(chr{i},uchr);   % Filter out which chromosome I want.
        idxp = idx{idx2};  % Select only the sites within the chromosome I want.

        distbp = pos(idxp) - start(i) + 60;

        idxprobe = find(abs(distbp) < 100);

        sumC = sum(countsC(idxp(idxprobe)));
        sumTot = sum(countsTot(idxp(idxprobe)));
        maxTot = max(countsTot(idxp(idxprobe)));
        
        if sumTot > 0
            methmat(i,fidx) = sumC/sumTot;
            countsmat(i,fidx) = sumTot;
            maxcountsmat(i,fidx) = maxTot;
        else
            methmat(i,fidx) = NaN;
            countsmat(i,fidx) = NaN;
            maxcountsmat(i,fidx) = NaN;
        end
        
        
        
    end  
end
    
kitid = start;

writematrix([kitid methmat],'turk_meth_143.csv');
writematrix([kitid countsmat],'turk_counts_143.csv');
writematrix([kitid maxcountsmat],'turk_counts_max_143.csv');





