function f_jordan_ensembleanalysis_2018(firing_rates)

peakframes=.99;
shuffsim = 1;
corrHIGHstate = 1; %the cell-cell correlation analysis can be done on either the continuous data (all timepoints; corrHIGHstate=0) or on the ensemble activations (corrHIGHstate=1).
corrsimil = 1;
shuffcor = 1;

CC = 1;
fil = 1;

% normalize
actpre = firing_rates;
actpre=actpre./max(actpre,[],2);

for framesimilarity=1
    disp('identifying ensemble activations');
    % find local peaks
    [zz, zzz] = if_find_peaks2(sum(actpre,1), 50);

    num_shuff_reps = 100;
    for shufflingactivity_definitions=1
        shufaz=zeros(num_frames, num_shuff_reps);
        %shufaz2 = zeros(num_frames, num_shuff_reps);
        for xx=1:num_shuff_reps
            tmp2=actpre-actpre;
            for n_cell=1:num_cells
                tmp2(n_cell,:) = circshift(actpre(n_cell,:),randsample(num_frames,1));
            end
            [~, shufaz(:,xx)] = if_find_peaks2(sum(tmp2,1), 50);
        end
    end
    cuto1=prctile(shufaz(:),peakframes*100);
    cuto2=101010;       
    acttmp=actpre(:,and(zz==1,and(zzz>cuto1,zzz<cuto2)));
    %simmat=zeros(size(acttmp,2),size(acttmp,2)); 
    disp(' complete  1/5');


    disp('comparing ensemble activations');
    simmat = similarity_index(acttmp',acttmp');
    simmat = simmat - diag(ones(size(simmat,1),1)*2);
    simmats{CC,fil}=simmat;%(sz:end,1:sz);
    disp('complete  2/5');

    num_upst = size(acttmp,2);
    for shuffsim=shuffsim
        disp('shuffling ensemble activations');
        if shuffsim>0
            mn=0;
            screeshuffs=[];
            for z=1:100
                tmp2=acttmp;
                % shuffle between cells
                if shuffsim==1
                    for t=1:size(tmp2,2)
                        a=randsample(num_cells, num_cells);
                        tmp2(:,t)=tmp2(a,t);
                    end
                end

                % shuffle between times
                if shuffsim==2
                    for n_cell=1:num_cells
                        tmp2(n_cell,:) = circshift(tmp2(n_cell,:),randsample(num_upst,1));
                    end
                end

                % shuffle cells and times
                if shuffsim==3             
                    for n_cell=1:num_cells
                        tmp2(n_cell,:) = circshift(tmp2(n_cell,:),randsample(num_upst,1));
                    end

                    for t=1:size(tmp2,2)
                        a=randsample(num_cells, num_cells);
                        tmp2(:,t)=tmp2(a,t);
                    end
                end

                tmpsimmat=zeros(floor(num_upst/10),floor(num_upst/10))-1;
                ab=randsample(5,1);
                gb=randsample(5,1)+5;
                ak1=1;
                for t1=ab:10:size(tmp2,2)
                    ak2=1;
                    for t2=gb:10:size(tmp2,2)
                        if t1~=t2
                            tmpsimmat(ak1,ak2)=dot(tmp2(:,t1),tmp2(:,t2))/(((norm(tmp2(:,t1)).^2)+(norm(tmp2(:,t2)).^2))/2);
                            %simmat(t1,t2)=min(min(corrcoef(acttmp(:,t1),acttmp(:,t2))));
                        end
                        ak2=ak2+1;
                    end
                    ak1=ak1+1;
                end
                mn=vertcat(mn,tmpsimmat(:));

            end
            tsmat=simmat;
            simmatZ=tsmat-tsmat;
            for st1=1:num_upst
                for st2=1:num_upst
                    % returns z score
                    simmatZ(st1,st2)=norminv(mean(tsmat(st1,st2)>mn));
                    if simmatZ(st1,st2)==-Inf
                        simmatZ(st1,st2)=norminv(1/size(mn,1));
                    elseif simmatZ(st1,st2)==Inf
                        simmatZ(st1,st2)=norminv((size(mn,1)-1)/size(mn,1));
                    end
                end
            end
            simmatsZ{CC,fil}=simmatZ;
            disp('complete  3/5');
        end
    end

    for custering=1
        disp('clustering')
        for doPREpca=1 %simplify/orthogonalize the data before clustering
            [coeffa,scorea,~,~,scree,mu] = pca(acttmp'); %spatial or "cell number reduction" PCA
            for t=1:size(scree,1) %only include components with greater than expected amount of explained varience for a single cell
                if scree(t)>((1/num_cells)*100)
                    pcs=t;
                end
            end
            clear scoreaa;
            % this gets scores without subtracting the mean mu
            scoreaa = coeffa(:,1:pcs)'*acttmp; 
            numcomps(fil,CC)=pcs;
        end %determines number of PCs and then does PCA

        for dokmeansshuffle=1 %determines number of clusters
            clear scree screeshuff
            for z=1:100 %make shufflescree
                for c=1:10
                    [~,~,sumd,~] = kmeans(scoreaa,c,'Distance','cosine');
                    scree(z,c)=sum(sumd(:));
                end
            end

            for z=1:100
                clear ac
                for t=1:pcs
                    ac(:,t)=scoreaa(randsample(size(scoreaa,1),size(scoreaa,1)),t);
                end
                for c=1:10
                    [~,~,sumd,~] = kmeans(ac,c,'Distance','cosine');
                    screeshuff(z,c)=sum(sumd(:));
                end
            end
            figure;
            ss=mean(screeshuff,3);
            sss=ss(:,:)-ss(:,:);
            sss(:,2:10)=(ss(:,2:10)-ss(:,1:9))./ss(:,1:9);
            mshuf=mean(sss,1);
            sshuff=std(sss,0,1);
            errorbar(1:10,mean(sss,1)',(std(sss,0,1))'); hold on
            ss2=mean(scree,1);
            clusterscree(1)=0;
            clusterscree(2:10)=(ss2(2:10)-ss2(1:9))./ss2(1:9); %actual "cluster
            plot(clusterscree,'k');
            clusses=1;
            for t=2:10
                if clusterscree(t)<(mshuf(t)-sshuff(2)*1.67)
                    clusses=t;
                end
            end
        end

        numensembles(fil,CC)=clusses;
        STATEidentifiers=kmeans(scoreaa(:,1:pcs),clusses);
        disp('complete 4/5');
    end
end


for cellcorrels=1

    disp('cell cell correlations');
    for corrHIGHstate=corrHIGHstate
        if corrHIGHstate==0
            acttmp=actpre;
        end
        docoactivation=0; %ignore this
        clear conmat clustmp clustmp2 corrs
        corrs=zeros(size(acttmp,1))-1;
        for c1=1:size(acttmp,1)
            for c2=1:size(acttmp,1)
                if c1~=c2
                    if corrsimil==0
                        corrs(c1,c2)=min(min(corrcoef(acttmp(c1,:),acttmp(c2,:))));
                    else
                        corrs(c1,c2)=dot(acttmp(c1,:),acttmp(c2,:))/(((norm(acttmp(c1,:)).^2)+(norm(acttmp(c2,:)).^2))/2);
                    end
                end
            end
        end
        cormats{CC,fil}=corrs;
        for shuffcor=shuffcor

            mn=0;
            screeshufs=[];
            for z=1:1000
                tmp2=acttmp;
                tmp1=horzcat(tmp2,tmp2,tmp2);
                for c=1:size(tmp2,1)
                    llk=randsample(size(tmp2,2),1);
                    tmp2(c,:)=tmp1(c,llk:llk+size(tmp2,2)-1);
                end
                if shuffcor==1
                    tmp2=acttmp;
                    for t=1:size(tmp2,2)
                        a=randsample(size(tmp2,1),size(tmp2,1));
                        tmp2(:,t)=tmp2(a,t);
                    end
                end
                if shuffcor==3
                    tmp3=tmp2;
                    for t=1:size(tmp2,2)
                        a=randsample(size(tmp2,1),size(tmp2,1));
                        tmp3(:,t)=tmp3(a,t);
                    end
                    tmp2=tmp3;
                end
                tmpcormat=zeros(size(tmp2,1))-1;
                for c1=1:size(acttmp,1)
                    for c2=1:size(acttmp,1)
                        if c1~=c2
                            if corrHIGHstate==1
                                tmpcormat(c1,c2)=min(min(corrcoef(tmp2(c1,:),tmp2(c2,:))));
                            elseif and(docoactivation==0,corrsimil==0)
                                %tmpcormat(c1,c2)=mean(xcorr(tmp2(c1,:),tmp2(c2,:),0,'coeff'));
                                tmpcormat(c1,c2)=min(min(corrcoef(tmp2(c1,:),tmp2(c2,:))));
                            end
                            if docoactivation==1
                                tmpcormat(c1,c2)=mean(tmp2(c1,:).*tmp2(c2,:))./sqrt(mean(tmp2(c1,:)).*mean(tmp2(c2,:)));
                            elseif corrsimil==1
                                tmpcormat(c1,c2)=dot(tmp2(c1,:),tmp2(c2,:))/(((norm(tmp2(c1,:)).^2)+(norm(tmp2(c2,:)).^2))/2);
                            end
                        end
                    end
                end
                mn=vertcat(mn,tmpcormat(tmpcormat>-1));
            end
            tsmat=cormats{CC,fil}; cormatZ=tsmat-tsmat;
            for st1=1:size(tsmat,1)
                for st2=1:size(tsmat,1)
                    cormatZ(st1,st2)=norminv(mean(tsmat(st1,st2)>mn));
                    if cormatZ(st1,st2)==-Inf
                        cormatZ(st1,st2)=norminv(1/size(mn,1));
                    elseif cormatZ(st1,st2)==Inf
                        cormatZ(st1,st2)=norminv((size(mn,1)-1)/size(mn,1));
                    end
                end
            end
            cormatsZ{CC,fil}=cormatZ;
        end
    end
    disp('complete 5/5');
end

figure; hold on;
plot(sum(actpre,1))
plot(ones(num_frames,1)*cuto1)

figure; imagesc(simmat)
figure; imagesc(simmatZ)


end

function [peaks_trace, peaks_vals] = if_find_peaks1(data, percentile)

num_frames = numel(data);
peaks_trace=data-data; %upstate idetnify
peaks_vals=data-data; %upstate values
local_median=prctile(data,percentile);
% find local peaks
for t=15:num_frames-15
    if and(data(1,t)>local_median,data(1,t)==max(data(1,t-1:t+1)))
        peaks_trace(1,t)=1;
        peaks_vals(1,t)=data(1,t);
    end
end

end

function [peaks_trace, peaks_vals] = if_find_peaks2(data, percentile)

num_frames = numel(data);
local_median=prctile(data,percentile);
peaks_trace = logical(((data - circshift(data,1))>0).*((data - circshift(data,-1))>0).*(data>local_median));
peaks_vals = zeros(1,num_frames);
peaks_vals(peaks_trace) = data(peaks_trace);

end