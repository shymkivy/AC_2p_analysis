function [dfofGRADpos, dfofGRADposN] = calc_ddfof(data, halos)
% this calculates derivative of dfof of calcium traces, with Jordans method

dfof = 1;

% extract the number of cells and the number of frames
total_frames = size(data,2);
num_cells = size(data,1);


data_halosub=data-halos;
% get minimum value of each ROI
mins=min(data_halosub,[],2);

% adjust baseline to 1
data_halosub_norm = data_halosub - mins + 1;


% dfof calculation (Jordan's code)
dfof2=zeros(size(data_halosub_norm));


for c=1:num_cells
    % take the min and max within a 20 bin window
    if dfof == 1
        [minimg, maximg] = minmaxfilt(data_halosub_norm(c,:),10,'minmax','same');
        baso=minimg+((maximg-minimg)/10);
        baso=smooth(baso,4);
        for t=21:size(data,2)
            a=(baso(t-20:t-1,1));
            m=prctile(a,20);
            dfof2(c,t)=(data_halosub_norm(c,t)-m);
        end
        %disp(c);
    elseif dfof == 0
        dfof2=data;
    end
    dfofall=dfof2;
end


% smooth the dfof signal and take the first derivative
% dfofFOOPall=run_fast_oopsi((dfofall));
dfofFOOPall=dfofall-dfofall; 
for c=1:num_cells    % size(data,1)
    dfofFOOPall(c,1:end-1)=diff(smooth(dfofall(c,:),10,'lowess'));
    %dfofFOOPall(c,:)=diff(smooth(dfofall(c,:),4));
end



% make positive and normalize
dfoftmp=dfofFOOPall;
for c=1:num_cells %size(dfofbinary,1)
    mk=max(dfoftmp(c,:));
    for t=1:total_frames %size(dfofbinary,2)
        if or(dfoftmp(c,t)==0,dfoftmp(c,t)>0)
            dfofGRADpos(c,t)=dfoftmp(c,t);
            dfofGRADposN(c,t)=dfoftmp(c,t)./mk;
        end
    end
end

    
%     
%     % plot individual traces, if you want
%     for ii = 1:20
%         figure;
%         %plot(data_halosub_norm(1,:));
%         % hold on;
%         subplot(3,1,1);
%         plot(dfofall(ii,:))
%         title('dfof');
%         subplot(3,1,2);
%         plot(smoothed_dfofall(ii,:));
%         title('smoothed dfof');
%         subplot(3,1,3);
%         plot(dfofGRADpos(ii,:));
%         title(sprintf('derivative of dfof, cell %d', ii));
%     end

end