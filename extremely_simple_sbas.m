% simple_sbas_gmtsar_geocoded.m -- script to read in GMTSAR interferograms,
% and then do SBAS on 'em
%
% based on my lics geotiff code simple_sbas_geotiff.m
%
% i get by with a little help from the internets (my friends?)
%
% gjf, 18-apr-2020
% geniusinaction.com
%
% modification history 
%
% 20-mar-2020  gjf  first stab at a aria reader
% 24-mar-2020  gjf  added interferogram info output
% 17-apr-2020  gjf  modified for radar geometry GMTSAR grd files
% 27-apr-2020  gjf  modified again for geocoded GMTSAR grd files
% 16-feb-2022  gjf  first attempt at a 'unified' code

%% Step 1: define some key parameters

% start from fresh
clear;

% define your interferogram list
ifglist='ifgdates_doy.txt';  % a text list of the interferograms you want
                             % to include in your analysis, two columns for
                             % date1 date2, in the date format used for 
                             % naming 

% some data parameters
trackname=166;
trackgeom='A';
wavl=55.5;   % use the units you want output in
signconv=-1; % 1 for range change, -1 for ground displacement

% coherence threshold (values for the GMTSAR variant are lower in general)
cohthresh=0.05;

% and reference pixel
reflon=-114.96513;
reflat=34.84444;
searchpx=3;   % number of pixels to average over, centered on ref pixel

% what is your data source?
dsource='gmtsar';    % 'gmtsar', 'lics' and 'aria' are all valid choices


%% Step 2: load in the data

% which data loader do we need to run here?

switch dsource
    case 'gmtsar'
        disp('GMTSAR-format interferograms selected')
        read_data_gmtsar
    case 'aria'
        error('ARIA-format interferograms are not implemented yet. Soz.')
    case 'lics'
        error('LiCS-format interferograms are not implemented yet. My bad.')
    otherwise
        error('naah mate, you''ve had a mare')
end


%% Step 3: some optional plotting

% coloraxis
cmin=-15;
cmax=15;

% plot all the interferograms! 

% nplots=ceil(nifgs/100);
% startlayer=0;
% 
% for j=1:nplots
% 
%    plotrows=ceil(sqrt(nifgs/nplots));
%  
%    endlayer=startlayer+plotrows^2;
%    
%    if (endlayer>nifgs)
%        endlayer=nifgs;
%    end
%    
%    figure;
%    for i=1:(endlayer-startlayer)
%        % now with added NaN mask
%        nanmask=ones(size(ifgdata(:,:,startlayer+i)));
%        nanmask(isnan(ifgdata(:,:,startlayer+i)))=0;
%     
%        subplot(plotrows,plotrows,i);
%        imagesc(ifgdata(:,:,startlayer+i)-nanmean(nanmean...
%            (ifgdata(:,:,startlayer+i))),'AlphaData',nanmask,...
%            [cmin cmax]); axis image
%        colormap jet
%        title(sprintf('%d-%d',ifgdates(i,1),ifgdates(i,2)))
%    end
%    colorbar('location','Manual','position',[0.93 0.1 0.02 0.81])
% 
%    startlayer=startlayer+plotrows^2;   
%    
% end

% some diagnostic stuffs
figure; spy(A)
figure; spy(N)
figure; imagesc(T); axis image; colorbar
figure; plot((t1t2(:,2)+t1t2(:,1))/2,n,'o'); datetick('x','yyyy-mm-dd')

% view dates and coherent fractions
%[datestr(t1t2(:,1),'yyyy-mm-dd') repmat(' ',[nifgs,1]) datestr(t1t2(:,2),'yyyy-mm-dd') repmat(' ',[nifgs,1]) num2str(n)]

%% Step 4: look at phase closures

% identify closed loops in the set of interferograms and calculate the
% closure phase... hopefully it will be zero(...?)

fprintf(1,'\nChecking for phase closures...\n')

cnt_loops=0;                  % closed loop counter
ifg_loops=zeros(nepochs,3);   % indices of closed loop interferograms
%ifgdata_masked=ifgdata;       % somewhere to put the output of this

for i=1:nepochs-1    
   testA=zeros(1,nepochs);  % make a row vector to match
   testA(i)=1;
   i_ifg1=find(ismember(A,testA,'rows')); % find row of A that matches it
   testA=zeros(1,nepochs);  % and again with a different ifg
   testA(i+1)=1;
   i_ifg2=find(ismember(A,testA,'rows')); % find row of A that matches that
   testA(i)=1;
   i_ifg3=find(ismember(A,testA,'rows')); % and once more
   
   if not(isempty(i_ifg1) || isempty(i_ifg2) || isempty(i_ifg3))
      fprintf(1,'closed loop between %d, %d and %d\n',i_ifg1,i_ifg2,i_ifg3)
      cnt_loops=cnt_loops+1;
      ifg_loops(cnt_loops,:)=[i_ifg1 i_ifg2 i_ifg3];
   end
end

ifg_loops(cnt_loops+1:end,:)=[];

% closed loop should be column1 + column2 - column3

for i=1:cnt_loops
   i_ifg1=ifg_loops(i,1); 
   i_ifg2=ifg_loops(i,2); 
   i_ifg3=ifg_loops(i,3); 

   % make the closed loop
   pha_closure=ifgdata(:,:,i_ifg1)+ifgdata(:,:,i_ifg2)-ifgdata(:,:,i_ifg3);
   % convert to number of cycles
   pha_closure_cycles=round(pha_closure/2/pi);
   
   
   % masks for plotting 
   nanmask1=ones(size(pha_closure));
   nanmask2=ones(size(pha_closure));
   nanmask3=ones(size(pha_closure));
   nanmaskloop=ones(size(pha_closure));
   nanmask1(isnan(ifgdata(:,:,i_ifg1)))=0;
   nanmask2(isnan(ifgdata(:,:,i_ifg2)))=0;
   nanmask3(isnan(ifgdata(:,:,i_ifg3)))=0;
   nanmaskloop(isnan(pha_closure))=0;
   
   % plot them if you want!
%    figure;
%    subplot(2,2,1); imagesc(ifgdata(:,:,i_ifg1),'AlphaData',nanmask1);...
%        axis image; colorbar; title(sprintf('ifg %d',i_ifg1));
%    subplot(2,2,2); imagesc(ifgdata(:,:,i_ifg2),'AlphaData',nanmask2);... 
%        axis image; colorbar; title(sprintf('ifg %d',i_ifg2));
%    subplot(2,2,3); imagesc(ifgdata(:,:,i_ifg3),'AlphaData',nanmask3);... 
%        axis image; colorbar; title(sprintf('ifg %d',i_ifg3));
%    subplot(2,2,4); imagesc(pha_closure_cycles,'AlphaData',nanmaskloop);... 
%        axis image; colorbar; title('closure phase');

   % make an NaN error mask for these interferograms...
   errormask=ones(size(pha_closure));
%   errormask(isnan(pha_closure))=nan;
%   errormask(pha_closure_cycles~=0)=nan;   % crudest option
   errormask(pha_closure>1.5)=nan;          % more conservative option
   errormask(pha_closure<-2)=nan;
   
   % ...and apply it
   ifgdata(:,:,i_ifg1)=ifgdata(:,:,i_ifg1).*errormask;
   ifgdata(:,:,i_ifg2)=ifgdata(:,:,i_ifg2).*errormask;   
   ifgdata(:,:,i_ifg3)=ifgdata(:,:,i_ifg3).*errormask;
   
end

%% Step 4a: some more optional plotting

% let's recalculate N

% N = ones(size(ifgdata(:,:,1)));
N=~isnan(sum(ifgdata,3));

% coloraxis
cmin=-15;
cmax=15;

% plot all the interferograms! 

% nplots=ceil(nifgs/100);
% startlayer=0;
% 
% for j=1:nplots
% 
%    plotrows=ceil(sqrt(nifgs/nplots));
%  
%    endlayer=startlayer+plotrows^2;
%    
%    if (endlayer>nifgs)
%        endlayer=nifgs;
%    end
%    
%    figure;
%    for i=1:(endlayer-startlayer)
%        % now with added NaN mask
%        nanmask=ones(size(ifgdata(:,:,startlayer+i)));
%        nanmask(isnan(ifgdata(:,:,startlayer+i)))=0;
%     
%        subplot(plotrows,plotrows,i);
%        imagesc(ifgdata(:,:,startlayer+i)-nanmean(nanmean...
%            (ifgdata(:,:,startlayer+i))),'AlphaData',nanmask,...
%            [cmin cmax]); axis image
%        colormap jet
%        title(sprintf('%d-%d',ifgdates(i,1),ifgdates(i,2)))
%    end
%    colorbar('location','Manual','position',[0.93 0.1 0.02 0.81])
% 
%    startlayer=startlayer+plotrows^2;   
%    
% end


%% Step 5: try some SBAS?

fprintf(1,'\nLet us try some SBAS...\n')

% use a truncated svd to invert A
% you need to test with different numbers of singular values to choose an 
% appropriate number for nsv
[U,S,V]=svd(A,0);
nsv=85;
struncinv=1./diag(S);
struncinv(nsv+1:end)=zeros(size(struncinv(nsv+1:end)));
Struncinv=diag(struncinv);
H=V*Struncinv*U';               % the generalized inverse!



% find all the complete time series
[r,c]=find(N==1);               % row and column indices of coherent pixels
ptidx=find(N==1);               % linear indices of the same
nts=length(r);

% find the averaged reference point time series

refts=zeros(nepochs,searchpx^2);
cnt=0;

%loop through pixels around the reference point
for i=refrow-refshft:refrow+refshft
    for j=refcol-refshft:refcol+refshft
      
       cnt=cnt+1;
       dref=squeeze(ifgdata(i,j,:));
       %mref=H*dref;   % if using svd
       mref=A\dref;    % if using backslash
       refts(:,cnt)=mref;
       
    end
end

mref=nanmean(refts')';              % mean reference time series

tic
% and loop over them
for i=1:nts

    d=squeeze(ifgdata(r(i),c(i),:));
    m=H*d;

    sbasinc(r(i),c(i),:)=m-mref;
    sbascum(r(i),c(i),:)=cumsum(m);
    
end
t=toc;
fprintf(1,'SBAS complete in %f s\n',t)

%figure; plot(tn,cumsum(m)); datetick('x','yyyy-mm-dd')

%cmin=min(min(min(sbascum)));
%cmax=max(max(max(sbascum)));
cmin=-15;
cmax=15;

% plot all the epochs! 
% plotrows=ceil(sqrt(nepochs));
% figure;
% for i=1:nepochs
%    subplot(plotrows,plotrows,i);
%    imagesc(sbascum(:,:,i),'AlphaData',N,[cmin cmax]); axis image; colormap jet
%    title(datestr(tn(i),'yyyy-mm-dd'))
% end
% colorbar('location','Manual','position',[0.93 0.1 0.02 0.81])
% figure;
% for i=1:nepochs
%    subplot(plotrows,plotrows,i);
%    imagesc(sbasinc(:,:,i),'AlphaData',N,[cmin cmax]); axis image; colormap jet
%    title(datestr(tn(i),'yyyy-mm-dd'))
% end
% colorbar('location','Manual','position',[0.93 0.1 0.02 0.81])


%% Step 6: export the time series

% coordinate sanity
if(min(lon)>180)
   lon=lon-360;
end

fprintf(1,'\nnow extracting time series for points of interest (pois)\n');

% plot individual time series for pixels of interest (pois)

poi=[-116.712089 33.733331; -116.209142 34.699427; -114.004952 32.722545;...
    -116.169785 32.616524; -116.569323 33.059163; -115.765846 33.516411;...
    -116.710253 34.192457; -116.67145 34.93683; -116.879623 34.462064;...
    -116.571644 35.320643; -114.599404 34.188934; -114.759349 32.497872;...
    -115.88373 32.417007; -114.453859155 35.183120004;...
    -115.307320823 32.4169609916; -115.484626567 33.1661580465;
    -114.347563602 34.4715911005; -116.54323095 33.8235198829;
    -115.366441305 34.2650047143; -115.564796754 34.8716557886;
    -115.860233429 33.4602590452; -114.630723358 33.5699305303;
    -115.878094547 32.5855684071; -114.468097034 32.9750945841];

poinames=["P472"; "P608"; "BMHL"; "P604"; "P598"; "P483"; "P504";... 
    "P623"; "P503"; "PJZX"; "P473"; "P502"; "RD10"; "RD11"; "RD12";...
    "RD13"; "RD14"; "RD15"; "RD16"; "RD17"; "RD18"; "RD19"; "P625";];

poinames=["DSSC"; "LDSW"; "P003"; "P066"; "P483"; "P504"; "P598";... 
    "P604"; "P606"; "P617"; "P623"; "P796"; "PJZX"; "THUM"; "RD10";...
    "RD11"; "RD14"; "RD15"; "RD16"; "RD17"; "RD18"; "RD19";...
    "RD22"; "RD24"];

% get the motion of the reference point (should be zero)
tsrefpoits=zeros(size(mref));
npoi=length(poi);

% vertices of the area of interest
xv=[min(lon); min(lon); max(lon); max(lon)];
yv=[min(lat); max(lat); max(lat); min(lat)];

for i=1:npoi

    % inpolygon test for the poi
    in=inpolygon(poi(i,1),poi(i,2),xv,yv);
    
    if(in==1)
    
        % find row and column for each poi (in data and metadata)
        [poirow,poicol]=latlon2pix(R,poi(i,2),poi(i,1));
        %[poimetarow,poimetacol]=latlon2pix(Rmeta,poi(i,2),poi(i,1));
        poirow=round(poirow);
        poicol=round(poicol);
        %poimetarow=round(poimetarow);
        %poimetacol=round(poimetacol);

        % extract the information for those locations
        
        %poiazi=azimuth(poimetarow,poimetacol);
        %poiinc=incidence(poimetarow,poimetacol);
        
        %poienu=squeeze(pointing_vector_components(poirow,poicol,:))*-1;  % LICS pointing vector is from ground to sat
        
        % find the mean poi time series (including neighbours)
        poinbts=zeros(nepochs,searchpx^2);  % poi and neighbours time series
        cnt=0;

        %loop through pixels around the reference point
        for j=poirow-refshft:poirow+refshft
            for k=poicol-refshft:poicol+refshft
                if(sbascum(j,k,end)~=0)
                    cnt=cnt+1;
                    poinbts(:,cnt)=squeeze(sbascum(j,k,:))*wavl/4/pi;
                end
            end
        end

        poits=nanmean(poinbts(:,1:cnt)')';              % mean poi time series
        
        % sanity check (for output data)
        if(nansum(poits)==0)
            fprintf(1,' poi %s has no nonzero time series\n',poinames{i}); 
        elseif(length(poits)~=nepochs)
            fprintf(1,' poi %s has incomplete time series\n',poinames{i}); 
        else
            % plot it
            figure; plot(tn,poits-tsrefpoits,'o'); datetick('x','yyyy-mm-dd')
                    title(sprintf('T%03d %s: %fE, %fN',trackname,poinames(i),...
                poi(i,1),poi(i,2)));
            % open an output file
            poioutfile=sprintf('T%03d_%s_%fE_%fN.csv',trackname,poinames(i),...
                poi(i,1),poi(i,2));
            poioutid=fopen(poioutfile,'wt');
            % write headers
            fprintf(poioutid,'%% Track %s\n',trackname);
            fprintf(poioutid,'%% Site %s\n',poinames(i));
            fprintf(poioutid,'%% Longitude %f\n',poi(i,1));
            fprintf(poioutid,'%% Latitude %f\n',poi(i,2));
            %fprintf(poioutid,'%% LOS_E %f\n',poienu(1));
            %fprintf(poioutid,'%% LOS_N %f\n',poienu(2));
            %fprintf(poioutid,'%% LOS_U %f\n',poienu(3));
            fprintf(poioutid,'%% date, datenum, displacement_m \n');
            fprintf(poioutid,'%s, %d, %f\n',datestr(t0n,'yyyy-mm-dd'),t0n,0); % t=0

            % and write the displacement time series
            for j=1:nepochs
                fprintf(poioutid,'%s, %d, %f\n',datestr(tn(j),'yyyy-mm-dd'),tn(j),poits(j)-tsrefpoits(j));   
            end

            fclose(poioutid);
        end
    else
        % if we're outside the region
        fprintf(1,' poi %s is not in the area of coverage\n',poinames{i});
    end
end


%% Step 7: export some displacements

for i=1:nepochs
    % extract date info for yer epoch
    epochvec=datevec(epochs(i));
    % make a file name
    outgrdfile=sprintf('%04d%02d%02d.grd',epochvec(1),epochvec(2),...
        epochvec(3));    
    % and write it out - couldn't do it without this function!
    % thanks to kelsey jordahl and the mathworks file exchange!
    % https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2
    grdwrite2(lon,lat,sbascum(:,:,i)*wavl/4/pi,outgrdfile);
    
end

%% Step 8: export velocities

tny=(tn-tn(1))/365.25;

velmap=zeros(nrows,ncols);

for i=1:nrows
    
    for j=1:ncols
        
        ts=squeeze(sbascum(i,j,:));
        
        if(mean(ts)~=0)
            
            G=[tny sin(2*pi*tny) cos(2*pi*tny) ones(size(ts))];
%           A=[tny sin(2*pi*tny) cos(2*pi*tny) sin(pi*tny) cos(pi*tny) ones(size(ts))];
%           A=[tny sin(2*pi*tny) cos(2*pi*tny) sin(pi*tny) cos(pi*tny) sin(2*pi*tny/3) cos(2*pi*tny/3) ones(size(ts))];

            m=G\ts;
            tsprime=G*m;
            
            velmap(i,j)=m(1);
   
        end       
    end
end

grdwrite2(lon,lat,velmap*wavl/4/pi,'velocities.grd');

%% Step 9: make a nice plot?

%figure; plot(tn,squeeze(sbascum(320,320,:)),'o'); datetick('x','yyyy-mm-dd')



% output cumulative displacements
%totdisp=sbascum(:,:,end);
%outdata=[LON(ptidx) LAT(ptidx) totdisp(ptidx)*wavl/4/pi];
%save total_displacements.dat outdata -ascii

%dlmwrite('total_displacements.txt',outdata);

%save sbascum.mat sbascum -v7.3

% this plays nice with zeros when plotting geographically
figure; 
t=geoshow(velmap*wavl/4/pi,R,'DisplayType','surface'); 
alpha(t,double(velmap~=0))
colorbar

