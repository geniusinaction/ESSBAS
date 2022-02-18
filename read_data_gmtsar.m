% read_data_gmtsar.m 
%
% read in geocoded GMTSAR interferograms for SBAS analysis
%
% gjf, 16-feb-2022
% geniusinaction.com


% filenames in each directory
unwgrd='unwrap_ll_dsamp.grd';
corgrd='corr_ll_dsamp.grd';

% and let's load in some things: first the interferogram list

fid=fopen(ifglist);
inlist=textscan(fid,'%s %s');
fclose(fid);

% loop through and extract the dates

nifgs=length(inlist{1});
ifgdates=zeros(nifgs,2);
date1n=zeros(nifgs,1);
date2n=zeros(nifgs,1);

for i=1:nifgs

    C=textscan(inlist{1}{i},'%4d%3d');
    yr1=double(C{1});
    doy1=double(C{2});
    D=textscan(inlist{2}{i},'%4d%3d');
    yr2=double(D{1});
    doy2=double(D{2});
    
    date1n(i)=datenum(yr1,1,doy1)+1;  % add 1 because GMTSAR needs it
    date2n(i)=datenum(yr2,1,doy2)+1;
    
    
    ifgdates(i,:)=[yr1*1000+doy1 yr2*1000+doy2];
    
end

% what do we have here? 
fprintf(1,'\n%d interferograms in total\n\n',nifgs);

% load the first file for a look
date1=ifgdates(1,1);  
date2=ifgdates(1,2);

% construct the file names
unwfile=sprintf('%d_%d/%s',date1,date2,unwgrd);
corfile=sprintf('%d_%d/%s',date1,date2,corgrd);

% read in the files
fprintf(1,'opening %s\n\n',unwfile);

% unwphase=ncread(unwfile,'/z')'; % interferometric phase
% phsigcoh=ncread(corfile,'/z')'; % correlation (GMTSAR convention)
% 
% lon=ncread(unwfile,'/lon'); % range (x) pixels
% lat=ncread(unwfile,'/lat'); % azimuth (y) pixels

[lon,lat,unwphase]=grdread2(unwfile);
[lon,lat,phsigcoh]=grdread2(corfile);

% to prevent sign convention issues later on
if(min(lon)>180)
   lon=lon-360;
end
%[RNG,AZI]=meshgrid(rngx,aziy);
[nrows,ncols]=size(unwphase);

% mask the worst pixels

maskunwphase=unwphase.*(phsigcoh>cohthresh);
maskunwphase(maskunwphase==0)=NaN;

% geocoding information
R=georefpostings([min(lat) max(lat)],[min(lon) max(lon)],...
    size(maskunwphase),'ColumnsStartFrom','South');
[X,Y]=meshgrid(1:nrows,1:ncols);
[LAT,LON]=pix2latlon(R,X',Y');

% extract reference point phase (mean of a window centered on the point)
[refrow,refcol]=latlon2pix(R,reflat,reflon);
refrow=round(refrow);
refcol=round(refcol);
refshft=floor(searchpx/2);
refpha=nanmean(nanmean(maskunwphase(refrow-refshft:refrow+refshft,...
    refcol-refshft:refcol+refshft)));
fprintf(1,'  reference point phase: %f\n',refpha);

% make a 3d matrix to put all the data in
ifgdata=zeros(nrows,ncols,nifgs);

% and put this first one in there, subtracting reference phase, and applying
% sign convention
ifgdata(1:nrows,1:ncols,1)=(maskunwphase-refpha)*signconv;

% also, retain the reference phase
refphase=zeros(nifgs,1);
refphase(1)=refpha;

% do some things with dates
epochs=unique([date1n date2n]);                                        % get all the unique datenums
t0n=epochs(1);                                                         % datenum for 1st epoch
tn=zeros(length(epochs)-1,1);                                          % time for each epoch
for i=2:length(epochs)
    tn(i-1)=epochs(i);                                                 % the other datenums
end


nepochs=length(tn);                                                    % number of time series epochs
tspan=zeros(nifgs,1);                                                  % time span of each interferogram
tspan(1)=date2n(1)-date1n(1);
t1t2=zeros(nifgs,2);                                                   % datenums of interferogram dates
t1t2(1,:)=[date1n(1) date2n(1)];

% make a template for the design matrix
A=zeros(nifgs,nepochs);
startepoch=1;
endepoch=find(tn==date2n(1));
A(1,startepoch:endepoch)=ones(size(startepoch:endepoch));              % add in the first row!

% make a pixel non-NaN mask and total time image
n=zeros(nifgs,1);                                                      % fraction of coherent pixels
nz=zeros(nifgs,1);                                                     % fraction of nonzero pixels
N=zeros(nrows,ncols);                                                  % always coherent (non-NaN) mask
T=zeros(nrows,ncols);                                                  % total timespan
n(1)=sum(sum(isnan(maskunwphase)==0))/nrows/ncols;
nz(1)=sum(sum(maskunwphase~=0))/nrows/ncols;
N=N+(isnan(maskunwphase)==0);
T=T+(isnan(maskunwphase)==0)*(date2n(1)-date1n(1));

% let's make an informative output file
finfo=fopen('interferogram_info.out','w');
fprintf(finfo,'%s_%s %d %f %f \n',inlist{1}{1},inlist{2}{1},date2n-date1n,n(1),nz(1));

% and finally, make a 3d matrix to put the sbas results in
sbasinc=zeros(nrows,ncols,nepochs);                                      % incremental displacements
sbascum=zeros(nrows,ncols,nepochs);                                      % cumulative displacements


% loop from the second image to the last

for i=2:nifgs
   
   % when is this? 
   date1=ifgdates(i,1); 
   date2=ifgdates(i,2);
   
   % construct the file names
   unwfile=sprintf('%d_%d/%s',date1,date2,unwgrd);
   corfile=sprintf('%d_%d/%s',date1,date2,corgrd);


   % make an appropriate entry in the design matrix
   if date1n(i)==t0n     % if this starts on the 'zero' date
       startepoch=1;  % then set the start to the first date
   else
       startepoch=(find(tn==date1n(i)))+1;    % displacements are incremental
   end
   endepoch=find(tn==date2n(i));
   A(i,startepoch:endepoch)=ones(size(startepoch:endepoch));
   tspan(i)=date2n(i)-date1n(i);
   t1t2(i,:)=[date1n(i) date2n(i)];

 
   % read in moar data!!!
   fprintf(1,'opening %s\n',unwfile);
   %unwphase=ncread(unwfile,'/z')'; % interferometric phase
   %phsigcoh=ncread(corfile,'/z')'; % correlation (GMTSAR convention)
   [lon,lat,unwphase]=grdread2(unwfile);
   [lon,lat,phsigcoh]=grdread2(corfile);

   
   %maskunwphase=unwphase.*(conncomp>0);
   %maskunwphase=unwphase.*(phsigcoh>cohthresh);
   maskunwphase=unwphase;
   maskunwphase(maskunwphase==0)=NaN;

   refpha=nanmean(nanmean(maskunwphase(refrow-refshft:refrow+refshft,...
      refcol-refshft:refcol+refshft)));
   fprintf('  reference point phase: %f\n',refpha);
   refphase(i)=refpha;
   ifgdata(1:nrows,1:ncols,i)=(maskunwphase-refpha)*signconv;

   % add to yer mask and image
   n(i)=sum(sum(isnan(maskunwphase)==0))/nrows/ncols;
   nz(i)=sum(sum(maskunwphase~=0))/nrows/ncols;
   N=N.*(isnan(maskunwphase)==0);
   T=T+(isnan(maskunwphase)==0)*(date2n(i)-date1n(i));
   
   % append to yer info file
   fprintf(finfo,'%s_%s %d %f %f \n',inlist{1}{i},inlist{2}{i},date2n-date1n,n(i),nz(i));

   
end

fclose(finfo);

disp('Done loading!');