function [totalpower,pow,F,lvar,ntrials]=getpower3_2(M,fs)
% return power from each trace
[ntrials,m]=size(M);

pow=[];
NW= 2.5;%10;%5%2.5%1.5;%2.5;
ntapers=2*NW-1;

for i=1:ntrials,

  
  data=M(i,find(~isnan(M(i,:))));
  
%  [P,F]=pmtm(data(1,:),NW,[],fs);


% the following does the same as [P,F]=pmtm(data(1,:),NW,[],fs);
% but it is copied here to keep info and do bootstrap on tapers
% to compute error bars
  N=length(data);
  nfft = max(256,2^nextpow2(N));
%   nfft=max(2^(nextpow2(N)+pad),N);

  data=data(:);
  k = min(round(2*NW),N);
  k = max(k-1,1);
  k = min(k,ntapers);
  [E,V]=dpss(N,NW,k);%tapers
  
  % Compute the windowed DFTs.
  Sk=abs(fft(E.*data(:,ones(1,k)),nfft)).^2;

  wt = ones(k,1);
  S = Sk*wt/k;

  w = 2*pi*(0:1/nfft:1-1/nfft);
  w=w(1:nfft/2+1);
  S=S(1:nfft/2+1);
  S=S(:);
  P=[S(1); 2*S(2:end-1); S(end)]./fs; %converts to single-sided and
%   normalize; Xq comment out 2017 11 22
  F=w.*fs./(2.*pi);
  F=F(:);
%    P=[S(1); 2*S(2:end-1); S(end)]./F;% XQ 2017 11/22
%   %up to here same as pmtm.m
  
  pow(i,:)=P';
end

% freq60=0;
% for i=1:length(F),
%   if freq60==0&F(i)>60, freq60=i; end
% end 

% I=1;
% while ~isempty(I),
%   clear lsp
  for i=1:ntrials,
    indices=setdiff(1:ntrials,i);
    lsp(i,:)=log10(1/(ntrials-1)*sum(pow(indices,:),1));
  end
  % now remove any single trial that when excluded modifies the
  % 60 Hz power by more than 3 standard deviations
%   st=std(lsp(:,freq60));
%   mn=mean(lsp(:,freq60));      
%   I=find(abs(lsp(:,freq60)-mn)>3*st);
%   indxs=setdiff(1:ntrials,I);
%   pow=pow(indxs,:);
%   [ntrials,m]=size(pow);
%   % continue until no trial affects on its own the power at 60 Hz 
%   % by more than 3 std's
% end


totalpower=mean(pow,1);
lvar=sqrt((ntrials-1)/2)*std(lsp,0,1); %jacknife std over trials (((ntrials-1)/2)*std(lsp,0,1)^2)

alpha=0.95;

jkCI=tinv((1-alpha)/2,ntrials-1)*lvar;

%dof=2*ntrials*k;

%theo1=chi2conf(alpha,dof);

%theo2=2*sqrt(2*(1/dof+1/dof^2));

%semilogy(F,totalpower,'b',F,totalpower.*10.^lvar,'r',F,totalpower.*10.^(- ...
%						  jkCI),'g',F,totalpower* ...
%	 theo1(2),'k',F,totalpower*exp(theo2),'m')

%keyboard

lvar=-jkCI;