function [Fro,rc]=pred_check(samples, breaks, datCounts, burnin, thin)
  R=length(samples);
  N=size(samples{1}.X,2);
  D=size(samples{1}.G,1);
  maxbins = size(datCounts,2);
  if nargin<6
    thin=2;
  end
  if nargin<5
    burnin=R/5;
  end
  sampleseq=(burnin+1):thin:R;
  samples=samples(sampleseq);
  S=length(sampleseq);
  Fro=zeros(1,S);
  rc=cell(1,S);

  %for each posterior sample
  for r=1:S
    %sample from predictive distribution
    Ypred=(samples{r}.Z.*samples{r}.G)*samples{r}.X + mvnrnd(zeros(N,D),diag(samples{r}.lambdae.^-1))';
    %get histograms using breaks for each variable with araryfun and pad with zeros as in datCounts.
    repCounts=arrayfun(@(vind) [histc(Ypred(vind,:),breaks{vind})(1:end-1),zeros(1,maxbins-length(breaks{vind})+1)],1:D,'UniformOutput',false);
    %convert from cell to matrix
    repCounts=reshape(cell2mat(repCounts),maxbins,D)';
    %compute error
    datCounts=double(datCounts);
    repCounts=double(repCounts);
    Frob=norm(datCounts-repCounts,'fro');
    rFro=norm(repCounts,'fro');
    dFro=norm(datCounts,'fro');
    minFro=abs(dFro-rFro);
    Fro(r)=(Frob-minFro)/(dFro+rFro-minFro);
    rc{r}=repCounts;
  end
end
