function [aicm,bicm,dic]=model_selection_criteria(results,N,burnin,thin)
  R=size(results,1);
  if nargin<4
    thin=2;
  end
  if nargin<3
    burnin=R/5;
  end
  sampleseq=(burnin+1):thin:R;
  postll=results(sampleseq,9);
  lbar=mean(postll);
  dhat=2*var(postll);
  lmax=2*lbar+dhat;
  aicm=4*lbar-lmax;
  bicm=lmax-dhat*log(N);
  dic =2*(3*lbar-lmax);
end
