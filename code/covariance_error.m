function err=covariance_error(samples, empCov, burnin, thin)
  R=length(samples);
  if nargin<4
    thin=2;
  end
  if nargin<3
    burnin=R/5;
  end
 
  sampleseq=(burnin+1):thin:R;
  samples=samples(sampleseq);
  S=length(samples);
  mse=zeros(1,S); 
  medse=zeros(1,S);
  mae=zeros(1,S);
  medae=zeros(1,S);
  rmse=zeros(1,S);
  nrmse=zeros(1,S);
  postCov=zeros([size(empCov),S]);
  Covrange=max(empCov(:))-min(empCov(:));

  for r=1:S
    Lambda=samples{r}.Z.*samples{r}.G;
    postCov(:,:,r)=Lambda*Lambda'+diag(samples{r}.lambdae.^-1);
    diff = empCov - postCov(:,:,r);
    sqdiff = diff .* diff;
    absdiff = abs(diff);
    mse(r) = mean(sqdiff(:));
    medse(r) = median(sqdiff(:));
    mae(r) = mean(absdiff(:));
    medae(r) = median(absdiff(:));
    rmse(r) = sqrt(mse(r));
    nrmse(r) = rmse(r)/Covrange;
  end

  err.medians=[median(mse),median(medse),median(mae),...
             median(medae),median(rmse),median(nrmse)];

  err.last=[mse(end),medse(end),mae(end),medae(end),rmse(end),nrmse(end)];

  postCov_mean=mean(postCov,3);
  diff = empCov - postCov_mean;
  sqdiff = diff .* diff;
  absdiff = abs(diff);
  msepm = mean(sqdiff(:));
  medsepm = median(sqdiff(:));
  maepm = mean(absdiff(:));
  medaepm = median(absdiff(:));
  rmsepm = sqrt(msepm);
  nrmsepm = rmsepm/Covrange;

  err.postmean=[msepm,medsepm,maepm,medaepm,rmsepm,nrmsepm];
end 
