library(rstan)
library(HDInterval)
library(DescTools)

#PA = Plate assay performed at the Ifremer quarantine lab in La Tremblade, France
#TB = Field trial performed in Tomales Bay, California

Data_PA = read.csv("Data_PA.csv")
Data_TB = read.csv("Data_TB.csv")

#Additive pedigree relationship matrix
A = as.matrix(read.csv("A.csv",row.names = 1,check.names = F))

#Only keep genotypes that have a pedigree, i.e., remove controls
Data_PA = Data_PA[-which(is.na(match(Data_PA$Family,rownames(A)))),]
A_PA = A[match(Data_PA$Family,rownames(A)),match(Data_PA$Family,rownames(A))]

Data_TB = Data_TB[-which(is.na(match(Data_TB$Family,rownames(A)))),]
A_TB = A[match(Data_TB$Family,rownames(A)),match(Data_TB$Family,rownames(A))]

#CPS = Cumulative percent survival
#PS = Percent survival
#D = Day

#Data setup
PA_CPS_D3_data = list(y = Data_PA$CPS_D3,n = length(Data_PA$CPS_D3),A = A_PA)
PA_CPS_D4_data = list(y = Data_PA$CPS_D4,n = length(Data_PA$CPS_D4),A = A_PA)
PA_CPS_D5_data = list(y = Data_PA$CPS_D5,n = length(Data_PA$CPS_D5),A = A_PA)
PA_CPS_D6_data = list(y = Data_PA$CPS_D6,n = length(Data_PA$CPS_D6),A = A_PA)
PA_CPS_D7_data = list(y = Data_PA$CPS_D7,n = length(Data_PA$CPS_D7),A = A_PA)
TB_PS_data = list(y = Data_TB$PS,n = length(Data_TB$PS),A = A_TB)

#Gaussian process regression model
Model = stan_model(model_code = "
data {
  int n;
  vector[n] y;
  matrix[n,n] A;
}

parameters {
  real<lower=0,upper=variance(y)> sigma2_a;
  real<lower=0,upper=variance(y)> sigma2_e;
}
model {
  matrix[n,n] K;
  K = sigma2_a*A + diag_matrix(rep_vector(sigma2_e,n));
  y ~ multi_normal(rep_vector(mean(y),n),K);
}
")

#Run models to estimate genetic and error variances
PA_CPS_D3_data_fit = sampling(Model, data = PA_CPS_D3_data, iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99))
PA_CPS_D4_data_fit = sampling(Model, data = PA_CPS_D4_data, iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99))
PA_CPS_D5_data_fit = sampling(Model, data = PA_CPS_D5_data, iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99))
PA_CPS_D6_data_fit = sampling(Model, data = PA_CPS_D6_data, iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99))
PA_CPS_D7_data_fit = sampling(Model, data = PA_CPS_D7_data, iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99))
TB_PS_data_fit = sampling(Model, data = TB_PS_data, iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99))

#To exactly replicate the results in the paper (MCMC requires sampling), load the saved models below instead of running the code above
PA_CPS_D3_data_fit = readRDS("PA_CPS_D3_data_fit.rds")
PA_CPS_D4_data_fit = readRDS("PA_CPS_D4_data_fit.rds")
PA_CPS_D5_data_fit = readRDS("PA_CPS_D5_data_fit.rds")
PA_CPS_D6_data_fit = readRDS("PA_CPS_D6_data_fit.rds")
PA_CPS_D7_data_fit = readRDS("PA_CPS_D7_data_fit.rds")
TB_PS_data_fit = readRDS("TB_PS_data_fit.rds")


#Extract genetic and error variance samples of all models
#Calculate heritability for the phenotypes from these samples

PA_CPS_D3_sigma2_a = extract(PA_CPS_D3_data_fit)$sigma2_a
PA_CPS_D3_sigma2_e = extract(PA_CPS_D3_data_fit)$sigma2_e
PA_CPS_D3_h2 = PA_CPS_D3_sigma2_a/(PA_CPS_D3_sigma2_a+PA_CPS_D3_sigma2_e/2)

PA_CPS_D4_sigma2_a = extract(PA_CPS_D4_data_fit)$sigma2_a
PA_CPS_D4_sigma2_e = extract(PA_CPS_D4_data_fit)$sigma2_e
PA_CPS_D4_h2 = PA_CPS_D4_sigma2_a/(PA_CPS_D4_sigma2_a+PA_CPS_D4_sigma2_e/2)

PA_CPS_D5_sigma2_a = extract(PA_CPS_D5_data_fit)$sigma2_a
PA_CPS_D5_sigma2_e = extract(PA_CPS_D5_data_fit)$sigma2_e
PA_CPS_D5_h2 = PA_CPS_D5_sigma2_a/(PA_CPS_D5_sigma2_a+PA_CPS_D5_sigma2_e/2)

PA_CPS_D6_sigma2_a = extract(PA_CPS_D6_data_fit)$sigma2_a
PA_CPS_D6_sigma2_e = extract(PA_CPS_D6_data_fit)$sigma2_e
PA_CPS_D6_h2 = PA_CPS_D6_sigma2_a/(PA_CPS_D6_sigma2_a+PA_CPS_D6_sigma2_e/2)

PA_CPS_D7_sigma2_a = extract(PA_CPS_D7_data_fit)$sigma2_a
PA_CPS_D7_sigma2_e = extract(PA_CPS_D7_data_fit)$sigma2_e
PA_CPS_D7_h2 = PA_CPS_D7_sigma2_a/(PA_CPS_D7_sigma2_a+PA_CPS_D7_sigma2_e/2)

TB_PS_sigma2_a = extract(TB_PS_data_fit)$sigma2_a
TB_PS_sigma2_e = extract(TB_PS_data_fit)$sigma2_e
TB_PS_h2 = TB_PS_sigma2_a/(TB_PS_sigma2_a+TB_PS_sigma2_e/4)


#Functions to calculate the posterior mode and lower and upper 95% highest posterior density interval limits
PosteriorMode = function(vals,min,max){
  density = density(vals,from=min,to=max)
  return(density$x[density$y==max(density$y)])
}

lHDI = function(vals,min,max){
  density = density(vals,from=min,to=max)
  return(hdi(density)[1])
}

uHDI = function(vals,min,max){
  density = density(vals,from=min,to=max)
  return(hdi(density)[2])
}

#Posterior mode and lower and upper 95% highest posterior density interval for the plate assay day 3 phenotype's heritability
#Replace PA_CPS_D3_h2 with PA_CPS_D4_h2, PA_CPS_D5_h2, PA_CPS_D6_h2, PA_CPS_D7_h2, or TB_PS_h2 to get the same statistics for the other phenotypes
PosteriorMode(PA_CPS_D3_h2,0,1)
lHDI(PA_CPS_D3_h2,0,1)
uHDI(PA_CPS_D3_h2,0,1)

#Posterior mode and lower and upper 95% highest posterior density interval for the plate assay day 3 phenotype's genetic variance
PosteriorMode(PA_CPS_D3_sigma2_a,0,var(PA_CPS_D3_data$y))
lHDI(PA_CPS_D3_sigma2_a,0,var(PA_CPS_D3_data$y))
uHDI(PA_CPS_D3_sigma2_a,0,var(PA_CPS_D3_data$y))

#Posterior mode and lower and upper 95% highest posterior density interval for the plate assay day 3 phenotype's error variance
PosteriorMode(PA_CPS_D3_sigma2_e,0,var(PA_CPS_D3_data$y))
lHDI(PA_CPS_D3_sigma2_e,0,var(PA_CPS_D3_data$y))
uHDI(PA_CPS_D3_sigma2_e,0,var(PA_CPS_D3_data$y))


#Relationship matrices between all 71 MBP families and those in the two experiments
A_All_PA = A[,match(Data_PA$Family,rownames(A))]
A_All_TB = A[,match(Data_TB$Family,rownames(A))]

#Breeding value (BV) estimation

y = PA_CPS_D3_data$y
sigma2_a = PosteriorMode(PA_CPS_D3_sigma2_a,0,var(y))
sigma2_e = PosteriorMode(PA_CPS_D3_sigma2_e,0,var(y))
K = sigma2_a*A_PA + sigma2_e*diag(nrow(A_PA))
PA_CPS_D3_BV = mean(y) + sigma2_a*A_All_PA %*% solve(K) %*% (y-mean(y))

y = PA_CPS_D4_data$y
sigma2_a = PosteriorMode(PA_CPS_D4_sigma2_a,0,var(y))
sigma2_e = PosteriorMode(PA_CPS_D4_sigma2_e,0,var(y))
K = sigma2_a*A_PA + sigma2_e*diag(nrow(A_PA))
PA_CPS_D4_BV = mean(y) + sigma2_a*A_All_PA %*% solve(K) %*% (y-mean(y))

y = PA_CPS_D5_data$y
sigma2_a = PosteriorMode(PA_CPS_D5_sigma2_a,0,var(y))
sigma2_e = PosteriorMode(PA_CPS_D5_sigma2_e,0,var(y))
K = sigma2_a*A_PA + sigma2_e*diag(nrow(A_PA))
PA_CPS_D5_BV = mean(y) + sigma2_a*A_All_PA %*% solve(K) %*% (y-mean(y))

y = PA_CPS_D6_data$y
sigma2_a = PosteriorMode(PA_CPS_D6_sigma2_a,0,var(y))
sigma2_e = PosteriorMode(PA_CPS_D6_sigma2_e,0,var(y))
K = sigma2_a*A_PA + sigma2_e*diag(nrow(A_PA))
PA_CPS_D6_BV = mean(y) + sigma2_a*A_All_PA %*% solve(K) %*% (y-mean(y))

y = PA_CPS_D7_data$y
sigma2_a = PosteriorMode(PA_CPS_D7_sigma2_a,0,var(y))
sigma2_e = PosteriorMode(PA_CPS_D7_sigma2_e,0,var(y))
K = sigma2_a*A_PA + sigma2_e*diag(nrow(A_PA))
PA_CPS_D7_BV = mean(y) + sigma2_a*A_All_PA %*% solve(K) %*% (y-mean(y))

y = TB_PS_data$y
sigma2_a = PosteriorMode(TB_PS_sigma2_a,0,var(y))
sigma2_e = PosteriorMode(TB_PS_sigma2_e,0,var(y))
K = sigma2_a*A_TB + sigma2_e*diag(nrow(A_TB))
TB_PS_BV = mean(y) + sigma2_a*A_All_TB %*% solve(K) %*% (y-mean(y))

#Genetic correlation and 95% confidence interval between the plate assay day 3 and 4 phenotypes
#Exchange PA_CPS_D3_BV and PA_CPS_D4_BV for the breeding values estimated using other phenotypes to obtain the other genetic correlations
CorCI(cor(PA_CPS_D3_BV[,1],PA_CPS_D4_BV[,1]),71)
