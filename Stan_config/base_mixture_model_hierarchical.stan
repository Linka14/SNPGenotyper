// Hierarchical base mixture model with one Stan call for all samples and loci combined.
// Input are allele counts per sample and locus pair (k = allele1 reads, d = all usable reads).
// - error rate per locus and per sample, partially pooled
// - beta-binomial overdispersion: phi_het per sample, phi_hom shared or per locus and/or sample
// - allele bias per locus (heterozygote allele1 fraction can differ from 0.5)
// - genotype prior: fixed per pair (pair_prior) or learned wolf ancestry per sample (use_ancestry = 1)

data {
  int<lower=1> N; // sample and locus pairs with enough reads
  int<lower=1> S; // samples
  int<lower=1> L; // loci
  array[N] int<lower=1, upper=S> sample_id;
  array[N] int<lower=1, upper=L> locus_id;
  array[N] int<lower=0> k; // reads for allele1
  array[N] int<lower=1> d; // all usable reads (allele1 + allele2)
  real<lower=0, upper=0.5> max_error; // upper cap on any error rate (0.10)
  real<lower=0, upper=0.5> het_min; // HT allele1 fraction cap [het_min, 1 - het_min] (0.2)
  array[N] simplex[3] pair_prior; // genotype prior of every pair: P(hom1), P(hom2), P(het)
  int<lower=0, upper=1> use_ancestry; // 1 = genotype prior from the learned wolf ancestry of each sample
  vector<lower=0, upper=1>[L] p_wolf; // allele1 frequency in wolves, floored (only used when use_ancestry = 1)
  vector<lower=0, upper=1>[L] p_dog; // allele1 frequency in dogs, floored (only used when use_ancestry = 1)
  real<lower=0> ancestry_prior_shape; // Beta(shape, shape) prior on wolf ancestry, 1 is uniform, 0.5 leans to pure wolf or pure dog
  int<lower=0, upper=1> phi_hom_by_locus; // 1 = HM overdispersion differs between loci (partially pooled)
  int<lower=0, upper=1> phi_hom_by_sample; // 1 = HM overdispersion differs between samples (partially pooled)
}

parameters { // error rate: logit scale, typical value + locus effect + sample effect (non-centred)
  real mu_error;
  real<lower=0> sigma_error_locus;
  real<lower=0> sigma_error_sample;
  vector[L] z_error_locus;
  vector[S] z_error_sample;

  // allele bias: HT allele1 fraction on logit scale, typical value + locus effect
  real mu_bias;
  real<lower=0> sigma_bias;
  vector[L] z_bias;

  // overdispersion: concentration of the beta-binomial (large = close to binomial)
  real mu_log_phi_het; // typical sample: phi_het = 10 + exp(mu_log_phi_het)
  real<lower=0> sigma_phi_het; // how much samples differ (log scale)
  vector[S] z_phi_het; // sample effects (non-centred)
  real mu_log_phi_hom; // typical pair: phi_hom = 10 + exp(mu_log_phi_hom)
  real<lower=0> sigma_phi_hom_locus; // how much loci differ (log scale)
  vector[phi_hom_by_locus ? L : 0] z_phi_hom_locus;
  real<lower=0> sigma_phi_hom_sample; // how much samples differ (log scale)
  vector[phi_hom_by_sample ? S : 0] z_phi_hom_sample;

  // wolf ancestry of every sample, only exists when use_ancestry = 1
  vector<lower=0, upper=1>[use_ancestry ? S : 0] wolf_ancestry;
}

transformed parameters {
  // expected allele1 share in a HT, kept between het_min and 1 - het_min
  vector[L] het_fraction = het_min + (1 - 2 * het_min) * inv_logit(mu_bias + sigma_bias * z_bias);
  // HT overdispersion of every sample, at least ~10 effective molecules
  vector[S] phi_het_sample = 10 + exp(mu_log_phi_het + sigma_phi_het * z_phi_het);
  vector[N] phi_hom_pair; // HM overdispersion of every pair: typical + locus effect + sample effect
  for (n in 1:N) {
    real log_phi = mu_log_phi_hom;
    if (phi_hom_by_locus) log_phi += sigma_phi_hom_locus * z_phi_hom_locus[locus_id[n]];
    if (phi_hom_by_sample) log_phi += sigma_phi_hom_sample * z_phi_hom_sample[sample_id[n]];
    phi_hom_pair[n] = 10 + exp(log_phi);
  }
  vector[N] error_rate; // expected wrong-allele share in a HM, capped at max_error
  for (n in 1:N) {
    error_rate[n] = max_error * inv_logit(mu_error + sigma_error_locus * z_error_locus[locus_id[n]]
                                                   + sigma_error_sample * z_error_sample[sample_id[n]]);
  }
  // genotype prior of every pair on the log scale: hom1, hom2, het
  array[N] vector[3] log_genotype_prior;
  for (n in 1:N) {
    if (use_ancestry) {
      real p = wolf_ancestry[sample_id[n]] * p_wolf[locus_id[n]] + (1 - wolf_ancestry[sample_id[n]]) * p_dog[locus_id[n]];
      log_genotype_prior[n] = [log(square(p)), log(square(1 - p)), log(2 * p * (1 - p))]';   // Hardy-Weinberg on the sample's own allele frequency
    } else {
      log_genotype_prior[n] = log(pair_prior[n]);
    }
  }
}

model { // hyperpriors: mild, the data overrule them quickly
  mu_error ~ normal(-2.3, 1); // typical error around 0.10 * inv_logit(-2.3) = 0.9 %
  sigma_error_locus ~ normal(0, 1);
  sigma_error_sample ~ normal(0, 1);
  z_error_locus ~ std_normal();
  z_error_sample ~ std_normal();
  mu_bias ~ normal(0, 0.5); // typical HT centred on 50/50
  sigma_bias ~ normal(0, 0.25); // loci differ only mildly in allele bias
  z_bias ~ std_normal();
  mu_log_phi_het ~ normal(log(20), 1); // typical sample is 10 + 20 = 30 effective molecules, as the earlier beta-binomial fit
  sigma_phi_het ~ normal(0, 0.5); // samples differ mildly, data can make it larger
  z_phi_het ~ std_normal();
  mu_log_phi_hom ~ normal(log(200), 1.5); // same centre and width as the old lognormal(log(200), 1.5) prior
  sigma_phi_hom_locus ~ normal(0, 1);
  sigma_phi_hom_sample ~ normal(0, 1);
  if (phi_hom_by_locus) z_phi_hom_locus ~ std_normal();
  if (phi_hom_by_sample) z_phi_hom_sample ~ std_normal();
  if (use_ancestry) wolf_ancestry ~ beta(ancestry_prior_shape, ancestry_prior_shape); // one per sample, no shared mean, so a rare hybrid is not pulled toward the wolves

  for (n in 1:N) {
    real e = error_rate[n];
    real h = het_fraction[locus_id[n]];
    real phi_het_n = phi_het_sample[sample_id[n]];
    real phi_hom = phi_hom_pair[n];
    real p_het = h * (1 - e) + (1 - h) * e; // read errors also apply inside a HT
    vector[3] lp;
    lp[1] = log_genotype_prior[n][1] + beta_binomial_lpmf(k[n] | d[n], (1 - e) * phi_hom, e * phi_hom); // hom1
    lp[2] = log_genotype_prior[n][2] + beta_binomial_lpmf(k[n] | d[n], e * phi_hom, (1 - e) * phi_hom); // hom2
    lp[3] = log_genotype_prior[n][3] + beta_binomial_lpmf(k[n] | d[n], p_het * phi_het_n, (1 - p_het) * phi_het_n); // het
    target += log_sum_exp(lp);
  }
}

generated quantities {
  array[N] vector[3] posterior;   // per pair: P(hom1), P(hom2), P(het)
  vector[L] error_locus = max_error * inv_logit(mu_error + sigma_error_locus * z_error_locus); // locus error for an average sample
  vector[S] error_sample = max_error * inv_logit(mu_error + sigma_error_sample * z_error_sample); // sample error for an average locus
  real phi_het = 10 + exp(mu_log_phi_het); // phi_het of a typical sample
  real phi_hom_typical = 10 + exp(mu_log_phi_hom); // phi_hom of a typical pair
  for (n in 1:N) {
    real e = error_rate[n];
    real h = het_fraction[locus_id[n]];
    real phi_het_n = phi_het_sample[sample_id[n]];
    real phi_hom = phi_hom_pair[n];
    real p_het = h * (1 - e) + (1 - h) * e;
    vector[3] lp;
    lp[1] = log_genotype_prior[n][1] + beta_binomial_lpmf(k[n] | d[n], (1 - e) * phi_hom, e * phi_hom);
    lp[2] = log_genotype_prior[n][2] + beta_binomial_lpmf(k[n] | d[n], e * phi_hom, (1 - e) * phi_hom);
    lp[3] = log_genotype_prior[n][3] + beta_binomial_lpmf(k[n] | d[n], p_het * phi_het_n, (1 - p_het) * phi_het_n);
    posterior[n] = softmax(lp);   // prior x likelihood, normalised (with a flat prior this is the normalised likelihood)
  }
}