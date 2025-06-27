import numpy as np
import os
import sys
import pdb
import time
import argparse
from bgen import BgenReader
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV, Lasso
from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.utils.validation import check_is_fitted





class NonZeroLassoCV_OLD(BaseEstimator, RegressorMixin):
	"""
	LassoCV that refuses to choose an alpha whose final model is all-zero.

	Parameters
	----------
	alphas : array-like, shape (n_alphas,), optional
		Grid of α’s to search.  If None, uses the same default grid as LassoCV.
	cv : int or cross-validation generator, default=5
	min_nonzero : int, default=1
		Minimum number of non-zero coefficients required for an α to be eligible.
	**lasso_kwargs
		Extra keyword arguments forwarded to both LassoCV and Lasso
		(e.g. `fit_intercept`, `max_iter`, `tol`, …).
	"""
	def __init__(self, alphas=None, cv=5, min_nonzero=1, **lasso_kwargs):
		self.alphas        = alphas
		self.cv            = cv
		self.min_nonzero   = min_nonzero
		self.lasso_kwargs  = lasso_kwargs

	# --------------------------------------------------------------------- #
	def fit(self, X, y):
		# 1) ordinary LassoCV to get CV errors for every α
		self._lasso_cv_ = LassoCV(alphas=self.alphas,
								  cv=self.cv,
								  **self.lasso_kwargs).fit(X, y)
		#print(np.sum(self._lasso_cv_.coef_))
		#print(self._lasso_cv_.alpha_)

		alphas     = self._lasso_cv_.alphas_
		mse_mean   = self._lasso_cv_.mse_path_.mean(axis=1)

		# 2) run a brief Lasso fit at each α on the full data
		nz_mask = []
		for a in alphas:
			coef = Lasso(alpha=a, **self.lasso_kwargs).fit(X, y).coef_
			nz_mask.append(np.count_nonzero(coef) >= self.min_nonzero)
		nz_mask = np.array(nz_mask)

		# 3) select the α with *lowest CV error* among those that passed
		if nz_mask.any():
			eligible_mse      = mse_mean.copy()
			eligible_mse[~nz_mask] = np.inf          # make them ineligible
			best_idx          = int(np.argmin(eligible_mse))
		else:                  # nobody passed → fall back to LassoCV’s choice
			best_idx = int(np.argmin(mse_mean))

		self.alpha_     = float(alphas[best_idx])

		# 4) final refit with the chosen α
		self._model_    = Lasso(alpha=self.alpha_, **self.lasso_kwargs).fit(X, y)
		self.coef_      = self._model_.coef_
		self.intercept_ = self._model_.intercept_
		return self

	# --------------------------------------------------------------------- #
	def predict(self, X):
		check_is_fitted(self, "_model_")
		return self._model_.predict(X)

	# compatible score() delegates to the underlying Lasso
	def score(self, X, y):
		check_is_fitted(self, "_model_")
		return self._model_.score(X, y)

	# convenience: same API as LassoCV
	@property
	def n_iter_(self):
		check_is_fitted(self, "_model_")
		return self._model_.n_iter_

class NonZeroLassoCV(BaseEstimator, RegressorMixin):
    """
    LassoCV that refuses to choose an alpha whose final model is all-zero.

    Parameters
    ----------
    alphas : array-like, shape (n_alphas,), optional
        Grid of α’s to search.  If None, uses the same default grid as LassoCV.
    cv : int or cross-validation generator, default=5
    min_nonzero : int, default=1
        Minimum number of non-zero coefficients required for an α to be eligible.
    **lasso_kwargs
        Extra keyword arguments forwarded to both LassoCV and Lasso
        (e.g. `fit_intercept`, `max_iter`, `tol`, …).
    """
    def __init__(self, alphas=None, cv=5, min_nonzero=1, **lasso_kwargs):
        self.alphas        = alphas
        self.cv            = cv
        self.min_nonzero   = min_nonzero
        self.lasso_kwargs  = lasso_kwargs

    # --------------------------------------------------------------------- #
    def fit(self, X, y):
        # 1) ordinary LassoCV to get CV errors for every α
        self._lasso_cv_ = LassoCV(alphas=self.alphas,
                                  cv=self.cv,
                                  **self.lasso_kwargs).fit(X, y)

        alphas     = self._lasso_cv_.alphas_
        mse_mean   = self._lasso_cv_.mse_path_.mean(axis=1)

        # Non-zero coef off the bat
        if np.sum(self._lasso_cv_.coef_ != 0.0) > 0:
            return self._lasso_cv_, True
        else:
            for alpha_index in np.argsort(mse_mean):
                alpha = alphas[alpha_index]
                lasso_fit = Lasso(alpha=alpha, **self.lasso_kwargs).fit(X, y)
                lasso_coef = lasso_fit.coef_
                if np.count_nonzero(lasso_coef) >= self.min_nonzero:
                    break
            self.alpha_ = float(alpha)
            self.coef_ = np.copy(lasso_coef)
            return self, False




def load_in_ref_alt_allele_arr(pvar_file):
	ref_alt_alleles = []
	f = open(pvar_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 5:
			print('assumptino erooror')
			pdb.set_trace()
		ref_allele = data[3]
		alt_allele = data[4]
		if ref_allele == alt_allele:
			print('assumptino eororor')
			pdb.set_trace()
		ref_alt_alleles.append((ref_allele, alt_allele))
	f.close()
	return ref_alt_alleles

def extract_array_of_chromosomes(chromosome_file):
	if chromosome_file == 'None':
		chrom_arr = np.arange(1,23)
	else:
		chrom_arr = []
		f = open(chromosome_file)
		for line in f:
			line = line.rstrip()
			chrom_arr.append(int(line))
		chrom_arr = np.asarray(chrom_arr)
		f.close()
	return chrom_arr


def extract_expression_file_sample_names(expr_file):
	f = open(expr_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		sample_names = np.asarray(data[3:])
		# Only need first line
		break
	f.close()
	return sample_names

def extract_indices_corresponding_to_genotyped_samples(genotype_sample_names, expr_sample_names):
	dicti1 = {}
	for ii,val in enumerate(genotype_sample_names):
		dicti1[val.split('_')[0]] = ii

	mapping = []
	for expr_sample_name in expr_sample_names:
		mapping.append(dicti1[expr_sample_name])

	return mapping


def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices):
	dosages = []

	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.alt_dosage
		#ma = var.minor_allele

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


def create_mapping_from_rsid_to_snp_index(G_obj_rsids):
	dicti = {}
	for ii, val in enumerate(G_obj_rsids):
		dicti[val] = ii
	return dicti

def get_non_zero_beta(G_cis_std, expr_vec, alpha_init, ens_id):
	factor_grid = [.5, .1, .01, .001, .0001]

	booler = False
	for factor in factor_grid:
		cur_alpha = alpha_init*factor
		lasso_obj = Lasso(alpha=cur_alpha, max_iter=100000).fit(G_cis_std, expr_vec)
		std_beta_hat = lasso_obj.coef_
		non_sparsity = np.sum(std_beta_hat!=0)/len(std_beta_hat)	
		if non_sparsity > 0.0:
			booler = True
			break
	if booler == False:
		print('Zero-gene ' + ens_id)

	return std_beta_hat, cur_alpha


#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()

parser.add_argument('--expr', default='None', type=str,
					help='Expression file')
parser.add_argument('--bgen', default='None', type=str,
					help='Bgen file')
parser.add_argument('--chromosome-file', default='None', type=str,
					help='Filename containing list of chromosomes to run (no header). If None, then defaults to 22 chromosomes')
parser.add_argument('--alpha', default=0.01, type=float,
					help='Lasso regularization parameter')
parser.add_argument('--cis-window', default=500000.0, type=float,
					help='BP region around gene to call cis-eqtls')
parser.add_argument('--cross-val-alpha', default=False, action='store_true',
					help='Boolean on whether or not to use cross validation to get alpha for a given gene')
parser.add_argument('--output', default=None, type=str,
					help='Output file stem to save data to')
args = parser.parse_args()
np.random.seed(1)


# Open output file handle
t = open(args.output,'w')
# Print header
t.write('GENE\tCHR\tSNP\tSNP_POS\tA1\tA2\tDOSAGE_EFFECT\tSPARSITY_PARAM\tNON_SPARSITY_LEVEL\tCV_boolean\n')


# Extract array of chromosomes
chromosome_arr = extract_array_of_chromosomes(args.chromosome_file)

# Extract sample names
expr_sample_names = extract_expression_file_sample_names(args.expr)


counter = 0
aa = []
t_old = time.time()
# Loop through chromosomes
for chrom_num in chromosome_arr:


	# Load in genotype data for this chromosome
	genotype_stem = args.bgen + str(chrom_num) 
	genotype_obj = BgenReader(genotype_stem + '.bgen')
	G_obj_pos = np.asarray(genotype_obj.positions())
	G_obj_rsids = np.asarray(genotype_obj.rsids())
	snp_integers = np.arange(len(G_obj_rsids))
	
	# Create mapping from rsid to snp index
	rsid_to_snp_index = create_mapping_from_rsid_to_snp_index(G_obj_rsids)


	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')
	ref_alt_alleles = np.asarray(ref_alt_alleles)

	# Get indices corresponding to genotyped samples
	genotype_sample_indices = extract_indices_corresponding_to_genotyped_samples(genotype_obj.samples, expr_sample_names)


	# Now loop through expression file (one gene at a time)
	f = open(args.expr)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()

		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields
		ens_id = data[0]
		line_chrom_num = int(data[1])
		gene_position = int(data[2])

		# Skip genes not on current chromosome
		if line_chrom_num != chrom_num:
			continue
		counter = counter + 1

		if np.mod(counter,100) == 0:
			print(counter)
			t_cur = time.time()
			print(t_cur-t_old)
			t_old = t_cur

		# Extract vector of gene expression 
		expr_vec = np.asarray(data[3:]).astype(float)

		# Extract cis snp matrix
		cis_start_pos = gene_position - args.cis_window
		cis_end_pos = gene_position + args.cis_window
		cis_snp_indices = (G_obj_pos >= cis_start_pos) & (G_obj_pos <= cis_end_pos)
		cis_rsids = G_obj_rsids[cis_snp_indices]
		cis_pos = G_obj_pos[cis_snp_indices]
		cis_ref_alt_alleles = ref_alt_alleles[cis_snp_indices, :]
		G_cis_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, snp_integers[cis_snp_indices])
		# Filter to relevent individuals
		G_cis_dosage = G_cis_dosage[genotype_sample_indices, :]

		# Center Genotype dosage
		col_means = G_cis_dosage.mean(axis=0)
		G_cis_centered = G_cis_dosage - col_means

		# Standardize genotype
		snp_sdevs = G_cis_centered.std(axis=0)
		# Deal with nans if they exist
		if np.sum(snp_sdevs == 0.0) > 0.0:
			valid_indices = snp_sdevs != 0.0
			G_cis_centered = G_cis_centered[:, valid_indices]
			G_cis_dosage = G_cis_dosage[:, valid_indices]
			cis_rsids = cis_rsids[valid_indices]
			cis_pos = cis_pos[valid_indices]
			cis_ref_alt_alleles = cis_ref_alt_alleles[valid_indices]
			snp_sdevs = snp_sdevs[valid_indices]
		# Do the standardization
		G_cis_std = G_cis_centered/snp_sdevs


		# Fit lasso model
		#print('#####')
		if args.cross_val_alpha:
			lasso_obj, booler = NonZeroLassoCV(cv=5,alphas=np.arange(.01, .5, .02), max_iter=50000).fit(G_cis_std, expr_vec)
			std_beta_hat = lasso_obj.coef_
			sparsity_param = np.copy(lasso_obj.alpha_)*1.0
		else:
			lasso_obj = Lasso(alpha=args.alpha, max_iter=100000).fit(G_cis_std, expr_vec)
			std_beta_hat = lasso_obj.coef_
			sparsity_param = np.copy(args.alpha)*1.0
			booler = 'NA'
		#print(np.mean(aa))

		non_sparsity_level = np.sum(std_beta_hat!=0)/len(std_beta_hat)
		if non_sparsity_level == 0.0:
			std_beta_hat, sparsity_param = get_non_zero_beta(G_cis_std, expr_vec, args.alpha, ens_id)
			non_sparsity_level = np.sum(std_beta_hat!=0)/len(std_beta_hat)

		# Convert to dosage based beta-hats
		dosage_beta_hat = np.copy(std_beta_hat)/snp_sdevs

		# Print to output
		for snp_iter, snp_rsid in enumerate(cis_rsids):
			# Skip snps that have no effect
			tmp_beta = dosage_beta_hat[snp_iter]
			if tmp_beta == 0.0:
				continue
			# Extract relevent snp information
			snp_pos = cis_pos[snp_iter]
			snp_a1 = cis_ref_alt_alleles[snp_iter, 1]
			snp_a2 = cis_ref_alt_alleles[snp_iter, 0]
			# print
			t.write(ens_id + '\t' + str(chrom_num) + '\t' + snp_rsid + '\t' + str(snp_pos) + '\t' + snp_a1 + '\t' + snp_a2 + '\t' + str(tmp_beta) + '\t' + str(sparsity_param) + '\t' + str(non_sparsity_level) + '\t' + str(booler) + '\n')
	f.close()
t.close()





















