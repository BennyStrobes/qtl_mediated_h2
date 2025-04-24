import numpy as np
import os
import sys
import pdb
import gzip


def get_anno_names(anno_file):
	if anno_file.endswith('gz'):
		f = gzip.open(anno_file)
	else:
		f = open(anno_file)

	for line in f:
		if anno_file.endswith('gz'):
			line = line.decode('utf-8').rstrip()
		else:
			line = line.rstrip()
		data = line.split('\t')
		anno_names = np.asarray(data[3:])
		break
	f.close()
	return anno_names


def extract_tau_and_tau_se_from_log_file(sldsc_log_file, anno_names):
	f = open(sldsc_log_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('Coefficients:'):
			coef = np.asarray(line.split()[1:]).astype(float)
			if len(coef) < len(anno_names):
				#line = f.next()
				line = next(f)
				data = np.asarray(line.split()).astype(float)
				coef = np.hstack([coef, data])
		if line.startswith('Coefficient SE:'):
			coef_se = np.asarray(line.split()[2:]).astype(float)
			if len(coef_se) < len(anno_names):
				#line = f.next()
				line = next(f)
				data = np.asarray(line.split()).astype(float)
				coef_se = np.hstack([coef_se, data])
	return coef, coef_se



def print_organized_summary_file(output_file_name, anno_names, tau, tau_se):
	tau_z = tau/tau_se
	t = open(output_file_name,'w')
	t.write('Annotation\ttau\ttau_se\ttau_z\n')
	for ii in range(len(anno_names)):
		t.write(anno_names[ii] + '\t' + str(tau[ii]) + '\t' + str(tau_se[ii]) + '\t' + str(tau_z[ii]) + '\n')
	t.close()

def load_in_mvec(anno_stem, suffix):
	for chrom_num in range(1,23):
		filer = anno_stem  + str(chrom_num) + '.l2' + suffix
		tmp = np.loadtxt(filer)
		if tmp.size == 1:
			tmp = np.asarray([tmp*1.0])
		if chrom_num == 1:
			counter = np.zeros(len(tmp))
		counter = counter + tmp
	return counter

def extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_vec):
	jacknifed_med_h2 = jacknifed_taus*m_vec


	jacknifed_geno_h2 = np.sum(jacknifed_med_h2[:,:eqtl_start_index],axis=1)
	jacknifed_expr_h2 = np.sum(jacknifed_med_h2[:,eqtl_start_index:],axis=1)

	jacknifed_expr_med_h2 = jacknifed_expr_h2/(jacknifed_expr_h2 + jacknifed_geno_h2)

	mean_expr_med_h2 = np.mean(jacknifed_expr_med_h2)

	diff = jacknifed_expr_med_h2 - mean_expr_med_h2
	num_jacknife_samples = jacknifed_expr_med_h2.shape[0]

	jacknifed_se = np.sqrt(np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples)
	
	return mean_expr_med_h2, jacknifed_se


def print_organized_h2_mediated(output_file_name, h2_med, h2_med_se):
	t = open(output_file_name,'w')
	t.write('h2_med\th2_med_se\n')
	t.write(str(h2_med) + '\t' + str(h2_med_se) + '\n')
	t.close()

#######################
# Command line args
#######################
tglr_output_stem = sys.argv[1]
variant_anno_stem = sys.argv[2]
expression_anno_stem = sys.argv[3]


example_variant_anno_file = variant_anno_stem + '21.l2.ldscore.gz'
variant_anno_names = get_anno_names(example_variant_anno_file)
example_gene_anno_file = expression_anno_stem + '21.l2.ldscore'
gene_anno_names = get_anno_names(example_gene_anno_file)

joint_anno_names = np.hstack((variant_anno_names, gene_anno_names))

tglr_log_file = tglr_output_stem + '.log'
tau, tau_se = extract_tau_and_tau_se_from_log_file(tglr_log_file, joint_anno_names)
tau_z = tau/tau_se


# Print to output
print_organized_summary_file(tglr_output_stem + '_organized_coef.txt', joint_anno_names, tau, tau_se)


# Total counts
variant_m_vec = load_in_mvec(variant_anno_stem, '.M')
variant_m_5_50_vec = load_in_mvec(variant_anno_stem, '.M_5_50')
gene_m_vec = load_in_mvec(expression_anno_stem, '.M')
gene_m_5_50_vec = load_in_mvec(expression_anno_stem, '.M_5_50')

m_vec = np.hstack((variant_m_vec, gene_m_vec))
m_5_50_vec = np.hstack((variant_m_5_50_vec, gene_m_5_50_vec))

# Print partioned h2 to output
print_organized_summary_file(tglr_output_stem + '_mediated_h2.txt', joint_anno_names, tau*m_vec, tau_se*m_vec)
print_organized_summary_file(tglr_output_stem + '_mediated_h2_5_50.txt', joint_anno_names, tau*m_5_50_vec, tau_se*m_5_50_vec)



# Load in jacknifed data
jacknifed_taus_file = tglr_output_stem + '.part_delete'
jacknifed_taus = np.loadtxt(jacknifed_taus_file)


# Jacknife expression mediated h2
eqtl_start_index = len(variant_anno_names)

h2_med, h2_med_se = extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_vec)
h2_5_50_med, h2_5_50_med_se = extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_5_50_vec)
print_organized_h2_mediated(tglr_output_stem + '_frac_h2_med.txt', h2_med, h2_med_se)
print_organized_h2_mediated(tglr_output_stem + '_frac_h2_med_5_50.txt', h2_5_50_med, h2_5_50_med_se)


# Delete unnecessary files
os.system('rm ' + tglr_output_stem + '.log')
os.system('rm ' + tglr_output_stem + '.delete')
os.system('rm ' + tglr_output_stem + '.part_delete')
os.system('rm ' + tglr_output_stem + '.intercept_delete')







