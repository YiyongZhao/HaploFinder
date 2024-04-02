import pandas as pd
from ete3 import PhyloTree
import sys
import random
from collections import defaultdict
import time
import subprocess

def generate_sps_voucher(sps_num:int) -> list:
	characters = [chr(i) for i in range(65, 91)] + [chr(i) for i in range(97, 123)] + [str(i) for i in range(10)]
	unique_strings = set()

	while len(unique_strings) < sps_num:
		unique_strings.add(''.join(random.sample(characters, 3)))

	return sorted(list(unique_strings))


def gene_id_transfer(gene2taxa_list:str) -> dict:
	gene2taxa_dic = read_and_return_dict(gene2taxa_list)
	taxa_list = list(set(gene2taxa_dic.values()))
	taxa2voucher_dic = dict(zip(taxa_list, generate_sps_voucher(len(taxa_list))))
	voucher2taxa_dic = {value: key for key, value in taxa2voucher_dic.items()}
	gene_count = {}

	for species in gene2taxa_dic.values():
		gene_count[species] = gene_count.get(species, 0) + 1

	new_gene_names = [f"{taxa2voucher_dic[species]}_{i}" for species, count in gene_count.items() for i in range(1, count + 1)]
	gene2new_named_gene_dic = dict(zip(gene2taxa_dic.keys(), new_gene_names))
	new_named_gene2gene_dic = {value: key for key, value in gene2new_named_gene_dic.items()}
	return gene2new_named_gene_dic,new_named_gene2gene_dic,taxa2voucher_dic
#gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("gene2taxa.list")

def read_and_return_dict(filename, separator="\t") -> dict:
	df=pd.read_csv(filename,sep=separator,header=None)
	return df.set_index([0])[1].to_dict()

def rename_input_tre(Phylo_t:object, gene2new_named_gene_dic:dict) -> object:
	Phylo_t1=Phylo_t.copy()
	for node in Phylo_t1.traverse():
		if node.name in gene2new_named_gene_dic.keys():
			node.name = gene2new_named_gene_dic[node.name]
	return Phylo_t1

def read_phylo_tree(tre_path:str) -> object:
	return PhyloTree(tre_path)
######################################################################################################################

def find_dup_node(Phylo_t:object)->list:#After searching for duplication events, the list of node names where duplication events occurred is as follows:
	events = Phylo_t.get_descendant_evol_events()
	dup_node_name_list = []
	for ev in events:
		if ev.etype == "D":
			i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
			events_node_name_list = i.split(',')
			common_ancestor_node_name = Phylo_t.get_common_ancestor(events_node_name_list)
			dup_node_name_list.append(common_ancestor_node_name)
	return dup_node_name_list


def generate_combinations(list1, list2):
	gene_pairs=[]
	for elem1 in list1:
		for elem2 in list2:
			gene_pairs.append((elem1,elem2))
	return gene_pairs

def filter_gd(gene_pair, paral_clade_list, new_named_gene2gene_dic):
	item1_set = set(gene_pair)

	for gd in paral_clade_list:
		if get_species_set(gd) ==set(taxa2voucher_dic.values()):
			gd_clade1, gd_clade2 = gd.get_children()
			gd_tips1 = [leaf for leaf in gd_clade1.get_leaf_names()]
			gd_tips2 = [leaf for leaf in gd_clade2.get_leaf_names()]
			set1 = set(gd_tips1)
			set2 = set(gd_tips2)
			total = set(gd_tips1 + gd_tips2)

			if item1_set <= total:
				if item1_set <= set1 or item1_set <= set2:
					return new_named_gene2gene_dic[gene_pair[0]] + '\t' + new_named_gene2gene_dic[gene_pair[1]] + '\t' + 'red'+'\n'
				else:
					return new_named_gene2gene_dic[gene_pair[0]] + '\t' + new_named_gene2gene_dic[gene_pair[1]] + '\t' + 'green'+'\n'

	return ''  # Return empty string if no match found

def get_species_list(node):
	if node is None:
		return []  # or any other appropriate default value
	return [leaf.name.split('_')[0] for leaf in node.iter_leaves()]

def get_species_set(Phylo_t:object)->set:
	return set(get_species_list(Phylo_t))

def get_gene_pairs(t, outgroup):
	sp_count = defaultdict(list)
	for i in t.get_leaf_names():
		sps = i.split('_')[0]
		if sps == outgroup:
			continue
		sp_count[sps].append(i)

	sp1_gene_list = list(sp_count.values())[0]
	sp2_gene_list = list(sp_count.values())[1]
	gene_pairs = generate_combinations(sp1_gene_list, sp2_gene_list)
	return gene_pairs

def main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic, taxa2voucher_dic, outgroup):
	with open('tree_label', 'w') as o:
		for tre_ID, tre_path in tre_dic.items():
			Phylo_t0 = read_phylo_tree(tre_path)
			Phylo_t0.resolve_polytomy(recursive=True)
			Phylo_t0.sort_descendants()
			Phylo_t1 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)
			if len(get_species_set(Phylo_t1)) == 1:
				continue
			sps_tol=get_species_set(Phylo_t1)

			if len(sps_tol)==2 and taxa2voucher_dic[outgroup] in sps_tol:
				continue

			paral_clade_list = find_dup_node(Phylo_t1)
			gene_pairs = get_gene_pairs(Phylo_t1, taxa2voucher_dic[outgroup])

			
			for item in gene_pairs:
				result = filter_gd(item, paral_clade_list, new_named_gene2gene_dic)
				if result:

					o.write(result)
					
		

if __name__ == "__main__":
	gf=sys.argv[1]
	imap=sys.argv[2]
	outgroup=sys.argv[3]
	t1=time.time()
	tre_dic=read_and_return_dict(gf)
	gene2new_named_gene_dic, new_named_gene2gene_dic ,taxa2voucher_dic= gene_id_transfer(imap)
	main(tre_dic,gene2new_named_gene_dic, new_named_gene2gene_dic,taxa2voucher_dic,outgroup)

	sort_command = f"sort -k3,3 -t$'\t' --stable -o tree_label tree_label"
	subprocess.run(sort_command, shell=True, check=True)
	t2=time.time()
	print("Label tree_gd took "+str(t2-t1)+" second") 


