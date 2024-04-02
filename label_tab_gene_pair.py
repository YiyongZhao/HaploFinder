#Usage :python label_tab_gene_pair.py DAOD.tab_2 phyto73.lineage-70.txt wang.txt
import concurrent.futures
import sys
import subprocess
import time

def get_tab(tab_filename):
    with open(tab_filename,'r') as file :
        tab = []
        count = 0
        for line in file:
            if count == 0:
                count += 1
                continue
            elements = line.strip().split('\t')
            if len(elements) == 3:
                tab.append(tuple(elements))
            else:
                print(f"Error: Line '{line}' does not contain three elements.")
                continue 
    return tab

def get_gd(gd_pair):
    with open(gd_pair,'r') as file :
        gd = []
        count = 0
        for line in file:
            if count == 0:
                count += 1
                continue
            elements = line.strip().split()
            if len(elements) == 5:
                gd.append(tuple(elements[3:]))
            else:
                print(f"Error: Line '{line}' does not contain three elements.")
                continue    
    return gd

def filter_gd(i, gd):
    gene_pair = set(list(i)[1:])
    result = '\t'.join(list(i)) + '\t' + 'NONE'  

    for j in gd:
        clade1, clade2 = j
        lst1 = [i.split('_')[1] for i in clade1.split(',')]
        lst2 = [i.split('_')[1] for i in clade2.split(',')]
        set1 = set(lst1)
        set2 = set(lst2)
        total = set(lst1 + lst2)

        if gene_pair <= total:
            if gene_pair <= set1 or gene_pair <= set2:
                result = '\t'.join(list(i)) + '\t' + 'red'
                break  
                
            else:
                result = '\t'.join(list(i)) + '\t' + 'blue'
                break  

    return result

def process_data_in_threads(tab, gd,outfile,num_threads: int=4):
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = []
        for i in tab:
            results.append(executor.submit(filter_gd, i, gd))
        
        with open(outfile, 'w') as o:
            for future in concurrent.futures.as_completed(results):
                s = future.result()
                o.write(s + '\n')


if __name__=="__main__":

    tab_filename=sys.argv[1]
    gd_pair=sys.argv[2]
    outfile=sys.argv[3]
    t1=time.time()
    print('Label tab_gene_pair is ready to begin')
    tab=get_tab(tab_filename)
    gd=get_gd(gd_pair)
    process_data_in_threads(tab, gd,outfile)

    sort_command = f"sort -k4,4 -t$'\t' --stable -o {outfile} {outfile}"
    subprocess.run(sort_command, shell=True, check=True)

    command=f'cut -f2- {outfile} > sort_{outfile}'
    subprocess.run(command, shell=True, check=True)

    t2=time.time()
    print("Label tab_gene_pairs took "+str(t2-t1)+" second") 
    

