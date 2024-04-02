import sys
import time
import subprocess

def process_blastp_result(input_file, num):
    col_count = {}
    processed_lines = []
    with open(input_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            col1 = cols[0]
            if col1 in col_count:
                col_count[col1] += 1
            else:
                col_count[col1] = 1

    with open(input_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            col1 = cols[0]
            if col_count[col1] == 1:  # If col1 has only one occurrence as key
                processed_lines.append(f"{cols[0]}\t{cols[1]}\tred\n")
            else:
                processed_lines.append(f"{cols[0]}\t{cols[1]}\tNONE\n")

    with open('blastp_label1', 'w') as o:
        o.writelines(processed_lines)



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python label_blastp.py input_file num")
        sys.exit(1)
    input_file = sys.argv[1]
    num = int(sys.argv[2])
    t1=time.time()

    process_blastp_result(input_file,num)

    sort_command = f"sort -k3,3 -t$'\t' --stable -o blastp_label1 blastp_label1"
    subprocess.run(sort_command, shell=True, check=True)

    t2=time.time()
    print("Label blastp took "+str(t2-t1)+" second") 
