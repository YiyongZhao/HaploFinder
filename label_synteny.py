import sys
import time

def process_synteny_result(input_file):
    with open(input_file, 'r') as file:
        processed_lines = [f"{cols[0]}\t{cols[2]}\torange\n" for line in file if not line.startswith('#') for cols in [line.strip().split('\t')]]

    with open('synteny_label', 'w') as o:
        o.writelines(processed_lines)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python label_synteny.py input_file")
        sys.exit(1)
    input_file = sys.argv[1]
    t1=time.time()

    process_synteny_result(input_file)

    t2=time.time()
    print("Label synteny took "+str(t2-t1)+" second") 
