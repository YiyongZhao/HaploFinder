import subprocess
import sys
# Function for running subprocess commands with error handling
def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        raise

def main():
    # Input arguments
    gff1 = sys.argv[1]
    gff2 = sys.argv[2]
    lens1 = sys.argv[3]
    lens2 = sys.argv[4]
    spe1 = sys.argv[5]
    spe2 = sys.argv[6]
    blastp = sys.argv[7]
    synteny = sys.argv[8]
    num = sys.argv[9]
    gf = sys.argv[10]
    imap = sys.argv[11]
    outgroup = sys.argv[12]

    # Step One: Run label_blastp.py
    step_one_command = f"python label_blastp.py {blastp} {num}"
    run_command(step_one_command)

    # Step Two: Run label_synteny.py
    step_two_command = f"python label_synteny.py {synteny}"
    run_command(step_two_command)

    # Step Three: Run label_tree_gd.py
    step_three_command = f"python label_tree_gd.py {gf} {imap} {outgroup}"
    run_command(step_three_command)

    # Combine labeled files into gd_pairs
    cat_command = f"cat blastp_label synteny_label tree_label > total_pairs"
    run_command(cat_command)

    # Step Four: Run dotplot.py for each labeled file
    labeled_files = ["blastp_label", "synteny_label", "tree_label"]
    for file in labeled_files:
        step_four_command = f"python dotplot.py {gff1} {gff2} {lens1} {lens2} {file} {spe1} {spe2}"
        run_command(step_four_command)

    # Final dotplot for gd_pairs
    final_step_four_command = f"python dotplot.py {gff1} {gff2} {lens1} {lens2} total_pairs {spe1} {spe2}"
    run_command(final_step_four_command)

if __name__ == "__main__":
    main()
