import os
import math
import subprocess
import numpy as np

BASE_DIR = "/data/home/bingbing/mywork/lamppost"
EXE_PATH = os.path.join(BASE_DIR, "main")
DATA_DIR = os.path.join(BASE_DIR, "data")
LOG_DIR = "/data/home/bingbing/slurm_logs"
SCRIPTS_DIR = os.path.join(BASE_DIR, "slurm_scripts")

Q_list = [-0.25, 0.0, 0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375]

for d in [DATA_DIR, LOG_DIR, SCRIPTS_DIR]:
    os.makedirs(d, exist_ok=True)

def generate_sbatch(q_idx, q_val, a_idx, a_val):
    r_plus = 1.0 + math.sqrt(max(0.0, 1.0 - a_val**2 - q_val**2))
    hmin0 = 1.099 * r_plus
    hmax0 = 50.0 * r_plus
    heights = np.power(np.arange(100) / 99.0,  2.5) * (hmax0- hmin0) + hmin0

    script_path = os.path.join(SCRIPTS_DIR, f"job_Q{q_idx}_a{a_idx}.sh")
    with open(script_path, "w") as f:
        f.write(f"#!/bin/bash\n")
        f.write(f"#SBATCH --job-name=Q{q_idx}_a{a_idx}\n")
        f.write(f"#SBATCH --partition=cpu\n")
        f.write(f"#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --ntasks=1\n")
        f.write(f"#SBATCH --cpus-per-task=1\n")
        f.write(f"#SBATCH --mem=4G\n")
        f.write(f"#SBATCH --time=5-12:00:00\n")
        f.write(f"#SBATCH --output={LOG_DIR}/slurm_%j_Q{q_idx}_a{a_idx}.out\n")
        f.write(f"#SBATCH --error={LOG_DIR}/slurm_%j_Q{q_idx}_a{a_idx}.err\n")
        f.write(f"#SBATCH --mail-user=bingbinghu@email.ncu.edu.cn\n")
        f.write(f"#SBATCH --mail-type=FAIL\n\n")
        
        f.write(f"cd {BASE_DIR}\n")
        f.write(f'echo "Job started in $(pwd) at $(date)"\n\n')

        # 每个 h 一行调用 main
        for h in heights:
            f.write(f"{EXE_PATH} {a_val:.6f} {h:.6f} {q_val:.6f}\n")

        f.write('\necho "Finished at $(date)"\n')

    return script_path

if __name__ == "__main__":
    for q_idx, q_val in enumerate(Q_list):
        spin_max = math.sqrt(max(0, 1.0 - q_val**2))
        spin_offsets = [0.0, 0.5*(spin_max-0.230), spin_max-0.230, spin_max-0.207,
                        spin_max-0.185, spin_max-0.165, spin_max-0.147, spin_max-0.129,
                        spin_max-0.112, spin_max-0.097, spin_max-0.081, spin_max-0.067,
                        spin_max-0.053, spin_max-0.0265, spin_max-0.014, spin_max-0.0018]
        spin_offsets = [max(0, min(1, a)) for a in spin_offsets]

        for a_idx, a_val in enumerate(spin_offsets):
            path = generate_sbatch(q_idx, q_val, a_idx, a_val)
            subprocess.run(["sbatch", path])

    print(f"任务已提交。总 job 数 = {len(Q_list) * len(spin_offsets)}")
