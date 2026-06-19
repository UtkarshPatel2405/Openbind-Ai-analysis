import os
import subprocess
import sys

def main():
    print("Setting up OpenBind Benchmarking Workspace...")
    
    # 1. Clone the OpenBind EV-A71 2A Benchmark Repository
    if not os.path.exists("EV-A71_2A_benchmark"):
        print("Cloning OpenBind EV-A71 2A Benchmark Repository...")
        subprocess.run(["git", "clone", "https://github.com/OpenBind-Consortium/EV-A71_2A_benchmark.git"], check=True)
    else:
        print("EV-A71_2A_benchmark repository already exists.")

    # 2. Clone the Rowan RBFE Benchmark Repository
    if not os.path.exists("openbind_ev_a71_rbfe"):
        print("Cloning Rowan RBFE Benchmark Repository...")
        subprocess.run(["git", "clone", "https://github.com/rowansci/openbind_ev_a71_rbfe.git"], check=True)
    else:
        print("openbind_ev_a71_rbfe repository already exists.")
        
    print("\nWorkspace setup complete!")
    print("To install dependencies, run: pip install pandas scipy matplotlib numpy requests")

if __name__ == "__main__":
    main()
