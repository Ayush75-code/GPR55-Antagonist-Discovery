"""
AutoDock Vina Setup & Installation
==================================
Sets up environment and installs AutoDock Vina for docking.
Run this first in Google Colab.

Run: Execute cells in Colab notebook
"""

import os
import subprocess


def mount_google_drive():
    """Mount Google Drive for persistent storage"""
    try:
        from google.colab import drive
        print("Mounting Google Drive...")
        drive.mount('/content/drive')
        print("✓ Google Drive mounted!")
        return True
    except ImportError:
        print("Not running in Colab - skipping drive mount")
        return False


def check_system_resources():
    """Display system information"""
    print("=" * 80)
    print("SYSTEM INFORMATION")
    print("=" * 80)
    
    # CPU info
    try:
        result = subprocess.run(['lscpu'], capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if any(x in line for x in ['CPU(s)', 'Model name', 'Thread', 'Core', 'Socket']):
                print(line)
    except:
        print("Could not get CPU info")
    
    print("")
    
    # Memory info
    try:
        result = subprocess.run(['free', '-h'], capture_output=True, text=True)
        print(result.stdout)
    except:
        print("Could not get memory info")


def install_vina():
    """Install AutoDock Vina and dependencies"""
    print("=" * 80)
    print("INSTALLING AUTODOCK VINA")
    print("=" * 80)
    
    commands = [
        "pip install -q pandas matplotlib seaborn numpy",
        "apt-get update -qq",
        "apt-get install -y -qq wget bc",
        "wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/local/bin/vina",
        "chmod +x /usr/local/bin/vina",
        "apt-get install -y -qq parallel"
    ]
    
    for cmd in commands:
        print(f"Running: {cmd[:50]}...")
        os.system(cmd)
    
    # Verify installation
    print("\nVerifying Vina installation:")
    os.system("vina --version")
    print("\n✓ AutoDock Vina installed successfully!")


def check_input_files(work_dir):
    """Check for required input files"""
    print("=" * 80)
    print("CHECKING INPUT FILES")
    print("=" * 80)
    print(f"Location: {work_dir}/\n")
    
    required_files = [
        "GPR55_receptor.pdbqt",
        "AM251_control.pdbqt",
        "pubchem_am251.pdbqt",
        "conf_target_P0.txt",
        "conf_target_P1.txt",
        "conf_target_P2.txt",
        "conf_target_P3.txt",
        "conf_target_P4.txt",
        "conf_target_P5.txt"
    ]
    
    missing_files = []
    found_files = []
    
    for file in required_files:
        filepath = os.path.join(work_dir, file)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            found_files.append(file)
            print(f"✓ {file:<30} ({size:,} bytes)")
        else:
            missing_files.append(file)
            print(f"✗ {file:<30} MISSING!")
    
    print(f"\nFound: {len(found_files)}/{len(required_files)} files")
    
    if missing_files:
        print(f"\n⚠️  ERROR: Missing {len(missing_files)} files!")
        print(f"\nPlease upload these files to {work_dir}/:")
        for f in missing_files:
            print(f"  - {f}")
        return False
    else:
        print("\n✓ All required files found!")
        return True


def main():
    # Configuration
    GDRIVE_BASE = "/content/drive/MyDrive/Major Project/MAIN DATA/CONTROL/Phase1_validation"
    
    # Setup steps
    mount_google_drive()
    check_system_resources()
    install_vina()
    
    if os.path.exists(GDRIVE_BASE):
        check_input_files(GDRIVE_BASE)
    else:
        print(f"\n⚠️  Working directory not found: {GDRIVE_BASE}")
        print("Please create this directory in Google Drive first.")


if __name__ == "__main__":
    main()
