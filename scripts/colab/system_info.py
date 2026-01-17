"""
System Hardware Detection for HPC/Cloud Environments
=====================================================
Detects CPU, RAM, GPU specs for optimal parallel processing configuration.
Run: python system_info.py
"""

import platform
import sys
import psutil

try:
    import GPUtil
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False

try:
    from tabulate import tabulate
    TABULATE_AVAILABLE = True
except ImportError:
    TABULATE_AVAILABLE = False


def get_system_info():
    """Collect comprehensive system information"""
    data = [
        ["Python Version", sys.version.split()[0]],
        ["Platform", platform.platform()],
        ["Processor", platform.processor()],
        ["Architecture", str(platform.architecture())],
        ["Machine Type", platform.machine()],
        ["System Name", platform.system()],
        ["Node Name", platform.node()],
        ["Release Version", platform.release()],
        ["CPU Cores (Physical)", psutil.cpu_count(logical=False)],
        ["CPU Cores (Logical)", psutil.cpu_count(logical=True)],
        ["CPU Frequency (Current)", f"{psutil.cpu_freq().current:.2f} MHz" if psutil.cpu_freq() else "N/A"],
        ["CPU Usage", f"{psutil.cpu_percent(interval=1):.2f}%"],
        ["RAM Total", f"{psutil.virtual_memory().total / (1024**3):.2f} GB"],
        ["RAM Available", f"{psutil.virtual_memory().available / (1024**3):.2f} GB"],
        ["Swap Memory", f"{psutil.swap_memory().total / (1024**3):.2f} GB"],
    ]

    # GPU info if available
    if GPU_AVAILABLE:
        gpus = GPUtil.getGPUs()
        if gpus:
            for i, gpu in enumerate(gpus):
                data.append([f"GPU {i} Model", gpu.name])
                data.append([f"GPU {i} Total Memory", f"{gpu.memoryTotal:.2f} MB"])
                data.append([f"GPU {i} Used Memory", f"{gpu.memoryUsed:.2f} MB"])
                data.append([f"GPU {i} Load", f"{gpu.load*100:.2f}%"])
        else:
            data.append(["GPU", "No GPU detected"])
    else:
        data.append(["GPU", "GPUtil not installed"])

    return data


def determine_system_class(logical_cores, total_ram_gb):
    """Determine system performance class for optimization"""
    if logical_cores >= 40 and total_ram_gb >= 150:
        return "HIGH-END WORKSTATION", "~100-140 mol/s", 0.90
    elif logical_cores >= 20 and total_ram_gb >= 60:
        return "HIGH-PERFORMANCE", "~50-80 mol/s", 0.85
    elif logical_cores >= 8 and total_ram_gb >= 16:
        return "STANDARD", "~20-40 mol/s", 0.80
    else:
        return "BASIC", "~5-15 mol/s", 0.70


def main():
    print("\n" + "=" * 80)
    print(" SYSTEM HARDWARE DETECTION")
    print("=" * 80)

    data = get_system_info()

    if TABULATE_AVAILABLE:
        print(tabulate(data, headers=["Detail", "Value"], tablefmt="grid"))
    else:
        for key, value in data:
            print(f"{key:30s}: {value}")

    # Determine optimal configuration
    logical_cores = psutil.cpu_count(logical=True)
    total_ram_gb = psutil.virtual_memory().total / (1024 ** 3)

    system_class, expected_rate, core_usage = determine_system_class(logical_cores, total_ram_gb)

    print("\n" + "=" * 80)
    print(" RECOMMENDED CONFIGURATION")
    print("=" * 80)
    print(f"  System Class: {system_class}")
    print(f"  Expected Processing Rate: {expected_rate}")
    print(f"  Recommended Worker Count: {int(logical_cores * core_usage)}")
    print(f"  Recommended Memory Limit: {int(total_ram_gb * 0.85)} GB")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
