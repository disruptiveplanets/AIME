# Trims the directory so that all -1 are moved to -0.
import os

for d in os.listdir():
    if os.path.isdir(d):
        for f in os.listdir(d):
            if "-1-" in f:
                remove_index = f.find("-1-")
                new_name = f[:remove_index] + "-0-" + f[remove_index+3:]
                print(f"Moving {f} to {new_name}")
                os.rename(f"{d}/{f}", f"{d}/{new_name}")

