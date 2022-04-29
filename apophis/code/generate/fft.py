import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) == 2:
    file_name = sys.argv[-1]
else:
    file_name = "2-params"

plt.figure(figsize=(8, 4))

f = open(f"{file_name}-resolved.dat", 'r')
xs = []
ys = []
zs = []
for line in f.readlines():
    if line == "": continue
    x, y, z = line.split(" ")
    xs.append(float(x))
    ys.append(float(y))
    zs.append(float(z))
time = np.arange(0, len(xs), 1) * 120/3600
xs = np.array(xs) * np.sin(np.pi * time / np.max(time))
ys = np.array(ys) * np.sin(np.pi * time / np.max(time))
zs = np.array(zs) * np.sin(np.pi * time / np.max(time))
period = 2 * np.sin(2 * np.pi /27.38547 * time) * np.sin(2 * np.pi / (264.178/2) * time) * np.sin(np.pi * time / np.max(time))
#period = 2 * np.pi / np.sqrt(xs*xs + ys*ys + zs*zs)

plt.plot(time, xs, label='x')
plt.plot(time, ys, label='y')
#plt.plot(time, zs, label='z')
#plt.plot(time, period, label='period')
plt.xlabel("Time (hr)")
plt.ylabel("spins")
plt.legend()
plt.tight_layout()
    
TRIGGER = 0.3
def get_peaks(freqs, data):
    data = data[np.argsort(freqs)] / np.max(data)
    above_trigger = False
    maxes = []
    for f, d in zip(freqs, data):
        if f < 0:
            continue
        if not above_trigger:
            if d > TRIGGER:
                above_trigger = True
                maxes.append((f, d))
        else:
            if d > maxes[-1][1]:
                maxes[-1] = (f, d)
            if d < TRIGGER:
                above_trigger = False
    maxes = sorted(maxes, key= lambda d: -d[1])
    return [f for (f, d) in maxes[:2]]

def get_envelopes(freqs, data):
    try:
        f1, f2 = get_peaks(freqs, data)
    except Exception:
        return [None]
    return 4 / (f1+f2), 2 / abs(f1-f2)

x_fft = np.abs(np.fft.fft(xs))
y_fft = np.abs(np.fft.fft(ys))
#z_fft = np.abs(np.fft.fft(zs))
p_fft = np.abs(np.fft.fft(period))
freqs = np.fft.fftfreq(len(xs), 120 / 3600)

plt.figure()
mask = (0.1 > np.abs(freqs)) & (np.abs(freqs) > 0.01)
print("\t\t\t", 27.38547, "\t\t", 264.178/2)
print("X sum & diff periods\t", '\t'.join([str(f) for f in get_envelopes(freqs[mask], x_fft[mask])]))
print("Y sum & diff periods\t", '\t'.join([str(f) for f in get_envelopes(freqs[mask], y_fft[mask])]))
#print("Z sum & diff periods\t", '\t'.join([str(f) for f in get_envelopes(freqs[mask], z_fft[mask])]))
print("P sum & diff periods\t", '\t'.join([str(f) for f in get_envelopes(freqs[mask], p_fft[mask])]))
plt.plot(1/freqs[mask], x_fft[mask] / np.max(x_fft[mask]))
plt.plot(1/freqs[mask], y_fft[mask] / np.max(y_fft[mask]))
#plt.plot(1/freqs[mask], z_fft[mask] / np.max(z_fft[mask]))
plt.plot(1/freqs[mask], p_fft[mask] / np.max(p_fft[mask]), color='k')

# plt.axvline(x=264.178 / 2)
# plt.axvline(x=-264.178 / 2)

plt.axvline(x=27.38547, linestyle='dotted', color='k')
plt.axvline(x=-27.38547, linestyle='dotted', color='k')

plt.xlabel("Frequency (1/hr)")
plt.ylabel("Power")
plt.tight_layout()
plt.show()
