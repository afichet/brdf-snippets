#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math

from utia import Utia

def main():
    brdf = Utia('data/m072_fabric107.bin')

    theta_h = math.radians(45)
    samples = 1024
    phi_h = [2 * math.pi * x / samples for x in range(0, samples)]

    reflectance_raw = [brdf.eval_raw(theta_h, phi, theta_h, phi)
                       for phi in phi_h]
    reflectance_interp = [brdf.eval_interpolated(theta_h, phi, theta_h, phi)
                          for phi in phi_h]

    plot_raw = plt.subplot(121, projection='polar')
    plot_raw.plot(phi_h + [phi_h[0]],
                  reflectance_raw + [reflectance_raw[0]])
    plot_raw.grid(True)
    plot_raw.set_title('Raw values')
    
    plot_interp = plt.subplot(122, projection='polar')
    plot_interp.plot(phi_h + [phi_h[0]],
                     reflectance_interp + [reflectance_interp[0]])
    plot_interp.grid(True)
    plot_interp.set_title('Interpolated values')
    
    plt.show()

if __name__ == '__main__':
    main()
