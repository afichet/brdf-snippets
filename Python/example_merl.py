#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt

from merl import Merl

def main():
    brdf = Merl('blue-metallic-paint.binary')
    samples = 1024

    theta_h = [math.pi/2 * x / samples for x in range(0, samples)]
    reflectance_raw_ndf = [brdf.eval_raw(th, 0, 0) for th in theta_h]
    reflectance_interp_ndf = [brdf.eval_interp(th, 0, 0) for th in theta_h]

    plot_raw_ndf = plt.subplot(121)
    plot_raw_ndf.plot(theta_h, reflectance_raw_ndf)
    plot_raw_ndf.grid(True)
    plot_raw_ndf.set_title("NDF Raw")

    plot_interp_ndf = plt.subplot(122)
    plot_interp_ndf.plot(theta_h, reflectance_interp_ndf)
    plot_interp_ndf.grid(True)
    plot_interp_ndf.set_title("NDF Interpolated")
    
    plt.show()
    
if __name__ == '__main__':
    main()