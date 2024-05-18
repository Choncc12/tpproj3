import build.Debug.pybind11module as lib
import numpy as np
import os

x = np.linspace(0, 10 * np.pi, 1000)
y = np.sin(x)

lib.plot_line(x,y)
#say_hello(2)

#square(100,10)

#sin(10000,1)

#cos(10000,100)

#test("sine-wave.wav")