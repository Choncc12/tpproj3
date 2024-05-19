import visual.Debug.pybind11module as lib
import numpy as np
import os


f = 10
N = 1000
t = np.linspace(0, 2 * np.pi, N)

x = np.square(f * t)
X= lib.DFT(x)
Y = np.abs(X)
#lib.plot_line(t,x)

lib.plot_line(np.arange(N),x)
lib.plot_line(np.arange(N),Y)
lib.plot_line(np.arange(N),lib.IDFT(X))

#say_hello(2)

#square(100,10)

#sin(10000,1)

#cos(10000,100)

#test("sine-wave.wav")