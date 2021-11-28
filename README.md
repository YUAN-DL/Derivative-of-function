# Use the finite difference method and the Fourier spectral method to find the first and second derivatives of functions.
**Compiler Environment: Linux**
### Libraries used:  

#### C++:
 - FFTW3   
 https://www.fftw.org/

#### Python:  
 - Matplotlib   
 https://matplotlib.org/
 - Numpy   
 https://numpy.org/
 
### Compile command：
```
g++ Func_derivative.cpp -o Func_derivative -lfftw3 -lm
```
###  Run command：
###### First, need to run the C++ program to output the test data to files
```
./Func_derivative
```
###### Then use the pyhton image library to draw the image
```
python3 test1_plot1.py
python3 test1_plot2.py
python3 test2_plot1.py
python3 test2_plot2.py
python3 test3_plot1.py
python3 test3_plot2.py
```

