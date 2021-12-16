# izh_neuralnet
Simulation codes for Izhikevich neural network 

# Prerequisites
* gcc 5.4.0  
* Intel Math Kernel Library 2018 version
* Python 3
* numpy, matplotlib library
* Jupyter Notebook to see the result easily

## Intel MKL library installation
* Need to adjust LD_LIBRARY_PATH
  * export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/lib/intel64

ref: https://www.intel.com/content/www/us/en/developer/articles/guide/intel-math-kernel-library-intel-mkl-2019-install-guide.html

If you cannot install intel mkl library, add macro on the top of the file "lib/izh.h" #define NO_MKL_LIB  
This would slow down the code runtime.

# Run
* make single  
  * Makefile would compile all dependent object files and create run_single.out
* create simulation information json file with same template in the simul_infos/info_single.json
  * "tag" is the prefix of teh output files
* ./run_single.out {simulation_info_json_file_name}
  * e.g. ./run_single.out ./simul_infos/info_single.json
* If the simulation done, output files would be created
  * e.g. if "tag": ./data/test
  * ./data/test_v.dat, ./data/test_u.dat, ./data/test_i.dat, ./data/test_r.dat, ./data/test_t.txt, ./data/test_env.txt
  * .dat files are binary files, so you need to run python code to read and visualize results -> refer result/single_cell_result.ipynb
  * 

# Example


