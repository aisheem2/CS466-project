# CS466-project

Hirschberg Algorithm for Local and Fitting Alignment 

# How to run the code
Run the script in the main.py file - and select the type of algorithm ("classic" and "space_efficient") using the `--type` argument. To decide if you want to run local or fitting alignment, choose the `--mode` argument and specify "local" and "fitting". To use one of the sample data texts, the argument `--data` with a number can be used (the default is "8"). Based on these choices, it runs the algorithm we have created in our space-efficient Hirscberg algorithm for fitting and local alignment in `space_efficient.py` or the classic versions in `classic.py` with the data chosen.

# What happens in our code
We consider local and fitting alignments for the classic algorithms and our space-efficient algorithms. The `classic.py` uses the fitting and local alignment algorithms worked on in class and the homework. These are used in benchmarking and to compare memory profiles for our results. In the `space_efficient.py`, we define our `find_end` and `find_start` functions. We then use the hirschberg algorithm for global alignment as defined in `hirschberg.py`, along with our start and ending positions found from the `find_end` and `find_start` functions for v and w, to get the fitting alignment from the function `space_efficient_fitting_align`. The same is done with local alignment, by using the `find_end`, `find_start`, and `hirschberg` functions and we get the final aligned strings in the `space_efficient_local_align` function. 
 
